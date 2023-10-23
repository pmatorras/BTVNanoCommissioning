import collections, numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    muSFs,
    puwei,
    btagSFs,
    load_jmefactory,
    JME_shifts,
    Roccor_shifts,
)
#from BTVNanoCommissioning.utils.testRC2 import getRCFile, applyRC
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_mvatightid
import sys

class NanoProcessor(processor.ProcessorABC):
    '''
    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False, roCorr=False, isSyst=False):
        self._year = year
        self._campaign = campaign
        self.isJERC = isJERC
        self.isCorr = isCorr
        self.roCorr = roCorr
        print(campaign, correction_config[campaign])
        self.lumiMask = load_lumi(self._campaign)#correction_config[self._campaign]["lumiMask"])
    '''
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ## Load corrections
        self.SF_map = load_SF(self._campaign)
        '''
        if isCorr:
            if "BTV" in correction_config[self._campaign].keys():
                self._deepjetc_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepJetC"
                )
                self._deepjetb_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepJetB"
                )
                self._deepcsvc_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepCSVC"
                )
                self._deepcsvb_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepCSVB"
                )
            if "PU" in correction_config[self._campaign].keys():
                self._pu = load_pu(
                    self._campaign, correction_config[self._campaign]["PU"]
                )
        '''
        roCorr =False
        if roCorr:
            print("im doing rocorr")
            rcFilename = getRCFile(self,'2022')
            print(rcFilename)
            #sys.exit()
        self.isSyst = isSyst
        self.lumiMask = load_lumi(self._campaign)

        ## Load corrections
        #if isCorr:
        self.SF_map = load_SF(self._campaign)
        #if isJERC:
        #    self._jet_factory = load_jmefactory(self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        events = missing_branch(events)
        shifts = []
        if "JME" in self.SF_map.keys():
            syst_JERC = True if self.isSyst != None else False
            if self.isSyst == "JERC_split":
                syst_JERC = "split"
            shifts = JME_shifts(
                shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
            )
        else:
            if "Run3" not in self._campaign:
                shifts = [
                    ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
                ]
            else:
                shifts = [
                    (
                        {
                            "Jet": events.Jet,
                            "MET": events.PuppiMET,
                            "Muon": events.Muon,
                        },
                        None,
                    )
                ]
        if "roccor" in self.SF_map.keys():
            shifts = Roccor_shifts(shifts, self.SF_map, events, isRealData, False)
        else:
            shifts[0][0]["Muon"] = events.Muon

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )
    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        self.isRealData = isRealData
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "QCDmuen")
        )

        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }
        #if self.isJERC:
        #    events = add_jec(events, self._campaign, self._jet_factory)
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        ## HLT
        triggers = ["BTagMu_AK4DiJet40_Mu5"]
        checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
        if ak.all(checkHLT == False):
            raise ValueError("HLT paths:", triggers, " are all invalid in", dataset)
        elif ak.any(checkHLT == False):
            print(np.array(triggers)[~checkHLT], " not exist in", dataset)
        trig_arrs = [
            events.HLT[_trig] for _trig in triggers if hasattr(events.HLT, _trig)
        ]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        
        #if self.roCorr:
        #    print("Applying rochester corrections")
        #    events.Muon = applyRC(self,events)
        ## Muon cuts
        muo = events.Muon[(events.Muon.pt > 5) & (abs(events.Muon.eta)<2.5) & mu_idiso(events, self._campaign)]
        req_muon= ak.num(muo.pt)==1
        
        ## Jet cuts
        events.Jet[
           (events.Jet.pt>50)
            & jet_id(events, self._campaign)
            & ak.all(
                events.Jet.metric_table(muo) <= 0.4, axis=2, mask_identity=True
            )
            & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1))
        ]        
        req_jets  = ak.num(events.Jet.pt) >= 2
        events.Jet = ak.pad_none(events.Jet, 1, axis=1)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(muo) <= 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1))
            )
            == 1,
        )

        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        event_level = ak.fill_none(
            req_muon & req_lumi & req_trig & req_jets, False
        )
        if len(events[event_level]) == 0:
            return {dataset: output}
        ####################
        # Selected objects #
        ####################
        sjets = events.Jet[event_level]
        smu   = events.Muon[event_level]
        njet  = ak.count(sjets.pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            spfcands = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx[event_level]
                ]
                .pFCandsIdx
            ]
        sel_jet = sjets[:, 0]
        ####################
        # Weight & Geninfo #
        ####################
        # create Weights object to save individual weights
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            # Load SFs
            if len(self.SF_map.keys()) > 0:
                syst_wei = (
                    True if self.isSyst != None else False
                )  # load systematic flag
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if "MUO" in self.SF_map.keys():
                    muSFs(
                        smu, self.SF_map, weights, syst_wei, False
                    )  # input selected muon
                if "BTV" in self.SF_map.keys():
                    # For BTV weight, you need to specify type
                    btagSFs(sjets, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVC", syst_wei)
        else:
            genflavor = ak.zeros_like(sjets.pt)

        # Systematics information (add name of systematics)
        if shift_name is None:  # weight variations
            systematics = ["nominal"] + list(weights.variations)
        else:  # resolution/ scale variation would use the shift_name
            systematics = [shift_name]
        exclude_btv = [
            "DeepCSVC",
            "DeepCSVB",
            "DeepJetB",
            "DeepJetB",
        ]  # exclude b-tag SFs for btag inputs

        ####################
        #  Fill histogram  #
        ####################
        for syst in systematics:
            if self.isSyst == None and syst != "nominal":
                break
            if self.noHist:
                break
            weight = (
                weights.weight()
                if syst == "nominal" or syst == shift_name
                else weights.weight(modifier=syst)
            )  # shift up/down for weight systematics

            # fill the histogram (check axis defintion in histogrammer and following the order)

            for histname, h in output.items():
                print ("histname",histname)
                if (
                    "Deep" in histname
                    and "btag" not in histname
                    and histname in events.Jet.fields
                ):
                    h.fill(
                        flatten(genflavor),
                        flatten(sjets[histname]),
                        weight=flatten(
                            ak.broadcast_arrays(weights.weight(), sjets["pt"])[0]
                        ),
                    )
                elif (
                    "PFCands" in events.fields
                    and "PFCands" in histname
                    and histname.split("_")[1] in events.PFCands.fields
                ):
                    h.fill(
                        flatten(ak.broadcast_arrays(genflavor[:, 0], spfcands["pt"])[0]),
                        flatten(spfcands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(weights.weight(), spfcands["pt"])[0]
                        ),
                    )

                elif "jet" in histname and "dr" not in histname and "njet" != histname:
                    for i in range(2):
                        print ("name", histname, i)
                        sel_jet = sjets[:, i]
                        print("check this", str(i), histname, len(genflavor[:, i]), len(sel_jet[histname.replace(f"jet{i}_", "")]))
                        if str(i) in histname:
                            try:
                                h.fill(
                                flatten(genflavor[:, i]),
                                flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                                weight=weight,
                                                )
                            except ValueError as ve:
                                print("There is a value error here:", str(i), histname, len(genflavor[:, i]), len(sel_jet[histname.replace(f"jet{i}_", "")]), len(weights.weight()), ve)

                                #exit()
                elif (
                    "btagDeep" in histname
                    and "0" in histname
                    and histname.replace("_0", "") in events.Jet.fields
                ):
                    h.fill(
                        flav=genflavor[:, 0],
                        syst="noSF",
                        discr=np.where(
                            sel_jet[histname.replace("_0", "")] < 0,
                            -0.2,
                            sel_jet[histname.replace("_0", "")],
                        ),
                        weight=weights.weight(),
                    )
                    if (
                        not isRealData
                        and self.isCorr
                        and "btag" in self.SF_map.keys()
                        and "_b" not in histname
                        and "_bb" not in histname
                        and "_lepb" not in histname
                    ):
                        for syst in disc_list[histname.replace("_0", "")][0].keys():
                            h.fill(
                                flav=genflavor[:, 0],
                                syst=syst,
                                discr=np.where(
                                    sel_jet[histname.replace("_0", "")] < 0,
                                    -0.2,
                                    sel_jet[histname.replace("_0", "")],
                                ),
                                weight=weights.weight()
                                * disc_list[histname.replace("_0", "")][0][syst],
                            )
            output["njet"].fill(njet, weight=weights.weight())
            output["dr_mumu"].fill(dr=snegmu.delta_r(sposmu), weight=weights.weight())
            output["z_pt"].fill(flatten(sz.pt), weight=weights.weight())
            output["z_eta"].fill(flatten(sz.eta), weight=weights.weight())
            output["z_phi"].fill(flatten(sz.phi), weight=weights.weight())
            output["z_mass"].fill(flatten(sz.mass), weight=weights.weight())
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev.Jet = sel_jet
            pruned_ev.Muon = smu
            pruned_ev["dilep"] = sposmu + snegmu
            pruned_ev["dilep", "pt"] = pruned_ev.dilep.pt
            pruned_ev["dilep", "eta"] = pruned_ev.dilep.eta
            pruned_ev["dilep", "phi"] = pruned_ev.dilep.phi
            pruned_ev["dilep", "mass"] = pruned_ev.dilep.mass
            if "PFCands" in events.fields:
                pruned_ev.PFCands = spfcands
            # Add custom variables
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )

            for i in range(2):
                output[f"dr_mujet{i}"].fill(
                    flav=flatten(genflavor[:, i]),
                    dr=flatten(smu.delta_r(sjets[:, i])),
                    weight=weights.weight(),
                )
        
            # Create a list of variables want to store. For objects from the PFNano file, specify as {object}_{variable}, wildcard option only accepted at the end of the string
            out_branch = np.setdiff1d(
                np.array(pruned_ev.fields), np.array(events.fields)
            )
            out_branch = np.delete(
                out_branch,
                np.where((out_branch == "dilep")),
            )

            for kin in ["pt", "eta", "phi", "mass", "pfRelIso04_all", "dxy", "dz"]:
                for obj in ["Muon", "Jet", "dilep"]:
                    if (obj != "Muon") and ("pfRelIso04_all" == kin or "d" in kin):
                        continue
                    out_branch = np.append(out_branch, [f"{obj}_{kin}"])
            out_branch = np.append(
                out_branch, ["Jet_btagDeep*", "Jet_DeepJet*", "PFCands_*"]
            )
            # write to root files
            os.system(f"mkdir -p {self.name}/{dataset}")
            with uproot.recreate(
                f"{self.name}/{dataset}/f{events.metadata['filename'].split('_')[-1].replace('.root','')}_{systematics[0]}_{int(events.metadata['entrystop']/self.chunksize)}.root"
            ) as fout:
                fout["Events"] = uproot_writeable(pruned_ev, include=out_branch)

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
