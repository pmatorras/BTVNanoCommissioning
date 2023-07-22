import collections, numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    muSFs,
    puwei,
    btagSFs,
    load_jmefactory,
)
from BTVNanoCommissioning.utils.testRC2 import getRCFile, applyRC
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_mvatightid
import sys

class NanoProcessor(processor.ProcessorABC):
    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False, roCorr=False, isSyst=False):
        self._year = year
        self._campaign = campaign
        self.isJERC = isJERC
        self.isCorr = isCorr
        self.roCorr = roCorr
        print(campaign, correction_config[campaign])
        self.lumiMask = load_lumi(self._campaign)#correction_config[self._campaign]["lumiMask"])

        ## Load corrections
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

        if roCorr:
            print("im doing rocooricos")
            rcFilename = getRCFile(self,'2022')
            print(rcFilename)
            #sys.exit()
        self.isSyst = isSyst
        self.lumiMask = load_lumi(self._campaign)

        ## Load corrections
        if isCorr:
            self.SF_map = load_SF(self._campaign)
        if isJERC:
            self._jet_factory = load_jmefactory(self._campaign)
        _hist_event_dict = histogrammer("QCDmuen")
        self.make_output = lambda: {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.make_output()
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        self.isRealData = isRealData
        events = missing_branch(events)

        if self.isJERC:
            events = add_jec(events, self._campaign, self._jet_factory)
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
        #triggers = ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"]
        #triggers = ["HLT_BTagMu_AK4DiJet40_Mu5"]
        triggers = ["BTagMu_AK4DiJet40_Mu5"]
        print("---------------->all triggers:", triggers)
        #sys.exit()
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

        if self.roCorr:
            print("Applying rochester corrections")
            events.Muon = applyRC(self,events)
        ## Muon cuts
        dilep_muo = events.Muon[(events.Muon.pt > 5) & (abs(events.Muon.eta)<2.5) & mu_idiso(events, self._campaign)]
        ## Electron cuts
        dilep_ele = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]
        ## dilepton
        '''
        req_dilep = (
            (ak.num(pos_dilep.pt) >= 1)
            & (ak.num(neg_dilep.pt) >= 1)
            & (ak.num(dilep_mu.charge) >= 2)
            & (ak.num(dilep_ele.charge) == 0)
        )
        pos_dilep = ak.pad_none(pos_dilep, 1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep, 1, axis=1)

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81) & (dilep_mass.mass < 101) & (dilep_mass.pt > 15)
        )
        '''
        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & ak.all(
                events.Jet.metric_table(dilep_muo) <= 0.4, axis=2, mask_identity=True
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 50
        event_jet = ak.pad_none(event_jet, 1, axis=1)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(dilep_muo) <= 0.4,
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
            req_lumi & req_trig & req_jets, False
        )
        if len(events[event_level]) == 0:
            return {dataset: output}
        ####################
        # Selected objects #
        ####################
        sjets = event_jet[event_level]
        njet = ak.count(sjets.pt, axis=1)
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
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
        if not isRealData and self.isCorr:
            if "PU" in self.SF_map.keys():
                weights.add(
                    "puweight", puwei(self.SF_map, events[event_level].Pileup.nTrueInt)
                )
            if "MUO" in self.SF_map.keys() or "EGM" in self.SF_map.keys():
                weights.add("lep1sf", muSFs(sposmu, self.SF_map, True))
                weights.add("lep2sf", muSFs(snegmu, self.SF_map, True))
        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            if self.isCorr and (
                "btag" in self.SF_map.keys() or "ctag" in self.SF_map.keys()
            ):
                jetsfs_c = collections.defaultdict(dict)
                jetsfs_b = collections.defaultdict(dict)
                csvsfs_c = collections.defaultdict(dict)
                csvsfs_b = collections.defaultdict(dict)

                ## for each jet
                if self.isCorr and (
                    "btag" in self.SF_map.keys() or "ctag" in self.SF_map.keys()
                ):
                    jetsfs_c[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepJetC")
                    jetsfs_b[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepJetB")
                    csvsfs_c[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepCSVC")
                    csvsfs_b[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepCSVB")
                    if self.isSyst:
                        for syst in [
                            "hf",
                            "lf",
                            "cferr1",
                            "cferr2",
                            "hfstat1",
                            "hfstat2",
                            "lfstats1",
                            "lfstats2",
                        ]:
                            jetsfs_c[0][f"SF_{syst}_up"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepJetC", f"up_{syst}"
                            )
                            jetsfs_c[0][f"SF_{syst}_dn"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepJetC", f"down_{syst}"
                            )
                            csvsfs_c[0][f"SF_{syst}_up"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepCSVC", f"up_{syst}"
                            )
                            csvsfs_c[0][f"SF_{syst}_dn"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepCSVC", f"down_{syst}"
                            )
                        csvsfs_b[0][f"SF_{syst}_up"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepCSVB", f"up"
                        )
                        csvsfs_b[0][f"SF_{syst}_dn"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepCSVB", f"down"
                        )
                        jetsfs_b[0][f"SF_{syst}_up"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepJetB", f"up"
                        )
                        jetsfs_b[0][f"SF_{syst}_dn"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepJetB", f"down"
                        )

                disc_list = {
                    "btagDeepB": csvsfs_b,
                    "btagDeepC": csvsfs_b,
                    "btagDeepFlavB": jetsfs_b,
                    "btagDeepFlavC": jetsfs_b,
                    "btagDeepCvL": csvsfs_c,
                    "btagDeepCvB": csvsfs_c,
                    "btagDeepFlavCvL": jetsfs_c,
                    "btagDeepFlavCvB": jetsfs_c,
                }
        ####################
        #  Fill histogram  #
        ####################
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
                print("inside jet", histname)
                for i in range(2):
                    sel_jet = sjets[:, i]
                    print( str(i), histname, len(genflavor[:, i]), len(sel_jet[histname.replace(f"jet{i}_", "")]))
                    if str(i) in histname:
                        try:
                            h.fill(
                            flatten(genflavor[:, i]),
                            flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                            weight=weights.weight(),
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
        for i in range(2):
            output[f"dr_mujet{i}"].fill(
                flav=flatten(genflavor[:, i]),
                dr=flatten(smu.delta_r(sjets[:, i])),
                weight=weights.weight(),
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
