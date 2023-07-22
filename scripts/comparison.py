import numpy as np
import argparse, os, arrow, glob
from coffea.util import load
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
import hist
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.helpers.definitions import definitions, axes_name
from BTVNanoCommissioning.utils.plot_utils import (
    plotratio,
    markers,
    autoranger,
    rebin_hist,
)

plt.style.use(hep.style.ROOT)
from BTVNanoCommissioning.helpers.xs_scaler import collate

bininfo = definitions()
parser = argparse.ArgumentParser(description="make comparison for different campaigns")
parser.add_argument(
    "-p",
    "--phase",
    required=True,
    choices=list(workflows.keys()),
    dest="phase",
    help="which phase space",
)
parser.add_argument(
    "-i",
    "--input1",
    required=True,
    type=str,
    help="input coffea files (str), splitted different files with ','. Wildcard option * available as well.",
)
parser.add_argument(
    "--input2",
    default = None,
    type=str,
    help="input coffea files for comparison (str), splitted different files with ','. Wildcard option * available as well.",
)


parser.add_argument("-r", "--ref", required=True, help="referance dataset")
parser.add_argument(
    "-c", "--compared", required=True, type=str, help="compared datasets, splitted by ,"
)
parser.add_argument(
    "--sepflav", action="store_true", help="seperate flavour(b/c/light)"
)
parser.add_argument("--log", action="store_true", help="log on y axis")
parser.add_argument(
    "--norm",
    action="store_true",
    help="compare shape, normalized yield to reference",
    default=False,
)
parser.add_argument(
    "-v",
    "--variable",
    type=str,
    help="variables to plot, splitted by ,. Wildcard option * available as well. Specifying `all` will run through all variables.",
)
parser.add_argument(
    "--xrange", type=str, default=None, help="custom x-range, --xrange xmin,xmax"
)
parser.add_argument(
    "--flow",
    type=str,
    default=None,
    help="str, optional {None, 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.",
)
parser.add_argument("--ext", type=str, default="", help="prefix name/btv name tag")
parser.add_argument("--com", default="13", type=str, help="sqrt(s) in TeV")
parser.add_argument(
    "--shortref",
    default="",
    type=str,
    help="short name for reference dataset for legend",
)
parser.add_argument(
    "--shortcomp",
    default="",
    type=str,
    help="short names for compared datasets for legend, split by ','",
)
parser.add_argument(
    "--autorebin",
    default=None,
    help="Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)",
)
parser.add_argument(
    "--xlabel",
    type=str,
    help="rename the label of variables to plot, splitted by ,.  Wildcard option * NOT available here",
)
parser.add_argument(
    "--ylabel", type=str, default=None, help="Modify y-axis label of plot"
)
parser.add_argument(
    "--splitOSSS",
    type=int,
    default=None,
    help="Only for W+c phase space, split opposite sign(1) and same sign events(-1), if not specified, the combined OS-SS phase space is used",
)
args = parser.parse_args()
<<<<<<< HEAD:plotting/comparison.py
def getOutput(args_input):
    output_ = {}
    print("args_input", args_input)
    if len(args_input.split(",")) > 1:
        output_ = {i: load(i) for i in args_input.split(",")}
    elif "*" in args_input:
        files = glob.glob(args_input)
        output_ = {i: load(i) for i in files}
    else:
        output_ = load(args_input)
    return output_
print ("input", args.input1)
print (args.input1.split(','))
output  = getOutput(args.input1)
#    output = {}
#    if len(args.input.split(",")) > 1:
#        output = {i: load(i) for i in args.input.split(",")}
#    elif "*" in args.input:
#        files = glob.glob(args.input)
#        output = {i: load(i) for i in files}
#    else:
#        output = load(args.input)

=======
output = {}
if len(args.input.split(",")) > 1:
    output = {i: load(i) for i in args.input.split(",")}
elif "*" in args.input:
    files = glob.glob(args.input)
    output = {i: load(i) for i in files}
else:
    output = {args.input: load(args.input)}
>>>>>>> upstream/master:scripts/comparison.py
mergemap = {}
time = arrow.now().format("YY_MM_DD")
if not os.path.isdir(f"plot/BTV/{args.phase}_{args.ext}_{time}/"):
    os.makedirs(f"plot/BTV/{args.phase}_{args.ext}_{time}/")

def createCollated(output1, output2, areSame=True):
    if not any(".coffea" in o for o in output1.keys()):
        mergemap[args.ref] = [m for m in output1.keys() if args.ref == m]
        for c in args.compared.split(","):
            mergemap[c] = [m for m in output2.keys() if c == m]
    else:
        reflist = []
        for f in output1.keys():
            reflist.extend([m for m in output1[f].keys() if args.ref == m])
        mergemap[args.ref] = reflist

        for c in args.compared.split(","):
            comparelist = []
            for f in output2.keys():
                comparelist.extend([m for m in output2[f].keys() if c == m])
            mergemap[c] = comparelist
    if areSame:
        return collate(output1, mergemap)
    else: return collate(output1, mergemap), collate(output2, mergemap)
    
#collated = collate(output, mergemap)
if args.input2:
    output2 = getOutput(args.input2)
    collated, compareCollated = createCollated(output, output2, False)
else:
    collated = createCollated(output, output)
# if not any(".coffea" in o for o in output.keys()):
#     mergemap[args.ref] = [m for m in output.keys() if args.ref == m]
#     for c in args.compared.split(","):
#         mergemap[c] = [m for m in output.keys() if c == m]
# else:
#     reflist = []
#     for f in output.keys():
#         reflist.extend([m for m in output[f].keys() if args.ref == m])
#     mergemap[args.ref] = reflist

<<<<<<< HEAD:plotting/comparison.py
#     for c in args.compared.split(","):
#         comparelist = []
#         for f in output.keys():
#             comparelist.extend([m for m in output[f].keys() if c == m])
#         mergemap[c] = comparelist
# collated = collate(output, mergemap)
#print ("collated", collated)
=======
    for c in args.compared.split(","):
        comparelist = []
        for f in output.keys():
            comparelist.extend([m for m in output[f].keys() if c == m])
        mergemap[c] = comparelist
print(output.keys(), mergemap.keys())
collated = collate(output, mergemap)
>>>>>>> upstream/master:scripts/comparison.py
### style settings
if "Run" in args.ref:
    hist_type = "errorbar"
    label = "Preliminary"
else:
    hist_type = "step"
    label = "Simulation Preliminary"
nj = 1
### input text settings
if "Wc" in args.phase:
    input_txt = "W+c"
    if args.splitOSSS == 1:
        input_txt = input_txt + " OS"
    elif args.splitOSSS == -1:
        input_txt = input_txt + " SS"
    else:
        input_txt = input_txt + " OS-SS"
elif "DY" in args.phase:
    input_txt = "DY+jets"
elif "semilep" in args.phase:
    input_txt = r"t$\bar{t}$ semileptonic"
    nj = 4
elif "dilep" in args.phase:
    input_txt = r"t$\bar{t}$ dileptonic"
    nj = 2
if "emctag" in args.phase:
    input_txt = input_txt + " (e$\mu$)"
elif "ectag" in args.phase:
    input_txt = input_txt + " (e)"
elif "ttdilep_sf" == args.phase:
    input_txt = input_txt + " (e$\mu$)"
else:
    input_txt = input_txt + " ($\mu$)"
if "ctag" in args.phase and "DY" not in args.phase:
    input_txt = input_txt + "\nw/ soft-$\mu$"
if args.shortref == "":
    args.shortref  = args.input1 if args.input2 else args.ref

if args.shortcomp == "":
    
    args.shortcomp = args.input2 if args.input2 else args.compared

if args.variable == "all":
    var_set = collated[args.ref].keys()
elif "*" in args.variable:
    var_set = []
    print(len(args.variable.split('*')), "\n")
    for var in collated[args.ref].keys():
        argsplit = args.variable.split('*')
        ishere = False 

        

        if len(argsplit[-1]) == 0 and argsplit[0] in var: var_set.append(var) 
        elif len(argsplit)>0 and  len(argsplit[-1]) >0:
            if argsplit[0] in var and argsplit[1] in var.split(argsplit[0])[1]: var_set.append(var)
    
else:
    var_set = args.variable.split(",")
<<<<<<< HEAD:plotting/comparison.py
#print ("var", var_set)
#print(len(args.variable.split('*')), args.variable.split('*'))
#exit()
def makePlotting(collated1, collated2, input2text):
    for index, discr in enumerate(var_set):
        allaxis = {}
        if "flav" in collated1[args.ref][discr].axes.name:
            allaxis["flav"] = sum
        if "syst" in collated1[args.ref][discr].axes.name:
            allaxis["syst"] = sum
        if "osss" in collated1[args.ref][discr].axes.name:  ## do dominal OS-SS
            if args.splitOSSS is None:  # OS-SS
                collated1[args.ref][discr] = (
                    collated1[args.ref][discr][{"osss": 0}]
                    + collated1[args.ref][discr][{"osss": 1}] * -1
                )
                for c in args.compared.split(","):
                    collated2[c][discr] = (
                        collated2[c][discr][{"osss": 0}]
                        + collated2[c][discr][{"osss": 1}] * -1
                    )
            elif args.splitOSSS == 1:
                allaxis["osss"] = 0  # opposite sign
            elif args.splitOSSS == -1:
                allaxis["osss"] = 1  # same sign
        do_xerr = False
        if args.autorebin is not None:
            if args.autorebin.isdigit():
                rebin = int(args.autorebin)
            else:
                rebin = np.array([float(i) for i in args.autorebin.split(",")])
                do_xerr = True
            collated["mc"][discr] = rebin_hist(
                collated["mc"][discr], collated["mc"][discr].axes[-1].name, rebin
            )
            collated["data"][discr] = rebin_hist(
                collated["data"][discr], collated["data"][discr].axes[-1].name, rebin
=======
for index, discr in enumerate(var_set):
    allaxis = {}
    if "flav" in collated[args.ref][discr].axes.name:
        allaxis["flav"] = sum
    if "syst" in collated[args.ref][discr].axes.name:
        allaxis["syst"] = "nominal"
    if "osss" in collated[args.ref][discr].axes.name:  ## do dominal OS-SS
        if args.splitOSSS is None:  # OS-SS
            collated[args.ref][discr] = (
                collated[args.ref][discr][{"osss": 0}]
                + collated[args.ref][discr][{"osss": 1}] * -1
>>>>>>> upstream/master:scripts/comparison.py
            )
        # FIXME: add wildcard option for xlabel
        # xlabel = (
        #     args.xlabel
        #     if args.xlabel is not None
        #     else collated["data"][discr].axes[-1].label
        # )  # Use label from stored hists
        ## FIXME: Set temporary fix for the x-axis
        if args.xlabel is not None:
            args.xlabel.split(",")[index]
        elif "DeepJet" in discr or "DeepCSV" in discr or "PFCands" in discr:
            xlabel = (
                bininfo[discr]["displayname"]
                + " ["
                + bininfo[discr]["inputVar_units"]
                + "]"
                if bininfo[discr]["inputVar_units"] is not None
                else bininfo[discr]["displayname"]
            )
        else:
            xlabel = axes_name(discr)

        text = args.ext
        if args.norm:
            text = args.ext + "\nNormalized to Ref."
            for c in args.compared.split(","):
                collated2[c][discr] = collated2[c][discr] * float(
                    np.sum(collated1[args.ref][discr][allaxis].values())
                    / np.sum(collated2[c][discr][allaxis].values())
                )

        if args.sepflav:  # split into 3 panels for different flavor
            fig, (ax, rax, rax2, rax3) = plt.subplots(
                4,
                1,
                figsize=(8, 8),
                gridspec_kw={"height_ratios": (3, 1, 1, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            ax.set_xlabel(None)
            hep.cms.label(label, com=args.com, data=True, loc=0, ax=ax)
            laxis = {"flav": 0}
            puaxis = {"flav": 1}
            caxis = {"flav": 2}
            baxis = {"flav": 3}

            if "syst" in collated1[args.ref][discr].axes.name:
                laxis["syst"] = "noSF"
                puaxis["syst"] = "noSF"
                caxis["syst"] = "noSF"
                baxis["syst"] = "noSF"

            hep.histplot(
                collated1[args.ref][discr][laxis] + collated1[args.ref][discr][puaxis],
                label=args.shortref + "-l",
                color="b",
                histtype=hist_type,
                yerr=True,
                ax=ax,
                xerr=do_xerr
                # flow=args.flow,
            )
            hep.histplot(
                collated1[args.ref][discr][caxis],
                label=args.shortref + "-c",
                color="g",
                histtype=hist_type,
                yerr=True,
                ax=ax,
                xerr=do_xerr
                # flow=args.flow,
            )
            hep.histplot(
                collated1[args.ref][discr][baxis],
                label=args.shortref + "-b",
                yerr=True,
                color="r",
                histtype=hist_type,
                ax=ax,
                xerr=do_xerr
                # flow=args.flow,
            )

            mindex = 0
            ax.legend(ncol=3, loc=1)
            for c, s in zip(args.compared.split(","), args.shortcomp.split(",")):
                hep.histplot(
                    collated2[c][discr][laxis] + collated2[c][discr][puaxis],
                    label=s + "-l",
                    color="b",
                    marker=markers[mindex + 1],
                    histtype="errorbar",
                    yerr=True,
                    ax=ax,
                    xerr=do_xerr
                    # flow=args.flow,
                )
                hep.histplot(
                    collated2[c][discr][caxis],
                    label=s + "-c",
                    color="g",
                    marker=markers[mindex + 1],
                    histtype="errorbar",
                    yerr=True,
                    ax=ax,
                    xerr=do_xerr
                    # flow=args.flow,
                )
                hep.histplot(
                    collated2[c][discr][baxis],
                    label=s + "-b",
                    yerr=True,
                    color="r",
                    marker=markers[mindex + 1],
                    histtype="errorbar",
                    ax=ax,
                    xerr=do_xerr
                    # flow=args.flow,
                )
                # comparison splitted by flavor
                rax = plotratio(
                    collated2[c][discr][laxis] + collated2[c][discr][puaxis],
                    collated1[args.ref][discr][laxis] + collated1[args.ref][discr][puaxis],
                    ax=rax,
                    denom_fill_opts=None,
                    error_opts={"color": "b", "marker": markers[mindex + 1]},
                    clear=False,
                    xerr=do_xerr
                    # flow=args.flow,
                )
                rax2 = plotratio(
                    collated2[c][discr][caxis],
                    collated1[args.ref][discr][caxis],
                    ax=rax2,
                    denom_fill_opts=None,
                    error_opts={"color": "g", "marker": markers[mindex + 1]},
                    clear=False,
                    xerr=do_xerr
                    # flow=args.flow,
                )
                rax3 = plotratio(
                    collated2[c][discr][baxis],
                    collated1[args.ref][discr][baxis],
                    ax=rax3,
                    denom_fill_opts=None,
                    error_opts={"color": "r", "marker": markers[mindex + 1]},
                    clear=False,
                    xerr=do_xerr
                    # flow=args.flow,
                )
                mindex += 1

            discrs = discr
            ax.set_xlabel("A.U.")
            rax.set_ylabel("udsg-jets")
            rax2.set_ylabel("c-jets")
            rax3.set_ylabel("b-jets")
            rax.set_ylim(0.5, 1.5)
            rax2.set_ylim(0.5, 1.5)
            rax3.set_ylim(0.5, 1.5)

            rax3.set_xlabel(xlabel)
            ax.legend()

            at = AnchoredText(input_txt + "\n" + text, loc=2, frameon=False)
            ax.add_artist(at)
            hep.mpl_magic(ax=ax)
            if args.log:
                ax.set_yscale("log")
            fig.savefig(
                f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_inclusive{discrs}.png"
            )
            fig.savefig(
                f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_inclusive{discrs}.pdf"
            )

        else:
            fig, ((ax), (rax)) = plt.subplots(
                2, 1, figsize=(12, 12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            hep.cms.label(label, com=args.com, data=True, loc=0, ax=ax)
            ax.set_xlabel(None)
            hep.histplot(
                collated1[args.ref][discr][allaxis],
                label=args.shortref + " (Ref)",
                histtype=hist_type,
                yerr=True,
                ax=ax,
                xerr=do_xerr
                # flow=args.flow,
            )
<<<<<<< HEAD:plotting/comparison.py
            for c, s in zip(args.compared.split(","), args.shortcomp.split(",")):
                hep.histplot(
                    collated2[c][discr][allaxis],
                    label=s,
                    histtype=hist_type,
                    yerr=True,
                    ax=ax,
                    xerr=do_xerr
                    # flow=args.flow,
                )
=======
        for i, c in enumerate(args.compared.split(",")):
            plotratio(
                collated[c][discr][allaxis],
                collated[args.ref][discr][allaxis],
                ax=rax,
                denom_fill_opts=None,
                error_opts={"color": ax.get_lines()[i * 2 + 1].get_color()},
                clear=False,
                xerr=do_xerr
                # flow=args.flow,
            )  ## No error band used
        alls = collated[args.ref][discr][allaxis]
        for c in args.compared.split(","):
            alls = collated[c][discr][allaxis] + alls
        xmin, xmax = autoranger(alls)
        rax.set_xlim(xmin, xmax)
        if args.xrange is not None:
            xmin, xmax = float(args.xrange.split(",")[0]), float(
                args.xrange.split(",")[1]
            )
        rax.set_xlabel(xlabel)
        ax.set_xlabel(None)
        ax.set_ylabel(args.ylabel)
        rax.set_ylabel("Other/Ref")
        ax.legend(loc=1)
        rax.set_ylim(0.0, 2.0)
>>>>>>> upstream/master:scripts/comparison.py

            for i, c in enumerate(args.compared.split(",")):
                plotratio(
                    collated2[c][discr][allaxis],
                    collated1[args.ref][discr][allaxis],
                    ax=rax,
                    denom_fill_opts=None,
                    error_opts={"color": ax.get_lines()[i + 1].get_color()},
                    clear=False,
                    xerr=do_xerr
                    # flow=args.flow,
                )  ## No error band used
            alls = collated1[args.ref][discr][allaxis]
            for c in args.compared.split(","):
                alls = collated2[c][discr][allaxis] + alls
            xmin, xmax = autoranger(alls)
            rax.set_xlim(xmin, xmax)
            if args.xrange is not None:
                xmin, xmax = float(args.xrange.split(",")[0]), float(
                    args.xrange.split(",")[1]
                )
            rax.set_xlabel(xlabel)
            ax.set_xlabel(None)
            ax.set_ylabel(args.ylabel)
            rax.set_ylabel("Other/Ref")
            ax.legend(loc=1)
            rax.set_ylim(0.0, 2.0)

            at = AnchoredText(input_txt + "\n" + args.ext, loc=2, frameon=False)
            ax.add_artist(at)
            hep.mpl_magic(ax=ax)
            ax.set_ylim(bottom=0)
            varext = "inclusive" if args.ref=="all" else args.ref.split('_')[0]
            logext = ""
            normtext = ""
            print (args.ref, args.ref.split('_')[0])
            if args.norm:
                normtext = "_norm"
            if args.log:
                ax.set_yscale("log")
                logext += "_log"
                ax.set_ylim(bottom=0.1)
                hep.mpl_magic(ax=ax)
            print(
                "creating:",
                f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_{varext}_{discr}{logext}{normtext}{input2text}.png",
            )
            fig.savefig(
                f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_{varext}_{discr}{logext}{normtext}{input2text}.pdf"
            )
            fig.savefig(
                f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_{varext}_{discr}{logext}{normtext}{input2text}.png"
            )


if args.input2:
    print("i should be comparing")
    input2text = '_'+'compareRC'# Might need to be generalised in the future
    print(args.input1, args.input2, input2text)
    makePlotting(collated,compareCollated, input2text)
else: 
    input2text = ''
    makePlotting(collated, collated, input2text)
