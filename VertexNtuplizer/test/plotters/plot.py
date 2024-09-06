import ROOT
ROOT.gROOT.SetBatch(1)
import time

from plotConfig import *
from plotFunctions import *

startScript = time.time()

file_TTnoPU = ROOT.TFile.Open(fname_TTnoPU, "READ")
file_TTPU200 = ROOT.TFile.Open(fname_TTPU200, "READ")
dir_TTnoPU = file_TTnoPU.Get("vertexNtuplizer")
dir_TTPU200 = file_TTPU200.Get("vertexNtuplizer")

n = {}
p = {}

refHist = "nGPs"
nevents_TTnoPU = dir_TTnoPU.Get(refHist).GetEntries()
nevents_TTPU200 = dir_TTPU200.Get(refHist).GetEntries()
ratio = nevents_TTnoPU / nevents_TTPU200
print("There are %d TTnoPU events and %d TTPU200 events. The ratio is %.2f." % (nevents_TTnoPU, nevents_TTPU200, ratio))

# Plot PF candidate times
start = time.time()
h1 = dir_TTnoPU.Get("pfCandidate_bs_tval").Clone()
h2 = dir_TTnoPU.Get("pfCandidate_bs_terr").Clone()
h3 = dir_TTnoPU.Get("pfCandidate_bs_tsig").Clone()
h4 = dir_TTnoPU.Get("pfCandidate_pv_tval").Clone()
h5 = dir_TTnoPU.Get("pfCandidate_pv_terr").Clone()
h6 = dir_TTnoPU.Get("pfCandidate_pv_tsig").Clone()
h7 = dir_TTPU200.Get("pfCandidate_bs_tval").Clone()
h8 = dir_TTPU200.Get("pfCandidate_bs_terr").Clone()
h9 = dir_TTPU200.Get("pfCandidate_bs_tsig").Clone()
h10 = dir_TTPU200.Get("pfCandidate_pv_tval").Clone()
h11 = dir_TTPU200.Get("pfCandidate_pv_terr").Clone()
h12 = dir_TTPU200.Get("pfCandidate_pv_tsig").Clone()
h1.Scale(1.0 / h1.Integral())
h2.Scale(1.0 / h2.Integral())
h3.Scale(1.0 / h3.Integral())
h4.Scale(1.0 / h4.Integral())
h5.Scale(1.0 / h5.Integral())
h6.Scale(1.0 / h6.Integral())
h7.Scale(1.0 / h7.Integral())
h8.Scale(1.0 / h8.Integral())
h9.Scale(1.0 / h9.Integral())
h10.Scale(1.0 / h10.Integral())
h11.Scale(1.0 / h11.Integral())
h12.Scale(1.0 / h12.Integral())
labels = [ttnoPU+", BS extrapolation", ttnoPU+", PV extrapolation", ttPU200+", BS extrapolation", ttPU200+", PV extrapolation"]
ratio_labels = ["BS extrapolation", "PV extrapolation"]
colors = [COLORS[0], COLORS[2], COLORS[0], COLORS[2]]
styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
hists = [h1, h4, h7, h10]
plot(hists, labels, colors, styles, "Time [ps]", "Normalized PF candidates", save_path+"Tracks", "tval")
plot_ratio(hists, labels, 2, ratio_labels, colors, styles, "Time [ps]", "Normalized PF candidates", "PU0/PU200", save_path+"Tracks", "tval")
hists = [h2, h5, h8, h11]
plot(hists, labels, colors, styles, "Time error [ps]", "Normalized PF candidates", save_path+"Tracks", "terr")
plot_ratio(hists, labels, 2, ratio_labels, colors, styles, "Time error [ps]", "Normalized PF candidates", "PU0/PU200", save_path+"Tracks", "terr")
hists = [h3, h6, h9, h12]
plot(hists, labels, colors, styles, "Time significance", "Normalized PF candidates", save_path+"Tracks", "tsig")
plot_ratio(hists, labels, 2, ratio_labels, colors, styles, "Time significance", "Normalized PF candidates", "PU0/PU200", save_path+"Tracks", "tsig")
end = time.time()
print("PF candidate plotting time = %.2fs" % (end-start))

# Plot GenVertex variables
for name in gv_names:
  for var in gv_vars:
    histName = name + "_" + var
    n[histName] = dir_TTnoPU.Get(histName).Clone()
    p[histName] = dir_TTPU200.Get(histName).Clone()
    p[histName].Scale(ratio)
    # p[histName].SetEntries(p[histName].GetEntries() * ratio)
start = time.time()
for iVar, var in enumerate(gv_vars):
  suffix = "_" + var
  hists = [n["gv"+suffix], n["gvB"+suffix], n["gvD"+suffix], p["gv"+suffix], p["gvB"+suffix], p["gvD"+suffix]]
  labels = [ttnoPU+", B & D mother", ttnoPU+", B mother", ttnoPU+", D mother", ttPU200+", B & D mother", ttPU200+", B mother", ttPU200+", D mother"]
  ratio_labels = ["B & D mother", "B mother", "D mother"]
  colors = [COLORS[0], COLORS[1], COLORS[2], COLORS[0], COLORS[1], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed, ROOT.kDashed]
  ylabel = "Gen Vertices"
  if "trk_" in var: ylabel = "Gen Vertex Daughters"
  plot(hists, labels, colors, styles, gv_var_labels[iVar], ylabel, save_path+"GenVertex", var)
  plot_ratio(hists, labels, 1, ratio_labels, colors, styles, gv_var_labels[iVar], ylabel, "PU0/PU200", save_path+"GenVertex", var)
end = time.time()
print("GenVertex plotting time = %.2fs" % (end-start))

# Plot SecondaryVertex variables
for name in sv_names:
  for var in sv_vars:
    histName = name + "_" + var
    getName = histName
    if "tval" in var or "terr" in var or "tsig" in var:
      getName = name + "MTDPV_" + var
    n[histName] = dir_TTnoPU.Get(getName).Clone()
    p[histName] = dir_TTPU200.Get(getName).Clone()
    p[histName].Scale(ratio)
    # p[histName].SetEntries(p[histName].GetEntries() * ratio)
start = time.time()
for iVar, var in enumerate(sv_vars):
  suffix = "_" + var
  hists = [n["svInclusive"+suffix], n["svSlimmed"+suffix], p["svInclusive"+suffix], p["svSlimmed"+suffix]]
  labels = [ttnoPU+", inclusive", ttnoPU+", slimmed", ttPU200+", inclusive", ttPU200+", slimmed"]
  ratio_labels = ["inclusive", "slimmed"]
  colors = [COLORS[0], COLORS[2], COLORS[0], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  ylabel = "Secondary Vertices"
  if "trk_" in var: ylabel = "Secondary Vertex Tracks"
  if "tval" in var or "terr" in var or "tsig" in var:
    hists = [n["svInclusive"+suffix], n["svMerged"+suffix], p["svInclusive"+suffix], p["svMerged"+suffix]]
  plot(hists, labels, colors, styles, sv_var_labels[iVar], ylabel, save_path+"SecondaryVertex", var)
  plot_ratio(hists, labels, 2, ratio_labels, colors, styles, sv_var_labels[iVar], ylabel, "PU0/PU200", save_path+"SecondaryVertex", var)
for name in sv_2d_names:
  for iVar, var in enumerate(sv_var2d):
    hist2DnoPU = dir_TTnoPU.Get(name+"_"+var).Clone()
    hist2DPU200 = dir_TTPU200.Get(name+"_"+var).Clone()
    hist2DPU200.Scale(ratio)
    plot_2D(hist2DnoPU, sv_var2d_labels[iVar][0], sv_var2d_labels[iVar][1], "Secondary Vertex Tracks", save_path+"SecondaryVertex", name+"_"+var+"_noPU")
    plot_2D(hist2DPU200, sv_var2d_labels[iVar][0], sv_var2d_labels[iVar][1], "Secondary Vertex Tracks", save_path+"SecondaryVertex", name+"_"+var+"_PU200")
end = time.time()
print("SecondaryVertex plotting time = %.2fs" % (end-start))

# Plot GenVertex and SecondaryVertex match quality
start = time.time()
for iVar, var in enumerate(vtx_vars):
  h1 = dir_TTnoPU.Get("gv_svInclusive_"+var).Clone()
  h2 = dir_TTnoPU.Get("gv_svSlimmed_"+var).Clone()
  h3 = dir_TTPU200.Get("gv_svInclusive_"+var).Clone()
  h4 = dir_TTPU200.Get("gv_svSlimmed_"+var).Clone()
  h3.Scale(ratio)
  h4.Scale(ratio)
  hists = [h1, h2, h3, h4]
  labels = [ttnoPU+", inclusive", ttnoPU+", slimmed", ttPU200+", inclusive", ttPU200+", slimmed"]
  colors = [COLORS[0], COLORS[2], COLORS[0], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  plot(hists, labels, colors, styles, vtx_var_labels[iVar], "Gen Vertices", save_path+"GenVertex", var)

  h1 = dir_TTnoPU.Get("gv_svSlimmed_"+var).Clone()
  h2 = dir_TTnoPU.Get("gvB_svSlimmed_"+var).Clone()
  h3 = dir_TTnoPU.Get("gvD_svSlimmed_"+var).Clone()
  h4 = dir_TTPU200.Get("gv_svSlimmed_"+var).Clone()
  h5 = dir_TTPU200.Get("gvB_svSlimmed_"+var).Clone()
  h6 = dir_TTPU200.Get("gvD_svSlimmed_"+var).Clone()
  h4.Scale(ratio)
  h5.Scale(ratio)
  h6.Scale(ratio)
  hists = [h1, h2, h3, h4, h5, h6]
  labels = [ttnoPU+", B & D mother", ttnoPU+", B mother", ttnoPU+", D mother", ttPU200+", B & D mother", ttPU200+", B mother", ttPU200+", D mother"]
  colors = [COLORS[0], COLORS[1], COLORS[2], COLORS[0], COLORS[1], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed, ROOT.kDashed]
  plot(hists, labels, colors, styles, vtx_var_labels[iVar], "Gen Vertices", save_path+"GenVertex", var+"_slimmed")
end = time.time()
print("Match quality plotting time = %.2fs" % (end-start))

# Plot GenVertex efficiencies
start = time.time()
for iVar, var in enumerate(efficiency_vars):
  gv_svInclusive_eff_noPU = dir_TTnoPU.Get("gv_svInclusive_"+var).Clone()
  gv_svInclusive_eff_noPU.Divide(dir_TTnoPU.Get("gv_"+var).Clone())
  gv_svSlimmed_eff_noPU = dir_TTnoPU.Get("gv_svSlimmed_"+var).Clone()
  gv_svSlimmed_eff_noPU.Divide(dir_TTnoPU.Get("gv_"+var).Clone())
  gv_svInclusive_eff_PU200 = dir_TTPU200.Get("gv_svInclusive_"+var).Clone()
  gv_svInclusive_eff_PU200.Divide(dir_TTPU200.Get("gv_"+var).Clone())
  gv_svSlimmed_eff_PU200 = dir_TTPU200.Get("gv_svSlimmed_"+var).Clone()
  gv_svSlimmed_eff_PU200.Divide(dir_TTPU200.Get("gv_"+var).Clone())
  hists = [gv_svInclusive_eff_noPU, gv_svSlimmed_eff_noPU, gv_svInclusive_eff_PU200, gv_svSlimmed_eff_PU200]
  labels = [ttnoPU+", inclusive", ttnoPU+", slimmed", ttPU200+", inclusive", ttPU200+", slimmed"]
  colors = [COLORS[0], COLORS[2], COLORS[0], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  plot(hists, labels, colors, styles, eff_var_labels[iVar], "Efficiency", save_path+"Efficiency", var)

  gvB_svSlimmed_eff_noPU = dir_TTnoPU.Get("gvB_svSlimmed_"+var).Clone()
  gvB_svSlimmed_eff_noPU.Divide(dir_TTnoPU.Get("gvB_"+var).Clone())
  gvD_svSlimmed_eff_noPU = dir_TTnoPU.Get("gvD_svSlimmed_"+var).Clone()
  gvD_svSlimmed_eff_noPU.Divide(dir_TTnoPU.Get("gvD_"+var).Clone())
  gvB_svSlimmed_eff_PU200 = dir_TTPU200.Get("gvB_svSlimmed_"+var).Clone()
  gvB_svSlimmed_eff_PU200.Divide(dir_TTPU200.Get("gvB_"+var).Clone())
  gvD_svSlimmed_eff_PU200 = dir_TTPU200.Get("gvD_svSlimmed_"+var).Clone()
  gvD_svSlimmed_eff_PU200.Divide(dir_TTPU200.Get("gvD_"+var).Clone())
  hists = [gv_svSlimmed_eff_noPU, gvB_svSlimmed_eff_noPU, gvD_svSlimmed_eff_noPU, gv_svSlimmed_eff_PU200, gvB_svSlimmed_eff_PU200, gvD_svSlimmed_eff_PU200]
  labels = [ttnoPU+", B & D mother", ttnoPU+", B mother", ttnoPU+", D mother", ttPU200+", B & D mother", ttPU200+", B mother", ttPU200+", D mother"]
  colors = [COLORS[0], COLORS[1], COLORS[2], COLORS[0], COLORS[1], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed, ROOT.kDashed]
  plot(hists, labels, colors, styles, eff_var_labels[iVar], "Efficiency", save_path+"Efficiency", var+"_slimmed")
end = time.time()
print("Efficiency plotting time = %.2fs" % (end-start))

endScript = time.time()
print("Total plotting time = %.2fs" % (endScript-startScript))
