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

# Maximum matching efficiencies
histMaxMatchedWithInclusiveNoPU = dir_TTnoPU.Get("maxMatchedWithInclusive")
histMaxMatchedWithSlimmedNoPU = dir_TTnoPU.Get("maxMatchedWithSlimmed")
histMaxMatchedWithInclusivePU200 = dir_TTPU200.Get("maxMatchedWithInclusive")
histMaxMatchedWithSlimmedPU200 = dir_TTPU200.Get("maxMatchedWithSlimmed")
histnGVNoPU = dir_TTnoPU.Get("ngv")
histnGVPU200 = dir_TTPU200.Get("ngv")

maxMatchedWithInclusiveNoPU = 0.0
maxMatchedWithSlimmedNoPU = 0.0
maxMatchedWithInclusivePU200 = 0.0
maxMatchedWithSlimmedPU200 = 0.0
nGVNoPU = 0.0
nGVPU200 = 0.0
for i in range(1, histMaxMatchedWithInclusiveNoPU.GetNbinsX()):
  maxMatchedWithInclusiveNoPU += histMaxMatchedWithInclusiveNoPU.GetBinContent(i)*(i-1)
for i in range(1, histMaxMatchedWithSlimmedNoPU.GetNbinsX()):
  maxMatchedWithSlimmedNoPU += histMaxMatchedWithSlimmedNoPU.GetBinContent(i)*(i-1)
for i in range(1, histMaxMatchedWithInclusivePU200.GetNbinsX()):
  maxMatchedWithInclusivePU200 += histMaxMatchedWithInclusivePU200.GetBinContent(i)*(i-1)
for i in range(1, histMaxMatchedWithSlimmedPU200.GetNbinsX()):
  maxMatchedWithSlimmedPU200 += histMaxMatchedWithSlimmedPU200.GetBinContent(i)*(i-1)
for i in range(1, histnGVNoPU.GetNbinsX()):
  nGVNoPU += histnGVNoPU.GetBinContent(i)*(i-1)
for i in range(1, histnGVPU200.GetNbinsX()):
  nGVPU200 += histnGVPU200.GetBinContent(i)*(i-1)
maxMatchedWithInclusiveNoPU /= nGVNoPU
maxMatchedWithSlimmedNoPU /= nGVNoPU
maxMatchedWithInclusivePU200 /= nGVPU200
maxMatchedWithSlimmedPU200 /= nGVPU200

print("Maximum possible matching efficiency with inclusive collection with no PU  = %.2f percent" % (maxMatchedWithInclusiveNoPU*100.0))
print("Maximum possible matching efficiency with slimmed collection with no PU    = %.2f percent" % (maxMatchedWithSlimmedNoPU*100.0))
print("Maximum possible matching efficiency with inclusive collection with PU 200 = %.2f percent" % (maxMatchedWithInclusivePU200*100.0))
print("Maximum possible matching efficiency with slimmed collection with PU 200   = %.2f percent" % (maxMatchedWithSlimmedPU200*100.0))

# Plot PF candidate times
start = time.time()
for mtdRegion in ["BTL", "ETL"]:
  h1 = dir_TTnoPU.Get("pfCandidate_bs_tval"+mtdRegion).Clone()
  h2 = dir_TTnoPU.Get("pfCandidate_bs_terr"+mtdRegion).Clone()
  h3 = dir_TTnoPU.Get("pfCandidate_bs_tsig"+mtdRegion).Clone()
  h4 = dir_TTnoPU.Get("pfCandidate_pv_tval"+mtdRegion).Clone()
  h5 = dir_TTnoPU.Get("pfCandidate_pv_terr"+mtdRegion).Clone()
  h6 = dir_TTnoPU.Get("pfCandidate_pv_tsig"+mtdRegion).Clone()
  h7 = dir_TTPU200.Get("pfCandidate_bs_tval"+mtdRegion).Clone()
  h8 = dir_TTPU200.Get("pfCandidate_bs_terr"+mtdRegion).Clone()
  h9 = dir_TTPU200.Get("pfCandidate_bs_tsig"+mtdRegion).Clone()
  h10 = dir_TTPU200.Get("pfCandidate_pv_tval"+mtdRegion).Clone()
  h11 = dir_TTPU200.Get("pfCandidate_pv_terr"+mtdRegion).Clone()
  h12 = dir_TTPU200.Get("pfCandidate_pv_tsig"+mtdRegion).Clone()
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
  colors = [CMS.p8.kOrange, CMS.p8.kBlue, CMS.p8.kOrange, CMS.p8.kBlue]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  markers = [26, 26, 32, 32]
  hists = [h1, h4, h7, h10]
  plot(hists, labels, colors, markers, styles, "Time [ps]", "Normalized PF candidates", save_path+"Tracks", "tval"+mtdRegion)
  plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, "Time [ps]", "Normalized PF candidates", "PU0/PU200", save_path+"Tracks", "tval"+mtdRegion)
  hists = [h2, h5, h8, h11]
  plot(hists, labels, colors, markers, styles, "Time error [ps]", "Normalized PF candidates", save_path+"Tracks", "terr"+mtdRegion)
  plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, "Time error [ps]", "Normalized PF candidates", "PU0/PU200", save_path+"Tracks", "terr"+mtdRegion)
  hists = [h3, h6, h9, h12]
  plot(hists, labels, colors, markers, styles, "Time significance", "Normalized PF candidates", save_path+"Tracks", "tsig"+mtdRegion)
  plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, "Time significance", "Normalized PF candidates", "PU0/PU200", save_path+"Tracks", "tsig"+mtdRegion)
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
  hists = [n["gvB"+suffix], n["gvD"+suffix], p["gvB"+suffix], p["gvD"+suffix]]
  labels = [ttnoPU+", B mother", ttnoPU+", D mother", ttPU200+", B mother", ttPU200+", D mother"]
  ratio_labels = ["B mother", "D mother"]
  colors = [COLORS[4], COLORS[3], COLORS[4], COLORS[3]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  markers = [26, 26, 32, 32]
  ylabel = "Gen Vertices"
  if "trk_" in var: ylabel = "Gen Vertex Daughters"
  plot(hists, labels, colors, markers, styles, gv_var_labels[iVar], ylabel, save_path+"GenVertex", var)
  plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, gv_var_labels[iVar], ylabel, "PU0/PU200", save_path+"GenVertex", var)
end = time.time()
print("GenVertex plotting time = %.2fs" % (end-start))

# Plot SecondaryVertex variables
for name in sv_names:
  for var in sv_vars:
    histName = name + "_" + var
    n[histName] = dir_TTnoPU.Get(histName).Clone()
    p[histName] = dir_TTPU200.Get(histName).Clone()
    if n[histName].Integral() != 0:
      n[histName].Scale(1.0 / n[histName].Integral())
    if p[histName].Integral() != 0:
      p[histName].Scale(1.0 / p[histName].Integral())
    # p[histName].Scale(ratio)
    # p[histName].SetEntries(p[histName].GetEntries() * ratio)
start = time.time()
for refPoint in ["BS", "PV"]:
  for iVar, var in enumerate(sv_vars):
    suffix = "MTD" + refPoint + "_" + var
    hists = [n["svInclusive"+suffix], n["svSlimmed"+suffix], p["svInclusive"+suffix], p["svSlimmed"+suffix]]
    labels = [ttnoPU+", inclusive", ttnoPU+", slimmed", ttPU200+", inclusive", ttPU200+", slimmed"]
    ratio_labels = ["inclusive", "slimmed"]
    colors = [COLORS[0], COLORS[2], COLORS[0], COLORS[2]]
    styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
    markers = [26, 26, 32, 32]
    ylabel = "Normalized Secondary Vertices"
    if "trk_" in var: ylabel = "Normalized Secondary Vertex Tracks"
    if "tval" in var or "terr" in var or "tsig" in var:
      hists = [n["svInclusive"+suffix], n["svSlimmed"+suffix], p["svInclusive"+suffix], p["svSlimmed"+suffix]]
    plot(hists, labels, colors, markers, styles, sv_var_labels[iVar], ylabel, save_path+"SecondaryVertex", refPoint+"_"+var)
    plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, sv_var_labels[iVar], ylabel, "PU0/PU200", save_path+"SecondaryVertex", refPoint+"_"+var)
  for name in sv_names:
    for iVar, var in enumerate(sv_var2d):
      print(name+"_"+var)
      hist2DnoPU = dir_TTnoPU.Get(name+"_"+var).Clone()
      hist2DPU200 = dir_TTPU200.Get(name+"_"+var).Clone()
      hist2DPU200.Scale(ratio)
      plot_2D(hist2DnoPU, sv_var2d_labels[iVar][0], sv_var2d_labels[iVar][1], "Secondary Vertex Tracks", save_path+"SecondaryVertex", refPoint+"_"+name+"_"+var+"_noPU")
      plot_2D(hist2DPU200, sv_var2d_labels[iVar][0], sv_var2d_labels[iVar][1], "Secondary Vertex Tracks", save_path+"SecondaryVertex", refPoint+"_"+name+"_"+var+"_PU200")
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
  ratio_labels = ["inclusive", "slimmed"]
  colors = [COLORS[0], COLORS[2], COLORS[0], COLORS[2]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  markers = [26, 26, 32, 32]
  plot(hists, labels, colors, markers, styles, vtx_var_labels[iVar], "Gen Vertices", save_path+"GenVertex", var)
  plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, vtx_var_labels[iVar], "Gen Vertices", "PU0/PU200", save_path+"GenVertex", var)

  h1 = dir_TTnoPU.Get("gvB_svSlimmed_"+var).Clone()
  h2 = dir_TTnoPU.Get("gvD_svSlimmed_"+var).Clone()
  h3 = dir_TTPU200.Get("gvB_svSlimmed_"+var).Clone()
  h4 = dir_TTPU200.Get("gvD_svSlimmed_"+var).Clone()
  h3.Scale(ratio)
  h4.Scale(ratio)
  hists = [h1, h2, h3, h4]
  labels = [ttnoPU+", B mother", ttnoPU+", D mother", ttPU200+", B mother", ttPU200+", D mother"]
  colors = [COLORS[4], COLORS[3], COLORS[4], COLORS[3]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  markers = [26, 26, 32, 32]
  plot(hists, labels, colors, markers, styles, vtx_var_labels[iVar], "Gen Vertices", save_path+"GenVertex", var+"_slimmed")
  plot_ratio(hists, labels, 2, ratio_labels, colors, markers, styles, vtx_var_labels[iVar], "Gen Vertices", "PU0/PU200", save_path+"GenVertex", var+"_slimmed")
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
  markers = [26, 26, 32, 32]
  plot(hists, labels, colors, markers, styles, eff_var_labels[iVar], "Efficiency", save_path+"Efficiency", var)

  gvB_svSlimmed_eff_noPU = dir_TTnoPU.Get("gvB_svSlimmed_"+var).Clone()
  gvB_svSlimmed_eff_noPU.Divide(dir_TTnoPU.Get("gvB_"+var).Clone())
  gvD_svSlimmed_eff_noPU = dir_TTnoPU.Get("gvD_svSlimmed_"+var).Clone()
  gvD_svSlimmed_eff_noPU.Divide(dir_TTnoPU.Get("gvD_"+var).Clone())
  gvB_svSlimmed_eff_PU200 = dir_TTPU200.Get("gvB_svSlimmed_"+var).Clone()
  gvB_svSlimmed_eff_PU200.Divide(dir_TTPU200.Get("gvB_"+var).Clone())
  gvD_svSlimmed_eff_PU200 = dir_TTPU200.Get("gvD_svSlimmed_"+var).Clone()
  gvD_svSlimmed_eff_PU200.Divide(dir_TTPU200.Get("gvD_"+var).Clone())
  hists = [gvB_svSlimmed_eff_noPU, gvD_svSlimmed_eff_noPU, gvB_svSlimmed_eff_PU200, gvD_svSlimmed_eff_PU200]
  labels = [ttnoPU+", B mother", ttnoPU+", D mother", ttPU200+", B mother", ttPU200+", D mother"]
  colors = [COLORS[4], COLORS[3], COLORS[4], COLORS[3]]
  styles = [ROOT.kSolid, ROOT.kSolid, ROOT.kDashed, ROOT.kDashed]
  markers = [26, 26, 32, 32]
  plot(hists, labels, colors, markers, styles, eff_var_labels[iVar], "Efficiency", save_path+"Efficiency", var+"_slimmed")
end = time.time()
print("Efficiency plotting time = %.2fs" % (end-start))

endScript = time.time()
print("Total plotting time = %.2fs" % (endScript-startScript))
