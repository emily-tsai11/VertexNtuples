import ROOT
ROOT.gROOT.SetBatch(1)
import sys
sys.path.insert(0, "../cmsstyle/src/cmsstyle")
import cmsstyle as CMS


style1 = "hist"
style2 = "p e x0"


def plot(hists, labels, colors, markers, styles, xlabel, ylabel, save_path, save_name, header="", log=False):

    CMS.SetLumi("")
    CMS.SetEnergy("14")
    CMS.SetExtraText("Work in Progress")

    left = hists[0].GetBinLowEdge(1)
    right = hists[0].GetBinLowEdge(hists[0].GetNbinsX() + 1)
    maxi = -1.0
    for hist in hists:
        temp = hist.Clone()
        # temp.Scale(1.0/hist.Integral())
        maxi = max(maxi, temp.GetMaximum())

    canv = CMS.cmsCanvas(save_name,
        left,
        right,
        0.01 if log else 0,
        50 * maxi if log else 1.6 * maxi,
        xlabel, ylabel, square=CMS.kRectangular, extraSpace=0.03, iPos=0)
    canv.SetLogy(log)
    leg = CMS.cmsLeg(0.3, 0.89 - 0.036 * 7, 0.89, 0.89, textSize=0.034)
    if len(header): CMS.cmsHeader(leg, header, textSize=0.034)

    final_hists = []
    for iHist, hist in enumerate(hists):
        h = hist.Clone()
        # h.Scale(1.0/h.Integral())
        CMS.cmsDraw(h, style1, mcolor=colors[iHist], fstyle=0, lwidth=3, lstyle=styles[iHist], marker=markers[iHist])
        final_hists.append(h)
        points = h.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[iHist], fstyle=0, lwidth=3, marker=markers[iHist])
        final_hists.append(points)
        hlabel = "%s, #mu=%.2f, #sigma=%.2f" % (labels[iHist], h.GetMean(), h.GetStdDev())
        if "pdgIdBin" in save_name:
            hlabel = labels[iHist]
        leg.AddEntry(h, hlabel, "PLE")

    if "pdgIdBin" in save_name:
        CMS.GetcmsCanvasHist(canv).GetXaxis().SetBinLabel(100, "B meson")
        CMS.GetcmsCanvasHist(canv).GetXaxis().SetBinLabel(300, "B baryon")
        CMS.GetcmsCanvasHist(canv).GetXaxis().SetBinLabel(500, "C meson")
        CMS.GetcmsCanvasHist(canv).GetXaxis().SetBinLabel(700, "C baryon")
        CMS.GetcmsCanvasHist(canv).GetXaxis().LabelsOption("h")

    canv.SetRightMargin(0.03)
    CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleOffset(1.2)
    CMS.GetcmsCanvasHist(canv).GetXaxis().SetLabelSize(0.04)
    CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleSize(0.05)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleOffset(1.25)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetLabelSize(0.04)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleSize(0.05)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetMaxDigits(2)
    if maxi < 100.0: CMS.GetcmsCanvasHist(canv).GetYaxis().SetNoExponent()
    if maxi > 100.0: ROOT.TGaxis.SetExponentOffset(-0.065, 0.015, "Y")

    CMS.SaveCanvas(canv, save_path+"/"+save_name+".png", close=False)
    CMS.SaveCanvas(canv, save_path+"/"+save_name+".pdf")


def plot_ratio(hists, labels, ratio_style, ratio_labels, colors, markers, styles, xlabel,
        ylabel, ratiolabel, save_path, save_name, header="", log=False, ratioLog=False):

    CMS.SetLumi("")
    CMS.SetEnergy("14")
    CMS.SetExtraText("Work in Progress")

    left = hists[0].GetBinLowEdge(1)
    right = hists[0].GetBinLowEdge(hists[0].GetNbinsX() + 1)
    maxi = -1.0
    for hist in hists:
        temp = hist.Clone()
        # temp.Scale(1.0/hist.Integral())
        maxi = max(maxi, temp.GetMaximum())

    ratioMin = 0.0
    ratioMax = 2.0
    if "SecondaryVertex" in save_path:
        if "pt" in save_name:
            ratioMax = 4.0
            log = True
        if "trk_terrETL" in save_name:
            ratioMin = 0.8e-1
            ratioMax = 80.0
            ratioLog = True
    if "Tracks" in save_path:
        if "terrETL" in save_name:
            ratioMin = 0.8e-1
            ratioMax = 20.0
            ratioLog = True
    dicanv = CMS.cmsDiCanvas(
        save_name + "_ratio",
        left,
        right,
        0.5 if log else 0,
        2500 * maxi if log else 1.6 * maxi,
        ratioMin, ratioMax,
        xlabel, ylabel,
        ratiolabel,
        square=CMS.kRectangular, extraSpace=0.01, iPos=0)
    topcanv = dicanv.cd(1)
    topcanv.SetLogy(log)
    leg = CMS.cmsLeg(0.3, 0.89 - 0.036 * 7, 0.89, 0.89, textSize=0.034)
    if len(header): CMS.cmsHeader(leg, header, textSize=0.034)

    final_hists = []
    final_points = []
    for iHist, hist in enumerate(hists):
        h = hist.Clone()
        # h.Scale(1.0/h.Integral())
        CMS.cmsDraw(h, style1, mcolor=colors[iHist], fstyle=0, lwidth=3, lstyle=styles[iHist], marker=markers[iHist])
        final_hists.append(h)
        points = h.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[iHist], fstyle=0, lwidth=3, marker=markers[iHist])
        final_points.append(points)
        hlabel = "%s, #mu=%.2f, #sigma=%.2f" % (labels[iHist], h.GetMean(), h.GetStdDev())
        if "pdgIdBin" in save_name:
            hlabel = labels[iHist]
        leg.AddEntry(h, hlabel, "PLE")

    CMS.GetcmsCanvasHist(topcanv).GetYaxis().SetTitleOffset(1.1)
    CMS.GetcmsCanvasHist(topcanv).GetYaxis().SetLabelSize(0.05)
    CMS.GetcmsCanvasHist(topcanv).GetYaxis().SetTitleSize(0.06)
    CMS.GetcmsCanvasHist(topcanv).GetYaxis().SetMaxDigits(2)
    if maxi < 100.0 and not log: CMS.GetcmsCanvasHist(topcanv).GetYaxis().SetNoExponent()
    if maxi > 100.0: ROOT.TGaxis.SetExponentOffset(-0.055, 0.015, "Y")

    botcanv = dicanv.cd(2)
    botcanv.SetLogy(ratioLog)
    final_ratio = []
    if ratio_style == 0:
        for iHist in range(len(final_hists)):
            if iHist == 0: continue # skip nominal case
            ratio = final_hists[iHist].Clone()
            ratio.Divide(final_hists[0])
            CMS.cmsDraw(ratio, style1, mcolor=colors[iHist], fstyle=0, lwidth=3, marker=markers[iHist])
            final_ratio.append(ratio)
    elif ratio_style == 1: # exactly 6 histograms
        legRatio = CMS.cmsLeg(0.89 - 0.25 * 3, 0.89 - 0.11, 0.89, 0.89, textSize=0.07, columns=3)
        # -------- #
        ratio = final_hists[0].Clone()
        ratio.Divide(final_hists[3])
        CMS.cmsDraw(ratio, style1, mcolor=colors[0], fstyle=0, lwidth=3, lstyle=styles[0], marker=markers[0])
        points = ratio.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[0], fstyle=0, lwidth=3, marker=markers[0])
        final_ratio.append(ratio)
        final_ratio.append(points)
        legRatio.AddEntry(ratio, ratio_labels[0], "PLE")
        # -------- #
        ratio = final_hists[1].Clone()
        ratio.Divide(final_hists[4])
        CMS.cmsDraw(ratio, style1, mcolor=colors[1], fstyle=0, lwidth=3, lstyle=styles[1], marker=markers[1])
        points = ratio.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[1], fstyle=0, lwidth=3, marker=markers[1])
        final_ratio.append(ratio)
        final_ratio.append(points)
        legRatio.AddEntry(ratio, ratio_labels[1], "PLE")
        # -------- #
        ratio = final_hists[2].Clone()
        ratio.Divide(final_hists[5])
        CMS.cmsDraw(ratio, style1, mcolor=colors[2], fstyle=0, lwidth=3, lstyle=styles[2], marker=markers[2])
        points = ratio.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[2], fstyle=0, lwidth=3, marker=markers[2])
        final_ratio.append(ratio)
        final_ratio.append(points)
        legRatio.AddEntry(ratio, ratio_labels[2], "PLE")
        # -------- #
    elif ratio_style == 2: # exactly 4 histograms
        legRatio = CMS.cmsLeg(0.89 - 0.25 * 2, 0.89 - 0.11, 0.89, 0.89, textSize=0.07, columns=3)
        # -------- #
        ratio = final_hists[0].Clone()
        ratio.Divide(final_hists[2])
        CMS.cmsDraw(ratio, style1, mcolor=colors[0], fstyle=0, lwidth=3, lstyle=styles[0], marker=markers[0])
        points = ratio.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[0], fstyle=0, lwidth=3, marker=markers[0])
        final_ratio.append(ratio)
        final_ratio.append(points)
        legRatio.AddEntry(ratio, ratio_labels[0], "PLE")
        # -------- #
        ratio = final_hists[1].Clone()
        ratio.Divide(final_hists[3])
        CMS.cmsDraw(ratio, style1, mcolor=colors[1], fstyle=0, lwidth=3, lstyle=styles[1], marker=markers[1])
        points = ratio.Clone()
        CMS.cmsDraw(points, style2, mcolor=colors[1], fstyle=0, lwidth=3, marker=markers[1])
        final_ratio.append(ratio)
        final_ratio.append(points)
        legRatio.AddEntry(ratio, ratio_labels[1], "PLE")
        # -------- #
    ref_line = ROOT.TLine(left, 1, right, 1)
    CMS.cmsDrawLine(ref_line, lcolor=ROOT.kBlack, lstyle=ROOT.kDotted)

    if "pdgIdBin" in save_name:
        CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetBinLabel(100, "B meson")
        CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetBinLabel(300, "B baryon")
        CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetBinLabel(500, "C meson")
        CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetBinLabel(700, "C baryon")
        CMS.GetcmsCanvasHist(botcanv).GetXaxis().LabelsOption("h")

    CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetTitleOffset(1.2)
    CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetLabelSize(0.1)
    CMS.GetcmsCanvasHist(botcanv).GetXaxis().SetTitleSize(0.11)
    CMS.GetcmsCanvasHist(botcanv).GetYaxis().SetTitleOffset(0.6)
    CMS.GetcmsCanvasHist(botcanv).GetYaxis().SetLabelSize(0.1)
    CMS.GetcmsCanvasHist(botcanv).GetYaxis().SetTitleSize(0.11)
    CMS.GetcmsCanvasHist(botcanv).GetYaxis().SetMaxDigits(2)
    CMS.GetcmsCanvasHist(botcanv).GetYaxis().SetNoExponent()

    CMS.SaveCanvas(dicanv, save_path+"/"+save_name+"_ratio.png", close=False)
    CMS.SaveCanvas(dicanv, save_path+"/"+save_name+"_ratio.pdf")

def plot_2D(hist, xlabel, ylabel, zlabel, save_path, save_name):

    CMS.SetLumi("")
    CMS.SetEnergy("14")
    CMS.SetExtraText("Work in Progress")

    leftX = hist.GetXaxis().GetBinLowEdge(1)
    rightX = hist.GetXaxis().GetBinLowEdge(hist.GetNbinsX() + 1)
    leftY = hist.GetYaxis().GetBinLowEdge(1)
    rightY = hist.GetYaxis().GetBinLowEdge(hist.GetNbinsY() + 1)

    canv = CMS.cmsCanvas(
        save_name,
        leftX, rightX, leftY, rightY,
        xlabel, ylabel,
        square=CMS.kRectangular, extraSpace=0.04, iPos=0,
        with_z_axis=True, scaleLumi=0.8,
    )

    hist.GetZaxis().SetTitle(zlabel)
    hist.Draw("same colz")

    canv.SetRightMargin(0.15)
    canv.SetLogz()
    CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleOffset(1.2)
    CMS.GetcmsCanvasHist(canv).GetXaxis().SetLabelSize(0.04)
    CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleSize(0.05)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleOffset(1.3)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetLabelSize(0.04)
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleSize(0.05)
    hist.GetZaxis().SetTitleOffset(1.1)
    hist.GetZaxis().SetLabelSize(0.04)
    hist.GetZaxis().SetTitleSize(0.05)
    hist.GetZaxis().SetMaxDigits(2)

    CMS.SetAlternative2DColor(hist, CMS.cmsStyle)
    CMS.UpdatePalettePosition(hist, canv)
    # CMS.UpdatePalettePosition(hist, Y2=0.9)

    CMS.SaveCanvas(canv, save_path+"/"+save_name+".png", close=False)
    CMS.SaveCanvas(canv, save_path+"/"+save_name+".pdf")
