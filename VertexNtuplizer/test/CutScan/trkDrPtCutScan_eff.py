import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")

samples = ["$t\\bar{t}$ no PU", "$t\\bar{t}$ PU 200"]
sampleNames = ["TTToHadronic_noPU", "TTToHadronic_PU200"]
counts = [2000, 299]

dRfix = 4.0
pTfix = 1.0

sort = "noSort"

savePath = "/eos/home-b/btvweb/www/Offline/Phase2/etsai/Phase2/CutScan/"
labels = ["GV-SV", "GV-SVt", "GVs-SV", "GVs-SVt", "GVn-SV", "GVn-SVt", "GVns-SV", "GVns-SVt"]

for iSample, sample in enumerate(samples):
  sampleName = sampleNames[iSample]
  nevents = counts[iSample]

  print("Plotting " + sampleName)

  fileName = sampleName + "_cutScan_" + sort + ".txt"
  saveName = savePath + sampleName + "_" + sort

  dR = []
  dReff = [[], [], [], [], [], [], [], []]
  dRdxy = [[], [], [], [], [], [], [], []]
  dRdxyerr = [[], [], [], [], [], [], [], []]
  dRd3d = [[], [], [], [], [], [], [], []]
  dRd3derr = [[], [], [], [], [], [], [], []]
  pT = []
  pTeff = [[], [], [], [], [], [], [], []]
  pTdxy = [[], [], [], [], [], [], [], []]
  pTdxyerr = [[], [], [], [], [], [], [], []]
  pTd3d = [[], [], [], [], [], [], [], []]
  pTd3derr = [[], [], [], [], [], [], [], []]

  with open(fileName) as f:
    line = f.readline().strip()
    while len(line) > 0:
      vals = line.split()
      if len(vals) == 3:
        if vals[0] == "trkMatchDrCut":
          dR.append(float(vals[2]))
        elif vals[0] == "trkMatchPtCut":
          pT.append(float(vals[2]))
      elif len(vals) == 6:
        if vals[0] == "gv_sv:":
          dReff[0].append(float(vals[1]))
          dRdxy[0].append(float(vals[2]))
          dRdxyerr[0].append(float(vals[3]))
          dRd3d[0].append(float(vals[4]))
          dRd3derr[0].append(float(vals[5]))
          pTeff[0].append(float(vals[1]))
          pTdxy[0].append(float(vals[2]))
          pTdxyerr[0].append(float(vals[3]))
          pTd3d[0].append(float(vals[4]))
          pTd3derr[0].append(float(vals[5]))
        elif vals[0] == "gv_svt:":
          dReff[1].append(float(vals[1]))
          dRdxy[1].append(float(vals[2]))
          dRdxyerr[1].append(float(vals[3]))
          dRd3d[1].append(float(vals[4]))
          dRd3derr[1].append(float(vals[5]))
          pTeff[1].append(float(vals[1]))
          pTdxy[1].append(float(vals[2]))
          pTdxyerr[1].append(float(vals[3]))
          pTd3d[1].append(float(vals[4]))
          pTd3derr[1].append(float(vals[5]))
        elif vals[0] == "gvs_sv:":
          dReff[2].append(float(vals[1]))
          dRdxy[2].append(float(vals[2]))
          dRdxyerr[2].append(float(vals[3]))
          dRd3d[2].append(float(vals[4]))
          dRd3derr[2].append(float(vals[5]))
          pTeff[2].append(float(vals[1]))
          pTdxy[2].append(float(vals[2]))
          pTdxyerr[2].append(float(vals[3]))
          pTd3d[2].append(float(vals[4]))
          pTd3derr[2].append(float(vals[5]))
        elif vals[0] == "gvs_svt:":
          dReff[3].append(float(vals[1]))
          dRdxy[3].append(float(vals[2]))
          dRdxyerr[3].append(float(vals[3]))
          dRd3d[3].append(float(vals[4]))
          dRd3derr[3].append(float(vals[5]))
          pTeff[3].append(float(vals[1]))
          pTdxy[3].append(float(vals[2]))
          pTdxyerr[3].append(float(vals[3]))
          pTd3d[3].append(float(vals[4]))
          pTd3derr[3].append(float(vals[5]))
        elif vals[0] == "gvn_sv:":
          dReff[4].append(float(vals[1]))
          dRdxy[4].append(float(vals[2]))
          dRdxyerr[4].append(float(vals[3]))
          dRd3d[4].append(float(vals[4]))
          dRd3derr[4].append(float(vals[5]))
          pTeff[4].append(float(vals[1]))
          pTdxy[4].append(float(vals[2]))
          pTdxyerr[4].append(float(vals[3]))
          pTd3d[4].append(float(vals[4]))
          pTd3derr[4].append(float(vals[5]))
        elif vals[0] == "gvn_svt:":
          dReff[5].append(float(vals[1]))
          dRdxy[5].append(float(vals[2]))
          dRdxyerr[5].append(float(vals[3]))
          dRd3d[5].append(float(vals[4]))
          dRd3derr[5].append(float(vals[5]))
          pTeff[5].append(float(vals[1]))
          pTdxy[5].append(float(vals[2]))
          pTdxyerr[5].append(float(vals[3]))
          pTd3d[5].append(float(vals[4]))
          pTd3derr[5].append(float(vals[5]))
        elif vals[0] == "gvns_sv:":
          dReff[6].append(float(vals[1]))
          dRdxy[6].append(float(vals[2]))
          dRdxyerr[6].append(float(vals[3]))
          dRd3d[6].append(float(vals[4]))
          dRd3derr[6].append(float(vals[5]))
          pTeff[6].append(float(vals[1]))
          pTdxy[6].append(float(vals[2]))
          pTdxyerr[6].append(float(vals[3]))
          pTd3d[6].append(float(vals[4]))
          pTd3derr[6].append(float(vals[5]))
        elif vals[0] == "gvns_svt:":
          dReff[7].append(float(vals[1]))
          dRdxy[7].append(float(vals[2]))
          dRdxyerr[7].append(float(vals[3]))
          dRd3d[7].append(float(vals[4]))
          dRd3derr[7].append(float(vals[5]))
          pTeff[7].append(float(vals[1]))
          pTdxy[7].append(float(vals[2]))
          pTdxyerr[7].append(float(vals[3]))
          pTd3d[7].append(float(vals[4]))
          pTd3derr[7].append(float(vals[5]))
      line = f.readline().strip()

  dR = dR[0:24]
  for i in range(len(dReff)):
    dReff[i] = dReff[i][0:24]
    dRdxy[i] = dRdxy[i][0:24]
    dRdxyerr[i] = dRdxyerr[i][0:24]
    dRd3d[i] = dRd3d[i][0:24]
    dRd3derr[i] = dRd3derr[i][0:24]

  pT = pT[24:39]
  for i in range(len(pTeff)):
    pTeff[i] = pTeff[i][24:39]
    pTdxy[i] = pTdxy[i][24:39]
    pTdxyerr[i] = pTdxyerr[i][24:39]
    pTd3d[i] = pTd3d[i][24:39]
    pTd3derr[i] = pTd3derr[i][24:39]

  # print(dR)
  # print(dReff[0])
  # print(dReff[1])
  # print(dReff[2])
  # print(dReff[3])
  # print(dReff[4])
  # print(dReff[5])
  # print(dReff[6])
  # print(dReff[7])

  # print(pT)
  # print(pTeff[0])
  # print(pTeff[1])
  # print(pTeff[2])
  # print(pTeff[3])
  # print(pTeff[4])
  # print(pTeff[5])
  # print(pTeff[6])
  # print(pTeff[7])

  plt.figure()
  for i in range(8):
    plt.plot(dR, dReff[i], marker=".", linestyle="-", markersize=15, label=labels[i])
  plt.legend()
  plt.xlim(-0.1, dRfix + 0.1)
  plt.ylim(0.0, 1.2)
  plt.xlabel("$\Delta R$ cut with $p_{T}$ ratio cut = %.1f" % pTfix)
  plt.ylabel("Matching efficiency")
  plt.grid(True, which="both")
  hep.style.use("CMS")
  hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
  plt.savefig(saveName + "_dReff.png", bbox_inches="tight")
  plt.savefig(saveName + "_dReff.pdf", bbox_inches="tight")
  plt.close()

  plt.figure()
  for i in range(8):
    plt.errorbar(dR, dRdxy[i], yerr=dRdxyerr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
  plt.legend()
  plt.xlim(-0.1, dRfix + 0.1)
  plt.ylim(0.0, 2.0)
  plt.xlabel("$\Delta R$ cut with $p_{T}$ ratio cut = %.1f" % pTfix)
  plt.ylabel("Match <$d_{xy}$> [cm]")
  plt.grid(True, which="both")
  hep.style.use("CMS")
  hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
  plt.savefig(saveName + "_dRdxy.png", bbox_inches="tight")
  plt.savefig(saveName + "_dRdxy.pdf", bbox_inches="tight")
  plt.close()

  plt.figure()
  for i in range(8):
    plt.errorbar(dR, dRd3d[i], yerr=dRd3derr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
  plt.legend()
  plt.xlim(-0.1, dRfix + 0.1)
  plt.ylim(0.0, 2.0)
  plt.xlabel("$\Delta R$ cut with $p_{T}$ ratio cut = %.1f" % pTfix)
  plt.ylabel("Match <$d_{3D}$> [cm]")
  plt.grid(True, which="both")
  hep.style.use("CMS")
  hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
  plt.savefig(saveName + "_dRd3d.png", bbox_inches="tight")
  plt.savefig(saveName + "_dRd3d.pdf", bbox_inches="tight")
  plt.close()

  plt.figure()
  for i in range(8):
    plt.plot(pT, pTeff[i], marker=".", linestyle="-", markersize=15, label=labels[i])
  plt.legend()
  plt.xlim(-0.1, pTfix + 0.1)
  plt.ylim(0.0, 1.2)
  plt.xlabel("$p_{T}$ ratio cut with $\Delta R$ cut = %.1f" % dRfix)
  plt.ylabel("Matching efficiency")
  plt.grid(True, which="both")
  hep.style.use("CMS")
  hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
  plt.savefig(saveName + "_pTeff.png", bbox_inches="tight")
  plt.savefig(saveName + "_pTeff.pdf", bbox_inches="tight")
  plt.close()

  plt.figure()
  for i in range(8):
    plt.errorbar(pT, pTdxy[i], yerr=pTdxyerr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
  plt.legend()
  plt.xlim(-0.1, pTfix + 0.1)
  plt.ylim(0.0, 2.0)
  plt.xlabel("$p_{T}$ ratio cut with $\Delta R$ cut = %.1f" % dRfix)
  plt.ylabel("Match <$d_{xy}$> [cm]")
  plt.grid(True, which="both")
  hep.style.use("CMS")
  hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
  plt.savefig(saveName + "_pTdxy.png", bbox_inches="tight")
  plt.savefig(saveName + "_pTdxy.pdf", bbox_inches="tight")
  plt.close()

  plt.figure()
  for i in range(8):
    plt.errorbar(pT, pTd3d[i], yerr=pTd3derr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
  plt.legend()
  plt.xlim(-0.1, pTfix + 0.1)
  plt.ylim(0.0, 2.0)
  plt.xlabel("$p_{T}$ ratio cut with $\Delta R$ cut = %.1f" % dRfix)
  plt.ylabel("Match <$d_{3D}$> [cm]")
  plt.grid(True, which="both")
  hep.style.use("CMS")
  hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
  plt.savefig(saveName + "_pTd3d.png", bbox_inches="tight")
  plt.savefig(saveName + "_pTd3d.pdf", bbox_inches="tight")
  plt.close()
