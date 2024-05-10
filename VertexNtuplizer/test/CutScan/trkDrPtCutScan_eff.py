import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")

nevents = 2000
sample = "$t\\bar{t}$ no PU"
sampleName = "TTToHadronic_noPU"

# nevents = 200
# sample = "$t\\bar{t}$ PU 200"
# sampleName = "TTToHadronic_PU200"

dRfix = 5.0
pTfix = 5.0

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

with open(sampleName + "_cutScan.txt") as f:
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

dR = dR[0:18]
for i in range(len(dReff)):
  dReff[i] = dReff[i][0:18]
  dRdxy[i] = dRdxy[i][0:18]
  dRdxyerr[i] = dRdxyerr[i][0:18]
  dRd3d[i] = dRd3d[i][0:18]
  dRd3derr[i] = dRd3derr[i][0:18]

pT = pT[18:36]
for i in range(len(pTeff)):
  pTeff[i] = pTeff[i][18:36]
  pTdxy[i] = pTdxy[i][18:36]
  pTdxyerr[i] = pTdxyerr[i][18:36]
  pTd3d[i] = pTd3d[i][18:36]
  pTd3derr[i] = pTd3derr[i][18:36]

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

labels = ["GV-SV", "GV-SVt", "GVs-SV", "GVs-SVt", "GVn-SV", "GVn-SVt", "GVns-SV", "GVns-SVt"]

plt.figure()
for i in range(8):
  plt.plot(dR, dReff[i], marker=".", linestyle="-", markersize=15, label=labels[i])
plt.legend()
plt.xlim(-0.1, pTfix + 0.1)
plt.ylim(0.0, 1.2)
plt.xlabel("$\Delta R$ cut [cm]")
plt.ylabel("Matching efficiency at $p_{T}$ ratio cut = %.1f" % pTfix)
hep.style.use("CMS")
hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
plt.savefig(sampleName + "_dReff.pdf", bbox_inches="tight")
plt.close()

plt.figure()
for i in range(8):
  plt.errorbar(dR, dRdxy[i], yerr=dRdxyerr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
plt.legend()
plt.xlim(-0.1, pTfix + 0.1)
plt.xlabel("$\Delta R$ cut [cm]")
plt.ylabel("Match $d_{xy}$ [cm] at $p_{T}$ ratio cut = %.1f" % pTfix)
hep.style.use("CMS")
hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
plt.savefig(sampleName + "_dRdxy.pdf", bbox_inches="tight")
plt.close()

plt.figure()
for i in range(8):
  plt.errorbar(dR, dRd3d[i], yerr=dRd3derr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
plt.legend()
plt.xlim(-0.1, pTfix + 0.1)
plt.xlabel("$\Delta R$ cut [cm]")
plt.ylabel("Match $d_{3D}$ [cm] at $p_{T}$ ratio cut = %.1f" % pTfix)
hep.style.use("CMS")
hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
plt.savefig(sampleName + "_dRd3d.pdf", bbox_inches="tight")
plt.close()

plt.figure()
for i in range(8):
  plt.plot(pT, pTeff[i], marker=".", linestyle="-", markersize=15, label=labels[i])
plt.legend()
plt.xlim(-0.1, dRfix + 0.1)
plt.ylim(0.0, 1.2)
plt.xlabel("$p_{T}$ ratio cut")
plt.ylabel("Matching efficiency at $\Delta R$ cut = %.1f cm" % dRfix)
hep.style.use("CMS")
hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
plt.savefig(sampleName + "_pTeff.pdf", bbox_inches="tight")
plt.close()

plt.figure()
for i in range(8):
  plt.errorbar(pT, pTdxy[i], yerr=pTdxyerr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
plt.legend()
plt.xlim(-0.1, dRfix + 0.1)
plt.xlabel("$p_{T}$ ratio cut")
plt.ylabel("Match $d_{xy}$ [cm] at $\Delta R$ cut = %.1f cm" % dRfix)
hep.style.use("CMS")
hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
plt.savefig(sampleName + "_pTdxy.pdf", bbox_inches="tight")
plt.close()

plt.figure()
for i in range(8):
  plt.errorbar(pT, pTd3d[i], yerr=pTd3derr[i], marker=".", linestyle="-", markersize=15, label=labels[i])
plt.legend()
plt.xlim(-0.1, dRfix + 0.1)
plt.xlabel("$p_{T}$ ratio cut")
plt.ylabel("Match $d_{3D}$ [cm] at $\Delta R$ cut = %.1f cm" % dRfix)
hep.style.use("CMS")
hep.cms.label("Private Work", data=False, rlabel="%s, %d events" % (sample, nevents), fontsize=22)
plt.savefig(sampleName + "_pTd3d.pdf", bbox_inches="tight")
plt.close()
