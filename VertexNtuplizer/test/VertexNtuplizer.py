# Parse job options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register("inputFiles",
  "file:/eos/user/e/etsai/workspace/TTToHadronic_noPU_CMSSW_13_1_0/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_noPU_slimmed.root",
  # "file:/eos/cms/store/user/etsai/TT_TuneCP5_14TeV-powheg-pythia8/TTToHadronic_noPU/240429_225553/0000/TTToHadronic_noPU_slimmed_1.root",
  # "file:/eos/cms/store/group/phys_btag/etsai/TT_TuneCP5_14TeV-powheg-pythia8/TTToHadronic_PU200/240429_231224/0000/TTToHadronic_PU200_slimmed_1.root",
  VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Input file(s) (default is ttbar no pileup)")
options.register("outputFile",  "histo", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Output file (w/o .root)")
options.register("maxEvents",   10,      VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,    "Maximum number of events")
options.register("reportEvery", 100,     VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,    "Report every N events")
options.register("run2",        False,   VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,   "Sets |eta| to max 2.5 for 2018 tracker coverage")
options.register("scanCuts",    False,   VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,   "Scan pT and dR cuts")
options.parseArguments()


# Poor man's solution because the "RECREATE" option is not working to delete the output ROOT file for some reason
import os
if os.path.exists(options.outputFile + ".root"):
  os.remove(options.outputFile + ".root")


# Define process
import FWCore.ParameterSet.Config as cms
process = cms.Process("VertexNtuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile + ".root"))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))


# Load analyzer
process.load("VertexNtuples.VertexNtuplizer.VertexNtuplizer_cfi")
if options.run2:
  process.vertexNtuplizer.etaMax = 2.5


# Also create a process with loose GV to SV matching cuts for comparison
process.vertexNtuplizerLoose = process.vertexNtuplizer.clone(
  recoTrkMatchDrCut=4.0,
  recoTrkMatchPtCut=1.0
)


# Run process
release = os.environ["CMSSW_VERSION"][6:]
print("Using release " + release)
process.p = cms.Path(process.vertexNtuplizer * process.vertexNtuplizerLoose)


if options.scanCuts:
  print("Scanning 24 dR cuts and 15 pT cuts with %d events." % options.maxEvents)

  pTfix = 1.0
  process.vertexNtuplizer1 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.05, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer2 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.06, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer3 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.07, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer4 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.08, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer5 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.09, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer6 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.1, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer7 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.2, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer8 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.3, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer9 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.4, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer10 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.5, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer11 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.6, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer12 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.7, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer13 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.8, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer14 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=0.9, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer15 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=1.0, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer16 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=1.1, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer17 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=1.3, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer18 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=1.5, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer19 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=1.7, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer20 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=2.0, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer21 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=2.5, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer22 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=3.0, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer23 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=3.5, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)
  process.vertexNtuplizer24 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=4.0, recoTrkMatchPtCut=pTfix, scanCuts=options.scanCuts)

  dRfix = 4.0
  process.vertexNtuplizer25 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.05, scanCuts=options.scanCuts)
  process.vertexNtuplizer26 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.06, scanCuts=options.scanCuts)
  process.vertexNtuplizer27 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.07, scanCuts=options.scanCuts)
  process.vertexNtuplizer28 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.08, scanCuts=options.scanCuts)
  process.vertexNtuplizer29 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.09, scanCuts=options.scanCuts)
  process.vertexNtuplizer30 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.1, scanCuts=options.scanCuts)
  process.vertexNtuplizer31 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.2, scanCuts=options.scanCuts)
  process.vertexNtuplizer32 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.3, scanCuts=options.scanCuts)
  process.vertexNtuplizer33 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.4, scanCuts=options.scanCuts)
  process.vertexNtuplizer34 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.5, scanCuts=options.scanCuts)
  process.vertexNtuplizer35 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.6, scanCuts=options.scanCuts)
  process.vertexNtuplizer36 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.7, scanCuts=options.scanCuts)
  process.vertexNtuplizer37 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.8, scanCuts=options.scanCuts)
  process.vertexNtuplizer38 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=0.9, scanCuts=options.scanCuts)
  process.vertexNtuplizer39 = process.vertexNtuplizer.clone(recoTrkMatchDrCut=dRfix, recoTrkMatchPtCut=1.0, scanCuts=options.scanCuts)

  process.p = cms.Path(
      process.vertexNtuplizer1*process.vertexNtuplizer2*process.vertexNtuplizer3*process.vertexNtuplizer4*
      process.vertexNtuplizer5*process.vertexNtuplizer6*process.vertexNtuplizer7*process.vertexNtuplizer8*
      process.vertexNtuplizer9*process.vertexNtuplizer10*process.vertexNtuplizer11*process.vertexNtuplizer12*
      process.vertexNtuplizer13*process.vertexNtuplizer14*process.vertexNtuplizer15*process.vertexNtuplizer16*
      process.vertexNtuplizer17*process.vertexNtuplizer18*process.vertexNtuplizer19*process.vertexNtuplizer20*
      process.vertexNtuplizer21*process.vertexNtuplizer22*process.vertexNtuplizer23*process.vertexNtuplizer24*
      process.vertexNtuplizer25*process.vertexNtuplizer26*process.vertexNtuplizer27*process.vertexNtuplizer28*
      process.vertexNtuplizer29*process.vertexNtuplizer30*process.vertexNtuplizer31*process.vertexNtuplizer32*
      process.vertexNtuplizer33*process.vertexNtuplizer34*process.vertexNtuplizer35*process.vertexNtuplizer36*
      process.vertexNtuplizer37*process.vertexNtuplizer38*process.vertexNtuplizer39)


# Dump process config
open("dumpVertexNtuplizer.py", "w").write(process.dumpPython())
