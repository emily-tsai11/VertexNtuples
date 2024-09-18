# Parse job options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register("inputFiles",
  # "file:/eos/user/e/etsai/workspace/CMSSW_13_1_0/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_noPU_slimmed.root",
  # "file:/eos/user/e/etsai/workspace/CMSSW_13_1_3/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_PU200_slimmed.root",
  "file:/eos/cms/store/group/phys_btag/etsai/TT_TuneCP5_14TeV-powheg-pythia8/TTToHadronic_noPU/240917_090552/0000/TTToHadronic_noPU_slimmed_1.root",
  # "file:/eos/cms/store/group/phys_btag/etsai/TT_TuneCP5_14TeV-powheg-pythia8/TTToHadronic_PU200/240917_090637/0000/TTToHadronic_PU200_slimmed_1.root",
  VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Input file(s) (default is ttbar no pileup)")
options.register("outputFile",  "histo", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Output file (w/o .root)")
options.register("maxEvents",   -1,      VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,    "Maximum number of events")
options.register("reportEvery", 100,     VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,    "Report every N events")
options.register("run2",        True,    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,   "Sets |eta| to max 2.5 for 2018 tracker coverage")
options.parseArguments()
# file_path = "/eos/cms/store/group/phys_btag/etsai/VertexNtuples/preliminary/" + options.outputFile + ".root"
file_path = options.outputFile + ".root" # for CRAB submission


# Poor man's solution because the "RECREATE" option is not working to delete the existing ROOT file for some reason
import os
if os.path.exists(file_path):
  os.remove(file_path)


# Define process
import FWCore.ParameterSet.Config as cms
process = cms.Process("VertexNtuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.TFileService = cms.Service("TFileService", fileName=cms.string(file_path))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))


# Load analyzer
process.load("VertexNtuples.VertexNtuplizer.VertexNtuplizer_cfi")
if options.run2:
  process.vertexNtuplizer.absEtaMax = 2.5


# Run process
release = os.environ["CMSSW_VERSION"][6:]
print("Using release " + release)
process.p = cms.Path(process.vertexNtuplizer)


# Dump process config
open("dumpVertexNtuplizer.py", "w").write(process.dumpPython())
