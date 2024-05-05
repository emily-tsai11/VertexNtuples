# Parse job options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register("inputFiles",
    "file:/eos/user/e/etsai/workspace/TTToHadronic_noPU_CMSSW_13_1_0/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_noPU_slimmed.root",
    VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Input file(s) (default is ttbar no pileup)")
options.register("outputFile",  "histo", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Output file (w/o .root)")
options.register("maxEvents",   -1,      VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,    "Maximum number of events")
options.register("reportEvery", 1,       VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,    "Report every N events")
options.register("run2",        False,   VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,   "Sets |eta| to max 2.5 for 2018 tracker coverage")
if hasattr(sys, "argv"):
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


# Run process
release = os.environ["CMSSW_VERSION"][6:]
print("Using release " + release)
process.p = cms.Path(process.vertexNtuplizer)
