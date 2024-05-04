import FWCore.ParameterSet.Config as cms

process = cms.Process("VertexNtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("VertexNtuples.VertexNtuplizer.VertexNtuplizer_cfi")

# Add arg parser

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "file:/eos/user/e/etsai/workspace/TTToHadronic_noPU_CMSSW_13_1_0/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_noPU_slimmed.root"
  )
)

process.p = cms.Path(process.vertexNtuplizer)
