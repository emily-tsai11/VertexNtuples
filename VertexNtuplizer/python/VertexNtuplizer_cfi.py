import FWCore.ParameterSet.Config as cms

vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
  tracks = cms.untracked.InputTag("generalTracks")
)
