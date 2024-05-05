import FWCore.ParameterSet.Config as cms

vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
  genParticles     = cms.untracked.InputTag("prunedGenParticles", "", "BTV"),
  simTracks        = cms.untracked.InputTag("g4SimHits", "", "SIM"),
  tracks           = cms.untracked.InputTag("generalTracks", "", "BTV"),

  # Kinematic cuts
  absEtaMax        = cms.untracked.double(3.0), # MTD coverage
  genMotherPtMin   = cms.untracked.double(10.0),
  genDaughterPtMin = cms.untracked.double(1.0),
  trkPtMin         = cms.untracked.double(0.8),
  jetPtMin         = cms.untracked.double(20.0),
  jetPtMax         = cms.untracked.double(1000.0),

  # Matching cuts
  trkMatchDrCut    = cms.untracked.double(0.05),
  trkMatchPtCut    = cms.untracked.double(0.1),
  vtxMatchFrac     = cms.untracked.double(0.66),
)
