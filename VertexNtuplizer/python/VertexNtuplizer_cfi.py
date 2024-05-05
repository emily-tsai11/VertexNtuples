import FWCore.ParameterSet.Config as cms

vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
  # Collections
  genParticles               = cms.untracked.InputTag("prunedGenParticles",             "",           "BTV"),
  simTracks                  = cms.untracked.InputTag("g4SimHits",                      "",           "SIM"),
  primaryVertices            = cms.untracked.InputTag("offlinePrimaryVertices",         "",           "BTV"),
  secondaryVertices          = cms.untracked.InputTag("inclusiveVertexFinder",          "",           "BTV"),
  secondaryVerticesMTDTiming = cms.untracked.InputTag("inclusiveVertexFinderMTDTiming", "",           "BTV"),
  IVFclusters                = cms.untracked.InputTag("inclusiveVertexFinder",          "nClusters",  "BTV"),
  IVFclustersMTDTiming       = cms.untracked.InputTag("inclusiveVertexFinderMTDTiming", "nClusters",  "BTV"),
  trackTimeValueMap          = cms.untracked.InputTag("tofPID",                         "t0",         "BTV"),
  trackTimeErrorMap          = cms.untracked.InputTag("tofPID",                         "sigmat0",    "BTV"),
  trackTimeQualityMap        = cms.untracked.InputTag("mtdTrackQualityMVA",             "mtdQualMVA", "BTV"),
  jets                       = cms.untracked.InputTag("slimmedJets",                    "",           "BTV"),
  genJets                    = cms.untracked.InputTag("slimmedGenJets",                 "",           "BTV"),
  genJetsFlavourInfo         = cms.untracked.InputTag("slimmedGenJetsFlavourInfos",     "",           "BTV"),

  # Kinematic cuts
  absEtaMax                  = cms.untracked.double(3.0), # MTD coverage
  genMotherPtMin             = cms.untracked.double(10.0),
  genDaughterPtMin           = cms.untracked.double(1.0),
  trkPtMin                   = cms.untracked.double(0.8),
  jetPtMin                   = cms.untracked.double(20.0),
  jetPtMax                   = cms.untracked.double(1000.0),

  # Matching cuts
  trkMatchDrCut              = cms.untracked.double(0.05),
  trkMatchPtCut              = cms.untracked.double(0.1),
  vtxMatchFrac               = cms.untracked.double(0.66),
  trkTimeQualityCut          = cms.untracked.double(0.5), # Recommended by MTD DPG (verbally)
  vtxMatchDrCut              = cms.untracked.double(0.1),
  jetMatchDrCut              = cms.untracked.double(0.1),
)
