import FWCore.ParameterSet.Config as cms

vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
  # Collections
  prunedGenParticles      = cms.untracked.InputTag("prunedGenParticles",             "",                    "BTV"),
  packedGenParticles      = cms.untracked.InputTag("packedGenParticles",             "",                    "BTV"),
  simTracks               = cms.untracked.InputTag("g4SimHits",                      "",                    "SIM"),
  # trackingParticles       = cms.untracked.InputTag("mix",                            "MergedTrackTruth",    "HLT"),
  trackingVertices        = cms.untracked.InputTag("mix",                            "MergedTrackTruth",    "HLT"),
  primaryVertices         = cms.untracked.InputTag("offlinePrimaryVertices",         "",                    "BTV"),

  nIVFClusters            = cms.untracked.InputTag("inclusiveVertexFinder",          "nClusters",           "BTV"),
  nIVFClustersMTDBS       = cms.untracked.InputTag("inclusiveVertexFinderMTDBS",     "nClusters",           "BTV"),
  # nIVFClustersMTDBS4      = cms.untracked.InputTag("inclusiveVertexFinderMTDTiming", "nClusters",           "BTV"),
  nIVFClustersMTDBS4      = cms.untracked.InputTag("inclusiveVertexFinderMTDBS4",    "nClusters",           "BTV"),
  nIVFClustersMTDPV       = cms.untracked.InputTag("inclusiveVertexFinderMTDPV",     "nClusters",           "BTV"),

  secondaryVertices       = cms.untracked.InputTag("inclusiveVertexFinder",          "",                    "BTV"), # without MTD information
  secondaryVerticesMTDBS  = cms.untracked.InputTag("inclusiveVertexFinderMTDBS",     "",                    "BTV"),
  # secondaryVerticesMTDBS4 = cms.untracked.InputTag("inclusiveVertexFinderMTDTiming", "",                    "BTV"),
  secondaryVerticesMTDBS4 = cms.untracked.InputTag("inclusiveVertexFinderMTDBS4",    "",                    "BTV"),
  secondaryVerticesMTDPV  = cms.untracked.InputTag("inclusiveVertexFinderMTDPV",     "",                    "BTV"),

  # trackTimeBSValueMap     = cms.untracked.InputTag("tofPID",                         "t0",                  "BTV"),
  # trackTimeBSErrorMap     = cms.untracked.InputTag("tofPID",                         "sigmat0",             "BTV"),
  # trackTimeBSQualityMap   = cms.untracked.InputTag("mtdTrackQualityMVA",             "mtdQualMVA",          "BTV"),
  # trackTimePVValueMap     = cms.untracked.InputTag("tofPID",                         "t0",                  "BTV"),
  # trackTimePVErrorMap     = cms.untracked.InputTag("tofPID",                         "sigmat0",             "BTV"),
  # trackTimePVQualityMap   = cms.untracked.InputTag("mtdTrackQualityMVA",             "mtdQualMVA",          "BTV"),

  trackTimeBSValueMap     = cms.untracked.InputTag("trackExtenderWithMTD",           "generalTrackt0",      "BTV"),
  trackTimeBSErrorMap     = cms.untracked.InputTag("trackExtenderWithMTD",           "generalTracksigmat0", "BTV"),
  trackTimeBSQualityMap   = cms.untracked.InputTag("mtdTrackQualityMVA",             "mtdQualMVA",          "BTV"),
  trackTimePVValueMap     = cms.untracked.InputTag("trackExtenderWithMTDPV",         "generalTrackt0",      "BTV"),
  trackTimePVErrorMap     = cms.untracked.InputTag("trackExtenderWithMTDPV",         "generalTracksigmat0", "BTV"),
  # trackTimePVQualityMap   = cms.untracked.InputTag("mtdTrackQualityMVA",             "mtdQualMVA",          "BTV"), # doesn't exist right now

  jets                    = cms.untracked.InputTag("slimmedJets",                    "",                    "BTV"),
  genJets                 = cms.untracked.InputTag("slimmedGenJets",                 "",                    "BTV"),
  genJetsFlavourInfo      = cms.untracked.InputTag("slimmedGenJetsFlavourInfos",     "",                    "BTV"),

  # Kinematic cuts
  absEtaMax               = cms.untracked.double(2.5), # 2.5: 2018 coverage, 3.0: MTD coverage
  genMotherPtMin          = cms.untracked.double(10.0),
  genDaughterPtMin        = cms.untracked.double(1.0),
  trkPtMin                = cms.untracked.double(0.8),
  jetPtMin                = cms.untracked.double(20.0),
  jetPtMax                = cms.untracked.double(1000.0),

  # Fit cuts
  trkTimeQualityCut       = cms.untracked.double(0.5), # Recommended by MTD DPG (verbally)
  svChi2dofMax            = cms.untracked.double(10.0),

  # Matching cuts
  simTrkMatchDrCut        = cms.untracked.double(0.05),
  simTrkMatchPtCut        = cms.untracked.double(0.1),
  recoTrkMatchDrCut       = cms.untracked.double(0.3),
  recoTrkMatchPtCut       = cms.untracked.double(0.2),
  vtxMatchFrac            = cms.untracked.double(0.66),
  jetRadius               = cms.untracked.double(0.4),
  jetMatchDrCut           = cms.untracked.double(0.1),
  scanCuts                = cms.untracked.bool(False),
)
