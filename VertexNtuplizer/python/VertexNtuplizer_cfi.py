import FWCore.ParameterSet.Config as cms

vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
  # Collections
  genParticles            = cms.untracked.InputTag("mergedGenParticles",            "",                    "BTV"),
  simTracks               = cms.untracked.InputTag("g4SimHits",                     "",                    "SIM"),
  generalTracks           = cms.untracked.InputTag("generalTracks",                 "",                    "BTV"),
  particleFlowCandidates  = cms.untracked.InputTag("particleFlow",                  "",                    "BTV"),
  trackTimeBSValueMap     = cms.untracked.InputTag("trackExtenderWithMTD",          "generalTrackt0",      "BTV"),
  trackTimeBSErrorMap     = cms.untracked.InputTag("trackExtenderWithMTD",          "generalTracksigmat0", "BTV"),
  trackTimeBSQualityMap   = cms.untracked.InputTag("mtdTrackQualityMVA",            "mtdQualMVA",          "BTV"),
  trackTimePVValueMap     = cms.untracked.InputTag("trackExtenderWithMTDPV",        "generalTrackt0",      "BTV"),
  trackTimePVErrorMap     = cms.untracked.InputTag("trackExtenderWithMTDPV",        "generalTracksigmat0", "BTV"),
  trackTimePVQualityMap   = cms.untracked.InputTag("mtdTrackQualityMVA",            "mtdQualMVA",          "BTV"), # Doesn't exist right now
  primaryVertices         = cms.untracked.InputTag("offlineSlimmedPrimaryVertices", "",                    "BTV"),
  nIVFClusters            = cms.untracked.InputTag("inclusiveVertexFinder",         "nClusters",           "BTV"),
  nIVFClustersMTDBS       = cms.untracked.InputTag("inclusiveVertexFinderMTDBS",    "nClusters",           "BTV"),
  nIVFClustersMTDBS4      = cms.untracked.InputTag("inclusiveVertexFinderMTDBS4",   "nClusters",           "BTV"),
  nIVFClustersMTDPV       = cms.untracked.InputTag("inclusiveVertexFinderMTDPV",    "nClusters",           "BTV"),
  secondaryVertices       = cms.untracked.InputTag("inclusiveVertexFinder",         "",                    "BTV"), # Without MTD information
  secondaryVerticesMTDBS  = cms.untracked.InputTag("inclusiveVertexFinderMTDBS",    "",                    "BTV"),
  secondaryVerticesMTDBS4 = cms.untracked.InputTag("inclusiveVertexFinderMTDBS4",   "",                    "BTV"),
  secondaryVerticesMTDPV  = cms.untracked.InputTag("inclusiveVertexFinderMTDPV",    "",                    "BTV"),
  slimmedCandSVs          = cms.untracked.InputTag("slimmedSecondaryVertices",      "",                    "BTV"),
  jets                    = cms.untracked.InputTag("slimmedJets",                   "",                    "BTV"),
  genJets                 = cms.untracked.InputTag("slimmedGenJets",                "",                    "BTV"),
  genJetsFlavourInfo      = cms.untracked.InputTag("slimmedGenJetsFlavourInfos",    "",                    "BTV"),

  # Kinematic cuts
  absEtaMax               = cms.untracked.double(3.0), # 2.5: 2018 coverage, 3.0: MTD coverage
  genMotherPtMin          = cms.untracked.double(10.0),
  genDaughterPtMin        = cms.untracked.double(0.8),
  trkPtMin                = cms.untracked.double(0.8),
  jetPtMin                = cms.untracked.double(20.0),
  jetPtMax                = cms.untracked.double(1000.0),

  # Fit cuts
  trkTimeQualityCut       = cms.untracked.double(0.5), # (Verbally) Recommended by the MTD DPG
  svChi2dofMax            = cms.untracked.double(10.0),

  # Matching cuts
  trkMatchDrCut           = cms.untracked.double(0.05),
  trkMatchPtCut           = cms.untracked.double(0.1),
  vtxMatchFrac            = cms.untracked.double(0.66),
  jetRadius               = cms.untracked.double(0.4),
  jetMatchDrCut           = cms.untracked.double(0.1),
)
