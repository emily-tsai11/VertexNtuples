import FWCore.ParameterSet.Config as cms

vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
  # Collections
  genParticles                    = cms.untracked.InputTag("mergedGenParticles",                       "",                    "BTV"),
  simTracks                       = cms.untracked.InputTag("g4SimHits",                                "",                    "SIM"),
  generalTracks                   = cms.untracked.InputTag("generalTracks",                            "",                    "BTV"),
  pfCandidates                    = cms.untracked.InputTag("particleFlow",                             "",                    "BTV"),
  trackT0FromBS                   = cms.untracked.InputTag("trackExtenderWithMTD",                     "generalTrackt0",      "BTV"),
  trackSigmaT0FromBS              = cms.untracked.InputTag("trackExtenderWithMTD",                     "generalTracksigmat0", "BTV"),
  # trackQualityFromBS              = cms.untracked.InputTag("mtdTrackQualityMVA",                       "mtdQualMVA",          "BTV"),
  trackT0FromPV                   = cms.untracked.InputTag("trackExtenderFromPointWithMTD",            "generalTrackt0",      "BTV"),
  trackSigmaT0FromPV              = cms.untracked.InputTag("trackExtenderFromPointWithMTD",            "generalTracksigmat0", "BTV"),
  # trackQualityFromPV              = cms.untracked.InputTag("mtdTrackQualityMVA",                       "mtdQualMVA",          "BTV"), # Doesn't exist right now
  primaryVertices                 = cms.untracked.InputTag("offlineSlimmedPrimaryVertices",            "",                    "BTV"),
  nIVFClusters                    = cms.untracked.InputTag("inclusiveCandidateVertexFinder",           "nClusters",           "BTV"),
  nIVFClustersMTDPV               = cms.untracked.InputTag("inclusiveCandidateVertexFinderMTDPV",      "nClusters",           "BTV"),
  inclusiveSecondaryVertices      = cms.untracked.InputTag("inclusiveCandidateVertexFinder",           "",                    "BTV"),
  inclusiveSecondaryVerticesMTDPV = cms.untracked.InputTag("inclusiveCandidateVertexFinderMTDPV",      "",                    "BTV"),
  mergedSecondaryVertices         = cms.untracked.InputTag("inclusiveCandidateSecondaryVertices",      "",                    "BTV"),
  mergedSecondaryVerticesMTDPV    = cms.untracked.InputTag("inclusiveCandidateSecondaryVerticesMTDPV", "",                    "BTV"),
  slimmedSecondaryVertices        = cms.untracked.InputTag("slimmedSecondaryVertices",                 "",                    "BTV"),
  slimmedSecondaryVerticesMTDPV   = cms.untracked.InputTag("slimmedSecondaryVerticesMTDPV",            "",                    "BTV"),
  jets                            = cms.untracked.InputTag("slimmedJets",                              "",                    "BTV"),
  genJets                         = cms.untracked.InputTag("slimmedGenJets",                           "",                    "BTV"),
  genJetsFlavourInfo              = cms.untracked.InputTag("slimmedGenJetsFlavourInfos",               "",                    "BTV"),

  # Kinematic cuts
  genMotherPtMin                  = cms.untracked.double(10.0),
  trkPtMin                        = cms.untracked.double(0.8), # also for daughters
  jetPtMin                        = cms.untracked.double(20.0),
  jetPtMax                        = cms.untracked.double(1000.0),
  absEtaMax                       = cms.untracked.double(3.0), # 2.5: 2018 coverage, 3.0: MTD coverage

  # Fit cuts
  trkTimeQualityCut               = cms.untracked.double(0.5), # (Verbally) Recommended by the MTD DPG
  svChi2dofMax                    = cms.untracked.double(10.0),

  # Matching cuts
  trkMatchDrCut                   = cms.untracked.double(0.05),
  trkMatchPtCut                   = cms.untracked.double(0.1),
  vtxMatchFrac                    = cms.untracked.double(0.66),
  jetRadius                       = cms.untracked.double(0.4),
  jetMatchDrCut                   = cms.untracked.double(0.1),
)
