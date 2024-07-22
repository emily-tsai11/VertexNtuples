import FWCore.ParameterSet.Config as cms

process = cms.Process("VertexNtuplizer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/user/e/etsai/workspace/TTToHadronic_noPU_CMSSW_13_1_0/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_noPU_slimmed.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.maxLuminosityBlocks = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

process.vertexNtuplizer = cms.EDAnalyzer("VertexNtuplizer",
    absEtaMax = cms.untracked.double(2.5),
    genDaughterPtMin = cms.untracked.double(0.8),
    genJets = cms.untracked.InputTag("slimmedGenJets","","BTV"),
    genJetsFlavourInfo = cms.untracked.InputTag("slimmedGenJetsFlavourInfos","","BTV"),
    genMotherPtMin = cms.untracked.double(10.0),
    genParticles = cms.untracked.InputTag("mergedGenParticles","","BTV"),
    generalTracks = cms.untracked.InputTag("generalTracks","","BTV"),
    jetMatchDrCut = cms.untracked.double(0.1),
    jetPtMax = cms.untracked.double(1000.0),
    jetPtMin = cms.untracked.double(20.0),
    jetRadius = cms.untracked.double(0.4),
    jets = cms.untracked.InputTag("slimmedJets","","BTV"),
    nIVFClusters = cms.untracked.InputTag("inclusiveVertexFinder","nClusters","BTV"),
    nIVFClustersMTDBS = cms.untracked.InputTag("inclusiveVertexFinderMTDBS","nClusters","BTV"),
    nIVFClustersMTDBS4 = cms.untracked.InputTag("inclusiveVertexFinderMTDBS4","nClusters","BTV"),
    nIVFClustersMTDPV = cms.untracked.InputTag("inclusiveVertexFinderMTDPV","nClusters","BTV"),
    particleFlowCandidates = cms.untracked.InputTag("particleFlow","","BTV"),
    primaryVertices = cms.untracked.InputTag("offlineSlimmedPrimaryVertices","","BTV"),
    secondaryVertices = cms.untracked.InputTag("inclusiveVertexFinder","","BTV"),
    secondaryVerticesMTDBS = cms.untracked.InputTag("inclusiveVertexFinderMTDBS","","BTV"),
    secondaryVerticesMTDBS4 = cms.untracked.InputTag("inclusiveVertexFinderMTDBS4","","BTV"),
    secondaryVerticesMTDPV = cms.untracked.InputTag("inclusiveVertexFinderMTDPV","","BTV"),
    simTracks = cms.untracked.InputTag("g4SimHits","","SIM"),
    slimmedCandSVs = cms.untracked.InputTag("slimmedSecondaryVertices","","BTV"),
    svChi2dofMax = cms.untracked.double(10.0),
    trackTimeBSErrorMap = cms.untracked.InputTag("trackExtenderWithMTD","generalTracksigmat0","BTV"),
    trackTimeBSQualityMap = cms.untracked.InputTag("mtdTrackQualityMVA","mtdQualMVA","BTV"),
    trackTimeBSValueMap = cms.untracked.InputTag("trackExtenderWithMTD","generalTrackt0","BTV"),
    trackTimePVErrorMap = cms.untracked.InputTag("trackExtenderWithMTDPV","generalTracksigmat0","BTV"),
    trackTimePVQualityMap = cms.untracked.InputTag("mtdTrackQualityMVA","mtdQualMVA","BTV"),
    trackTimePVValueMap = cms.untracked.InputTag("trackExtenderWithMTDPV","generalTrackt0","BTV"),
    trkMatchDrCut = cms.untracked.double(0.05),
    trkMatchPtCut = cms.untracked.double(0.1),
    trkPtMin = cms.untracked.double(0.8),
    trkTimeQualityCut = cms.untracked.double(0.5),
    vtxMatchFrac = cms.untracked.double(0.66)
)


process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(100)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        enable = cms.untracked.bool(True),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.untracked.bool(False),
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.untracked.string('WARNING'),
        threshold = cms.untracked.string('INFO'),
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(False),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.optional.untracked.bool,
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.optional.untracked.string,
        threshold = cms.optional.untracked.string,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    debugModules = cms.untracked.vstring(),
    default = cms.untracked.PSet(
        limit = cms.optional.untracked.int32,
        lineLength = cms.untracked.int32(80),
        noLineBreaks = cms.untracked.bool(False),
        noTimeStamps = cms.untracked.bool(False),
        reportEvery = cms.untracked.int32(1),
        statisticsThreshold = cms.untracked.string('INFO'),
        threshold = cms.untracked.string('INFO'),
        timespan = cms.optional.untracked.int32,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    files = cms.untracked.PSet(
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            enableStatistics = cms.untracked.bool(False),
            extension = cms.optional.untracked.string,
            filename = cms.optional.untracked.string,
            lineLength = cms.optional.untracked.int32,
            noLineBreaks = cms.optional.untracked.bool,
            noTimeStamps = cms.optional.untracked.bool,
            output = cms.optional.untracked.string,
            resetStatistics = cms.untracked.bool(False),
            statisticsThreshold = cms.optional.untracked.string,
            threshold = cms.optional.untracked.string,
            allowAnyLabel_=cms.optional.untracked.PSetTemplate(
                limit = cms.optional.untracked.int32,
                reportEvery = cms.untracked.int32(1),
                timespan = cms.optional.untracked.int32
            )
        )
    ),
    suppressDebug = cms.untracked.vstring(),
    suppressFwkInfo = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    allowAnyLabel_=cms.optional.untracked.PSetTemplate(
        limit = cms.optional.untracked.int32,
        reportEvery = cms.untracked.int32(1),
        timespan = cms.optional.untracked.int32
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('/eos/cms/store/group/phys_btag/etsai/VertexNtuples/preliminary/histo.root')
)


process.p = cms.Path(process.vertexNtuplizer)


