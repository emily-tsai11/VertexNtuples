from CRABClient.UserUtilities import config, getUsername
config = config()

sample = "TTToHadronic_noPU"
# sample = "TTToHadronic_PU200"

config.General.requestName       = sample
config.General.transferOutputs   = True
config.General.instance          = "prod"

config.JobType.pluginName        = "Analysis"
config.JobType.psetName          = "VertexNtuplizer.py"
config.JobType.maxMemoryMB       = 3000
config.JobType.maxJobRuntimeMin  = 60
# config.JobType.pyCfgParams       = ["run2=True"]

config.Data.splitting            = "FileBased"
config.Data.unitsPerJob          = 1
config.Data.outLFNDirBase        = "/store/group/phys_btag/%s" % getUsername()
config.Data.publication          = False
config.Data.outputDatasetTag     = sample
config.Data.userInputFiles       = open(sample + ".list").readlines()
config.Data.outputPrimaryDataset = "VertexNtuples"

config.Site.storageSite          = "T2_CH_CERN"
config.Site.whitelist            = ["T2_CH_CERN"]
