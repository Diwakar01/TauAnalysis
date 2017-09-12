from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DYJetsToLL_M-50_13TeV_Fall15_25ns_miniAODv2_v3'
config.General.workArea = 'DYJetsLL'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'producePFTauIdEffNtupleFromMiniAOD_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/bluj-DYJetsToLL_M-50_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1-35a6577e52abd393294705164f2fba84/USER'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/bluj-DYJetsToLL_M-50_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1p1-35a6577e52abd393294705164f2fba84/USER'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
#config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/anayak/TauIDMVATraining2015/Validation/DYJetsLL_Fall15_17Feb2016/'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
