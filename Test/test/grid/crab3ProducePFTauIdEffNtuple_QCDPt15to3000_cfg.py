from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'QCD_Pt-15to3000_Flat_13TeV_Fall15_25ns_miniAODv2_User_v1'
config.General.workArea = 'QCD_Pt-15to3000'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'producePFTauIdEffNtupleFromMiniAOD_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDataset = '/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/bluj-QCD_Pt-15to3000_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1-35a6577e52abd393294705164f2fba84/USER'
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/anayak/TauIDMVATraining2015/Validation/QCDFlat_Fall15_17Feb2016/'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
