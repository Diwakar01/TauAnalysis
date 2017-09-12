from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'QCD_Pt-15to7000_Flat_13TeV_Fall15_25ns_miniAODv2_v3'
config.General.workArea = 'QCD_Pt-15to7000'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'producePFTauIdEffNtupleFromMiniAOD_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_magnetOn_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/anayak/TauIDMVATraining2015/Validation/QCDFlat_Fall15_17Feb2016/'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
