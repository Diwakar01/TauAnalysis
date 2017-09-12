import FWCore.ParameterSet.Config as cms

process = cms.Process("producePFTauIdEffNtupleFromMiniAOD")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('76X_mcRun2_asymptotic_v12')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

flist = []
#signal
#for ii in range(1,15):                                                                                              
#    flist.append('/store/user/bluj/76XMiniAODv2PFTau/v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1/160426_121904/0000/miniAOD-prod_PAT_' + str(ii) + '.root') 

#for ii in range(1,12):
#    flist.append('/store/cmst3/user/ytakahas/TauIDPAS/Zp1_run2/miniAOD-prod_PAT_' + str(ii) + '.root')
#for ii in range(1,7):
#    if ii in [3]: continue
#    flist.append('/store/cmst3/user/ytakahas/TauIDPAS/Zp4p5_run2/miniAOD-prod_PAT_' + str(ii) + '.root')
#for ii in range(1,32):
#    flist.append('/store/cmst3/user/ytakahas/TauIDPAS/DY_v1p1_run2/miniAOD-prod_PAT_' + str(ii) + '.root')
#for ii in range(1,12):
#    flist.append('/store/cmst3/user/ytakahas/TauIDPAS/GGH_run2/miniAOD-prod_PAT_' + str(ii) + '.root') 

#bkg
for ii in range(1,52):
    if ii in [3, 9, 12, 18, 33, 36, 48]: continue
    flist.append('/store/cmst3/user/ytakahas/TauIDPAS/QCD_run2/miniAOD-prod_PAT_' + str(ii) + '.root')

#bkg
#for ii in range(0,17):                                                                                              
#    if ii in [0, 3, 9, 12]: continue 
#    flist.append('/store/user/bluj/76XMiniAODv2PFTau/v1/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to3000_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1/160426_122846/0000/miniAOD-prod_PAT_' +str(ii)+'.root')



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/002ABFCA-A0B9-E511-B9BA-0CC47A57CD6A.root'
        #'root:://cms-xrd-global.cern.ch///store/relval/CMSSW_7_6_2/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/76X_mcRun2_asymptotic_v12-v1/00000/549E843B-979C-E511-A22A-0025905A6092.root',
        #'root:://cms-xrd-global.cern.ch///store/relval/CMSSW_7_6_2/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/76X_mcRun2_asymptotic_v12-v1/00000/7488793C-979C-E511-90DB-0026189438A7.root'
        #'/store/user/bluj/76XMiniAODv2PFTau/v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1/160426_121904/0000/miniAOD-prod_PAT_1.root'
        #'/store/user/bluj/76XMiniAODv2PFTau/v1/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to3000_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1/160426_122846/0000/miniAOD-prod_PAT_1.root'
        #'/store/cmst3/user/ytakahas/TauIDPAS/Zp2_run2/miniAOD-prod_PAT_1.root'
        flist
    ),
##    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
##    inputCommands=cms.untracked.vstring(
##        'keep *',
##        'drop patTaus_*_*_*',
##        'drop *PFTau*_*_*_*'
##    )
)
 
#####################################################
  
process.producePFTauIdEffNtupleFromMiniAODSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("prunedGenParticles"),
    maxEventsToPrint = cms.untracked.int32(1)
)
process.producePFTauIdEffNtupleFromMiniAODSequence += process.printGenParticleList
#--------------------------------------------------------------------------------

process.elecVeto  = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt > 15. && abs(eta) < 2.5 && pfIsolationVariables.sumChargedHadronPt < 0.2 && full5x5_sigmaIetaIeta < 0.0532 && abs(deltaEtaSuperClusterTrackAtVtx) < 0.0152 && abs(deltaPhiSuperClusterTrackAtVtx) < 0.237 && hcalOverEcal < 0.181 && passConversionVeto > 0"),
    filter = cms.bool(False)
)
process.producePFTauIdEffNtupleFromMiniAODSequence += process.elecVeto

process.pfTauIdEffNtupleFromMiniAODProducer = cms.EDProducer("PFTauIdEffMiniAODNtupleProducer",
    srcGenParticles = cms.InputTag('prunedGenParticles'),
    srcGenJets = cms.InputTag('slimmedGenJets'),
    srcRecTaus = cms.InputTag('slimmedTaus'),
    srcRecJets = cms.InputTag('slimmedJets'),                                                   
    srcRecVetoElectrons = cms.InputTag("elecVeto"),
    srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    srcWeights = cms.VInputTag()                                               
)
process.producePFTauIdEffNtupleFromMiniAODSequence += process.pfTauIdEffNtupleFromMiniAODProducer

process.p = cms.Path(process.producePFTauIdEffNtupleFromMiniAODSequence)
print process.p

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("pfTauIdEffNtupleFromMiniAOD_Fall15.root")
)

processDumpFile = open('producePFTauIdEffNtupleFromMiniAOD.dump', 'w')
print >> processDumpFile, process.dumpPython()




