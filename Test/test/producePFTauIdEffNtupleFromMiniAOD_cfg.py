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
process.GlobalTag.globaltag = cms.string('92X_upgrade2017_realistic_v10')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

flist = []
#signal
#for ii in range(1,15):                                                                                              
#    flist.append('/store/user/bluj/76XMiniAODv2PFTau/v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1/160426_121904/0000/miniAOD-prod_PAT_' + str(ii) + '.root') 
#bkg
#for ii in range(0,17):                                                                                              
#    if ii in [0, 3, 9, 12]: continue 
#    flist.append('/store/user/bluj/76XMiniAODv2PFTau/v1/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to3000_76X_mcRun2_asymptotic_v13_MiniAODv2_PFTau_v1/160426_122846/0000/miniAOD-prod_PAT_' +str(ii)+'.root')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #flist
        '/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/0A8639B4-1D94-E711-9C68-02163E0135C6.root'
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




