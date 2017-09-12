import FWCore.ParameterSet.Config as cms

process = cms.Process("producePFTauIdEffNtuple2")

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

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://xrootd.ba.infn.it//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        #'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/AODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/125A6B71-C56A-E411-9D2B-0025907609BE.root'
        #'file:/nfs/dust/cms/user/anayak/CMS/Ntuple_Phys14TauId/PickEvents_JetPt400/pickevents_merged_QCD_Phys14.root'
        #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02B3B2D1-BC13-E511-A895-008CFA110B10.root'
        'root://cms-xrd-global.cern.ch///store/mc/RunIIFall15DR76/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/003220F0-CDA8-E511-817B-5065F3816291.root'
        #'file:copyEvent/pickevents_merged.root'
    ),
##    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
##    inputCommands=cms.untracked.vstring(
##        'keep *',
##        'drop patTaus_*_*_*',
##        'drop *PFTau*_*_*_*'
##    )
)
 
#####################################################
  
process.producePFTauIdEffNtuple2Sequence = cms.Sequence()

#--------------------------------------------------------------------------------
# print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(1)
)
process.producePFTauIdEffNtuple2Sequence += process.printGenParticleList
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# rerun tau reconstruction with latest tags

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
#process.loadRecoTauTagMVAsFromPrepDB2.connect = cms.string('sqlite_file:RecoTauTag_MVAs_2015Sep23.db')
process.producePFTauIdEffNtuple2Sequence += process.PFTau
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

process.producePFTauIdEffNtuple2Sequence += process.selectPrimaryVertex
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# compute event weights for pile-up reweighting
# (Summer'12 MC to 2012 run A data)

##from TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi import vertexMultiplicityReweight
##process.vertexMultiplicityReweight3d2012RunABC = vertexMultiplicityReweight.clone(
##    inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs190456to208686_Mu17_Mu8.root"),
##    type = cms.string("gen3d"),
##    mcPeriod = cms.string("Summer12_S10")
##)
##process.producePFTauIdEffNtuple2Sequence += process.vertexMultiplicityReweight3d2012RunABC
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")

# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute' )
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    jetSource = cms.InputTag('ak4PFJets'),
    jetCorrections = ( 'AK4PF', jetCorrections, "" ),
    outputModules = []
)

# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
   
process.producePFTauIdEffNtuple2Sequence += process.patDefaultSequence
#--------------------------------------------------------------------------------
simpleCutsWP95 = \
    "(userFloat('nHits') <= 999 " + \
    "&& ( (isEB && userFloat('sihih') < 0.01 && userFloat('dPhi') < 0.8 " + \
    "&&            userFloat('dEta') < 0.007 && userFloat('HoE') < 0.15) " + \
    "||   " + \
    "     (isEE && userFloat('sihih') < 0.03 && userFloat('dPhi') < 0.7 " + \
    "&&            userFloat('dEta') < 0.01 && userFloat('HoE') < 999) ))"

process.elecVeto  = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("selectedPatElectrons"),
    cut = cms.string(simpleCutsWP95 + " && pt > 15. && abs(eta) < 2.5 && userFloat('PFRelIsoDB04v3') < 0.3 && abs(userFloat('dzWrtPV')) < 0.2"),
    filter = cms.bool(False)
)
process.producePFTauIdEffNtuple2Sequence += process.elecVeto

process.pfTauIdEffNtuple2Producer = cms.EDProducer("PFTauIdEffNtupleProducer2",
    srcGenParticles = cms.InputTag('genParticles'),
    srcGenJets = cms.InputTag('ak4GenJets'),                                                
    srcRecTaus = cms.InputTag('patTaus'),
    srcRecJets = cms.InputTag('patJets'),                                                   
    srcRecVetoElectrons = cms.InputTag("elecVeto"),
    #srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    srcTracks = cms.InputTag('generalTracks'),
    srcBeamSpot = cms.InputTag('offlineBeamSpot'),
    ##srcWeights = cms.VInputTag('vertexMultiplicityReweight3d2012RunABC')
    srcWeights = cms.VInputTag()                                               
)
process.producePFTauIdEffNtuple2Sequence += process.pfTauIdEffNtuple2Producer

process.p = cms.Path(process.producePFTauIdEffNtuple2Sequence)
print process.p

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("pfTauIdEffNtuple2_Fall15.root")
)

processDumpFile = open('producePFTauIdEffNtuple2.dump', 'w')
print >> processDumpFile, process.dumpPython()




