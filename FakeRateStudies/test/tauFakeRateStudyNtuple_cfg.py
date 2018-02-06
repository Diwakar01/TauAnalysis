import FWCore.ParameterSet.Config as cms

isData = True
#isData = False
runOnData=isData #data/MC switch

process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
if runOnData:
  process.GlobalTag.globaltag = '92X_dataRun2_Prompt_v9'
else:
  process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/data/Run2017D/JetHT/MINIAOD/PromptReco-v1/000/302/031/00000/4634BAB9-2F8F-E711-86BA-02163E014673.root'
    )
)

process.slimmedPatTriggerUnpacked = cms.EDProducer('PATTriggerObjectStandAloneUnpacker',
                                                   patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
                                                   triggerResults = cms.InputTag('TriggerResults::HLT'),
                                                   unpackFilterLabels = cms.bool(True)
                                                   )
                                                   
process.tauNtuple = cms.EDAnalyzer('FakeRateStudyNtupleProducer',
                                   PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                   MetCollectionTag = cms.InputTag('slimmedMETs'),
                                   MuonCollectionTag = cms.InputTag('slimmedMuons'),
                                   PFJetCollectionTag = cms.InputTag('slimmedJets'),
                                   TriggerObjectCollectionTag = cms.InputTag('slimmedPatTriggerUnpacked'),
                                   pileupTag = cms.InputTag('slimmedAddPileupInfo'),
                                   triggerFiltersTag = cms.InputTag('slimmedPatTrigger::filterLabels::PAT'),
                                   hltprocess = cms.InputTag('TriggerResults::HLT'),
                                   triggerPaths = cms.vstring('HLT_PFJet60_v',
                                                              'HLT_PFJet140_v',
                                                              'HLT_PFJet500_v',
                                                              'HLT_IsoMu24_eta2p1_v',
                                                              'HLT_IsoMu24_v'
                                                              ),
                                   triggerFilters = cms.vstring('hltSinglePFJet60',
                                                                'hltSinglePFJet140',
                                                                'hltSinglePFJet500',
                                                                'hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07',
                                                                'hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07'
                                                                )
                                   )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("tauNtuple.root")
)
process.p = cms.Path(process.slimmedPatTriggerUnpacked + process.tauNtuple)
