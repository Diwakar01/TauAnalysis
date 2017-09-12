import FWCore.ParameterSet.Config as cms

process = cms.Process('dumpTriggerInfo')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/user/veelken/CMSSW_5_3_x/skims/data_TauPlusX_2012runD_AOD_1_1_29z.root'                        
        '/store/user/veelken/CMSSW_5_3_x/skims/simQCDmuEnrichedPt470to600_AOD_1_1_A8V.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.dumpTriggerInfo = cms.EDAnalyzer("DumpTriggerInfo",
    srcL1GtReadoutRecord = cms.InputTag('gtDigis::RECO'),
    srcHLTresults = cms.InputTag('TriggerResults::HLT')
)

process.p = cms.Path(process.dumpTriggerInfo)
