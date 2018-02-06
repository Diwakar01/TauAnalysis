// -*- C++ -*-
//
// Package:    FakeRateStudyNtupleProducer
// Class:      FakeRateStudyNtupleProducer
// 
/**\class FakeRateStudyNtupleProducer FakeRateStudyNtupleProducer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arun Nayak
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
//
// class declaration
//
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h" 
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <iostream>
#include <vector>
#include "TTree.h"
#include <string>


class FakeRateStudyNtupleProducer : public edm::EDAnalyzer {
   public:
      explicit FakeRateStudyNtupleProducer(const edm::ParameterSet&);
      ~FakeRateStudyNtupleProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      bool isICHEPMuon(const reco::Muon&);

      // ----------member data ---------------------------
  //edm::InputTag hlt_;
  edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
  edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;
  edm::EDGetTokenT<reco::VertexCollection> PVToken_;
  edm::EDGetTokenT<pat::JetCollection> PFJetCollectionToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken_;
  edm::EDGetTokenT<edm::TriggerResults> hlt_;
  edm::EDGetTokenT<std::vector<std::string> > triggerFiltersToken_;
  std::vector<std::string> triggerPaths_;
  std::vector<std::string> triggerFilters_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >pileupToken_;

  TTree* tree_;
  unsigned long run_,event_,lumi_;
  float metPx_, metPy_, metPt_, metPhi_;
  Float_t npuIT, npuOOT, npuTrue;
  UInt_t  nPV_;
  Float_t PVx_;
  Float_t PVy_;
  Float_t PVz_;

  //pat jets
  std::vector<float> pfjetPt_;
  std::vector<float> pfjetEta_;
  std::vector<float> pfjetPhi_;
  //hlt jets
  std::vector<float> hltjetPt_;
  std::vector<float> hltjetEta_;
  std::vector<float> hltjetPhi_;
  std::vector<std::string> hltjetHLTpath_;
  std::vector<std::string> hltjetHLTfilter_;
  //HLT info
  std::vector<std::string> triggerResults_;
  std::vector<float> hltObjPt_;
  std::vector<float> hltObjEta_;
  std::vector<float> hltObjPhi_;
  std::vector<int> hltObjType_;
  std::vector<std::string> hltObjFilter_;

  // pat muons
  std::vector<float> muonPx_;
  std::vector<float> muonPy_;
  std::vector<float> muonPz_;
  std::vector<float> muonPt_;
  std::vector<float> muonEta_;
  std::vector<float> muonPhi_;
  std::vector<float> muonCharge_;
  std::vector<float> muonR04SumChargedHadronPt_;
  std::vector<float> muonR04SumChargedParticlePt_;
  std::vector<float> muonR04SumNeutralHadronEt_;
  std::vector<float> muonR04SumPhotonEt_;
  std::vector<float> muonR04SumPUPt_;
  std::vector<bool> muonIsPF_;
  std::vector<bool> muonIsGlobal_;
  std::vector<bool> muonIsTracker_;
  std::vector<bool> muonIsICHEPMedium_;
  std::vector<float> muonDz_;
  std::vector<float> muonDxy_;
  //std::vector<float> muonNormChi2_;

};

