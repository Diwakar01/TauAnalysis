#include "FakeRateStudyNtupleProducer.h"
#include <iostream>

FakeRateStudyNtupleProducer::FakeRateStudyNtupleProducer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  PVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag"));
  MetCollectionToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCollectionTag"));
  MuonCollectionToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("MuonCollectionTag"));
  PFJetCollectionToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("PFJetCollectionTag"));
  TriggerObjectCollectionToken_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerObjectCollectionTag"));
  triggerPaths_ = iConfig.getParameter<std::vector<std::string> >("triggerPaths");
  triggerFilters_ = iConfig.getParameter<std::vector<std::string> >("triggerFilters");
  triggerFiltersToken_ = consumes<std::vector<std::string> >(iConfig.getParameter<edm::InputTag>("triggerFiltersTag"));
  hlt_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltprocess"));
  pileupToken_ = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupTag"));
}


FakeRateStudyNtupleProducer::~FakeRateStudyNtupleProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FakeRateStudyNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ///event info
  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.getLuminosityBlock().luminosityBlock();
  
  //---------- pu -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  npuIT = 0.; npuOOT = 0.; npuTrue = 0.;
  if (!iEvent.isRealData()) {
    iEvent.getByToken(pileupToken_,PupInfo);
    if(PupInfo.isValid()){
      for(std::vector<PileupSummaryInfo>::const_iterator pileup = PupInfo->begin(); pileup != PupInfo->end(); ++pileup){
	if (pileup->getBunchCrossing() == 0) {
	  npuIT = pileup->getPU_NumInteractions();
	}
	else { 
	  npuOOT += pileup->getPU_NumInteractions();
	}
	npuTrue = pileup->getTrueNumInteractions();
      }
    }
  }


  math::XYZPoint pv_position = math::XYZPoint(0.,0.,0.);
  edm::Handle<reco::VertexCollection> Vertex;
  iEvent.getByToken(PVToken_, Vertex);
  nPV_ = 0;
  if(Vertex.isValid()) {
    for(unsigned i = 0 ; i < Vertex->size(); i++) {
      nPV_++;
      if(i == 0) {
	PVx_ = (*Vertex)[i].x();
	PVy_ = (*Vertex)[i].y();
	PVz_ = (*Vertex)[i].z();
	pv_position = (*Vertex)[i].position();
      }
    }
  }

  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(hlt_,hltresults);
  const edm::TriggerNames& TrigNames_ = iEvent.triggerNames(*hltresults);
  const int ntrigs = hltresults->size();
  //Get hlt jets
  edm::Handle<pat::TriggerObjectStandAloneCollection> hltobjects;
  iEvent.getByToken(TriggerObjectCollectionToken_,hltobjects);
  //store trigger paths that are passed in the event
  triggerResults_.clear();
  hltObjPt_.clear(); hltObjEta_.clear(); hltObjPhi_.clear();
  hltObjType_.clear(); hltObjFilter_.clear();
  for(int itrig=0; itrig<ntrigs; itrig++)
    {
      if(!hltresults->wasrun(itrig) )continue;
      std::string trigName = TrigNames_.triggerName(itrig);
      if(triggerPaths_.size() > 0){
        for(size_t ip = 0; ip < triggerPaths_.size(); ip++){
          if(trigName.find(triggerPaths_[ip]) != std::string::npos){
            if(hltresults->accept(itrig)){
              triggerResults_.push_back(trigName);
	      
	      //get HLT candidates passing that trigger
	      for(pat::TriggerObjectStandAloneCollection::const_iterator iobj = hltobjects->begin(); iobj != hltobjects->end(); ++iobj) {
		if( (trigName.find("Jet") != std::string::npos && iobj->hasTriggerObjectType(trigger::TriggerJet) ) ||
		    (trigName.find("Mu") != std::string::npos && iobj->hasTriggerObjectType(trigger::TriggerMuon) )){
		  if(iobj->hasPathName(trigName) && iobj->hasFilterLabel(triggerFilters_[ip])){
		    hltObjPt_.push_back(iobj->p4().pt());
		    hltObjEta_.push_back(iobj->p4().eta());
		    hltObjPhi_.push_back(iobj->p4().phi());
		    if(iobj->hasTriggerObjectType(trigger::TriggerMuon))hltObjType_.push_back(trigger::TriggerMuon);
		    else if(iobj->hasTriggerObjectType(trigger::TriggerJet))hltObjType_.push_back(trigger::TriggerJet);
		    
		    //get the last filter name for this object
		    if(iobj->hasPathLastFilterAccepted()){
		    std::vector< std::string > filters = iobj->filterLabels();
		    if(filters.size() > 0)hltObjFilter_.push_back(filters[filters.size() - 1]);
		    }

		  }
		}
	      }

            }
          }
        }
      }
    }

  //Get offline Reco/PAT Jets
  edm::Handle<pat::JetCollection> offlinepfjets;
  iEvent.getByToken(PFJetCollectionToken_,offlinepfjets);
  pfjetPt_.clear(); pfjetEta_.clear(); pfjetPhi_.clear();
  hltjetHLTpath_.clear(); hltjetHLTfilter_.clear();
  for(pat::JetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) {  
    if (ipf->pt() > 20){
      pfjetPt_.push_back(ipf->pt());
      pfjetEta_.push_back(ipf->eta()); 
      pfjetPhi_.push_back(ipf->phi()); 
    }
  }
  /*
  //Get hlt jets
  edm::Handle<pat::TriggerObjectStandAloneCollection> hltobjects;
  iEvent.getByToken(TriggerObjectCollectionToken_,hltobjects);
  hltjetPt_.clear(); hltjetEta_.clear(); hltjetPhi_.clear();
  for(pat::TriggerObjectStandAloneCollection::const_iterator iobj = hltobjects->begin(); iobj != hltobjects->end(); ++iobj) {  
    if (iobj->p4().pt() > 40 && iobj->hasTriggerObjectType(trigger::TriggerJet)){
      hltjetPt_.push_back(iobj->p4().pt());
      hltjetEta_.push_back(iobj->p4().eta()); 
      hltjetPhi_.push_back(iobj->p4().phi());
    }
  }
  */

  edm::Handle<pat::METCollection> METs;
  iEvent.getByToken(MetCollectionToken_, METs);
  
  metPt_ = 0.; metPhi_ = -999.;
  metPx_ = 0.; metPy_ = 0.;
  if(METs.isValid() && METs->size() > 0){
    metPt_ = (*METs)[0].pt();
    metPhi_ = (*METs)[0].phi();
    metPx_ = (*METs)[0].px();
    metPy_ = (*METs)[0].py();
  }
  
  muonPx_.clear(); muonPy_.clear(); muonPz_.clear(); 
  muonPt_.clear(); muonEta_.clear(); muonPhi_.clear(); muonCharge_.clear();
  muonR04SumChargedHadronPt_.clear(); muonR04SumChargedParticlePt_.clear();
  muonR04SumNeutralHadronEt_.clear(); muonR04SumPhotonEt_.clear();
  muonR04SumPUPt_.clear(); muonIsPF_.clear();
  muonIsGlobal_.clear(); muonIsTracker_.clear(); muonIsICHEPMedium_.clear();
  muonDz_.clear(); muonDxy_.clear(); //muonNormChi2_.clear();
  
  edm::Handle<pat::MuonCollection> Muons;
  iEvent.getByToken(MuonCollectionToken_, Muons);
  if(Muons.isValid())
    {
      for(unsigned i = 0 ; i < Muons->size() ; i++){
	if ((*Muons)[i].pt() < 10.) continue; //require a minimim cut
	muonPx_.push_back((*Muons)[i].px());
	muonPy_.push_back((*Muons)[i].py());
	muonPz_.push_back((*Muons)[i].pz());
	muonPt_.push_back((*Muons)[i].pt());
	muonEta_.push_back((*Muons)[i].eta());
	muonPhi_.push_back((*Muons)[i].phi());
	muonCharge_.push_back((*Muons)[i].charge());
	muonR04SumChargedHadronPt_.push_back((*Muons)[i].pfIsolationR04().sumChargedHadronPt);
	muonR04SumChargedParticlePt_.push_back((*Muons)[i].pfIsolationR04().sumChargedParticlePt);
	muonR04SumNeutralHadronEt_.push_back((*Muons)[i].pfIsolationR04().sumNeutralHadronEt);
	muonR04SumPhotonEt_.push_back((*Muons)[i].pfIsolationR04().sumPhotonEt);
	muonR04SumPUPt_.push_back((*Muons)[i].pfIsolationR04().sumPUPt);
	muonIsPF_.push_back((*Muons)[i].isPFMuon());
	muonIsGlobal_.push_back((*Muons)[i].isGlobalMuon());
	muonIsTracker_.push_back((*Muons)[i].isTrackerMuon());
	muonIsICHEPMedium_.push_back(isICHEPMuon((*Muons)[i]));
	reco::TrackRef bestTrack  = (*Muons)[i].muonBestTrack();
	if (bestTrack.isNonnull()) {
	  muonDxy_.push_back(bestTrack->dxy(pv_position));
	  muonDz_.push_back(bestTrack->dz(pv_position));
	}
	else{
	  muonDxy_.push_back(-9999.);
	  muonDz_.push_back(-9999.);
	}
      }
    }

  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
FakeRateStudyNtupleProducer::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");

  tree_->Branch("run", &run_, "run/l");
  tree_->Branch("event", &event_, "event/l");
  tree_->Branch("lumi", &lumi_, "lumi/l");

  tree_->Branch("npuIT", &npuIT, "npuIT/F");
  tree_->Branch("npuOOT", &npuOOT, "npuOOT/F");
  tree_->Branch("npuTrue", &npuTrue, "npuTrue/F");

  tree_->Branch("nPV", &nPV_, "nPV/i"); 
  tree_->Branch("PVx", &PVx_, "PVx/F");
  tree_->Branch("PVy", &PVy_, "PVy/F");
  tree_->Branch("PVz", &PVz_, "PVz/F");

  tree_->Branch("pfjetPt", "std::vector<float>", &pfjetPt_);
  tree_->Branch("pfjetEta", "std::vector<float>", &pfjetEta_);
  tree_->Branch("pfjetPhi", "std::vector<float>", &pfjetPhi_);
  /*
  tree_->Branch("hltjetPt", "std::vector<float>", &hltjetPt_);
  tree_->Branch("hltjetEta", "std::vector<float>", &hltjetEta_);
  tree_->Branch("hltjetPhi", "std::vector<float>", &hltjetPhi_);
  tree_->Branch("hltjetHLTpath", "std::vector<std::string>", &hltjetHLTpath_);
  tree_->Branch("hltjetHLTfilter", "std::vector<std::string>", &hltjetHLTfilter_);
  */
  tree_->Branch("metPx", &metPx_, "metPx/F");
  tree_->Branch("metPy", &metPy_, "metPy/F");
  tree_->Branch("metPt", &metPt_, "metPt/F");
  tree_->Branch("metPhi", &metPhi_, "metPhi/F");
  
  tree_->Branch("muonPx", "std::vector<float>", &muonPx_);
  tree_->Branch("muonPy", "std::vector<float>",&muonPy_);
  tree_->Branch("muonPz", "std::vector<float>",&muonPz_);
  tree_->Branch("muonPt", "std::vector<float>",&muonPt_);
  tree_->Branch("muonEta", "std::vector<float>",&muonEta_);
  tree_->Branch("muonPhi", "std::vector<float>",&muonPhi_);
  tree_->Branch("muonCharge", "std::vector<float>",&muonCharge_);
  tree_->Branch("muonR04SumChargedHadronPt", "std::vector<float>",&muonR04SumChargedHadronPt_);
  tree_->Branch("muonR04SumChargedParticlePt", "std::vector<float>",&muonR04SumChargedParticlePt_);
  tree_->Branch("muonR04SumNeutralHadronEt", "std::vector<float>",&muonR04SumNeutralHadronEt_);
  tree_->Branch("muonR04SumPhotonEt", "std::vector<float>",&muonR04SumPhotonEt_);
  tree_->Branch("muonR04SumPUPt", "std::vector<float>",&muonR04SumPUPt_);
  tree_->Branch("muonIsPF", "std::vector<bool>", &muonIsPF_);
  tree_->Branch("muonIsGlobal", "std::vector<bool>", &muonIsGlobal_);
  tree_->Branch("muonIsTracker", "std::vector<bool>", &muonIsTracker_);
  tree_->Branch("muonIsICHEPMedium", "std::vector<bool>", &muonIsICHEPMedium_);
  tree_->Branch("muonDz", "std::vector<float>", &muonDz_);
  tree_->Branch("muonDxy", "std::vector<float>",&muonDxy_);
  //tree_->Branch("muonNormChi2", "std::vector<float>",&muonNormChi2_);

  tree_->Branch("triggerResults", "std::vector<std::string>", &triggerResults_);
  tree_->Branch("hltObjPt", "std::vector<float>", &hltObjPt_);
  tree_->Branch("hltObjEta", "std::vector<float>", &hltObjEta_);
  tree_->Branch("hltObjPhi", "std::vector<float>", &hltObjPhi_);
  tree_->Branch("hltObjType", "std::vector<int>", &hltObjType_);
  tree_->Branch("hltObjFilter", "std::vector<std::string>", &hltObjFilter_);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FakeRateStudyNtupleProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
FakeRateStudyNtupleProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FakeRateStudyNtupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FakeRateStudyNtupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FakeRateStudyNtupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FakeRateStudyNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool FakeRateStudyNtupleProducer::isICHEPMuon(const reco::Muon& recoMu){
  bool goodGlob = recoMu.isGlobalMuon() && 
    recoMu.globalTrack()->normalizedChi2() < 3 && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.49 && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}
 
//define this as a plug-in
DEFINE_FWK_MODULE(FakeRateStudyNtupleProducer);
