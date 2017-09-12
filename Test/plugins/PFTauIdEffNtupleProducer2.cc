#include "TauAnalysis/Test/plugins/PFTauIdEffNtupleProducer2.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

int PFTauIdEffNtupleProducer2::verbosity_ = 0;

PFTauIdEffNtupleProducer2::PFTauIdEffNtupleProducer2(const edm::ParameterSet& cfg) 
  : moduleLabel_(cfg.getParameter<std::string>("@module_label")),
    ntuple_(0)
{
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  srcGenJets_ = cfg.getParameter<edm::InputTag>("srcGenJets");
  srcRecVetoElectrons_ = cfg.getParameter<edm::InputTag>("srcRecVetoElectrons");
  srcRecTaus_ = cfg.getParameter<edm::InputTag>("srcRecTaus");
  srcRecJets_ = cfg.getParameter<edm::InputTag>("srcRecJets");
  srcVertices_ = cfg.getParameter<edm::InputTag>("srcVertices");
  srcBeamSpot_ = cfg.getParameter<edm::InputTag>("srcBeamSpot");
  srcTracks_ = cfg.getParameter<edm::InputTag>("srcTracks");

  srcWeights_ = cfg.getParameter<vInputTag>("srcWeights");

  tauIdDiscriminators_.push_back("decayModeFindingNewDMs");
  tauIdDiscriminators_.push_back("decayModeFindingOldDMs");
  tauIdDiscriminators_.push_back("decayModeFinding");
  //tauIdDiscriminators_.push_back("byVLooseCombinedIsolationDeltaBetaCorr");
  //tauIdDiscriminators_.push_back("byLooseCombinedIsolationDeltaBetaCorr");
  //tauIdDiscriminators_.push_back("byMediumCombinedIsolationDeltaBetaCorr");
  //tauIdDiscriminators_.push_back("byTightCombinedIsolationDeltaBetaCorr");
  tauIdDiscriminators_.push_back("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  tauIdDiscriminators_.push_back("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  tauIdDiscriminators_.push_back("byTightCombinedIsolationDeltaBetaCorr3Hits");
  tauIdDiscriminators_.push_back("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  tauIdDiscriminators_.push_back("byLoosePileupWeightedIsolation3Hits");
  tauIdDiscriminators_.push_back("byMediumPileupWeightedIsolation3Hits");
  tauIdDiscriminators_.push_back("byTightPileupWeightedIsolation3Hits");
  tauIdDiscriminators_.push_back("byPileupWeightedIsolationRaw3Hits");
  tauIdDiscriminators_.push_back("byPhotonPtSumOutsideSignalCone");
  tauIdDiscriminators_.push_back("chargedIsoPtSum");
  tauIdDiscriminators_.push_back("neutralIsoPtSum");
  tauIdDiscriminators_.push_back("puCorrPtSum");
  tauIdDiscriminators_.push_back("neutralIsoPtSumWeight");
  tauIdDiscriminators_.push_back("photonPtSumOutsideSignalCone");
  tauIdDiscriminators_.push_back("footprintCorrection");
  ////tauIdDiscriminators_.push_back("byIsolationMVA3oldDMwoLTraw");
  ////tauIdDiscriminators_.push_back("byVLooseIsolationMVA3oldDMwoLT");
  ////tauIdDiscriminators_.push_back("byLooseIsolationMVA3oldDMwoLT");
  ////tauIdDiscriminators_.push_back("byMediumIsolationMVA3oldDMwoLT");
  ////tauIdDiscriminators_.push_back("byTightIsolationMVA3oldDMwoLT");
  ////tauIdDiscriminators_.push_back("byVTightIsolationMVA3oldDMwoLT");
  ////tauIdDiscriminators_.push_back("byVVTightIsolationMVA3oldDMwoLT");
  //tauIdDiscriminators_.push_back("byIsolationMVA3oldDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationMVA3oldDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationMVA3oldDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationMVA3oldDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationMVA3oldDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationMVA3oldDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationMVA3oldDMwLT");
  ////tauIdDiscriminators_.push_back("byIsolationMVA3newDMwoLTraw");
  ////tauIdDiscriminators_.push_back("byVLooseIsolationMVA3newDMwoLT");
  ////tauIdDiscriminators_.push_back("byLooseIsolationMVA3newDMwoLT");
  ////tauIdDiscriminators_.push_back("byMediumIsolationMVA3newDMwoLT");
  ////tauIdDiscriminators_.push_back("byTightIsolationMVA3newDMwoLT");
  ////tauIdDiscriminators_.push_back("byVTightIsolationMVA3newDMwoLT");
  ////tauIdDiscriminators_.push_back("byVVTightIsolationMVA3newDMwoLT");
  //tauIdDiscriminators_.push_back("byIsolationMVA3newDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationMVA3newDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationMVA3newDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationMVA3newDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationMVA3newDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationMVA3newDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationMVA3newDMwLT");
  //tauIdDiscriminators_.push_back("byIsolationRun2MVA1isoDBoldDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationRun2MVA1isoDBoldDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationRun2MVA1isoDBoldDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationRun2MVA1isoDBoldDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationRun2MVA1isoDBoldDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationRun2MVA1isoDBoldDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationRun2MVA1isoDBoldDMwLT");
  //tauIdDiscriminators_.push_back("byIsolationRun2MVA1isoDBnewDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationRun2MVA1isoDBnewDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationRun2MVA1isoDBnewDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationRun2MVA1isoDBnewDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationRun2MVA1isoDBnewDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationRun2MVA1isoDBnewDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationRun2MVA1isoDBnewDMwLT");
  //tauIdDiscriminators_.push_back("byIsolationRun2MVA1isoPWoldDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationRun2MVA1isoPWoldDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationRun2MVA1isoPWoldDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationRun2MVA1isoPWoldDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationRun2MVA1isoPWoldDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationRun2MVA1isoPWoldDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationRun2MVA1isoPWoldDMwLT");
  //tauIdDiscriminators_.push_back("byIsolationRun2MVA1isoPWnewDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationRun2MVA1isoPWnewDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationRun2MVA1isoPWnewDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationRun2MVA1isoPWnewDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationRun2MVA1isoPWnewDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationRun2MVA1isoPWnewDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationRun2MVA1isoPWnewDMwLT");
  //tauIdDiscriminators_.push_back("chargedIsoPtSumR03");
  //tauIdDiscriminators_.push_back("neutralIsoPtSumR03");
  //tauIdDiscriminators_.push_back("puCorrPtSumR03");
  //tauIdDiscriminators_.push_back("neutralIsoPtSumWeightR03");
  //tauIdDiscriminators_.push_back("footprintCorrectionR03");
  //tauIdDiscriminators_.push_back("photonPtSumOutsideSignalConeR03");
  //tauIdDiscriminators_.push_back("byIsolationRun2MVA1isoDBR03oldDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationRun2MVA1isoDBR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationRun2MVA1isoDBR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationRun2MVA1isoDBR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationRun2MVA1isoDBR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationRun2MVA1isoDBR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationRun2MVA1isoDBR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byIsolationRun2MVA1isoPWR03oldDMwLTraw");
  //tauIdDiscriminators_.push_back("byVLooseIsolationRun2MVA1isoPWR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byLooseIsolationRun2MVA1isoPWR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byMediumIsolationRun2MVA1isoPWR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byTightIsolationRun2MVA1isoPWR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byVTightIsolationRun2MVA1isoPWR03oldDMwLT");
  //tauIdDiscriminators_.push_back("byVVTightIsolationRun2MVA1isoPWR03oldDMwLT");
  //tauIdDiscriminators_.push_back("againstElectronLoose");
  //tauIdDiscriminators_.push_back("againstElectronMedium");
  //tauIdDiscriminators_.push_back("againstElectronTight");
  tauIdDiscriminators_.push_back("againstElectronMVA5raw");
  tauIdDiscriminators_.push_back("againstElectronMVA5category");
  tauIdDiscriminators_.push_back("againstElectronVLooseMVA5");
  tauIdDiscriminators_.push_back("againstElectronLooseMVA5");
  tauIdDiscriminators_.push_back("againstElectronMediumMVA5");
  tauIdDiscriminators_.push_back("againstElectronTightMVA5");
  tauIdDiscriminators_.push_back("againstElectronVTightMVA5");
  //tauIdDiscriminators_.push_back("againstElectronMVA3raw");
  //tauIdDiscriminators_.push_back("againstElectronMVA3category");
  //tauIdDiscriminators_.push_back("againstElectronLooseMVA3");
  //tauIdDiscriminators_.push_back("againstElectronMediumMVA3");
  //tauIdDiscriminators_.push_back("againstElectronTightMVA3");
  //tauIdDiscriminators_.push_back("againstElectronVTightMVA3");
  //tauIdDiscriminators_.push_back("againstElectronDeadECAL");
  //tauIdDiscriminators_.push_back("againstMuonLoose");
  //tauIdDiscriminators_.push_back("againstMuonMedium");
  //tauIdDiscriminators_.push_back("againstMuonTight");
  //tauIdDiscriminators_.push_back("againstMuonLoose2");
  //tauIdDiscriminators_.push_back("againstMuonMedium2");
  //tauIdDiscriminators_.push_back("againstMuonTight2");
  tauIdDiscriminators_.push_back("againstMuonLoose3");
  tauIdDiscriminators_.push_back("againstMuonTight3");
  //tauIdDiscriminators_.push_back("againstMuonMVAraw");
  //tauIdDiscriminators_.push_back("againstMuonLooseMVA");
  //tauIdDiscriminators_.push_back("againstMuonMediumMVA");
  //tauIdDiscriminators_.push_back("againstMuonTightMVA");

  minGenVisPt_ = 15.;
  minRecJetPt_ = 20.;
}

PFTauIdEffNtupleProducer2::~PFTauIdEffNtupleProducer2()
{
// nothing to be done yet...
}

void PFTauIdEffNtupleProducer2::beginJob()
{
//--- create TTree
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("pfTauIdEffNtuple", "pfTauIdEffNtuple");

//--- add branches 
  addBranchL("run");
  addBranchL("ls");
  addBranchL("event");
  addBranchF("genEvtWeight");

  addBranch_EnPxPyPz("recTau");
  addBranch_EnPxPyPz("recTauAlternate");
  addBranchI("recTauDecayMode");
  addBranch_EnPxPyPz("recRefJet");
  addBranchF("recRefJetZ");
  addBranch_EnPxPyPz("leadPFCand");
  addBranch_EnPxPyPz("leadPFChargedHadrCand");
  addBranchF("Tau_Energy_ECAL");
  addBranchF("Tau_Energy_HCAL");
  addBranchF("Tau_HoP");
  addBranchF("Tau_etaImpact_ECAL");
  addBranchF("Tau_phiImpact_ECAL");
  addBranchF("Tau_leadChHadDxy");
  addBranchI("Tau_leadChHadPixelHits");
  addBranchF("Tau_leadChHadDzWrtPV");
  addBranchF("Tau_MatchedVtxSumTkPt2");
  addBranchF("Tau_MatchedVtxZ");
  addBranchF("Tau_MatchedVtxRho");
  addBranchI("Tau_MatchedVtxNdof");
  addBranchI("Tau_MatchedVtxIsValid");
  addBranchF("Tau_MatchedVtxNChi2");
  addBranchI("Tau_MatchedVtxTkSize");
  addBranchF("Tau_leadChHadDzVtxMatch");
  addBranchI("Tau_leadChHadTkAlgo");
  addBranchI("Tau_leadChHadTkType");
  addBranchF("Tau_dzWrtPV");
  addBranchF("Tau_z");
  addBranchIV("Tau_sigChHadPixelHits");
  addBranchIV("Tau_sigChHadTkAlgo");
  for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
    addBranchF(*tauIdDiscriminator);
  }
  addBranchF("Muon_CaloCompatibility");
  addBranchF("Muon_SegmentCompatibility");
  addBranchI("Muon_StationMask");
  addBranchI("Muon_NumMatches");
  addBranchI("Muon_NumHitsDT_station0");
  addBranchI("Muon_NumHitsDT_station1");
  addBranchI("Muon_NumHitsDT_station2");
  addBranchI("Muon_NumHitsDT_station3");
  addBranchI("Muon_NumHitsCSC_station0");
  addBranchI("Muon_NumHitsCSC_station1");
  addBranchI("Muon_NumHitsCSC_station2");
  addBranchI("Muon_NumHitsCSC_station3");
  addBranchI("Muon_NumHitsRPC_station0");
  addBranchI("Muon_NumHitsRPC_station1");
  addBranchI("Muon_NumHitsRPC_station2");
  addBranchI("Muon_NumHitsRPC_station3");
  addBranchF("ElectronVetoDeltaR");
  addBranchI("ElectronVetoMatch");
  addBranch_EnPxPyPz("genTau");
  addBranchI("genTauDecayMode");
  addBranchI("genTauMatch");
  addBranchF("genTauDeltaR");
  addBranch_EnPxPyPz("genTauLeadChHad");
  addBranchI("genTauLchMother");
  addBranchF("genTauLchPhPt");
  addBranchI("genTauNDaughters");
  addBranchI("genTauNDaPhotons");
  addBranchF("genTauDecayLength");
  addBranch_EnPxPyPz("genElectron");
  addBranchI("genElectronMatch");
  addBranchF("genElectronDeltaR");
  addBranch_EnPxPyPz("genMuon");
  addBranchI("genMuonMatch");
  addBranchF("genMuonDeltaR");
  addBranch_EnPxPyPz("genQuarkOrGluon");
  addBranchI("genQuarkOrGluonMatch");
  addBranchF("genQuarkOrGluonDeltaR");
  addBranchI("genQuarkOrGluonPdgId");
  addBranch_EnPxPyPz("genJet");
  addBranchI("genJetMatch");
  addBranchF("genJetDeltaR");
  addBranchF("genJetZ");
  addBranchI("genNIsoChHad");
  addBranchF("genPtIsoChHad");
  addBranchF("genPtIsoNeuHad");
  addBranchF("genPtIsoPhoton");
  addBranch_EnPxPyPz("recJet");
  addBranch_EnPxPyPz("recJetRaw");
  addBranchI("recJetMatch");
  addBranchF("recJetDeltaR");
  addBranchF("recJetZ");
  addBranchI("ntrkiso");
  addBranchI("ntrkiso_pt1p0");
  addBranchI("ntrkiso_pt1p5");
  addBranchF("pttrkiso");
  addBranchF("pttrkiso_pt1p0");
  addBranchF("pttrkiso_pt1p5");
  addBranchI("ntrkiso_nopf");
  addBranchF("pttrkiso_nopf");
  addBranchI("ntrkiso_pfSigGamma");
  addBranchF("pttrkiso_pfSigGamma");
  addBranchI("ntrkiso_pfIsoGamma");
  addBranchF("pttrkiso_pfIsoGamma");
  addBranchIV("isotrkalgo");
  addBranchIV("isotrkpixelhits");
  addBranchIV("isotrktrackhits");
  addBranchFV("isotrkpt");
  addBranchFV("isotrketa");
  addBranchFV("isotrkphi");
  addBranchI("numVertices");
  addBranchF("vertexZ");
  addBranchF("evtWeight");
}

const pat::Tau* PFTauIdEffNtupleProducer2::findMatchingRecTau(const pat::TauCollection& recTaus, const reco::Candidate::LorentzVector& genParticleP4)
{
  const pat::Tau* recTau_matched = 0;
  
  double genTauDeltaR = 9.9;
  for ( pat::TauCollection::const_iterator recTau = recTaus.begin();
	recTau != recTaus.end(); ++recTau ) {
    
    if ( recTau->pt() < 15. ) continue;
    
    double dR = deltaR(recTau->p4(), genParticleP4);
    if ( dR < 0.3 && dR < genTauDeltaR ) {
      genTauDeltaR = dR;
      recTau_matched = &(*recTau);
    }
  }

  return recTau_matched;
}

const pat::Electron* PFTauIdEffNtupleProducer2::findMatchingElectronVeto(const pat::Tau& recTau, const pat::ElectronCollection& recVetoElectrons)
{
  const pat::Electron* recElectronVeto_matched = 0;
  
  double EleVetoDeltaR = 9.9;
  for ( pat::ElectronCollection::const_iterator recElectron = recVetoElectrons.begin();
	recElectron != recVetoElectrons.end(); ++recElectron ) {
    
    if ( recElectron->pt() < 10. ) continue;
    
    double dR = deltaR(recElectron->p4(), recTau);
    if ( dR < 0.3 && dR < EleVetoDeltaR ) {
      EleVetoDeltaR = dR;
      recElectronVeto_matched = &(*recElectron);
    }
  }

  return recElectronVeto_matched;
}

namespace
{
  void compECALeta_and_phiImpact(const pat::Tau* recTau_matched, double& etaImpact, double& phiImpact)
  {
    etaImpact = -99.;
    phiImpact = -99.;

    if ( recTau_matched ) {
      double sumEtaTimesEnergy = 0.;
      double sumPhiTimesEnergy = 0.;
      double sumEnergy         = 0.;
      std::vector<reco::PFCandidatePtr> signalPFCands = recTau_matched->signalPFCands();
      for ( std::vector<reco::PFCandidatePtr>::const_iterator signalPFCand = signalPFCands.begin();
	    signalPFCand != signalPFCands.end(); ++signalPFCand ) {
	double energy = (*signalPFCand)->energy();
	sumEtaTimesEnergy += (energy*(*signalPFCand)->positionAtECALEntrance().eta());
	sumPhiTimesEnergy += (energy*(*signalPFCand)->positionAtECALEntrance().phi());
	sumEnergy         += energy;
      }
      if ( sumEnergy > 0. ) {
	etaImpact = sumEtaTimesEnergy/sumEnergy;
	phiImpact = sumPhiTimesEnergy/sumEnergy;
      }
    }
  }
}


double compSumTrackPt2(const reco::Vertex& vertex, double maxTrackWeight)
{
  double sumTrackPt2 = 0.;
  for ( reco::Vertex::trackRef_iterator track = vertex.tracks_begin();
	track != vertex.tracks_end(); ++track ) {
    if ( !(track->isNonnull() && track->isAvailable()) ) continue;
    double trackPt = (*track)->pt();
    if ( maxTrackWeight > 0. && trackPt > maxTrackWeight ) trackPt = maxTrackWeight;
    sumTrackPt2 += (trackPt*trackPt);
  }
  return sumTrackPt2;
}

double compGenTauDecayDistance(const reco::GenParticle* genTau)
{
  double genTauDecayDistance = -1.;
  std::vector<const reco::GenParticle*> genDaughters;
  findDaughters(genTau, genDaughters, 1);
  int numDaughters = genDaughters.size();
  for ( int iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    if ( TMath::Abs(genDaughters.at(iDaughter)->pdgId()) == 16 )
      genTauDecayDistance = TMath::Sqrt((genTau->vertex() - genDaughters.at(iDaughter)->vertex()).Mag2());
  }
  return genTauDecayDistance;
}

void PFTauIdEffNtupleProducer2::setRecTauValues(const pat::Tau* recTau_matched, const pat::ElectronCollection& recVetoElectrons, const reco::VertexCollection& vertices, const reco::TrackCollection& recoTracks, const reco::BeamSpot& beamSpot)
{
  //reco::Vertex* vtx = 0;
  //if(vertices.size() > 0) vtx = (&(vertices[0]));
  if ( recTau_matched ) {
    //std::cout << "rec Tau Vis Pt "<<recTau_matched->pt()<<" eta "<<recTau_matched->eta()<<" phi "<<recTau_matched->phi()<<std::endl;
    setValue_EnPxPyPz("recTau", recTau_matched->p4());
    setValue_EnPxPyPz("recTauAlternate", recTau_matched->alternatLorentzVect());
    setValueI("recTauDecayMode", recTau_matched->decayMode());
    setValue_EnPxPyPz("recRefJet", recTau_matched->pfJetRef()->p4());
    //compute the Vertex.Z of refJet
    double recRefJetZ = 0, recRefJetE = 0;
    reco::TrackRefVector recRefJetTracks = recTau_matched->pfJetRef()->getTrackRefs();
    if(recRefJetTracks.size() <= 0)recRefJetZ = -999.;
    else{
      for(size_t it = 0; it < recRefJetTracks.size(); it++){
	recRefJetZ += (recRefJetTracks)[it]->vz() * (recRefJetTracks)[it]->p();
	recRefJetE += (recRefJetTracks)[it]->p();
      }
    }
    if(recRefJetE != 0) recRefJetZ = recRefJetZ/recRefJetE;
    setValueF("recRefJetZ", recRefJetZ);
    if ( recTau_matched->leadPFCand().isNonnull() ) setValue_EnPxPyPz("leadPFCand", recTau_matched->leadPFCand()->p4());
    else setValue_EnPxPyPz("leadPFCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    int Tau_leadChHadTkAlgo_ = -99, Tau_leadChHadTkType_ = -99;
    if ( recTau_matched->leadPFChargedHadrCand().isNonnull() ){
      setValue_EnPxPyPz("leadPFChargedHadrCand", recTau_matched->leadPFChargedHadrCand()->p4());
      setValueF("Tau_Energy_ECAL", recTau_matched->leadPFChargedHadrCand()->ecalEnergy());
      setValueF("Tau_Energy_HCAL", recTau_matched->leadPFChargedHadrCand()->hcalEnergy());
      setValueF("Tau_HoP", recTau_matched->leadPFChargedHadrCand()->hcalEnergy()/recTau_matched->leadPFChargedHadrCand()->p());
      if(recTau_matched->leadPFChargedHadrCand()->trackRef().isNonnull() ){
	setValueF("Tau_leadChHadDxy", recTau_matched->leadPFChargedHadrCand()->trackRef()->dxy(beamSpot.position()));
	setValueI("Tau_leadChHadPixelHits", recTau_matched->leadPFChargedHadrCand()->trackRef()->hitPattern().numberOfValidPixelHits());
	if(vertices.size() > 0 && vertices[0].isValid())setValueF("Tau_leadChHadDzWrtPV", recTau_matched->leadPFChargedHadrCand()->trackRef()->dz(vertices[0].position()));
	else setValueF("Tau_leadChHadDzWrtPV", -999.);
	Tau_leadChHadTkAlgo_ = recTau_matched->leadPFChargedHadrCand()->trackRef()->algo();
	Tau_leadChHadTkType_ = 1;
      }
      else if(recTau_matched->leadPFChargedHadrCand()->gsfTrackRef().isNonnull() ){
        setValueF("Tau_leadChHadDxy", recTau_matched->leadPFChargedHadrCand()->gsfTrackRef()->dxy(beamSpot.position()));
        setValueI("Tau_leadChHadPixelHits", recTau_matched->leadPFChargedHadrCand()->gsfTrackRef()->hitPattern().numberOfValidPixelHits());
        if(vertices.size() > 0 && vertices[0].isValid())setValueF("Tau_leadChHadDzWrtPV", recTau_matched->leadPFChargedHadrCand()->gsfTrackRef()->dz(vertices[0].position()));
        else setValueF("Tau_leadChHadDzWrtPV", -999.);
	Tau_leadChHadTkAlgo_ = recTau_matched->leadPFChargedHadrCand()->gsfTrackRef()->algo();
	Tau_leadChHadTkType_ = 2;
      }
    }
    else setValue_EnPxPyPz("leadPFChargedHadrCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValueI("Tau_leadChHadTkAlgo", Tau_leadChHadTkAlgo_);
    setValueI("Tau_leadChHadTkType", Tau_leadChHadTkType_);
    double etaImpact_ECAL, phiImpact_ECAL;
    compECALeta_and_phiImpact(recTau_matched, etaImpact_ECAL, phiImpact_ECAL);
    setValueF("Tau_etaImpact_ECAL", etaImpact_ECAL);
    setValueF("Tau_phiImpact_ECAL", phiImpact_ECAL);
    double dzWrtPV = (vertices.size() > 0 && vertices[0].isValid()) ? fabs( recTau_matched->vertex().z() - vertices[0].position().z() ) : -99;
    setValueF("Tau_dzWrtPV", dzWrtPV);
    setValueF("Tau_z", recTau_matched->vertex().z());
    double sumTrackPt2 = -999., vtxZ = -999., vtxRho = -999.; 
    int vtxIsValid = -999, vtxNdof = 0, vtxTrackSize = 0;
    double vtxNChi2 = 999.;
    double dZ_vtxMatch = 1.e+3;
    
    for(reco::VertexCollection::const_iterator iv = vertices.begin(); iv != vertices.end(); iv++){
      if(iv->z() == recTau_matched->vertex().z()){
	sumTrackPt2 = compSumTrackPt2((*iv), 0.);
	vtxZ = iv->z();
	vtxRho = iv->position().Rho();
	vtxIsValid = int(iv->isValid());
	vtxNdof = iv->ndof();
	vtxNChi2 = iv->normalizedChi2();
	vtxTrackSize = iv->tracksSize();
      }
      //std::cout<<" vertex (x, y, z) "<<iv->x()<<" "<<iv->y()<<" "<<iv->z()<<" trk size "<<iv->tracksSize()<<std::endl; 
      for(reco::Vertex::trackRef_iterator it = iv->tracks_begin(); it != iv->tracks_end(); it++){
	//const reco::TrackBaseRef track = reco::TrackBaseRef(*it);
	//if(track == reco::TrackBaseRef(recTau_matched->leadPFChargedHadrCand()->trackRef()) ||
	//track == reco::TrackBaseRef(recTau_matched->leadPFChargedHadrCand()->gsfTrackRef())){
	if(it->isNonnull() && it->isAvailable() && recTau_matched->leadPFChargedHadrCand().isNonnull()){
	  //std::cout<<" Has lead ch. hadron "<<std::endl;
	  if((recTau_matched->leadPFChargedHadrCand()->trackRef().isNonnull() && 
	     recTau_matched->leadPFChargedHadrCand()->trackRef().isAvailable() &&
	      it->get() == recTau_matched->leadPFChargedHadrCand()->trackRef().get()) ||
	     (recTau_matched->leadPFChargedHadrCand()->gsfTrackRef().isNonnull() &&
	      recTau_matched->leadPFChargedHadrCand()->gsfTrackRef().isAvailable() &&
	      it->get() == recTau_matched->leadPFChargedHadrCand()->gsfTrackRef().get())){
	    double tkdZ = (*it)->dz((iv->position()));
	    //std::cout<<" matches to KF or Gsf track with dZ = "<<tkdZ<<std::endl;
	    if(fabs(tkdZ) < fabs(dZ_vtxMatch)){
	      dZ_vtxMatch = tkdZ;
	    }
	  }
	  /*
	  if((recTau_matched->leadPFChargedHadrCand()->trackRef().isNonnull() &&
	      recTau_matched->leadPFChargedHadrCand()->trackRef().isAvailable() &&
	      deltaR((*it)->momentum(), recTau_matched->leadPFChargedHadrCand()->trackRef()->momentum()) < 0.01) ||
	     (recTau_matched->leadPFChargedHadrCand()->gsfTrackRef().isNonnull() &&
              recTau_matched->leadPFChargedHadrCand()->gsfTrackRef().isAvailable() &&
	      deltaR((*it)->momentum(), recTau_matched->leadPFChargedHadrCand()->gsfTrackRef()->momentum())< 0.01)){
	    double tkdZ = (*it)->dz((iv->position()));
	    std::cout<<" matches to KF or Gsf track using deltaR with dZ = "<<tkdZ<<std::endl;
	  }
	  */
	}
      }
    }
    setValueF("Tau_MatchedVtxSumTkPt2", sumTrackPt2);
    setValueF("Tau_MatchedVtxZ", vtxZ);
    setValueF("Tau_MatchedVtxRho", vtxRho);
    setValueI("Tau_MatchedVtxIsValid", vtxIsValid);
    setValueI("Tau_MatchedVtxNdof", vtxNdof);
    setValueF("Tau_MatchedVtxNChi2", vtxNChi2);
    setValueI("Tau_MatchedVtxTkSize", vtxTrackSize);
    setValueF("Tau_leadChHadDzVtxMatch", dZ_vtxMatch);
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
      setValueF(*tauIdDiscriminator, recTau_matched->tauID(*tauIdDiscriminator));
    }
    
    double muon_CaloCompatibility    = 0.;
    double muon_SegmentCompatibility = 0.;
    int    muon_StationMask          = 0;
    int    muon_NumMatches           = 0;
    int    muon_NumHitsDT_station0   = 0;
    int    muon_NumHitsDT_station1   = 0;
    int    muon_NumHitsDT_station2   = 0;
    int    muon_NumHitsDT_station3   = 0;
    int    muon_NumHitsCSC_station0  = 0;
    int    muon_NumHitsCSC_station1  = 0;
    int    muon_NumHitsCSC_station2  = 0;
    int    muon_NumHitsCSC_station3  = 0;
    int    muon_NumHitsRPC_station0  = 0;
    int    muon_NumHitsRPC_station1  = 0;
    int    muon_NumHitsRPC_station2  = 0;
    int    muon_NumHitsRPC_station3  = 0;
    if ( recTau_matched->leadPFChargedHadrCand().isNonnull() ) {
      reco::MuonRef muonRef = recTau_matched->leadPFChargedHadrCand()->muonRef();      
      if ( muonRef.isNonnull() ) {
	muon_CaloCompatibility = muonRef->caloCompatibility();
	muon_SegmentCompatibility = muon::segmentCompatibility(*muonRef);
	muon_StationMask = muonRef->stationMask(reco::Muon::NoArbitration);
	muon_NumMatches = muonRef->numberOfMatches(reco::Muon::NoArbitration);
	if ( muonRef->outerTrack().isNonnull() ) {
	  const reco::HitPattern& muonHitPattern = muonRef->outerTrack()->hitPattern();
	  for ( int iHit = 0; iHit < muonHitPattern.numberOfValidHits(); ++iHit ) {
	    uint32_t hit = muonHitPattern.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, iHit);
	    if ( hit == 0 ) break;	    
	    if ( muonHitPattern.muonHitFilter(hit) && muonHitPattern.getHitType(hit) == TrackingRecHit::valid ) {
	      int muonStation = muonHitPattern.getMuonStation(hit);
	      if ( muonHitPattern.muonDTHitFilter(hit) ) {
		if      ( muonStation == 1 ) ++muon_NumHitsDT_station0;
		else if ( muonStation == 2 ) ++muon_NumHitsDT_station1;
		else if ( muonStation == 3 ) ++muon_NumHitsDT_station2;
		else if ( muonStation == 4 ) ++muon_NumHitsDT_station3;
	      }
	      if ( muonHitPattern.muonCSCHitFilter(hit) ) {
		if      ( muonStation == 1 ) ++muon_NumHitsCSC_station0;
		else if ( muonStation == 2 ) ++muon_NumHitsCSC_station1;
		else if ( muonStation == 3 ) ++muon_NumHitsCSC_station2;
		else if ( muonStation == 4 ) ++muon_NumHitsCSC_station3;
	      }
	      if ( muonHitPattern.muonRPCHitFilter(hit) ) {
		if      ( muonStation == 1 ) ++muon_NumHitsRPC_station0;
		else if ( muonStation == 2 ) ++muon_NumHitsRPC_station1;
		else if ( muonStation == 3 ) ++muon_NumHitsRPC_station2;
		else if ( muonStation == 4 ) ++muon_NumHitsRPC_station3;
	      }
	    }
	  }
	}
      }
    }
    setValueF("Muon_CaloCompatibility", muon_CaloCompatibility);
    setValueF("Muon_SegmentCompatibility", muon_SegmentCompatibility);
    setValueI("Muon_StationMask", muon_StationMask);
    setValueI("Muon_NumMatches", muon_NumMatches);
    setValueI("Muon_NumHitsDT_station0", muon_NumHitsDT_station0);
    setValueI("Muon_NumHitsDT_station1", muon_NumHitsDT_station1);
    setValueI("Muon_NumHitsDT_station2", muon_NumHitsDT_station2);
    setValueI("Muon_NumHitsDT_station3", muon_NumHitsDT_station3);
    setValueI("Muon_NumHitsCSC_station0", muon_NumHitsCSC_station0);
    setValueI("Muon_NumHitsCSC_station1", muon_NumHitsCSC_station1);
    setValueI("Muon_NumHitsCSC_station2", muon_NumHitsCSC_station2);
    setValueI("Muon_NumHitsCSC_station3", muon_NumHitsCSC_station3);
    setValueI("Muon_NumHitsRPC_station0", muon_NumHitsRPC_station0);
    setValueI("Muon_NumHitsRPC_station1", muon_NumHitsRPC_station1);
    setValueI("Muon_NumHitsRPC_station2", muon_NumHitsRPC_station2);
    setValueI("Muon_NumHitsRPC_station3", muon_NumHitsRPC_station3);

    const pat::Electron* ElectronVeto_matched = findMatchingElectronVeto(*recTau_matched, recVetoElectrons);
    double ElectronVetoDeltaR = ( ElectronVeto_matched ) ? deltaR(ElectronVeto_matched->p4(), recTau_matched->p4()) : 9.9;
    bool ElectronVetoMatch = (ElectronVetoDeltaR < 0.3);
    setValueF("ElectronVetoDeltaR", ElectronVetoDeltaR);
    setValueI("ElectronVetoMatch", ElectronVetoMatch);

    //Get signal tracks from Taus                                                                                                                                                
    int ntracks = 0, ntracks_pt1p0 = 0, ntracks_pt1p5 = 0;
    double pttracks = 0, pttracks_pt1p0 = 0, pttracks_pt1p5 = 0;                                                                                                          \
    int ntracks_nopf = 0, ntracks_pfSigGamma = 0, ntracks_pfIsoGamma = 0;
    double pttracks_nopf = 0, pttracks_pfSigGamma = 0, pttracks_pfIsoGamma = 0;
    std::vector<Int_t> isotrkalgo_; isotrkalgo_.clear();
    std::vector<Int_t> isotrkpixelhits_; isotrkpixelhits_.clear();
    std::vector<Int_t> isotrktrackhits_; isotrktrackhits_.clear();
    std::vector<Float_t> isotrkpt_; isotrkpt_.clear();
    std::vector<Float_t> isotrketa_; isotrketa_.clear();
    std::vector<Float_t> isotrkphi_; isotrkphi_.clear();

    std::vector<reco::TrackBaseRef> SignalTracks;
    SignalTracks.clear();
    std::vector<Int_t> sigtrkalgo_; sigtrkalgo_.clear();
    std::vector<Int_t> sigtrkpixelhits_; sigtrkpixelhits_.clear();
    if(recTau_matched->signalPFChargedHadrCands().size() > 0){
      std::vector<reco::PFCandidatePtr> cands = recTau_matched->signalPFChargedHadrCands();
      for (std::vector<reco::PFCandidatePtr>::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
	if(!iter->isNonnull())continue;
	if(iter->get()->trackRef().isNonnull()){
	  SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
	  sigtrkalgo_.push_back((*iter)->trackRef()->algo());
	  sigtrkpixelhits_.push_back((*iter)->trackRef()->hitPattern().numberOfValidPixelHits());
	}
	else if(iter->get()->gsfTrackRef().isNonnull()){
	  SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));
	  sigtrkalgo_.push_back((*iter)->gsfTrackRef()->algo());
          sigtrkpixelhits_.push_back((*iter)->gsfTrackRef()->hitPattern().numberOfValidPixelHits());
	}
      }
    }
    setValueIV("Tau_sigChHadPixelHits", sigtrkpixelhits_);
    setValueIV("Tau_sigChHadTkAlgo", sigtrkalgo_);

    std::vector<reco::TrackBaseRef> allTracks;
    allTracks.clear();
    if(recTau_matched->signalPFCands().size() > 0){
      std::vector<reco::PFCandidatePtr> cands = recTau_matched->signalPFCands();
      for (std::vector<reco::PFCandidatePtr>::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
        if(!iter->isNonnull())continue;
        if(iter->get()->trackRef().isNonnull()) allTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
        else if(iter->get()->gsfTrackRef().isNonnull()){allTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
      }
    }
    if(recTau_matched->isolationPFCands().size() > 0){
      std::vector<reco::PFCandidatePtr> cands = recTau_matched->isolationPFCands();
      for (std::vector<reco::PFCandidatePtr>::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
        if(!iter->isNonnull())continue;
        if(iter->get()->trackRef().isNonnull()) allTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
        else if(iter->get()->gsfTrackRef().isNonnull()){allTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
      }
    }

    std::vector<reco::TrackBaseRef> sigGammaTracks;
    sigGammaTracks.clear();
    if(recTau_matched->signalPFGammaCands().size() > 0){
      std::vector<reco::PFCandidatePtr> cands = recTau_matched->signalPFGammaCands();
      for (std::vector<reco::PFCandidatePtr>::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
        if(!iter->isNonnull())continue;
        if(iter->get()->trackRef().isNonnull()) sigGammaTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
        else if(iter->get()->gsfTrackRef().isNonnull()){sigGammaTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
      }
    }

    std::vector<reco::TrackBaseRef> isoGammaTracks;
    isoGammaTracks.clear();
    if(recTau_matched->isolationPFGammaCands().size() > 0){
      std::vector<reco::PFCandidatePtr> cands = recTau_matched->isolationPFGammaCands();
      for (std::vector<reco::PFCandidatePtr>::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
        if(!iter->isNonnull())continue;
        if(iter->get()->trackRef().isNonnull()) isoGammaTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
        else if(iter->get()->gsfTrackRef().isNonnull()){isoGammaTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
      }
    }

    //std::cout<<"size  of tracks "<<allTracks.size()<<" "<<sigGammaTracks.size()<<" "<<isoGammaTracks.size()<<" "<<SignalTracks.size()<<std::endl; 
    unsigned int idx = 0;
    for (reco::TrackCollection::const_iterator trk = recoTracks.begin(); trk != recoTracks.end(); ++trk, idx++) {
      if(deltaR(recTau_matched->p4(), trk->momentum()) < 0.5){

	// apply quality cuts
	bool track_quality = false;
	double tk_dz = TMath::Abs(trk->dz(recTau_matched->vertex()));
	double tk_dxy = TMath::Abs(trk->dxy(recTau_matched->vertex()));
	int tk_hits = trk->hitPattern().numberOfValidHits();
	int tk_nChi2 = trk->normalizedChi2();
	
	if ( tk_dz < 0.2 && tk_dxy < 0.03 && tk_hits >= 3 && tk_nChi2 < 100 && trk->pt() > 0.5) track_quality = true;

	if(!track_quality) continue;
	
	//check it is not signal track
	reco::TrackRef tmpRef(&recoTracks, idx);
	reco::TrackRef tmpRefForBase=tmpRef;
	bool isSigTrk = false;
	for (unsigned int sigTrk = 0; sigTrk < SignalTracks.size(); sigTrk++) {
	  //if (reco::TrackBaseRef(tmpRefForBase)==SignalTracks.at(sigTrk)){ 
	  if(SignalTracks[sigTrk]->pt() == trk->pt() && SignalTracks[sigTrk]->eta() == trk->eta() && SignalTracks[sigTrk]->phi() == trk->phi()){
	    isSigTrk = true;
	    //std::cout<<" is a signal track "<<std::endl;
	  }
	}
	if(!isSigTrk){
	  if(trk->pt() > 0.5){
	    ntracks++;
	    pttracks += trk->pt();
	  }
	  if(trk->pt() > 1.0){
	    ntracks_pt1p0++;
	    pttracks_pt1p0 += trk->pt();
	  }
	  if(trk->pt() > 1.5){
	    ntracks_pt1p5++;
	    pttracks_pt1p5 += trk->pt();
	  }
	  
	  //check if no PF cand is found
	  bool isPFCand = false;
	  for (unsigned int allTrk = 0; allTrk < allTracks.size(); allTrk++) {
	    //if (reco::TrackBaseRef(tmpRefForBase)==allTracks.at(allTrk)){
	    if(allTracks[allTrk]->pt() == trk->pt() && allTracks[allTrk]->eta() == trk->eta() && allTracks[allTrk]->phi() == trk->phi()){
	      isPFCand = true;
	      //std::cout<<" is a pf cand "<<std::endl;
	    }
	  }
	  if(!isPFCand){
	    ntracks_nopf++;
	    pttracks_nopf += trk->pt();
	  }

	  //check if PF gamma cand is created
	  //in signal cone
	  bool isSigPFGamma = false;
	  for (unsigned int gammaTrk = 0; gammaTrk < sigGammaTracks.size(); gammaTrk++) {
            //if (reco::TrackBaseRef(tmpRefForBase)==sigGammaTracks.at(gammaTrk)){
	    if(sigGammaTracks[gammaTrk]->pt() == trk->pt() && sigGammaTracks[gammaTrk]->eta() == trk->eta() && sigGammaTracks[gammaTrk]->phi() == trk->phi()){
	      isSigPFGamma = true;
	      //std::cout<<" is a sig ph cand "<<std::endl;
	    }
          }
	  if(isSigPFGamma){
            ntracks_pfSigGamma++;
            pttracks_pfSigGamma += trk->pt();
          }

	  //in isolation cone
	  bool isIsoPFGamma = false;
          for (unsigned int gammaTrk = 0; gammaTrk < isoGammaTracks.size(); gammaTrk++) {
            //if (reco::TrackBaseRef(tmpRefForBase)==isoGammaTracks.at(gammaTrk)){
	    if(isoGammaTracks[gammaTrk]->pt() == trk->pt() && isoGammaTracks[gammaTrk]->eta() == trk->eta() && isoGammaTracks[gammaTrk]->phi() == trk->phi()){
	      isIsoPFGamma = true;
	      //std::cout<<" is a iso ph cand "<<std::endl;
	    }
          }
          if(isIsoPFGamma){
            ntracks_pfIsoGamma++;
            pttracks_pfIsoGamma += trk->pt();
          }
	 
	  //Fill individual tracks
	  isotrkalgo_.push_back(trk->algo());
	  isotrkpixelhits_.push_back(trk->hitPattern().numberOfValidPixelHits());
	  isotrktrackhits_.push_back(trk->hitPattern().numberOfValidTrackerHits());
	  isotrkpt_.push_back(trk->pt());
	  isotrketa_.push_back(trk->eta());
	  isotrkphi_.push_back(trk->phi());

	}
      }
    }
    setValueI("ntrkiso", ntracks);
    setValueI("ntrkiso_pt1p0", ntracks_pt1p0);
    setValueI("ntrkiso_pt1p5", ntracks_pt1p5);
    setValueF("pttrkiso", pttracks);
    setValueF("pttrkiso_pt1p0", pttracks_pt1p0);
    setValueF("pttrkiso_pt1p5", pttracks_pt1p5);
    setValueI("ntrkiso_nopf", ntracks_nopf);
    setValueF("pttrkiso_nopf", pttracks_nopf);
    setValueI("ntrkiso_pfSigGamma", ntracks_pfSigGamma);
    setValueF("pttrkiso_pfSigGamma", pttracks_pfSigGamma);
    setValueI("ntrkiso_pfIsoGamma", ntracks_pfIsoGamma);
    setValueF("pttrkiso_pfIsoGamma", pttracks_pfIsoGamma);
    setValueIV("isotrkalgo", isotrkalgo_);
    setValueIV("isotrkpixelhits", isotrkpixelhits_);
    setValueIV("isotrktrackhits", isotrktrackhits_);
    setValueFV("isotrkpt", isotrkpt_);
    setValueFV("isotrketa", isotrketa_);
    setValueFV("isotrkphi", isotrkphi_);
  } else {
    setValue_EnPxPyPz("recTau", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValueI("recTauDecayMode", -1);
    setValue_EnPxPyPz("recRefJet", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValueF("recRefJetZ", -999.);
    setValue_EnPxPyPz("leadPFCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValue_EnPxPyPz("leadPFChargedHadrCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValueF("Tau_Energy_ECAL", 0.);
    setValueF("Tau_Energy_HCAL", 0.);
    setValueF("Tau_HoP", 0.);
    setValueF("Tau_etaImpact_ECAL", 0.);
    setValueF("Tau_phiImpact_ECAL", 0.);    
    setValueF("Tau_leadChHadDxy", -999.);
    setValueF("Tau_leadChHadPixelHits", 0.);
    setValueF("Tau_leadChHadDzWrtPV", -999.);
    setValueF("Tau_MatchedVtxSumTkPt2", -999.);
    setValueF("Tau_MatchedVtxZ", -999.);
    setValueF("Tau_MatchedVtxRho", -999.);
    setValueI("Tau_MatchedVtxIsValid", -999);
    setValueI("Tau_MatchedVtxNdof", 0);
    setValueF("Tau_MatchedVtxNChi2", 999.);
    setValueI("Tau_MatchedVtxTkSize", 0);
    setValueF("Tau_leadChHadDzVtxMatch", 1.e+3);
    setValueF("Tau_dzWrtPV", 99.);
    setValueF("Tau_z", -999.);
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
      setValueF(*tauIdDiscriminator, 0.);
    }
    setValueF("ElectronVetoDeltaR", 0.);
    setValueI("ElectronVetoMatch", 0.);
    setValueI("ntrkiso", -99);
    setValueI("ntrkiso_pt1p0", -99);
    setValueI("ntrkiso_pt1p5", -99);
    setValueF("pttrkiso", -99.);
    setValueF("pttrkiso_pt1p0", -99.);
    setValueF("pttrkiso_pt1p5", -99.);
    setValueI("ntrkiso_nopf", -99);
    setValueF("pttrkiso_nopf", -99.);
    setValueI("ntrkiso_pfSigGamma", -99);
    setValueF("pttrkiso_pfSigGamma", -99.);
    setValueI("ntrkiso_pfIsoGamma", -99);
    setValueF("pttrkiso_pfIsoGamma", -99.);
  }
}
/*
void PFTauIdEffNtupleProducer2::setRecTauExtraValues(int ntrkiso_pt0p5, int ntrkiso_pt1p0, int ntrkiso_pt1p5, double pttrkiso_pt0p5, double pttrkiso_pt1p0, double pttrkiso_pt1p5){
  setValueI("ntrkiso_pt0p5", ntrkiso_pt0p5);
  setValueI("ntrkiso_pt1p0", ntrkiso_pt1p0);
  setValueI("ntrkiso_pt1p5", ntrkiso_pt1p5);
  setValueF("pttrkiso_pt0p5", pttrkiso_pt0p5);
  setValueF("pttrkiso_pt1p0", pttrkiso_pt1p0);
  setValueF("pttrkiso_pt1p5", pttrkiso_pt1p5);
}
*/
void PFTauIdEffNtupleProducer2::setGenTauMatchValues(bool genTauMatch, const reco::Candidate::LorentzVector& genTauP4, int genTauDecayMode, double genTauDeltaR, const reco::Candidate::LorentzVector& genTauLeadChHadP4, double genTauDecayLength)
{
  setValue_EnPxPyPz("genTau", genTauP4);
  setValueI("genTauDecayMode", genTauDecayMode);
  setValueI("genTauMatch", genTauMatch);
  setValueF("genTauDeltaR", genTauDeltaR);
  setValue_EnPxPyPz("genTauLeadChHad", genTauLeadChHadP4);
  setValueF("genTauDecayLength", genTauDecayLength);
}
void PFTauIdEffNtupleProducer2::setGenTauExtraValues(int genTauLchMother, double genTauLchPhPt, int genTauNDaughters, int genTauNDaPhotons)
{
  setValueI("genTauLchMother", genTauLchMother);
  setValueF("genTauLchPhPt", genTauLchPhPt);
  setValueI("genTauNDaughters", genTauNDaughters);
  setValueI("genTauNDaPhotons", genTauNDaPhotons);
} 

void PFTauIdEffNtupleProducer2::setGenElectronMatchValues(bool genElectronMatch, const reco::Candidate::LorentzVector& genElectronP4, double genElectronDeltaR)
{
  setValue_EnPxPyPz("genElectron", genElectronP4);
  setValueI("genElectronMatch", genElectronMatch);
  setValueF("genElectronDeltaR", genElectronDeltaR);
}

void PFTauIdEffNtupleProducer2::setGenMuonMatchValues(bool genMuonMatch, const reco::Candidate::LorentzVector& genMuonP4, double genMuonDeltaR)
{
  setValue_EnPxPyPz("genMuon", genMuonP4);
  setValueI("genMuonMatch", genMuonMatch);
  setValueF("genMuonDeltaR", genMuonDeltaR);
}

void PFTauIdEffNtupleProducer2::setGenQuarkOrGluonMatchValues(bool genQuarkOrGluonMatch, const reco::Candidate::LorentzVector& genQuarkOrGluonP4, double genQuarkOrGluonDeltaR, int genQuarkOrGluonPdgId)
{
  setValue_EnPxPyPz("genQuarkOrGluon", genQuarkOrGluonP4);
  setValueI("genQuarkOrGluonMatch", genQuarkOrGluonMatch);
  setValueF("genQuarkOrGluonDeltaR", genQuarkOrGluonDeltaR);
  setValueI("genQuarkOrGluonPdgId", genQuarkOrGluonPdgId);
}

void PFTauIdEffNtupleProducer2::setGenJetMatchValues(bool genJetMatch, const reco::Candidate::LorentzVector& genJetP4, double genJetDeltaR, double genJetZ, int nIsoGenChHad, float ptIsoGenChHad, float ptIsoGenNeuHad, float ptIsoGenPhoton)
{
  setValue_EnPxPyPz("genJet", genJetP4);
  setValueI("genJetMatch", genJetMatch);
  setValueF("genJetDeltaR", genJetDeltaR);
  setValueF("genJetZ", genJetZ);
  setValueI("genNIsoChHad", nIsoGenChHad);
  setValueF("genPtIsoChHad", ptIsoGenChHad);
  setValueF("genPtIsoNeuHad", ptIsoGenNeuHad);
  setValueF("genPtIsoPhoton", ptIsoGenPhoton);
}

void PFTauIdEffNtupleProducer2::setRecJetMatchValues(bool recJetMatch, const reco::Candidate::LorentzVector& recJetP4, const reco::Candidate::LorentzVector& recJetRawP4, double recJetDeltaR, double recJetZ)
{
  setValue_EnPxPyPz("recJet", recJetP4);
  setValue_EnPxPyPz("recJetRaw", recJetRawP4);
  setValueI("recJetMatch", recJetMatch);
  setValueF("recJetDeltaR", recJetDeltaR);
  setValueF("recJetZ", recJetZ);
}

void PFTauIdEffNtupleProducer2::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  setValueL("run" ,evt.run());
  setValueL("ls", evt.luminosityBlock());
  setValueL("event", evt.eventAuxiliary().event());
  //std::cout<<" run "<<evt.run()<<" ls "<<evt.luminosityBlock()<<" event "<<evt.eventAuxiliary().event()<<std::endl;
  
  assert(ntuple_);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);

  std::vector<int> pdgIdsGenTau;
  pdgIdsGenTau.push_back(-15);
  pdgIdsGenTau.push_back(+15);

  std::vector<int> pdgIdsGenElectron;
  pdgIdsGenElectron.push_back(-11);
  pdgIdsGenElectron.push_back(+11);

  std::vector<int> pdgIdsGenMuon;
  pdgIdsGenMuon.push_back(-13);
  pdgIdsGenMuon.push_back(+13);

  std::vector<int> pdgIdsGenQuarkOrGluon;
  pdgIdsGenQuarkOrGluon.push_back(-21);
  //pdgIdsGenQuarkOrGluon.push_back(-6); //not top quark
  pdgIdsGenQuarkOrGluon.push_back(-5);
  pdgIdsGenQuarkOrGluon.push_back(-4);
  pdgIdsGenQuarkOrGluon.push_back(-3);
  pdgIdsGenQuarkOrGluon.push_back(-2);
  pdgIdsGenQuarkOrGluon.push_back(-1);
  pdgIdsGenQuarkOrGluon.push_back(+1);
  pdgIdsGenQuarkOrGluon.push_back(+2);
  pdgIdsGenQuarkOrGluon.push_back(+3);
  pdgIdsGenQuarkOrGluon.push_back(+4);
  pdgIdsGenQuarkOrGluon.push_back(+5);
  //pdgIdsGenQuarkOrGluon.push_back(+6);
  pdgIdsGenQuarkOrGluon.push_back(+21);
  
  edm::Handle<pat::ElectronCollection> recVetoElectrons;
  evt.getByLabel(srcRecVetoElectrons_, recVetoElectrons);

  edm::Handle<pat::TauCollection> recTaus;
  evt.getByLabel(srcRecTaus_, recTaus);

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(srcVertices_, vertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByLabel(srcBeamSpot_,beamSpot);
  
  edm::Handle<reco::TrackCollection> tracks;
  evt.getByLabel(srcTracks_,tracks);

  double evtWeight = 1.0;
  for ( vInputTag::const_iterator srcWeight = srcWeights_.begin();
	srcWeight != srcWeights_.end(); ++srcWeight ) {
    edm::Handle<double> weight;
    evt.getByLabel(*srcWeight, weight);
    evtWeight *= (*weight);
  }
  
  //weight from MC@NLO 
  double weightevt = 1;
  try{
    edm::Handle<GenEventInfoProduct> genEvt;
    evt.getByLabel("generator",genEvt);
    weightevt = genEvt->weight();
    //std::cout<<" mc@nlo weight "<<weightevt<<std::endl;
  }
  catch(std::exception &e){ std::cerr << e.what();}

  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	genParticle != genParticles->end(); ++genParticle ) {

    if ( genParticle->pt() < minGenVisPt_ ) continue;

    unsigned numHypotheses = 0;

    bool isGenTau = false;
    for ( std::vector<int>::const_iterator pdgIdGenTau = pdgIdsGenTau.begin();
	  pdgIdGenTau != pdgIdsGenTau.end(); ++pdgIdGenTau ) {
      if ( genParticle->status() == 2 && genParticle->pdgId() == (*pdgIdGenTau) ) isGenTau = true;
    }
    if ( isGenTau ) {
      reco::Candidate::LorentzVector genTauP4 = getVisMomentum(&(*genParticle));
      std::string genTauDecayMode_string = getGenTauDecayMode(&(*genParticle));
      int genTauDecayMode_int = -1;
      if      ( genTauDecayMode_string == "oneProng0Pi0"    ) genTauDecayMode_int = reco::PFTau::kOneProng0PiZero;
      else if ( genTauDecayMode_string == "oneProng1Pi0"    ) genTauDecayMode_int = reco::PFTau::kOneProng1PiZero;
      else if ( genTauDecayMode_string == "oneProng2Pi0"    ) genTauDecayMode_int = reco::PFTau::kOneProng2PiZero;
      else if ( genTauDecayMode_string == "threeProng0Pi0"  ) genTauDecayMode_int = reco::PFTau::kThreeProng0PiZero;
      else if ( genTauDecayMode_string == "threeProng1Pi0"  ) genTauDecayMode_int = reco::PFTau::kThreeProng1PiZero;
      //else if ( genTauDecayMode_string == "oneProngOther"   ||
      //	  genTauDecayMode_string == "threeProngOther" ||
      //	  genTauDecayMode_string == "rare"            ) genTauDecayMode_int = reco::PFTau::kRareDecayMode;
      if ( genTauDecayMode_int != -1 && genTauP4.pt() > minGenVisPt_ ) {
	const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, genTauP4);
	double genTauDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genTauP4) : 9.9;
	bool genTauMatch = (genTauDeltaR < 0.3);
	int genTauDecayMode = genTauDecayMode_int;
	reco::Candidate::LorentzVector genTauLeadChHadP4 = getLeadChHadMomentum(&(*genParticle));
	double genTauDecayLength_ = compGenTauDecayDistance(&(*genParticle));
	//std::cout << "gen Tau Vis Pt "<<genTauP4.pt()<<" eta "<<genTauP4.eta()<<" phi "<<genTauP4.phi()<<std::endl;    
	//const reco::Vertex* pv = 0;
	//if(vertices->size() > 0)pv = &((*vertices)[0]);
	setRecTauValues(recTau_matched, *recVetoElectrons, (*vertices), *tracks, *beamSpot);
	//setRecTauExtraValues();
	setGenTauMatchValues(genTauMatch, genTauP4, genTauDecayMode, genTauDeltaR, genTauLeadChHadP4, genTauDecayLength_);
	setGenElectronMatchValues(false);
	setGenMuonMatchValues(false);
	setGenJetMatchValues(false);

	std::vector<const reco::GenParticle*> genTauStableDaughters;
	findDaughters((&(*genParticle)), genTauStableDaughters, 1);
	//std::cout << "gen Tau Vis Pt "<<genTauP4.pt()<<" eta "<<genTauP4.eta()<<" phi "<<genTauP4.phi()<<std::endl;
	int nDaughters = genTauStableDaughters.size();
	int nDPhotons = 0;
	const reco::GenParticle* leadChParticle = 0;
	float lchPt = 0; float lchDPhPt = 0;
	for ( std::vector<const reco::GenParticle*>::const_iterator daughter = genTauStableDaughters.begin();
	      daughter != genTauStableDaughters.end(); ++daughter ) {
	  //std::cout << "daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","                                                                 
	  //<< " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
	  if((*daughter)->pdgId() == 22){
	    nDPhotons++;
	    if((*daughter)->mother()->pdgId() != 111) lchDPhPt += (*daughter)->pt();
	  }
	  if ( (*daughter)->status() == 1 && (*daughter)->charge() != 0 ) {
	    if((*daughter)->pt() > lchPt){
	      lchPt = (*daughter)->pt();
	      leadChParticle = (&(*(*daughter)));
	    }
	  }
	}
	int lchMother = 0;
	if(leadChParticle){
	  lchMother = leadChParticle->mother()->pdgId();
	  //if(leadChParticle->mother()->pdgId() == leadChParticle->pdgId()){
	  //  //std::cout<<"mother Id "<<leadChParticle->mother()->pdgId()<<" pt "<<leadChParticle->mother()->pt()<<" lch "<<leadChParticle->pdgId()<< "pt "<<leadChParticle->pt()<<std::endl;
	  //  lchDPhPt = (leadChParticle->mother()->p4() - leadChParticle->p4()).pt();
	  //}
	}
	setGenTauExtraValues(lchMother, lchDPhPt, nDaughters, nDPhotons);

	++numHypotheses;
      }
    }

    bool isGenElectron = false;
    for ( std::vector<int>::const_iterator pdgIdGenElectron = pdgIdsGenElectron.begin();
	  pdgIdGenElectron != pdgIdsGenElectron.end(); ++pdgIdGenElectron ) {
      if ( genParticle->status() == 1 && genParticle->pdgId() == (*pdgIdGenElectron) ) isGenElectron = true;
    }
    if ( isGenElectron && genParticle->pt() > minGenVisPt_ ) {
      reco::Candidate::LorentzVector genElectronP4 = genParticle->p4();
      const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, genElectronP4);
      double genElectronDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genElectronP4) : 9.9;
      bool genElectronMatch = (genElectronDeltaR < 0.3);
      //const reco::Vertex* pv = 0;
      //if(vertices->size() > 0)pv = &((*vertices)[0]);
      setRecTauValues(recTau_matched, *recVetoElectrons, (*vertices), *tracks, *beamSpot);
      //setRecTauExtraValues();
      setGenElectronMatchValues(genElectronMatch, genElectronP4, genElectronDeltaR);
      setGenTauMatchValues(false);
      setGenMuonMatchValues(false);
      setGenTauExtraValues();
      setGenJetMatchValues(false);
      ++numHypotheses;
    }
  
    bool isGenMuon = false;
    for ( std::vector<int>::const_iterator pdgIdGenMuon = pdgIdsGenMuon.begin();
	  pdgIdGenMuon != pdgIdsGenMuon.end(); ++pdgIdGenMuon ) {
      if ( genParticle->status() == 1 && genParticle->pdgId() == (*pdgIdGenMuon) ) isGenMuon = true;
    }
    if ( isGenMuon && genParticle->pt() > minGenVisPt_ ) {
      reco::Candidate::LorentzVector genMuonP4 = genParticle->p4();
      const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, genMuonP4);
      double genMuonDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genMuonP4) : 9.9;
      bool genMuonMatch = (genMuonDeltaR < 0.3);
      //const reco::Vertex* pv = 0;
      //if(vertices->size() > 0)pv = &((*vertices)[0]);
      setRecTauValues(recTau_matched, *recVetoElectrons, (*vertices), *tracks, *beamSpot);
      //setRecTauExtraValues();
      setGenMuonMatchValues(genMuonMatch, genMuonP4, genMuonDeltaR);
      setGenTauMatchValues(false);
      setGenElectronMatchValues(false);
      setGenTauExtraValues();
      setGenJetMatchValues(false);
      ++numHypotheses;
    }

    if ( numHypotheses > 1 ) 
      edm::LogWarning("PFTauIdEffNtupleProducer2::analyze")
	<< " Matching between reconstructed PFTau and generator level tau-jets, electrons and is ambiguous --> skipping !!";
    if ( numHypotheses != 1 ) continue;
    
    setValueI("numVertices", vertices->size());
    double vertexZ_ = (vertices->size() > 0) ? (*vertices)[0].z() : -99.;
    setValueF("vertexZ", vertexZ_);
    setValueF("evtWeight", evtWeight);
    setValueF("genEvtWeight", weightevt);
    
//--- fill all computed quantities into TTree
    assert(ntuple_);
    ntuple_->Fill();
  }

  typedef edm::View<reco::Jet> JetView;
  edm::Handle<JetView> genJets;
  evt.getByLabel(srcGenJets_, genJets);

  for ( JetView::const_iterator genJet = genJets->begin();
	genJet != genJets->end(); ++genJet ) {
    
    if ( genJet->pt() < minGenVisPt_ ) continue;

    const reco::GenParticle* bestGenParticleMatch = 0;
    for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	  genParticle != genParticles->end(); ++genParticle ) {
      if ( deltaR(genParticle->p4(), genJet->p4()) < 0.3 ) {
	bool isGenQuarkOrGluon = false;
	for ( std::vector<int>::const_iterator pdgIdGenQuarkOrGluon = pdgIdsGenQuarkOrGluon.begin();
	      pdgIdGenQuarkOrGluon != pdgIdsGenQuarkOrGluon.end(); ++pdgIdGenQuarkOrGluon ) {
	  if ( genParticle->pdgId() == (*pdgIdGenQuarkOrGluon) ) isGenQuarkOrGluon = true;
	}
	if ( isGenQuarkOrGluon && (!bestGenParticleMatch || genParticle->energy() > bestGenParticleMatch->energy()) ) bestGenParticleMatch = &(*genParticle);
      }
    }
    if ( bestGenParticleMatch && genJet->pt() > minGenVisPt_ && bestGenParticleMatch->pt() > minGenVisPt_ && TMath::Abs(bestGenParticleMatch->eta()) < 3.0 ) {
      reco::Candidate::LorentzVector genQuarkOrGluonP4 = bestGenParticleMatch->p4();
      const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, genQuarkOrGluonP4);
      double genQuarkOrGluonDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genQuarkOrGluonP4) : 9.9;
      bool genQuarkOrGluonMatch = (genQuarkOrGluonDeltaR < 0.3);
      //const reco::Vertex* pv = 0;
      //if(vertices->size() > 0)pv = &((*vertices)[0]);
      setRecTauValues(recTau_matched, *recVetoElectrons, (*vertices), *tracks, *beamSpot);
      setGenQuarkOrGluonMatchValues(genQuarkOrGluonMatch, genQuarkOrGluonP4, genQuarkOrGluonDeltaR, bestGenParticleMatch->pdgId());
      reco::Candidate::LorentzVector genJetP4 = genJet->p4();
      double genJetDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genJetP4) : 9.9;
      bool genJetMatch = (genJetDeltaR < 0.3);

      //compute genJet.vZ from constituents
      //std::vector <const GenParticle*> jetConstituents = genJet->getGenConstituents();
      reco::Jet::Constituents jetConstituents = genJet->getJetConstituents();
      
      double genJetZ =0, genJetE =0;
      for(size_t ip = 0; ip < jetConstituents.size(); ip++){
	if(TMath::Abs(jetConstituents[ip]->pdgId()) == 12 || TMath::Abs(jetConstituents[ip]->pdgId()) == 14 ||
	   TMath::Abs(jetConstituents[ip]->pdgId()) == 16) continue;
	genJetZ += (jetConstituents[ip]->vz() * jetConstituents[ip]->energy());
	genJetE += jetConstituents[ip]->energy();
      }
      if(genJetE != 0) genJetZ = genJetZ/genJetE;

      //Compute isolation at generator level:
      int nGenChHad = 0;
      float ptGenChHad = 0; float ptGenNeuHad = 0; float ptGenPhoton = 0;
      if(recTau_matched){
	double DeltaR_sig = std::max(0.05, std::min(0.10, 3.0/recTau_matched->pt()));
	double DeltaR_iso  = 0.5;
	for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	      genParticle != genParticles->end(); ++genParticle ) {
	  if(genParticle->status() != 1 || genParticle->pt() < 0.5) continue;
	  
	  double deltaR_wrtTau = deltaR(genParticle->p4(), recTau_matched->p4());
	  if(deltaR_wrtTau > DeltaR_sig && deltaR_wrtTau < DeltaR_iso){
	    if(genParticle->charge() != 0){
	      nGenChHad++;
	      ptGenChHad += genParticle->pt();
	    }
	    else if(genParticle->charge() == 0){
	      ptGenNeuHad += genParticle->pt();
	      if(TMath::Abs(genParticle->pdgId()) == 22)ptGenPhoton += genParticle->pt();
	    }
	  }
	}
      }
      
      setGenJetMatchValues(genJetMatch, genJetP4, genJetDeltaR, genJetZ, nGenChHad, ptGenChHad, ptGenNeuHad, ptGenPhoton);
      setGenTauMatchValues(false);
      setGenElectronMatchValues(false);
      setGenMuonMatchValues(false);
      setGenTauExtraValues();
    }

    setValueI("numVertices", vertices->size());
    double vertexZ_ = (vertices->size() > 0) ? (*vertices)[0].z() : -99.;
    setValueF("vertexZ", vertexZ_);
    setValueF("evtWeight", evtWeight);

//--- fill all computed quantities into TTree
    assert(ntuple_);
    ntuple_->Fill();
  }

  edm::Handle<pat::JetCollection> recJets;
  evt.getByLabel(srcRecJets_, recJets);

  for ( pat::JetCollection::const_iterator recJet = recJets->begin();
	recJet != recJets->end(); ++recJet ) {
    
    if ( recJet->pt() < minRecJetPt_ ) continue;
    
    reco::Candidate::LorentzVector recJetP4 = recJet->p4();
    const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, recJetP4);
    double recJetDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), recJetP4) : 9.9;
    bool recJetMatch = (recJetDeltaR < 0.3);
    //const reco::Vertex* pv = 0;
    //if(vertices->size() > 0)pv = &((*vertices)[0]);
    setRecTauValues(recTau_matched, *recVetoElectrons, (*vertices), *tracks, *beamSpot);
    setGenQuarkOrGluonMatchValues(false);
    setGenJetMatchValues(false);
    setGenTauMatchValues(false);
    setGenElectronMatchValues(false);
    setGenMuonMatchValues(false);
    setGenTauExtraValues();
    reco::Candidate::LorentzVector recJetRawP4 = recJet->correctedJet("Uncorrected").p4();

    //compute the Vertex.Z of refJet    
    double recJetZ =0, recJetE =0;
    //reco::TrackRefVector recJetTracks = recJet->originalObject()->getTrackRefs(); //associatedTracks();
    //if(recJetTracks.size() <= 0)recJetZ =-999.;
    const std::vector<reco::PFCandidatePtr> pfConstituents = recJet->getPFConstituents ();
    for(size_t it = 0; it < pfConstituents.size(); it++){
      if(pfConstituents[it]->trackRef().isNonnull()){
	recJetZ += (pfConstituents[it]->trackRef()->vz() * pfConstituents[it]->trackRef()->p());
	recJetE += pfConstituents[it]->trackRef()->p();
      } 
    }
    recJetZ = (recJetE != 0) ? recJetZ/recJetE : -999.;

    setRecJetMatchValues(recJetMatch, recJetP4, recJetRawP4, recJetDeltaR, recJetZ);

    setValueI("numVertices", vertices->size());
    double vertexZ_ = (vertices->size() > 0) ? (*vertices)[0].z() : -99.;
    setValueF("vertexZ", vertexZ_);
    setValueF("evtWeight", evtWeight);

    //--- fill all computed quantities into TTree                                                                                                                                      
    assert(ntuple_);
    ntuple_->Fill();
  }

}

void PFTauIdEffNtupleProducer2::addBranchF(const std::string& name) 
{
  assert(branches_.count(name) == 0);
  std::string name_and_format = name + "/F";
  ntuple_->Branch(name.c_str(), &branches_[name].valueF_, name_and_format.c_str());
}

void PFTauIdEffNtupleProducer2::addBranchI(const std::string& name) 
{
  assert(branches_.count(name) == 0);
  std::string name_and_format = name + "/I";
  ntuple_->Branch(name.c_str(), &branches_[name].valueI_, name_and_format.c_str());
}

void PFTauIdEffNtupleProducer2::addBranchL(const std::string& name) 
{
  assert(branches_.count(name) == 0);
  std::string name_and_format = name + "/L";
  ntuple_->Branch(name.c_str(), &branches_[name].valueL_, name_and_format.c_str());
}

void PFTauIdEffNtupleProducer2::addBranchIV(const std::string& name)
{
  assert(branches_.count(name) == 0);
  ntuple_->Branch(name.c_str(), "std::vector<Int_t>", &branches_[name].valueIV_);
}

void PFTauIdEffNtupleProducer2::addBranchFV(const std::string& name)
{
  assert(branches_.count(name) == 0);
  ntuple_->Branch(name.c_str(), "std::vector<Float_t>", &branches_[name].valueFV_);
}

void PFTauIdEffNtupleProducer2::printBranches(std::ostream& stream)
{
  stream << "<PFTauIdEffNtupleProducer2::printBranches>:" << std::endl;
  stream << " registered branches for module = " << moduleLabel_ << std::endl;
  for ( branchMap::const_iterator branch = branches_.begin();
	branch != branches_.end(); ++branch ) {
    stream << " " << branch->first << std::endl;
  }
  stream << std::endl;
}

void PFTauIdEffNtupleProducer2::setValueF(const std::string& name, double value) 
{
  if ( verbosity_ ) std::cout << "branch = " << name << ": value = " << value << std::endl;
  branchMap::iterator branch = branches_.find(name);
  if ( branch != branches_.end() ) {
    branch->second.valueF_ = value;
  } else {
    throw cms::Exception("InvalidParameter") 
      << "No branch with name = " << name << " defined !!\n";
  }
}

void PFTauIdEffNtupleProducer2::setValueI(const std::string& name, int value) 
{
  if ( verbosity_ ) std::cout << "branch = " << name << ": value = " << value << std::endl;
  branchMap::iterator branch = branches_.find(name);
  if ( branch != branches_.end() ) {
    branch->second.valueI_ = value;
  } else {
    throw cms::Exception("InvalidParameter") 
      << "No branch with name = " << name << " defined !!\n";
  }
}

void PFTauIdEffNtupleProducer2::setValueL(const std::string& name, long value) 
{
  if ( verbosity_ ) std::cout << "branch = " << name << ": value = " << value << std::endl;
  branchMap::iterator branch = branches_.find(name);
  if ( branch != branches_.end() ) {
    branch->second.valueL_ = value;
  } else {
    throw cms::Exception("InvalidParameter") 
      << "No branch with name = " << name << " defined !!\n";
  }
}

void PFTauIdEffNtupleProducer2::setValueIV(const std::string& name, std::vector<Int_t> value)
{
  if ( verbosity_ ) std::cout << "branch = " << name << ": value = " << value.size() << std::endl;
  branchMap::iterator branch = branches_.find(name);
  if ( branch != branches_.end() ) {
    branch->second.valueIV_ = value;
  } else {
    throw cms::Exception("InvalidParameter")
      << "No branch with name = " << name << " defined !!\n";
  }
}

void PFTauIdEffNtupleProducer2::setValueFV(const std::string& name, std::vector<Float_t> value)
{
  if ( verbosity_ ) std::cout << "branch = " << name << ": value = " << value.size() << std::endl;
  branchMap::iterator branch = branches_.find(name);
  if ( branch != branches_.end() ) {
    branch->second.valueFV_ = value;
  } else {
    throw cms::Exception("InvalidParameter")
      << "No branch with name = " << name << " defined !!\n";
  }
}

//
//-------------------------------------------------------------------------------
//

void PFTauIdEffNtupleProducer2::addBranch_EnPxPyPz(const std::string& name) 
{
  addBranchF(std::string(name).append("En"));
  addBranchF(std::string(name).append("P"));
  addBranchF(std::string(name).append("Px"));
  addBranchF(std::string(name).append("Py"));
  addBranchF(std::string(name).append("Pz"));
  addBranchF(std::string(name).append("M"));
  addBranchF(std::string(name).append("Eta"));
  addBranchF(std::string(name).append("Phi"));
  addBranchF(std::string(name).append("Pt"));
}
//
//-------------------------------------------------------------------------------
//

void PFTauIdEffNtupleProducer2::setValue_EnPxPyPz(const std::string& name, const reco::Candidate::LorentzVector& p4)
{
  setValueF(std::string(name).append("En"), p4.E());
  setValueF(std::string(name).append("P"), p4.P());
  setValueF(std::string(name).append("Px"), p4.px());
  setValueF(std::string(name).append("Py"), p4.py());
  setValueF(std::string(name).append("Pz"), p4.pz());
  setValueF(std::string(name).append("M"), p4.M());
  setValueF(std::string(name).append("Eta"), p4.eta());
  setValueF(std::string(name).append("Phi"), p4.phi());
  setValueF(std::string(name).append("Pt"), p4.pt());
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PFTauIdEffNtupleProducer2);
