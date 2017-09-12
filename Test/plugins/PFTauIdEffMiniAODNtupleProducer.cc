#include "TauAnalysis/Test/plugins/PFTauIdEffMiniAODNtupleProducer.h"

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

int PFTauIdEffMiniAODNtupleProducer::verbosity_ = 0;

PFTauIdEffMiniAODNtupleProducer::PFTauIdEffMiniAODNtupleProducer(const edm::ParameterSet& cfg) 
  : moduleLabel_(cfg.getParameter<std::string>("@module_label")),
    ntuple_(0)
{
  srcGenParticles_ = consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("srcGenParticles"));
  srcGenJets_ = consumes<edm::View<reco::Jet>>(cfg.getParameter<edm::InputTag>("srcGenJets"));
  srcRecVetoElectrons_ = consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("srcRecVetoElectrons"));
  srcRecTaus_ = consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("srcRecTaus"));
  srcRecJets_ = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("srcRecJets"));
  srcVertices_ = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("srcVertices"));

  srcWeights_ = cfg.getParameter<vInputTag>("srcWeights");

  

  tauIdDiscriminators_.push_back("decayModeFindingNewDMs");
  //tauIdDiscriminators_.push_back("decayModeFindingOldDMs");
  tauIdDiscriminators_.push_back("decayModeFinding");
  tauIdDiscriminators_.push_back("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  tauIdDiscriminators_.push_back("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  tauIdDiscriminators_.push_back("byTightCombinedIsolationDeltaBetaCorr3Hits");
  tauIdDiscriminators_.push_back("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03");
  tauIdDiscriminators_.push_back("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03");
  tauIdDiscriminators_.push_back("byTightCombinedIsolationDeltaBetaCorr3HitsdR03");
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
  tauIdDiscriminators_.push_back("byIsolationMVA3oldDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVA3oldDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVA3oldDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVA3oldDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVA3oldDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVA3oldDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVA3oldDMwLT");
  tauIdDiscriminators_.push_back("byIsolationMVA3newDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVA3newDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVA3newDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVA3newDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVA3newDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVA3newDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVA3newDMwLT");
  tauIdDiscriminators_.push_back("byIsolationMVArun2v1DBoldDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVArun2v1DBoldDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVArun2v1DBoldDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVArun2v1DBoldDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVArun2v1DBoldDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVArun2v1DBoldDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVArun2v1DBoldDMwLT");
  tauIdDiscriminators_.push_back("byIsolationMVArun2v1DBnewDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVArun2v1DBnewDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVArun2v1DBnewDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVArun2v1DBnewDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVArun2v1DBnewDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVArun2v1DBnewDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVArun2v1DBnewDMwLT");
  tauIdDiscriminators_.push_back("byIsolationMVArun2v1PWoldDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVArun2v1PWoldDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVArun2v1PWoldDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVArun2v1PWoldDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVArun2v1PWoldDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVArun2v1PWoldDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVArun2v1PWoldDMwLT");
  tauIdDiscriminators_.push_back("byIsolationMVArun2v1PWnewDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVArun2v1PWnewDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVArun2v1PWnewDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVArun2v1PWnewDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVArun2v1PWnewDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVArun2v1PWnewDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVArun2v1PWnewDMwLT");
  tauIdDiscriminators_.push_back("chargedIsoPtSumdR03");
  tauIdDiscriminators_.push_back("neutralIsoPtSumdR03");
  tauIdDiscriminators_.push_back("puCorrPtSumdR03");
  tauIdDiscriminators_.push_back("neutralIsoPtSumWeightdR03");
  tauIdDiscriminators_.push_back("footprintCorrectiondR03");
  tauIdDiscriminators_.push_back("photonPtSumOutsideSignalConedR03");
  tauIdDiscriminators_.push_back("byIsolationMVArun2v1DBdR03oldDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVArun2v1DBdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVArun2v1DBdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVArun2v1DBdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVArun2v1DBdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byIsolationMVArun2v1PWdR03oldDMwLTraw");
  tauIdDiscriminators_.push_back("byVLooseIsolationMVArun2v1PWdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byLooseIsolationMVArun2v1PWdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byMediumIsolationMVArun2v1PWdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byTightIsolationMVArun2v1PWdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byVTightIsolationMVArun2v1PWdR03oldDMwLT");
  tauIdDiscriminators_.push_back("byVVTightIsolationMVArun2v1PWdR03oldDMwLT");
  tauIdDiscriminators_.push_back("againstElectronMVA6raw");
  tauIdDiscriminators_.push_back("againstElectronMVA6category");
  tauIdDiscriminators_.push_back("againstElectronVLooseMVA6");
  tauIdDiscriminators_.push_back("againstElectronLooseMVA6");
  tauIdDiscriminators_.push_back("againstElectronMediumMVA6");
  tauIdDiscriminators_.push_back("againstElectronTightMVA6");
  tauIdDiscriminators_.push_back("againstElectronVTightMVA6");
  tauIdDiscriminators_.push_back("againstMuonLoose3");
  tauIdDiscriminators_.push_back("againstMuonTight3");

  minGenVisPt_ = 15.;
  minRecJetPt_ = 20.;
}

PFTauIdEffMiniAODNtupleProducer::~PFTauIdEffMiniAODNtupleProducer()
{
// nothing to be done yet...
}

void PFTauIdEffMiniAODNtupleProducer::beginJob()
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
  addBranch_EnPxPyPz("leadPFCand");
  addBranch_EnPxPyPz("leadPFChargedHadrCand");
  addBranchF("Tau_dzWrtPV");
  addBranchF("Tau_z");
  for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
    addBranchF(*tauIdDiscriminator);
  }
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
  addBranchF("recJetNHEn");
  addBranchF("recJetNEMEn");
  addBranchF("recJetCHEn");
  addBranchF("recJetCEMEn");
  addBranchI("recJetNHM");
  addBranchI("recJetCHM");
  addBranchI("numVertices");
  addBranchF("vertexZ");
  addBranchF("evtWeight");
}

const pat::Tau* PFTauIdEffMiniAODNtupleProducer::findMatchingRecTau(const pat::TauCollection& recTaus, const reco::Candidate::LorentzVector& genParticleP4)
{
  const pat::Tau* recTau_matched = 0;
  
  double genTauDeltaR = 9.9;
  for ( pat::TauCollection::const_iterator recTau = recTaus.begin();
	recTau != recTaus.end(); ++recTau ) {
    
    //if ( recTau->pt() < 15. ) continue;
    
    double dR = deltaR(recTau->p4(), genParticleP4);
    if ( dR < 0.3 && dR < genTauDeltaR ) {
      genTauDeltaR = dR;
      recTau_matched = &(*recTau);
    }
  }

  return recTau_matched;
}

const pat::Electron* PFTauIdEffMiniAODNtupleProducer::findMatchingElectronVeto(const pat::Tau& recTau, const pat::ElectronCollection& recVetoElectrons)
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
/*
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
*/
double  PFTauIdEffMiniAODNtupleProducer::compGenTauDecayDistance(const reco::GenParticle* genTau)
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

void PFTauIdEffMiniAODNtupleProducer::setRecTauValues(const pat::Tau* recTau_matched, const pat::ElectronCollection& recVetoElectrons)
{

  if ( recTau_matched ) {
    setValue_EnPxPyPz("recTau", recTau_matched->p4());
    setValue_EnPxPyPz("recTauAlternate", recTau_matched->alternatLorentzVect());
    setValueI("recTauDecayMode", recTau_matched->decayMode());
    if ( recTau_matched->leadCand().isNonnull() ) setValue_EnPxPyPz("leadPFCand", recTau_matched->leadCand()->p4());
    else setValue_EnPxPyPz("leadPFCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    if ( recTau_matched->leadChargedHadrCand().isNonnull() ){
      setValue_EnPxPyPz("leadPFChargedHadrCand", recTau_matched->leadChargedHadrCand()->p4());
    }
    else setValue_EnPxPyPz("leadPFChargedHadrCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(recTau_matched->leadChargedHadrCand().get());
    double dzWrtPV = (packedLeadTauCand != 0) ? packedLeadTauCand->dz() : 999.;
    setValueF("Tau_dzWrtPV", dzWrtPV);
    setValueF("Tau_z", 0.);
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
      setValueF(*tauIdDiscriminator, recTau_matched->tauID(*tauIdDiscriminator));
    }

    const pat::Electron* ElectronVeto_matched = findMatchingElectronVeto(*recTau_matched, recVetoElectrons);
    double ElectronVetoDeltaR = ( ElectronVeto_matched ) ? deltaR(ElectronVeto_matched->p4(), recTau_matched->p4()) : 9.9;
    bool ElectronVetoMatch = (ElectronVetoDeltaR < 0.3);
    setValueF("ElectronVetoDeltaR", ElectronVetoDeltaR);
    setValueI("ElectronVetoMatch", ElectronVetoMatch);
  } else {
    setValue_EnPxPyPz("recTau", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValueI("recTauDecayMode", -1);
    setValue_EnPxPyPz("leadPFCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValue_EnPxPyPz("leadPFChargedHadrCand", reco::Candidate::LorentzVector(0.,0.,0.,0.));
    setValueF("Tau_dzWrtPV", 99.);
    setValueF("Tau_z", -999.);
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
      setValueF(*tauIdDiscriminator, 0.);
    }
    setValueF("ElectronVetoDeltaR", 0.);
    setValueI("ElectronVetoMatch", 0.);
  }
}

void PFTauIdEffMiniAODNtupleProducer::setGenTauMatchValues(bool genTauMatch, const reco::Candidate::LorentzVector& genTauP4, int genTauDecayMode, double genTauDeltaR, const reco::Candidate::LorentzVector& genTauLeadChHadP4, double genTauDecayLength)
{
  setValue_EnPxPyPz("genTau", genTauP4);
  setValueI("genTauDecayMode", genTauDecayMode);
  setValueI("genTauMatch", genTauMatch);
  setValueF("genTauDeltaR", genTauDeltaR);
  setValue_EnPxPyPz("genTauLeadChHad", genTauLeadChHadP4);
  setValueF("genTauDecayLength", genTauDecayLength);
}
void PFTauIdEffMiniAODNtupleProducer::setGenTauExtraValues(int genTauLchMother, double genTauLchPhPt, int genTauNDaughters, int genTauNDaPhotons)
{
  setValueI("genTauLchMother", genTauLchMother);
  setValueF("genTauLchPhPt", genTauLchPhPt);
  setValueI("genTauNDaughters", genTauNDaughters);
  setValueI("genTauNDaPhotons", genTauNDaPhotons);
} 

void PFTauIdEffMiniAODNtupleProducer::setGenElectronMatchValues(bool genElectronMatch, const reco::Candidate::LorentzVector& genElectronP4, double genElectronDeltaR)
{
  setValue_EnPxPyPz("genElectron", genElectronP4);
  setValueI("genElectronMatch", genElectronMatch);
  setValueF("genElectronDeltaR", genElectronDeltaR);
}

void PFTauIdEffMiniAODNtupleProducer::setGenMuonMatchValues(bool genMuonMatch, const reco::Candidate::LorentzVector& genMuonP4, double genMuonDeltaR)
{
  setValue_EnPxPyPz("genMuon", genMuonP4);
  setValueI("genMuonMatch", genMuonMatch);
  setValueF("genMuonDeltaR", genMuonDeltaR);
}

void PFTauIdEffMiniAODNtupleProducer::setGenQuarkOrGluonMatchValues(bool genQuarkOrGluonMatch, const reco::Candidate::LorentzVector& genQuarkOrGluonP4, double genQuarkOrGluonDeltaR, int genQuarkOrGluonPdgId)
{
  setValue_EnPxPyPz("genQuarkOrGluon", genQuarkOrGluonP4);
  setValueI("genQuarkOrGluonMatch", genQuarkOrGluonMatch);
  setValueF("genQuarkOrGluonDeltaR", genQuarkOrGluonDeltaR);
  setValueI("genQuarkOrGluonPdgId", genQuarkOrGluonPdgId);
}

void PFTauIdEffMiniAODNtupleProducer::setGenJetMatchValues(bool genJetMatch, const reco::Candidate::LorentzVector& genJetP4, double genJetDeltaR)
{
  setValue_EnPxPyPz("genJet", genJetP4);
  setValueI("genJetMatch", genJetMatch);
  setValueF("genJetDeltaR", genJetDeltaR);
}

void PFTauIdEffMiniAODNtupleProducer::setRecJetMatchValues(bool recJetMatch, const reco::Candidate::LorentzVector& recJetP4, const reco::Candidate::LorentzVector& recJetRawP4, double recJetDeltaR)
{
  setValue_EnPxPyPz("recJet", recJetP4);
  setValue_EnPxPyPz("recJetRaw", recJetRawP4);
  setValueI("recJetMatch", recJetMatch);
  setValueF("recJetDeltaR", recJetDeltaR);
}
void PFTauIdEffMiniAODNtupleProducer::setRecJetIDValues(float recJetNHEn, float recJetNEMEn, float recJetCHEn, float recJetCEMEn, float recJetNHM, float recJetCHM)
{
  setValueF("recJetNHEn", recJetNHEn);
  setValueF("recJetNEMEn", recJetNEMEn);
  setValueF("recJetCHEn", recJetCHEn);
  setValueF("recJetCEMEn", recJetCEMEn);
  setValueI("recJetNHM", recJetNHM);
  setValueI("recJetCHM", recJetCHM);
}

void PFTauIdEffMiniAODNtupleProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  setValueL("run" ,evt.run());
  setValueL("ls", evt.luminosityBlock());
  setValueL("event", evt.eventAuxiliary().event());
  //std::cout<<" run "<<evt.run()<<" ls "<<evt.luminosityBlock()<<" event "<<evt.eventAuxiliary().event()<<std::endl;
  
  assert(ntuple_);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(srcGenParticles_, genParticles);

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
  evt.getByToken(srcRecVetoElectrons_, recVetoElectrons);

  edm::Handle<pat::TauCollection> recTaus;
  evt.getByToken(srcRecTaus_, recTaus);

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(srcVertices_, vertices);

  double evtWeight = 1.0;
  //for ( vInputTag::const_iterator srcWeight = srcWeights_.begin();
  //	srcWeight != srcWeights_.end(); ++srcWeight ) {
  //  edm::Handle<double> weight;
  //  evt.getByLabel(*srcWeight, weight);
  //  evtWeight *= (*weight);
  //}
  
  //weight from MC@NLO 
  double weightevt = 1;
  //try{
  //  edm::Handle<GenEventInfoProduct> genEvt;
  //  evt.getByLabel("generator",genEvt);
  //  weightevt = genEvt->weight();
  //  //std::cout<<" mc@nlo weight "<<weightevt<<std::endl;
  //}
  //catch(std::exception &e){ std::cerr << e.what();}

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
      //if ( genTauDecayMode_int != -1 && genTauP4.pt() > minGenVisPt_ ) {
      if (genTauP4.pt() > minGenVisPt_ ) {
	const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, genTauP4);
	double genTauDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genTauP4) : 9.9;
	bool genTauMatch = (genTauDeltaR < 0.3);
	int genTauDecayMode = genTauDecayMode_int;
	reco::Candidate::LorentzVector genTauLeadChHadP4 = getLeadChHadMomentum(&(*genParticle));
	double genTauDecayLength_ = compGenTauDecayDistance(&(*genParticle));
	setRecTauValues(recTau_matched, *recVetoElectrons);
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
      setRecTauValues(recTau_matched, *recVetoElectrons);
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
      setRecTauValues(recTau_matched, *recVetoElectrons);
      //setRecTauExtraValues();
      setGenMuonMatchValues(genMuonMatch, genMuonP4, genMuonDeltaR);
      setGenTauMatchValues(false);
      setGenElectronMatchValues(false);
      setGenTauExtraValues();
      setGenJetMatchValues(false);
      ++numHypotheses;
    }

    if ( numHypotheses > 1 ) 
      edm::LogWarning("PFTauIdEffMiniAODNtupleProducer::analyze")
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
  evt.getByToken(srcGenJets_, genJets);

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
      setRecTauValues(recTau_matched, *recVetoElectrons);
      setGenQuarkOrGluonMatchValues(genQuarkOrGluonMatch, genQuarkOrGluonP4, genQuarkOrGluonDeltaR, bestGenParticleMatch->pdgId());
      reco::Candidate::LorentzVector genJetP4 = genJet->p4();
      double genJetDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genJetP4) : 9.9;
      bool genJetMatch = (genJetDeltaR < 0.3);

      setGenJetMatchValues(genJetMatch, genJetP4, genJetDeltaR);
      setGenTauMatchValues(false);
      setGenElectronMatchValues(false);
      setGenMuonMatchValues(false);
      setGenTauExtraValues();

      setValueI("numVertices", vertices->size());
      double vertexZ_ = (vertices->size() > 0) ? (*vertices)[0].z() : -99.;
      setValueF("vertexZ", vertexZ_);
      setValueF("evtWeight", evtWeight);
      
      //--- fill all computed quantities into TTree
      assert(ntuple_);
      ntuple_->Fill();
    }
  }

  edm::Handle<pat::JetCollection> recJets;
  evt.getByToken(srcRecJets_, recJets);

  for ( pat::JetCollection::const_iterator recJet = recJets->begin();
	recJet != recJets->end(); ++recJet ) {
    
    if ( recJet->pt() < minRecJetPt_ ) continue;
    
    reco::Candidate::LorentzVector recJetP4 = recJet->p4();
    const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, recJetP4);
    double recJetDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), recJetP4) : 9.9;
    bool recJetMatch = (recJetDeltaR < 0.3);
    setRecTauValues(recTau_matched, *recVetoElectrons);
    setGenQuarkOrGluonMatchValues(false);
    setGenJetMatchValues(false);
    setGenTauMatchValues(false);
    setGenElectronMatchValues(false);
    setGenMuonMatchValues(false);
    setGenTauExtraValues();
    reco::Candidate::LorentzVector recJetRawP4 = recJet->correctedJet("Uncorrected").p4();

    setRecJetMatchValues(recJetMatch, recJetP4, recJetRawP4, recJetDeltaR);
    
    float recJetNHEn_ = recJet->neutralHadronEnergy();
    float recJetNEMEn_ = recJet->neutralEmEnergy();
    float recJetCHEn_ = recJet->chargedHadronEnergy();
    float recJetCEMEn_ = recJet->chargedEmEnergy();
    int recJetNHM_ = recJet->neutralMultiplicity();
    int recJetCHM_ = recJet->chargedMultiplicity();
    setRecJetIDValues(recJetNHEn_, recJetNEMEn_, recJetCHEn_, recJetCEMEn_,
		      recJetNHM_, recJetCHM_);

    setValueI("numVertices", vertices->size());
    double vertexZ_ = (vertices->size() > 0) ? (*vertices)[0].z() : -99.;
    setValueF("vertexZ", vertexZ_);
    setValueF("evtWeight", evtWeight);

    //--- fill all computed quantities into TTree                                                                                                                                      
    assert(ntuple_);
    ntuple_->Fill();
  }
  
}

void PFTauIdEffMiniAODNtupleProducer::addBranchF(const std::string& name) 
{
  assert(branches_.count(name) == 0);
  std::string name_and_format = name + "/F";
  ntuple_->Branch(name.c_str(), &branches_[name].valueF_, name_and_format.c_str());
}

void PFTauIdEffMiniAODNtupleProducer::addBranchI(const std::string& name) 
{
  assert(branches_.count(name) == 0);
  std::string name_and_format = name + "/I";
  ntuple_->Branch(name.c_str(), &branches_[name].valueI_, name_and_format.c_str());
}

void PFTauIdEffMiniAODNtupleProducer::addBranchL(const std::string& name) 
{
  assert(branches_.count(name) == 0);
  std::string name_and_format = name + "/L";
  ntuple_->Branch(name.c_str(), &branches_[name].valueL_, name_and_format.c_str());
}

void PFTauIdEffMiniAODNtupleProducer::addBranchIV(const std::string& name)
{
  assert(branches_.count(name) == 0);
  ntuple_->Branch(name.c_str(), "std::vector<Int_t>", &branches_[name].valueIV_);
}

void PFTauIdEffMiniAODNtupleProducer::addBranchFV(const std::string& name)
{
  assert(branches_.count(name) == 0);
  ntuple_->Branch(name.c_str(), "std::vector<Float_t>", &branches_[name].valueFV_);
}

void PFTauIdEffMiniAODNtupleProducer::printBranches(std::ostream& stream)
{
  stream << "<PFTauIdEffMiniAODNtupleProducer::printBranches>:" << std::endl;
  stream << " registered branches for module = " << moduleLabel_ << std::endl;
  for ( branchMap::const_iterator branch = branches_.begin();
	branch != branches_.end(); ++branch ) {
    stream << " " << branch->first << std::endl;
  }
  stream << std::endl;
}

void PFTauIdEffMiniAODNtupleProducer::setValueF(const std::string& name, double value) 
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

void PFTauIdEffMiniAODNtupleProducer::setValueI(const std::string& name, int value) 
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

void PFTauIdEffMiniAODNtupleProducer::setValueL(const std::string& name, long value) 
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

void PFTauIdEffMiniAODNtupleProducer::setValueIV(const std::string& name, std::vector<Int_t> value)
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

void PFTauIdEffMiniAODNtupleProducer::setValueFV(const std::string& name, std::vector<Float_t> value)
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

void PFTauIdEffMiniAODNtupleProducer::addBranch_EnPxPyPz(const std::string& name) 
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

void PFTauIdEffMiniAODNtupleProducer::setValue_EnPxPyPz(const std::string& name, const reco::Candidate::LorentzVector& p4)
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

DEFINE_FWK_MODULE(PFTauIdEffMiniAODNtupleProducer);
