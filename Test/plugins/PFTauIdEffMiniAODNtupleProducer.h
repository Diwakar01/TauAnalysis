#ifndef TauAnalysis_Test_PFTauIdEffMiniAODNtupleProducer_h  
#define TauAnalysis_Test_PFTauIdEffMiniAODNtupleProducer_h

/** \class PFTauIdEffMiniAODNtupleProducer
 *
 * Produce an Ntuple of various quantities useful 
 * to check tau id. efficiencies and e/mu -> tau fake-rates
 *
 * \author Arun Nayak
 *
 * \version $Revision: 1.3 $
 *
 * $Id: PFTauIdEffMiniAODNtupleProducer.h,v 1.3 2012/03/08 10:31:49 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include <TTree.h>

#include <map>
#include <string>
#include <vector>
#include <ostream>

class PFTauIdEffMiniAODNtupleProducer : public edm::EDProducer 
{
 public:
  
  PFTauIdEffMiniAODNtupleProducer(const edm::ParameterSet&);
  ~PFTauIdEffMiniAODNtupleProducer();

  void produce(edm::Event&, const edm::EventSetup&);
  void beginJob();

 private:

  void addBranchF(const std::string&);
  void addBranchI(const std::string&);
  void addBranchL(const std::string&);
  void addBranchIV(const std::string&);
  void addBranchFV(const std::string&);

  void addBranch_EnPxPyPz(const std::string&);

  void printBranches(std::ostream&);

  void setValueF(const std::string&, double);
  void setValueI(const std::string&, int);
  void setValueL(const std::string&, long);
  void setValueFV(const std::string&, std::vector<Float_t>);
  void setValueIV(const std::string&, std::vector<Int_t>);

  void setValue_EnPxPyPz(const std::string&, const reco::Candidate::LorentzVector&);

  const pat::Tau* findMatchingRecTau(const pat::TauCollection&, const reco::Candidate::LorentzVector&);
  const pat::Electron* findMatchingElectronVeto(const pat::Tau&, const pat::ElectronCollection&);
  void setRecTauValues(const pat::Tau*,const pat::ElectronCollection&);
  void setGenTauMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), int = -1, double = 9.9, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = -99.);
  void setGenTauExtraValues(int genTauLchMother=0, double genTauLchPhPt=0, int genTauNDaughters=0, int genTauNDaPhotons=0);
  void setGenElectronMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9);
  void setGenMuonMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9);
  void setGenQuarkOrGluonMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9, int = 0);
  void setGenJetMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9);
  void setRecJetMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9);
  void setRecJetIDValues(float recJetNHEn = 0, float recJetNEMEn = 0, float recJetCHEn = 0, float recJetCEMEn = 0, float recJetNHM = 0, float recJetCHM = 0);
  double compGenTauDecayDistance(const reco::GenParticle* = 0);
  
  std::string moduleLabel_;

  //edm::InputTag srcGenParticles_;
  //edm::InputTag srcGenJets_;
  //edm::InputTag srcRecVetoElectrons_;
  //edm::InputTag srcRecTaus_;
  //edm::InputTag srcRecJets_;
  //edm::InputTag srcVertices_;

  edm::EDGetTokenT<reco::GenParticleCollection> srcGenParticles_;
  edm::EDGetTokenT<edm::View<reco::Jet>> srcGenJets_;
  edm::EDGetTokenT<pat::ElectronCollection> srcRecVetoElectrons_;
  edm::EDGetTokenT<pat::TauCollection> srcRecTaus_;
  edm::EDGetTokenT<pat::JetCollection> srcRecJets_;
  edm::EDGetTokenT<reco::VertexCollection> srcVertices_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights_;

  typedef std::vector<std::string> vstring;
  vstring tauIdDiscriminators_;

  double minGenVisPt_;
  double minRecJetPt_;

  struct branchEntryType
  {
    Float_t valueF_;
    Int_t valueI_;
    Long_t valueL_;
    std::vector<Int_t> valueIV_;
    std::vector<Float_t> valueFV_;
  };

  typedef std::map<std::string, branchEntryType> branchMap; // key = branch name
  branchMap branches_;

  TTree* ntuple_;

  static int verbosity_;
};

#endif


