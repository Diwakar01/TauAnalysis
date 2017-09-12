#ifndef TauAnalysis_Test_PFTauIdEffNtupleProducer2_h  
#define TauAnalysis_Test_PFTauIdEffNtupleProducer2_h

/** \class PFTauIdEffNtupleNtupleProducer2
 *
 * Produce an Ntuple of various quantities useful 
 * to check tau id. efficiencies and e/mu -> tau fake-rates
 *
 * \author Christian Veelken, LLR
 *
 * \version $Revision: 1.3 $
 *
 * $Id: PFTauIdEffNtupleProducer2.h,v 1.3 2012/03/08 10:31:49 veelken Exp $
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

#include <TTree.h>

#include <map>
#include <string>
#include <vector>
#include <ostream>

class PFTauIdEffNtupleProducer2 : public edm::EDProducer 
{
 public:
  
  PFTauIdEffNtupleProducer2(const edm::ParameterSet&);
  ~PFTauIdEffNtupleProducer2();

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
  void setRecTauValues(const pat::Tau*,const pat::ElectronCollection&, const reco::VertexCollection&, const reco::TrackCollection&, const reco::BeamSpot&);
  //void setRecTauExtraValues(int ntrkiso_pt0p5=-99, int ntrkiso_pt1p0=-99, int ntrkiso_pt1p5=-99, double pttrkiso_pt0p5=-99., double pttrkiso_pt1p0=-99., double pttrkiso_pt1p5=-99.);
  void setGenTauMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), int = -1, double = 9.9, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = -99.);
  void setGenTauExtraValues(int genTauLchMother=0, double genTauLchPhPt=0, int genTauNDaughters=0, int genTauNDaPhotons=0);
  void setGenElectronMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9);
  void setGenMuonMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9);
  void setGenQuarkOrGluonMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9, int = 0);
  void setGenJetMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9, double = -999, int = -99, float = 99., float = 99., float = 99.);
  void setRecJetMatchValues(bool, const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), const reco::Candidate::LorentzVector& = reco::Candidate::LorentzVector(0.,0.,0.,0.), double = 9.9, double = -999);

  std::string moduleLabel_;

  edm::InputTag srcGenParticles_;
  edm::InputTag srcGenJets_;
  edm::InputTag srcRecVetoElectrons_;
  edm::InputTag srcRecTaus_;
  edm::InputTag srcRecJets_;
  edm::InputTag srcVertices_;
  edm::InputTag srcTracks_;
  edm::InputTag srcBeamSpot_;

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


