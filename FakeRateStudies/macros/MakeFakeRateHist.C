#include <vector>
#include "TString.h"
#include <string>

#define  PI 3.14159265358979312e+00

float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
  float deta = fabs(eta_1-eta_2);
  float dphi = fabs(phi_1-phi_2);
  if(dphi > PI)
    dphi = 2*PI -dphi;
  return abs(sqrt(deta*deta+dphi*dphi));
}

float pileupwght(TH1F* hin,int nPu){
  
  int b = hin->FindBin(nPu);
  float wght = hin->GetBinContent(b);
  return wght;
}

void make(bool isData){
  
  TFile* file_pu = new TFile("MyDataPileupHistogram.root");
  TH1F* h_pu = (TH1F*)file_pu->Get("pileup");
  h_pu->Scale(1/h_pu->Integral());

  std::vector<TString> filenames; filenames.clear();
  if(isData)filenames.push_back("JetHT2017.root");
  else{
    //filenames.push_back("QCD_Pt_30to50_TuneCP5_13TeV.root");
    filenames.push_back("QCD_Pt_50to80_TuneCP5_13TeV_94X.root");
    filenames.push_back("QCD_Pt-80to120_TuneCP5_13TeV.root");
    filenames.push_back("QCD_Pt-120to170_94X.root");
    filenames.push_back("QCD_Pt-170to300_TuneCP5.root");
    filenames.push_back("QCD_Pt-300to470_TuneCP5.root");
    filenames.push_back("QCD_Pt-470to600_TuneCP5_13TeV.root");
    filenames.push_back("QCD_Pt-600to800_TuneCP5_13TeV.root");
    //filenames.push_back("QCD_Pt-800to1000_94X.root");
    //filenames.push_back("QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8.root");
  }
  
  std::vector<double> xsections; xsections.clear();
  if(!isData){
    //xsections.push_back(140932000);
    xsections.push_back(19204300);
    xsections.push_back(2762530);
    xsections.push_back(471100);
    xsections.push_back(117276);
    xsections.push_back(7823);
    xsections.push_back(648.2);
    xsections.push_back(186.9);
    //xsections.push_back(32.293);
  }
  Float_t Lumi = 41481.0; //in pb^-1

  TString outFileName;
  if(isData)outFileName = TString("histograms_JetHT2017.root");
  else outFileName= TString("histograms_QCD_MC2017.root");
  
  TFile* fout = new TFile(outFileName, "recreate");

  //Define histograms
  //double pTbins[19] ={20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
  double pTbins[12] ={100,120,140,160,180,200,250,300,400,500,750,1000};
  TH1F *Denom_jetPt= new TH1F("Denom_jetPt"," ",11,pTbins);Denom_jetPt->Sumw2();
  TH1F *NumeDM_Pt = new TH1F("NumeDM_Pt"," ",11,pTbins);NumeDM_Pt->Sumw2();
  TH1F *LNume_Pt = new TH1F("LNume_Pt"," ",11,pTbins);LNume_Pt->Sumw2();
  TH1F *MVALNume_Pt = new TH1F("MVALNume_Pt"," ",11,pTbins);MVALNume_Pt->Sumw2();
  TH1F *MNume_Pt = new TH1F("MNume_Pt"," ",11,pTbins);MNume_Pt->Sumw2();
  TH1F *MVAMNume_Pt = new TH1F("MVAMNume_Pt"," ",11,pTbins);MVAMNume_Pt->Sumw2();
  TH1F *TNume_Pt = new TH1F("TNume_Pt"," ",11,pTbins);TNume_Pt->Sumw2();
  TH1F *MVATNume_Pt = new TH1F("MVATNume_Pt"," ",11,pTbins);MVATNume_Pt->Sumw2();
  TH1F *MVAVTNume_Pt = new TH1F("MVAVTNume_Pt"," ",11,pTbins);MVAVTNume_Pt->Sumw2();
  TH1F *MVAVLNume_Pt = new TH1F("MVAVLNume_Pt"," ",11,pTbins);MVAVLNume_Pt->Sumw2();
  TH1F *Denom_jetEta = new TH1F("Denom_jetEta","",20,-2.5,2.5);Denom_jetEta->Sumw2();
  TH1F *NumeDM_Eta = new TH1F("NumeDM_Eta"," ",20,-2.5,2.5);NumeDM_Eta->Sumw2();
  TH1F *LNume_Eta = new TH1F("LNume_Eta"," ",20,-2.5,2.5);LNume_Eta->Sumw2();
  TH1F *MNume_Eta = new TH1F("MNume_Eta"," ",20,-2.5,2.5);MNume_Eta->Sumw2();
  TH1F *TNume_Eta = new TH1F("TNume_Eta"," ",20,-2.5,2.5);TNume_Eta->Sumw2();
  TH1F *MVALNume_Eta = new TH1F("MVALNume_Eta"," ",20,-2.5,2.5);MVALNume_Eta->Sumw2();
  TH1F *MVAVLNume_Eta = new TH1F("MVAVLNume_Eta"," ",20,-2.5,2.5);MVAVLNume_Eta->Sumw2();
  TH1F *MVAMNume_Eta = new TH1F("MVAMNume_Eta"," ",20,-2.5,2.5);MVAMNume_Eta->Sumw2();
  TH1F *MVATNume_Eta = new TH1F("MVATNume_Eta"," ",20,-2.5,2.5);MVATNume_Eta->Sumw2();
  TH1F *MVAVTNume_Eta = new TH1F("MVAVTNume_Eta"," ",20,-2.5,2.5);MVAVTNume_Eta->Sumw2();

  //Loop over input files
  for(size_t ifile = 0; ifile < filenames.size(); ifile++){
    std::cout<<" reading file: "<<filenames[ifile]<<std::endl;

    TFile* inFile = new TFile(filenames[ifile]);
    TString treeName = "tauNtuple/tree";
    TTree *tree_jet = (TTree*)inFile->Get(treeName.Data());

    TH1F* h_p = (TH1F*)h_pu->Clone("h_p");
    if(isData){
      TH1F *hMC_pu = new TH1F("hMC_pu","",1001,0,100);
      tree_jet->Draw("npuIT>>hMC_pu","npuIT > 5");
      hMC_pu->Scale(1/(hMC_pu->Integral()));
      hMC_pu->Draw();
      h_p->Divide(hMC_pu);
    }

    std::vector<float> *trig_jetPt = new std::vector<float>();
    std::vector<float> *trig_jetEta = new std::vector<float>();
    std::vector<float> *trig_jetPhi = new std::vector<float>();
    std::vector<float> *PF_jetEta = new std::vector<float>();
    std::vector<float> *PF_jetPhi = new std::vector<float>();
    std::vector<float> *PF_jetPt = new std::vector<float>();
    std::vector<float> *TauPt_ = new std::vector<float>();
    std::vector<float> *TauEta_= new std::vector<float>();
    std::vector<float> *TauPhi_= new std::vector<float>();
    std::vector<float> *DM_find= new std::vector<float>();
    std::vector<float> *loose_iso= new std::vector<float>();
    std::vector<float> *medium_iso= new std::vector<float>();
    std::vector<float> *tight_iso= new std::vector<float>();
    std::vector<float> *MVADB_find= new std::vector<float>();
    std::vector<float> *MVAVloose_iso= new std::vector<float>();
    std::vector<float> *MVAloose_iso= new std::vector<float>();
    std::vector<float> *MVAmedium_iso= new std::vector<float>();
    std::vector<float> *MVAtight_iso= new std::vector<float>();
    std::vector<float> *MVAVtight_iso= new std::vector<float>();
    std::vector<float> *MVAVVtight_iso= new std::vector<float>();
    std::vector<string> *triggerResults_ = new std::vector<string>();

    std::vector<float> *agstEleTightMVA_ = new std::vector<float>();
    std::vector<float> *agstMuTight_ = new std::vector<float>();

    tree_jet->SetBranchAddress("hltObjPt",&trig_jetPt);
    tree_jet->SetBranchAddress("hltObjEta",&trig_jetEta);
    tree_jet->SetBranchAddress("hltObjPhi",&trig_jetPhi);
    tree_jet->SetBranchAddress("pfjetEta",&PF_jetEta);
    tree_jet->SetBranchAddress("pfjetPt",&PF_jetPt);
    tree_jet->SetBranchAddress("pfjetPhi",&PF_jetPhi);
    tree_jet->SetBranchAddress("tauPt",&TauPt_);
    tree_jet->SetBranchAddress("tauPhi",&TauPhi_);
    tree_jet->SetBranchAddress("tauEta",&TauEta_);
    tree_jet->SetBranchAddress("tauIdDisc_DMfi",&DM_find);
    tree_jet->SetBranchAddress("tauIdDisc_byLooseDb3", &loose_iso);
    tree_jet->SetBranchAddress("tauIdDisc_byMediumDb3",&medium_iso);
    tree_jet->SetBranchAddress("tauIdDisc_byTightDb3",&tight_iso);
    tree_jet->SetBranchAddress("tauIdDisc_MVADBoldDM",&MVADB_find);
    tree_jet->SetBranchAddress("tauIdDisc_VLooseMVADBold",&MVAVloose_iso);
    tree_jet->SetBranchAddress("tauIdDisc_LooseMVADBold",&MVAloose_iso);
    tree_jet->SetBranchAddress("tauIdDisc_MediumMVADBold",&MVAmedium_iso);
    tree_jet->SetBranchAddress("tauIdDisc_TightMVADBold",&MVAtight_iso);
    tree_jet->SetBranchAddress("tauIdDisc_VTightMVADBold",&MVAVtight_iso);
    tree_jet->SetBranchAddress("tauIdDisc_VVTightMVADBold",&MVAVVtight_iso);
    tree_jet->SetBranchAddress("triggerResults", &triggerResults_);
    tree_jet->SetBranchAddress("tauIdDisc_agstEleTightMVA", &agstEleTightMVA_);
    tree_jet->SetBranchAddress("tauIdDisc_agstMuTight", &agstMuTight_);

    //Loop over entries
    Int_t nEntries = tree_jet->GetEntries();
    for(size_t iEvnt=0; iEvnt < nEntries; iEvnt++){
      tree_jet->GetEntry(iEvnt);

      //Calculate event weights
      float Weight = 1.0;
      if(!isData && pileupwght(h_p,iEvnt) != 0){
	Weight = xsections[ifile]*Lumi*pileupwght(h_p,iEvnt)/nEntries; 
      }

      // cout<<Weight<<endl;
      //Check trigger selection
      bool passPFJet500 = false;
      for(size_t i = 0; i < triggerResults_->size(); i++){
	if((*triggerResults_)[i].find("HLT_PFJet500_v") != std::string::npos)passPFJet500 = true;
      }
      if(!passPFJet500) continue;  //Firing Jet500 trigger

      //Find jet that fires the trigger
      int tx = 0; int index = -1;
      bool  jet_match = false;
      for(size_t ijet = 0; ijet < PF_jetPt->size(); ijet++){
	bool  jet_match = false;

	for(size_t tjet = 0; tjet < trig_jetPt->size(); tjet++){
	  if((*trig_jetPt)[tjet] > 500){
	    float dR_jet = deltaR((*PF_jetEta)[ijet],(*trig_jetEta)[tjet],(*PF_jetPhi)[ijet],(*trig_jetPhi)[tjet]);
	    if(dR_jet < 0.4){
	      jet_match = true;
	    }
	  }
	}
	if(jet_match){
	  tx++; //count the number of jet matchings
	  index = ijet; //store this index to exclude this jet
	}
      }

      if(index < 0) continue; //exclude event if no jet matched to trigger

      //Exclude jet that fires the trigger   
      //If more than one jet fire the trigger use all jets. 
      for(int ijet = 0; ijet < PF_jetPt->size(); ijet++){
	if(tx > 1 || ijet!=index){
	  if((*PF_jetPt)[ijet] < 100 || fabs((*PF_jetEta)[ijet]) > 2.3 ) continue; //require jet to be within |eta| < 2.3
	  
	  Denom_jetPt->Fill((*PF_jetPt)[ijet], Weight);
	  Denom_jetEta->Fill((*PF_jetEta)[ijet], Weight);
	  
	  float dRmax = 1.0; int itau = -1;
	  for(size_t kk = 0; kk < TauPt_->size(); kk++){
	    float dR_tau = deltaR((*PF_jetEta)[ijet],(*TauEta_)[kk],(*PF_jetPhi)[ijet],(*TauPhi_)[kk]);
	    if(dR_tau < dRmax){
	      dRmax = dR_tau;
	      itau = kk;
	    }
	  }
	  
	  if(dRmax < 0.3){
	    
	    if((*DM_find)[itau] > 0.5 && (*agstEleTightMVA_)[itau] > 0.5 && (*agstMuTight_)[itau] > 0.5){
	      NumeDM_Pt->Fill((*PF_jetPt)[ijet], Weight);
	      NumeDM_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      
	      if((*loose_iso)[itau] > 0.5){
		LNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		LNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*medium_iso)[itau] > 0.5){
		MNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		MNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*tight_iso)[itau] > 0.5){
		TNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		TNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*MVAVloose_iso)[itau] > 0.5){
		MVAVLNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		MVAVLNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*MVAloose_iso)[itau] > 0.5){
		MVALNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		MVALNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*MVAmedium_iso)[itau] > 0.5){
		MVAMNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		MVAMNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*MVAtight_iso)[itau] > 0.5){
		MVATNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		MVATNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }
	      if((*MVAVtight_iso)[itau] > 0.5){
		MVAVTNume_Pt->Fill((*PF_jetPt)[ijet], Weight);
		MVAVTNume_Eta->Fill((*PF_jetEta)[ijet], Weight);
	      }

	      
	    }// end if tau DM 
	  }//End tau loop
	}
	
      }//End pf jet loop
      
    }// End Event loop

    //delete all pointers
    delete trig_jetPt;
    delete trig_jetEta;
    delete trig_jetPhi;
    delete PF_jetEta;
    delete PF_jetPhi;
    delete PF_jetPt;
    delete TauPt_;
    delete TauEta_;
    delete TauPhi_;
    delete DM_find;
    delete loose_iso;
    delete medium_iso;
    delete tight_iso;
    delete MVADB_find;
    delete MVAVloose_iso;
    delete MVAloose_iso;
    delete MVAmedium_iso;
    delete MVAtight_iso;
    delete MVAVtight_iso;
    delete MVAVVtight_iso;
    delete triggerResults_;

    delete h_p;
    delete tree_jet;
    inFile->Close();
    delete inFile;

  }// End loop over files
  fout->cd();
  Denom_jetPt->Write();
  NumeDM_Pt->Write();
  LNume_Pt->Write();
  MNume_Pt->Write();
  TNume_Pt->Write();
  MVAVLNume_Pt->Write();
  MVALNume_Pt->Write();
  MVAMNume_Pt->Write();
  MVATNume_Pt->Write();
  MVAVTNume_Pt->Write();
  Denom_jetEta->Write();
  NumeDM_Eta->Write();
  LNume_Eta->Write();
  MNume_Eta->Write();
  TNume_Eta->Write();
  MVAVLNume_Eta->Write();
  MVALNume_Eta->Write();
  MVAMNume_Eta->Write();
  MVATNume_Eta->Write();
  MVAVTNume_Eta->Write();
  
  fout->Close();
}
