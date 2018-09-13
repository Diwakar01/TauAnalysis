void LegendSettings(TLegend *leg, int ncolumn){
  leg->SetNColumns(ncolumn);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.055);
  leg->SetTextFont(42);
}

void makePlots(TString denom, TString num, TString xTitle, TString label_string, TString outFileName)
{
  TFile fData("histograms_JetHT2017.root");
  TH1F* h_denom_data = (TH1F*) fData.Get(denom);
  TH1F* h_num_data   = (TH1F*) fData.Get(num);

  TFile fMC("histograms_QCD_MC2017.root");
  TH1F* h_denom_mc = (TH1F*) fMC.Get(denom);
  TH1F* h_num_mc   = (TH1F*) fMC.Get(num);

  TH1F* h_fr_data = (TH1F*)h_num_data->Clone();
  h_fr_data->Divide(h_denom_data);

  TH1F* h_fr_mc = (TH1F*)h_num_mc->Clone();
  h_fr_mc->Divide(h_denom_mc);

  TCanvas* c1 = new TCanvas();
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  gStyle->SetOptStat(0);

  TPad *pad1 = new TPad("pad1", "",0.05,0.22,0.96,1.0);
  pad1->Draw();
  pad1->cd();
  //pad1->Range(-25,-0.1375,225,1.2375);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetLogy();
  
  TH1F* histogram_base;
  if(denom.Contains("Pt")) histogram_base = new TH1F("histogram_base", "", 100, 0., 500.);
  else if(denom.Contains("Eta"))histogram_base = new TH1F("histogram_base", "", 50, -2.5, 2.5);
  histogram_base->SetTitle("");
  histogram_base->SetStats(false);
  histogram_base->SetMinimum(0.00001);
  histogram_base->SetMaximum(1.0);
  histogram_base->GetXaxis()->SetTitle(xTitle);
  histogram_base->GetYaxis()->SetTitle("jet #rightarrow #tau MisID probability");
  histogram_base->GetYaxis()->SetTitleSize(0.05);
  histogram_base->GetXaxis()->SetTitleSize(0.05);
  histogram_base->Draw("hist");


  h_fr_data->SetMarkerStyle(20);
  h_fr_data->SetMarkerColor(kBlack);
  h_fr_data->SetMarkerSize(1.0);
  h_fr_mc->SetMarkerStyle(24);
  h_fr_mc->SetMarkerColor(kBlack);
  h_fr_mc->SetMarkerSize(1.0);

  h_fr_data->Draw("same");
  h_fr_mc->Draw("same");

  TLegend *leg1 = new TLegend(0.6,0.7,0.9,0.8);
  LegendSettings(leg1, 1);
  leg1->AddEntry(h_fr_data, "Data", "lp");
  leg1->AddEntry(h_fr_mc, "MC", "lp");
  leg1->Draw();

  TPaveText* label = new TPaveText(0.2, 0.8, 0.55, 0.85, "brNDC");
  label->AddText(label_string);
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(0.05);
  label->SetTextFont(42);
  label->Draw();

  pad1->Modified();
  c1->cd();

  TPad* pad2 = new TPad("pad2", "",0.05,0.02,0.96,0.22);
  pad2->Draw();
  pad2->cd();
  //pad2->Range(-25,-0.2687085,225,0.2687085);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  pad2->SetFrameBorderMode(0);
  pad2->SetFrameBorderMode(0);
  pad2->SetGridy();

  TH1F* h_fr_ratio = (TH1F*)h_fr_data->Clone();
  h_fr_ratio->Divide(h_fr_mc);
  h_fr_ratio->GetYaxis()->SetTitle("Data/MC");
  h_fr_ratio->GetYaxis()->SetLabelSize(0.1);
  h_fr_ratio->GetYaxis()->SetTitleSize(0.15);
  h_fr_ratio->GetYaxis()->SetTitleOffset(0.3);
  h_fr_ratio->Draw();

  pad2->Modified();
  c1->cd();
  c1->Modified();
  c1->cd();

  c1->SaveAs("plots/FakeRate_"+outFileName+".png");
  //c1->SaveAs("plots/FakeRate_"+outFileName+".pdf");
  //c1->SaveAs("plots/FakeRate_"+outFileName+".jpg");
}

void makePlotsAll()
{

  makePlots("Denom_jetPt", "NumeDM_Pt", "Jet p_{T} (GeV)", "Decay Mode", "DecayMode_pt");
  makePlots("Denom_jetEta", "NumeDM_Eta", "Jet #eta", "Decay Mode", "DecayMode_eta");
  
  makePlots("Denom_jetPt", "LNume_Pt", "Jet p_{T} (GeV)", "Cut Loose", "Cut_Loose_pt");
  makePlots("Denom_jetEta", "LNume_Eta", "Jet #eta", "Cut Loose", "Cut_Loose_eta");

  makePlots("Denom_jetPt", "MNume_Pt", "Jet p_{T} (GeV)", "Cut Medium", "Cut_Medium_pt");
  makePlots("Denom_jetEta", "MNume_Eta", "Jet #eta", "Cut Medium", "Cut_Medium_eta");

  makePlots("Denom_jetPt", "TNume_Pt", "Jet p_{T} (GeV)", "Cut Tight", "Cut_Tight_pt");
  makePlots("Denom_jetEta", "TNume_Eta", "Jet #eta", "Cut Tight", "Cut_Tight_eta");

  makePlots("Denom_jetPt", "MVALNume_Pt", "Jet p_{T} (GeV)", "MVA Loose", "MVA_Loose_pt");
  makePlots("Denom_jetEta", "MVALNume_Eta", "Jet #eta", "MVA Loose", "MVA_Loose_eta");

  makePlots("Denom_jetPt", "MVAVLNume_Pt", "Jet p_{T} (GeV)", "MVA VLoose", "MVA_Loose_pt");
  makePlots("Denom_jetEta", "MVAVLNume_Eta", "Jet #eta", "MVA VLoose", "MVA_Loose_eta");

  makePlots("Denom_jetPt", "MVAMNume_Pt", "Jet p_{T} (GeV)", "MVA Medium", "MVA_Medium_pt");
  makePlots("Denom_jetEta", "MVAMNume_Eta", "Jet #eta", "MVA Medium", "MVA_Medium_eta");

  makePlots("Denom_jetPt", "MVATNume_Pt", "Jet p_{T} (GeV)", "MVA Tight", "MVA_Tight_pt");
  makePlots("Denom_jetEta", "MVATNume_Eta", "Jet #eta", "MVA Tight", "MVA_Tight_eta");

  makePlots("Denom_jetPt", "MVAVTNume_Pt", "Jet p_{T} (GeV)", "MVA VTight", "MVA_VTight_pt");
  makePlots("Denom_jetEta", "MVAVTNume_Eta", "Jet #eta", "MVA VTight", "MVA_VTight_eta");
}
