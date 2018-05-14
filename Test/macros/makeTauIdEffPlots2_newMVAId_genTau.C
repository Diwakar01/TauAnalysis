#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>

#include <vector>
#include <iostream>
#include <iomanip>

#include "assert.h"

using namespace std;

TH1* getHistogram(TFile* inputFile, const TString& dqmDirectory, const TString& meName)
{  
  TString histogramName = dqmDirectory;
  if ( !histogramName.EndsWith("/") ) histogramName.Append("/");
  histogramName.Append(meName);
  std::cout << "histogramName = " << histogramName.Data() << std::endl; 

  TH1* histogram = (TH1*)inputFile->Get(histogramName.Data());
  std::cout << "histogram = " << histogram << std::endl; 

  if ( !histogram->GetSumw2N()     ) histogram->Sumw2();
  if (  histogram->Integral() > 0. ) histogram->Scale(1./histogram->Integral());

  return histogram;
}

TGraphAsymmErrors* makeGraph(TH1* histogramNumerator, TH1* histogramDenominator)
{
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(histogramNumerator, histogramDenominator, "w");
  return graph;
}

struct plotEntryType
{
  plotEntryType(const TString& name)
    : name_(name),
      graphTauPtS_(0),
      graphTauPtL_(0),
      graphTauEta_(0),
      graphTauPhi_(0),
      graphTauLeadChHadPt_(0),
      graphTauLeadChHadEta_(0),
      graphNumVertices_(0),
      graphNumVerticesLbin_(0),
      graphUnbinned_(0)
  {
    float pTbin_[] = {0., 20., 40., 60., 80., 100., 120., 140., 160., 200., 250., 300., 350., 450., 600., 1000.};
    float pTbinS_[] = {0., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100.};
    //numeratorTauPtS_             = bookHistogram("numeratorTauPtS",             name, 20,       0.,         100.);
    numeratorTauPtS_             = bookHistogram("numeratorTauPtS",             name, 11, pTbinS_);
    //numeratorTauPtL_             = bookHistogram("numeratorTauPtL",             name, 20,       0.,        1000.);
    numeratorTauPtL_             = bookHistogram("numeratorTauPtL",             name, 15, pTbin_);
    //numeratorTauEta_             = bookHistogram("numeratorTauEta",             name, 46,     -2.3,         +2.3);
    numeratorTauEta_             = bookHistogram("numeratorTauEta",             name, 23,     -2.3,         +2.3);
    numeratorTauEtaImpactECAL_   = bookHistogram("numeratorTauEtaImpactECAL",   name, 460,     -2.3,         +2.3);
    numeratorTauPhi_             = bookHistogram("numeratorTauPhi",             name, 360, -TMath::Pi(), +TMath::Pi());
    numeratorTauLeadChHadPt_     = bookHistogram("numeratorTauLeadChHadPt",     name, 20,       0.,         100.);
    numeratorTauLeadChHadEta_    = bookHistogram("numeratorTauLeadChHadEta",    name, 460,     -2.3,         +2.3);
    numeratorTauNeuHadPt_        = bookHistogram("numeratorTauNeuHadPt",     name, 20,       0.,         100.);
    numeratorTauDecayMode_       = bookHistogram("numeratorTauDecayMode",       name, 15,      -0.5,        +14.5);
    numeratorCategoryRel_        = bookHistogram("numeratorCategoryRel",        name, 21,      -1.5,        +19.5);
    numeratorCategoryAbs_        = bookHistogram("numeratorCategoryAbs",        name, 21,      -1.5,        +19.5);
    numeratorNumVertices_        = bookHistogram("numeratorNumVertices",        name, 35,      -0.5,        +34.5);
    numeratorNumVerticesLbin_    = bookHistogram("numeratorNumVerticesLbin",    name, 3,      0.,        +30.);
    numeratorUnbinned_           = bookHistogram("numeratorUnbinned",           name,  1,      -0.5,         +0.5);
    //denominatorTauPtS_           = bookHistogram("denominatorTauPtS",           name, 20,       0.,         100.);
    denominatorTauPtS_           = bookHistogram("denominatorTauPtS",           name, 11, pTbinS_);
    //denominatorTauPtL_           = bookHistogram("denominatorTauPtL",           name, 20,       0.,        1000.);
    denominatorTauPtL_           = bookHistogram("denominatorTauPtL",           name, 15, pTbin_);
    //denominatorTauEta_           = bookHistogram("denominatorTauEta",           name, 46,     -2.3,         +2.3);
    denominatorTauEta_           = bookHistogram("denominatorTauEta",           name, 23,     -2.3,         +2.3);
    denominatorTauEtaImpactECAL_ = bookHistogram("denominatorTauEtaImpactECAL", name, 460,     -2.3,         +2.3);
    denominatorTauPhi_           = bookHistogram("denominatorTauPhi",           name, 360, -TMath::Pi(), +TMath::Pi());
    denominatorTauLeadChHadPt_   = bookHistogram("denominatorTauLeadChHadPt",   name, 20,       0.,         100.);
    denominatorTauLeadChHadEta_  = bookHistogram("denominatorTauLeadChHadEta",  name, 460,     -2.3,         +2.3);
    denominatorTauNeuHadPt_   = bookHistogram("denominatorTauNeuHadPt",   name, 20,       0.,         100.);
    denominatorTauDecayMode_     = bookHistogram("denominatorTauDecayMode",     name, 15,      -0.5,        +14.5);
    denominatorCategoryRel_      = bookHistogram("denominatorCategoryRel",      name, 21,      -1.5,        +19.5);
    denominatorCategoryAbs_      = bookHistogram("denominatorCategoryAbs",      name, 21,      -1.5,        +19.5);
    denominatorNumVertices_      = bookHistogram("denominatorNumVertices",      name, 35,      -0.5,        +34.5);
    denominatorNumVerticesLbin_  = bookHistogram("denominatorNumVerticesLbin",      name, 3,      0.,        +30.);
    denominatorUnbinned_         = bookHistogram("denominatorUnbinned",         name,  1,      -0.5,         +0.5);
  }
  ~plotEntryType()
  {
    delete numeratorTauPtS_;
    delete numeratorTauPtL_;
    delete numeratorTauEta_;
    delete numeratorTauEtaImpactECAL_;
    delete numeratorTauPhi_;
    delete numeratorTauLeadChHadPt_;
    delete numeratorTauLeadChHadEta_;
    delete numeratorTauNeuHadPt_;
    delete numeratorTauDecayMode_;
    delete numeratorCategoryRel_;
    delete numeratorCategoryAbs_;
    delete numeratorNumVertices_;
    delete numeratorNumVerticesLbin_;
    delete numeratorUnbinned_;
    delete denominatorTauPtS_;
    delete denominatorTauPtL_;
    delete denominatorTauEta_;
    delete denominatorTauEtaImpactECAL_;   
    delete denominatorTauPhi_;
    delete denominatorTauLeadChHadPt_;
    delete denominatorTauLeadChHadEta_;
    delete denominatorTauNeuHadPt_;
    delete denominatorTauDecayMode_;
    delete denominatorCategoryRel_;
    delete denominatorCategoryAbs_;
    delete denominatorNumVertices_;
    delete denominatorNumVerticesLbin_;
    delete denominatorUnbinned_;
    delete graphTauPtS_;
    delete graphTauPtL_;
    delete graphTauEta_;
    delete graphTauPhi_;
    delete graphNumVertices_;
    delete graphNumVerticesLbin_;
    delete graphCategoryRel_;
    delete graphCategoryAbs_;
    delete graphUnbinned_;
  }
  TH1* bookHistogram(const TString& name1, const TString& name2, Int_t numBinsX, Float_t xMin, Float_t xMax)
  {
    TString histogramName = Form("%s_%s", name1.Data(), name2.Data());
    TH1* histogram = new TH1F(histogramName.Data(), histogramName.Data(), numBinsX, xMin, xMax);
    return histogram;
  }
  TH1* bookHistogram(const TString& name1, const TString& name2, Int_t numBinsX, Float_t* xBins)
  {
    TString histogramName = Form("%s_%s", name1.Data(), name2.Data());
    TH1* histogram = new TH1F(histogramName.Data(), histogramName.Data(), numBinsX, xBins);
    return histogram;
  }
  void fillHistograms(bool passesNumerator, bool passesDenominator, 
		      Float_t tauPt, Float_t tauEta, Float_t tauEtaImpactECAL, Float_t tauPhi, Float_t tauleadChHadPt, Float_t tauleadChHadEta, Float_t tauNeutralHadPt,
		      Int_t tauDecayMode, Int_t category, Float_t numVertices,
		      Float_t evtWeight)
  {
    if ( passesDenominator ) {
      denominatorTauPtS_->Fill(tauPt, evtWeight);
      denominatorTauPtL_->Fill(tauPt, evtWeight);
      denominatorTauEta_->Fill(tauEta, evtWeight);
      denominatorTauEtaImpactECAL_->Fill(tauEtaImpactECAL, evtWeight);
      denominatorTauPhi_->Fill(tauPhi, evtWeight);
      denominatorTauLeadChHadPt_->Fill(tauleadChHadPt, evtWeight);
      denominatorTauLeadChHadEta_->Fill(tauleadChHadEta, evtWeight);
      denominatorTauNeuHadPt_->Fill(tauNeutralHadPt, evtWeight);
      denominatorTauDecayMode_->Fill(tauDecayMode, evtWeight);
      denominatorCategoryRel_->Fill(category, evtWeight);
      for ( int i = -1; i < 20; ++i ) {
	denominatorCategoryAbs_->Fill(i, evtWeight);
      }
      denominatorNumVertices_->Fill(numVertices, evtWeight);
      denominatorNumVerticesLbin_->Fill(numVertices, evtWeight);
      denominatorUnbinned_->Fill(0., evtWeight);
      if ( passesNumerator ) {
	numeratorTauPtS_->Fill(tauPt, evtWeight);
	numeratorTauPtL_->Fill(tauPt, evtWeight);
	numeratorTauEta_->Fill(tauEta, evtWeight);
	numeratorTauEtaImpactECAL_->Fill(tauEtaImpactECAL, evtWeight);
	numeratorTauPhi_->Fill(tauPhi, evtWeight);
	numeratorTauLeadChHadPt_->Fill(tauleadChHadPt, evtWeight);
	numeratorTauLeadChHadEta_->Fill(tauleadChHadEta, evtWeight);
	numeratorTauNeuHadPt_->Fill(tauNeutralHadPt, evtWeight);
	numeratorTauDecayMode_->Fill(tauDecayMode, evtWeight);
	numeratorCategoryRel_->Fill(category, evtWeight);
	numeratorCategoryAbs_->Fill(category, evtWeight);
	numeratorNumVertices_->Fill(numVertices, evtWeight);
	numeratorNumVerticesLbin_->Fill(numVertices, evtWeight);
	numeratorUnbinned_->Fill(0., evtWeight);
      }
    }
  }
  void makeGraphs()
  {    
    graphTauPtS_           = makeGraph(numeratorTauPtS_,           denominatorTauPtS_);
    graphTauPtL_           = makeGraph(numeratorTauPtL_,           denominatorTauPtL_);
    graphTauEta_           = makeGraph(numeratorTauEta_,           denominatorTauEta_);
    graphTauEtaImpactECAL_ = makeGraph(numeratorTauEtaImpactECAL_, denominatorTauEtaImpactECAL_);
    graphTauPhi_           = makeGraph(numeratorTauPhi_,           denominatorTauPhi_);
    graphTauLeadChHadPt_   = makeGraph(numeratorTauLeadChHadPt_,   denominatorTauLeadChHadPt_);
    graphTauLeadChHadEta_  = makeGraph(numeratorTauLeadChHadEta_,   denominatorTauLeadChHadEta_);
    graphTauNeuHadPt_      = makeGraph(numeratorTauNeuHadPt_,   denominatorTauNeuHadPt_);
    graphTauDecayMode_     = makeGraph(numeratorTauDecayMode_,     denominatorTauDecayMode_);
    graphCategoryRel_      = makeGraph(numeratorCategoryRel_,      denominatorCategoryRel_);
    graphCategoryAbs_      = makeGraph(numeratorCategoryAbs_,      denominatorCategoryAbs_);
    graphNumVertices_      = makeGraph(numeratorNumVertices_,      denominatorNumVertices_);    
    graphNumVerticesLbin_  = makeGraph(numeratorNumVerticesLbin_,      denominatorNumVerticesLbin_);
    graphUnbinned_         = makeGraph(numeratorUnbinned_,         denominatorUnbinned_); 
    std::cout << "<makeGraphs>:" << std::endl;
    Double_t x, y;
    graphUnbinned_->GetPoint(0, x, y);
    std::cout << " name = " << name_.Data() << ": numerator = " << numeratorUnbinned_->Integral() << "," 
	      << " denominator = " << denominatorUnbinned_->Integral() << " --> efficiency/fake-rate = " << y << std::endl;
  }
  TString name_;
  TH1* numeratorTauPtS_;
  TH1* numeratorTauPtL_;
  TH1* numeratorTauEta_;
  TH1* numeratorTauEtaImpactECAL_;
  TH1* numeratorTauPhi_;
  TH1* numeratorTauLeadChHadPt_;
  TH1* numeratorTauLeadChHadEta_;
  TH1* numeratorTauNeuHadPt_;
  TH1* numeratorTauDecayMode_;
  TH1* numeratorCategoryRel_;
  TH1* numeratorCategoryAbs_;
  TH1* numeratorNumVertices_;
  TH1* numeratorNumVerticesLbin_;
  TH1* numeratorUnbinned_;
  TH1* denominatorTauPtS_;
  TH1* denominatorTauPtL_;
  TH1* denominatorTauEta_;
  TH1* denominatorTauEtaImpactECAL_;
  TH1* denominatorTauPhi_;
  TH1* denominatorTauLeadChHadPt_;
  TH1* denominatorTauLeadChHadEta_;
  TH1* denominatorTauNeuHadPt_;
  TH1* denominatorTauDecayMode_;
  TH1* denominatorCategoryRel_;
  TH1* denominatorCategoryAbs_;
  TH1* denominatorNumVertices_;
  TH1* denominatorNumVerticesLbin_;
  TH1* denominatorUnbinned_;
  TGraphAsymmErrors* graphTauPtS_;
  TGraphAsymmErrors* graphTauPtL_;
  TGraphAsymmErrors* graphTauEta_;
  TGraphAsymmErrors* graphTauEtaImpactECAL_;
  TGraphAsymmErrors* graphTauPhi_;
  TGraphAsymmErrors* graphTauLeadChHadPt_;
  TGraphAsymmErrors* graphTauLeadChHadEta_;
  TGraphAsymmErrors* graphTauNeuHadPt_;
  TGraphAsymmErrors* graphTauDecayMode_;
  TGraphAsymmErrors* graphCategoryRel_;
  TGraphAsymmErrors* graphCategoryAbs_;
  TGraphAsymmErrors* graphNumVertices_;
  TGraphAsymmErrors* graphNumVerticesLbin_;
  TGraphAsymmErrors* graphUnbinned_;
};

void showGraphs(double canvasSizeX, double canvasSizeY,
                TGraph* graph1, const std::string& legendEntry1,
                TGraph* graph2, const std::string& legendEntry2,
                TGraph* graph3, const std::string& legendEntry3,
                TGraph* graph4, const std::string& legendEntry4,
		TGraph* graph5, const std::string& legendEntry5,
		TGraph* graph6, const std::string& legendEntry6,
                double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
                std::vector<std::string>& labelTextLines, double labelTextSize,
                double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
                double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
                bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                const std::string& outputFileName)
{
  //std::cout << "<showGraphs>:" << std::endl;
  //std::cout << " graph1 = " << graph1;
  //if ( graph1 ) std::cout << " (" << graph1->GetName() << ")";
  //std::cout << std::endl;
  //std::cout << " graph2 = " << graph2;
  //if ( graph2 ) std::cout << " (" << graph2->GetName() << ")";
  //std::cout << std::endl;
  assert(graph1 && graph2);

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  canvas->SetLogy(useLogScale);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 100, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  int colors[6] = { 1, 2, 4, 6, 8, 9 };
  int markerStyles[6] = { 20, 21, 22, 23, 33, 24 };

  graph1->SetLineColor(colors[0]);
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->Draw("p");

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->Draw("p");
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->Draw("p");
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->Draw("p");
  }
  
  if ( graph5 ) {
    graph5->SetLineColor(colors[4]);
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->Draw("p");
  }

  if ( graph6 ) {
    graph6->SetLineColor(colors[5]);
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->Draw("p");
  }

  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);
  legend->AddEntry(graph1, legendEntry1.data(), "pl");
  if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "pl");
  if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "pl");
  if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "pl");
  if ( graph5 ) legend->AddEntry(graph5, legendEntry5.data(), "pl");
  if ( graph6 ) legend->AddEntry(graph6, legendEntry6.data(), "pl");
  legend->Draw();

  TPaveText* label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
  for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
        labelTextLine != labelTextLines.end(); ++labelTextLine ) {
    label->AddText(labelTextLine->data());
  }
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(labelTextSize);
  label->Draw();

  TPaveText* label_cms = new TPaveText(0.14, 0.905, 0.88, 0.965, "brNDC");
  label_cms->AddText("CMS Simulation 2015, #sqrt{s} = 13 TeV");
  label_cms->SetFillColor(10);
  label_cms->SetBorderSize(0);
  label_cms->SetTextColor(1);
  label_cms->SetTextAlign(12);
  label_cms->SetTextSize(0.045);
  label_cms->Draw();

  canvas->Update();
  std::string outputFileName_plot = "plots_Sep2016_genJet_v4/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".eps").data());
  canvas->Print(std::string(outputFileName_plot).append(".C").data());
  
  delete dummyHistogram;
  delete label;
  delete legend;
  delete canvas;  
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		plotEntryType* plots1, const std::string& legendEntry1,
		plotEntryType* plots2, const std::string& legendEntry2,
		plotEntryType* plots3, const std::string& legendEntry3,
		plotEntryType* plots4, const std::string& legendEntry4,
		plotEntryType* plots5, const std::string& legendEntry5,
		plotEntryType* plots6, const std::string& legendEntry6,
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
                const std::string& labelTextLine, double labelTextSize,
                double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
                bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                const std::string& outputFileName)
{
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_base = std::string(outputFileName, 0, idx);

  std::vector<std::string> labelTextLines;  
  labelTextLines.push_back(labelTextLine);

  std::string outputFileName_TauPtS = std::string(outputFileName_base).append("_PtS");
  TGraph* graph1_TauPtS = ( plots1 ) ? plots1->graphTauPtS_ : 0;
  TGraph* graph2_TauPtS = ( plots2 ) ? plots2->graphTauPtS_ : 0;
  TGraph* graph3_TauPtS = ( plots3 ) ? plots3->graphTauPtS_ : 0;
  TGraph* graph4_TauPtS = ( plots4 ) ? plots4->graphTauPtS_ : 0;
  TGraph* graph5_TauPtS = ( plots5 ) ? plots5->graphTauPtS_ : 0;
  TGraph* graph6_TauPtS = ( plots6 ) ? plots6->graphTauPtS_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauPtS, legendEntry1,
             graph2_TauPtS, legendEntry2,
             graph3_TauPtS, legendEntry3,
             graph4_TauPtS, legendEntry4,
	     graph5_TauPtS, legendEntry5,
	     graph6_TauPtS, legendEntry6,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             0., 100., "P_{T}^{#tau_{h}} / GeV", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauPtS);
  std::string outputFileName_TauPtL = std::string(outputFileName_base).append("_PtL");
  TGraph* graph1_TauPtL = ( plots1 ) ? plots1->graphTauPtL_ : 0;
  TGraph* graph2_TauPtL = ( plots2 ) ? plots2->graphTauPtL_ : 0;
  TGraph* graph3_TauPtL = ( plots3 ) ? plots3->graphTauPtL_ : 0;
  TGraph* graph4_TauPtL = ( plots4 ) ? plots4->graphTauPtL_ : 0;
  TGraph* graph5_TauPtL = ( plots5 ) ? plots5->graphTauPtL_ : 0;
  TGraph* graph6_TauPtL = ( plots6 ) ? plots6->graphTauPtL_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauPtL, legendEntry1,
             graph2_TauPtL, legendEntry2,
             graph3_TauPtL, legendEntry3,
             graph4_TauPtL, legendEntry4,
	     graph5_TauPtL, legendEntry5,
	     graph6_TauPtL, legendEntry6,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             0., 1000., "P_{T}^{#tau_{h}} / GeV", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauPtL);

  std::string outputFileName_TauEta = std::string(outputFileName_base).append("_Eta");
  TGraph* graph1_TauEta = ( plots1 ) ? plots1->graphTauEta_ : 0;
  TGraph* graph2_TauEta = ( plots2 ) ? plots2->graphTauEta_ : 0;
  TGraph* graph3_TauEta = ( plots3 ) ? plots3->graphTauEta_ : 0;
  TGraph* graph4_TauEta = ( plots4 ) ? plots4->graphTauEta_ : 0;
  TGraph* graph5_TauEta = ( plots5 ) ? plots5->graphTauEta_ : 0;
  TGraph* graph6_TauEta = ( plots6 ) ? plots6->graphTauEta_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauEta, legendEntry1,
             graph2_TauEta, legendEntry2,
             graph3_TauEta, legendEntry3,
             graph4_TauEta, legendEntry4,
	     graph5_TauEta, legendEntry5,
	     graph6_TauEta, legendEntry6,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -2.3, +2.3, "#eta_{#tau_{h}}", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauEta);
  /*
  std::string outputFileName_TauEtaImpactECAL = std::string(outputFileName_base).append("_EtaImpactECAL");
  TGraph* graph1_TauEtaImpactECAL = ( plots1 ) ? plots1->graphTauEtaImpactECAL_ : 0;
  TGraph* graph2_TauEtaImpactECAL = ( plots2 ) ? plots2->graphTauEtaImpactECAL_ : 0;
  TGraph* graph3_TauEtaImpactECAL = ( plots3 ) ? plots3->graphTauEtaImpactECAL_ : 0;
  TGraph* graph4_TauEtaImpactECAL = ( plots4 ) ? plots4->graphTauEtaImpactECAL_ : 0;
  TGraph* graph5_TauEtaImpactECAL = ( plots5 ) ? plots5->graphTauEtaImpactECAL_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauEtaImpactECAL, legendEntry1,
             graph2_TauEtaImpactECAL, legendEntry2,
             graph3_TauEtaImpactECAL, legendEntry3,
             graph4_TauEtaImpactECAL, legendEntry4,
	     graph5_TauEtaImpactECAL, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -2.3, +2.3, "#eta_{impact}^{ECAL}", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauEtaImpactECAL);

  std::string outputFileName_TauPhi = std::string(outputFileName_base).append("_Phi");
  TGraph* graph1_TauPhi = ( plots1 ) ? plots1->graphTauPhi_ : 0;
  TGraph* graph2_TauPhi = ( plots2 ) ? plots2->graphTauPhi_ : 0;
  TGraph* graph3_TauPhi = ( plots3 ) ? plots3->graphTauPhi_ : 0;
  TGraph* graph4_TauPhi = ( plots4 ) ? plots4->graphTauPhi_ : 0;
  TGraph* graph5_TauPhi = ( plots5 ) ? plots5->graphTauPhi_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauPhi, legendEntry1,
             graph2_TauPhi, legendEntry2,
             graph3_TauPhi, legendEntry3,
             graph4_TauPhi, legendEntry4,
	     graph5_TauPhi, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -TMath::Pi(), +TMath::Pi(), "#phi_{#tauJet} / Rad", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauPhi);
  */
  //std::string outputFileName_TauLeadChHadPt = std::string(outputFileName_base).append("_LeadChHadPt");
  //TGraph* graph1_TauLeadChHadPt = ( plots1 ) ? plots1->graphTauLeadChHadPt_ : 0;
  //TGraph* graph2_TauLeadChHadPt = ( plots2 ) ? plots2->graphTauLeadChHadPt_ : 0;
  //TGraph* graph3_TauLeadChHadPt = ( plots3 ) ? plots3->graphTauLeadChHadPt_ : 0;
  //TGraph* graph4_TauLeadChHadPt = ( plots4 ) ? plots4->graphTauLeadChHadPt_ : 0;
  //TGraph* graph5_TauLeadChHadPt = ( plots5 ) ? plots5->graphTauLeadChHadPt_ : 0;
  //showGraphs(canvasSizeX, canvasSizeY,
  //           graph1_TauLeadChHadPt, legendEntry1,
  //           graph2_TauLeadChHadPt, legendEntry2,
  //           graph3_TauLeadChHadPt, legendEntry3,
  //           graph4_TauLeadChHadPt, legendEntry4,
  //           graph5_TauLeadChHadPt, legendEntry5,
  //           legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY,
  //           labelTextLines, labelTextSize,
  //           labelPosX, labelPosY, labelSizeX, labelSizeY,
  //           0., 100., "P_{T}^{lead.Ch.Had} / GeV", 1.3,
  //           useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
  //           outputFileName_TauLeadChHadPt);
  //
  //std::string outputFileName_TauLeadChHadEta = std::string(outputFileName_base).append("_LeadChHadEta");
  //TGraph* graph1_TauLeadChHadEta = ( plots1 ) ? plots1->graphTauLeadChHadEta_ : 0;
  //TGraph* graph2_TauLeadChHadEta = ( plots2 ) ? plots2->graphTauLeadChHadEta_ : 0;
  //TGraph* graph3_TauLeadChHadEta = ( plots3 ) ? plots3->graphTauLeadChHadEta_ : 0;
  //TGraph* graph4_TauLeadChHadEta = ( plots4 ) ? plots4->graphTauLeadChHadEta_ : 0;
  //TGraph* graph5_TauLeadChHadEta = ( plots5 ) ? plots5->graphTauLeadChHadEta_ : 0;
  //showGraphs(canvasSizeX, canvasSizeY,
  //           graph1_TauLeadChHadEta, legendEntry1,
  //           graph2_TauLeadChHadEta, legendEntry2,
  //           graph3_TauLeadChHadEta, legendEntry3,
  //           graph4_TauLeadChHadEta, legendEntry4,
  //           graph5_TauLeadChHadEta, legendEntry5,
  //           legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY,
  //           labelTextLines, labelTextSize,
  //           labelPosX, labelPosY, labelSizeX, labelSizeY,
  //           -2.3, +2.3, "#eta_{lead.Ch.Had}", 1.3,
  //           useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
  //           outputFileName_TauLeadChHadEta);
  //
  //std::string outputFileName_TauNeuHadPt = std::string(outputFileName_base).append("_NeutralHadPt");
  //TGraph* graph1_TauNeuHadPt = ( plots1 ) ? plots1->graphTauNeuHadPt_ : 0;
  //TGraph* graph2_TauNeuHadPt = ( plots2 ) ? plots2->graphTauNeuHadPt_ : 0;
  //TGraph* graph3_TauNeuHadPt = ( plots3 ) ? plots3->graphTauNeuHadPt_ : 0;
  //TGraph* graph4_TauNeuHadPt = ( plots4 ) ? plots4->graphTauNeuHadPt_ : 0;
  //TGraph* graph5_TauNeuHadPt = ( plots5 ) ? plots5->graphTauNeuHadPt_ : 0;
  //showGraphs(canvasSizeX, canvasSizeY,
  //           graph1_TauNeuHadPt, legendEntry1,
  //           graph2_TauNeuHadPt, legendEntry2,
  //           graph3_TauNeuHadPt, legendEntry3,
  //           graph4_TauNeuHadPt, legendEntry4,
  //           graph5_TauNeuHadPt, legendEntry5,
  //           legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY,
  //           labelTextLines, labelTextSize,
  //           labelPosX, labelPosY, labelSizeX, labelSizeY,
  //           0., 100., "P_{T}^{Neutral Had} / GeV", 1.3,
  //           useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
  //           outputFileName_TauNeuHadPt);
  /*
  std::string outputFileName_TauDecayMode = std::string(outputFileName_base).append("_DecayMode");
  TGraph* graph1_TauDecayMode = ( plots1 ) ? plots1->graphTauDecayMode_ : 0;
  TGraph* graph2_TauDecayMode = ( plots2 ) ? plots2->graphTauDecayMode_ : 0;
  TGraph* graph3_TauDecayMode = ( plots3 ) ? plots3->graphTauDecayMode_ : 0;
  TGraph* graph4_TauDecayMode = ( plots4 ) ? plots4->graphTauDecayMode_ : 0;
  TGraph* graph5_TauDecayMode = ( plots5 ) ? plots5->graphTauDecayMode_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_TauDecayMode, legendEntry1,
             graph2_TauDecayMode, legendEntry2,
             graph3_TauDecayMode, legendEntry3,
             graph4_TauDecayMode, legendEntry4,
	     graph5_TauDecayMode, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -0.5, +14.5, "Decay mode", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_TauDecayMode);

  std::string outputFileName_CategoryRel = std::string(outputFileName_base).append("_CategoryRel");
  TGraph* graph1_CategoryRel = ( plots1 ) ? plots1->graphCategoryRel_ : 0;
  TGraph* graph2_CategoryRel = ( plots2 ) ? plots2->graphCategoryRel_ : 0;
  TGraph* graph3_CategoryRel = ( plots3 ) ? plots3->graphCategoryRel_ : 0;
  TGraph* graph4_CategoryRel = ( plots4 ) ? plots4->graphCategoryRel_ : 0;
  TGraph* graph5_CategoryRel = ( plots5 ) ? plots5->graphCategoryRel_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_CategoryRel, legendEntry1,
             graph2_CategoryRel, legendEntry2,
             graph3_CategoryRel, legendEntry3,
             graph4_CategoryRel, legendEntry4,
	     graph5_CategoryRel, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -1.5, +19.5, "Category", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_CategoryRel);
  std::string outputFileName_CategoryAbs = std::string(outputFileName_base).append("_CategoryAbs");
  TGraph* graph1_CategoryAbs = ( plots1 ) ? plots1->graphCategoryAbs_ : 0;
  TGraph* graph2_CategoryAbs = ( plots2 ) ? plots2->graphCategoryAbs_ : 0;
  TGraph* graph3_CategoryAbs = ( plots3 ) ? plots3->graphCategoryAbs_ : 0;
  TGraph* graph4_CategoryAbs = ( plots4 ) ? plots4->graphCategoryAbs_ : 0;
  TGraph* graph5_CategoryAbs = ( plots5 ) ? plots5->graphCategoryAbs_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_CategoryAbs, legendEntry1,
             graph2_CategoryAbs, legendEntry2,
             graph3_CategoryAbs, legendEntry3,
             graph4_CategoryAbs, legendEntry4,
	     graph5_CategoryAbs, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -1.5, +19.5, "Category", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_CategoryAbs);
  */
  std::string outputFileName_NumVertices = std::string(outputFileName_base).append("_NumVertices");
  TGraph* graph1_NumVertices = ( plots1 ) ? plots1->graphNumVertices_ : 0;
  TGraph* graph2_NumVertices = ( plots2 ) ? plots2->graphNumVertices_ : 0;
  TGraph* graph3_NumVertices = ( plots3 ) ? plots3->graphNumVertices_ : 0;
  TGraph* graph4_NumVertices = ( plots4 ) ? plots4->graphNumVertices_ : 0;
  TGraph* graph5_NumVertices = ( plots5 ) ? plots5->graphNumVertices_ : 0;
  TGraph* graph6_NumVertices = ( plots6 ) ? plots6->graphNumVertices_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_NumVertices, legendEntry1,
             graph2_NumVertices, legendEntry2,
             graph3_NumVertices, legendEntry3,
             graph4_NumVertices, legendEntry4,
	     graph5_NumVertices, legendEntry5,
	     graph6_NumVertices, legendEntry6,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -0.5, 34.5, "Num. reconstructed Vertices", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_NumVertices);
  /*
  std::string outputFileName_NumVerticesLbin = std::string(outputFileName_base).append("_NumVerticesLbin");
  TGraph* graph1_NumVerticesLbin = ( plots1 ) ? plots1->graphNumVerticesLbin_ : 0;
  TGraph* graph2_NumVerticesLbin = ( plots2 ) ? plots2->graphNumVerticesLbin_ : 0;
  TGraph* graph3_NumVerticesLbin = ( plots3 ) ? plots3->graphNumVerticesLbin_ : 0;
  TGraph* graph4_NumVerticesLbin = ( plots4 ) ? plots4->graphNumVerticesLbin_ : 0;
  TGraph* graph5_NumVerticesLbin = ( plots5 ) ? plots5->graphNumVerticesLbin_ : 0;
  showGraphs(canvasSizeX, canvasSizeY,
             graph1_NumVerticesLbin, legendEntry1,
             graph2_NumVerticesLbin, legendEntry2,
             graph3_NumVerticesLbin, legendEntry3,
             graph4_NumVerticesLbin, legendEntry4,
             graph5_NumVerticesLbin, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY,
             labelTextLines, labelTextSize,
             labelPosX, labelPosY, labelSizeX, labelSizeY,
             -0.5, 34.5, "Num. reconstructed Vertices", 1.3,
             useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_NumVerticesLbin);
  */
}

double square(double x)
{
  return x*x;
}

TGraphAsymmErrors* compRatioGraph(const std::string& ratioGraphName, const TGraph* numerator, const TGraph* denominator)
{
  assert(numerator->GetN() == denominator->GetN());
  int nPoints = numerator->GetN();

  TGraphAsymmErrors* graphRatio = new TGraphAsymmErrors(nPoints);
  graphRatio->SetName(ratioGraphName.data());

  for ( int iPoint = 0; iPoint < nPoints; ++iPoint ){
    double x_numerator, y_numerator;
    numerator->GetPoint(iPoint, x_numerator, y_numerator);
    double xErrUp_numerator = 0.;
    double xErrDown_numerator = 0.;
    double yErrUp_numerator = 0.;
    double yErrDown_numerator = 0.;
    if ( dynamic_cast<const TGraphAsymmErrors*>(numerator) ) {
      const TGraphAsymmErrors* numerator_asymmerrors = dynamic_cast<const TGraphAsymmErrors*>(numerator);
      xErrUp_numerator = numerator_asymmerrors->GetErrorXhigh(iPoint);
      xErrDown_numerator = numerator_asymmerrors->GetErrorXlow(iPoint);
      yErrUp_numerator = numerator_asymmerrors->GetErrorYhigh(iPoint);
      yErrDown_numerator = numerator_asymmerrors->GetErrorYlow(iPoint);
    } else if ( dynamic_cast<const TGraphErrors*>(numerator) ) {
      const TGraphErrors* numerator_errors = dynamic_cast<const TGraphErrors*>(numerator);
      xErrUp_numerator = numerator_errors->GetErrorX(iPoint);
      xErrDown_numerator = xErrUp_numerator;
      yErrUp_numerator = numerator_errors->GetErrorY(iPoint);
      yErrDown_numerator = yErrUp_numerator;
    }

    double x_denominator, y_denominator;
    denominator->GetPoint(iPoint, x_denominator, y_denominator);
    assert(x_denominator == x_numerator);
    double xErrUp_denominator = 0.;
    double xErrDown_denominator = 0.;
    double yErrUp_denominator = 0.;
    double yErrDown_denominator = 0.;
    if ( dynamic_cast<const TGraphAsymmErrors*>(denominator) ) {
      const TGraphAsymmErrors* denominator_asymmerrors = dynamic_cast<const TGraphAsymmErrors*>(denominator);
      xErrUp_denominator = denominator_asymmerrors->GetErrorXhigh(iPoint);
      xErrDown_denominator = denominator_asymmerrors->GetErrorXlow(iPoint);
      yErrUp_denominator = denominator_asymmerrors->GetErrorYhigh(iPoint);
      yErrDown_denominator = denominator_asymmerrors->GetErrorYlow(iPoint);
    } else if ( dynamic_cast<const TGraphErrors*>(denominator) ) {
      const TGraphErrors* denominator_errors = dynamic_cast<const TGraphErrors*>(denominator);
      xErrUp_denominator = denominator_errors->GetErrorX(iPoint);
      xErrDown_denominator = xErrUp_denominator;
      yErrUp_denominator = denominator_errors->GetErrorY(iPoint);
      yErrDown_denominator = yErrUp_denominator;
    }

    double x_ratio = x_numerator;
    double y_ratio = ( y_denominator > 0. ) ? (y_numerator/y_denominator) : 0.;
    double xErrUp_ratio = TMath::Max(xErrUp_numerator, xErrUp_denominator);
    double xErrDown_ratio = TMath::Max(xErrDown_numerator, xErrDown_denominator);
    double yErr2Up_ratio = 0.;
    if ( y_numerator   ) yErr2Up_ratio += square(yErrUp_numerator/y_numerator);
    if ( y_denominator ) yErr2Up_ratio += square(yErrDown_denominator/y_numerator);
    double yErrUp_ratio = TMath::Sqrt(yErr2Up_ratio)*y_ratio;
    double yErr2Down_ratio = 0.;
    if ( y_numerator   ) yErr2Down_ratio += square(yErrDown_numerator/y_numerator);
    if ( y_denominator ) yErr2Down_ratio += square(yErrUp_denominator/y_numerator);
    double yErrDown_ratio = TMath::Sqrt(yErr2Down_ratio)*y_ratio;

    graphRatio->SetPoint(iPoint, x_ratio, y_ratio - 1.);
    graphRatio->SetPointError(iPoint, xErrDown_ratio, xErrUp_ratio, yErrDown_ratio, yErrUp_ratio);
  }
  
  graphRatio->SetLineColor(numerator->GetLineColor());
  graphRatio->SetLineWidth(numerator->GetLineWidth());
  graphRatio->SetMarkerColor(numerator->GetMarkerColor());
  graphRatio->SetMarkerStyle(numerator->GetMarkerStyle());
  graphRatio->SetMarkerSize(numerator->GetMarkerSize());

  return graphRatio;
}
/*
void showGraphs_SoverB(double canvasSizeX, double canvasSizeY,
		       TGraph* graph1_S, TGraph* graph1_B, const std::string& legendEntry1,
		       TGraph* graph2_S, TGraph* graph2_B, const std::string& legendEntry2,
		       TGraph* graph3_S, TGraph* graph3_B, const std::string& legendEntry3,
		       TGraph* graph4_S, TGraph* graph4_B, const std::string& legendEntry4,
		       TGraph* graph5_S, TGraph* graph5_B, const std::string& legendEntry5,
		       double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		       const std::string& labelTextLine, double labelTextSize,
		       double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		       bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		       const std::string& outputFileName)
{
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_base = std::string(outputFileName, 0, idx);

  std::vector<std::string> labelTextLines;  
  labelTextLines.push_back(labelTextLine);

  std::string outputFileName_Category = std::string(outputFileName_base).append("_Category");
  TGraph* graph1_SoverB = ( graph1_S && graph1_B ) ? compRatioGraph(Form("%s_div_%s", graph1_S->GetName(), graph1_B->GetName()), graph1_S, graph1_B) : 0;
  TGraph* graph2_SoverB = ( graph2_S && graph2_B ) ? compRatioGraph(Form("%s_div_%s", graph2_S->GetName(), graph2_B->GetName()), graph2_S, graph2_B) : 0;
  TGraph* graph3_SoverB = ( graph3_S && graph3_B ) ? compRatioGraph(Form("%s_div_%s", graph3_S->GetName(), graph3_B->GetName()), graph3_S, graph3_B) : 0;
  TGraph* graph4_SoverB = ( graph4_S && graph4_B ) ? compRatioGraph(Form("%s_div_%s", graph4_S->GetName(), graph4_B->GetName()), graph4_S, graph4_B) : 0;
  TGraph* graph5_SoverB = ( graph5_S && graph5_B ) ? compRatioGraph(Form("%s_div_%s", graph5_S->GetName(), graph5_B->GetName()), graph5_S, graph5_B) : 0;
  showGraphs(canvasSizeX, canvasSizeY,
	     graph1_SoverB, legendEntry1,
             graph2_SoverB, legendEntry2,
             graph3_SoverB, legendEntry3,
             graph4_SoverB, legendEntry4,
	     graph5_SoverB, legendEntry5,
             legendTextSize, legendPosX, legendPosY, legendSizeX, legendSizeY, 
             labelTextLines, labelTextSize,
	     labelPosX, labelPosY, labelSizeX, labelSizeY,
             -1.5, +19.5, "Category", 1.3,
	     useLogScale, yMin, yMax, yAxisTitle, yAxisOffset,
             outputFileName_Category);
}
*/
TH1* rebinHistogram(const TH1* histogram, unsigned numBinsMin_rebinned, double xMin, double xMax, bool normalize)
{
  TH1* histogram_rebinned = 0;

  if ( histogram ) {
    unsigned numBins = histogram->GetNbinsX();
    unsigned numBins_withinRange = 0;
    for ( unsigned iBin = 1; iBin <= numBins; ++iBin ) {
      double binCenter = histogram->GetBinCenter(iBin);
      if ( binCenter >= xMin && binCenter <= xMax ) ++numBins_withinRange;
    }

    std::cout << "histogram = " << histogram->GetName() << ":" 
              << " numBins(" << xMin << ".." << "xMax) = " << numBins_withinRange << ", integral = " << histogram->Integral() << std::endl;
    
    unsigned numBins_rebinned = numBins_withinRange;

    for ( int combineNumBins = 5; combineNumBins >= 2; --combineNumBins ) {
      if ( numBins_withinRange >= (combineNumBins*numBinsMin_rebinned) && (numBins % combineNumBins) == 0 ) {
        numBins_rebinned /= combineNumBins;
        numBins_withinRange /= combineNumBins;
      }
    }

    std::string histogramName_rebinned = std::string(histogram->GetName()).append("_rebinned");
    histogram_rebinned = new TH1D(histogramName_rebinned.data(), histogram->GetTitle(), numBins_rebinned, xMin, xMax);
    if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();

    TAxis* xAxis = histogram_rebinned->GetXaxis();
      
    unsigned iBin = 1;
    for ( unsigned iBin_rebinned = 1; iBin_rebinned <= numBins_rebinned; ++iBin_rebinned ) {
      double binContent_rebinned = 0.;
      double binError2_rebinned = 0.;

      double xMin_rebinnedBin = xAxis->GetBinLowEdge(iBin_rebinned);
      double xMax_rebinnedBin = xAxis->GetBinUpEdge(iBin_rebinned);

      while ( histogram->GetBinCenter(iBin) < xMin_rebinnedBin ) {
	++iBin;
      }

      while ( histogram->GetBinCenter(iBin) >= xMin_rebinnedBin && histogram->GetBinCenter(iBin) < xMax_rebinnedBin ) {
	binContent_rebinned += histogram->GetBinContent(iBin);
	binError2_rebinned += square(histogram->GetBinError(iBin));
	++iBin;
      }

      histogram_rebinned->SetBinContent(iBin_rebinned, binContent_rebinned);
      histogram_rebinned->SetBinError(iBin_rebinned, TMath::Sqrt(binError2_rebinned));
    }

    if ( normalize ) {
      if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();
      histogram_rebinned->Scale(1./histogram_rebinned->Integral());
    }

    std::cout << "histogram(rebinned) = " << histogram_rebinned->GetName() << ":" 
              << " numBins = " << histogram_rebinned->GetNbinsX() << ", integral = " << histogram_rebinned->Integral() << std::endl;
  }

  return histogram_rebinned;
}

TH1* compRatioHistogram(const std::string& ratioHistogramName, const TH1* numerator, const TH1* denominator)
{
  assert(numerator->GetDimension() == denominator->GetDimension());
  assert(numerator->GetNbinsX() == denominator->GetNbinsX());

  TH1* histogramRatio = (TH1*)numerator->Clone(ratioHistogramName.data());
  histogramRatio->Divide(denominator);

  int nBins = histogramRatio->GetNbinsX();
  for ( int iBin = 1; iBin <= nBins; ++iBin ){
    double binContent = histogramRatio->GetBinContent(iBin);
    histogramRatio->SetBinContent(iBin, binContent - 1.);
  }

  histogramRatio->SetLineColor(numerator->GetLineColor());
  histogramRatio->SetLineWidth(numerator->GetLineWidth());
  histogramRatio->SetMarkerColor(numerator->GetMarkerColor());
  histogramRatio->SetMarkerStyle(numerator->GetMarkerStyle());
  histogramRatio->SetMarkerSize(numerator->GetMarkerSize());

  return histogramRatio;
}

void showDistribution(double canvasSizeX, double canvasSizeY,
		      TH1* histogram_ref, const std::string& legendEntry_ref,
		      TH1* histogram2, const std::string& legendEntry2,
		      TH1* histogram3, const std::string& legendEntry3,
		      double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		      double xMin, double xMax, unsigned numBinsMin_rebinned, const std::string& xAxisTitle, double xAxisOffset,
		      bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		      const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.04);
  topPad->SetLeftMargin(0.15);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.05);
  topPad->SetLogy(useLogScale);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  int colors[3] = { 2, 8, 4, };
  int markerStyles[3] = { 34, 20, 21 };
  int markerSizes[3] = { 1, 1, 1 };

  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);

  TH1* histogram_ref_rebinned = rebinHistogram(histogram_ref, numBinsMin_rebinned, xMin, xMax, false);
  histogram_ref_rebinned->SetTitle("");
  histogram_ref_rebinned->SetStats(false);
  histogram_ref_rebinned->SetMinimum(yMin);
  histogram_ref_rebinned->SetMaximum(yMax);
  histogram_ref_rebinned->SetLineColor(colors[0]);
  histogram_ref_rebinned->SetLineWidth(2);
  histogram_ref_rebinned->SetMarkerColor(colors[0]);
  histogram_ref_rebinned->SetMarkerStyle(markerStyles[0]);
  histogram_ref_rebinned->SetMarkerSize(markerSizes[0]);
  histogram_ref_rebinned->Draw("e1p");
  legend->AddEntry(histogram_ref_rebinned, legendEntry_ref.data(), "p");

  TAxis* xAxis_top = histogram_ref_rebinned->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogram_ref_rebinned->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);

  TH1* histogram2_rebinned = 0;
  if ( histogram2 ) {
    histogram2_rebinned = rebinHistogram(histogram2, numBinsMin_rebinned, xMin, xMax, false);
    histogram2_rebinned->SetLineColor(colors[1]);
    histogram2_rebinned->SetLineWidth(2);
    histogram2_rebinned->SetMarkerColor(colors[1]);
    histogram2_rebinned->SetMarkerStyle(markerStyles[1]);
    histogram2_rebinned->SetMarkerSize(markerSizes[1]);
    histogram2_rebinned->Draw("e1psame");
    legend->AddEntry(histogram2_rebinned, legendEntry2.data(), "p");
  }

  TH1* histogram3_rebinned = 0;
  if ( histogram3 ) {
    histogram3_rebinned = rebinHistogram(histogram3, numBinsMin_rebinned, xMin, xMax, false);
    histogram3_rebinned->SetLineColor(colors[2]);
    histogram3_rebinned->SetLineWidth(2);
    histogram3_rebinned->SetMarkerColor(colors[2]);
    histogram3_rebinned->SetMarkerStyle(markerStyles[2]);
    histogram3_rebinned->SetMarkerSize(markerSizes[2]);
    histogram3_rebinned->Draw("e1psame");
    legend->AddEntry(histogram3_rebinned, legendEntry3.data(), "p");
  }

  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogram2_div_ref = 0;
  if ( histogram2 ) {
    std::string histogramName2_div_ref = std::string(histogram2->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram2_div_ref = compRatioHistogram(histogramName2_div_ref, histogram2_rebinned, histogram_ref_rebinned);
    histogram2_div_ref->SetTitle("");
    histogram2_div_ref->SetStats(false);
    histogram2_div_ref->SetMinimum(-0.50);
    histogram2_div_ref->SetMaximum(+0.50);

    TAxis* xAxis_bottom = histogram2_div_ref->GetXaxis();
    xAxis_bottom->SetTitle(xAxis_top->GetTitle());
    xAxis_bottom->SetLabelColor(1);
    xAxis_bottom->SetTitleColor(1);
    xAxis_bottom->SetTitleOffset(1.20);
    xAxis_bottom->SetTitleSize(0.08);
    xAxis_bottom->SetLabelOffset(0.02);
    xAxis_bottom->SetLabelSize(0.08);
    xAxis_bottom->SetTickLength(0.055);
    
    TAxis* yAxis_bottom = histogram2_div_ref->GetYaxis();
    yAxis_bottom->SetTitle("#frac{Embedding - Z/#gamma^{*} #rightarrow #tau #tau}{Z/#gamma^{*} #rightarrow #tau #tau}");
    yAxis_bottom->SetTitleOffset(0.70);
    yAxis_bottom->SetNdivisions(505);
    yAxis_bottom->CenterTitle();
    yAxis_bottom->SetTitleSize(0.08);
    yAxis_bottom->SetLabelSize(0.08);
    yAxis_bottom->SetTickLength(0.04);  
  
    histogram2_div_ref->Draw("e1p");
  }

  TH1* histogram3_div_ref = 0;
  if ( histogram3 ) {
    std::string histogramName3_div_ref = std::string(histogram3->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram3_div_ref = compRatioHistogram(histogramName3_div_ref, histogram3_rebinned, histogram_ref_rebinned);
    histogram3_div_ref->SetTitle("");
    histogram3_div_ref->SetStats(false);
    histogram3_div_ref->SetMinimum(-0.50);
    histogram3_div_ref->SetMaximum(+0.50);

    TAxis* xAxis_bottom = histogram3_div_ref->GetXaxis();
    xAxis_bottom->SetTitle(xAxis_top->GetTitle());
    xAxis_bottom->SetLabelColor(1);
    xAxis_bottom->SetTitleColor(1);
    xAxis_bottom->SetTitleOffset(1.20);
    xAxis_bottom->SetTitleSize(0.08);
    xAxis_bottom->SetLabelOffset(0.02);
    xAxis_bottom->SetLabelSize(0.08);
    xAxis_bottom->SetTickLength(0.055);
    
    TAxis* yAxis_bottom = histogram3_div_ref->GetYaxis();
    yAxis_bottom->SetTitle("#frac{Embedding - Z/#gamma^{*} #rightarrow #tau #tau}{Z/#gamma^{*} #rightarrow #tau #tau}");
    yAxis_bottom->SetTitleOffset(0.70);
    yAxis_bottom->SetNdivisions(505);
    yAxis_bottom->CenterTitle();
    yAxis_bottom->SetTitleSize(0.08);
    yAxis_bottom->SetLabelSize(0.08);
    yAxis_bottom->SetTickLength(0.04);  

    if ( histogram2 ) histogram3_div_ref->Draw("e1psame");
    else histogram3_div_ref->Draw("e1p");
  }

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete legend;
  delete histogram2_div_ref;
  delete histogram3_div_ref;
  delete topPad;
  delete bottomPad;
  delete canvas;  
}

enum { kOldTags, kNewTags };

std::map<std::string, std::map<std::string, plotEntryType*> > makePlots(const TString& inputFileName, const std::string& name, int mode)
{
  TFile* inputFile = TFile::Open(inputFileName.Data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName.Data() << " !!" << std::endl;
    assert(0);
  }

  int maxEvents = -1;
  //int maxEvents = 100000;

  TString treeName = "pfTauIdEffNtupleFromMiniAODProducer/pfTauIdEffNtuple";
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get(treeName.Data()));
  if ( !tree ) {
    std::cerr << "Failed to find tree = " << treeName.Data() << " in input file = " << inputFileName.Data() << " !!" << std::endl;
    assert(0);
  }

  Float_t recTauEn, recTauPx, recTauPy, recTauPz;
  Int_t   recTauDecayMode;
  Float_t leadPFCandEn, leadPFCandPx, leadPFCandPy, leadPFCandPz;
  Float_t leadPFChargedHadrCandEn, leadPFChargedHadrCandPx, leadPFChargedHadrCandPy, leadPFChargedHadrCandPz;
  //Float_t recTauEtaImpactECAL, recTauPhiImpactECAL;
  Float_t genTauEn, genTauPx, genTauPy, genTauPz;
  Int_t   genTauDecayMode;
  Int_t   genTauMatch;
  Float_t genTauDeltaR;
  Float_t genTauLeadChHadEn, genTauLeadChHadPx, genTauLeadChHadPy, genTauLeadChHadPz;
  Int_t   genTauNDaPhotons;
  Float_t genElectronEn, genElectronPx, genElectronPy, genElectronPz;
  Int_t   genElectronMatch;
  Float_t genElectronDeltaR;
  Float_t genMuonEn, genMuonPx, genMuonPy, genMuonPz;
  Int_t   genMuonMatch;
  Float_t genMuonDeltaR;
  Float_t genQuarkOrGluonEn, genQuarkOrGluonPx, genQuarkOrGluonPy, genQuarkOrGluonPz;
  Int_t   genQuarkOrGluonMatch;
  Float_t genQuarkOrGluonDeltaR;
  Float_t recJetEn, recJetPx, recJetPy, recJetPz;
  Int_t   recJetMatch;
  Float_t recJetDeltaR;
  Float_t recJetRawPt;
  Int_t   ElectronVetoMatch;
  Int_t   numVertices;
  Float_t evtWeight, evtWeight_dummy;
  Float_t decayModeFindingNewDMs;
  Float_t decayModeFindingOldDMs;
  Float_t byLooseCombinedIsolationDeltaBetaCorr3Hits;
  Float_t byMediumCombinedIsolationDeltaBetaCorr3Hits;
  Float_t byTightCombinedIsolationDeltaBetaCorr3Hits;
  Float_t byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
  Float_t byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
  Float_t byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
  Float_t byLoosePileupWeightedIsolation3Hits;
  Float_t byMediumPileupWeightedIsolation3Hits;
  Float_t byTightPileupWeightedIsolation3Hits;
  Float_t chargedIsoPtSum;
  Float_t neutralIsoPtSum;
  Float_t puCorrPtSum;
  Float_t againstElectronMVA6raw;
  Float_t againstElectronMVA6category;
  Float_t againstElectronVLooseMVA6;
  Float_t againstElectronLooseMVA6;
  Float_t againstElectronMediumMVA6;
  Float_t againstElectronTightMVA6;
  Float_t againstElectronVTightMVA6;
  Float_t againstMuonLoose3;
  Float_t againstMuonTight3;
  Float_t byVLooseIsolationMVArun2v1DBoldDMwLT;
  Float_t byLooseIsolationMVArun2v1DBoldDMwLT;
  Float_t byMediumIsolationMVArun2v1DBoldDMwLT;
  Float_t byTightIsolationMVArun2v1DBoldDMwLT;
  Float_t byVTightIsolationMVArun2v1DBoldDMwLT;
  Float_t byVVTightIsolationMVArun2v1DBoldDMwLT;
  Float_t byVLooseIsolationMVArun2v1DBnewDMwLT;
  Float_t byLooseIsolationMVArun2v1DBnewDMwLT;
  Float_t byMediumIsolationMVArun2v1DBnewDMwLT;
  Float_t byTightIsolationMVArun2v1DBnewDMwLT;
  Float_t byVTightIsolationMVArun2v1DBnewDMwLT;
  Float_t byVVTightIsolationMVArun2v1DBnewDMwLT;
  Float_t byVLooseIsolationMVArun2v1PWoldDMwLT;
  Float_t byLooseIsolationMVArun2v1PWoldDMwLT;
  Float_t byMediumIsolationMVArun2v1PWoldDMwLT;
  Float_t byTightIsolationMVArun2v1PWoldDMwLT;
  Float_t byVTightIsolationMVArun2v1PWoldDMwLT;
  Float_t byVVTightIsolationMVArun2v1PWoldDMwLT;
  Float_t byVLooseIsolationMVArun2v1PWnewDMwLT;
  Float_t byLooseIsolationMVArun2v1PWnewDMwLT;
  Float_t byMediumIsolationMVArun2v1PWnewDMwLT;
  Float_t byTightIsolationMVArun2v1PWnewDMwLT;
  Float_t byVTightIsolationMVArun2v1PWnewDMwLT;
  Float_t byVVTightIsolationMVArun2v1PWnewDMwLT;
  Float_t byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
  Float_t byLooseIsolationMVArun2v1DBdR03oldDMwLT;
  Float_t byMediumIsolationMVArun2v1DBdR03oldDMwLT;
  Float_t byTightIsolationMVArun2v1DBdR03oldDMwLT;
  Float_t byVTightIsolationMVArun2v1DBdR03oldDMwLT;
  Float_t byVVTightIsolationMVArun2v1DBdR03oldDMwLT;
  Float_t byVLooseIsolationMVArun2v1PWdR03oldDMwLT;
  Float_t byLooseIsolationMVArun2v1PWdR03oldDMwLT;
  Float_t byMediumIsolationMVArun2v1PWdR03oldDMwLT;
  Float_t byTightIsolationMVArun2v1PWdR03oldDMwLT;
  Float_t byVTightIsolationMVArun2v1PWdR03oldDMwLT;
  Float_t byVVTightIsolationMVArun2v1PWdR03oldDMwLT;

  Long64_t run, ls, event;

  tree->SetBranchAddress("recTauEn", &recTauEn);
  tree->SetBranchAddress("recTauPx", &recTauPx);
  tree->SetBranchAddress("recTauPy", &recTauPy);
  tree->SetBranchAddress("recTauPz", &recTauPz);
  tree->SetBranchAddress("recTauDecayMode", &recTauDecayMode);
  tree->SetBranchAddress("leadPFCandEn", &leadPFCandEn);
  tree->SetBranchAddress("leadPFCandPx", &leadPFCandPx);
  tree->SetBranchAddress("leadPFCandPy", &leadPFCandPy);
  tree->SetBranchAddress("leadPFCandPz", &leadPFCandPz);
  tree->SetBranchAddress("leadPFChargedHadrCandEn", &leadPFChargedHadrCandEn);
  tree->SetBranchAddress("leadPFChargedHadrCandPx", &leadPFChargedHadrCandPx); 
  tree->SetBranchAddress("leadPFChargedHadrCandPy", &leadPFChargedHadrCandPy);
  tree->SetBranchAddress("leadPFChargedHadrCandPz", &leadPFChargedHadrCandPz);
  //tree->SetBranchAddress("Tau_etaImpact_ECAL", &recTauEtaImpactECAL);
  //tree->SetBranchAddress("Tau_phiImpact_ECAL", &recTauPhiImpactECAL);
  tree->SetBranchAddress("genTauEn", &genTauEn);
  tree->SetBranchAddress("genTauPx", &genTauPx);
  tree->SetBranchAddress("genTauPy", &genTauPy);
  tree->SetBranchAddress("genTauPz", &genTauPz);
  tree->SetBranchAddress("genTauDecayMode", &genTauDecayMode);
  tree->SetBranchAddress("genTauMatch", &genTauMatch);
  tree->SetBranchAddress("genTauDeltaR", &genTauDeltaR);
  tree->SetBranchAddress("genTauLeadChHadEn", &genTauLeadChHadEn);
  tree->SetBranchAddress("genTauLeadChHadPx", &genTauLeadChHadPx);
  tree->SetBranchAddress("genTauLeadChHadPy", &genTauLeadChHadPy);
  tree->SetBranchAddress("genTauLeadChHadPz", &genTauLeadChHadPz);
  tree->SetBranchAddress("genTauNDaPhotons", &genTauNDaPhotons);
  tree->SetBranchAddress("genElectronEn", &genElectronEn);
  tree->SetBranchAddress("genElectronPx", &genElectronPx);
  tree->SetBranchAddress("genElectronPy", &genElectronPy);
  tree->SetBranchAddress("genElectronPz", &genElectronPz);
  tree->SetBranchAddress("genElectronMatch", &genElectronMatch);
  tree->SetBranchAddress("genElectronDeltaR", &genElectronDeltaR);
  tree->SetBranchAddress("genMuonEn", &genMuonEn);
  tree->SetBranchAddress("genMuonPx", &genMuonPx);
  tree->SetBranchAddress("genMuonPy", &genMuonPy);
  tree->SetBranchAddress("genMuonPz", &genMuonPz);
  tree->SetBranchAddress("genMuonMatch", &genMuonMatch);
  tree->SetBranchAddress("genMuonDeltaR", &genMuonDeltaR);
  tree->SetBranchAddress("genQuarkOrGluonEn", &genQuarkOrGluonEn);
  tree->SetBranchAddress("genQuarkOrGluonPx", &genQuarkOrGluonPx);
  tree->SetBranchAddress("genQuarkOrGluonPy", &genQuarkOrGluonPy);
  tree->SetBranchAddress("genQuarkOrGluonPz", &genQuarkOrGluonPz);
  tree->SetBranchAddress("genQuarkOrGluonMatch", &genQuarkOrGluonMatch);
  tree->SetBranchAddress("genQuarkOrGluonDeltaR", &genQuarkOrGluonDeltaR);
  tree->SetBranchAddress("recJetEn", &recJetEn);
  tree->SetBranchAddress("recJetPx", &recJetPx);
  tree->SetBranchAddress("recJetPy", &recJetPy);
  tree->SetBranchAddress("recJetPz", &recJetPz);
  tree->SetBranchAddress("recJetMatch", &recJetMatch);
  tree->SetBranchAddress("recJetDeltaR", &recJetDeltaR);
  tree->SetBranchAddress("recJetRawPt", &recJetRawPt);
  tree->SetBranchAddress("ElectronVetoMatch", &ElectronVetoMatch);
  tree->SetBranchAddress("numVertices", &numVertices);
  tree->SetBranchAddress("evtWeight", &evtWeight);
  if ( mode == kOldTags ) {
    tree->SetBranchAddress("decayModeFinding", &decayModeFindingOldDMs);
  } else if ( mode == kNewTags ) {
    tree->SetBranchAddress("decayModeFindingNewDMs", &decayModeFindingNewDMs);
    tree->SetBranchAddress("decayModeFinding", &decayModeFindingOldDMs);
  } else assert(0);
  tree->SetBranchAddress("byLooseCombinedIsolationDeltaBetaCorr3Hits", &byLooseCombinedIsolationDeltaBetaCorr3Hits);
  tree->SetBranchAddress("byMediumCombinedIsolationDeltaBetaCorr3Hits", &byMediumCombinedIsolationDeltaBetaCorr3Hits);
  tree->SetBranchAddress("byTightCombinedIsolationDeltaBetaCorr3Hits", &byTightCombinedIsolationDeltaBetaCorr3Hits);
  tree->SetBranchAddress("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03", &byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
  tree->SetBranchAddress("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03", &byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
  tree->SetBranchAddress("byTightCombinedIsolationDeltaBetaCorr3HitsdR03", &byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
  tree->SetBranchAddress("byLoosePileupWeightedIsolation3Hits", &byLoosePileupWeightedIsolation3Hits);
  tree->SetBranchAddress("byMediumPileupWeightedIsolation3Hits", &byMediumPileupWeightedIsolation3Hits);
  tree->SetBranchAddress("byTightPileupWeightedIsolation3Hits", &byTightPileupWeightedIsolation3Hits);
  if ( mode == kNewTags ) {
    tree->SetBranchAddress("chargedIsoPtSum", &chargedIsoPtSum);
    tree->SetBranchAddress("neutralIsoPtSum", &neutralIsoPtSum);
    tree->SetBranchAddress("puCorrPtSum", &puCorrPtSum);
    tree->SetBranchAddress("byVLooseIsolationMVArun2v1DBoldDMwLT", &byVLooseIsolationMVArun2v1DBoldDMwLT);
    tree->SetBranchAddress("byLooseIsolationMVArun2v1DBoldDMwLT", &byLooseIsolationMVArun2v1DBoldDMwLT);
    tree->SetBranchAddress("byMediumIsolationMVArun2v1DBoldDMwLT", &byMediumIsolationMVArun2v1DBoldDMwLT);
    tree->SetBranchAddress("byTightIsolationMVArun2v1DBoldDMwLT", &byTightIsolationMVArun2v1DBoldDMwLT);
    tree->SetBranchAddress("byVTightIsolationMVArun2v1DBoldDMwLT", &byVTightIsolationMVArun2v1DBoldDMwLT);
    tree->SetBranchAddress("byVVTightIsolationMVArun2v1DBoldDMwLT", &byVVTightIsolationMVArun2v1DBoldDMwLT);
    tree->SetBranchAddress("byVLooseIsolationMVArun2v1DBnewDMwLT", &byVLooseIsolationMVArun2v1DBnewDMwLT);
    tree->SetBranchAddress("byLooseIsolationMVArun2v1DBnewDMwLT", &byLooseIsolationMVArun2v1DBnewDMwLT);
    tree->SetBranchAddress("byMediumIsolationMVArun2v1DBnewDMwLT", &byMediumIsolationMVArun2v1DBnewDMwLT);
    tree->SetBranchAddress("byTightIsolationMVArun2v1DBnewDMwLT", &byTightIsolationMVArun2v1DBnewDMwLT);
    tree->SetBranchAddress("byVTightIsolationMVArun2v1DBnewDMwLT", &byVTightIsolationMVArun2v1DBnewDMwLT);
    tree->SetBranchAddress("byVVTightIsolationMVArun2v1DBnewDMwLT", &byVVTightIsolationMVArun2v1DBnewDMwLT);
    tree->SetBranchAddress("byVLooseIsolationMVArun2v1PWoldDMwLT", &byVLooseIsolationMVArun2v1PWoldDMwLT);
    tree->SetBranchAddress("byLooseIsolationMVArun2v1PWoldDMwLT", &byLooseIsolationMVArun2v1PWoldDMwLT);
    tree->SetBranchAddress("byMediumIsolationMVArun2v1PWoldDMwLT", &byMediumIsolationMVArun2v1PWoldDMwLT);
    tree->SetBranchAddress("byTightIsolationMVArun2v1PWoldDMwLT", &byTightIsolationMVArun2v1PWoldDMwLT);
    tree->SetBranchAddress("byVTightIsolationMVArun2v1PWoldDMwLT", &byVTightIsolationMVArun2v1PWoldDMwLT);
    tree->SetBranchAddress("byVVTightIsolationMVArun2v1PWoldDMwLT", &byVVTightIsolationMVArun2v1PWoldDMwLT);
    tree->SetBranchAddress("byVLooseIsolationMVArun2v1PWnewDMwLT", &byVLooseIsolationMVArun2v1PWnewDMwLT);
    tree->SetBranchAddress("byLooseIsolationMVArun2v1PWnewDMwLT", &byLooseIsolationMVArun2v1PWnewDMwLT);
    tree->SetBranchAddress("byMediumIsolationMVArun2v1PWnewDMwLT", &byMediumIsolationMVArun2v1PWnewDMwLT);
    tree->SetBranchAddress("byTightIsolationMVArun2v1PWnewDMwLT", &byTightIsolationMVArun2v1PWnewDMwLT);
    tree->SetBranchAddress("byVTightIsolationMVArun2v1PWnewDMwLT", &byVTightIsolationMVArun2v1PWnewDMwLT);
    tree->SetBranchAddress("byVVTightIsolationMVArun2v1PWnewDMwLT", &byVVTightIsolationMVArun2v1PWnewDMwLT);
    tree->SetBranchAddress("byVLooseIsolationMVArun2v1DBdR03oldDMwLT", &byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
    tree->SetBranchAddress("byLooseIsolationMVArun2v1DBdR03oldDMwLT", &byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    tree->SetBranchAddress("byMediumIsolationMVArun2v1DBdR03oldDMwLT", &byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    tree->SetBranchAddress("byTightIsolationMVArun2v1DBdR03oldDMwLT", &byTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->SetBranchAddress("byVTightIsolationMVArun2v1DBdR03oldDMwLT", &byVTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->SetBranchAddress("byVVTightIsolationMVArun2v1DBdR03oldDMwLT", &byVVTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->SetBranchAddress("byVLooseIsolationMVArun2v1PWdR03oldDMwLT", &byVLooseIsolationMVArun2v1PWdR03oldDMwLT);
    tree->SetBranchAddress("byLooseIsolationMVArun2v1PWdR03oldDMwLT", &byLooseIsolationMVArun2v1PWdR03oldDMwLT);
    tree->SetBranchAddress("byMediumIsolationMVArun2v1PWdR03oldDMwLT", &byMediumIsolationMVArun2v1PWdR03oldDMwLT);
    tree->SetBranchAddress("byTightIsolationMVArun2v1PWdR03oldDMwLT", &byTightIsolationMVArun2v1PWdR03oldDMwLT);
    tree->SetBranchAddress("byVTightIsolationMVArun2v1PWdR03oldDMwLT", &byVTightIsolationMVArun2v1PWdR03oldDMwLT);
    tree->SetBranchAddress("byVVTightIsolationMVArun2v1PWdR03oldDMwLT", &byVVTightIsolationMVArun2v1PWdR03oldDMwLT);
  }
  if ( mode == kNewTags ) {
    tree->SetBranchAddress("againstElectronMVA6raw", &againstElectronMVA6raw);
    tree->SetBranchAddress("againstElectronMVA6category", &againstElectronMVA6category);
    tree->SetBranchAddress("againstElectronVLooseMVA6", &againstElectronVLooseMVA6);
    tree->SetBranchAddress("againstElectronLooseMVA6", &againstElectronLooseMVA6);
    tree->SetBranchAddress("againstElectronMediumMVA6", &againstElectronMediumMVA6);
    tree->SetBranchAddress("againstElectronTightMVA6", &againstElectronTightMVA6);
    tree->SetBranchAddress("againstElectronVTightMVA6", &againstElectronVTightMVA6);
  }
  if ( mode == kNewTags ) {
    tree->SetBranchAddress("againstMuonLoose3", &againstMuonLoose3);
    tree->SetBranchAddress("againstMuonTight3", &againstMuonTight3);
  }
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("ls", &ls);
  tree->SetBranchAddress("event", &event);

  std::map<std::string, std::map<std::string, plotEntryType*> > plotEntryMap; // key = particle type, discriminator

  plotEntryType* plotsTau_recoTauJetmatching                    = new plotEntryType(Form("%s_tau_recoTauJetmatching", name.data()));
  plotEntryType* plotsTau_byDecayModeFindingNewDMs              = new plotEntryType(Form("%s_tau_byDecayModeFindingNewDMs", name.data()));
  plotEntryType* plotsTau_byDecayModeFindingOldDMs              = new plotEntryType(Form("%s_tau_byDecayModeFindingOldDMs", name.data()));
  plotEntryType* plotsTau_byCombinedIsolation3HitsLoose         = new plotEntryType(Form("%s_tau_byCombinedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsTau_byCombinedIsolation3HitsMedium        = new plotEntryType(Form("%s_tau_byCombinedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsTau_byCombinedIsolation3HitsTight         = new plotEntryType(Form("%s_tau_byCombinedIsolation3HitsTight", name.data()));
  plotEntryType* plotsTau_byCombinedIsolationNewDM3HitsLoose         = new plotEntryType(Form("%s_tau_byCombinedIsolationNewDM3HitsLoose", name.data()));
  plotEntryType* plotsTau_byCombinedIsolationNewDM3HitsMedium        = new plotEntryType(Form("%s_tau_byCombinedIsolationNewDM3HitsMedium", name.data()));
  plotEntryType* plotsTau_byCombinedIsolationNewDM3HitsTight         = new plotEntryType(Form("%s_tau_byCombinedIsolationNewDM3HitsTight", name.data()));
  plotEntryType* plotsTau_byCombinedIsolation3HitsdR03Loose         = new plotEntryType(Form("%s_tau_byCombinedIsolation3HitsdR03Loose", name.data()));
  plotEntryType* plotsTau_byCombinedIsolation3HitsdR03Medium        = new plotEntryType(Form("%s_tau_byCombinedIsolation3HitsdR03Medium", name.data()));
  plotEntryType* plotsTau_byCombinedIsolation3HitsdR03Tight         = new plotEntryType(Form("%s_tau_byCombinedIsolation3HitsdR03Tight", name.data()));
  plotEntryType* plotsTau_byChargedIsolationLoose               = new plotEntryType(Form("%s_tau_byChargedIsolationLoose", name.data()));
  plotEntryType* plotsTau_byChargedIsolationMedium              = new plotEntryType(Form("%s_tau_byChargedIsolationMedium", name.data()));
  plotEntryType* plotsTau_byChargedIsolationTight               = new plotEntryType(Form("%s_tau_byChargedIsolationTight", name.data()));
  
  plotEntryType* plotsTau_byPileupWeightedIsolation3HitsLoose         = new plotEntryType(Form("%s_tau_byPileupWeightedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsTau_byPileupWeightedIsolation3HitsMedium        = new plotEntryType(Form("%s_tau_byPileupWeightedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsTau_byPileupWeightedIsolation3HitsTight         = new plotEntryType(Form("%s_tau_byPileupWeightedIsolation3HitsTight", name.data()));
  plotEntryType* plotsTau_byPileupWeightedIsolationNewDM3HitsLoose         = new plotEntryType(Form("%s_tau_byPileupWeightedIsolationNewDM3HitsLoose", name.data()));
  plotEntryType* plotsTau_byPileupWeightedIsolationNewDM3HitsMedium        = new plotEntryType(Form("%s_tau_byPileupWeightedIsolationNewDM3HitsMedium", name.data()));
  plotEntryType* plotsTau_byPileupWeightedIsolationNewDM3HitsTight         = new plotEntryType(Form("%s_tau_byPileupWeightedIsolationNewDM3HitsTight", name.data()));
  plotEntryType* plotsTau_byVLooseIsolationMVArun2v1DBoldDMwLT         = new plotEntryType(Form("%s_tau_byVLooseIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsTau_byLooseIsolationMVArun2v1DBoldDMwLT          = new plotEntryType(Form("%s_tau_byLooseIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsTau_byMediumIsolationMVArun2v1DBoldDMwLT         = new plotEntryType(Form("%s_tau_byMediumIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsTau_byTightIsolationMVArun2v1DBoldDMwLT          = new plotEntryType(Form("%s_tau_byTightIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsTau_byVTightIsolationMVArun2v1DBoldDMwLT         = new plotEntryType(Form("%s_tau_byVTightIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsTau_byVVTightIsolationMVArun2v1DBoldDMwLT        = new plotEntryType(Form("%s_tau_byVVTightIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsTau_byVLooseIsolationMVArun2v1DBnewDMwLT         = new plotEntryType(Form("%s_tau_byVLooseIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsTau_byLooseIsolationMVArun2v1DBnewDMwLT          = new plotEntryType(Form("%s_tau_byLooseIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsTau_byMediumIsolationMVArun2v1DBnewDMwLT         = new plotEntryType(Form("%s_tau_byMediumIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsTau_byTightIsolationMVArun2v1DBnewDMwLT          = new plotEntryType(Form("%s_tau_byTightIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsTau_byVTightIsolationMVArun2v1DBnewDMwLT         = new plotEntryType(Form("%s_tau_byVTightIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsTau_byVVTightIsolationMVArun2v1DBnewDMwLT        = new plotEntryType(Form("%s_tau_byVVTightIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsTau_byVLooseIsolationMVArun2v1PWoldDMwLT         = new plotEntryType(Form("%s_tau_byVLooseIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsTau_byLooseIsolationMVArun2v1PWoldDMwLT          = new plotEntryType(Form("%s_tau_byLooseIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsTau_byMediumIsolationMVArun2v1PWoldDMwLT         = new plotEntryType(Form("%s_tau_byMediumIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsTau_byTightIsolationMVArun2v1PWoldDMwLT          = new plotEntryType(Form("%s_tau_byTightIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsTau_byVTightIsolationMVArun2v1PWoldDMwLT         = new plotEntryType(Form("%s_tau_byVTightIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsTau_byVVTightIsolationMVArun2v1PWoldDMwLT        = new plotEntryType(Form("%s_tau_byVVTightIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsTau_byVLooseIsolationMVArun2v1PWnewDMwLT         = new plotEntryType(Form("%s_tau_byVLooseIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsTau_byLooseIsolationMVArun2v1PWnewDMwLT          = new plotEntryType(Form("%s_tau_byLooseIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsTau_byMediumIsolationMVArun2v1PWnewDMwLT         = new plotEntryType(Form("%s_tau_byMediumIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsTau_byTightIsolationMVArun2v1PWnewDMwLT          = new plotEntryType(Form("%s_tau_byTightIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsTau_byVTightIsolationMVArun2v1PWnewDMwLT         = new plotEntryType(Form("%s_tau_byVTightIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsTau_byVVTightIsolationMVArun2v1PWnewDMwLT        = new plotEntryType(Form("%s_tau_byVVTightIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsTau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT         = new plotEntryType(Form("%s_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT          = new plotEntryType(Form("%s_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT         = new plotEntryType(Form("%s_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byTightIsolationMVArun2v1DBdR03oldDMwLT          = new plotEntryType(Form("%s_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byVTightIsolationMVArun2v1DBdR03oldDMwLT         = new plotEntryType(Form("%s_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT        = new plotEntryType(Form("%s_tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT         = new plotEntryType(Form("%s_tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byLooseIsolationMVArun2v1PWdR03oldDMwLT          = new plotEntryType(Form("%s_tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byMediumIsolationMVArun2v1PWdR03oldDMwLT         = new plotEntryType(Form("%s_tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byTightIsolationMVArun2v1PWdR03oldDMwLT          = new plotEntryType(Form("%s_tau_byTightIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byVTightIsolationMVArun2v1PWdR03oldDMwLT         = new plotEntryType(Form("%s_tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsTau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT        = new plotEntryType(Form("%s_tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryMap["Tau"]["recoTauJetmatching"]                     = plotsTau_recoTauJetmatching;
  plotEntryMap["Tau"]["byDecayModeFindingNewDMs"]               = plotsTau_byDecayModeFindingNewDMs;
  plotEntryMap["Tau"]["byDecayModeFindingOldDMs"]               = plotsTau_byDecayModeFindingOldDMs;
  plotEntryMap["Tau"]["byCombinedIsolation3HitsLoose"]          = plotsTau_byCombinedIsolation3HitsLoose;
  plotEntryMap["Tau"]["byCombinedIsolation3HitsMedium"]         = plotsTau_byCombinedIsolation3HitsMedium;
  plotEntryMap["Tau"]["byCombinedIsolation3HitsTight"]          = plotsTau_byCombinedIsolation3HitsTight;
  plotEntryMap["Tau"]["byCombinedIsolationNewDM3HitsLoose"]          = plotsTau_byCombinedIsolationNewDM3HitsLoose;
  plotEntryMap["Tau"]["byCombinedIsolationNewDM3HitsMedium"]         = plotsTau_byCombinedIsolationNewDM3HitsMedium;
  plotEntryMap["Tau"]["byCombinedIsolationNewDM3HitsTight"]          = plotsTau_byCombinedIsolationNewDM3HitsTight;
  plotEntryMap["Tau"]["byCombinedIsolation3HitsdR03Loose"]          = plotsTau_byCombinedIsolation3HitsdR03Loose;
  plotEntryMap["Tau"]["byCombinedIsolation3HitsdR03Medium"]         = plotsTau_byCombinedIsolation3HitsdR03Medium;
  plotEntryMap["Tau"]["byCombinedIsolation3HitsdR03Tight"]          = plotsTau_byCombinedIsolation3HitsdR03Tight;
  plotEntryMap["Tau"]["byChargedIsolationLoose"]                = plotsTau_byChargedIsolationLoose;
  plotEntryMap["Tau"]["byChargedIsolationMedium"]               = plotsTau_byChargedIsolationMedium;
  plotEntryMap["Tau"]["byChargedIsolationTight"]                = plotsTau_byChargedIsolationTight;
  plotEntryMap["Tau"]["byPileupWeightedIsolation3HitsLoose"]          = plotsTau_byPileupWeightedIsolation3HitsLoose;
  plotEntryMap["Tau"]["byPileupWeightedIsolation3HitsMedium"]         = plotsTau_byPileupWeightedIsolation3HitsMedium;
  plotEntryMap["Tau"]["byPileupWeightedIsolation3HitsTight"]          = plotsTau_byPileupWeightedIsolation3HitsTight;
  plotEntryMap["Tau"]["byPileupWeightedIsolationNewDM3HitsLoose"]          = plotsTau_byPileupWeightedIsolationNewDM3HitsLoose;
  plotEntryMap["Tau"]["byPileupWeightedIsolationNewDM3HitsMedium"]         = plotsTau_byPileupWeightedIsolationNewDM3HitsMedium;
  plotEntryMap["Tau"]["byPileupWeightedIsolationNewDM3HitsTight"]          = plotsTau_byPileupWeightedIsolationNewDM3HitsTight;
  plotEntryMap["Tau"]["byVLooseIsolationMVArun2v1DBoldDMwLT"]          = plotsTau_byVLooseIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Tau"]["byLooseIsolationMVArun2v1DBoldDMwLT"]           = plotsTau_byLooseIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Tau"]["byMediumIsolationMVArun2v1DBoldDMwLT"]          = plotsTau_byMediumIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Tau"]["byTightIsolationMVArun2v1DBoldDMwLT"]           = plotsTau_byTightIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Tau"]["byVTightIsolationMVArun2v1DBoldDMwLT"]          = plotsTau_byVTightIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Tau"]["byVVTightIsolationMVArun2v1DBoldDMwLT"]          = plotsTau_byVVTightIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Tau"]["byLooseIsolationMVArun2v1DBnewDMwLT"]           = plotsTau_byLooseIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Tau"]["byMediumIsolationMVArun2v1DBnewDMwLT"]          = plotsTau_byMediumIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Tau"]["byTightIsolationMVArun2v1DBnewDMwLT"]           = plotsTau_byTightIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Tau"]["byVTightIsolationMVArun2v1DBnewDMwLT"]          = plotsTau_byVTightIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Tau"]["byLooseIsolationMVArun2v1PWoldDMwLT"]           = plotsTau_byLooseIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Tau"]["byMediumIsolationMVArun2v1PWoldDMwLT"]          = plotsTau_byMediumIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Tau"]["byTightIsolationMVArun2v1PWoldDMwLT"]           = plotsTau_byTightIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Tau"]["byVTightIsolationMVArun2v1PWoldDMwLT"]          = plotsTau_byVTightIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Tau"]["byLooseIsolationMVArun2v1PWnewDMwLT"]           = plotsTau_byLooseIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Tau"]["byMediumIsolationMVArun2v1PWnewDMwLT"]          = plotsTau_byMediumIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Tau"]["byTightIsolationMVArun2v1PWnewDMwLT"]           = plotsTau_byTightIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Tau"]["byVTightIsolationMVArun2v1PWnewDMwLT"]          = plotsTau_byVTightIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Tau"]["byLooseIsolationMVArun2v1DBdR03oldDMwLT"]           = plotsTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Tau"]["byMediumIsolationMVArun2v1DBdR03oldDMwLT"]          = plotsTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Tau"]["byTightIsolationMVArun2v1DBdR03oldDMwLT"]           = plotsTau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Tau"]["byVTightIsolationMVArun2v1DBdR03oldDMwLT"]          = plotsTau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Tau"]["byLooseIsolationMVArun2v1PWdR03oldDMwLT"]           = plotsTau_byLooseIsolationMVArun2v1PWdR03oldDMwLT;
  plotEntryMap["Tau"]["byMediumIsolationMVArun2v1PWdR03oldDMwLT"]          = plotsTau_byMediumIsolationMVArun2v1PWdR03oldDMwLT;
  plotEntryMap["Tau"]["byTightIsolationMVArun2v1PWdR03oldDMwLT"]           = plotsTau_byTightIsolationMVArun2v1PWdR03oldDMwLT;
  plotEntryMap["Tau"]["byVTightIsolationMVArun2v1PWdR03oldDMwLT"]          = plotsTau_byVTightIsolationMVArun2v1PWdR03oldDMwLT;


  /*
  plotEntryType* plotsTau_genDM0_recoTauJetmatching                    = new plotEntryType(Form("%s_tau_gendm0_recoTauJetmatching", name.data()));
  plotEntryType* plotsTau_genDM0_byDecayModeFindingNewDMs              = new plotEntryType(Form("%s_tau_gendm0_byDecayModeFindingNewDMs", name.data()));
  plotEntryType* plotsTau_genDM0_byDecayModeFindingOldDMs              = new plotEntryType(Form("%s_tau_gendm0_byDecayModeFindingOldDMs", name.data()));
  plotEntryType* plotsTau_genDM0_byCombinedIsolation8HitsLoose         = new plotEntryType(Form("%s_tau_gendm0_byCombinedIsolation8HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM0_byCombinedIsolation8HitsMedium        = new plotEntryType(Form("%s_tau_gendm0_byCombinedIsolation8HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM0_byCombinedIsolation8HitsTight         = new plotEntryType(Form("%s_tau_gendm0_byCombinedIsolation8HitsTight", name.data()));
  plotEntryType* plotsTau_genDM0_byCombinedIsolation3HitsLoose         = new plotEntryType(Form("%s_tau_gendm0_byCombinedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM0_byCombinedIsolation3HitsMedium        = new plotEntryType(Form("%s_tau_gendm0_byCombinedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM0_byCombinedIsolation3HitsTight         = new plotEntryType(Form("%s_tau_gendm0_byCombinedIsolation3HitsTight", name.data()));
  plotEntryType* plotsTau_genDM0_byChargedIsolationLoose               = new plotEntryType(Form("%s_tau_gendm0_byChargedIsolationLoose", name.data()));
  plotEntryType* plotsTau_genDM0_byChargedIsolationMedium              = new plotEntryType(Form("%s_tau_gendm0_byChargedIsolationMedium", name.data()));
  plotEntryType* plotsTau_genDM0_byChargedIsolationTight               = new plotEntryType(Form("%s_tau_gendm0_byChargedIsolationTight", name.data()));
  plotEntryMap["Tau_genDM0"]["recoTauJetmatching"]                     = plotsTau_genDM0_recoTauJetmatching;
  plotEntryMap["Tau_genDM0"]["byDecayModeFindingNewDMs"]               = plotsTau_genDM0_byDecayModeFindingNewDMs;
  plotEntryMap["Tau_genDM0"]["byDecayModeFindingOldDMs"]               = plotsTau_genDM0_byDecayModeFindingOldDMs;
  plotEntryMap["Tau_genDM0"]["byCombinedIsolation8HitsLoose"]          = plotsTau_genDM0_byCombinedIsolation8HitsLoose;
  plotEntryMap["Tau_genDM0"]["byCombinedIsolation8HitsMedium"]         = plotsTau_genDM0_byCombinedIsolation8HitsMedium;
  plotEntryMap["Tau_genDM0"]["byCombinedIsolation8HitsTight"]          = plotsTau_genDM0_byCombinedIsolation8HitsTight;
  plotEntryMap["Tau_genDM0"]["byCombinedIsolation3HitsLoose"]          = plotsTau_genDM0_byCombinedIsolation3HitsLoose;
  plotEntryMap["Tau_genDM0"]["byCombinedIsolation3HitsMedium"]         = plotsTau_genDM0_byCombinedIsolation3HitsMedium;
  plotEntryMap["Tau_genDM0"]["byCombinedIsolation3HitsTight"]          = plotsTau_genDM0_byCombinedIsolation3HitsTight;
  plotEntryMap["Tau_genDM0"]["byChargedIsolationLoose"]                = plotsTau_genDM0_byChargedIsolationLoose;
  plotEntryMap["Tau_genDM0"]["byChargedIsolationMedium"]               = plotsTau_genDM0_byChargedIsolationMedium;
  plotEntryMap["Tau_genDM0"]["byChargedIsolationTight"]                = plotsTau_genDM0_byChargedIsolationTight;

  plotEntryType* plotsTau_genDM1_recoTauJetmatching                    = new plotEntryType(Form("%s_tau_gendm1_recoTauJetmatching", name.data()));
  plotEntryType* plotsTau_genDM1_byDecayModeFindingNewDMs              = new plotEntryType(Form("%s_tau_gendm1_byDecayModeFindingNewDMs", name.data()));
  plotEntryType* plotsTau_genDM1_byDecayModeFindingOldDMs              = new plotEntryType(Form("%s_tau_gendm1_byDecayModeFindingOldDMs", name.data()));
  plotEntryType* plotsTau_genDM1_byCombinedIsolation8HitsLoose         = new plotEntryType(Form("%s_tau_gendm1_byCombinedIsolation8HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM1_byCombinedIsolation8HitsMedium        = new plotEntryType(Form("%s_tau_gendm1_byCombinedIsolation8HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM1_byCombinedIsolation8HitsTight         = new plotEntryType(Form("%s_tau_gendm1_byCombinedIsolation8HitsTight", name.data()));
  plotEntryType* plotsTau_genDM1_byCombinedIsolation3HitsLoose         = new plotEntryType(Form("%s_tau_gendm1_byCombinedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM1_byCombinedIsolation3HitsMedium        = new plotEntryType(Form("%s_tau_gendm1_byCombinedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM1_byCombinedIsolation3HitsTight         = new plotEntryType(Form("%s_tau_gendm1_byCombinedIsolation3HitsTight", name.data()));
  plotEntryType* plotsTau_genDM1_byChargedIsolationLoose               = new plotEntryType(Form("%s_tau_gendm1_byChargedIsolationLoose", name.data()));
  plotEntryType* plotsTau_genDM1_byChargedIsolationMedium              = new plotEntryType(Form("%s_tau_gendm1_byChargedIsolationMedium", name.data()));
  plotEntryType* plotsTau_genDM1_byChargedIsolationTight               = new plotEntryType(Form("%s_tau_gendm1_byChargedIsolationTight", name.data()));
  plotEntryMap["Tau_genDM1"]["recoTauJetmatching"]                     = plotsTau_genDM1_recoTauJetmatching;
  plotEntryMap["Tau_genDM1"]["byDecayModeFindingNewDMs"]               = plotsTau_genDM1_byDecayModeFindingNewDMs;
  plotEntryMap["Tau_genDM1"]["byDecayModeFindingOldDMs"]               = plotsTau_genDM1_byDecayModeFindingOldDMs;
  plotEntryMap["Tau_genDM1"]["byCombinedIsolation8HitsLoose"]          = plotsTau_genDM1_byCombinedIsolation8HitsLoose;
  plotEntryMap["Tau_genDM1"]["byCombinedIsolation8HitsMedium"]         = plotsTau_genDM1_byCombinedIsolation8HitsMedium;
  plotEntryMap["Tau_genDM1"]["byCombinedIsolation8HitsTight"]          = plotsTau_genDM1_byCombinedIsolation8HitsTight;
  plotEntryMap["Tau_genDM1"]["byCombinedIsolation3HitsLoose"]          = plotsTau_genDM1_byCombinedIsolation3HitsLoose;
  plotEntryMap["Tau_genDM1"]["byCombinedIsolation3HitsMedium"]         = plotsTau_genDM1_byCombinedIsolation3HitsMedium;
  plotEntryMap["Tau_genDM1"]["byCombinedIsolation3HitsTight"]          = plotsTau_genDM1_byCombinedIsolation3HitsTight;
  plotEntryMap["Tau_genDM1"]["byChargedIsolationLoose"]                = plotsTau_genDM1_byChargedIsolationLoose;
  plotEntryMap["Tau_genDM1"]["byChargedIsolationMedium"]               = plotsTau_genDM1_byChargedIsolationMedium;
  plotEntryMap["Tau_genDM1"]["byChargedIsolationTight"]                = plotsTau_genDM1_byChargedIsolationTight;

  plotEntryType* plotsTau_genDM2_recoTauJetmatching                    = new plotEntryType(Form("%s_tau_gendm2_recoTauJetmatching", name.data()));
  plotEntryType* plotsTau_genDM2_byDecayModeFindingNewDMs              = new plotEntryType(Form("%s_tau_gendm2_byDecayModeFindingNewDMs", name.data()));
  plotEntryType* plotsTau_genDM2_byDecayModeFindingOldDMs              = new plotEntryType(Form("%s_tau_gendm2_byDecayModeFindingOldDMs", name.data()));
  plotEntryType* plotsTau_genDM2_byCombinedIsolation8HitsLoose         = new plotEntryType(Form("%s_tau_gendm2_byCombinedIsolation8HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM2_byCombinedIsolation8HitsMedium        = new plotEntryType(Form("%s_tau_gendm2_byCombinedIsolation8HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM2_byCombinedIsolation8HitsTight         = new plotEntryType(Form("%s_tau_gendm2_byCombinedIsolation8HitsTight", name.data()));
  plotEntryType* plotsTau_genDM2_byCombinedIsolation3HitsLoose         = new plotEntryType(Form("%s_tau_gendm2_byCombinedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM2_byCombinedIsolation3HitsMedium        = new plotEntryType(Form("%s_tau_gendm2_byCombinedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM2_byCombinedIsolation3HitsTight         = new plotEntryType(Form("%s_tau_gendm2_byCombinedIsolation3HitsTight", name.data()));
  plotEntryType* plotsTau_genDM2_byChargedIsolationLoose               = new plotEntryType(Form("%s_tau_gendm2_byChargedIsolationLoose", name.data()));
  plotEntryType* plotsTau_genDM2_byChargedIsolationMedium              = new plotEntryType(Form("%s_tau_gendm2_byChargedIsolationMedium", name.data()));
  plotEntryType* plotsTau_genDM2_byChargedIsolationTight               = new plotEntryType(Form("%s_tau_gendm2_byChargedIsolationTight", name.data()));
  plotEntryMap["Tau_genDM2"]["recoTauJetmatching"]                     = plotsTau_genDM2_recoTauJetmatching;
  plotEntryMap["Tau_genDM2"]["byDecayModeFindingNewDMs"]               = plotsTau_genDM2_byDecayModeFindingNewDMs;
  plotEntryMap["Tau_genDM2"]["byDecayModeFindingOldDMs"]               = plotsTau_genDM2_byDecayModeFindingOldDMs;
  plotEntryMap["Tau_genDM2"]["byCombinedIsolation8HitsLoose"]          = plotsTau_genDM2_byCombinedIsolation8HitsLoose;
  plotEntryMap["Tau_genDM2"]["byCombinedIsolation8HitsMedium"]         = plotsTau_genDM2_byCombinedIsolation8HitsMedium;
  plotEntryMap["Tau_genDM2"]["byCombinedIsolation8HitsTight"]          = plotsTau_genDM2_byCombinedIsolation8HitsTight;
  plotEntryMap["Tau_genDM2"]["byCombinedIsolation3HitsLoose"]          = plotsTau_genDM2_byCombinedIsolation3HitsLoose;
  plotEntryMap["Tau_genDM2"]["byCombinedIsolation3HitsMedium"]         = plotsTau_genDM2_byCombinedIsolation3HitsMedium;
  plotEntryMap["Tau_genDM2"]["byCombinedIsolation3HitsTight"]          = plotsTau_genDM2_byCombinedIsolation3HitsTight;
  plotEntryMap["Tau_genDM2"]["byChargedIsolationLoose"]                = plotsTau_genDM2_byChargedIsolationLoose;
  plotEntryMap["Tau_genDM2"]["byChargedIsolationMedium"]               = plotsTau_genDM2_byChargedIsolationMedium;
  plotEntryMap["Tau_genDM2"]["byChargedIsolationTight"]                = plotsTau_genDM2_byChargedIsolationTight;

  plotEntryType* plotsTau_genDM10_recoTauJetmatching                    = new plotEntryType(Form("%s_tau_gendm10_recoTauJetmatching", name.data()));
  plotEntryType* plotsTau_genDM10_byDecayModeFindingNewDMs              = new plotEntryType(Form("%s_tau_gendm10_byDecayModeFindingNewDMs", name.data()));
  plotEntryType* plotsTau_genDM10_byDecayModeFindingOldDMs              = new plotEntryType(Form("%s_tau_gendm10_byDecayModeFindingOldDMs", name.data()));
  plotEntryType* plotsTau_genDM10_byCombinedIsolation8HitsLoose         = new plotEntryType(Form("%s_tau_gendm10_byCombinedIsolation8HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM10_byCombinedIsolation8HitsMedium        = new plotEntryType(Form("%s_tau_gendm10_byCombinedIsolation8HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM10_byCombinedIsolation8HitsTight         = new plotEntryType(Form("%s_tau_gendm10_byCombinedIsolation8HitsTight", name.data()));
  plotEntryType* plotsTau_genDM10_byCombinedIsolation3HitsLoose         = new plotEntryType(Form("%s_tau_gendm10_byCombinedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsTau_genDM10_byCombinedIsolation3HitsMedium        = new plotEntryType(Form("%s_tau_gendm10_byCombinedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsTau_genDM10_byCombinedIsolation3HitsTight         = new plotEntryType(Form("%s_tau_gendm10_byCombinedIsolation3HitsTight", name.data()));
  plotEntryType* plotsTau_genDM10_byChargedIsolationLoose               = new plotEntryType(Form("%s_tau_gendm10_byChargedIsolationLoose", name.data()));
  plotEntryType* plotsTau_genDM10_byChargedIsolationMedium              = new plotEntryType(Form("%s_tau_gendm10_byChargedIsolationMedium", name.data()));
  plotEntryType* plotsTau_genDM10_byChargedIsolationTight               = new plotEntryType(Form("%s_tau_gendm10_byChargedIsolationTight", name.data()));
  plotEntryMap["Tau_genDM10"]["recoTauJetmatching"]                     = plotsTau_genDM10_recoTauJetmatching;
  plotEntryMap["Tau_genDM10"]["byDecayModeFindingNewDMs"]               = plotsTau_genDM10_byDecayModeFindingNewDMs;
  plotEntryMap["Tau_genDM10"]["byDecayModeFindingOldDMs"]               = plotsTau_genDM10_byDecayModeFindingOldDMs;
  plotEntryMap["Tau_genDM10"]["byCombinedIsolation8HitsLoose"]          = plotsTau_genDM10_byCombinedIsolation8HitsLoose;
  plotEntryMap["Tau_genDM10"]["byCombinedIsolation8HitsMedium"]         = plotsTau_genDM10_byCombinedIsolation8HitsMedium;
  plotEntryMap["Tau_genDM10"]["byCombinedIsolation8HitsTight"]          = plotsTau_genDM10_byCombinedIsolation8HitsTight;
  plotEntryMap["Tau_genDM10"]["byCombinedIsolation3HitsLoose"]          = plotsTau_genDM10_byCombinedIsolation3HitsLoose;
  plotEntryMap["Tau_genDM10"]["byCombinedIsolation3HitsMedium"]         = plotsTau_genDM10_byCombinedIsolation3HitsMedium;
  plotEntryMap["Tau_genDM10"]["byCombinedIsolation3HitsTight"]          = plotsTau_genDM10_byCombinedIsolation3HitsTight;
  plotEntryMap["Tau_genDM10"]["byChargedIsolationLoose"]                = plotsTau_genDM10_byChargedIsolationLoose;
  plotEntryMap["Tau_genDM10"]["byChargedIsolationMedium"]               = plotsTau_genDM10_byChargedIsolationMedium;
  plotEntryMap["Tau_genDM10"]["byChargedIsolationTight"]                = plotsTau_genDM10_byChargedIsolationTight;
  */
  plotEntryType* plotsTau_againstElectronVLooseMVA6             = new plotEntryType(Form("%s_tau_againstElectronVLooseMVA6", name.data()));
  plotEntryType* plotsTau_againstElectronLooseMVA6              = new plotEntryType(Form("%s_tau_againstElectronLooseMVA6", name.data()));
  plotEntryType* plotsTau_againstElectronMediumMVA6             = new plotEntryType(Form("%s_tau_againstElectronMediumMVA6", name.data()));
  plotEntryType* plotsTau_againstElectronTightMVA6              = new plotEntryType(Form("%s_tau_againstElectronTightMVA6", name.data()));
  plotEntryType* plotsTau_againstElectronVTightMVA6             = new plotEntryType(Form("%s_tau_againstElectronVTightMVA6", name.data()));
  plotEntryMap["Tau"]["againstElectronLooseMVA6"]               = plotsTau_againstElectronLooseMVA6;
  plotEntryMap["Tau"]["againstElectronMediumMVA6"]              = plotsTau_againstElectronMediumMVA6;
  plotEntryMap["Tau"]["againstElectronTightMVA6"]               = plotsTau_againstElectronTightMVA6;
  plotEntryMap["Tau"]["againstElectronVTightMVA6"]              = plotsTau_againstElectronVTightMVA6;

  plotEntryType* plotsTau_againstMuonLoose3                     = new plotEntryType(Form("%s_tau_againstMuonLoose3", name.data()));
  plotEntryType* plotsTau_againstMuonTight3                     = new plotEntryType(Form("%s_tau_againstMuonTight3", name.data()));
  plotEntryMap["Tau"]["againstMuonLoose3"]                      = plotsTau_againstMuonLoose3;
  plotEntryMap["Tau"]["againstMuonTight3"]                      = plotsTau_againstMuonTight3;

  plotEntryType* plotsElectron_againstElectronVLooseMVA6        = new plotEntryType(Form("%s_e_againstElectronVLooseMVA6", name.data()));
  plotEntryType* plotsElectron_againstElectronLooseMVA6         = new plotEntryType(Form("%s_e_againstElectronLooseMVA6", name.data()));
  plotEntryType* plotsElectron_againstElectronMediumMVA6        = new plotEntryType(Form("%s_e_againstElectronMediumMVA6", name.data()));
  plotEntryType* plotsElectron_againstElectronTightMVA6         = new plotEntryType(Form("%s_e_againstElectronTightMVA6", name.data()));
  plotEntryType* plotsElectron_againstElectronVTightMVA6        = new plotEntryType(Form("%s_e_againstElectronVTightMVA6", name.data()));
  plotEntryMap["Electron"]["againstElectronLooseMVA6"]          = plotsElectron_againstElectronLooseMVA6;
  plotEntryMap["Electron"]["againstElectronMediumMVA6"]         = plotsElectron_againstElectronMediumMVA6;
  plotEntryMap["Electron"]["againstElectronTightMVA6"]          = plotsElectron_againstElectronTightMVA6;
  plotEntryMap["Electron"]["againstElectronVTightMVA6"]         = plotsElectron_againstElectronVTightMVA6;
  
  plotEntryType* plotsMuon_againstMuonLoose3                    = new plotEntryType(Form("%s_mu_againstMuonLoose3", name.data()));
  plotEntryType* plotsMuon_againstMuonTight3                    = new plotEntryType(Form("%s_mu_againstMuonTight3", name.data()));
  plotEntryMap["Muon"]["againstMuonLoose3"]                     = plotsMuon_againstMuonLoose3;
  plotEntryMap["Muon"]["againstMuonTight3"]                     = plotsMuon_againstMuonTight3;
  
  plotEntryType* plotsJet_byDecayModeFindingNewDMs              = new plotEntryType(Form("%s_jet_byDecayModeFindingNewDMs", name.data()));
  plotEntryType* plotsJet_byDecayModeFindingOldDMs              = new plotEntryType(Form("%s_jet_byDecayModeFindingOldDMs", name.data()));
  plotEntryType* plotsJet_byCombinedIsolation3HitsLoose         = new plotEntryType(Form("%s_jet_byCombinedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsJet_byCombinedIsolation3HitsMedium        = new plotEntryType(Form("%s_jet_byCombinedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsJet_byCombinedIsolation3HitsTight         = new plotEntryType(Form("%s_jet_byCombinedIsolation3HitsTight", name.data()));
  plotEntryType* plotsJet_byCombinedIsolationNewDM3HitsLoose         = new plotEntryType(Form("%s_jet_byCombinedIsolationNewDM3HitsLoose", name.data()));
  plotEntryType* plotsJet_byCombinedIsolationNewDM3HitsMedium        = new plotEntryType(Form("%s_jet_byCombinedIsolationNewDM3HitsMedium", name.data()));
  plotEntryType* plotsJet_byCombinedIsolationNewDM3HitsTight         = new plotEntryType(Form("%s_jet_byCombinedIsolationNewDM3HitsTight", name.data()));
  plotEntryType* plotsJet_byCombinedIsolation3HitsdR03Loose         = new plotEntryType(Form("%s_jet_byCombinedIsolation3HitsdR03Loose", name.data()));
  plotEntryType* plotsJet_byCombinedIsolation3HitsdR03Medium        = new plotEntryType(Form("%s_jet_byCombinedIsolation3HitsdR03Medium", name.data()));
  plotEntryType* plotsJet_byCombinedIsolation3HitsdR03Tight         = new plotEntryType(Form("%s_jet_byCombinedIsolation3HitsdR03Tight", name.data()));
  plotEntryType* plotsJet_byChargedIsolationLoose               = new plotEntryType(Form("%s_jet_byChargedIsolationLoose", name.data()));
  plotEntryType* plotsJet_byChargedIsolationMedium              = new plotEntryType(Form("%s_jet_byChargedIsolationMedium", name.data()));
  plotEntryType* plotsJet_byChargedIsolationTight               = new plotEntryType(Form("%s_jet_byChargedIsolationTight", name.data()));

  plotEntryType* plotsJet_byPileupWeightedIsolation3HitsLoose         = new plotEntryType(Form("%s_jet_byPileupWeightedIsolation3HitsLoose", name.data()));
  plotEntryType* plotsJet_byPileupWeightedIsolation3HitsMedium        = new plotEntryType(Form("%s_jet_byPileupWeightedIsolation3HitsMedium", name.data()));
  plotEntryType* plotsJet_byPileupWeightedIsolation3HitsTight         = new plotEntryType(Form("%s_jet_byPileupWeightedIsolation3HitsTight", name.data()));
  plotEntryType* plotsJet_byPileupWeightedIsolationNewDM3HitsLoose         = new plotEntryType(Form("%s_jet_byPileupWeightedIsolationNewDM3HitsLoose", name.data()));
  plotEntryType* plotsJet_byPileupWeightedIsolationNewDM3HitsMedium        = new plotEntryType(Form("%s_jet_byPileupWeightedIsolationNewDM3HitsMedium", name.data()));
  plotEntryType* plotsJet_byPileupWeightedIsolationNewDM3HitsTight         = new plotEntryType(Form("%s_jet_byPileupWeightedIsolationNewDM3HitsTight", name.data()));
  plotEntryType* plotsJet_byVLooseIsolationMVArun2v1DBoldDMwLT         = new plotEntryType(Form("%s_jet_byVLooseIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsJet_byLooseIsolationMVArun2v1DBoldDMwLT          = new plotEntryType(Form("%s_jet_byLooseIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsJet_byMediumIsolationMVArun2v1DBoldDMwLT         = new plotEntryType(Form("%s_jet_byMediumIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsJet_byTightIsolationMVArun2v1DBoldDMwLT          = new plotEntryType(Form("%s_jet_byTightIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsJet_byVTightIsolationMVArun2v1DBoldDMwLT         = new plotEntryType(Form("%s_jet_byVTightIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsJet_byVVTightIsolationMVArun2v1DBoldDMwLT        = new plotEntryType(Form("%s_jet_byVVTightIsolationMVArun2v1DBoldDMwLT", name.data()));
  plotEntryType* plotsJet_byVLooseIsolationMVArun2v1DBnewDMwLT         = new plotEntryType(Form("%s_jet_byVLooseIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsJet_byLooseIsolationMVArun2v1DBnewDMwLT          = new plotEntryType(Form("%s_jet_byLooseIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsJet_byMediumIsolationMVArun2v1DBnewDMwLT         = new plotEntryType(Form("%s_jet_byMediumIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsJet_byTightIsolationMVArun2v1DBnewDMwLT          = new plotEntryType(Form("%s_jet_byTightIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsJet_byVTightIsolationMVArun2v1DBnewDMwLT         = new plotEntryType(Form("%s_jet_byVTightIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsJet_byVVTightIsolationMVArun2v1DBnewDMwLT        = new plotEntryType(Form("%s_jet_byVVTightIsolationMVArun2v1DBnewDMwLT", name.data()));
  plotEntryType* plotsJet_byVLooseIsolationMVArun2v1PWoldDMwLT         = new plotEntryType(Form("%s_jet_byVLooseIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsJet_byLooseIsolationMVArun2v1PWoldDMwLT          = new plotEntryType(Form("%s_jet_byLooseIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsJet_byMediumIsolationMVArun2v1PWoldDMwLT         = new plotEntryType(Form("%s_jet_byMediumIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsJet_byTightIsolationMVArun2v1PWoldDMwLT          = new plotEntryType(Form("%s_jet_byTightIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsJet_byVTightIsolationMVArun2v1PWoldDMwLT         = new plotEntryType(Form("%s_jet_byVTightIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsJet_byVVTightIsolationMVArun2v1PWoldDMwLT        = new plotEntryType(Form("%s_jet_byVVTightIsolationMVArun2v1PWoldDMwLT", name.data()));
  plotEntryType* plotsJet_byVLooseIsolationMVArun2v1PWnewDMwLT         = new plotEntryType(Form("%s_jet_byVLooseIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsJet_byLooseIsolationMVArun2v1PWnewDMwLT          = new plotEntryType(Form("%s_jet_byLooseIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsJet_byMediumIsolationMVArun2v1PWnewDMwLT         = new plotEntryType(Form("%s_jet_byMediumIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsJet_byTightIsolationMVArun2v1PWnewDMwLT          = new plotEntryType(Form("%s_jet_byTightIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsJet_byVTightIsolationMVArun2v1PWnewDMwLT         = new plotEntryType(Form("%s_jet_byVTightIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsJet_byVVTightIsolationMVArun2v1PWnewDMwLT        = new plotEntryType(Form("%s_jet_byVVTightIsolationMVArun2v1PWnewDMwLT", name.data()));
  plotEntryType* plotsJet_byVLooseIsolationMVArun2v1DBdR03oldDMwLT         = new plotEntryType(Form("%s_jet_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byLooseIsolationMVArun2v1DBdR03oldDMwLT          = new plotEntryType(Form("%s_jet_byLooseIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byMediumIsolationMVArun2v1DBdR03oldDMwLT         = new plotEntryType(Form("%s_jet_byMediumIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byTightIsolationMVArun2v1DBdR03oldDMwLT          = new plotEntryType(Form("%s_jet_byTightIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byVTightIsolationMVArun2v1DBdR03oldDMwLT         = new plotEntryType(Form("%s_jet_byVTightIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byVVTightIsolationMVArun2v1DBdR03oldDMwLT        = new plotEntryType(Form("%s_jet_byVVTightIsolationMVArun2v1DBdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byVLooseIsolationMVArun2v1PWdR03oldDMwLT         = new plotEntryType(Form("%s_jet_byVLooseIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byLooseIsolationMVArun2v1PWdR03oldDMwLT          = new plotEntryType(Form("%s_jet_byLooseIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byMediumIsolationMVArun2v1PWdR03oldDMwLT         = new plotEntryType(Form("%s_jet_byMediumIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byTightIsolationMVArun2v1PWdR03oldDMwLT          = new plotEntryType(Form("%s_jet_byTightIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byVTightIsolationMVArun2v1PWdR03oldDMwLT         = new plotEntryType(Form("%s_jet_byVTightIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryType* plotsJet_byVVTightIsolationMVArun2v1PWdR03oldDMwLT        = new plotEntryType(Form("%s_jet_byVVTightIsolationMVArun2v1PWdR03oldDMwLT", name.data()));
  plotEntryMap["Jet"]["byDecayModeFindingNewDMs"]               = plotsJet_byDecayModeFindingNewDMs;
  plotEntryMap["Jet"]["byDecayModeFindingOldDMs"]               = plotsJet_byDecayModeFindingOldDMs;
  plotEntryMap["Jet"]["byCombinedIsolation3HitsLoose"]          = plotsJet_byCombinedIsolation3HitsLoose;
  plotEntryMap["Jet"]["byCombinedIsolation3HitsMedium"]         = plotsJet_byCombinedIsolation3HitsMedium;
  plotEntryMap["Jet"]["byCombinedIsolation3HitsTight"]          = plotsJet_byCombinedIsolation3HitsTight;
  plotEntryMap["Jet"]["byCombinedIsolationNewDM3HitsLoose"]          = plotsJet_byCombinedIsolationNewDM3HitsLoose;
  plotEntryMap["Jet"]["byCombinedIsolationNewDM3HitsMedium"]         = plotsJet_byCombinedIsolationNewDM3HitsMedium;
  plotEntryMap["Jet"]["byCombinedIsolationNewDM3HitsTight"]          = plotsJet_byCombinedIsolationNewDM3HitsTight;
  plotEntryMap["Jet"]["byCombinedIsolation3HitsdR03Loose"]          = plotsJet_byCombinedIsolation3HitsdR03Loose;
  plotEntryMap["Jet"]["byCombinedIsolation3HitsdR03Medium"]         = plotsJet_byCombinedIsolation3HitsdR03Medium;
  plotEntryMap["Jet"]["byCombinedIsolation3HitsdR03Tight"]          = plotsJet_byCombinedIsolation3HitsdR03Tight;
  plotEntryMap["Jet"]["byChargedIsolationLoose"]                = plotsJet_byChargedIsolationLoose;
  plotEntryMap["Jet"]["byChargedIsolationMedium"]               = plotsJet_byChargedIsolationMedium;
  plotEntryMap["Jet"]["byChargedIsolationTight"]                = plotsJet_byChargedIsolationTight;
  plotEntryMap["Jet"]["byPileupWeightedIsolation3HitsLoose"]          = plotsJet_byPileupWeightedIsolation3HitsLoose;
  plotEntryMap["Jet"]["byPileupWeightedIsolation3HitsMedium"]         = plotsJet_byPileupWeightedIsolation3HitsMedium;
  plotEntryMap["Jet"]["byPileupWeightedIsolation3HitsTight"]          = plotsJet_byPileupWeightedIsolation3HitsTight;
  plotEntryMap["Jet"]["byPileupWeightedIsolationNewDM3HitsLoose"]          = plotsJet_byPileupWeightedIsolationNewDM3HitsLoose;
  plotEntryMap["Jet"]["byPileupWeightedIsolationNewDM3HitsMedium"]         = plotsJet_byPileupWeightedIsolationNewDM3HitsMedium;
  plotEntryMap["Jet"]["byPileupWeightedIsolationNewDM3HitsTight"]          = plotsJet_byPileupWeightedIsolationNewDM3HitsTight;
  plotEntryMap["Jet"]["byLooseIsolationMVArun2v1DBoldDMwLT"]           = plotsJet_byLooseIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Jet"]["byMediumIsolationMVArun2v1DBoldDMwLT"]          = plotsJet_byMediumIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Jet"]["byTightIsolationMVArun2v1DBoldDMwLT"]           = plotsJet_byTightIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Jet"]["byVTightIsolationMVArun2v1DBoldDMwLT"]          = plotsJet_byVTightIsolationMVArun2v1DBoldDMwLT;
  plotEntryMap["Jet"]["byLooseIsolationMVArun2v1DBnewDMwLT"]           = plotsJet_byLooseIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Jet"]["byMediumIsolationMVArun2v1DBnewDMwLT"]          = plotsJet_byMediumIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Jet"]["byTightIsolationMVArun2v1DBnewDMwLT"]           = plotsJet_byTightIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Jet"]["byVTightIsolationMVArun2v1DBnewDMwLT"]          = plotsJet_byVTightIsolationMVArun2v1DBnewDMwLT;
  plotEntryMap["Jet"]["byLooseIsolationMVArun2v1PWoldDMwLT"]           = plotsJet_byLooseIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Jet"]["byMediumIsolationMVArun2v1PWoldDMwLT"]          = plotsJet_byMediumIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Jet"]["byTightIsolationMVArun2v1PWoldDMwLT"]           = plotsJet_byTightIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Jet"]["byVTightIsolationMVArun2v1PWoldDMwLT"]          = plotsJet_byVTightIsolationMVArun2v1PWoldDMwLT;
  plotEntryMap["Jet"]["byLooseIsolationMVArun2v1PWnewDMwLT"]           = plotsJet_byLooseIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Jet"]["byMediumIsolationMVArun2v1PWnewDMwLT"]          = plotsJet_byMediumIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Jet"]["byTightIsolationMVArun2v1PWnewDMwLT"]           = plotsJet_byTightIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Jet"]["byVTightIsolationMVArun2v1PWnewDMwLT"]          = plotsJet_byVTightIsolationMVArun2v1PWnewDMwLT;
  plotEntryMap["Jet"]["byLooseIsolationMVArun2v1DBdR03oldDMwLT"]           = plotsJet_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Jet"]["byMediumIsolationMVArun2v1DBdR03oldDMwLT"]          = plotsJet_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Jet"]["byTightIsolationMVArun2v1DBdR03oldDMwLT"]           = plotsJet_byTightIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Jet"]["byVTightIsolationMVArun2v1DBdR03oldDMwLT"]          = plotsJet_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
  plotEntryMap["Jet"]["byLooseIsolationMVArun2v1PWdR03oldDMwLT"]           = plotsJet_byLooseIsolationMVArun2v1PWdR03oldDMwLT;
  plotEntryMap["Jet"]["byMediumIsolationMVArun2v1PWdR03oldDMwLT"]          = plotsJet_byMediumIsolationMVArun2v1PWdR03oldDMwLT;
  plotEntryMap["Jet"]["byTightIsolationMVArun2v1PWdR03oldDMwLT"]           = plotsJet_byTightIsolationMVArun2v1PWdR03oldDMwLT;
  plotEntryMap["Jet"]["byVTightIsolationMVArun2v1PWdR03oldDMwLT"]          = plotsJet_byVTightIsolationMVArun2v1PWdR03oldDMwLT;

  int numElectrons = 0;
  int numMuons     = 0;
  int numTaus      = 0;
  int numJets      = 0;

  int numEntries = tree->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries && (iEntry < maxEvents || maxEvents == -1); ++iEntry ) {
    if ( iEntry > 0 && (iEntry % 100000) == 0 ) {
      std::cout << "processing Event " << iEntry 
		<< " (numTaus = " << numTaus << ", numElectrons = " << numElectrons << ", numMuons = " << numMuons << ", numJets = " << numJets << ")" << std::endl;
    }

    tree->GetEntry(iEntry);

    TLorentzVector recTauP4(recTauPx, recTauPy, recTauPz, recTauEn);
    Float_t recTauPt  = recTauP4.Pt();
    Float_t recTauEta = ( recTauPt > 1.e-3 ) ? recTauP4.Eta() : 9.9;
    Float_t recTauPhi = recTauP4.Phi();
    
    TLorentzVector leadPFCandP4(leadPFCandPx, leadPFCandPy, leadPFCandPz, leadPFCandEn);
    TLorentzVector leadPFChargedHadrCandP4(leadPFChargedHadrCandPx, leadPFChargedHadrCandPy, leadPFChargedHadrCandPz, leadPFChargedHadrCandEn);
    Float_t leadPFChargedHadrCandPt = leadPFChargedHadrCandP4.Pt();
    Float_t leadPFChargedHadrCandEta = ( leadPFChargedHadrCandPt > 1.e-3 ) ? leadPFChargedHadrCandP4.Eta() : 9.9;
    Float_t tauNeutralHadPt  = (recTauP4 - leadPFChargedHadrCandP4).Pt();

    if ( mode == kOldTags ) {
      decayModeFindingNewDMs = decayModeFindingOldDMs;
    }

    TLorentzVector genTauP4(genTauPx, genTauPy, genTauPz, genTauEn);
    Float_t genTauPt  = genTauP4.Pt();
    Float_t genTauEta = ( genTauPt > 1.e-3 ) ? genTauP4.Eta() : 9.9;
    Float_t genTauPhi = genTauP4.Phi();

    TLorentzVector genTauLeadChHadP4(genTauLeadChHadPx, genTauLeadChHadPy, genTauLeadChHadPz, genTauLeadChHadEn);
    Float_t genTauLeadChHadPt  = genTauLeadChHadP4.Pt();
    Float_t genTauLeadChHadEta = ( genTauLeadChHadPt > 1.e-3 ) ? genTauLeadChHadP4.Eta() : 9.9;
    //Float_t genTauLeadChHadPhi = genTauLeadChHadP4.Phi();
    Float_t genTauNeuHadPt  = (genTauP4 - genTauLeadChHadP4).Pt();

    TLorentzVector genElectronP4(genElectronPx, genElectronPy, genElectronPz, genElectronEn);
    TLorentzVector genMuonP4(genMuonPx, genMuonPy, genMuonPz, genMuonEn);

    //TLorentzVector genQuarkOrGluonP4(genQuarkOrGluonPx, genQuarkOrGluonPy, genQuarkOrGluonPz, genQuarkOrGluonEn);
    //Float_t genQuarkOrGluonPt  = genQuarkOrGluonP4.Pt();
    //Float_t genQuarkOrGluonEta = ( genQuarkOrGluonPt > 1.e-3 ) ? genQuarkOrGluonP4.Eta() : 9.9;
    //Float_t genQuarkOrGluonPhi = genQuarkOrGluonP4.Phi();
    
    TLorentzVector recJetP4(recJetPx, recJetPy, recJetPz, recJetEn);
    Float_t recJetPt  = recJetP4.Pt();
    Float_t recJetEta = ( recJetPt > 1.e-3 ) ? recJetP4.Eta() : 9.9;
    Float_t recJetPhi = recJetP4.Phi();

    evtWeight_dummy = 1.;
    /*
    if ( recTauPt > 20. && 
	 TMath::Abs(recTauEta) < 2.3 && 
	 //decayModeFindingNewDMs > 0.5 && 
	 decayModeFindingOldDMs > 0.5 && 
	 byLooseCombinedIsolationDeltaBetaCorr > 0.5 ) {
      if ( genTauMatch > 0.5 && genTauDeltaR < 0.3 && 
	   (genTauDecayMode ==  0 || 
	    genTauDecayMode ==  1 || 
	    genTauDecayMode ==  2 || 
	    genTauDecayMode == 10) ) {	
	if ( decayModeFindingOldDMs > 0.5 ) {
  	  plotsTau_againstElectronLoose->fillHistograms(
	    againstElectronLoose > 0.5, true, recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	  plotsTau_againstElectronMedium->fillHistograms(
            againstElectronMedium > 0.5, true, recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
          plotsTau_againstElectronTight->fillHistograms(
            againstElectronTight > 0.5, true, recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	  plotsTau_againstElectronLooseMVA3->fillHistograms(
            againstElectronLooseMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
          plotsTau_againstElectronMediumMVA3->fillHistograms(
            againstElectronMediumMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
          plotsTau_againstElectronTightMVA3->fillHistograms(
            againstElectronTightMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
          plotsTau_againstElectronVTightMVA3->fillHistograms(
            againstElectronVTightMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
	  plotsTau_againstElectronVLooseMVA6->fillHistograms(
            againstElectronVLooseMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
	  plotsTau_againstElectronLooseMVA6->fillHistograms(
            againstElectronLooseMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
          plotsTau_againstElectronMediumMVA6->fillHistograms(
            againstElectronMediumMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
          plotsTau_againstElectronTightMVA6->fillHistograms(
            againstElectronTightMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
          plotsTau_againstElectronVTightMVA6->fillHistograms(
            againstElectronVTightMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
        }

        plotsTau_againstMuonLoose->fillHistograms(
	  againstMuonLoose > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonMedium->fillHistograms(
          againstMuonMedium > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonTight->fillHistograms(
          againstMuonTight > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsTau_againstMuonLoose2->fillHistograms(
	  againstMuonLoose2 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonMedium2->fillHistograms(
          againstMuonMedium2 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonTight2->fillHistograms(
          againstMuonTight2 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonLoose3->fillHistograms(
	  againstMuonLoose3 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonTight3->fillHistograms(
          againstMuonTight3 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonLooseMVA->fillHistograms(
	  againstMuonLooseMVA > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonMediumMVA->fillHistograms(
          againstMuonMediumMVA > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsTau_againstMuonTightMVA->fillHistograms(
          againstMuonTightMVA > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);

        ++numTaus;
      }
      if ( genElectronMatch > 0.5 && genElectronDeltaR < 0.3  ) {
	if ( decayModeFindingOldDMs > 0.5 ) {
	  plotsElectron_againstElectronLoose->fillHistograms(
            againstElectronLoose > 0.5, true, recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	  plotsElectron_againstElectronMedium->fillHistograms(
            againstElectronMedium > 0.5, true, recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
          plotsElectron_againstElectronTight->fillHistograms(
	    againstElectronTight > 0.5, true, recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	  plotsElectron_againstElectronLooseMVA3->fillHistograms(
            againstElectronLooseMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
	  plotsElectron_againstElectronMediumMVA3->fillHistograms(
            againstElectronMediumMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
          plotsElectron_againstElectronTightMVA3->fillHistograms(
	    againstElectronTightMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
          plotsElectron_againstElectronVTightMVA3->fillHistograms(
	    againstElectronVTightMVA3 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA3category), numVertices, evtWeight_dummy);
	  plotsElectron_againstElectronVLooseMVA6->fillHistograms(
            againstElectronVLooseMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
	  plotsElectron_againstElectronLooseMVA6->fillHistograms(
            againstElectronLooseMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
	  plotsElectron_againstElectronMediumMVA6->fillHistograms(
            againstElectronMediumMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
          plotsElectron_againstElectronTightMVA6->fillHistograms(
	    againstElectronTightMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
          plotsElectron_againstElectronVTightMVA6->fillHistograms(
	    againstElectronVTightMVA6 > 0.5, true, 
	    recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, TMath::Nint(againstElectronMVA6category), numVertices, evtWeight_dummy);
	}

	++numElectrons;
      }
      if ( genMuonMatch > 0.5 && genMuonDeltaR < 0.3  ) {
	plotsMuon_againstMuonLoose->fillHistograms(
	  againstMuonLoose > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsMuon_againstMuonMedium->fillHistograms(
          againstMuonMedium > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsMuon_againstMuonTight->fillHistograms(
          againstMuonTight > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsMuon_againstMuonLoose2->fillHistograms(
	  againstMuonLoose2 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsMuon_againstMuonMedium2->fillHistograms(
          againstMuonMedium2 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsMuon_againstMuonTight2->fillHistograms(
          againstMuonTight2 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsMuon_againstMuonLoose3->fillHistograms(
	  againstMuonLoose3 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsMuon_againstMuonTight3->fillHistograms(
          againstMuonTight3 > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsMuon_againstMuonLooseMVA->fillHistograms(
	  againstMuonLooseMVA > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
	plotsMuon_againstMuonMediumMVA->fillHistograms(
          againstMuonMediumMVA > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
        plotsMuon_againstMuonTightMVA->fillHistograms(
          againstMuonTightMVA > 0.5, true, 
	  recTauPt, recTauEta, recTauEtaImpactECAL, recTauPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);

	++numMuons;
      }
    }
    */
    if ( genTauPt > 20. && TMath::Abs(genTauEta) < 2.3 &&
    	 ((genTauDecayMode ==  0 ) ||
    	  (genTauDecayMode ==  1 ) ||
    	  (genTauDecayMode ==  2 ) ||
    	  (genTauDecayMode == 10 ) ) //&&
	 //genTauMatch > 0.5 && genTauDeltaR < 0.3 && recTauPt > 18.) {
	 ){
      bool recTauPreselection_passed = 
       (genTauMatch > 0.5 && genTauDeltaR < 0.3 && 
	 recTauPt > 20. && 
	 TMath::Abs(recTauEta) < 2.3);
      bool recTauOldDecayMode_passed = decayModeFindingOldDMs > 0.5 &&
	(recTauDecayMode == -1 || recTauDecayMode == 0 || recTauDecayMode == 1 || recTauDecayMode == 2 || recTauDecayMode == 10);
      bool recTauNewDecayMode_passed = decayModeFindingNewDMs > 0.5 &&
        (recTauDecayMode == -1 || recTauDecayMode == 0 || recTauDecayMode == 1 || recTauDecayMode == 2 || recTauDecayMode == 5 || recTauDecayMode == 6 || recTauDecayMode == 10);
      //bool recTauNewDecayMode_passed = recTauOldDecayMode_passed;

      double genTauEtaImpactECAL = genTauEta;
      
      plotsTau_recoTauJetmatching->fillHistograms(
        recTauPreselection_passed, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt, genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byDecayModeFindingNewDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingNewDMs > 0.5, true,  
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byDecayModeFindingOldDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingOldDMs > 0.5, true,  
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      /*
      plotsTau_genDM0_recoTauJetmatching->fillHistograms(
        recTauPreselection_passed && genTauDecayMode ==  0, genTauDecayMode ==  0,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM0_byDecayModeFindingNewDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingNewDMs > 0.5 && genTauDecayMode ==  0, genTauDecayMode ==  0,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM0_byDecayModeFindingOldDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingOldDMs > 0.5 && genTauDecayMode ==  0, genTauDecayMode ==  0,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM1_recoTauJetmatching->fillHistograms(
	recTauPreselection_passed && genTauDecayMode ==  1, genTauDecayMode ==  1,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM1_byDecayModeFindingNewDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingNewDMs > 0.5 && genTauDecayMode ==  1, genTauDecayMode ==  1,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM1_byDecayModeFindingOldDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingOldDMs > 0.5 && genTauDecayMode ==  1, genTauDecayMode ==  1,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM2_recoTauJetmatching->fillHistograms(
	recTauPreselection_passed && genTauDecayMode ==  2, genTauDecayMode ==  2,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM2_byDecayModeFindingNewDMs->fillHistograms(
        recTauPreselection_passed && decayModeFindingNewDMs > 0.5 && genTauDecayMode ==  2, genTauDecayMode ==  2,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM2_byDecayModeFindingOldDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingOldDMs > 0.5 && genTauDecayMode ==  2, genTauDecayMode ==  2,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM10_recoTauJetmatching->fillHistograms(
	recTauPreselection_passed && genTauDecayMode ==  10, genTauDecayMode ==  10,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM10_byDecayModeFindingNewDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingNewDMs > 0.5 && genTauDecayMode ==  10, genTauDecayMode ==  10,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM10_byDecayModeFindingOldDMs->fillHistograms(
	recTauPreselection_passed && decayModeFindingOldDMs > 0.5 && genTauDecayMode ==  10, genTauDecayMode ==  10,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      */

      plotsTau_byCombinedIsolation3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,  
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byCombinedIsolation3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,  
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byCombinedIsolation3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,  
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byCombinedIsolationNewDM3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauNewDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byCombinedIsolationNewDM3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauNewDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byCombinedIsolationNewDM3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauNewDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byCombinedIsolation3HitsdR03Loose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byCombinedIsolation3HitsdR03Medium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byCombinedIsolation3HitsdR03Tight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3HitsdR03 > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);


      plotsTau_byPileupWeightedIsolation3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLoosePileupWeightedIsolation3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byPileupWeightedIsolation3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumPileupWeightedIsolation3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byPileupWeightedIsolation3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightPileupWeightedIsolation3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byPileupWeightedIsolationNewDM3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauNewDecayMode_passed && byLoosePileupWeightedIsolation3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byPileupWeightedIsolationNewDM3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauNewDecayMode_passed && byMediumPileupWeightedIsolation3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byPileupWeightedIsolationNewDM3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauNewDecayMode_passed && byTightPileupWeightedIsolation3Hits > 0.5, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);


      /*
      plotsTau_genDM0_byCombinedIsolation3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  0, genTauDecayMode ==  0 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM0_byCombinedIsolation3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  0, genTauDecayMode ==  0 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM0_byCombinedIsolation3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  0, genTauDecayMode ==  0 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM1_byCombinedIsolation3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  1, genTauDecayMode ==  1 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM1_byCombinedIsolation3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  1, genTauDecayMode ==  1 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM1_byCombinedIsolation3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  1, genTauDecayMode ==  1 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM2_byCombinedIsolation3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  2, genTauDecayMode ==  2 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM2_byCombinedIsolation3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  2, genTauDecayMode ==  2 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM2_byCombinedIsolation3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  2, genTauDecayMode ==  2 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM10_byCombinedIsolation3HitsLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  10, genTauDecayMode ==  10 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM10_byCombinedIsolation3HitsMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  10, genTauDecayMode ==  10 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM10_byCombinedIsolation3HitsTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5 && genTauDecayMode ==  10, genTauDecayMode ==  10 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      */
      plotsTau_byChargedIsolationLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 2.0, true,  
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byChargedIsolationMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 1.0, true,  
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byChargedIsolationTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 0.8, true,  
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      /*
      plotsTau_genDM0_byChargedIsolationLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 2.0 && genTauDecayMode ==  0, genTauDecayMode ==  0 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM0_byChargedIsolationMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 1.0 && genTauDecayMode ==  0, genTauDecayMode ==  0 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM0_byChargedIsolationTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 0.8 && genTauDecayMode ==  0, genTauDecayMode ==  0 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM1_byChargedIsolationLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 2.0 && genTauDecayMode ==  1, genTauDecayMode ==  1 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM1_byChargedIsolationMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 1.0 && genTauDecayMode ==  1, genTauDecayMode ==  1 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM1_byChargedIsolationTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 0.8 && genTauDecayMode ==  1, genTauDecayMode ==  1 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM2_byChargedIsolationLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 2.0 && genTauDecayMode ==  2, genTauDecayMode ==  2 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM2_byChargedIsolationMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 1.0 && genTauDecayMode ==  2, genTauDecayMode ==  2 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM2_byChargedIsolationTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 0.8 && genTauDecayMode ==  2, genTauDecayMode ==  2 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_genDM10_byChargedIsolationLoose->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 2.0 && genTauDecayMode ==  10, genTauDecayMode ==  10 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM10_byChargedIsolationMedium->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 1.0 && genTauDecayMode ==  10, genTauDecayMode ==  10 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_genDM10_byChargedIsolationTight->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && chargedIsoPtSum < 0.8 && genTauDecayMode ==  10, genTauDecayMode ==  10 && recTauPreselection_passed && recTauOldDecayMode_passed,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      */

      plotsTau_byVLooseIsolationMVArun2v1DBoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVLooseIsolationMVArun2v1DBoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byLooseIsolationMVArun2v1DBoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byLooseIsolationMVArun2v1DBoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byMediumIsolationMVArun2v1DBoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byMediumIsolationMVArun2v1DBoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byTightIsolationMVArun2v1DBoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byTightIsolationMVArun2v1DBoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVTightIsolationMVArun2v1DBoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVTightIsolationMVArun2v1DBoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVVTightIsolationMVArun2v1DBoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVVTightIsolationMVArun2v1DBoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byVLooseIsolationMVArun2v1DBnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byVLooseIsolationMVArun2v1DBnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byLooseIsolationMVArun2v1DBnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byLooseIsolationMVArun2v1DBnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byMediumIsolationMVArun2v1DBnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byMediumIsolationMVArun2v1DBnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byTightIsolationMVArun2v1DBnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byTightIsolationMVArun2v1DBnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVTightIsolationMVArun2v1DBnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byVTightIsolationMVArun2v1DBnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVVTightIsolationMVArun2v1DBnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byVVTightIsolationMVArun2v1DBnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byVLooseIsolationMVArun2v1PWoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVLooseIsolationMVArun2v1PWoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byLooseIsolationMVArun2v1PWoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byLooseIsolationMVArun2v1PWoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byMediumIsolationMVArun2v1PWoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byMediumIsolationMVArun2v1PWoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byTightIsolationMVArun2v1PWoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byTightIsolationMVArun2v1PWoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVTightIsolationMVArun2v1PWoldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVTightIsolationMVArun2v1PWoldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy); 
      //plotsTau_byVVTightIsolationMVArun2v1PWoldDMwLT->fillHistograms(
      //  recTauPreselection_passed && recTauOldDecayMode_passed && byVVTightIsolationMVArun2v1PWoldDMwLT, true,
      //  genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byVLooseIsolationMVArun2v1PWnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byVLooseIsolationMVArun2v1PWnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byLooseIsolationMVArun2v1PWnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byLooseIsolationMVArun2v1PWnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byMediumIsolationMVArun2v1PWnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byMediumIsolationMVArun2v1PWnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byTightIsolationMVArun2v1PWnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byTightIsolationMVArun2v1PWnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVTightIsolationMVArun2v1PWnewDMwLT->fillHistograms(
        recTauPreselection_passed && recTauNewDecayMode_passed && byVTightIsolationMVArun2v1PWnewDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      //plotsTau_byVVTightIsolationMVArun2v1PWnewDMwLT->fillHistograms(
      //  recTauPreselection_passed && recTauNewDecayMode_passed && byVVTightIsolationMVArun2v1PWnewDMwLT, true,
      //  genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVLooseIsolationMVArun2v1DBdR03oldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byLooseIsolationMVArun2v1DBdR03oldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byMediumIsolationMVArun2v1DBdR03oldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byTightIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byTightIsolationMVArun2v1DBdR03oldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVTightIsolationMVArun2v1DBdR03oldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
        recTauPreselection_passed && recTauOldDecayMode_passed && byVVTightIsolationMVArun2v1DBdR03oldDMwLT, true,
        genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      plotsTau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
	 recTauPreselection_passed && recTauOldDecayMode_passed && byVLooseIsolationMVArun2v1PWdR03oldDMwLT, true,
	 genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byLooseIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byLooseIsolationMVArun2v1PWdR03oldDMwLT, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byMediumIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
	 recTauPreselection_passed && recTauOldDecayMode_passed && byMediumIsolationMVArun2v1PWdR03oldDMwLT, true,
	 genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byTightIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
	recTauPreselection_passed && recTauOldDecayMode_passed && byTightIsolationMVArun2v1PWdR03oldDMwLT, true,
	genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      plotsTau_byVTightIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
	 recTauPreselection_passed && recTauOldDecayMode_passed && byVTightIsolationMVArun2v1PWdR03oldDMwLT, true,
	 genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);
      //plotsTau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
      //  recTauPreselection_passed && recTauOldDecayMode_passed && byVVTightIsolationMVArun2v1PWdR03oldDMwLT, true,
      //  genTauPt, genTauEta, genTauEtaImpactECAL, genTauPhi, genTauLeadChHadPt, genTauLeadChHadEta, genTauNeuHadPt,genTauDecayMode, 0, numVertices, evtWeight_dummy);

      ++numTaus;
    }
    

//    //if ( genQuarkOrGluonMatch > 0.5 && genQuarkOrGluonDeltaR < 0.3 && genQuarkOrGluonPt > 15. && TMath::Abs(genQuarkOrGluonEta) < 2.3 ) {
//    if ( recJetPt > 20. && TMath::Abs(recJetEta) < 2.3  && ElectronVetoMatch <= 0) { 
// 
//      bool recJetPreselection_passed = 
//	(recTauPt > 20. && TMath::Abs(recTauEta) < 2.3 && recJetMatch > 0.5 && recJetDeltaR < 0.3);
//      bool recJetOldDecayMode_passed = decayModeFindingOldDMs > 0.5 &&
//	(recTauDecayMode == -1 || recTauDecayMode == 0 || recTauDecayMode == 1 || recTauDecayMode == 2 || recTauDecayMode == 10);
//      bool recJetNewDecayMode_passed = decayModeFindingNewDMs > 0.5 &&
//        (recTauDecayMode == -1 || recTauDecayMode == 0 || recTauDecayMode == 1 || recTauDecayMode == 2 || recTauDecayMode == 5 || recTauDecayMode == 6 || recTauDecayMode == 10);
//      //bool recJetNewDecayMode_passed = recJetOldDecayMode_passed;
//
//      double recJetEtaImpactECAL = recJetEta;
//     
//      plotsJet_byDecayModeFindingNewDMs->fillHistograms(
//	recJetPreselection_passed && decayModeFindingNewDMs > 0.5, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byDecayModeFindingOldDMs->fillHistograms(
//	recJetPreselection_passed && decayModeFindingOldDMs > 0.5, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byCombinedIsolation3HitsLoose->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byCombinedIsolation3HitsMedium->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byCombinedIsolation3HitsTight->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,  
//        recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byCombinedIsolationNewDM3HitsLoose->fillHistograms(
//	recJetPreselection_passed && recJetNewDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byCombinedIsolationNewDM3HitsMedium->fillHistograms(
//	recJetPreselection_passed && recJetNewDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byCombinedIsolationNewDM3HitsTight->fillHistograms(
//	recJetPreselection_passed && recJetNewDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byCombinedIsolation3HitsdR03Loose->fillHistograms(
//        recJetPreselection_passed && recJetOldDecayMode_passed && byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byCombinedIsolation3HitsdR03Medium->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byCombinedIsolation3HitsdR03Tight->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byTightCombinedIsolationDeltaBetaCorr3HitsdR03 > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      
//      plotsJet_byPileupWeightedIsolation3HitsLoose->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byLoosePileupWeightedIsolation3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byPileupWeightedIsolation3HitsMedium->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byMediumPileupWeightedIsolation3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byPileupWeightedIsolation3HitsTight->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && byTightPileupWeightedIsolation3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byPileupWeightedIsolationNewDM3HitsLoose->fillHistograms(
//	recJetPreselection_passed && recJetNewDecayMode_passed && byLoosePileupWeightedIsolation3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byPileupWeightedIsolationNewDM3HitsMedium->fillHistograms(
//	recJetPreselection_passed && recJetNewDecayMode_passed && byMediumPileupWeightedIsolation3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byPileupWeightedIsolationNewDM3HitsTight->fillHistograms(
//	recJetPreselection_passed && recJetNewDecayMode_passed && byTightPileupWeightedIsolation3Hits > 0.5, true,
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byChargedIsolationLoose->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && chargedIsoPtSum < 2.0, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byChargedIsolationMedium->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && chargedIsoPtSum < 1.0, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byChargedIsolationTight->fillHistograms(
//	recJetPreselection_passed && recJetOldDecayMode_passed && chargedIsoPtSum < 0.8, true,  
//	recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byVLooseIsolationMVArun2v1DBoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVLooseIsolationMVArun2v1DBoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byLooseIsolationMVArun2v1DBoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byLooseIsolationMVArun2v1DBoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byMediumIsolationMVArun2v1DBoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byMediumIsolationMVArun2v1DBoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byTightIsolationMVArun2v1DBoldDMwLT->fillHistograms(
//         recJetPreselection_passed && recJetOldDecayMode_passed && byTightIsolationMVArun2v1DBoldDMwLT, true,
//         recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byVTightIsolationMVArun2v1DBoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVTightIsolationMVArun2v1DBoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      //plotsJet_byVVTightIsolationMVArun2v1DBoldDMwLT->fillHistograms(
//      // recJetPreselection_passed && recJetOldDecayMode_passed && byVVTightIsolationMVArun2v1DBoldDMwLT, true,
//      // recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byVLooseIsolationMVArun2v1DBnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byVLooseIsolationMVArun2v1DBnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byLooseIsolationMVArun2v1DBnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byLooseIsolationMVArun2v1DBnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byMediumIsolationMVArun2v1DBnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byMediumIsolationMVArun2v1DBnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byTightIsolationMVArun2v1DBnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byTightIsolationMVArun2v1DBnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byVTightIsolationMVArun2v1DBnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byVTightIsolationMVArun2v1DBnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//   //plotsJet_byVVTightIsolationMVArun2v1DBnewDMwLT->fillHistograms(
//   //	 recJetPreselection_passed && recJetNewDecayMode_passed && byVVTightIsolationMVArun2v1DBnewDMwLT, true,
//   //	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byVLooseIsolationMVArun2v1PWoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVLooseIsolationMVArun2v1PWoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byLooseIsolationMVArun2v1PWoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byLooseIsolationMVArun2v1PWoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byMediumIsolationMVArun2v1PWoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byMediumIsolationMVArun2v1PWoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byTightIsolationMVArun2v1PWoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byTightIsolationMVArun2v1PWoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byVTightIsolationMVArun2v1PWoldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVTightIsolationMVArun2v1PWoldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
////plotsJet_byVVTightIsolationMVArun2v1PWoldDMwLT->fillHistograms(
////	 recJetPreselection_passed && recJetOldDecayMode_passed && byVVTightIsolationMVArun2v1PWoldDMwLT, true,
////	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byVLooseIsolationMVArun2v1PWnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byVLooseIsolationMVArun2v1PWnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byLooseIsolationMVArun2v1PWnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byLooseIsolationMVArun2v1PWnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byMediumIsolationMVArun2v1PWnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byMediumIsolationMVArun2v1PWnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byTightIsolationMVArun2v1PWnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byTightIsolationMVArun2v1PWnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byVTightIsolationMVArun2v1PWnewDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetNewDecayMode_passed && byVTightIsolationMVArun2v1PWnewDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
////plotsJet_byVVTightIsolationMVArun2v1PWnewDMwLT->fillHistograms(
////	 recJetPreselection_passed && recJetNewDecayMode_passed && byVVTightIsolationMVArun2v1PWnewDMwLT, true,
////	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
//         recJetPreselection_passed && recJetOldDecayMode_passed && byVLooseIsolationMVArun2v1DBdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byLooseIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byLooseIsolationMVArun2v1DBdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byMediumIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byMediumIsolationMVArun2v1DBdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byTightIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byTightIsolationMVArun2v1DBdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byVTightIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVTightIsolationMVArun2v1DBdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
////      plotsJet_byVVTightIsolationMVArun2v1DBdR03oldDMwLT->fillHistograms(
////	 recJetPreselection_passed && recJetOldDecayMode_passed && byVVTightIsolationMVArun2v1DBdR03oldDMwLT, true,
////	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//
//      plotsJet_byVLooseIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVLooseIsolationMVArun2v1PWdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byLooseIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byLooseIsolationMVArun2v1PWdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byMediumIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byMediumIsolationMVArun2v1PWdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byTightIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byTightIsolationMVArun2v1PWdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      plotsJet_byVTightIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
//	 recJetPreselection_passed && recJetOldDecayMode_passed && byVTightIsolationMVArun2v1PWdR03oldDMwLT, true,
//	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
////      plotsJet_byVVTightIsolationMVArun2v1PWdR03oldDMwLT->fillHistograms(
////	 recJetPreselection_passed && recJetOldDecayMode_passed && byVVTightIsolationMVArun2v1PWdR03oldDMwLT, true,
////	 recJetPt, recJetEta, recJetEtaImpactECAL, recJetPhi, leadPFChargedHadrCandPt, leadPFChargedHadrCandEta, tauNeutralHadPt, recTauDecayMode, 0, numVertices, evtWeight_dummy);
//      
//      ++numJets;
//    }


  }

  delete inputFile;

  std::cout << "numElectrons = " << numElectrons << std::endl;
  std::cout << "numMuons = " << numMuons << std::endl;
  std::cout << "numTaus = " << numTaus << std::endl;
  std::cout << "numJets = " << numJets << std::endl;
  
  assert((numElectrons + numMuons + numTaus + numJets) > 0);

  plotsTau_againstElectronVLooseMVA6->makeGraphs();
  plotsTau_againstElectronLooseMVA6->makeGraphs();
  plotsTau_againstElectronMediumMVA6->makeGraphs();
  plotsTau_againstElectronTightMVA6->makeGraphs(); 
  plotsTau_againstElectronVTightMVA6->makeGraphs(); 

  plotsTau_againstMuonLoose3->makeGraphs();
  plotsTau_againstMuonTight3->makeGraphs();

  plotsElectron_againstElectronVLooseMVA6->makeGraphs();
  plotsElectron_againstElectronLooseMVA6->makeGraphs();
  plotsElectron_againstElectronMediumMVA6->makeGraphs();
  plotsElectron_againstElectronTightMVA6->makeGraphs();
  plotsElectron_againstElectronVTightMVA6->makeGraphs();

  plotsMuon_againstMuonLoose3->makeGraphs();
  plotsMuon_againstMuonTight3->makeGraphs();

  plotsTau_recoTauJetmatching->makeGraphs();
  plotsTau_byDecayModeFindingNewDMs->makeGraphs();   
  plotsTau_byDecayModeFindingOldDMs->makeGraphs();   
  plotsTau_byCombinedIsolation3HitsLoose->makeGraphs();   
  plotsTau_byCombinedIsolation3HitsMedium->makeGraphs(); 
  plotsTau_byCombinedIsolation3HitsTight->makeGraphs(); 
  plotsTau_byCombinedIsolationNewDM3HitsLoose->makeGraphs();
  plotsTau_byCombinedIsolationNewDM3HitsMedium->makeGraphs();
  plotsTau_byCombinedIsolationNewDM3HitsTight->makeGraphs();
  plotsTau_byCombinedIsolation3HitsdR03Loose->makeGraphs();
  plotsTau_byCombinedIsolation3HitsdR03Medium->makeGraphs();
  plotsTau_byCombinedIsolation3HitsdR03Tight->makeGraphs();
  plotsTau_byChargedIsolationLoose->makeGraphs(); 
  plotsTau_byChargedIsolationMedium->makeGraphs(); 
  plotsTau_byChargedIsolationTight->makeGraphs(); 
  plotsTau_byPileupWeightedIsolation3HitsLoose->makeGraphs();
  plotsTau_byPileupWeightedIsolation3HitsMedium->makeGraphs();
  plotsTau_byPileupWeightedIsolation3HitsTight->makeGraphs();
  plotsTau_byPileupWeightedIsolationNewDM3HitsLoose->makeGraphs();
  plotsTau_byPileupWeightedIsolationNewDM3HitsMedium->makeGraphs();
  plotsTau_byPileupWeightedIsolationNewDM3HitsTight->makeGraphs();
  /*
  plotsTau_genDM0_recoTauJetmatching->makeGraphs();
  plotsTau_genDM0_byDecayModeFindingNewDMs->makeGraphs();
  plotsTau_genDM0_byDecayModeFindingOldDMs->makeGraphs();
  plotsTau_genDM0_byCombinedIsolation3HitsLoose->makeGraphs();
  plotsTau_genDM0_byCombinedIsolation3HitsMedium->makeGraphs();
  plotsTau_genDM0_byCombinedIsolation3HitsTight->makeGraphs();
  plotsTau_genDM0_byChargedIsolationLoose->makeGraphs();
  plotsTau_genDM0_byChargedIsolationMedium->makeGraphs();
  plotsTau_genDM0_byChargedIsolationTight->makeGraphs();

  plotsTau_genDM1_recoTauJetmatching->makeGraphs();
  plotsTau_genDM1_byDecayModeFindingNewDMs->makeGraphs();
  plotsTau_genDM1_byDecayModeFindingOldDMs->makeGraphs();
  plotsTau_genDM1_byChargedIsolationLoose->makeGraphs();
  plotsTau_genDM1_byChargedIsolationMedium->makeGraphs();
  plotsTau_genDM1_byChargedIsolationTight->makeGraphs();

  plotsTau_genDM2_recoTauJetmatching->makeGraphs();
  plotsTau_genDM2_byDecayModeFindingNewDMs->makeGraphs();
  plotsTau_genDM2_byDecayModeFindingOldDMs->makeGraphs();
  plotsTau_genDM2_byCombinedIsolation3HitsLoose->makeGraphs();
  plotsTau_genDM2_byCombinedIsolation3HitsMedium->makeGraphs();
  plotsTau_genDM2_byCombinedIsolation3HitsTight->makeGraphs();
  plotsTau_genDM2_byChargedIsolationLoose->makeGraphs();
  plotsTau_genDM2_byChargedIsolationMedium->makeGraphs();
  plotsTau_genDM2_byChargedIsolationTight->makeGraphs();

  plotsTau_genDM10_recoTauJetmatching->makeGraphs();
  plotsTau_genDM10_byDecayModeFindingNewDMs->makeGraphs();
  plotsTau_genDM10_byDecayModeFindingOldDMs->makeGraphs();
  plotsTau_genDM10_byCombinedIsolation3HitsLoose->makeGraphs();
  plotsTau_genDM10_byCombinedIsolation3HitsMedium->makeGraphs();
  plotsTau_genDM10_byCombinedIsolation3HitsTight->makeGraphs();
  plotsTau_genDM10_byChargedIsolationLoose->makeGraphs();
  plotsTau_genDM10_byChargedIsolationMedium->makeGraphs();
  plotsTau_genDM10_byChargedIsolationTight->makeGraphs();
  */

  plotsTau_byVLooseIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsTau_byLooseIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsTau_byMediumIsolationMVArun2v1DBoldDMwLT ->makeGraphs();
  plotsTau_byTightIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsTau_byVTightIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsTau_byVVTightIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsTau_byVLooseIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsTau_byLooseIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsTau_byMediumIsolationMVArun2v1DBnewDMwLT ->makeGraphs();
  plotsTau_byTightIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsTau_byVTightIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsTau_byVVTightIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  
  plotsTau_byVLooseIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsTau_byLooseIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsTau_byMediumIsolationMVArun2v1PWoldDMwLT ->makeGraphs();
  plotsTau_byTightIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsTau_byVTightIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsTau_byVVTightIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsTau_byVLooseIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsTau_byLooseIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsTau_byMediumIsolationMVArun2v1PWnewDMwLT ->makeGraphs();
  plotsTau_byTightIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsTau_byVTightIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsTau_byVVTightIsolationMVArun2v1PWnewDMwLT->makeGraphs();

  plotsTau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT ->makeGraphs();
  plotsTau_byTightIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsTau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsTau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsTau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsTau_byLooseIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsTau_byMediumIsolationMVArun2v1PWdR03oldDMwLT ->makeGraphs();
  plotsTau_byTightIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsTau_byVTightIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsTau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();

  plotsJet_byDecayModeFindingNewDMs->makeGraphs();   
  plotsJet_byDecayModeFindingOldDMs->makeGraphs();   
  plotsJet_byCombinedIsolation3HitsLoose->makeGraphs();   
  plotsJet_byCombinedIsolation3HitsMedium->makeGraphs(); 
  plotsJet_byCombinedIsolation3HitsTight->makeGraphs(); 
  plotsJet_byCombinedIsolationNewDM3HitsLoose->makeGraphs();
  plotsJet_byCombinedIsolationNewDM3HitsMedium->makeGraphs();
  plotsJet_byCombinedIsolationNewDM3HitsTight->makeGraphs();
  plotsJet_byCombinedIsolation3HitsdR03Loose->makeGraphs();
  plotsJet_byCombinedIsolation3HitsdR03Medium->makeGraphs();
  plotsJet_byCombinedIsolation3HitsdR03Tight->makeGraphs();
  plotsJet_byChargedIsolationLoose->makeGraphs(); 
  plotsJet_byChargedIsolationMedium->makeGraphs(); 
  plotsJet_byChargedIsolationTight->makeGraphs(); 
  plotsJet_byPileupWeightedIsolation3HitsLoose->makeGraphs();
  plotsJet_byPileupWeightedIsolation3HitsMedium->makeGraphs();
  plotsJet_byPileupWeightedIsolation3HitsTight->makeGraphs();
  plotsJet_byPileupWeightedIsolationNewDM3HitsLoose->makeGraphs();
  plotsJet_byPileupWeightedIsolationNewDM3HitsMedium->makeGraphs();
  plotsJet_byPileupWeightedIsolationNewDM3HitsTight->makeGraphs();

  plotsJet_byVLooseIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsJet_byLooseIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsJet_byMediumIsolationMVArun2v1DBoldDMwLT ->makeGraphs();
  plotsJet_byTightIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsJet_byVTightIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsJet_byVVTightIsolationMVArun2v1DBoldDMwLT->makeGraphs();
  plotsJet_byVLooseIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsJet_byLooseIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsJet_byMediumIsolationMVArun2v1DBnewDMwLT ->makeGraphs();
  plotsJet_byTightIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsJet_byVTightIsolationMVArun2v1DBnewDMwLT->makeGraphs();
  plotsJet_byVVTightIsolationMVArun2v1DBnewDMwLT->makeGraphs();

  plotsJet_byVLooseIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsJet_byLooseIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsJet_byMediumIsolationMVArun2v1PWoldDMwLT ->makeGraphs();
  plotsJet_byTightIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsJet_byVTightIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsJet_byVVTightIsolationMVArun2v1PWoldDMwLT->makeGraphs();
  plotsJet_byVLooseIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsJet_byLooseIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsJet_byMediumIsolationMVArun2v1PWnewDMwLT ->makeGraphs();
  plotsJet_byTightIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsJet_byVTightIsolationMVArun2v1PWnewDMwLT->makeGraphs();
  plotsJet_byVVTightIsolationMVArun2v1PWnewDMwLT->makeGraphs();

  plotsJet_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsJet_byLooseIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsJet_byMediumIsolationMVArun2v1DBdR03oldDMwLT ->makeGraphs();
  plotsJet_byTightIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsJet_byVTightIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsJet_byVVTightIsolationMVArun2v1DBdR03oldDMwLT->makeGraphs();
  plotsJet_byVLooseIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsJet_byLooseIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsJet_byMediumIsolationMVArun2v1PWdR03oldDMwLT ->makeGraphs();
  plotsJet_byTightIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsJet_byVTightIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();
  plotsJet_byVVTightIsolationMVArun2v1PWdR03oldDMwLT->makeGraphs();

  return plotEntryMap;
}

void makeTauIdEffPlots2_newMVAId()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  //TString inputFileName_newTauId = "/nfs/dust/cms/user/anayak/CMS/Ntuple_Phys14TauId/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/";
  //TString inputFileName_newTauId = "/nfs/dust/cms/user/anayak/CMS/Ntuple_Spring15TauID/Validation/Fall15/";
  //inputFileName_newTauId.Append("DYJetsToLL_M-50_13TeV_Fall15_25ns_miniAODv2_v1.root");
  //inputFileName_newTauId.Append("ZprimeToTauTau_M1000_13TeV_Fall15_25ns_miniAODv2_v1.root");
  TString inputFileName_newTauId = "/nfs/dust/cms/user/anayak/CMS/OnSLC6/CMSSW_763_tauId/src/TauAnalysis/Test/test/";
  inputFileName_newTauId.Append("pfTauIdEffNtupleFromMiniAOD_Fall15_DY_v1p1_run2.root");
  //inputFileName_newTauId.Append("pfTauIdEffNtupleFromMiniAOD_Fall15_Zp2_run2.root");

  std::map<std::string, std::map<std::string, plotEntryType*> > plots_newTauId = makePlots(inputFileName_newTauId, "newTauId", kNewTags);

  //TString inputFileName_Phys14 = "/nfs/dust/cms/user/anayak/CMS/Ntuple_Phys14TauId/5thDec2014/";
  //inputFileName_Phys14.Append("pfTauIdEffNtuple2_Phys14_All_QCD_V5.root");

  //std::map<std::string, std::map<std::string, plotEntryType*> > plots_Phys14 = makePlots(inputFileName_Phys14, "Phys14", kNewTags);
  
  //std::string legenEntry_newTauId = "CMSSW 5_3_x";
  //std::string legenEntry_Phys14 = "Phys14";
  std::string legenEntry_vloose = "VLoose";
  std::string legenEntry_loose = "Loose";
  std::string legenEntry_medium = "Medium";
  std::string legenEntry_tight = "Tight";
  std::string legenEntry_vtight = "VTight";
  std::string legenEntry_vvtight = "VVTight";
  /*
  showGraphs(800, 600,
             plots_newTauId["Tau"]["recoTauJetmatching"], "Reco Tau Match",
             plots_newTauId["Tau"]["byDecayModeFindingNewDMs"], "New DM",
             plots_newTauId["Tau"]["byDecayModeFindingOldDMs"], "Old DM",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_decayMode.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byChargedIsolationLoose"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byChargedIsolationMedium"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byChargedIsolationTight"], legenEntry_tight.data(),
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_chargeIso.pdf");
  */
  showGraphs(800, 600,
             plots_newTauId["Tau"]["byCombinedIsolation3HitsLoose"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byCombinedIsolation3HitsMedium"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byCombinedIsolation3HitsTight"], legenEntry_tight.data(),
             0, "",
             0, "",
	     0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0.2, 0.9, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_combIso3Hits.pdf");
  /*
  showGraphs(800, 600,
             plots_newTauId["Tau"]["byCombinedIsolationNewDM3HitsLoose"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byCombinedIsolationNewDM3HitsMedium"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byCombinedIsolationNewDM3HitsTight"], legenEntry_tight.data(),
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_combIso3HitsNewDM.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byCombinedIsolation3HitsdR03Loose"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byCombinedIsolation3HitsdR03Medium"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byCombinedIsolation3HitsdR03Tight"], legenEntry_tight.data(),
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_combIso3HitsdR03.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byPileupWeightedIsolation3HitsLoose"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byPileupWeightedIsolation3HitsMedium"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byPileupWeightedIsolation3HitsTight"], legenEntry_tight.data(),
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_pwIso3Hits.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byPileupWeightedIsolationNewDM3HitsLoose"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byPileupWeightedIsolationNewDM3HitsMedium"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byPileupWeightedIsolationNewDM3HitsTight"], legenEntry_tight.data(),
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_pwIso3HitsNewDM.pdf");
  */
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byVLooseIsolationMVArun2v1DBoldDMwLT"], legenEntry_vloose.data(),
             plots_newTauId["Tau"]["byLooseIsolationMVArun2v1DBoldDMwLT"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byMediumIsolationMVArun2v1DBoldDMwLT"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byTightIsolationMVArun2v1DBoldDMwLT"], legenEntry_tight.data(),
             plots_newTauId["Tau"]["byVTightIsolationMVArun2v1DBoldDMwLT"], legenEntry_vtight.data(),
	     plots_newTauId["Tau"]["byVVTightIsolationMVArun2v1DBoldDMwLT"], legenEntry_vvtight.data(),
             0.045, 0.18, 0.63, 0.28, 0.25,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0.1, 1.1, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_run2MVA1isoDBoldDMwLT.pdf");
  /*
  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVArun2v1DBnewDMwLT"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byMediumIsolationMVArun2v1DBnewDMwLT"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byTightIsolationMVArun2v1DBnewDMwLT"], legenEntry_tight.data(),
             plots_newTauId["Tau"]["byVTightIsolationMVArun2v1DBnewDMwLT"], legenEntry_vtight.data(),
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_run2MVA1isoDBnewDMwLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVArun2v1PWoldDMwLT"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byMediumIsolationMVArun2v1PWoldDMwLT"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byTightIsolationMVArun2v1PWoldDMwLT"], legenEntry_tight.data(),
             plots_newTauId["Tau"]["byVTightIsolationMVArun2v1PWoldDMwLT"], legenEntry_vtight.data(),
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_run2MVA1isoPWoldDMwLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVArun2v1PWnewDMwLT"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byMediumIsolationMVArun2v1PWnewDMwLT"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byTightIsolationMVArun2v1PWnewDMwLT"], legenEntry_tight.data(),
             plots_newTauId["Tau"]["byVTightIsolationMVArun2v1PWnewDMwLT"], legenEntry_vtight.data(),
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_run2MVA1isoPWnewDMwLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byMediumIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byTightIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_tight.data(),
             plots_newTauId["Tau"]["byVTightIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_vtight.data(),
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_run2MVA1isoDBdR03oldDMwLT.pdf");
  
  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_loose.data(),
             plots_newTauId["Tau"]["byMediumIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_medium.data(),
             plots_newTauId["Tau"]["byTightIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_tight.data(),
             plots_newTauId["Tau"]["byVTightIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_vtight.data(),
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.3, "Efficiency", 1.2,
             "makeTauIdEffPlots_newMVAId_Taus_run2MVA1isoPWdR03oldDMwLT.pdf");
  */
  /*
  showGraphs(800, 600,
             plots_newTauId["Tau"]["recoTauJetmatching"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["recoTauJetmatching"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_recoTauJetmatching.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["recoTauJetmatching"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["recoTauJetmatching"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_recoTauJetmatching.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["recoTauJetmatching"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["recoTauJetmatching"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_recoTauJetmatching.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["recoTauJetmatching"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["recoTauJetmatching"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_recoTauJetmatching.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["recoTauJetmatching"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["recoTauJetmatching"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_recoTauJetmatching.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byDecayModeFindingNewDMs"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["byDecayModeFindingNewDMs"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_decayModeFindingNewDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byDecayModeFindingNewDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byDecayModeFindingNewDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_decayModeFindingNewDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byDecayModeFindingNewDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byDecayModeFindingNewDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_decayModeFindingNewDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byDecayModeFindingNewDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byDecayModeFindingNewDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_decayModeFindingNewDMs.pdf");
  
  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byDecayModeFindingNewDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byDecayModeFindingNewDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_decayModeFindingNewDMs.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byDecayModeFindingOldDMs"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["byDecayModeFindingOldDMs"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_decayModeFindingOldDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byDecayModeFindingOldDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byDecayModeFindingOldDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_decayModeFindingOldDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byDecayModeFindingOldDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byDecayModeFindingOldDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_decayModeFindingOldDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byDecayModeFindingOldDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byDecayModeFindingOldDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_decayModeFindingOldDMs.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byDecayModeFindingOldDMs"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byDecayModeFindingOldDMs"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_decayModeFindingOldDMs.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byCombinedIsolation8HitsLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["byCombinedIsolation8HitsLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_combIso8HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byCombinedIsolation8HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byCombinedIsolation8HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_combIso8HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byCombinedIsolation8HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byCombinedIsolation8HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_combIso8HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byCombinedIsolation8HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byCombinedIsolation8HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_combIso8HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byCombinedIsolation8HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byCombinedIsolation8HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_combIso8HitsLoose.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau_genDM0"]["byCombinedIsolation8HitsMedium"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau_genDM0"]["byCombinedIsolation8HitsMedium"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_combIso8HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byCombinedIsolation8HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byCombinedIsolation8HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_combIso8HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byCombinedIsolation8HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byCombinedIsolation8HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_combIso8HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byCombinedIsolation8HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byCombinedIsolation8HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_combIso8HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byCombinedIsolation8HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byCombinedIsolation8HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_combIso8HitsMedium.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byCombinedIsolation8HitsTight"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["byCombinedIsolation8HitsTight"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_combIso8HitsTight.pdf");
  
  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byCombinedIsolation8HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byCombinedIsolation8HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_combIso8HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byCombinedIsolation8HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byCombinedIsolation8HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_combIso8HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byCombinedIsolation8HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byCombinedIsolation8HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_combIso8HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byCombinedIsolation8HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byCombinedIsolation8HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_combIso8HitsTight.pdf");
  
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byCombinedIsolation3HitsLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["byCombinedIsolation3HitsLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_combIso3HitsLoose.pdf");
	     
  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byCombinedIsolation3HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byCombinedIsolation3HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_combIso3HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byCombinedIsolation3HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byCombinedIsolation3HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_combIso3HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byCombinedIsolation3HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byCombinedIsolation3HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_combIso3HitsLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byCombinedIsolation3HitsLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byCombinedIsolation3HitsLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_combIso3HitsLoose.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byCombinedIsolation3HitsMedium"], legenEntry_newTauId.data(),  
	     plots_Phys14["Tau"]["byCombinedIsolation3HitsMedium"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_combIso3HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byCombinedIsolation3HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byCombinedIsolation3HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_combIso3HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byCombinedIsolation3HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byCombinedIsolation3HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_combIso3HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byCombinedIsolation3HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byCombinedIsolation3HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_combIso3HitsMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byCombinedIsolation3HitsMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byCombinedIsolation3HitsMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_combIso3HitsMedium.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byCombinedIsolation3HitsTight"], legenEntry_newTauId.data(),  
	     plots_Phys14["Tau"]["byCombinedIsolation3HitsTight"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_combIso3HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byCombinedIsolation3HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byCombinedIsolation3HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_combIso3HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byCombinedIsolation3HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byCombinedIsolation3HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_combIso3HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byCombinedIsolation3HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byCombinedIsolation3HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_combIso3HitsTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byCombinedIsolation3HitsTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byCombinedIsolation3HitsTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_combIso3HitsTight.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byChargedIsolationLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["byChargedIsolationLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_chargedIsoLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byChargedIsolationLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byChargedIsolationLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_chargedIsoLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byChargedIsolationLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byChargedIsolationLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_chargedIsoLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byChargedIsolationLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byChargedIsolationLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_chargedIsoLoose.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byChargedIsolationLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byChargedIsolationLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_chargedIsoLoose.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byChargedIsolationMedium"], legenEntry_newTauId.data(),  
	     plots_Phys14["Tau"]["byChargedIsolationMedium"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_chargedIsoMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byChargedIsolationMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byChargedIsolationMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_chargedIsoMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byChargedIsolationMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byChargedIsolationMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_chargedIsoMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byChargedIsolationMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byChargedIsolationMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_chargedIsoMedium.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byChargedIsolationMedium"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau_genDM10"]["byChargedIsolationMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_chargedIsoMedium.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["byChargedIsolationTight"], legenEntry_newTauId.data(),  
	     plots_Phys14["Tau"]["byChargedIsolationTight"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_chargedIsoTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM0"]["byChargedIsolationTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM0"]["byChargedIsolationTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM0_chargedIsoTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM1"]["byChargedIsolationTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM1"]["byChargedIsolationTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM1_chargedIsoTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM2"]["byChargedIsolationTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM2"]["byChargedIsolationTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM2_chargedIsoTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau_genDM10"]["byChargedIsolationTight"], legenEntry_newTauId.data(),
             plots_Phys14["Tau_genDM10"]["byChargedIsolationTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_genDM10_chargedIsoTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVA3oldDMwoLT"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byLooseIsolationMVA3oldDMwoLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_mva3IsoLooseoldDMwoLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byMediumIsolationMVA3oldDMwoLT"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byMediumIsolationMVA3oldDMwoLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_mva3IsoMediumoldDMwoLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byTightIsolationMVA3oldDMwoLT"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byTightIsolationMVA3oldDMwoLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_mva3IsoTightoldDMwoLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byLooseIsolationMVA3oldDMwLT"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byLooseIsolationMVA3oldDMwLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_mva3IsoLooseoldDMwLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byMediumIsolationMVA3oldDMwLT"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byMediumIsolationMVA3oldDMwLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_mva3IsoMediumoldDMwLT.pdf");

  showGraphs(800, 600,
             plots_newTauId["Tau"]["byTightIsolationMVA3oldDMwLT"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["byTightIsolationMVA3oldDMwLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_mva3IsoTightoldDMwLT.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstElectronLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstElectronLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronLoose.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstElectronLooseMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstElectronLooseMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronLooseMVA6.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstElectronMediumMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstElectronMediumMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronMediumMVA6.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstElectronTightMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstElectronTightMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronTightMVA6.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstElectronVTightMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstElectronVTightMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronVTightMVA6.pdf");


  showGraphs(800, 600,                                                                                                                                                             
	     plots_newTauId["Tau"]["againstElectronLoose"], legenEntry_newTauId.data(),                                                                                        
	     plots_Phys14["Tau"]["againstElectronLoose"], legenEntry_Phys14.data(),                                                                                           
	     0, "",                                                                                                                                                                
	     0, "",                                                                                                                                                                
	     0, "",                                                                                                                                                                
	     0.045, 0.18, 0.67, 0.28, 0.21,                                                                                                                                        
	     "Z #rightarrow #tau#tau", 0.055,                                                                                                                                      
	     0.74, 0.16, 0.15, 0.06, 
	     false, 0., 1.5, "Efficiency", 1.2,                                                                                                                                    
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronLoose.pdf");                                                                                                 

  showGraphs(800, 600,                                                                                                                                                             
	     plots_newTauId["Tau"]["againstElectronMedium"], legenEntry_newTauId.data(),                                                                                       
	     plots_Phys14["Tau"]["againstElectronMedium"], legenEntry_Phys14.data(),                                                                                          
	     0, "",                                                                                                                                                              
	     0, "",                                                                                                                                                             
	     0, "",                                                                                                                                                            
	     0.045, 0.18, 0.67, 0.28, 0.21,                                                                                                                                   
	     "Z #rightarrow #tau#tau", 0.055,                                                                                                                                
	     0.74, 0.16, 0.15, 0.06,                                                                                                                                        
	     false, 0., 1.5, "Efficiency", 1.2,                                                                                                                            
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronMedium.pdf");                                                                                       

  showGraphs(800, 600,                                                                                                                                                             
	     plots_newTauId["Tau"]["againstElectronTight"], legenEntry_newTauId.data(),                                                                                      
	     plots_Phys14["Tau"]["againstElectronTight"], legenEntry_Phys14.data(),                                                                                         
	     0, "",                                                                                                                                                            
	     0, "",                                                                                                                                                           
	     0, "",                                                                                                                                                          
	     0.045, 0.18, 0.67, 0.28, 0.21,                                                                                                                                 
	     "Z #rightarrow #tau#tau", 0.055,                                                                                                                              
	     0.74, 0.16, 0.15, 0.06,                                                                                                                                      
	     false, 0., 1.5, "Efficiency", 1.2,                                                                                                                          
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstElectronTight.pdf");
  
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstMuonLoose3"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstMuonLoose3"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstMuonLoose3.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Tau"]["againstMuonTight3"], legenEntry_newTauId.data(),
	     plots_Phys14["Tau"]["againstMuonTight3"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #tau#tau", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     false, 0., 1.5, "Efficiency", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Taus_againstMuonTight3.pdf");
  showGraphs(800, 600,
             plots_newTauId["Tau"]["againstMuonLooseMVA"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["againstMuonLooseMVA"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_againstMuonLooseMVA.pdf");
  showGraphs(800, 600,
             plots_newTauId["Tau"]["againstMuonTightMVA"], legenEntry_newTauId.data(),
             plots_Phys14["Tau"]["againstMuonTightMVA"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #tau#tau", 0.055,
             0.74, 0.16, 0.15, 0.06,
             false, 0., 1.5, "Efficiency", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Taus_againstMuonTightMVA.pdf");
  */

//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["recoTauJetmatching"], "Reco Tau Match",
//             plots_newTauId["Jet"]["byDecayModeFindingNewDMs"], "New DM",
//             plots_newTauId["Jet"]["byDecayModeFindingOldDMs"], "Old DM",
//             0, "",
//             0, "",
//	     0.045, 0.58, 0.67, 0.28, 0.21,
//	     "QCD", 0.055,  
//	     0.74, 0.16, 0.15, 0.06,
//	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_decayMode.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byChargedIsolationLoose"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byChargedIsolationMedium"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byChargedIsolationTight"], legenEntry_tight.data(),
//             0, "",
//             0, "",
//	     0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_chargeIso.pdf");
//
//  showGraphs(800, 600,                                                     
//	     plots_newTauId["Jet"]["byCombinedIsolation3HitsLoose"], legenEntry_loose.data(),
//	     plots_newTauId["Jet"]["byCombinedIsolation3HitsMedium"], legenEntry_medium.data(),
//	     plots_newTauId["Jet"]["byCombinedIsolation3HitsTight"], legenEntry_tight.data(),
//	     0, "",
//	     0, "",
//	     0.045, 0.58, 0.67, 0.28, 0.21,
//	     "QCD", 0.055,
//	     0.74, 0.16, 0.15, 0.06,
//	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//	     "makeTauIdEffPlots_newMVAId_Jets_combIso3Hits.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byCombinedIsolationNewDM3HitsLoose"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byCombinedIsolationNewDM3HitsMedium"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byCombinedIsolationNewDM3HitsTight"], legenEntry_tight.data(),
//             0, "",
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_combIso3HitsNewDM.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byCombinedIsolation3HitsdR03Loose"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byCombinedIsolation3HitsdR03Medium"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byCombinedIsolation3HitsdR03Tight"], legenEntry_tight.data(),
//             0, "",
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_combIso3HitsdR03.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byPileupWeightedIsolation3HitsLoose"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byPileupWeightedIsolation3HitsMedium"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byPileupWeightedIsolation3HitsTight"], legenEntry_tight.data(),
//             0, "",
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_pwIso3Hits.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byPileupWeightedIsolationNewDM3HitsLoose"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byPileupWeightedIsolationNewDM3HitsMedium"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byPileupWeightedIsolationNewDM3HitsTight"], legenEntry_tight.data(),
//             0, "",
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_pwIso3HitsNewDM.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byLooseIsolationMVArun2v1DBoldDMwLT"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byMediumIsolationMVArun2v1DBoldDMwLT"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byTightIsolationMVArun2v1DBoldDMwLT"], legenEntry_tight.data(),
//             plots_newTauId["Jet"]["byVTightIsolationMVArun2v1DBoldDMwLT"], legenEntry_vtight.data(),
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_run2MVA1isoDBoldDMwLT.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byLooseIsolationMVArun2v1DBnewDMwLT"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byMediumIsolationMVArun2v1DBnewDMwLT"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byTightIsolationMVArun2v1DBnewDMwLT"], legenEntry_tight.data(),
//             plots_newTauId["Jet"]["byVTightIsolationMVArun2v1DBnewDMwLT"], legenEntry_vtight.data(),
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_run2MVA1isoDBnewDMwLT.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byLooseIsolationMVArun2v1PWoldDMwLT"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byMediumIsolationMVArun2v1PWoldDMwLT"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byTightIsolationMVArun2v1PWoldDMwLT"], legenEntry_tight.data(),
//             plots_newTauId["Jet"]["byVTightIsolationMVArun2v1PWoldDMwLT"], legenEntry_vtight.data(),
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_run2MVA1isoPWoldDMwLT.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byLooseIsolationMVArun2v1PWnewDMwLT"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byMediumIsolationMVArun2v1PWnewDMwLT"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byTightIsolationMVArun2v1PWnewDMwLT"], legenEntry_tight.data(),
//             plots_newTauId["Jet"]["byVTightIsolationMVArun2v1PWnewDMwLT"], legenEntry_vtight.data(),
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_run2MVA1isoPWnewDMwLT.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byLooseIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byMediumIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byTightIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_tight.data(),
//             plots_newTauId["Jet"]["byVTightIsolationMVArun2v1DBdR03oldDMwLT"], legenEntry_vtight.data(),
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_run2MVA1isoDBdR03oldDMwLT.pdf");
//
//  showGraphs(800, 600,
//             plots_newTauId["Jet"]["byLooseIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_loose.data(),
//             plots_newTauId["Jet"]["byMediumIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_medium.data(),
//             plots_newTauId["Jet"]["byTightIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_tight.data(),
//             plots_newTauId["Jet"]["byVTightIsolationMVArun2v1PWdR03oldDMwLT"], legenEntry_vtight.data(),
//             0, "",
//             0.045, 0.58, 0.67, 0.28, 0.21,
//             "QCD", 0.055,
//             0.74, 0.16, 0.15, 0.06,
//             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
//             "makeTauIdEffPlots_newMVAId_Jets_run2MVA1isoPWdR03oldDMwLT.pdf");

  /*
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byDecayModeFindingNewDMs"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byDecayModeFindingNewDMs"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_decayModeFindingNewDMs.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byDecayModeFindingOldDMs"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byDecayModeFindingOldDMs"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_decayModeFindingOldDMs.pdf");
  
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byCombinedIsolation8HitsLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byCombinedIsolation8HitsLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_combIso8HitsLoose.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byCombinedIsolation8HitsMedium"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byCombinedIsolation8HitsMedium"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_combIso8HitsMedium.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byCombinedIsolation8HitsTight"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byCombinedIsolation8HitsTight"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_combIso8HitsTight.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byCombinedIsolation3HitsLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byCombinedIsolation3HitsLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_combIso3HitsLoose.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byCombinedIsolation3HitsMedium"], legenEntry_newTauId.data(),  
	     plots_Phys14["Jet"]["byCombinedIsolation3HitsMedium"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_combIso3HitsMedium.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byCombinedIsolation3HitsTight"], legenEntry_newTauId.data(),  
	     plots_Phys14["Jet"]["byCombinedIsolation3HitsTight"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_combIso3HitsTight.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byChargedIsolationLoose"], legenEntry_newTauId.data(),
	     plots_Phys14["Jet"]["byChargedIsolationLoose"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_chargedIsoLoose.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byChargedIsolationMedium"], legenEntry_newTauId.data(),  
	     plots_Phys14["Jet"]["byChargedIsolationMedium"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_chargedIsoMedium.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Jet"]["byChargedIsolationTight"], legenEntry_newTauId.data(),  
	     plots_Phys14["Jet"]["byChargedIsolationTight"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.58, 0.67, 0.28, 0.21,
	     "QCD", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Jets_chargedIsoTight.pdf");

  showGraphs(800, 600,
             plots_newTauId["Jet"]["byChargedIsolationUser1"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byChargedIsolationUser1"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_chargedIsoUser1.pdf");

  showGraphs(800, 600,
             plots_newTauId["Jet"]["byChargedIsolationUser2"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byChargedIsolationUser2"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_chargedIsoUser2.pdf");

  showGraphs(800, 600,
             plots_newTauId["Jet"]["byChargedIsolationUser3"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byChargedIsolationUser3"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_chargedIsoUser3.pdf");

  showGraphs(800, 600,
             plots_newTauId["Jet"]["byLooseIsolationMVA3oldDMwoLT"], legenEntry_newTauId.data(),                                                                                 
             plots_Phys14["Jet"]["byLooseIsolationMVA3oldDMwoLT"], legenEntry_Phys14.data(),                                                                                     
             0, "",                                                                                                                                                              
             0, "",                                                                                                                                                              
             0, "",                                                                                                                                                              
             0.045, 0.58, 0.67, 0.28, 0.21,                                                                                                                                      
             "QCD", 0.055,                                                                                                                                    
	     0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,                                                                                                                                  
             "makeTauIdEffPlots_oldVsNewTauId_Jets_mva3IsoLooseoldDMwoLT.pdf");                                                                                                  
  showGraphs(800, 600,                                                                                                                                                           
             plots_newTauId["Jet"]["byMediumIsolationMVA3oldDMwoLT"], legenEntry_newTauId.data(),                                                                                
             plots_Phys14["Jet"]["byMediumIsolationMVA3oldDMwoLT"], legenEntry_Phys14.data(),                                                                                    
             0, "",                                                                                                                                                              
             0, "",                                                                                                                                                              
             0, "",                                                                                                                                                              
             0.045, 0.58, 0.67, 0.28, 0.21,                                                                                                                                      
	     "QCD", 0.055,                                                                                                                                    
             0.74, 0.16, 0.15, 0.06,                                                                                                                                             
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,                                                                                                                                 
             "makeTauIdEffPlots_oldVsNewTauId_Jets_mva3IsoMediumoldDMwoLT.pdf");                                                                                                 
  showGraphs(800, 600,
             plots_newTauId["Jet"]["byTightIsolationMVA3oldDMwoLT"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byTightIsolationMVA3oldDMwoLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_mva3IsoTightoldDMwoLT.pdf");
  showGraphs(800, 600,
             plots_newTauId["Jet"]["byLooseIsolationMVA3oldDMwLT"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byLooseIsolationMVA3oldDMwLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_mva3IsoLooseoldDMwLT.pdf");
  showGraphs(800, 600,
             plots_newTauId["Jet"]["byMediumIsolationMVA3oldDMwLT"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byMediumIsolationMVA3oldDMwLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_mva3IsoMediumoldDMwLT.pdf");
  showGraphs(800, 600,
             plots_newTauId["Jet"]["byTightIsolationMVA3oldDMwLT"], legenEntry_newTauId.data(),
             plots_Phys14["Jet"]["byTightIsolationMVA3oldDMwLT"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.58, 0.67, 0.28, 0.21,
             "QCD", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Jets_mva3IsoTightoldDMwLT.pdf");
  */
  /*
  showGraphs(800, 600,
	     plots_newTauId["Electron"]["againstElectronLooseMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Electron"]["againstElectronLooseMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow ee", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronLooseMVA6.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Electron"]["againstElectronMediumMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Electron"]["againstElectronMediumMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow ee", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronMediumMVA6.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Electron"]["againstElectronTightMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Electron"]["againstElectronTightMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow ee", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronTightMVA6.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Electron"]["againstElectronVTightMVA6"], legenEntry_newTauId.data(),
	     plots_Phys14["Electron"]["againstElectronVTightMVA6"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow ee", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronVTightMVA6.pdf");
  
  showGraphs(800, 600,
             plots_newTauId["Electron"]["againstElectronLoose"], legenEntry_newTauId.data(),
             plots_Phys14["Electron"]["againstElectronLoose"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow ee", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronLoose.pdf");
  showGraphs(800, 600,
             plots_newTauId["Electron"]["againstElectronMedium"], legenEntry_newTauId.data(),
             plots_Phys14["Electron"]["againstElectronMedium"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow ee", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronMedium.pdf");
  showGraphs(800, 600,
             plots_newTauId["Electron"]["againstElectronTight"], legenEntry_newTauId.data(),
             plots_Phys14["Electron"]["againstElectronTight"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow ee", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Electrons_againstElectronTight.pdf");

  showGraphs(800, 600,
	     plots_newTauId["Muon"]["againstMuonLoose3"], legenEntry_newTauId.data(),
	     plots_Phys14["Muon"]["againstMuonLoose3"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #mu#mu", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Muons_againstMuonLoose3.pdf");
  showGraphs(800, 600,
	     plots_newTauId["Muon"]["againstMuonTight3"], legenEntry_newTauId.data(),
	     plots_Phys14["Muon"]["againstMuonTight3"], legenEntry_Phys14.data(),
	     0, "",
	     0, "",
	     0, "",
	     0.045, 0.18, 0.67, 0.28, 0.21,
	     "Z #rightarrow #mu#mu", 0.055, 
	     0.74, 0.16, 0.15, 0.06,
	     true, 1.e-4, 1.e0, "Fake-rate", 1.2,
	     "makeTauIdEffPlots_oldVsNewTauId_Muons_againstMuonTight3.pdf");

  showGraphs(800, 600,
             plots_newTauId["Muon"]["againstMuonLooseMVA"], legenEntry_newTauId.data(),
             plots_Phys14["Muon"]["againstMuonLooseMVA"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #mu#mu", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Muons_againstMuonLooseMVA.pdf");
  showGraphs(800, 600,
             plots_newTauId["Muon"]["againstMuonTightMVA"], legenEntry_newTauId.data(),
             plots_Phys14["Muon"]["againstMuonTightMVA"], legenEntry_Phys14.data(),
             0, "",
             0, "",
             0, "",
             0.045, 0.18, 0.67, 0.28, 0.21,
             "Z #rightarrow #mu#mu", 0.055,
             0.74, 0.16, 0.15, 0.06,
             true, 1.e-4, 1.e0, "Fake-rate", 1.2,
             "makeTauIdEffPlots_oldVsNewTauId_Muons_againstMuonTightMVA.pdf");
  */
}
