#ifndef TauAnalysis_Test_DumpTriggerInfo_h
#define TauAnalysis_Test_DumpTriggerInfo_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/MonitorElement.h"

class DumpTriggerInfo : public edm::EDAnalyzer 
{
 public:

  explicit DumpTriggerInfo(const edm::ParameterSet&);
  ~DumpTriggerInfo() {}
    
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:

  edm::InputTag srcL1GtReadoutRecord_; 
  //edm::InputTag srcL1GtObjectMaps_;
  edm::InputTag srcHLTresults_; 
};

#endif
