#include "TNMRPulse.h"
#include "TrolleyAnalysisManager.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

int main(){
  // Need to add reading in arguments to set the analysis parameters
  // This would be great within Art and FHICL

  TrolleyAnalysisManager *manager = new TrolleyAnalysisManager;

  manager->GetPulse()->SetClockFrequencyInMHz(61.7);
  manager->SetInputFile(string{"data/data-nmr-Run1.dat"});

  int i = 0;

  while(manager->ReadNextPulse() > 0){
    cout << ++i << endl;
    shared_ptr<TGraphErrors> g = manager->GetPulse()->GetRawGraph();
    
    // Start analysis, will add this to the TrolleyAnalysisManager
  }


  return 1;
}
