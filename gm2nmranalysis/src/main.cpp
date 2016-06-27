#include "TNMRPulse.h"
#include "TrolleyAnalysisManager.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

int main(){
  TrolleyAnalysisManager *manager = new TrolleyAnalysisManager;

  manager->GetPulse()->SetClockFrequencyInMHz(61.7);
  manager->SetInputFile(string{"data/data-nmr-Run1.dat"});

  manager->ReadNextPulse();

  shared_ptr<TGraphErrors> g = manager->GetPulse()->GetRawGraph();
  TCanvas *c = new TCanvas();

  g->Draw("AP");
  c->SaveAs("c1.png");

  return 1;
}
