#include "TNMRPulse.h"
#include "TTrolleyNMRPulse.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

int main(){
  TTrolleyNMRPulse *p = new TTrolleyNMRPulse;

  p->SetClockFrequencyInMHz(61.7);
  p->SetInputFile(string{"data/data-nmr-Run1.dat"});

  p->ReadNextPulse();

  shared_ptr<TGraphErrors> g = p->GetRawGraph();
  TCanvas *c = new TCanvas();

  g->Draw("AP");
  c->SaveAs("c1.png");

  return 1;
}
