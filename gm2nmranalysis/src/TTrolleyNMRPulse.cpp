#include "TTrolleyNMRPulse.h"

TTrolleyNMRPulse::TTrolleyNMRPulse(const TString& Name, const TString& Title) : TNMRPulse(Name, Title)
{
}

TTrolleyNMRPulse::TTrolleyNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy) : TNMRPulse(Name, Title, fx, fy)
{
}

TTrolleyNMRPulse::~TTrolleyNMRPulse(){
  
}

void TTrolleyNMRPulse::SetFrequency(float f){
  frequency = f;
}

void TTrolleyNMRPulse::SetInputFile(string input){
  inputFileName.assign(input);
}
