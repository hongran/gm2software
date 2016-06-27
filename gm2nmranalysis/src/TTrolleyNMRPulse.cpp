#include "TTrolleyNMRPulse.h"

TTrolleyNMRPulse::TTrolleyNMRPulse(const TString& Name, const TString& Title) : TNMRPulse(Name, Title)
{
  Initialize();
}


TTrolleyNMRPulse::TTrolleyNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy, vector<double> fyerr) 
  : TNMRPulse(Name, Title, fx, fy, fyerr){
  Initialize();
}


TTrolleyNMRPulse::~TTrolleyNMRPulse(){
  
}


void TTrolleyNMRPulse::Initialize(){
  ClockFrequency = 0.;
  DigitizerFrequency = 0.;
}


void TTrolleyNMRPulse::SetClockFrequencyInMHz(float f){
  SetClockFrequency(f*1E6);
}


void TTrolleyNMRPulse::SetClockFrequency(float f){
  ClockFrequency = f;
  DigitizerFrequency = ClockFrequency / 62;
}


