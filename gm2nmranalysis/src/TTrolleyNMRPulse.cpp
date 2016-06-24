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
  
  inputFileName = "";
  inputFile = NULL;
}


void TTrolleyNMRPulse::SetClockFrequencyInMHz(float f){
  SetClockFrequency(f*1E6);
}


void TTrolleyNMRPulse::SetClockFrequency(float f){
  ClockFrequency = f;
  DigitizerFrequency = ClockFrequency / 62;
}


void TTrolleyNMRPulse::SetInputFile(string input){
  inputFileName.assign(input);
}


int TTrolleyNMRPulse::ReadNextPulse(){
  Packet_Header header;
  unsigned short Nsamples = 0;
  short sample;
  int j = 0;
  
  if(!inputFile){ 
    inputFile = fopen(inputFileName.c_str(), "r");
  }
  
  if(!inputFile){
    cerr << "TTrolleyNMRPulse::ReadNextPulse(): Couldn't open file " << inputFileName.c_str() << endl;
    return -1;
  }

  fX.clear();
  fY.clear();
  
  while(1){
    if(feof(inputFile)) return -2;
    if(fread(&header, sizeof(header), 1, inputFile)){
      
      // Let's look for the very specific start pattern for a new header 
      if(header.fixed_value1 == 0x80007FFF && 
	 header.fixed_value2 == 0x80007FFF &&
	 header.fixed_value3 == 0x80007FFF &&
	 header.fixed_value4 == 0xAAAA){
	
	cout << "New header found" << endl;

	// We have found a header, so we can decode the data
	Nsamples = header.nmr_samples;
	if(Nsamples == 0) return -3;
	
	for(int i=0; i<Nsamples; ++i){
	  fread(&sample, sizeof(sample), 1, inputFile);
	  SetPoint(i, 1./DigitizerFrequency * (double)i, (double)sample, 1.);
	}
	
	break;
      }
      else{ // Still looking for the header word, so let's spool back by 2 bytes 
	fseek(inputFile, -sizeof(header)+2, SEEK_CUR);
      }
    }
  }
  
  return 1;
}
