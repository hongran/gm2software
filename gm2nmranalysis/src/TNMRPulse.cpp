#include "TNMRPulse.h"

TNMRPulse::TNMRPulse(const TString& Name, const TString& Title) : TNamed(Name, Title)
{
  //Initialize attributes
  fNPoints = 0;		      //Number of points
}

TNMRPulse::TNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy) : TNamed(Name, Title)
{
  //Initialize data array
  if (fX.size()!=fY.size()){
    cerr <<"Error in constructing TNMRPulse. Dimensions do not fit."<<endl;
    return;
  }
  fX = fx;
  fY = fy;
  fNPoints = fX.size();
}

TNMRPulse::~TNMRPulse(){
  
}

void TNMRPulse::SetPoint(const int i, const double x, const double y){
  if (i<fX.size() && i<fY.size()){
    fX[i] = x;
    fY[i] = y;
  }else{
    fX.resize(i+1,-1.0);
    fY.resize(i+1,-1.0);
    fX[i] = x;
    fY[i] = y;
    fNPoints = i+1;
  }
}
