#include "TNMRPulse.h"

TNMRPulse::TNMRPulse(const TString& Name, const TString& Title) : TNamed(Name, Title)
{
  //Initialize attributes
  fNPoints = 0;		      //Number of points
}

TNMRPulse::TNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy, vector<double> fyerr) : TNamed(Name, Title)
{
  //Initialize data array
  if (fx.size()!=fy.size() || (fyerr.size()>0 && fyerr.size()!=fx.size())){
    cerr <<"Error in constructing TNMRPulse. Dimensions do not fit."<<endl;
    return;
  }
  fX = fx;
  fY = fy;
  
  if(fyerr.size()>0){
    fYErr = fyerr;
  }
  else{
    fYErr.resize(fX.size(), 0.);
  }
  
  fNPoints = fX.size();
}

TNMRPulse::~TNMRPulse(){
  
}

void TNMRPulse::SetPoint(const int i, const double x, const double y, const double yerr){
  if (i<fX.size() && i<fY.size()){
    fX[i] = x;
    fY[i] = y;
    fYErr[i] = yerr;
  }else{
    fX.resize(i+1,-1.0);
    fY.resize(i+1,-1.0);
    fYErr.resize(i+1,-1.0);
    fX[i] = x;
    fY[i] = y;
    fYErr[i] = yerr;
    fNPoints = i+1;
  }
}


shared_ptr<TGraphErrors> TNMRPulse::GetRawGraph() const
{
  auto graph_ptr = make_shared<TGraphErrors>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i],fY[i]);
    graph_ptr->SetPointError(i,0.,fYErr[i]);
  }
  graph_ptr->SetName("g"+fName);
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}
