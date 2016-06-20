/*
 * =====================================================================================
 *
 *       Filename:  gm2_TBarcode.cxx
 *
 *    Description:  gm2_TBarcode classes function definitions.
 *
 *        Version:  1.0
 *        Created:  03/24/2016 10:21:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ran Hong 
 *   Organization:  ANL
 *
 * =====================================================================================
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include "gm2_TBarcode.h"
#include "TString.h"
#include "TVirtualPad.h"

using namespace std;

/**********************************************************************/
//Base class TBarcode
/**********************************************************************/
gm2_TBarcode::gm2_TBarcode(const TString& Name, const TString& Title) : TNamed(Name,Title)
{
  //Initialize attributes
  fNPoints = 0;		      //Number of points
  fNHighLevels = 0;           //Number of high levels
  fNLowLevels = 0;            //Number of low levels
  fNExtrema = 0;              //Number o extrema
  fThreshold = 0.5;           //Threshold
  //Flags
  LogicLevelConverted = false;   //Whether logic levels are converted
  ExtremaFound = false;		//Whether Extrema are found
}

/**********************************************************************/
gm2_TBarcode::gm2_TBarcode(const TString& Name, const TString& Title , vector<double> fx, vector<double> fy) : TNamed(Name,Title)
{
  //Initialize attributes
  fNHighLevels = 0;           //Number of high levels
  fNLowLevels = 0;            //Number of low levels
  fNExtrema = 0;              //Number o extrema
  fThreshold = 0.5;           //Threshold
  //Initialize data array
  if (fX.size()!=fY.size()){
    cerr <<"Error in constructing TBarcode. Dimensions do not fit."<<endl;
    return;
  }
  fX = fx;
  fY = fy;
  fNPoints = fX.size();
  //Flags
  LogicLevelConverted = false;   //Whether logic levels are converted
  ExtremaFound = false;		//Whether Extrema are found
}

/**********************************************************************/
gm2_TBarcode::~gm2_TBarcode()
{}

/**********************************************************************/
void gm2_TBarcode::SetPoint(const int i, const double x,const double y)
{
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

/**********************************************************************/
int gm2_TBarcode::FindHalfRise(int low, int high)
{
  double Difference = fY[high]-fY[low];
  int i=low;
  for (i=low;i<=high;i++){
    if (fY[i]-fY[low]>Difference/2.0)break;
  }
  return i;
}

/**********************************************************************/
int gm2_TBarcode::FindHalfFall(int high, int low)
{
  double Difference = fY[high]-fY[low];
  int i=high;
  for (i=high;i<=low;i++){
    if (fY[i]-fY[low]<Difference/2.0)break;
  }
  return i;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TBarcode::GetRawGraph() const
{
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i],fY[i]);
  }
  graph_ptr->SetName("g"+fName);
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TBarcode::GetLogicLevelGraph() const
{
  if (!LogicLevelConverted){
    cout << "Logic levels for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i],fLogicLevels[i]);
  }

  graph_ptr->SetName("g"+fName+"_LogicLevels");
  graph_ptr->SetTitle(fTitle+"_LogicLevels");
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TBarcode::GetExtremaGraph(TString Option) const
{
  if (!ExtremaFound){
    cout << "Extrema for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNExtrema);
  for (int i=0;i<fNExtrema;i++){
    if (Option.CompareTo("VsIndex")==0){
      graph_ptr->SetPoint(i,i,fY[fExtremaList[i]]);
    }else if (Option.CompareTo("VsX")==0){
      graph_ptr->SetPoint(i,fX[fExtremaList[i]],fY[fExtremaList[i]]);
    }else{
      cout << "Option "<<Option<<" is not identified. Set to VsX."<<endl;
      graph_ptr->SetPoint(i,fX[i],fY[i]);
    }
  }
  graph_ptr->SetName("g"+fName+"_Extrema");
  graph_ptr->SetTitle(fTitle+"_Extrema");
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TBarcode::GetIntervalGraph() const
{
  if (!ExtremaFound){
    cout << "Extrema for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNExtrema-1);
  for (int i=0;i<fNExtrema-1;i++){
    graph_ptr->SetPoint(i,fX[fExtremaList[i]],fX[fExtremaList[i+1]]-fX[fExtremaList[i]]);
  }
  graph_ptr->SetName("g"+fName+"_Interval");
  graph_ptr->SetTitle(fTitle+"_Interval");
  return graph_ptr;
}

/**********************************************************************/
//Derived class TRegBarcode
/**********************************************************************/
gm2_TRegBarcode::gm2_TRegBarcode(const TString& Name, const TString& Title):gm2_TBarcode(Name,Title)
{}

/**********************************************************************/
gm2_TRegBarcode::gm2_TRegBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy):gm2_TBarcode(Name,Title,fx,fy)
{}

/**********************************************************************/
gm2_TRegBarcode::~gm2_TRegBarcode()
{}

/**********************************************************************/
shared_ptr<TGraph> gm2_TRegBarcode::GetAverageGraph() const
{
  if (!ExtremaFound){
    cout << "Extrema for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i],fAverage[i]);
  }

  graph_ptr->SetName("g"+fName+"_Average");
  graph_ptr->SetTitle(fTitle+"_Average");
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TRegBarcode::GetContrastGraph() const
{
  if (!ExtremaFound){
    cout << "Extrema for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i],fContrast[i]);
  }

  graph_ptr->SetName("g"+fName+"_Contrast");
  graph_ptr->SetTitle(fTitle+"_Contrast");
  return graph_ptr;
}

/**********************************************************************/
int gm2_TRegBarcode::FindExtrema()
{
  auto N = fX.size();
  int MaxN = 0;
  int MinN = 0;
  fNExtrema = 0;
  fMaxList.clear();
  fMinList.clear();

  //Read out values
  double y1,y2;

  //Contrast threshold method
  //To determine a minimum or maximum, the depth of the peak of trough needs to be larger than the Threshold
  int MaxIndex = 0;
  double LocalMax = -1E6;
  int MinIndex = 0;
  double LocalMin = 1E6;

  //Find initial trend
  bool LastWasMax = false;
  if(N > 4){
    y1 = fY[0];
    y2 = fY[4];
    if(y1 > y2){
      LastWasMax = true;
    }else{
      LastWasMax = false;
    }
  }else{
    cout <<"Barcode length is too short."<<endl;
    return -1;
  }

  for(int i=0; i<N-1; i++){
    y1 = fY[i];
    y2 = fY[i+1];

    if(y1 < LocalMin && LastWasMax){
      MinIndex = i;
      LocalMin = y1;
    }
    if(y1 > LocalMax && !LastWasMax){
      MaxIndex = i;
      LocalMax = y1;
    }

    if(y1 > LocalMin && y2 > LocalMin+fThreshold && y2 > y1 && LastWasMax){
      //	cout <<"min "<<MinIndex<<" "<<LocalMin<<endl;
      fMinList.push_back(MinIndex);
      fExtremaList.push_back(MinIndex);
      MinN++;
      LocalMin = 1E6;
      LastWasMax = false;
    }
    if(y1 < LocalMax && y2 < LocalMax-fThreshold && y2 < y1 && !LastWasMax){
      //	cout <<"max "<<MaxIndex<<" "<<LocalMax<<endl;
      fMaxList.push_back(MaxIndex);
      fExtremaList.push_back(MaxIndex);
      MaxN++;
      LocalMax = -1E6;
      LastWasMax = true;
    }
  }
  fNExtrema = MaxN + MinN;
  ExtremaFound = true;

  //Calculate Average and Contrast
  fAverage.resize(fNPoints,-1);
  fContrast.resize(fNPoints,-1);
  for (int i=0;i<fNExtrema-1;i++){
    double average_val = (fY[fExtremaList[i]]+fY[fExtremaList[i+1]])/2.0;
    if (average_val<=0){
      cerr<<"Error! Average Barcode signal is zero!"<<endl;
      return -1;
    }
    double contrast_val = abs(fY[fExtremaList[i]]-fY[fExtremaList[i+1]])/average_val/2.0;
    for (int j=fExtremaList[i];j<fExtremaList[i+1];j++){
      fAverage[j] = average_val;
      fContrast[j] = contrast_val;
    }
  }
  //handle the beginning section and ending section
  for (int j=0;j<fExtremaList[0];j++){
    fAverage[j] = fAverage[fExtremaList[0]];
    fContrast[j] = fContrast[fExtremaList[0]];
  }
  for (int j=fExtremaList[fNExtrema-1];j<N;j++){
    fAverage[j] = fAverage[fExtremaList[fNExtrema-2]];
    fContrast[j] = fContrast[fExtremaList[fNExtrema-2]];
  }
  return 0;
}

/**********************************************************************/
int gm2_TRegBarcode::ConvertToLogic()
{
  fLogicLevels.clear();
  fLogicLevels.resize(fNPoints,-1); //Initialize to -1 for safety check
  fNHighLevels = 0;
  fNLowLevels = 0;
  //If extrema are not found, find them first
  if (!ExtremaFound){
    FindExtrema();
  }

  //Initialize to -1, for error checking after processing the data
  for (int i=0;i<fNPoints;i++){
    fLogicLevels[i]=-1;
  }

  //Convert analog to lotic levels
  int Level = 0;
  int IndexLow = 0;
  int IndexHigh = 0;
  //Treat first incomplete section
  IndexLow = 0;
//  IndexHigh = (fMaxList[0]+fMinList[0])/2;
  if (fExtremaList[0] == fMinList[0]){
    IndexHigh = FindHalfRise(fExtremaList[0],fExtremaList[1]);
    Level = 0;
    fNLowLevels++;
  }else if (fExtremaList[0] == fMaxList[0]){
    IndexHigh = FindHalfFall(fExtremaList[0],fExtremaList[1]);
    Level = 1;
    fNHighLevels++;
  }else{
    cout << "Encountered mismatching of Max/Min list and Extrema list."<<endl;
    return -1;
  }
  for (int i=IndexLow;i<IndexHigh;i++){
    fLogicLevels[i]=Level;
  }
  //Alternating Level
  if (Level==0)Level=1;
  else if (Level==1)Level=0;
  //Treat normal section
  for (int i=1;i<fNExtrema-1;i++){
    if (Level==0){
      IndexLow = FindHalfFall(fExtremaList[i-1],fExtremaList[i]);
      IndexHigh = FindHalfRise(fExtremaList[i],fExtremaList[i+1]);
      fNLowLevels++;
    }else if(Level==1){
      IndexLow = FindHalfRise(fExtremaList[i-1],fExtremaList[i]);
      IndexHigh = FindHalfFall(fExtremaList[i],fExtremaList[i+1]);
      fNHighLevels++;
    }
    if (IndexHigh>=fNPoints){
      cout << "Error! fLogicLevels out of range!"<<endl;
      return -1;
    }
    for (int i=IndexLow;i<IndexHigh;i++){
      fLogicLevels[i]=Level;
    }
    if (Level==0)Level=1;
    else if (Level==1)Level=0;
  }
  //Treat last step
  IndexHigh = fNPoints-1;
  if (Level==0){
    IndexLow = FindHalfFall(fExtremaList[fNExtrema-2],fExtremaList[fNExtrema-1]);
    fNLowLevels++;
  }else if(Level==1){
    IndexLow = FindHalfRise(fExtremaList[fNExtrema-2],fExtremaList[fNExtrema-1]);
    fNHighLevels++;
  }
  for (int i=IndexLow;i<=IndexHigh;i++){
    fLogicLevels[i]=Level;
  }

  //Testing and throwing out warnings and errors
  if (fNExtrema!=(fNHighLevels+fNLowLevels)){
    cout <<"Warning! fNExtrema != fNHighLevels+fNLowLevels. fNExtrema="<<fNExtrema<<" fNHighLevels+fNLowLevels="<<fNHighLevels+fNLowLevels<<endl;
  }
  /*
  for (int i=0;i<fNPoints;i++){
    if (fLogicLevels[i]==-1){
      cout <<"Error! fLogicLevels is not completely filled at "<<i<<endl;
      return -1;
    }
  }*/
  //Try using lambda
  auto count = count_if(begin(fExtremaList), end(fExtremaList), [](int i){ return i==-1; } );
  if (count>0){
    cout <<"Error! fLogicLevels is not completely filled"<<endl;
    return -1;
  }
  LogicLevelConverted = true;

  return 0;
}

/**********************************************************************/
//Derived class TAbsBarcode
/**********************************************************************/
gm2_TAbsBarcode::gm2_TAbsBarcode(const TString& Name, const TString& Title):gm2_TBarcode(Name,Title)
{}

/**********************************************************************/
gm2_TAbsBarcode::gm2_TAbsBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy):gm2_TBarcode(Name,Title,fx,fy)
{}

/**********************************************************************/
gm2_TAbsBarcode::~gm2_TAbsBarcode()
{}

/**********************************************************************/
int gm2_TAbsBarcode::FindExtrema()
{
  //Not yet developed
  return 0;
}

/**********************************************************************/
int gm2_TAbsBarcode::ConvertToLogic()
{
  //Not yet developed
  return 0;
}

/**********************************************************************/
int gm2_TAbsBarcode::Decode(const gm2_TRegBarcode& RefReg)
{
  //Not yet developed
  return 0;
}

