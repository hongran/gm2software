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
  fLogicLevelScale = 1.0;     //By default, logic-1 is "1"
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
void gm2_TBarcode::SetDirection(const vector<int> ExternalDirectionList)
{
  fDirectionList = ExternalDirectionList;
  if (ExternalDirectionList.size()!=fNPoints){
    cout << "Warning! Direction list has a different size. Auto resized to fNpoints."<<endl;
    fDirectionList.resize(fNPoints,0);
  }
  DirectionSet = true;
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
shared_ptr<TGraph> gm2_TBarcode::GetRawGraph(double shift) const
{
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i]+shift,fY[i]);
  }
  graph_ptr->SetName("g"+fName);
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TBarcode::GetLogicLevelGraph(double shift) const
{
  if (!LogicLevelConverted){
    cout << "Logic levels for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(2*fNExtrema);
  for (int i=0;i<fNExtrema;i++){
    graph_ptr->SetPoint(2*i,fX[fLogicLevels[i].LEdge]+shift,fLogicLevels[i].Level*fLogicLevelScale);
    graph_ptr->SetPoint(2*i+1,fX[fLogicLevels[i].REdge]+shift,fLogicLevels[i].Level*fLogicLevelScale);
  }

  graph_ptr->SetName("g"+fName+"_LogicLevels");
  graph_ptr->SetTitle(fTitle+"_LogicLevels");
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TBarcode::GetExtremaGraph(TString Option,double shift) const
{
  if (!ExtremaFound){
    cout << "Extrema for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNExtrema);
  for (int i=0;i<fNExtrema;i++){
    if (Option.CompareTo("VsIndex")==0){
      graph_ptr->SetPoint(i,i+shift,fY[fExtremaList[i]]);
    }else if (Option.CompareTo("VsX")==0){
      graph_ptr->SetPoint(i,fX[fExtremaList[i]]+shift,fY[fExtremaList[i]]);
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
shared_ptr<TGraph> gm2_TBarcode::GetDirectionGraph(double shift) const
{
  if (!DirectionSet){
    cout << "Direction is not set."<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    graph_ptr->SetPoint(i,fX[i]+shift,fDirectionList[i]);
  }

  graph_ptr->SetName("g"+fName+"_Direction");
  graph_ptr->SetTitle(fTitle+"_Direction");
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
shared_ptr<TGraph> gm2_TBarcode::GetLevelWidthGraph() const
{
  if (!LogicLevelConverted){
    cout <<"Not yet converted to Logic levels."<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNExtrema-1);
  for (int i=0;i<fNExtrema;i++){
    graph_ptr->SetPoint(i,(fX[fLogicLevels[i].REdge]+fX[fLogicLevels[i].LEdge])/2.0,fX[fLogicLevels[i].REdge]-fX[fLogicLevels[i].LEdge]);
  }
  graph_ptr->SetName("g"+fName+"_LevelWidths");
  graph_ptr->SetTitle(fTitle+"_LevelWidths");
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TH1D> gm2_TBarcode::GetLevelWidthHist(string LevelSelection)const
{
  if (!LogicLevelConverted){
    cout <<"Not yet converted to Logic levels."<<endl;
    return nullptr;
  }
  bool IncludeHigh{false},IncludeLow{false};
  if (LevelSelection.compare("High")==0)IncludeHigh=true;
  else if (LevelSelection.compare("Low")==0)IncludeLow=true;
  else if (LevelSelection.compare("Both")==0){
    IncludeHigh=true;
    IncludeLow=true;
  }else{
    cout <<"Non-recognized option "<<LevelSelection<<". Set to Both."<<endl;
    LevelSelection="Both";
    IncludeHigh=true;
    IncludeLow=true;
  }
  auto Name = string{"Hist_"+fName+LevelSelection+"_LevelWidths"};
  double Upperlimit = 5*(fX[fLogicLevels[4].REdge]-fX[fLogicLevels[4].LEdge]);
  auto hist_ptr = make_shared<TH1D>(Name.c_str(),Name.c_str(),100,0,Upperlimit);
  for (int i=0;i<fNExtrema;i++){
    if(IncludeHigh && (fLogicLevels[i].Level==1))hist_ptr->Fill(fX[fLogicLevels[i].REdge]-fX[fLogicLevels[i].LEdge]);
    if(IncludeLow && (fLogicLevels[i].Level==0))hist_ptr->Fill(fX[fLogicLevels[i].REdge]-fX[fLogicLevels[i].LEdge]);
  }
  return hist_ptr;
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
shared_ptr<TGraph> gm2_TRegBarcode::GetVelocityGraph() const
{
  if (!ExtremaFound){
    cout << "Extrema for Barcode "<<fName<<" are not constructed!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNExtrema-1);
  for (int i=0;i<fNExtrema-1;i++){
    double delta_T = fX[fExtremaList[i+1]]-fX[fExtremaList[i]];//ms
    double delta_X = 2.0; //mm
    graph_ptr->SetPoint(i,fX[fExtremaList[i]],delta_X/delta_T*10000.0); // mm/s
  }
  graph_ptr->SetName("g"+fName+"_Velocity");
  graph_ptr->SetTitle(fTitle+"_Velocity");
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
int gm2_TBarcode::ConvertToLogic()
{
  //If extrema are not found, find them first
  if (!ExtremaFound){
    FindExtrema();
  }
  //Initialize
  fLogicLevels.clear();
  fLogicLevels.resize(fNExtrema,gm2Barcode::LogicLevel{-1,0,0}); //Initialize to -1 for safety check
  fNHighLevels = 0;
  fNLowLevels = 0;

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
  fLogicLevels[0] = gm2Barcode::LogicLevel{Level,IndexLow,IndexHigh};
/*  for (int i=IndexLow;i<IndexHigh;i++){
    fLogicLevels[i]=Level;
  }*/
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
    fLogicLevels[i] = gm2Barcode::LogicLevel{Level,IndexLow,IndexHigh};
    /*
    for (int i=IndexLow;i<IndexHigh;i++){
      fLogicLevels[i]=Level;
    }*/
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
  fLogicLevels[fNExtrema-1] = gm2Barcode::LogicLevel{Level,IndexLow,IndexHigh};
/*  for (int i=IndexLow;i<=IndexHigh;i++){
    fLogicLevels[i]=Level;
  }*/

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
  auto count = count_if(begin(fLogicLevels), end(fLogicLevels), [](gm2Barcode::LogicLevel i){ return i.Level==-1; } );
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
{
  fNSegments = 0;
  fSegmented = false;
}

/**********************************************************************/
gm2_TAbsBarcode::gm2_TAbsBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy):gm2_TBarcode(Name,Title,fx,fy)
{
  fNSegments = 0;
  fSegmented = false;
}

/**********************************************************************/
gm2_TAbsBarcode::~gm2_TAbsBarcode()
{}

/**********************************************************************/
/*int gm2_TAbsBarcode::FindExtrema()
{
  //Not yet developed
  return 0;
}*/
int gm2_TAbsBarcode::FindExtrema()
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
  return 0;
}

/**********************************************************************/
/*int gm2_TAbsBarcode::ConvertToLogic()
{
  //Not yet developed
  return 0;
}*/

/**********************************************************************/
int gm2_TAbsBarcode::ChopSegments(const gm2_TRegBarcode& RefReg)
{
  //Clean up the previous segmentations
  if (fSegmented){
    fSegmentList.clear();
    fAuxList.clear();
    fSegmented=false;
  }
  //Check if the logic levels are converted
  if (!LogicLevelConverted){
    ConvertToLogic();
  }
  if (!RefReg.IfLogicLevelConverted()){
    cout <<"Error! Reference RegBarcode is not yet converted to Logic levels."<<endl;
    return -1;
  }
  fAuxList.resize(fNExtrema);

  //Get Logic level list from RefReg
  auto RefLevels = RefReg.GetLogicLevels();
  auto NRefLevels = RefLevels.size();
  int i{0};
  int j{0};
  while (i<fNExtrema && j<NRefLevels){
    auto LeftBound = fLogicLevels[i].LEdge;
    auto RightBound = fLogicLevels[i].REdge;
    while (j<NRefLevels){
      if (RefLevels[j].REdge<LeftBound){
	j++;
	continue;
      }
      else{
	fAuxList[i].push_back(RefLevels[j]);
	if (RefLevels[j].REdge<RightBound){
	  j++;
	  continue;
	}
	else{
	  break;
	}
      }
    }
    i++;
  }
  
  //Constructs segments
  i=0;
  while (i<fNExtrema){
    int NLevelsRef = fAuxList[i].size();
    if (NLevelsRef>=16 && NLevelsRef<=20){
      fSegmentList.push_back(gm2Barcode::AbsBarcodeSegment{false,-1,1,NLevelsRef,vector<int>{i},fAuxList[i]});
      i++;
      continue;
    }else{
      vector<int> TempAbsList;
      vector<gm2Barcode::LogicLevel> TempRegList;
      while (i<fNExtrema){
	int NLevelsRef = fAuxList[i].size();
	if (NLevelsRef>=16 && NLevelsRef<=20)break;
	TempAbsList.push_back(i);
	for (j=0;j<fAuxList[i].size();j++){
	  if (TempRegList.size()!=0){
	    if ((TempRegList.back().LEdge==fAuxList[i][j].LEdge) && (TempRegList.back().REdge==fAuxList[i][j].REdge)){
	      continue;
	    }
	  }
	  TempRegList.push_back(fAuxList[i][j]);
	}
	i++;
      }
      fSegmentList.push_back(gm2Barcode::AbsBarcodeSegment{true,-1,static_cast<int>(TempAbsList.size()),static_cast<int>(TempRegList.size()),TempAbsList,TempRegList});
      cout << "Temp Reg List Size: " << TempRegList.size() << endl;
    }
  }

  fNSegments = fSegmentList.size();
  fSegmented=true;
  return 0;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TAbsBarcode::GetAbsWidthGraph() const
{
  if (!fSegmented){
    cout << "Abs Barcode "<<fName<<" is not segmented!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNExtrema);
  for (int i=0;i<fNExtrema;i++){
    graph_ptr->SetPoint(i,(fX[fLogicLevels[i].REdge]+fX[fLogicLevels[i].LEdge])/2.0,fAuxList[i].size());
  }
  graph_ptr->SetName("g"+fName+"_AbsWidths");
  graph_ptr->SetTitle(fTitle+"_AbsWidths");
  return graph_ptr;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TAbsBarcode::GetAbsSegWidthGraph() const
{
  if (!fSegmented){
    cout << "Abs Barcode "<<fName<<" is not segmented!"<<endl;
    return nullptr;
  }
  auto graph_ptr = make_shared<TGraph>(fNSegments);
  for (int i=0;i<fNSegments;i++){
    graph_ptr->SetPoint(i,i,fSegmentList[i].NRegLevel);
  }
  graph_ptr->SetName("g"+fName+"_AbsSegWidths");
  graph_ptr->SetTitle(fTitle+"_AbsSegWidths");
  return graph_ptr;
}

/**********************************************************************/
int gm2_TAbsBarcode::Decode()
{
  //Not yet developed
  for (int i=0; i<fNSegments; ++i){
  	if(!fSegmentList[i].IsCodeRegion){
  		continue;
  	}
  	auto AbsIndexList = fSegmentList[i].fLevelIndexList;
  	auto RegLevelList = fSegmentList[i].fRegLevelList;
  	RegLevelList.erase(RegLevelList.begin()); //Remove first element
  	RegLevelList.pop_back();                  //Remove last element
  	vector<int> DigitMap(RegLevelList.size(),-1); //Temporarily save the digits from the absolute barcode
  	
  	for (int j=0; j<RegLevelList.size(); ++j){
  		auto RegLevelUnit = RegLevelList[j];
  		for (int k=0; k<AbsIndexList.size(); ++k){
  			auto AbsLevelIndex = AbsIndexList[k];
  			auto AbsLevelUnit = fLogicLevels[AbsLevelIndex];
  			
  			if (RegLevelUnit.REdge <= AbsLevelUnit.LEdge){
  				continue;
  			}else if(AbsLevelUnit.REdge <= RegLevelUnit.LEdge){
  				continue;
  			}else if(RegLevelUnit.LEdge <= AbsLevelUnit.LEdge && AbsLevelUnit.REdge <= RegLevelUnit.REdge){
  				DigitMap[j] = AbsLevelUnit.Level;
  				break;
  			}else if(RegLevelUnit.LEdge >= AbsLevelUnit.LEdge && AbsLevelUnit.REdge >= RegLevelUnit.REdge){
  				DigitMap[j] = AbsLevelUnit.Level;
  				break;
  			}else if(RegLevelUnit.LEdge < AbsLevelUnit.LEdge && RegLevelUnit.REdge > AbsLevelUnit.LEdge && RegLevelUnit.REdge < AbsLevelUnit.REdge){
  				double OverlapRatio = (fX[RegLevelUnit.REdge]-fX[AbsLevelUnit.LEdge])/(fX[RegLevelUnit.REdge]-fX[RegLevelUnit.LEdge]);
  				if (OverlapRatio < 0.5){
  					continue;
  				}
  				if (OverlapRatio >= 0.5){
  					DigitMap[j] = AbsLevelUnit.Level;
  					break;
  				}
  			}else if(RegLevelUnit.LEdge > AbsLevelUnit.LEdge && RegLevelUnit.LEdge < AbsLevelUnit.REdge && RegLevelUnit.REdge > AbsLevelUnit.REdge){
  			  	double OverlapRatio = (fX[AbsLevelUnit.REdge]-fX[RegLevelUnit.LEdge])/(fX[RegLevelUnit.REdge]-fX[RegLevelUnit.LEdge]);
  				if (OverlapRatio < 0.5){
  					continue;
  				}
  				if (OverlapRatio >= 0.5){
  					DigitMap[j] = AbsLevelUnit.Level;
  					break;
  				}
  			}else{
  				cout << "Error! in Decode Function an exception case occured." << endl;
  				cout << "Regular Left: " << fX[RegLevelUnit.LEdge] << endl;
  				cout << "Regular Right: " << fX[RegLevelUnit.REdge] << endl;
  				cout << "Absolute Left: " << fX[AbsLevelUnit.LEdge] << endl;
  				cout << "Absolute Right: " << fX[AbsLevelUnit.REdge] << endl;
  				return -1;
  			}
  			
  		}
  	}
  	DigitMap.erase(DigitMap.begin()); //Delete first element
  	DigitMap.pop_back();              //Delete last element
  	int DigitMapSize = DigitMap.size();
  	cout << "Digit Map Size: " << DigitMapSize << " Reg Level List Size: " << RegLevelList.size() << endl;
/*  	vector<int> pos_vec;
  	for(unsigned int m=0; m<DigitMap.size(); m++){
  		cout << DigitMap[m] << " ";
  		if(DigitMap[m] == 0){
  			continue;
  		}else if (DigitMap[m] == 1){
  			pos_vec.push_back(m);
  		}
  	}*/
  	
  	//Determine Direction
  	int LDir = fDirectionList[fLogicLevels[fSegmentList[i].fLevelIndexList.front()].LEdge];
  	int RDir = fDirectionList[fLogicLevels[fSegmentList[i].fLevelIndexList.back()].REdge];
  	double TimeStamp = fX[fLogicLevels[fSegmentList[i].fLevelIndexList.front()].LEdge];
  	cout <<"LeftDirection= "<<LDir<<endl;
  	cout <<"RightDirection= "<<RDir<<endl;
  	//Make complement, change direction if moving CW
  	if (LDir!=RDir){
  		continue;
  	}else if(LDir==0 || RDir==0){
  		continue;
  	}else if(LDir==-1 && RDir==-1){
  		for(int n=0; n<DigitMapSize/2; n++){
  			int temp = DigitMap[n];
  			DigitMap[n] = DigitMap[DigitMapSize-1-n];
  			DigitMap[DigitMapSize-1-n] = temp;
  		}
  	}
  	for(int n=0; n<DigitMapSize; n++){
  		if(DigitMap[n] == 0){
  			DigitMap[n] = 1;
  		}else if(DigitMap[n] == 1){
  			DigitMap[n] = 0;
  		}else{
  			cout <<"Error! Missing one or more interpretations of digits in Abs Barcode."<<endl;
  			return -1;
  		}
  	}
  	//Convert binary to digital
  	if (DigitMapSize!=10){//In bad region, end region or something else
  		continue;
  	}
  	int DecimalNumber = 0;
  	for(int j=0; j<DigitMapSize; j++){
  		cout <<DigitMap[j];
  		DecimalNumber += DigitMap[j] << j;
  	} 
  	cout << "Decimal Number: " << DecimalNumber << endl;
  	cout << "TimeStamp "<< TimeStamp << endl << endl;
  }
  
  return 0;
}

