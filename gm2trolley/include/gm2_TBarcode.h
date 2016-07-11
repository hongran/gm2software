/*
 * =====================================================================================
 *       Filename:  gm2_TBarcode.h
 *
 *    Description:  Barcode classes for the GM2 experiment
 *
 *        Version:  1.0
 *        Created:  03/24/2016 10:21:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ran Hong 
 *   Organization:  ANL
 This class is inherited from 
 the ROOT TNamed class. 
 
 * =====================================================================================
*/
#include <TNamed.h>
#include <TGraph.h>
#include <TString.h>
#include <map>
#include <memory>
#include <TH1.h>

using namespace std;

namespace gm2Barcode{
  class LogicLevel{
    public:
      int Level;
      int LEdge;	//Index of left edge
      int REdge;	//Index of right edge
  };
  class AbsBarcodeSegment{
    public:
      int RegionType;		//Code region (1) or non-code region(0) or abnormal region(-1)
      int Code;			//Code represented by the binary number
      int NAbsLevel;		//Number of Abs levels
      int NRegLevel;		//Number of Reg levels
      vector<int> fLevelIndexList;	//List of Abs logic level indicies
      vector<LogicLevel> fRegLevelList;	//List of Reg logic levels
  };
  enum AbnormalStatus{
    None,
    ZeroVelocity,
    StartRegion,
    EndRegion,
    ChamberJunction
  };
}

class gm2_TBarcode : public TNamed
{
  //Attributes
  protected:
    vector<double> fX;	//time stamp
    vector<double> fY;  //barcode readout value
    vector<int> fDirection;	//Direction indicators;
    vector<gm2Barcode::AbnormalStatus> fAbnormal;
    vector<gm2Barcode::LogicLevel> fLogicLevels;	//logic levels converted from the analog signal
    vector<int> fMaxList;	//List of maximum time stamps
    vector<int> fMinList;	//List of minimum time stamps
    vector<int> fExtremaList;	//List of extremum time stamps
    int fNPoints;		//Number of points
    int fNHighLevels;		//Number of high levels
    int fNLowLevels;		//Number of low levels
    int fNExtrema;		//Number o extrema
    double fThreshold;		//Threshold
    double fLogicLevelScale;	//Value of logic-1 for plotting graph
    double fTransitionThreshold;//Threshold for determining the transition between high and low level
    //Flags
    bool LogicLevelConverted;	//Whether logic levels are converted
    bool ExtremaFound;		//Whether Extrema are found
    bool DirectionSet;		//Whether direction list is set

    //Protected methods
    int FindHalfRise(int low, int high);	//Find half-way rising edge
    int FindHalfFall(int high, int low);	//Find half-way falling edge

  public:
    gm2_TBarcode(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    gm2_TBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
    ~gm2_TBarcode();
    //Set Methods
    void SetThreshold(const double val){fThreshold = val;}
    void SetTransitionThreshold(const double val);
    void SetPoint(const int i, const double x,const double y);
    void SetLogicLevelScale(const double Level){fLogicLevelScale = Level;}
    void SetDirection(const vector<int> ExternalDirection);
    //Get Methods
    bool IfLogicLevelConverted() const{return LogicLevelConverted;}
    bool IfExtremaFound() const{return ExtremaFound;}
    bool IfDirectionSet() const{return DirectionSet;}
    int GetNPoints() const{return fNPoints;}
    int GetNHighLevels() const{return fNHighLevels;}
    int GetNLowLevels() const{return fNLowLevels;}
    int GetNExtrema() const{return fNExtrema;}
    double GetThreshold() const{return fThreshold;}
    vector<double> GetX() const{return fX;}
    vector<double> GetY() const{return fY;}
    vector<int> GetExtremaList() const{return fExtremaList;}
    vector<gm2Barcode::LogicLevel> GetLogicLevels() const{return fLogicLevels;}
    shared_ptr<TGraph> GetRawGraph(double shift=0) const;
    shared_ptr<TGraph> GetLogicLevelGraph(double shift=0) const;
    shared_ptr<TGraph> GetExtremaGraph(TString Option,double shift=0) const;
    shared_ptr<TGraph> GetDirectionGraph(double shift=0) const;
    shared_ptr<TGraph> GetIntervalGraph()const;
    shared_ptr<TGraph> GetLevelWidthGraph()const;
    shared_ptr<TH1D> GetLevelWidthHist(string LevelSelection)const;
    void Smooth(int PointsAveraged = 5);
    
    //virtual method for determining the extrema
    virtual int FindExtrema() = 0;
    //virtual int ConvertToLogic() = 0;
    int ConvertToLogic();
};

class gm2_TRegBarcode : public gm2_TBarcode
{
  protected:
    vector<double> fAverage;
    vector<double> fContrast;
    vector<int> fBigGapList;
  public:
    gm2_TRegBarcode(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    gm2_TRegBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
    ~gm2_TRegBarcode();
    //Get functions
    vector<double> GetAverage() const{return fAverage;}
    vector<double> GetContrast() const{return fContrast;}
    shared_ptr<TGraph> GetAverageGraph() const;
    shared_ptr<TGraph> GetContrastGraph() const;
    shared_ptr<TGraph> GetVelocityGraph()const;
    //method for determining the extrema
    int FindExtrema();
    int FindBigGaps();
};

class gm2_TAbsBarcode : public gm2_TBarcode
{
  protected:
    int fNSegments;
    vector<vector<gm2Barcode::LogicLevel>> fAuxList;	//Auxilliary list of RegLevels associated with each AbsLevel
    vector<gm2Barcode::AbsBarcodeSegment> fSegmentList;
    bool fSegmented;
  public:
    gm2_TAbsBarcode(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    gm2_TAbsBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
    ~gm2_TAbsBarcode();
    //Get Methods
    shared_ptr<TGraph> GetAbsWidthGraph()const;
    shared_ptr<TGraph> GetAbsSegWidthGraph()const;
    //method for determining the extrema
    int FindExtrema();
    int ChopSegments(const gm2_TRegBarcode& RefReg);
    int Decode();
};

