/*
 * =====================================================================================
 *       Filename:  gm2_TBarcode.cxx
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

using namespace std;

class gm2_TBarcode : public TNamed
{
  //Attributes
  protected:
    vector<double> fX;	//time stamp
    vector<double> fY;  //barcode readout value
    vector<int> fLogicLevels;	//logic levels converted from the analog signal
    vector<int> fMaxList;	//List of maximum time stamps
    vector<int> fMinList;	//List of minimum time stamps
    vector<int> fExtremaList;	//List of extremum time stamps
    int fNPoints;		//Number of points
    int fNHighLevels;		//Number of high levels
    int fNLowLevels;		//Number of low levels
    int fNExtrema;		//Number o extrema
    double fThreshold;		//Threshold
    //Flags
    bool LogicLevelConverted;	//Whether logic levels are converted
    bool ExtremaFound;		//Whether Extrema are found

    //Protected methods
    int FindHalfRise(int low, int high);	//Find half-way rising edge
    int FindHalfFall(int high, int low);	//Find half-way falling edge

  public:
    gm2_TBarcode(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    gm2_TBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
    ~gm2_TBarcode();
    //Set Methods
    void SetThreshold(double val){fThreshold = val;}
    void SetPoint(const int i, const double x,const double y);
    //Get Methods
    int GetNPoints() const{return fNPoints;}
    int GetNHighLevels() const{return fNHighLevels;}
    int GetNLowLevels() const{return fNLowLevels;}
    int GetNExtrema() const{return fNExtrema;}
    double GetThreshold() const{return fThreshold;}
    vector<double> GetX() const{return fX;}
    vector<double> GetY() const{return fY;}
    vector<int> GetExtremaList() const{return fExtremaList;}
    vector<int> GetLogicLevels() const{return fLogicLevels;}
    shared_ptr<TGraph> GetRawGraph() const;
    shared_ptr<TGraph> GetLogicLevelGraph() const;
    shared_ptr<TGraph> GetExtremaGraph(TString Option) const;
    shared_ptr<TGraph> GetIntervalGraph()const;
    
    //virtual method for determining the extrema
    virtual int FindExtrema() = 0;
    virtual int ConvertToLogic() = 0;
};

class gm2_TRegBarcode : public gm2_TBarcode
{
  protected:
    vector<double> fAverage;
    vector<double> fContrast;
  public:
    gm2_TRegBarcode(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    gm2_TRegBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
    ~gm2_TRegBarcode();
    //Get functions
    vector<double> GetAverage() const{return fAverage;}
    vector<double> GetContrast() const{return fContrast;}
    shared_ptr<TGraph> GetAverageGraph() const;
    shared_ptr<TGraph> GetContrastGraph() const;
    //method for determining the extrema
    int FindExtrema();
    //Convert analog signal to logic levels
    int ConvertToLogic();
};

class gm2_TAbsBarcode : public gm2_TBarcode
{
  protected:
    int fNGroupIndex = 0;
    vector<int> fGroupIndexList;
    vector<int> fGroupStartList;
  public:
    gm2_TAbsBarcode(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    gm2_TAbsBarcode(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
    ~gm2_TAbsBarcode();
    //method for determining the extrema
    int FindExtrema();
    //Convert analog signal to logic levels
    int ConvertToLogic();
    int Decode(const gm2_TRegBarcode& RefReg);
};

