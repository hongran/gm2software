/*
 * =====================================================================================
 *       Filename:  gm2_TEncoder.h
 *
 *    Description:  Encoder readout classes for the GM2 experiment
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

class gm2_TEncoder : public TNamed
{
  //Attributes
  protected:
    vector<double> fTimeStamp;  //time stamp
    vector<double> fPosition1;  //Encoder1 position value
    vector<double> fPosition2;  //Encoder2 position value
    vector<double> fVelocity1;	//Encoder1 velocity value
    vector<double> fVelocity2;	//Encoder2 velocity value
    vector<double> fTension1;	//Encoder1 tension value
    vector<double> fTension2;	//Encoder2 tension value
    vector<double> fControl1;	//Encoder1 control voltage value
    vector<double> fControl2;	//Encoder2 control voltage value

    int fNPoints;	//Number of samples

  public:
    gm2_TEncoder(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
    ~gm2_TEncoder();
    //Set Methods
    void SetPoint(const int i,const double Time,const double P1,const double P2,const double V1,const double V2,const double T1,const double T2,const double C1,const double C2);
    //Get Methods
    int GetNPoints() const{return fNPoints;}
    vector<double> GetTimeLine() const{return fTimeStamp;}
    vector<double> GetPosition(int channel) const;
    vector<double> GetVelocity(int channel) const;
    vector<double> GetTension(int channel) const;
    vector<double> GetControlVoltage(int channel) const;

    shared_ptr<TGraph> GetPositionGraph(int channel,double shift=0) const;
    shared_ptr<TGraph> GetVelocityGraph(int channel,double shift=0) const;
    shared_ptr<TGraph> GetTensionGraph(int channel,double shift=0) const;
    shared_ptr<TGraph> GetControlVoltageGraph(int channel,double shift=0) const;
};
