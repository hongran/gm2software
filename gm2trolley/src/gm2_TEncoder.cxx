/*
 * =====================================================================================
 *
 *       Filename:  gm2_TEncoder.cxx
 *
 *    Description:  gm2_TEncoder classes function definitions.
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
#include <string>
#include <sstream>
#include "gm2_TEncoder.h"
#include "TString.h"
#include "TVirtualPad.h"

using namespace std;

/**********************************************************************/
//TEncoder class
/**********************************************************************/

/**********************************************************************/
gm2_TEncoder::gm2_TEncoder(const TString& Name, const TString& Title) : TNamed(Name,Title)
{
  //Initialize attributes
  fNPoints = 0;               //Number of points
}

/**********************************************************************/
gm2_TEncoder::~gm2_TEncoder()
{}

/**********************************************************************/
void gm2_TEncoder::SetPoint(const int i,const double Time,const double P1,const double P2,const double V1,const double V2,const double T1,const double T2,const double C1,const double C2)
{
  if (i>=fNPoints){
    fTimeStamp.resize(i+1,-1.0);
    fPosition1.resize(i+1,-1.0);
    fPosition2.resize(i+1,-1.0);
    fVelocity1.resize(i+1,-1.0);
    fVelocity2.resize(i+1,-1.0);
    fTension1.resize(i+1,-1.0);
    fTension2.resize(i+1,-1.0);
    fControl1.resize(i+1,-1.0);
    fControl2.resize(i+1,-1.0);
    fNPoints = i+1;
  }
  fTimeStamp[i] = Time;
  fPosition1[i] = P1;
  fPosition2[i] = P2;
  fVelocity1[i] = V1;
  fVelocity2[i] = V2;
  fTension1[i] = T1;
  fTension2[i] = T2;
  fControl1[i] = C1;
  fControl2[i] = C2;
}

/**********************************************************************/
vector<double> gm2_TEncoder::GetPosition(int channel) const
{
  if (channel==1)return fPosition1;
  else if (channel==2)return fPosition2;
  else{
    cout <<"Incorrect channel number, set to default channel 1."<<endl;
  }
  return fPosition1;
}

/**********************************************************************/
vector<double> gm2_TEncoder::GetVelocity(int channel) const
{
  if (channel==1)return fVelocity1;
  else if (channel==2)return fVelocity2;
  else{
    cout <<"Incorrect channel number, set to default channel 1."<<endl;
  }
  return fVelocity1;
}

/**********************************************************************/
vector<double> gm2_TEncoder::GetTension(int channel) const
{
  if (channel==1)return fTension1;
  else if (channel==2)return fTension2;
  else{
    cout <<"Incorrect channel number, set to default channel 1."<<endl;
  }
  return fTension1;
}

/**********************************************************************/
vector<double> gm2_TEncoder::GetControlVoltage(int channel) const
{
  if (channel==1)return fControl1;
  else if (channel==2)return fControl2;
  else{
    cout <<"Incorrect channel number, set to default channel 1."<<endl;
  }
  return fControl1;
}

/**********************************************************************/
shared_ptr<TGraph> gm2_TEncoder::GetPositionGraph(int channel,double shift) const
{
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    double Y{0};
    if (channel==1)Y=fPosition1[i];
    else if (channel==2)Y=fPosition2[i];
    else{
      cout <<"Incorrect channel number, set to default channel 1."<<endl;
      Y=fPosition1[i];
    }
    graph_ptr->SetPoint(i,fTimeStamp[i]+shift,Y);
  }
  stringstream ss;
  ss<<channel;
  graph_ptr->SetName("g"+fName+"Position"+ss.str());
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}
/**********************************************************************/
shared_ptr<TGraph> gm2_TEncoder::GetVelocityGraph(int channel,double shift) const
{
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    double Y{0};
    if (channel==1)Y=fVelocity1[i];
    else if (channel==2)Y=fVelocity2[i];
    else{
      cout <<"Incorrect channel number, set to default channel 1."<<endl;
      Y=fVelocity1[i];
    }
    graph_ptr->SetPoint(i,fTimeStamp[i]+shift,Y);
  }
  stringstream ss;
  ss<<channel;
  graph_ptr->SetName("g"+fName+"Velocity"+ss.str());
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}
/**********************************************************************/
shared_ptr<TGraph> gm2_TEncoder::GetTensionGraph(int channel,double shift) const
{
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    double Y{0};
    if (channel==1)Y=fTension1[i];
    else if (channel==2)Y=fTension2[i];
    else{
      cout <<"Incorrect channel number, set to default channel 1."<<endl;
      Y=fTension1[i];
    }
    graph_ptr->SetPoint(i,fTimeStamp[i]+shift,Y);
  }
  stringstream ss;
  ss<<channel;
  graph_ptr->SetName("g"+fName+"Tension"+ss.str());
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}
/**********************************************************************/
shared_ptr<TGraph> gm2_TEncoder::GetControlVoltageGraph(int channel,double shift) const
{
  auto graph_ptr = make_shared<TGraph>(fNPoints);
  for (int i=0;i<fNPoints;i++){
    double Y{0};
    if (channel==1)Y=fControl1[i];
    else if (channel==2)Y=fControl2[i];
    else{
      cout <<"Incorrect channel number, set to default channel 1."<<endl;
      Y=fControl1[i];
    }
    graph_ptr->SetPoint(i,fTimeStamp[i]+shift,Y);
  }
  stringstream ss;
  ss<<channel;
  graph_ptr->SetName("g"+fName+"ControlVoltage"+ss.str());
  graph_ptr->SetTitle(fTitle);
  return graph_ptr;
}
/**********************************************************************/
