/*
 * =====================================================================================
 *
 *       Filename:  gm2_test_Barcode_main.cxx
 *
 *    Description:  Stand-alone program that tests the GM2_TBarcode class
 *
 *        Version:  1.0
 *        Created:  03/24/2016 10:34:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ran Hong 
 *   Organization:  ANL
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include "gm2_TBarcode.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace std;

int main(int argc,char** argv)
{
  auto ChPos = make_shared<gm2_TRegBarcode>("Pos","Position");
  auto ChDir = make_shared<gm2_TRegBarcode>("Dir","Direction determination");
  auto ChAbs = make_shared<gm2_TAbsBarcode>("Abs","Absolute position");
  double Ch1,Ch2,Ch3,Ch4,Ch5,Ch6;
  double Tension1,Tension2,V1,V2,P1,P2,O1,O2;
  double Time;
  string Device;

  vector<double> PosList1,PosList2;
  vector<double> Velocity1,Velocity2;

  //Read from file
  string FileName = string{argv[1]}+".txt";
  ifstream filein;
  filein.open(FileName.c_str(),ios::in);
  char Buffer[100000];
  int iGalil{0};
  int iBarcode{0};
  double Time0G{0};
  double Time0B{0};
  while(!filein.eof()){
//    filein>>T;
//    filein.ignore(10,',');
//    filein>>Ch1;
//    filein.ignore(10,',');
//    filein>>Ch2;
//    filein.ignore(10,',');
//    filein>>Ch3;
//    filein.ignore(10,',');
//    filein>>Ch4;
//    filein.ignore(10,'\n');
    filein.getline(Buffer,sizeof(Buffer));
    stringstream iss{string{Buffer}};
    iss>>Device;
    if (Device.compare("Galil")==0){
      iss>>Time>>Tension1>>Tension2>>P1>>P2>>V1>>V2>>O1>>O2;
      if (iGalil==0)Time0G=Time;
      iGalil++;
    }
    if (Device.compare("Barcode")==0){
      iss>>Time>>Ch1>>Ch2>>Ch3>>Ch4>>Ch5>>Ch6;
      if (iBarcode==0)Time0B=Time;
      ChPos->SetPoint(iBarcode,Time-Time0B,Ch1);
      ChAbs->SetPoint(iBarcode,Time-Time0B,Ch2);
      ChDir->SetPoint(iBarcode,Time-Time0B,Ch3);
/*      ChPos->SetPoint(iBarcode,Time-Time0B,Ch1);
      ChAbs->SetPoint(iBarcode,Time-Time0B,Ch2);
      ChDir->SetPoint(iBarcode,Time-Time0B,Ch3);
      */
      iBarcode++;
    }
    /*
    ChPos->SetPoint(i,P1,Ch1);
    ChAbs->SetPoint(i,P1,Ch2);
    ChDir->SetPoint(i,P1,Ch3);
    */
    PosList1.push_back(P1);
    PosList2.push_back(P2);
    Velocity1.push_back(V1);
    Velocity2.push_back(V2);
  }
  ChPos->SetThreshold(0.01);
  ChPos->FindExtrema();
  ChPos->ConvertToLogic();
  ChDir->SetThreshold(0.01);
  ChDir->ConvertToLogic();
  ChAbs->SetThreshold(0.06);
  ChAbs->FindExtrema();
  ChAbs->ConvertToLogic();
  ChAbs->ChopSegments(*ChPos);
  cout << ChPos->GetNLowLevels()<<" " << ChPos->GetNHighLevels()<<endl;
  cout << ChDir->GetNLowLevels()<<" " << ChDir->GetNHighLevels()<<endl;
  cout << ChAbs->GetNLowLevels()<<" " << ChAbs->GetNHighLevels()<<endl;
  //Set Logic Level Scale
  ChPos->SetLogicLevelScale(0.4);
  ChDir->SetLogicLevelScale(0.4);
  ChAbs->SetLogicLevelScale(0.4);

  auto NExtremaPos = ChPos->GetNExtrema();
  auto PosExtremaList = ChPos->GetExtremaList();
  auto gPosCorrelation = make_shared<TGraph>(NExtremaPos);
  for (int i=0;i<NExtremaPos;i++){
    gPosCorrelation->SetPoint(i,PosList1[PosExtremaList[i]],i*4.0);
//    gPosCorrelation->SetPoint(i,PosList1[PosExtremaList[i]],PosExtremaList[i]);
  }
  gPosCorrelation->SetName("PosCorrelation");
  gPosCorrelation->SetTitle("Position correlation to encoder");


  auto gPosRaw = ChPos->GetRawGraph(250);
  auto gPosLog = ChPos->GetLogicLevelGraph();
  auto gPosAve = ChPos->GetAverageGraph();
  auto gPosCon = ChPos->GetContrastGraph();
  auto gPosExt = ChPos->GetExtremaGraph("VsX");

  auto gDirRaw = ChDir->GetRawGraph();
  auto gDirLog = ChDir->GetLogicLevelGraph();
  auto gDirAve = ChDir->GetAverageGraph();
  auto gDirCon = ChDir->GetContrastGraph();
  auto gDirExt = ChDir->GetExtremaGraph("VsX");

  auto gAbsRaw = ChAbs->GetRawGraph(250);
  auto gAbsExt = ChAbs->GetExtremaGraph("VsX");
  auto gAbsLog = ChAbs->GetLogicLevelGraph();
  auto gAbsWidth = ChAbs->GetAbsWidthGraph();
  auto gAbsSegWidth = ChAbs->GetAbsSegWidthGraph();

  //Velocity
  auto gPosInterval = ChPos->GetIntervalGraph();
  auto gPosLevelWidth= ChPos->GetLevelWidthGraph();
  auto gDirInterval = ChDir->GetIntervalGraph();
  auto gAbsLevelWidth = ChAbs->GetLevelWidthGraph();

  auto NPoints = Velocity1.size();
  auto gVelocity1 = make_shared<TGraph>(NPoints);
  auto gVelocity2 = make_shared<TGraph>(NPoints);
  auto gPosition1 = make_shared<TGraph>(NPoints);
  auto gPosition2 = make_shared<TGraph>(NPoints);
  for (int j=0;j<NPoints;j++){
//    gVelocity1->SetPoint(j,j,Velocity1[j]);
//    gVelocity2->SetPoint(j,j,Velocity2[j]);
    gVelocity1->SetPoint(j,j,Velocity1[j]);
    gVelocity2->SetPoint(j,j,Velocity2[j]);
    gPosition1->SetPoint(j,j,PosList1[j]);
    gPosition2->SetPoint(j,j,PosList2[j]);
  }
  
  //Calculate position resolution
  auto NInterval = gPosInterval->GetN();
  double average{0};
  double average2{0};
  int k{0};
  for (int i=NInterval/10;i<NInterval*4/10;i++){
    double x,y;
    gPosInterval->GetPoint(i,x,y);
    average+=y;
    average2+=y*y;
    k++;
  }
  average/=k;
  average2/=k;
  cout << "average interval = "<<average<<endl;
  cout << "RMS = "<<sqrt(average2-average*average)<<endl;


  TCanvas c1;
  gPosRaw->Draw("APL");
  gDirRaw->Draw("samePL");
  gAbsRaw->Draw("samePL");
  gPosRaw->SetLineColor(kBlack);
  gDirRaw->SetLineColor(kBlue);
  gAbsRaw->SetLineColor(kRed);


  TCanvas c2("c2");
  gPosRaw->Draw("APL");
  gPosExt->Draw("sameP");
  gPosLog->Draw("sameL");
  gPosExt->SetMarkerStyle(20);
  gPosExt->SetMarkerSize(1);
  gPosExt->SetMarkerColor(kGreen);

  TCanvas c3("c3");
  gDirRaw->Draw("APL");
  gDirExt->Draw("sameP");
  gDirExt->SetMarkerStyle(20);
  gDirExt->SetMarkerSize(1);
  gDirExt->SetMarkerColor(kRed);

  TCanvas c4("c4");
  gAbsRaw->Draw("APL");
  gAbsRaw->SetLineColor(kGreen);
  gAbsExt->Draw("sameP");
  gAbsExt->SetMarkerStyle(20);
  gAbsExt->SetMarkerSize(1);
  gAbsExt->SetMarkerColor(kRed);
  gAbsLog->Draw("sameL");
  gAbsLog->SetLineColor(kBlue);

  TCanvas c5("c5");
  gAbsLog->Draw("APL");
  gAbsLog->SetLineColor(kBlue);
  gAbsLog->SetLineWidth(4);
  gPosLog->Draw("sameL");
  gPosLog->SetLineWidth(2);
  gPosLog->SetLineColor(kRed);

  TFile* output = new TFile((string{argv[1]}+".root").c_str(),"recreate");
  ChPos->Write();
  ChDir->Write();
  ChAbs->Write();
  gPosRaw->Write();
  gPosLog->Write();
  gPosAve->Write();
  gPosCon->Write();
  gPosExt->Write();

  gDirRaw->Write();
  gDirLog->Write();
  gDirAve->Write();
  gDirCon->Write();
  gDirExt->Write();

  gAbsRaw->Write();
  gAbsLog->Write();
  gAbsExt->Write();
  gAbsWidth->Write();
  gAbsSegWidth->Write();

  gPosInterval->Write();
  gDirInterval->Write();
  gAbsLevelWidth->Write();
  gPosLevelWidth->Write();
  gPosCorrelation->Write();
  gVelocity1->Write();
  gVelocity2->Write();
  gPosition1->Write();
  gPosition2->Write();

  c1.Write();
  c2.Write();
  c3.Write();
  c4.Write();
  c5.Write();
  output->Close();
  return 0;
}
