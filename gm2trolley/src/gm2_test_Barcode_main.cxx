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
#include "gm2_TEncoder.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace std;

int main(int argc,char** argv)
{
  auto ChPos = make_shared<gm2_TRegBarcode>("Pos","Position");
  auto ChDir = make_shared<gm2_TRegBarcode>("Dir","Direction determination");
  auto ChAbs = make_shared<gm2_TAbsBarcode>("Abs","Absolute position");
  auto ChGalil = make_shared<gm2_TEncoder>("Galil","Galil Encoder Information");
  double Ch1,Ch2,Ch3,Ch4,Ch5,Ch6;
  double Tension1,Tension2,V1,V2,P1,P2,O1,O2;
  double Time;
  string Device;

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
      ChGalil->SetPoint(iGalil,Time-Time0G,P1,P2,V1,V2,Tension1,Tension2,O1,O2);
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
  }
  ChPos->Smooth();
  ChPos->SetThreshold(0.01);
  ChPos->FindExtrema();
  ChPos->FindBigGaps();
  ChPos->ConvertToLogic();
  ChDir->SetThreshold(0.01);
  ChDir->ConvertToLogic();
  ChAbs->Smooth();
  ChAbs->SetThreshold(0.05);
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

//  auto NExtremaPos = ChPos->GetNExtrema();
  auto PosExtremaList = ChPos->GetExtremaList();
  /*
  auto gPosCorrelation = make_shared<TGraph>(NExtremaPos);
  for (int i=0;i<NExtremaPos;i++){
    gPosCorrelation->SetPoint(i,PosList1[PosExtremaList[i]],i*4.0);
//    gPosCorrelation->SetPoint(i,PosList1[PosExtremaList[i]],PosExtremaList[i]);
  }
  gPosCorrelation->SetName("PosCorrelation");
  gPosCorrelation->SetTitle("Position correlation to encoder");
*/

  auto gPosRaw = ChPos->GetRawGraph();
  auto gPosLog = ChPos->GetLogicLevelGraph();
  auto gPosAve = ChPos->GetAverageGraph();
  auto gPosCon = ChPos->GetContrastGraph();
  auto gPosExt = ChPos->GetExtremaGraph("VsX");
  auto gPosVel = ChPos->GetVelocityGraph();
  auto hPosHWidth = ChPos->GetLevelWidthHist("High");
  auto hPosLWidth = ChPos->GetLevelWidthHist("Low");

  auto gDirRaw = ChDir->GetRawGraph();
  auto gDirLog = ChDir->GetLogicLevelGraph();
  auto gDirAve = ChDir->GetAverageGraph();
  auto gDirCon = ChDir->GetContrastGraph();
  auto gDirExt = ChDir->GetExtremaGraph("VsX");

  auto gAbsRaw = ChAbs->GetRawGraph();
  auto gAbsExt = ChAbs->GetExtremaGraph("VsX");
  auto gAbsLog = ChAbs->GetLogicLevelGraph();
  auto gAbsWidth = ChAbs->GetAbsWidthGraph();
  auto gAbsSegWidth = ChAbs->GetAbsSegWidthGraph();

  //Determine direction
/*  auto PosLevelList = ChPos->GetLogicLevels();
  auto DirLevelList = ChDir->GetLogicLevels();
  auto NLevels = PosLevelList.size();
  auto NLevelsDir = DirLevelList.size();
  vector<int> AuxDirectionList(NLevels,0);
  for (int i=0;i<NLevels;i++){
    for (int j=i-2;j<=i+2;j++){
      if (j>=0 && j<=NLevelsDir){
	if (DirLevelList[j].REdge>PosLevelList[i].LEdge && DirLevelList[j].REdge<PosLevelList[i].REdge){
	  if (DirLevelList[j].Level==PosLevelList[i].Level){
	    AuxDirectionList[i]=-1;
	  }else{
	    AuxDirectionList[i]=1;
	  }
	  break;
	}
      }
    } 
  }
  vector<int> DirectionList(ChPos->GetNPoints(),0);
  for (int i=0;i<NLevels;i++){
    for (int j=PosLevelList[i].LEdge;j<=PosLevelList[i].REdge;j++){
      DirectionList[j]=AuxDirectionList[i];
    }
  }
  */
  //Determine direction using Galil
  double VThreshold = 5.0;
  auto GalilTimeLine = ChGalil->GetTimeLine();
  auto GalilVelocity = ChGalil->GetVelocity(1);
  auto BarcodeTimeLine = ChPos->GetX();
  auto NVPoints = GalilVelocity.size();
  auto NPosPoints = ChPos->GetNPoints();
  vector<int> DirectionList(NPosPoints,0);
  int j=0;
  for (int i=0;i<NVPoints;i++){
    double end = GalilTimeLine[i];
    while (BarcodeTimeLine[j]<end){
      if (GalilVelocity[i]>VThreshold) DirectionList[j]=1;
      else if (GalilVelocity[i]<-VThreshold) DirectionList[j]=-1;
      else DirectionList[j]=0;
      j++;
    }
    if (j==NPosPoints) break;
  }
  ChPos->SetDirection(DirectionList);
  ChAbs->SetDirection(DirectionList);

  auto gAbsDirection = ChAbs->GetDirectionGraph();

  //Velocity
  auto gPosInterval = ChPos->GetIntervalGraph();
  auto gPosLevelWidth= ChPos->GetLevelWidthGraph();
  auto gDirInterval = ChDir->GetIntervalGraph();
  auto gAbsLevelWidth = ChAbs->GetLevelWidthGraph();

//  auto NGalilPoints = ChGalil->GetNPoints();
  auto gVelocity1 = ChGalil->GetVelocityGraph(1);
  auto gVelocity2 = ChGalil->GetVelocityGraph(2);
  auto gPosition1 = ChGalil->GetPositionGraph(1);
  auto gPosition2 = ChGalil->GetPositionGraph(2);
  auto gTension1= ChGalil->GetTensionGraph(1);
  auto gTension2= ChGalil->GetTensionGraph(2);
  auto gControl1= ChGalil->GetControlVoltageGraph(1);
  auto gControl2= ChGalil->GetControlVoltageGraph(2);
  
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
  
  ChAbs->Decode();

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
  gPosVel->Write();

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
  gAbsDirection->Write();

  gPosInterval->Write();
  gDirInterval->Write();
  gAbsLevelWidth->Write();
  gPosLevelWidth->Write();
//  gPosCorrelation->Write();
  gVelocity1->Write();
  gVelocity2->Write();
  gPosition1->Write();
  gPosition2->Write();
  gTension1->Write();
  gTension2->Write();
  gControl1->Write();
  gControl2->Write();

  hPosHWidth->Write();
  hPosLWidth->Write();

  c1.Write();
  c2.Write();
  c3.Write();
  c4.Write();
  c5.Write();
  output->Close();
  return 0;
}
