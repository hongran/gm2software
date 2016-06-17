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
  double Ch1,Ch2,Ch3;
  double Tension1,Tension2,V1,V2,P1,P2,O1,O2;

  //Read from file
  string FileName = string{argv[1]}+".txt";
  ifstream filein;
  filein.open(FileName.c_str(),ios::in);
  filein.ignore(200,'\n');
  int i=0;
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
    filein>>Tension1>>Tension2>>V1>>V2>>P1>>P2>>O1>>O2>>Ch1>>Ch2>>Ch3;
    ChPos->SetPoint(i,i,Ch1);
    ChAbs->SetPoint(i,i,Ch2);
    ChDir->SetPoint(i,i,Ch3);
    i++;
  }
  ChPos->SetThreshold(0.02);
  ChPos->FindExtrema();
  ChPos->ConvertToLogic();
  ChDir->SetThreshold(0.02);
  ChDir->ConvertToLogic();
  ChAbs->SetThreshold(0.02);
  ChAbs->ConvertToLogic();
  cout << ChPos->GetNLowLevels()<<" " << ChPos->GetNHighLevels()<<endl;
  cout << ChDir->GetNLowLevels()<<" " << ChDir->GetNHighLevels()<<endl;
  cout << ChAbs->GetNLowLevels()<<" " << ChAbs->GetNHighLevels()<<endl;

  auto gPosRaw = ChPos->GetRawGraph();
  auto gPosLog = ChPos->GetLogicLevelGraph();
  auto gPosAve = ChPos->GetAverageGraph();
  auto gPosCon = ChPos->GetContrastGraph();
  auto gPosExt = ChPos->GetExtremaGraph("VsX");

  auto gDirRaw = ChDir->GetRawGraph();
  auto gDirLog = ChDir->GetLogicLevelGraph();
  auto gDirAve = ChDir->GetAverageGraph();
  auto gDirCon = ChDir->GetContrastGraph();
  auto gDirExt = ChDir->GetExtremaGraph("VsX");

  auto gAbsRaw = ChAbs->GetRawGraph();
//  auto gAbsLog = ChAbs->GetLogicLevelGraph();
//  auto gAbsExt = ChAbs->GetExtremaGraph("VsX");

  TCanvas c1;
  gPosRaw->Draw("APL");
  gDirRaw->Draw("samePL");
  gAbsRaw->Draw("samePL");
  gPosRaw->SetLineColor(kBlack);
  gDirRaw->SetLineColor(kBlue);
  gAbsRaw->SetLineColor(kRed);


  TCanvas c2;
  gPosRaw->Draw("APL");
  gPosExt->Draw("sameP");
  gPosExt->SetMarkerStyle(20);
  gPosExt->SetMarkerSize(1);
  gPosExt->SetMarkerColor(kRed);

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

  c1.Write();
  c2.Write();
  output->Close();
  return 0;
}
