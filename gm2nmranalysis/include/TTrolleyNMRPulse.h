/*
 * =====================================================================================
 *       Filename:  TTrolleyNMRPulse.h
 *
 *    Description:  Class for reading and analyzing NMR data from the trolley
 *
 *        Version:  1.0
 *        Created:  06/23/2016
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peter Winter
 *   Organization:  ANL
 This class is inherited from the class TNMRPulse. 
 
 * =====================================================================================
*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <TString.h>
#include "TNMRPulse.h"
#include "TrolleyHeader.h"

using namespace std;

class TTrolleyNMRPulse : public TNMRPulse {
 public:
  TTrolleyNMRPulse(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
  TTrolleyNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
  ~TTrolleyNMRPulse();
  
  // Set methods
  void ReadNextPulse();
  void SetFrequency(float f);
  void SetInputFile(string input);

  // Get methods
  
 private:
  string inputFileName;
  FILE *inputFile;
  float frequency;   // Digitizer frequency
};
