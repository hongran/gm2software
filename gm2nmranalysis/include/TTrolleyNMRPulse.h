#ifndef __TTrolleyNMRPulse_h 
#define __TTrolleyNMRPulse_h 
/*
 * =====================================================================================
 *       Filename:  TTrolleyNMRPulse.h
 *
 *    Description:  Class for reading and storing NMR data from the trolley
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
  TTrolleyNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy, vector<double> fyerr);
  ~TTrolleyNMRPulse();
  
  // Set methods
  void SetClockFrequency(float f);  // Set the trolley RF clock frequency in Hz
  void SetClockFrequencyInMHz(float f);  // Set the trolley RF clock frequency providing MHz

  // Get methods
  float GetClockFrequency() const { return ClockFrequency; };  // Get the trolley RF clock frequency in Hz
  float GetDigitizerClockFrequency() const { return DigitizerFrequency; };  // Get the trolley digitizer frequency in Hz

 private:
  float ClockFrequency;   // RF clock frequency
  float DigitizerFrequency;   // Digitizer frequency derived from clock frequency

  // Private methods
  void Initialize();
};

#endif
