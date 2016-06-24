#ifndef __AnalysisManager_h
#define __AnalysisManager_h

/*
 * =====================================================================================
 *       Filename:  AnalysisManager.h
 *
 *    Description:  AbstractAnalysis manager template class for analyzing NMR pulses
 *
 *        Version:  1.0
 *        Created:  06/24/2016
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peter Winter
 *   Organization:  ANL
 * =====================================================================================
 */

#include <iostream>
#include <stdio.h>
#include <TString.h>
#include "TNamed.h"
#include "TNMRPulse.h"

using namespace std;

template <class Pulse> class AnalysisManager {
 public:
  virtual void Initialize() = 0;

  // Set methods
  virtual void SetPulse(shared_ptr<Pulse> p) = 0;

  // Get methods
  virtual int ReadNextPulse() = 0;
  virtual shared_ptr<Pulse> GetPulse() const = 0;
  
 protected:
  shared_ptr<Pulse> nmrPulse;
};

#endif 
