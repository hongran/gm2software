#ifndef __AnalysisManager_h
#define __AnalysisManager_h

/*
 * =====================================================================================
 *       Filename:  AnalysisManager.h
 *
 *    Description:  Abstract analysis manager template class for analyzing NMR pulses
 *
 *        Version:  1.0
 *        Created:  06/24/2016
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peter Winter
 *   Organization:  ANL
 * =====================================================================================
 */

#include <memory>

using namespace std;

template <class Pulse, class AnalysisAlgorithm> class AnalysisManager {
 public:
  // Set methods
  virtual void SetPulse(shared_ptr<Pulse> p) = 0;
  virtual void SetVerbosityLevel(int verbose) = 0;

  // Get methods
  virtual int ReadNextPulse() = 0;
  virtual shared_ptr<Pulse> GetPulse() const = 0;
  
 protected:
  shared_ptr<Pulse> nmrPulse;

  AnalysisAlgorithm algorithm;

  int verbosityLevel;

  virtual void Initialize() = 0;
};

#endif 
