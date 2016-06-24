#ifndef __TrolleyAnalysisManager_h
#define __TrolleyAnalysisManager_h

/*
 * =====================================================================================
 *       Filename:  TrolleyAnalysisManager.h
 *
 *    Description:  Implementation of AnalysisManager class for the trolley NMR data
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

#include "AnalysisManager.h"
#include "TTrolleyNMRPulse.h"

using namespace std;

class TrolleyAnalysisManager : public AnalysisManager<TTrolleyNMRPulse>{
 public:
  TrolleyAnalysisManager();
  ~TrolleyAnalysisManager();

  void Initialize();

  // Set methods
  virtual void SetPulse(shared_ptr<TTrolleyNMRPulse> p);
  void SetInputFile(string input);
  virtual void SetVerbosityLevel(int verbose){ verbosityLevel = verbose; };

  // Get methods
  virtual int ReadNextPulse();
  virtual shared_ptr<TTrolleyNMRPulse> GetPulse() const { return nmrPulse; }

 protected:
  string inputFileName;
  FILE *inputFile;
};

#endif 
