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
 *       Compiler:  g++
 *
 *         Author:  Peter Winter
 *   Organization:  ANL
 * =====================================================================================
 */

#include "AnalysisManager.h"
#include "TTrolleyNMRPulse.h"
#include "TrolleyAnalysisAlgorithm.h"

using namespace std;

class TrolleyAnalysisManager : public AnalysisManager<TTrolleyNMRPulse, 
  TrolleyAnalysisAlgorithm>{
 public:
  TrolleyAnalysisManager();
  ~TrolleyAnalysisManager();

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

  virtual void Initialize();
};

#endif 
