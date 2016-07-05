#ifndef __TrolleyAnalysisAlgorithm_h
#define __TrolleyAnalysisAlgorithm_h

/*
 * =====================================================================================
 *       Filename:  AnalysisManager.h
 *
 *    Description:  Abstract analysis algorithm class for analyzing trolley NMR pulses
 *
 *        Version:  1.0
 *        Created:  06/25/2016
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peter Winter
 *   Organization:  ANL
 * =====================================================================================
 */

class TrolleyAnalysisAlgorithm {
 public:
  TrolleyAnalysisAlgorithm();
  ~TrolleyAnalysisAlgorithm();

  // Set methods
  void SetFitRange(float min, float max);

  // Get methods
  void GetFitRange(float &min, float &max) const;

 protected:
  bool useFitRange;
  double fitMin;
  double fitMax;

  void Initialize();
};

#endif
