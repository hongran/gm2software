#include "TrolleyAnalysisAlgorithm.h"

TrolleyAnalysisAlgorithm::TrolleyAnalysisAlgorithm(){
  Initialize();
}


TrolleyAnalysisAlgorithm::~TrolleyAnalysisAlgorithm(){

}


void TrolleyAnalysisAlgorithm::Initialize(){
  useFitRange = false;
  fitMin = 0.;
  fitMax = 0.;
}


void TrolleyAnalysisAlgorithm::SetFitRange(float min, float max){
  useFitRange = true;
  fitMin = min;
  fitMax = max;
}


void TrolleyAnalysisAlgorithm::GetFitRange(float &min, float &max) const {
  min = fitMin;
  max = fitMax;
}
