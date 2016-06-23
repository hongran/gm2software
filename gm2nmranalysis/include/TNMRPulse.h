/*
 * =====================================================================================
 *       Filename:  TNMRPulse.h
 *
 *    Description:  Class for reading and analyzing NMR data
 *
 *        Version:  1.0
 *        Created:  06/23/2016
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peter Winter
 *   Organization:  ANL
 This class is inherited from 
 the ROOT TNamed class. 
 
 * =====================================================================================
*/

#include <iostream>
#include <stdio.h>
#include <TString.h>
#include "TNamed.h"

using namespace std;

class TNMRPulse : public TNamed {
 public:
  TNMRPulse(const TString& Name = TString{"EmptyName"}, const TString& Title = TString{"EmptyTitle"});
  TNMRPulse(const TString& Name, const TString& Title, vector<double> fx, vector<double> fy);
  ~TNMRPulse();
  
  // Set methods
  void SetPoint(const int i, const double x, const double y);

  // Get methods
  int GetNPoints() const { return fNPoints; }
  vector<double> GetX() const { return fX; }
  vector<double> GetY() const { return fY; }

 private:
  vector<double> fX;
  vector<double> fY;
  int fNPoints;
};
