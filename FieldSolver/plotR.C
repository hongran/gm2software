{
  ifstream input;
  input.open("output52R.txt",ios::in);
  double phi,B_read,B_cal;
  TGraph *gB_read = new TGraph();
  gB_read->SetName("gB_read");
  gB_read->SetMarkerStyle(20);
  gB_read->SetMarkerSize(0.3);
  gB_read->SetMarkerColor(kBlue);
  gB_read->SetLineColor(kBlue);
  TGraph *gB_cal = new TGraph();
  gB_cal->SetName("gB_cal");
  gB_cal->SetMarkerStyle(21);
  gB_cal->SetMarkerSize(0.3);
  gB_cal->SetMarkerColor(kRed);
  gB_cal->SetLineColor(kRed);
  for (int i=0;i<7842;i++){
    input>>phi>>B_read>>B_cal;
    gB_read->SetPoint(i,phi,B_read);
    gB_cal->SetPoint(i,phi,B_cal);
  }
  TCanvas c1;
  gB_cal->Draw("APL");
  gB_cal->SetTitle("B field measured vs fit");
//  gB_cal->GetYaxis()->SetRangeUser(61.78,61.805);
  gB_cal->GetYaxis()->SetTitle("Frequency Readout [MHz]");
  gB_cal->GetYaxis()->SetTitleOffset(1.5);
  gB_cal->GetXaxis()->SetTitle("Phi [rad]");

  TLegend* Lg = new TLegend(0.7,0.75,0.9,0.9);
  Lg->AddEntry(gB_cal,"B fit","pl");
  Lg->Draw("same");

  c1.SaveAs("BFieldR.pdf");
}
