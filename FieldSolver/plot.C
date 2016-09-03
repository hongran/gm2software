{
  ifstream input;
  input.open("output52full.txt",ios::in);
  double phi,B_read,B_cal,B_R,B_Z,B_Phi;
  int Nprb;
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
  TGraph *gB_res = new TGraph();
  gB_res->SetName("gB_res");
  gB_res->SetMarkerStyle(21);
  gB_res->SetMarkerSize(0.3);
  gB_res->SetMarkerColor(kBlue);
  gB_res->SetLineColor(kBlue);
  double B0=0;
  double B0_cal=0;
  for (int i=0;i<196050;i++){
    input>>phi>>Nprb>>B_read>>B_Z>>B_R>>B_Phi>>B_cal;
    B0+=B_read;
    B0_cal+=B_cal;
  }
  input.close();
  B0/=196050.0;
  B0_cal/=196050.0;
  cout <<"Average B = "<<B0<<endl;
  cout <<"Average B_cal = "<<B0_cal<<endl;
  cout <<"Diff B average = "<<(B0_cal-B0)/B0*1000000<<endl;
  int iG=0;
  double RMS1=0;
  double RMS2=0;
  input.open("output52full.txt",ios::in);
  for (int i=0;i<196050;i++){
    input>>phi>>Nprb>>B_read>>B_Z>>B_Phi>>B_R>>B_cal;
    RMS1+=pow(B_read-B0,2.0);
    RMS2+=pow(B_cal-B_read,2.0);
    if (Nprb==1){
      gB_read->SetPoint(iG,phi,B_read);
      gB_cal->SetPoint(iG,phi,B_cal);
      gB_res->SetPoint(iG,phi,(B_read-B_cal)/B0);
      iG++;
    }
  }
  RMS1/=196050.0;
  RMS2/=196050.0;
  RMS1=sqrt(RMS1);
  RMS2=sqrt(RMS2);
  RMS1/=B0;
  RMS2/=B0;
  RMS1*=1000000;
  RMS2*=1000000;
  cout << "Fluctuation = "<<RMS1<<" ppm"<<endl;
  cout << "Fit Residual = "<<RMS2<<" ppm"<<endl;
  TCanvas c1;
  gB_read->Draw("APL");
  gB_read->SetTitle("B field measured vs fit");
  gB_read->GetYaxis()->SetRangeUser(61.78,61.805);
  gB_read->GetYaxis()->SetTitle("Frequency Readout [MHz]");
  gB_read->GetYaxis()->SetTitleOffset(1.5);
  gB_read->GetXaxis()->SetTitle("Phi [rad]");
  gB_cal->Draw("samePL");

  TLegend* Lg = new TLegend(0.7,0.75,0.9,0.9);
  Lg->AddEntry(gB_read,"B measured","pl");
  Lg->AddEntry(gB_cal,"B fit","pl");
  Lg->Draw("same");

  c1.SaveAs("BField.pdf");
  
  TCanvas c2;
  gB_res->Draw("APL");
}
