{
  ifstream input;
  input.open("output52full.txt",ios::in);
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
  TGraph *gB_Phi = new TGraph();
  gB_Phi->SetName("gB_Phi");
  gB_Phi->SetMarkerStyle(21);
  gB_Phi->SetMarkerSize(0.3);
  gB_Phi->SetMarkerColor(kRed);
  gB_Phi->SetLineColor(kRed);
  TGraph *gB_R = new TGraph();
  gB_R->SetName("gB_R");
  gB_R->SetMarkerStyle(21);
  gB_R->SetMarkerSize(0.3);
  gB_R->SetMarkerColor(kRed);
  gB_R->SetLineColor(kRed);
  int Nprb;
  double BR,BPhi,BZ;
  int iGraph=0;
  for (int i=0;i<196050;i++){
    input>>phi>>Nprb>>B_read>>BZ>>BR>>BPhi>>B_cal;
    if (Nprb==20){
      BPhi/=61.74;
      BR/=61.74;
      BPhi*=1000000;
      BR*=1000000;
      B_read-=61.74;
      B_cal-=61.74;
      B_read/=61.74;
      B_cal/=61.74;
      B_read*=1000000;
      B_cal*=1000000;

      gB_read->SetPoint(iGraph,phi,B_read);
      gB_cal->SetPoint(iGraph,phi,B_cal);
      gB_Phi->SetPoint(iGraph,phi,BPhi);
      gB_R->SetPoint(iGraph,phi,BR);
    }
    if (Nprb==20)iGraph++;
  }
  input.close();

  TFile fout("Plots.root","recreate");
  gB_read->Write();
  gB_cal->Write();
  gB_Phi->Write();
  gB_R->Write();
  fout.Close();
  TCanvas c1;
  gB_read->Draw("APL");
  gB_read->SetTitle("B field measured vs fit");
//  gB_read->GetYaxis()->SetRangeUser(-200,200);
  gB_read->GetYaxis()->SetTitle("ppm");
  gB_read->GetYaxis()->SetTitleOffset(1.5);
  gB_read->GetXaxis()->SetTitle("Phi [rad]");
  gB_cal->Draw("samePL");

  TLegend* Lg = new TLegend(0.7,0.75,0.9,0.9);
  Lg->AddEntry(gB_read,"B measured","pl");
  Lg->AddEntry(gB_cal,"B fit","pl");
  Lg->Draw("same");

  c1.SaveAs("BField.pdf");

  TCanvas c2;
  gB_Phi->Draw("APL");
  gB_Phi->SetTitle("B_Phi field fit");
  gB_Phi->GetYaxis()->SetTitle("ppm");
  gB_Phi->GetYaxis()->SetTitleOffset(1.5);
  gB_Phi->GetXaxis()->SetTitle("Phi [rad]");

  c2.SaveAs("BFieldPhi.pdf");

  TCanvas c3;
  gB_R->Draw("APL");
  gB_R->SetTitle("B_R field fit");
  gB_R->GetYaxis()->SetTitle("ppm");
  gB_R->GetYaxis()->SetTitleOffset(1.5);
  gB_R->GetXaxis()->SetTitle("Phi [rad]");

  c3.SaveAs("BFieldR.pdf");
}
