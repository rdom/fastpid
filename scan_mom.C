#include "prttools.C"
#include "DrcPidFast.cxx"

void scan_mom(int pdg = 211, double mom = 6) {

  DrcPidFast pid;


  TH1F *hPi = new TH1F("hPi", ";#Delta Ln;entries", 200, -150, 150);
  TH1F *hK = new TH1F("hK", ";#Delta Ln;entries", 200, -150, 150);

  TH1F *hPi10 = new TH1F("hPi10", ";#Delta Ln;entries", 200, -150, 150);
  TH1F *hK10 = new TH1F("hK10", ";#Delta Ln;entries", 200, -150, 150);
  
  TGraphAsymmErrors *gmom[2];

  DrcPidInfo info;

  for (int s = 0; s < 2; s++) {
    gmom[s] = new TGraphAsymmErrors();
    int id = 0;
    for (double m = 1.5; m <= 6.1; m += 0.1) {
      TVector3 momentum(0, 0, m);
      momentum.RotateX(30 / 180. * TMath::Pi());

      for (int i = 0; i < 2000; i++) {
        info = pid.GetInfo(pdg, momentum, 0.5, (s == 0) ? 17 : 10);
        hPi->Fill(info.sigma[2]);
        hK->Fill(info.sigma[3]);
      }

      auto r1 = hPi->Fit("gaus", "SQ0");
      double m1 = r1->Parameter(1);
      double s1 = r1->Parameter(2);
      double dm1 = r1->ParError(1);
      double ds1 = r1->ParError(2);

      auto r2 = hK->Fit("gaus", "SQ0");
      double m2 = r2->Parameter(1);
      double s2 = r2->Parameter(2);
      double dm2 = r2->ParError(1);
      double ds2 = r2->ParError(2);

      double sep = (fabs(m1 - m2)) / (0.5 * (s1 + s2));
      double e1, e2, e3, e4;
      e1 = 2 / (s1 + s2) * dm1;
      e2 = -2 / (s1 + s2) * dm2;
      e3 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds1;
      e4 = -2 * (m1 - m2) / ((s1 + s2) * (s1 + s2)) * ds2;
      double esep = sqrt(e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4) + 0.1;

      gmom[s]->SetPoint(id, m, sep);
      gmom[s]->SetPointEYhigh(id, esep);
      gmom[s]->SetPointEYlow(id, esep);
      id++;

      std::cout << "separation " << sep << " +/- " << esep << std::endl;
      hPi->Reset();
      hK->Reset();
    }
  }

  prt_canvasAdd("mom_pik", 800, 400);

  TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.87);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  for (int s = 0; s < 2; s++) {

    gmom[s]->SetMarkerStyle(20);
    gmom[s]->SetMarkerSize(0.75);
    gmom[s]->GetXaxis()->SetTitle("momentum [GeV/c]");
    gmom[s]->GetYaxis()->SetTitle("separation [s.d.]");

    gmom[s]->GetXaxis()->SetLimits(1.3, 6.5);
    gmom[s]->GetXaxis()->SetRangeUser(1.3, 6.5);
    gmom[s]->GetYaxis()->SetRangeUser(0, 20);

    // gmom[s]->GetXaxis()->SetLimits(0.5, 2.1);
    // gmom[s]->GetXaxis()->SetRangeUser(0.5, 2.1);
    // gmom[s]->GetYaxis()->SetRangeUser(0, 10);

    if (s == 0) {
      gmom[s]->SetLineColor(kGray+2);
      gmom[s]->SetMarkerColor(kBlack);      
      gmom[s]->Draw("APL");
      leg->AddEntry(gmom[s], "17 mm", "lp");
    } else {
      gmom[s]->SetLineColor(kRed);
      gmom[s]->SetMarkerColor(kRed+1);
      gmom[s]->Draw("PL");
      leg->AddEntry(gmom[s], "10 mm", "lp");
    }

  }
  leg->Draw();


  prt_savepath = "data/scan_mom";
  prt_canvasSave(2);
}
