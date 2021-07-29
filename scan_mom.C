#include "prttools.C"
#include "DrcPidFast.cxx"

void scan_mom(int pdg = 211, double mom = 6) {

  DrcPidFast pid(0);


  TH1F *hPi = new TH1F("hPi", ";#Delta Ln;entries", 200, -50, 50);
  TH1F *hK = new TH1F("hK", ";#Delta Ln;entries", 200, -50, 50);
  TGraphAsymmErrors *gmom = new TGraphAsymmErrors();

  DrcPidInfo info;
  
  int id = 0;
  for (double m = 0.5; m <= 3.1; m += 0.1) {
    TVector3 momentum(0, 0, m);
    momentum.RotateX(30 / 180. * TMath::Pi());

    for (int i = 0; i < 1000; i++) {
      info = pid.GetInfo(pdg, momentum, 0.5);
      hPi->Fill(info.sigma[2]);
      hK->Fill(info.sigma[0]);
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
    
    gmom->SetPoint(id, m, sep);
    gmom->SetPointEYhigh(id, esep);
    gmom->SetPointEYlow(id, esep);
    id++;
    
    std::cout << "separation " << sep << " +/- " << esep << std::endl;
    hPi->Reset();
    hK->Reset();
  }

  prt_canvasAdd("mom", 800, 500);
  gmom->SetMarkerStyle(20);
  gmom->SetMarkerSize(0.75);
  gmom->GetXaxis()->SetTitle("momentum [GeV/c]");
  gmom->GetYaxis()->SetTitle("separation [s.d.]");
  gmom->GetXaxis()->SetLimits(0.3, 3.2);
  gmom->GetXaxis()->SetRangeUser(0.3, 3.2);
  gmom->GetYaxis()->SetRangeUser(0, 14);
  
  gmom->Draw("APL");
  
  prt_savepath = "data/scan_mom";
  prt_canvasSave(2);
}
