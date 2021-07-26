#include "DrcPidFast.cxx"

void example(int pdg = 211, double mom = 6) {

  DrcPidFast pid;

  TVector3 momentum(0, 0, mom);
  momentum.RotateX(20 / 180. * TMath::Pi());

  TH1F *hPi = new TH1F("hPi", ";#Delta Ln;entries", 200, -10, 10);
  TH1F *hK = new TH1F("hK", ";#Delta Ln;entries", 200, -10, 10);

  DrcPidInfo info;
  for (int i = 0; i < 1000; i++) {
    info = pid.GetInfo(pdg, momentum, 0.5);
    hPi->Fill(info.sigma[2]);
    hK->Fill(info.sigma[3]);
  }

  auto r1 = hPi->Fit("gaus", "S");
  double m1 = r1->Parameter(1);
  double s1 = r1->Parameter(2);

  auto r2 = hK->Fit("gaus", "S");
  double m2 = r2->Parameter(1);
  double s2 = r2->Parameter(2);

  double sep = (fabs(m1 - m2)) / (0.5 * (s1 + s2));

  std::cout<<"separation "<<sep<<std::endl;  
  
  hPi->Draw();
  hPi->SetLineColor(kBlue);
  hK->SetLineColor(kRed);
  hK->Draw("same");
}
