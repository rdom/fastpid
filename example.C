#include "DrcPidFast.h"

void example(int pdg=211,double mom=6){

  DrcPidFast pid;
  
  TVector3 momentum(0,0,mom);
  momentum.RotateX(25/180.*TMath::Pi());

  TH1F *hPi = new TH1F("hPi","",200,-10,10);
  TH1F *hK = new TH1F("hK","",200,-10,10);

  DrcPidInfo info;
  for(int i=0; i<1000; i++){
    info = pid.GetInfo(pdg,momentum,0);
    hPi->Fill(info.sigma[2]);
    hK->Fill(info.sigma[3]);
  }

  hPi->Draw();
  hPi->SetLineColor(kBlue);
  hK->SetLineColor(kRed);
  hK->Draw("same");

}
