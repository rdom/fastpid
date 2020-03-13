//#include "../prttools/prttools.C"
#include "prttools.C"
#include "DrcPidFast.h"

void plot_map(int pdg=211){

  prt_setRootPalette(1);
  TH2F *hsep = new TH2F("hsep","#pi/K separation power [s.d.] ;polar angle [deg];momentum [GeV/c]",129,24.5,153.5,50,0.1,10.1);
  TH1F *hPi = new TH1F("hPi","hPi",500,-100,100);
  TH1F *hK = new TH1F("hK","hK",500,-100,100);
  double m1=0,m2=0,s1=0,s2=0,sep=0;
  DrcPidFast pid;    
  DrcPidInfo info;
  TF1 *ff;
  //TCanvas *cc = new TCanvas("cc","cc",800,500);
  for(double m=0.6; m<=10; m=m+0.2){
    if(m>3){
      delete hK;
      delete hPi;
      hK = new TH1F("hK","hK",200,-20,20);
      hPi = new TH1F("hPi","hPi",200,-20,20);
    }
    if(m<1){
      delete hK;
      delete hPi;
      hK = new TH1F("hK","hK",1000,-500,500);
      hPi = new TH1F("hPi","hPi",1000,-500,500);
    }
    std::cout<<"m "<<m<<std::endl;
    
    for(int t=25; t<=153; t++){
      TVector3 mom(0,0,m);
      mom.RotateX(t/180.*TMath::Pi());
      for(int i=0; i<1000; i++){
	info = pid.GetInfo(pdg,mom,0.5);
	hPi->Fill(info.sigma[2]);
	hK->Fill(info.sigma[3]);    
      }
      
      {
	hK->Fit("gaus","Q0");
	ff = hK->GetFunction("gaus");
	if(ff) m1=ff->GetParameter(1);
	if(ff) s1=ff->GetParameter(2);
	hPi->Fit("gaus","Q0");
	ff = hPi->GetFunction("gaus");
	if(ff) m2=ff->GetParameter(1);
	if(ff) s2=ff->GetParameter(2);

	if(s1+s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
      }
      
      // hPi->Draw();
      // hPi->SetLineColor(kBlue);
      // hK->SetLineColor(kRed);
      // hK->Draw("same");
      // cc->Update();
      // cc->WaitPrimitive();
      
      hsep->Fill(t,m,sep);

      hPi->Reset();
      hK->Reset();
    }
  }
  
  // prt_canvasAdd("diff",800,500);
  // hPi->Draw();
  // hPi->SetLineColor(kBlue);
  // hK->SetLineColor(kRed);
  // hK->Draw("same");

  gStyle->SetOptStat(0);
  prt_canvasAdd("map_sep",800,500);
  gPad->SetLogz();
  hsep->Draw("colz");
  hsep->GetYaxis()->SetRangeUser(0.6,10);
  hsep->SetMaximum(10);
  hsep->SetMinimum(1);


  // draw ctr map
  prt_canvasAdd("map_ctr",800,500);
  pid.GetTrrMap()->SetTitle("Cherenkov track resolution [mrad]");
  pid.GetTrrMap()->Draw("colz");
  
  prt_savepath ="data/plot_map";
  prt_canvasSave(2);  
}
