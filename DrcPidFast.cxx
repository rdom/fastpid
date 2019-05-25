#include "DrcPidFast.h"

DrcPidFast::DrcPidFast(){

  fMass[0]=0.000511;
  fMass[1]=0.105658;
  fMass[2]=0.139570;
  fMass[3]=0.49368;
  fMass[4]=0.938272;
  
  // read Cherenkov track resolution map
  TFile* file = TFile::Open("ctr_map.root");
  fTrrMap = new TH2F();
  file->GetObject("htrr", fTrrMap);  
}

DrcPidInfo DrcPidFast::GetInfo(int pdg,TVector3 mom, double err){
  
  const int max = 5;
  DrcPidInfo info;
  int pid = get_pid(pdg);
  double m = mom.Mag();
  double theta = mom.Theta()*180/TMath::Pi();
  
  // set default values
  for(int i=0; i<max; i++){
    info.probability[i]=0.25;
    info.sigma[i]=100;
  }
  
  // check range
  if(m>10) m = 10;
  if(theta<25 || theta>140) return info;
      
  int bin = fTrrMap->FindBin(theta,m);
  double trr = fTrrMap->GetBinContent(bin);
  double ctr = sqrt(trr*trr+err*err)*0.001;
  
  // 1.46907 - fused silica
  double true_cangle = acos(sqrt(m*m + fMass[pid]*fMass[pid])/m/1.46907);  
  true_cangle += fRand.Gaus(0,ctr);

  // return default values if momentum below Cherenkov threshold (true_cangle is NaN)
  if(true_cangle != true_cangle) return info;
  
  double cangle,sum=0,fsum=0;
  double delta[max]={0}, probability[max]={0};

  for(int i=0; i<max; i++){
    cangle = acos(sqrt(m*m + fMass[i]*fMass[i])/m/1.46907);
    if(cangle != cangle) continue;
    delta[i] = fabs(cangle-true_cangle);
    sum += delta[i];
    info.sigma[i]=(cangle-true_cangle)/ctr;
  }

  // normalization
  for(int i=0; i<max; i++){
    if(delta[i]>0) info.probability[i] = sum/delta[i];
    fsum += info.probability[i];
  }
  for(int i=0; i<max; i++) info.probability[i] /= fsum;

  return info;
}

int DrcPidFast::get_pid(int pdg){
  int pid=0;
  if(pdg==11)   pid=0; //e
  if(pdg==13)   pid=1; //mu
  if(pdg==211)  pid=2; //pi
  if(pdg==321)  pid=3; //K
  if(pdg==2212) pid=4; //p
  return pid;
}
