#include "DrcPidFast.h"

DrcPidFast::DrcPidFast() {

  fMass[0] = 0.000511;
  fMass[1] = 0.105658;
  fMass[2] = 0.139570;
  fMass[3] = 0.49368;
  fMass[4] = 0.938272;

  // read Cherenkov track resolution map
  ReadMap("ctr_map_p1_0.95.root");

  // multiple scattering for 17 mm thick radiator at 30 deg
  fMs_mom = new TF1("", "expo(0)+expo(2)+expo(4)");
  fMs_mom->SetParameters(9.39815e-01, -1.48243e-01, 4.69733e+00, -4.33960e+00, 2.19745e+00,
                        -9.68617e-01);

  fMs_thickness_17 = new TF1("", "pol1");
  fMs_thickness_17->SetParameters(4.2, 0.034); // 17 mm raidator bar

  fMs_thickness_10 = new TF1("", "pol1");
  fMs_thickness_10->SetParameters(3.6, 0.021); // 10 mm radiator bar
  fMs_thickness_max = fMs_thickness_17->Eval(70);
}

void DrcPidFast::ReadMap(TString name) {
  TFile *file = TFile::Open(name);
  fTrrMap = new TH2F();
  file->GetObject("htrr", fTrrMap);
}

DrcPidInfo DrcPidFast::GetInfo(int pdg, TVector3 mom, double track_err, int thickness) {
  double p = mom.Mag();
  double theta = mom.Theta() * TMath::RadToDeg();
  return GetInfo(pdg, p, theta, track_err, thickness);
}

DrcPidInfo DrcPidFast::GetInfo(int pdg, double p, double theta, double track_err, int thickness) {

  // double theta = 2.0*atan(exp(-eta))*TMath::RadToDeg();

  const int max = 5;
  DrcPidInfo info;
  int pid = get_pid(pdg);

  // set default values
  for (int i = 0; i < max; i++) {
    info.probability[i] = 0.25;
    info.sigma[i] = 100;
  }
  info.cangle = 0;
  info.cctr = 0;

  // check range
  if (theta < 19.99 || theta > 160.01){
    std::cout<<"theta out of [20,160] deg range: "<<theta<<std::endl;    
  }

  double ms_mom_err = fMs_mom->Eval(p); // vector deviation after radiator 

  double alpha = (theta < 90) ? 90 - theta : theta - 90;

  double ms_thick_frac;
  if (thickness == 10) ms_thick_frac = fMs_thickness_10->Eval(alpha) / fMs_thickness_max;
  else ms_thick_frac = fMs_thickness_17->Eval(alpha) / fMs_thickness_max;
  
  // 0.28 for averaging direction vector over the radiator thickness
  double ms_err = 0.22 * ms_mom_err * ms_thick_frac;
  // double ms_err = 0.45 * ms_mom_err * ms_thick_frac;

  // ctr map is for theta = [25,153] and p = [0,10] GeV/c
  if (theta < 25) theta = 25;
  if (theta > 153) theta = 153;
  if (p > 10) p = 10;
  
  int bin = fTrrMap->FindBin(theta, p);  
  double ctr = fTrrMap->GetBinContent(bin);             // Cherenkov track resolution [mrad]
  double cctr = sqrt(ctr * ctr + track_err * track_err + ms_err * ms_err) *
                0.001; // combined Cherenkov track resolution[rad]

  // 1.46907 - fused silica
  double true_cangle = acos(sqrt(p * p + fMass[pid] * fMass[pid]) / p / 1.46907);
  true_cangle += gRandom->Gaus(0, cctr);

  // return default values if momentum below Cherenkov threshold (true_cangle is NaN)
  if (true_cangle != true_cangle) return info;

  double cangle, sum = 0, fsum = 0;
  double delta[max] = {0}, probability[max] = {0};

  for (int i = 0; i < max; i++) {
    cangle = acos(sqrt(p * p + fMass[i] * fMass[i]) / p / 1.46907);
    if (cangle != cangle) continue;
    delta[i] = fabs(cangle - true_cangle);
    sum += delta[i];
    info.sigma[i] = (cangle - true_cangle) / cctr;
    if (i == pid) info.cangle = cangle;
  }
  // normalization
  for (int i = 0; i < max; i++) {
    if (delta[i] > 0) info.probability[i] = sum / delta[i];
    fsum += info.probability[i];
  }
  for (int i = 0; i < max; i++) info.probability[i] /= fsum;
  info.cctr = cctr;

  return info;
}

int DrcPidFast::get_pid(int pdg) {
  int pid = 0;
  if (pdg == 11) pid = 0;   // e
  if (pdg == 13) pid = 1;   // mu
  if (pdg == 211) pid = 2;  // pi
  if (pdg == 321) pid = 3;  // K
  if (pdg == 2212) pid = 4; // p
  return pid;
}
