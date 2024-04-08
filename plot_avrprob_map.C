#if defined(__ACLIC__)
#include "DrcPidFast.h"
#include "../prttools/PrtTools.h"
R__LOAD_LIBRARY(../prtdirc/build/libPrt.so)
#else
#include "DrcPidFast.h"
R__LOAD_LIBRARY(../prtdirc/build/libPrt.so)
#endif


void plot_avrprob_map(int pid = 3, int charge = 1, double tphi = -1) {

  PrtTools t;
  t.set_palette(1);
  int pdg = t.pdg(pid);
  gStyle->SetOptStat(0);

  const int nh = 5;
  TString lnames[] = {"e", "#mu", "#pi", "K", "p"};
  TString names[] = {"e", "mu", "pi", "K", "p"};

  const int npphi = 28;
  TH2F *hphimap[9][npphi];
  double bins[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                   0.9, 0.95, 1,   1.25, 1.5, 1.75, 2,   2.25, 2.5, 2.75, 3,   4,    5,   6};
  int pdglist[] = {-11, -13, 211, 321, 2212};

  if (charge == -1) {
    for (int i = 0; i < nh; i++) {
      names[i] += "-";
      lnames[i] += "^{-}";
      pdglist[i] *= -1;
    }
  }

  TH1F *hlh[nh];
  TH2F *heffmap[nh];
  for (int h = 0; h < nh; h++) {
    hlh[h] = new TH1F(names[h], ";L(h_{i}) - L(#pi^{+}); entries / max", 2000, -300, 300);
    heffmap[h] = new TH2F("heffmap" + names[h], names[pid] + "/" + names[h] +
			  " efficiency map ;polar angle [deg];momentum [GeV/c]",
                          135, 25, 160, 51, 0, 10);
  }

  double mean[nh] = {0}, sigma[nh] = {0} , sep[nh] = {0}, eff[nh] = {0};
  DrcPidFast pidf;
  DrcPidInfo info;
  TF1 *ff;
  TString lut = "";
  // TCanvas *cc = new TCanvas("cc", "cc", 1200, 600);

  // "average" track resolution vs momentum
  TF1 *ftrack_res = new TF1("tr", "exp(2.37-0.863*x) + 0.6", 0.2, 10);
  // separation pover vs efficiency assuming perfect Gaussians
  TF1 *fsep_eff = new TF1("eff", "50+19.38*x+1.14*TMath::Sq(x)-1.758*pow(x,3)+0.3586*pow(x,4)-0.03067*pow(x,5)+0.0009852*pow(x,6)",0,10);

 
  // read phi maps
  auto f = new TFile("allphimaps.root");
  
  for (int h = 0; h < nh; h++) {
    for (int pm = 0; pm < npphi; pm++) {
      TString n = Form("phimap_%d_%2.2f", pdglist[h], bins[pm]);
      if(charge == -1) n.ReplaceAll("-2212","2212"); // flip map for antiprotons
      n.ReplaceAll(".00", "");
      n = n.Strip(TString::EStripType::kTrailing, '0');
      hphimap[h][pm] = (TH2F *)f->Get(n);
      
      if(charge == -1 && h == 4){ // flip map for antiprotons
	TH2F * flip =(TH2F*)(hphimap[h][pm]->Clone());	
	int xx = hphimap[h][pm]->GetNbinsX();
	int yy = hphimap[h][pm]->GetNbinsY();
	for (int i = 1; i <= xx; i++) {
	  for (int j = 1; j <= yy; j++) {
	    int b = hphimap[h][pm]->GetBin(i, j);
	    double val = hphimap[h][pm]->GetBinContent(b);
	    flip->SetBinContent(flip->GetBin(i, yy-j), val);
	  }
	}
	hphimap[h][pm] = flip;
      }
      // hphimap[h][pm]->Draw("colz");
      // gPad->Update();
      // gPad->WaitPrimitive();
    }
  }

  for (double m = 0.2; m <= 10; m += 0.2) {
    std::cout << "m " << m << std::endl;
 
    // find momentum bin for phi map
    int pphibin = 0;
    for (int i = 0; i < npphi; i++) {
      if (m <= bins[i]) {
        pphibin = i;
        break;
      }
      pphibin = i;
    }
    if(pphibin == 22) pphibin = 23; // 2.5 GeV/c map missing

    for (int theta = 25; theta <= 160; theta++) {
      // int theta = 30;
      for (double phi = 0; phi <= 30; phi += 0.5) {
        if(tphi > 0) phi = tphi;
        double nphotons = hphimap[pid][pphibin]->GetBinContent(hphimap[pid][pphibin]->FindBin(theta, phi));

        TVector3 mom(0, 0, m);
        mom.RotateX(theta / 180. * TMath::Pi());
        for (int i = 0; i < 2000; i++) {
          double tr = ftrack_res->Eval(m);
          info = pidf.GetInfo(pdg, mom, 0.5, 0);
          for (int h = 0; h < nh; h++) {
            hlh[h]->Fill(info.sigma[h]);
          }
        }

        for (int h = 0; h < nh; h++) {
          hlh[h]->Fit("gaus", "Q0");
          ff = hlh[h]->GetFunction("gaus");
          if (ff) mean[h] = ff->GetParameter(1);
          if (ff) sigma[h] = ff->GetParameter(2);
        }

        double sum = 0;
        for (int h = 0; h < nh; h++) {
          if (sigma[pid] > 0) sep[h] = (fabs(mean[pid] - mean[h])) / (0.5 * (sigma[pid] + sigma[h]));
          eff[h] = 1 - 0.01 * fsep_eff->Eval(sep[h]);
          if (sep[h] > 8) eff[h] = 1 - 0.999999;
        }

        for (int h = 0; h < nh; h++) {
	  if (h == 1) continue; // no muons
          nphotons = hphimap[h][pphibin]->GetBinContent(hphimap[h][pphibin]->FindBin(theta, phi));
          // account for acceptance and threshold
          if (nphotons < 5 || nphotons > 500) eff[h] = 0;
          sum += eff[h];
        }

        // std::cout << "S mom " << m << " " << sep[0] << " " << sep[1] << " " << sep[2] << " "
        //           << sep[3] << " " << sep[4] << std::endl;
        // std::cout << "E mom " << m << " " << eff[0] << " " << eff[1] << " " << eff[2] << " "
        //           << eff[3] << " " << eff[4] << std::endl;

        // normalization
        if (sum > 0) {
          for (int h = 0; h < nh; h++) eff[h] /= sum;
        }

        // std::cout << "N mom " << m << " " << eff[0] << " " << eff[1] << " " << eff[2] << " "
        //           << eff[3] << " " << eff[4] << std::endl;
        // std::cout << "-------------------- " << std::endl;

        // t.normalize_to(hlh, nh, 1);
        // for (int h = 0; h < nh; h++) {
        // 	// hlh[h]->GetXaxis()->SetRangeUser(-100, 100);
        // 	hlh[h]->SetTitle(Form("likelihoods for #pi @ p = %2.1f GeV/c",m));
        //   hlh[h]->SetLineColor(t.color(h));
        //   hlh[h]->Draw((h == 0) ? "hist" : "hist same");
        // }
        // cc->Update();
        // cc->WaitPrimitive();

        for (int h = 0; h < nh; h++) {
          heffmap[h]->Fill(theta, m, eff[h]);
          hlh[h]->Reset();
        }



        lut += Form("%d %d %2.2f %2.2f %2.2f", t.pdg(pid), charge, m, (double)theta, phi);
        for (int h = 0; h < nh; h++) {
	  if (h == 1) continue;
          lut += Form(" %3.4f", eff[h]);
        }
	lut += "\n";
        // std::cout << "lut " << lut << std::endl;
	
	if (tphi > 0) break;	
      }
    }
  }

  t.write_string("hpdirc_"+names[pid] +".lut",lut);
  gStyle->SetOptStat(0);
  for (int h = 0; h < nh; h++) {
    if (h == 1) continue; // no muons
    t.add_canvas("map_eff_" + names[pid] + "_" + names[h] + Form("_phi_%2.1f", tphi), 800, 500);
    heffmap[h]->Draw("colz");
    heffmap[h]->GetYaxis()->SetRangeUser(0, 10);
    heffmap[h]->SetMaximum(1);
    heffmap[h]->SetMinimum(0);
  }  

  t.save_canvas("data/plot_avrprob_map_n", 0);
}
