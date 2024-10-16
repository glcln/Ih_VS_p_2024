#include <iostream>
#include <TFile.h>
#include <TH2F.h>

std::vector <TH1D*> ProjectTH2_eachBinX(const TH2F* h2)
{
  std::vector <TH1D*> h1s;
  for (int i=1; i<=h2->GetNbinsX(); ++i)
  {
    TH1D *h1 = h2->ProjectionY(Form("py_bin%d", i), i, i);
    h1s.push_back(h1);
  }

  return h1s;
}

void Return_K_C_fit(double &K, double &C, std::string filename)
{
  std::ifstream infile(filename);
  infile >> K >> C;
  infile.close();

  return;
}

TF1 *Fit_MassParametrization(double mass, double p_start, double p_end, double C, double K)
{
  TF1 *fit = new TF1("fit", Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  fit->SetParameter(0, K);
  fit->FixParameter(1, C);

  return fit;
}

void Display_TH2_Fit(TH2F *h2, TF1* fit_original_pion, TF1* fit_original_proton, TF1 *fit_pion, TF1 *fit_kaon, TF1 *fit_proton, TF1 *fit_deuteron, std::string filename_PDF, std::string filename_C, std::string filename_ROOT, std::string filename_PNG)
{
  gROOT->SetBatch(kTRUE);

  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("Track momentum [GeV/c]");
  h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  h2->GetZaxis()->SetTitle("Number of tracks");
  h2->SetTitle(" ");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetRangeUser(0, 14);

  TCanvas *c = new TCanvas("c","c",2500,1773);
  c->SetRightMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  h2->Draw("colz");
  c->SetLogz();
  fit_pion->Draw("same");
  fit_kaon->Draw("same");
  fit_proton->Draw("same");
  fit_deuteron->Draw("same");
  fit_original_pion->SetLineColor(kBlack);
  fit_original_proton->SetLineColor(kBlack);
  fit_original_pion->Draw("same");
  fit_original_proton->Draw("same");

  TLatex * tex = new TLatex(0.85,0.915,"2024 (13.6 TeV)");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.10,0.915,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.08);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.22,0.96,"Preliminary");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(52);
  tex->SetTextSize(0.0608);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.12,0.3,"#pi");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.115,0.86,"K");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.15,0.86,"p");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.215,0.86,"D");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  TLegend *leg = new TLegend(0.57,0.7,0.82,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(fit_original_pion,"Fit to reference data","l");
  leg->AddEntry(fit_kaon,"Extrapolation","l");
  leg->Draw();


  c->SaveAs(filename_PDF.c_str());
  c->SaveAs(filename_C.c_str());
  c->SaveAs(filename_ROOT.c_str());
  c->SaveAs(filename_PNG.c_str());

  delete c;

  return;
}

void Display_TH2_NoFit(TH2F *h2, std::string filename_PDF, std::string filename_C, std::string filename_ROOT, std::string filename_PNG, bool charge)
{
  gROOT->SetBatch(kTRUE);

  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("Track momentum [GeV/c]");
  if (charge) h2->GetXaxis()->SetTitle("charge sign x Track momentum [GeV/c]");
  h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  h2->GetZaxis()->SetTitle("Number of tracks");
  h2->SetTitle(" ");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetRangeUser(0, 14);

  TCanvas *c = new TCanvas("c","c",2500,1773);
  c->SetRightMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  h2->Draw("colz");
  c->SetLogz();

  TLatex * tex = new TLatex(0.85,0.915,"2024 (13.6 TeV)");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.10,0.915,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.08);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.22,0.96,"Preliminary");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(52);
  tex->SetTextSize(0.0608);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();

  c->SaveAs(filename_PDF.c_str());
  c->SaveAs(filename_C.c_str());
  c->SaveAs(filename_ROOT.c_str());
  c->SaveAs(filename_PNG.c_str());

  delete c;

  return;
}


void dEdx_VS_p_Display()
{
    // Initialisation
    TFile *file_C = new TFile("ROOT_histograms/analysis_2024_C.root", "READ");
    TFile *file_D = new TFile("ROOT_histograms/analysis_2024_D.root", "READ");
    
    // Récupérer le TH2F
    TH2F *dEdX0stripVsP_lowp_C = (TH2F*)file_C->Get("dEdX0stripVsP_lowp");
    TH2F *dEdX0stripVsP_charge_C = (TH2F*)file_C->Get("dEdX0stripVsP_charge");
    TH2F *dEdX0stripVsP_lowp_D = (TH2F*)file_D->Get("dEdX0stripVsP_lowp");
    TH2F *dEdX0stripVsP_charge_D = (TH2F*)file_D->Get("dEdX0stripVsP_charge");

    TH2F *dEdX0stripVsP_lowp_BOTH = (TH2F*)dEdX0stripVsP_lowp_C->Clone("dEdX0stripVsP_lowp_BOTH");
    dEdX0stripVsP_lowp_BOTH->Add(dEdX0stripVsP_lowp_D);
    
    // Récupérer le Fit
    double K_Cperiod = -1, K_Dperiod = -1;
    double C_Cperiod = -1, C_Dperiod = -1;
    double K_BOTH = -1, C_BOTH = -1;
    Return_K_C_fit(K_Cperiod, C_Cperiod, "Results/K_C_fit_C_file.txt");
    Return_K_C_fit(K_Dperiod, C_Dperiod, "Results/K_C_fit_D_file.txt");
    Return_K_C_fit(K_BOTH, C_BOTH, "Results/K_C_fit_BOTH_file.txt");

    TF1* Fit_pion_original_Cperiod = Fit_MassParametrization(0.13957, 2.5, 5, C_Cperiod, K_Cperiod);
    TF1* Fit_proton_original_Cperiod = Fit_MassParametrization(0.93827, 0.575, 1, C_Cperiod, K_Cperiod);
    TF1* Fit_pion_original_Dperiod = Fit_MassParametrization(0.13957, 2.5, 5, C_Dperiod, K_Dperiod);
    TF1* Fit_proton_original_Dperiod = Fit_MassParametrization(0.93827, 0.57, 1, C_Dperiod, K_Dperiod);
    TF1* Fit_pion_original_BOTH = Fit_MassParametrization(0.13957, 2.5, 5, C_BOTH, K_BOTH);
    TF1* Fit_proton_original_BOTH = Fit_MassParametrization(0.93827, 0.575, 1, C_BOTH, K_BOTH);
    
    TF1* Fit_pion_Cperiod = Fit_MassParametrization(0.13957, 0.2, 5, C_Cperiod, K_Cperiod);
    TF1* Fit_kaon_Cperiod = Fit_MassParametrization(0.49368, 0.2, 5, C_Cperiod, K_Cperiod);
    TF1* Fit_proton_Cperiod = Fit_MassParametrization(0.93827, 0.2, 5, C_Cperiod, K_Cperiod);
    TF1* Fit_deuteron_Cperiod = Fit_MassParametrization(1.87561, 0.2, 5, C_Cperiod, K_Cperiod);
    TF1* Fit_pion_Dperiod = Fit_MassParametrization(0.13957, 0.2, 5, C_Dperiod, K_Dperiod);
    TF1* Fit_kaon_Dperiod = Fit_MassParametrization(0.49368, 0.2, 5, C_Dperiod, K_Dperiod);
    TF1* Fit_proton_Dperiod = Fit_MassParametrization(0.93827, 0.2, 5, C_Dperiod, K_Dperiod);
    TF1* Fit_deuteron_Dperiod = Fit_MassParametrization(1.87561, 0.2, 5, C_Dperiod, K_Dperiod);
    TF1* Fit_pion_BOTH = Fit_MassParametrization(0.13957, 0.19, 5, C_BOTH, K_BOTH);
    TF1* Fit_kaon_BOTH = Fit_MassParametrization(0.49368, 0.2, 5, C_BOTH, K_BOTH);
    TF1* Fit_proton_BOTH = Fit_MassParametrization(0.93827, 0.19, 5, C_BOTH, K_BOTH);
    TF1* Fit_deuteron_BOTH = Fit_MassParametrization(1.87561, 0.90, 5, C_BOTH, K_BOTH);

    // Affichage
    Display_TH2_NoFit(dEdX0stripVsP_charge_C, "Results/dEdX0stripVsP_charge_C.pdf", "dResults/EdX0stripVsP_charge_C.C", "Results/dEdX0stripVsP_charge_C.root", "Results/dEdX0stripVsP_charge_C.png", true);
    Display_TH2_NoFit(dEdX0stripVsP_lowp_C, "Results/dEdX0stripVsP_lowp_C.pdf", "Results/dEdX0stripVsP_lowp_C.C", "Results/dEdX0stripVsP_lowp_C.root", "Results/dEdX0stripVsP_lowp_C.png", false);
    Display_TH2_Fit(dEdX0stripVsP_lowp_C, Fit_pion_original_Cperiod, Fit_proton_original_Cperiod, Fit_pion_Cperiod, Fit_kaon_Cperiod, Fit_proton_Cperiod, Fit_deuteron_Cperiod, "Results/dEdX0stripVsP_lowp_FIT_C.pdf", "Results/dEdX0stripVsP_lowp_FIT_C.C", "Results/dEdX0stripVsP_lowp_FIT_C.root", "Results/dEdX0stripVsP_lowp_FIT_C.png");
   
    Display_TH2_NoFit(dEdX0stripVsP_charge_D, "Results/dEdX0stripVsP_charge_D.pdf", "Results/dEdX0stripVsP_charge_D.C", "Results/dEdX0stripVsP_charge_D.root", "Results/dEdX0stripVsP_charge_D.png", true);
    Display_TH2_NoFit(dEdX0stripVsP_lowp_D, "Results/dEdX0stripVsP_lowp_D.pdf", "Results/dEdX0stripVsP_lowp_D.C", "Results/dEdX0stripVsP_lowp_C.root", "Results/dEdX0stripVsP_lowp_C.png", false);
    Display_TH2_Fit(dEdX0stripVsP_lowp_D, Fit_pion_original_Dperiod, Fit_proton_original_Dperiod, Fit_pion_Dperiod, Fit_kaon_Dperiod, Fit_proton_Dperiod, Fit_deuteron_Dperiod, "Results/dEdX0stripVsP_lowp_FIT_D.pdf", "Results/dEdX0stripVsP_lowp_FIT_D.C", "Results/dEdX0stripVsP_lowp_FIT_C.root", "Results/dEdX0stripVsP_lowp_FIT_C.png");

    Display_TH2_NoFit(dEdX0stripVsP_charge_D, "Results/dEdX0stripVsP_charge_BOTH.pdf", "Results/dEdX0stripVsP_charge_BOTH.C", "Results/dEdX0stripVsP_charge_BOTH.root", "Results/dEdX0stripVsP_charge_BOTH.png", true);
    Display_TH2_NoFit(dEdX0stripVsP_lowp_BOTH, "Results/dEdX0stripVsP_lowp_BOTH.pdf", "Results/dEdX0stripVsP_lowp_BOTH.C", "Results/dEdX0stripVsP_lowp_BOTH.root", "Results/dEdX0stripVsP_lowp_BOTH.png", false);
    Display_TH2_Fit(dEdX0stripVsP_lowp_BOTH, Fit_pion_original_BOTH, Fit_proton_original_BOTH, Fit_pion_BOTH, Fit_kaon_BOTH, Fit_proton_BOTH, Fit_deuteron_BOTH, "Results/dEdX0stripVsP_lowp_FIT_BOTH.pdf", "Results/dEdX0stripVsP_lowp_FIT_BOTH.C", "Results/dEdX0stripVsP_lowp_FIT_BOTH.root", "Results/dEdX0stripVsP_lowp_FIT_BOTH.png");


    // Eta distribution and mass reconstruction
    TFile *file_eta = new TFile("ROOT_histograms/EtaDistribution.root", "READ");
    TFile *save_eta_mass = new TFile("ROOT_histograms/save_eta_mass.root", "RECREATE");

    TH1D *eta_distribution = (TH1D*)file_eta->Get("eta_distribution");
    TH2D *mass_reconstruction = (TH2D*)file_eta->Get("mass_reconstruction");

    mass_reconstruction->GetXaxis()->SetTitle("Track momentum [GeV/c]");
    mass_reconstruction->GetYaxis()->SetTitle("Reconstructed mass [GeV/c^{2}]");
    mass_reconstruction->SetTitle(" ");
    eta_distribution->GetXaxis()->SetTitle("#eta");

    save_eta_mass->cd();
    eta_distribution->Write();
    mass_reconstruction->Write();
    save_eta_mass->Close();


}
