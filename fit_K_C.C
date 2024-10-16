#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include <iostream>

void SaveCanvasFit (TFile *file, TString name, TF1 *ffit, TH1D *h1, TH2F *h2, bool logscale, bool ifTH2, std::string legend_text)
{
  gROOT->SetBatch(kTRUE);

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  gPad->SetGrid();
  gStyle->SetOptStat(0);
  if (logscale) c->SetLogy();
  gPad->SetGrid(0, 0);
  h1->SetMaximum(h1->GetMaximum()*1.2);
  h1->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h1->GetYaxis()->SetTitle("Entries");
  ffit->SetLineColor(kRed);
  
  if (ifTH2)
  {
    h2->GetXaxis()->SetTitle("p [GeV/c]");
    h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
    h2->GetZaxis()->SetTitle("Entries");
    c->SetRightMargin(0.15);
    h2->Draw("COLZ");
    c->SetLogz();
  }
  h1->Draw("SAME E0");
  ffit->Draw("SAME");
  
  TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->AddEntry(ffit, legend_text.c_str(), "l");
  legend->Draw();
  
  file->cd();
  c->Write();

  cout << "Canvas " << name << " saved in " << file->GetName() << endl;

  delete c;

  return;
}

TF1 *FitC (TH1D *h1, double p_start, double p_end, double &C, double &Cerr)
{
  TF1 *fit_C = new TF1("fit_C", "pol0", p_start, p_end);
  h1->Fit("fit_C", "RM0Q");

  C = fit_C->GetParameter(0);
  Cerr = fit_C->GetParError(0);

  return fit_C;
}

TF1 *FitK (TH1D *h1, double mass, double p_start, double p_end, double C, double &K, double &Kerr)
{
  TF1 *fit_K = new TF1("fit_K", Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  fit_K->FixParameter(1, C);
  h1->Fit("fit_K", "RM0Q");

  K = fit_K->GetParameter(0);
  Kerr = fit_K->GetParError(0);

  return fit_K;
}

double langaufun(double *x, double *par)
{
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit1peak(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF,TFitResultPtr &fitResultPtr)
{
   int i;
   char FunName[100];
   sprintf(FunName,"Fitfcn_%s",his->GetName());
   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);

   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MPV","Area","GSigma");
   for (i=0; i<4; i++) ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);

   fitResultPtr = his->Fit(FunName,"RB0SLQ");
   ffit->GetParameters(fitparams);
   for (i=0; i<4; i++) fiterrors[i] = ffit->GetParError(i);
  
   ChiSqr[0] = ffit->GetChisquare();
   NDF[0] = ffit->GetNDF();
   return (ffit);
}

TF1 *langaufit2peaks(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF,TFitResultPtr &fitResultPtr)
{
  int i;
  char FunName[100];
  sprintf(FunName,"Fitfcn_%s",his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName, [fitparams](double *x, double *p) {
  return langaufun(x, p) + langaufun(x, &p[4]);
  }, fitrange[0], fitrange[1], 8);

  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width1","MP1","Area1","GSigma1","Width2","MP2","Area2","GSigma2");
  for (i=0; i<8; i++) ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);

  fitResultPtr = his->Fit(FunName,"RB0SQ");
  ffit->GetParameters(fitparams);
  for (i=0; i<8; i++) fiterrors[i] = ffit->GetParError(i);
  
  ChiSqr[0] = ffit->GetChisquare();
  NDF[0] = ffit->GetNDF();

  return (ffit);
}

std::vector <TH1D*> ProjectTH2_eachBinX(const TH2F* h2, bool rebin)
{
  std::vector <TH1D*> h1s;
  for (int i=1; i<=h2->GetNbinsX(); ++i)
  {
    TH1D *h1 = h2->ProjectionY(Form("py_bin%d", i), i, i);
    
    if (rebin) h1 = (TH1D*)h1->Rebin(2, h1->GetName());

    h1s.push_back(h1);
  }

  return h1s;
}

void Compare_Fit(TH1D *h1, double p_start, double p_end, std::string filename)
{
  gROOT->SetBatch(kTRUE);

  TF1 *fit_gaus = new TF1("Gaussian fit", "gaus", p_start, p_end);
  fit_gaus->SetParameters(h1->Integral(), h1->GetMean(), h1->GetStdDev());
  fit_gaus->SetParNames("Normalization","Mean","StdDev");
  h1->Fit("Gaussian fit", "RM0Q");

  TF1 *fit_langaus = new TF1("Langaus fit", langaufun, p_start, p_end, 4);
  fit_langaus->SetParameters(h1->GetStdDev()/2, h1->GetMean(), h1->Integral()*0.08, h1->GetStdDev()*0.1);
  fit_langaus->SetParLimits(0, 0, 0.6);
  fit_langaus->SetParLimits(1, h1->GetMean()-0.5, h1->GetMean()+0.5);
  fit_langaus->SetParLimits(2, h1->Integral()*0.01, h1->Integral()*0.25);
  fit_langaus->SetParLimits(3, 0, h1->GetStdDev());
  fit_langaus->SetParNames("Width","MPV","Area","GSigma");
  h1->Fit("Langaus fit", "RM0Q");

  TCanvas *c = new TCanvas("Comparison of different fits", "Comparison of different fits", 2500, 1773);
  c->cd();
  gPad->SetGrid();
  gStyle->SetOptStat(0);
  h1->SetMaximum(h1->GetMaximum()*1.2);
  h1->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h1->GetYaxis()->SetTitle("Entries");
  h1->GetYaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetRangeUser(0, 7);
  h1->Draw("SAME E0");
  fit_gaus->SetLineColor(kRed);
  fit_gaus->Draw("SAME");
  fit_langaus->SetLineColor(kBlue);
  fit_langaus->Draw("SAME");
  c->SetLogy();
  gPad->SetGrid(0, 0);

  TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->AddEntry(fit_gaus, Form("Gaussian fit : #mu=%3.2f", fit_gaus->GetParameter(1)), "l");
  legend->AddEntry(fit_langaus, Form("Langaus fit : MPV=%3.2f", fit_langaus->GetParameter(1)), "l");
  legend->Draw();

  c->SaveAs(filename.c_str());
}

void fit_K_C()
{
  // Choose file
  bool file_C = false;
  bool file_D = false;
  bool bothfile = true;

  char fileName[100];
  char fileName2[100];
  char SaveFileName[100];
  char SaveFileROOT[100];
  if (file_C)
  {
    sprintf(fileName,"analysis_2024_C.root");
    sprintf(fileName2, "analysis_2024_C.root");
    sprintf(SaveFileName,"Results/K_C_fit_C_file.txt");
    sprintf(SaveFileROOT,"Results/K_C_fit_C_file.root");
  }
  if (file_D)
  {
    sprintf(fileName,"analysis_2024_D.root");
    sprintf(fileName2, "analysis_2024_D.root");
    sprintf(SaveFileName,"Results/K_C_fit_D_file.txt");
    sprintf(SaveFileROOT,"Results/K_C_fit_D_file.root");
  }
  if (bothfile)
  {
    sprintf(fileName, "analysis_2024_C.root");
    sprintf(fileName2, "analysis_2024_D.root");
    sprintf(SaveFileName,"Results/K_C_fit_BOTH_file.txt");
    sprintf(SaveFileROOT,"Results/K_C_fit_BOTH_file.root");
  }

  TFile *file = new TFile(fileName, "READ");
  TFile *file2 = new TFile(fileName2, "READ");
  ofstream Save_K_C(SaveFileName, std::ofstream::out);
  TFile *ROOT_output = new TFile(SaveFileROOT, "RECREATE");


  TH2F *dEdX0stripVsP_lowp = (TH2F*)file->Get("dEdX0stripVsP_lowp");
  TH2F *dEdX0stripVsP_lowp2 = (TH2F*)file2->Get("dEdX0stripVsP_lowp");

  if (bothfile) dEdX0stripVsP_lowp->Add(dEdX0stripVsP_lowp2);

  TH2F *CLONE_Rebin_dEdX0stripVsP_lowp = (TH2F*)dEdX0stripVsP_lowp->Clone("CLONE_Rebin_dEdX0stripVsP_lowp");
  std::vector <TH1D*> dEdx_lowP_C = ProjectTH2_eachBinX(dEdX0stripVsP_lowp, false);
  std::vector <TH1D*> dEdx_lowP_C_CLONE = ProjectTH2_eachBinX(CLONE_Rebin_dEdX0stripVsP_lowp, true);

  TH1D* histo_pion = new TH1D("histo_pion", "histo_pion", dEdX0stripVsP_lowp->GetNbinsX(), 0, dEdX0stripVsP_lowp->GetXaxis()->GetXmax());
  TH1D* histo_kaon = new TH1D("histo_kaon", "histo_kaon", dEdX0stripVsP_lowp->GetNbinsX(), 0, dEdX0stripVsP_lowp->GetXaxis()->GetXmax());
  TH1D* histo_proton = new TH1D("histo_proton", "histo_proton", dEdX0stripVsP_lowp->GetNbinsX(), 0, dEdX0stripVsP_lowp->GetXaxis()->GetXmax());

  int bound_2peak_1peak = 20;
  int bound_2peak = 35;
  int bound_1peak = 60;
  int bound_end = dEdx_lowP_C.size();

  std::vector <double> table_pion;
  std::vector <double> table_kaon;
  std::vector <double> table_proton;
  std::vector <double> table_p;
  for (int i=0; i<bound_2peak_1peak; i++) {table_pion.push_back(0); table_kaon.push_back(0); table_proton.push_back(0);}
  for (int i=0; i<dEdX0stripVsP_lowp->GetNbinsX(); i++) table_p.push_back(dEdX0stripVsP_lowp->GetXaxis()->GetBinCenter(i+1));

  // Iterative fit ---- START
  double startvalues_1peak[4] = {0.2, 9, 2*110000*0.40, 0.2};
  double parlimitslo_1peak[4] = {0.05, 8, 2*110000*0.005, 0};
  double parlimitshi_1peak[4] = {0.6, 10, 2*110000*1.5, 0.5};
  double fitrange_1peak[2] = {8, 14};
  double fitparams_1peak[4];
  double fiterrors_1peak[4];
  double ChiSqr_1peak;
  int NDF_1peak;
  TFitResultPtr fitResultPtr_1peak;

  double startvalues_1peak_pk[8] = {0.1, 3.2, 1000000*0.25, 0.2, 0.1, 4.5, 100000*0.25, 0.5};
  double parlimitslo_1peak_pk[8] = {0, 3, 1000000*0.01, 0, 0, 4.3, 100000*0.01, 0};
  double parlimitshi_1peak_pk[8] = {0.6, 3.3, 1000000, 0.5, 0.6, 5.6, 100000, 0.8};
  double fitrange_1peak_pk[2] = {1, 7};
  double fitparams_1peak_pk[8];
  double fiterrors_1peak_pk[8];
  double ChiSqr_1peak_pk;
  int NDF_1peak_pk;
  TFitResultPtr fitResultPtr_1peak_pk;

  double fitrange_2peak[2] = {1.5, 7};
  double startvalues_2peak[8] = {0.00379035, 3.07077, 1000000, 0.31263, 0.6, 5.03734, 100000, 1.17822e-07};
  double parlimitslo_2peak[8] = {0, 2.8, 1000000*0.1, 0, 0, 5.7, 1000000*0.05, 0};
  double parlimitshi_2peak[8] = {0.6, 3.2, 10000000, 2, 0.6, 6.2, 1000000, 2};
  double fitparams_2peak[8];
  double fiterrors_2peak[8];
  double ChiSqr_2peak;
  int NDF_2peak;
  TFitResultPtr fitResultPtr_2peak;

  // FIT
  for (int i=bound_2peak_1peak; i<bound_2peak; i++)
  {  
    if (i <= bound_2peak_1peak+2) table_proton.push_back(0);
    if (i > bound_2peak_1peak+2)
    {
      TF1 *fit = langaufit1peak(dEdx_lowP_C_CLONE[i], fitrange_1peak, startvalues_1peak, parlimitslo_1peak, parlimitshi_1peak, fitparams_1peak, fiterrors_1peak, &ChiSqr_1peak, &NDF_1peak, fitResultPtr_1peak);
    
      table_proton.push_back(fitparams_1peak[1]);
      histo_proton->SetBinContent(i+1, fitparams_1peak[1]);
      histo_proton->SetBinError(i+1, fiterrors_1peak[1]);

      //SaveCanvasFit(ROOT_output, Form("fit_%d", i), fit, dEdx_lowP_C_CLONE[i], dEdX0stripVsP_lowp, true, false);
      
      if (i > 30 && i < 32) { parlimitslo_1peak[0] = fitparams_1peak[0]*0.5; parlimitshi_1peak[0] = fitparams_1peak[0]*1.5; startvalues_1peak[0] = fitparams_1peak[0];}

      fitrange_1peak[0] -= 0.25; 
      fitrange_1peak[1] -= 0.5;

      for (int i = 1; i < 4; ++i)
      {
        startvalues_1peak[i] = fitparams_1peak[i];
        parlimitslo_1peak[i] = fitparams_1peak[i]*0.3;
        parlimitshi_1peak[i] = fitparams_1peak[i]*1.7;
      }

      startvalues_1peak[2] = fitparams_1peak[2];
      parlimitslo_1peak[2] = fitparams_1peak[2]*0.5;
      parlimitshi_1peak[2] = fitparams_1peak[2]*5;
      
      if (i > 32) { fitrange_1peak[0] = 5.5; fitrange_1peak[1] = 7.4;}
    }
    
    TF1 *fit_pk = langaufit2peaks(dEdx_lowP_C[i], fitrange_1peak_pk, startvalues_1peak_pk, parlimitslo_1peak_pk, parlimitshi_1peak_pk, fitparams_1peak_pk, fiterrors_1peak_pk, &ChiSqr_1peak_pk, &NDF_1peak_pk, fitResultPtr_1peak_pk);
    table_pion.push_back(fitparams_1peak_pk[1]);
    table_kaon.push_back(fitparams_1peak_pk[5]);
    histo_pion->SetBinContent(i, fitparams_1peak_pk[1]);
    histo_pion->SetBinError(i, fiterrors_1peak_pk[1]);
    histo_kaon->SetBinContent(i, fitparams_1peak_pk[5]);
    histo_kaon->SetBinError(i, fiterrors_1peak_pk[5]);

    for (int i = 0; i < 8; ++i)
    {
      startvalues_1peak_pk[i] = fitparams_1peak_pk[i];
      parlimitslo_1peak_pk[i] = fitparams_1peak_pk[i]*0.8;
      parlimitshi_1peak_pk[i] = fitparams_1peak_pk[i]*1.2;
    }
  }


  for (int i=bound_2peak; i<bound_1peak; ++i)
  {
    TF1 *fit = langaufit2peaks(dEdx_lowP_C[i], fitrange_2peak, startvalues_2peak, parlimitslo_2peak, parlimitshi_2peak, fitparams_2peak, fiterrors_2peak, &ChiSqr_2peak, &NDF_2peak, fitResultPtr_2peak);
    table_pion.push_back(fitparams_2peak[1]);
    table_kaon.push_back(fitparams_2peak[1]);
    table_proton.push_back(fitparams_2peak[5]);
    histo_pion->SetBinContent(i+1, fitparams_2peak[1]);
    histo_pion->SetBinError(i+1, fiterrors_2peak[1]);
    histo_kaon->SetBinContent(i+1, fitparams_2peak[1]);
    histo_kaon->SetBinError(i+1, fiterrors_2peak[1]);
    histo_proton->SetBinContent(i+1, fitparams_2peak[5]);
    histo_proton->SetBinError(i+1, fiterrors_2peak[5]);

    if (i < 55)
    {
      if (i == 43) continue;
      for (int i = 0; i < 8; ++i)
      {
        if (i==1 || i==5)
        {
          startvalues_2peak[i] = fitparams_2peak[i];
          parlimitslo_2peak[i] = fitparams_2peak[i]*0.9;
          parlimitshi_2peak[i] = fitparams_2peak[i]*1.1;
        }
      }
    }
    
    if (i < 55) fitrange_2peak[1] -= 0.1;

    //SaveCanvasFit(ROOT_output, Form("fit_%d", i), fit, dEdx_lowP_C[i], dEdX0stripVsP_lowp, true, false);
  }


  for (int i=bound_1peak; i<bound_end; ++i)
  {
    double mean = dEdx_lowP_C[i]->GetMean();
    double stddev = dEdx_lowP_C[i]->GetStdDev();
    double integral = dEdx_lowP_C[i]->Integral();
    double fitrange_1peakend[2] = {mean-2*stddev, mean+2*stddev};
    double startvalues_1peakend[4] = {stddev/2, mean, integral*0.08, stddev*0.1};
    double parlimitslo_1peakend[4] = {0, mean-0.5, integral*0.01, 0};
    double parlimitshi_1peakend[4] = {0.6, mean+0.5, integral*0.25, stddev};
    double fitparams_1peakend[4];
    double fiterrors_1peakend[4];
    double ChiSqr_1peakend;
    int NDF_1peakend;
    TFitResultPtr fitResultPtr_1peakend;

    TF1 *fit = langaufit1peak(dEdx_lowP_C[i], fitrange_1peakend, startvalues_1peakend, parlimitslo_1peakend, parlimitshi_1peakend, fitparams_1peakend, fiterrors_1peakend, &ChiSqr_1peakend, &NDF_1peakend, fitResultPtr_1peakend);
    table_pion.push_back(fitparams_1peakend[1]);
    table_kaon.push_back(fitparams_1peakend[1]);
    table_proton.push_back(fitparams_1peakend[1]);
    histo_pion->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_pion->SetBinError(i+1, fiterrors_1peakend[1]);
    histo_kaon->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_kaon->SetBinError(i+1, fiterrors_1peakend[1]);
    histo_proton->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_proton->SetBinError(i+1, fiterrors_1peakend[1]);

    //SaveCanvasFit(ROOT_output, Form("fit_%d", i), fit, dEdx_lowP_C[i], dEdX0stripVsP_lowp, true, false);
  }


  TGraph *graph_pion = new TGraph(table_p.size(), &table_p[0], &table_pion[0]);
  TGraph *graph_kaon = new TGraph(table_p.size(), &table_p[0], &table_kaon[0]);
  TGraph *graph_proton = new TGraph(table_p.size(), &table_p[0], &table_proton[0]);

  gROOT->SetBatch(kTRUE);
  TCanvas *c_test = new TCanvas("c_test","c_test",2500,1500);
  c_test->cd();
  graph_pion->SetMarkerColor(kRed);
  graph_kaon->SetMarkerColor(kBlue);
  graph_proton->SetMarkerColor(kGreen);
  dEdX0stripVsP_lowp->Draw("COLZ");
  graph_pion->Draw("SAME *");
  graph_kaon->Draw("SAME *");
  graph_proton->Draw("SAME *");

  double C = -1, K = -1;
  double Cerr = -1, Kerr = -1;
  TF1 *fit_C = FitC(histo_pion, 2.5, 5, C, Cerr);
  TF1 *fit_K = FitK(histo_proton, 0.938, 0.56, 1, C, K, Kerr);

  cout << "C = " << C << " pm " << Cerr << "; K = " << K << " pm " << Kerr << " ; Saved in " << SaveFileName << endl;
  SaveCanvasFit(ROOT_output, "C_fit", fit_C, histo_pion, dEdX0stripVsP_lowp, false, true, "Fit to reference data");
  SaveCanvasFit(ROOT_output, "K_fit", fit_K, histo_proton, dEdX0stripVsP_lowp, false, true, "Fit to reference data");

  Save_K_C << K << " " << C << endl;

  Compare_Fit(dEdx_lowP_C[93], 2.6, 4, "Results/CompareFit.root");

  ROOT_output->cd();
  c_test->Write();
  histo_pion->Write();
  histo_kaon->Write();
  histo_proton->Write();
  ROOT_output->Close();
}
