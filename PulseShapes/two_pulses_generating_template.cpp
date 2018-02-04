#include <iostream>
#include "Math/Interpolator.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <TMath.h>
#include <string>
#include <TF1.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;
const int nsamples = 150;
const int nboards = 5;
const int nchannels_per_board = 5;
const int nchannels = nboards*nchannels_per_board;
const int impinged_channel = 12; //Change me each new crystal!
const int impinged_channel2 = 2;

const float beam_energy = 100.0; //If want to do two different energy crystals can do a beam_energy2


//Need to duplicate these for the second pulse and that should fix the problem.
float temp_amp[nchannels];
float temp_time[nchannels];
float EA_X;
float EA_Y;

float WF_time[nsamples*nchannels];
float WF_val[nsamples*nchannels];
float X[2];
float Y[2];

//Need to duplicate these for the second pulse and that should fix the problem.
float temp_amp2[nchannels];
float temp_time2[nchannels];
float EA_X2;
float EA_Y2;

float WF_time2[nsamples*nchannels];
float WF_val2[nsamples*nchannels];
float X2[2];
float Y2[2];

void two_pulses_generating_template(){
  TCanvas *canvas = new TCanvas("generating template","generating template");
  TFile *file = new TFile("/afs/cern.ch/user/k/kyee/work/pulsetest/vfe_55_WC7416.root");
  //TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_55_C3_100.root");

  TFile *file2 = new TFile("/afs/cern.ch/user/k/kyee/work/pulsetest/vfe_55_WC7290.root");

  TTree *mtree = (TTree*) file->Get("h4");
  mtree->SetBranchAddress("X",X);
  mtree->SetBranchAddress("Y",Y);
  mtree->SetBranchAddress("WF_time",WF_time);
  mtree->SetBranchAddress("WF_val",WF_val);

  TTree *mtree2 = (TTree*) file2->Get("h4");
  mtree2->SetBranchAddress("X",X2);
  mtree2->SetBranchAddress("Y",Y2);
  mtree2->SetBranchAddress("WF_time",WF_time2);
  mtree2->SetBranchAddress("WF_val",WF_val2);

 TFile *template_recos = new TFile(" /afs/cern.ch/user/k/kyee/work/pulsetest/template_recos_D3_100.root ");
//  TFile *template_recos = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/template_recos_new_runs_C3_100.root");

 TFile *template_recos2 = new TFile(" /afs/cern.ch/user/k/kyee/work/pulsetest/template_recos_C3_100.root ");

  TTree *template_tree = (TTree*) template_recos->Get("template_tree");
  template_tree->SetBranchAddress("temp_amp",temp_amp);
  template_tree->SetBranchAddress("temp_time",temp_time);
  template_tree->SetBranchAddress("EA_X",&EA_X);
  template_tree->SetBranchAddress("EA_Y",&EA_Y);
  TF1 *g = new TF1("g","gaus",0,1000);
  TH1F *waveform = new TH1F("waveform - D3 - 100 GeV - 6*6 mm^2 at center","waveform - D3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);
  TProfile *mean_waveform = new TProfile("mean waveform - D3 - 100 GeV - 6*6 mm^2 at center","mean waveform - D3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);

  TTree *template_tree2 = (TTree*) template_recos2->Get("template_tree");
  template_tree2->SetBranchAddress("temp_amp",temp_amp2);
  template_tree2->SetBranchAddress("temp_time",temp_time2);
  template_tree2->SetBranchAddress("EA_X",&EA_X2);
  template_tree2->SetBranchAddress("EA_Y",&EA_Y2);
  TF1 *g2 = new TF1("g2","gaus",0,1000);
  TH1F *waveform2 = new TH1F("waveform - C3 - 100 GeV - 6*6 mm^2 at center","waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);
  TProfile *mean_waveform2 = new TProfile("mean waveform - C3 - 100 GeV - 6*6 mm^2 at center","mean waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);

// First Crystal

  const Long64_t nentries = mtree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0,ientry;
  float dummy_amp[nentries];
  float dummy_time[nentries];
  bool dummy_constraint[nentries];
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = template_tree->LoadTree(jentry);
    nb = template_tree->GetEntry(jentry);
    dummy_amp[jentry] = temp_amp[impinged_channel];
    dummy_time[jentry] = temp_time[impinged_channel];
    if (TMath::Abs(EA_X)<2&&TMath::Abs(EA_Y)<2&&temp_amp[2]>35&&temp_amp[13]>35&&temp_amp[11]>35&&temp_amp[7]>35&&dummy_amp[jentry] > 3000.0*beam_energy/100.0) dummy_constraint[jentry] = true; //NOTE: Please also change me with every crystal! //Note: I'm also worried the EA_X and Y cuts are too much. Also check each reconstructed run to make sure I did, indeed, recenter X and Y each time, otherwise this won't work.
    else dummy_constraint[jentry] = false;
  }
  
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);   
    nbytes += nb;
    if (TMath::Abs(X[0]+6)<10&&TMath::Abs(Y[0]-7)<10){                                                 //To be tuned
    if (dummy_constraint[jentry] == true){                                                 //To be tuned
      waveform->Reset();
      for (int i = 0;i < nsamples;i++) waveform->SetBinContent(i+1,WF_val[i+impinged_channel*nsamples]);
//      for (int i = 0;i < nsamples;i++) mean_waveform->Fill(WF_time[i+impinged_channel*nsamples]-dummy_time[jentry],WF_val[i+impinged_channel*nsamples]/dummy_amp[jentry]); //comment this line and uncomment the rest of lines in the loop if it's the first time you are generating templates and they haven't been used for fits (basically there's no temp_reco file)
      if (waveform->GetMaximum() < 2000.0*beam_energy/100.0) continue;                          //vice versa if otherwise
      g->SetParameters(1.2*waveform->GetMaximum(),6.25*(waveform->GetMaximumBin()-1),6);
      g->SetParLimits(0,1500.0*beam_energy/100.0,10000);
      g->SetParLimits(1,240,320);
      waveform->Fit("g","Q");
      waveform->Fit("g","Q","",-18+g->GetParameter(1),18+g->GetParameter(1));
      if (g->GetParameter(1) < 260 || g->GetParameter(1) > 280) continue;
      if (g->GetParameter(0) < 2300*beam_energy/100.0) continue;
      for (int i = 0;i < nsamples;i++) mean_waveform->Fill(270+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1),WF_val[i+impinged_channel*nsamples]); // */
    }
    }
  }
  
  //Next Crystal

 const Long64_t nentries2 = mtree2->GetEntriesFast();
  Long64_t nbytes2 = 0, nb2 = 0,ientry2;
  float dummy_amp2[nentries2];
  float dummy_time2[nentries2];
  bool dummy_constraint2[nentries2];
  for (Long64_t jentry=0; jentry<nentries2;jentry++){
    ientry2 = template_tree2->LoadTree(jentry);
    nb2 = template_tree2->GetEntry(jentry);
    dummy_amp2[jentry] = temp_amp2[impinged_channel2];
    dummy_time2[jentry] = temp_time2[impinged_channel2];
    if (TMath::Abs(EA_X2)<2&&TMath::Abs(EA_Y2)<2&&temp_amp2[1]>35&&temp_amp2[3]>35&&temp_amp2[12]>35&&temp_amp2[17]>35&&dummy_amp2[jentry] > 3000.0*beam_energy/100.0) dummy_constraint2[jentry] = true; //NOTE: Please also change me with every crystal! //Note: I'm also worried the EA_X and Y cuts are too much. Also check each reconstructed run to make sure I did, indeed, recenter X and Y each time, otherwise this won't work.
    else dummy_constraint2[jentry] = false;
  }
  
  for (Long64_t jentry=0; jentry<nentries2;jentry++){
    Long64_t ientry2 = mtree2->LoadTree(jentry);
    nb2 = mtree2->GetEntry(jentry);   
    nbytes2 += nb2;
    if (TMath::Abs(X2[0]+6)<10&&TMath::Abs(Y2[0]-7)<10){                                                 //To be tuned
    if (dummy_constraint2[jentry] == true){                                                 //To be tuned
      waveform2->Reset();
      for (int i = 0;i < nsamples;i++) waveform2->SetBinContent(i+1,WF_val2[i+impinged_channel2*nsamples]);
//      for (int i = 0;i < nsamples;i++) mean_waveform->Fill(WF_time[i+impinged_channel*nsamples]-dummy_time[jentry],WF_val2[i+impinged_channel*nsamples]/dummy_amp2[jentry]); //comment this line and uncomment the rest of lines in the loop if it's the first time you are generating templates and they haven't been used for fits (basically there's no temp_reco file)
      if (waveform2->GetMaximum() < 2000.0*beam_energy/100.0) continue;                          //vice versa if otherwise
      g2->SetParameters(1.2*waveform2->GetMaximum(),6.25*(waveform2->GetMaximumBin()-1),6);
      g2->SetParLimits(0,1500.0*beam_energy/100.0,10000);
      g2->SetParLimits(1,240,320);
      waveform2->Fit("g2","Q");
      waveform2->Fit("g2","Q","",-18+g2->GetParameter(1),18+g2->GetParameter(1));
      if (g2->GetParameter(1) < 260 || g2->GetParameter(1) > 280) continue;
      if (g2->GetParameter(0) < 2300*beam_energy/100.0) continue;
      for (int i = 0;i < nsamples;i++) mean_waveform2->Fill(270+WF_time2[i+impinged_channel2*nsamples]-g2->GetParameter(1),WF_val2[i+impinged_channel2*nsamples]); // */
    }
    }
  }

  //End Crystal


  TAxis *xaxis = mean_waveform->GetXaxis();
  const int nTempBins = 100*nsamples;
  TH1F *interpolated_mean_waveform = new TH1F("interpolated mean waveform - D3 - 100 GeV - 6*6 mm^2 at center","interpolated mean waveform - D3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nTempBins,-0.125,937.375);
  TAxis *xaxis2 = interpolated_mean_waveform->GetXaxis();
  double x[nsamples+2],y[nsamples+2];
  x[0] = xaxis->GetBinCenter(1) - 6.25;
  y[0] = 0;
  for (int i = 0;i < nsamples;i++){
    //x[i] = 6.25 * i;
    x[i+1] = xaxis->GetBinCenter(i+1);
    y[i+1] = mean_waveform->GetBinContent(i+1);
  }
  x[nsamples+1] = x[nsamples]+6.25;
  y[nsamples+1] = 0;
  //ROOT::Math::Interpolator inter(nsamples+1, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator inter(nsamples+2, ROOT::Math::Interpolation::kAKIMA);
  inter.SetData(nsamples+2,x,y);
  for (int i = 0;i < nTempBins;i++) interpolated_mean_waveform->SetBinContent(i+1,inter.Eval(xaxis2->GetBinCenter(i+1)));
  cout << "time of maximum of template:    " << xaxis2->GetBinCenter(interpolated_mean_waveform->GetMaximumBin()) << endl;
  cout << "maximum of template:    " << interpolated_mean_waveform->GetBinContent(interpolated_mean_waveform->GetMaximumBin()) << endl;
  int bin1 = interpolated_mean_waveform->FindFirstBinAbove(interpolated_mean_waveform->GetMaximum()/2);
  int bin2 = interpolated_mean_waveform->FindLastBinAbove(interpolated_mean_waveform->GetMaximum()/2);
  cout << "FWHM of template:    " << interpolated_mean_waveform->GetBinCenter(bin2) - interpolated_mean_waveform->GetBinCenter(bin1) << endl;
  mean_waveform->Draw();
  interpolated_mean_waveform->SetLineColor(4);
  interpolated_mean_waveform->Draw("same");

  //Next Crystal


//const int nTempBins = 100*nsamples;
  TAxis *xaxis_ = mean_waveform2->GetXaxis();
  TH1F *interpolated_mean_waveform2 = new TH1F("interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center","interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nTempBins,-0.125,937.375);
  TAxis *xaxis2_ = interpolated_mean_waveform2->GetXaxis();
  double x2[nsamples+2],y2[nsamples+2];
  x2[0] = xaxis_->GetBinCenter(1) - 6.25;
  y2[0] = 0;
  for (int i = 0;i < nsamples;i++){
    //x[i] = 6.25 * i;
    x2[i+1] = xaxis_->GetBinCenter(i+1);
    y2[i+1] = mean_waveform2->GetBinContent(i+1);
  }
  x2[nsamples+1] = x2[nsamples]+6.25;
  y2[nsamples+1] = 0;
  //ROOT::Math::Interpolator inter(nsamples+1, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator inter2(nsamples+2, ROOT::Math::Interpolation::kAKIMA);
  inter2.SetData(nsamples+2,x2,y2);
  for (int i = 0;i < nTempBins;i++) interpolated_mean_waveform2->SetBinContent(i+1,inter2.Eval(xaxis2_->GetBinCenter(i+1)));
  cout << "time of maximum of template:    " << xaxis2_->GetBinCenter(interpolated_mean_waveform2->GetMaximumBin()) << endl;
  cout << "maximum of template:    " << interpolated_mean_waveform2->GetBinContent(interpolated_mean_waveform2->GetMaximumBin()) << endl;
  int bin1_ = interpolated_mean_waveform2->FindFirstBinAbove(interpolated_mean_waveform2->GetMaximum()/2);
  int bin2_ = interpolated_mean_waveform2->FindLastBinAbove(interpolated_mean_waveform2->GetMaximum()/2);
  cout << "FWHM of template:    " << interpolated_mean_waveform2->GetBinCenter(bin2_) - interpolated_mean_waveform2->GetBinCenter(bin1_) << endl;
  mean_waveform2->Draw("same");
  interpolated_mean_waveform2->SetLineColor(3);
  interpolated_mean_waveform2->Draw("same");


  //End Crystal



  TFile *template_file = new TFile("/afs/cern.ch/user/k/kyee/work/pulsetest/testing.root","recreate");
  template_file->cd();
  mean_waveform->Write();
  mean_waveform2->Write();
  template_file->Close();

  /*TFile *template_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/template_July2017_C3_100.root","recreate");
  template_file->cd();
  mean_waveform->Write();
  template_file->Close();*/
}
