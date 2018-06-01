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

//Authored by Arash Jofrehei. Modified by Kelsey Yee.

using namespace std;
const int nsamples = 150;
const int nboards = 3;
const int nchannels_per_board = 5;
const int nchannels = nboards*nchannels_per_board;
const int impinged_channel = 1;
const float beam_energy = 150.0;

float amp_max[nchannels];
float time_max[nchannels];
float maximum[nchannels];

float WF_time[nsamples*nchannels];
float WF_val[nsamples*nchannels];
//float WF_ch[nsamples*nchannels];
float X[2];
float Y[2];
float hodoX[2];
float hodoY[2];
float xcenter = -2;
float ycenter = -5;

void Apr_generating_template(){
  TCanvas *canvas = new TCanvas("generating template","generating template");
  TFile *file = new TFile("apr_35_10280.root");
  //TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_55_C3_100.root");
  TTree *mtree = (TTree*) file->Get("h4");
  mtree->SetBranchAddress("WF_time",WF_time);
  mtree->SetBranchAddress("WF_val",WF_val);
  TTree *dtree = (TTree*) file->Get("digi");
  dtree->SetBranchAddress("amp_max",amp_max);
  dtree->SetBranchAddress("time_max",time_max);
  dtree->SetBranchAddress("maximum",maximum);
  
  TTree *htree = (TTree*) file->Get("hodo");
  htree->SetBranchAddress("X",hodoX);
  htree->SetBranchAddress("Y",hodoY);
  cout << "Checkpoint 1" << endl;
  TF1 *g = new TF1("g","gaus",0,1000);
  TH1F *waveform = new TH1F("waveform - C2 - 150 GeV - 6*6 mm^2 at center","waveform - C2 - 150 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);
  TProfile *mean_waveform = new TProfile("mean waveform - C2 - 150 GeV - 6*6 mm^2 at center","mean waveform - C2 - 150 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);


  const Long64_t nentries = htree->GetEntriesFast();
  cout << "nentries = " << nentries << endl;
  Long64_t nbytes = 0, nb = 0,ientry;
  float dummy_amp[nentries];
  float dummy_time[nentries];
  bool dummy_constraint[nentries];

  
  float hodo_dummy_X0[nentries];
  float hodo_dummy_Y0[nentries];
  float hodo_dummy_X1[nentries];
  float hodo_dummy_Y1[nentries];
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = htree->LoadTree(jentry);
    nb = htree->GetEntry(jentry);
    hodo_dummy_X0[jentry] = hodoX[0];
    hodo_dummy_Y0[jentry] = hodoY[0];
    hodo_dummy_X1[jentry] = hodoX[1];
    hodo_dummy_Y1[jentry] = hodoY[1];
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = dtree->LoadTree(jentry);
    nb = dtree->GetEntry(jentry);
    dummy_amp[jentry] = maximum[impinged_channel];
    dummy_time[jentry] = time_max[impinged_channel];
   if (TMath::Abs(hodo_dummy_X0[jentry] - xcenter)<5&&TMath::Abs(hodo_dummy_Y0[jentry] - ycenter)<5&& maximum[0]>35&&maximum[13]>35&&maximum[8]>35&&maximum[2]>35&&dummy_amp[jentry] > 4000.0*beam_energy/150.0) {dummy_constraint[jentry] = true; cout << "max " << dummy_amp[jentry] << endl; }
    else dummy_constraint[jentry] = false;
   if(dummy_amp[jentry] > 9000) dummy_constraint[jentry] = false;
  }
  cout << "No segvi yet" << endl;
  bool test = 1;
//  const Long64_t nentries2 = mtree->GetEntriesFast();
//  cout << "nentries2 = " << nentries2 << endl; 

  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);    
    nbytes += nb;

    X[0] = hodo_dummy_X0[jentry];
    Y[0] = hodo_dummy_Y0[jentry];
    X[1] = hodo_dummy_X1[jentry];
    Y[1] = hodo_dummy_Y1[jentry];

    if (TMath::Abs(X[0]- xcenter)<5&&TMath::Abs(Y[0]- ycenter)<5){                                                 //To be tuned
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
      if (g->GetParameter(1) < 260 || g->GetParameter(1) > 290) continue;
      if (g->GetParameter(0) < 2300*beam_energy/100.0) continue;
	//cout << "Fit parameter time?" << g->GetParameter(1) << endl;
      //for (int i = 0;i < nsamples;i++) mean_waveform->Fill(277+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1),WF_val[i+impinged_channel*nsamples]); //cout << "WF_time fill " << 277+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1) << " WF_val fill " << WF_val[i+impinged_channel*nsamples] << endl;} // */ //This does in fact work...
	for (int i = 0;i < nsamples;i++) //Find the average time_max of each of the three populations. 
	{
	    if(g->GetParameter(1) >= 280 && g->GetParameter(1) < 284) mean_waveform->Fill(281.4+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1),WF_val[i+impinged_channel*nsamples]);
	    if(g->GetParameter(1) >= 273 && g->GetParameter(1) < 280) mean_waveform->Fill(276.4+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1),WF_val[i+impinged_channel*nsamples]);
	    if(g->GetParameter(1) >= 268 && g->GetParameter(1) < 273) mean_waveform->Fill(271.1+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1),WF_val[i+impinged_channel*nsamples]);
	}
    }
    }
  }
   //mean_waveform->Draw();

  TAxis *xaxis = mean_waveform->GetXaxis();
  const int nTempBins = 100*nsamples;
  TH1F *interpolated_mean_waveform = new TH1F("interpolated mean waveform - C2 - 150 GeV - 6*6 mm^2 at center","interpolated mean waveform - C2 - 150 GeV - 6*6 mm^2 at center;time(ns)",nTempBins,-0.125,937.375);
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

  TFile *template_file = new TFile("template_new_10280.root","recreate");
  template_file->cd();
  interpolated_mean_waveform->Write();
  mean_waveform->Write();
  template_file->Close();

  /*TFile *template_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/template_July2017_C3_100.root","recreate");
  template_file->cd();
  mean_waveform->Write();
  template_file->Close();*/
}
