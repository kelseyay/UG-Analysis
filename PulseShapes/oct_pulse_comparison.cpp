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

#include <TGraph.h>

using namespace std;

const int nsamples = 150;
const int nboards = 5;
const int nchannels_per_board = 5;
const int nchannels = nboards*nchannels_per_board;
const int nTempBins = 100*nsamples;
float tmax, tmax2;
float amp_max, amp_max2;

void oct_pulse_comparison(){

  TCanvas *canvas = new TCanvas("Pulse Shape Comparison","Pulse Shape Comparison");
  TFile *template_file = new TFile("/afs/cern.ch/user/k/kyee/work/template_new_9900.root");
  TProfile *mean_waveform = (TProfile*) template_file->Get("mean waveform - C3 - 50 GeV - 6*6 mm^2 at center"); //Note: This title MUST align with what the template generator set it to
  TProfile *interpolated_mean_waveform = (TProfile*) template_file->Get("interpolated mean waveform - C3 - 50 GeV - 6*6 mm^2 at center");

  TFile *template_file2 = new TFile("/afs/cern.ch/user/k/kyee/work/template_new_9906.root");
  TProfile *mean_waveform2 = (TProfile*) template_file2->Get("mean waveform - C3 - 100 GeV - 6*6 mm^2 at center");
  TProfile *interpolated_mean_waveform2 = (TProfile*) template_file2->Get("interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center");

  TAxis *xaxis_ = interpolated_mean_waveform->GetXaxis();
  TAxis *xaxis2_ = interpolated_mean_waveform2->GetXaxis();
  cout << "No segvi yet" << endl;
  tmax = xaxis_->GetBinCenter(interpolated_mean_waveform->GetMaximumBin());
  amp_max = interpolated_mean_waveform->GetBinContent(interpolated_mean_waveform->GetMaximumBin());
//  cout << "tmax of C3 is " << tmax << endl;
  tmax2 = xaxis2_->GetBinCenter(interpolated_mean_waveform2->GetMaximumBin());
//  cout << "tmax of D3 is " << tmax2 << endl; //Both are correct
  amp_max2 = interpolated_mean_waveform2->GetBinContent(interpolated_mean_waveform2->GetMaximumBin());

  double x[nsamples+1],y[nsamples+1];
  for (int i = 0;i < nsamples+1;i++){
    x[i] = 6.25 * i + 2.65625 - tmax; //I can't remember why I did this... I think it's to center the peaks on 0
    if (i == nsamples) y[i] = 0;
    else y[i] = mean_waveform->GetBinContent(i+1) / amp_max;
    if(y[i] > 0.95) cout << "x value is " << x[i] << endl;
  }

  double x2[nsamples+1],y2[nsamples+1];
  for (int i = 0;i < nsamples+1;i++){
    x2[i] = 6.25 * i + 2.65625 - tmax2;
    if (i == nsamples) y2[i] = 0;
    else y2[i] = mean_waveform2->GetBinContent(i+1) / amp_max2;
  }

  double xave[nsamples+1],ydiff[nsamples+1];
  for (int i = 0;i < nsamples+1;i++){
    xave[i] = ((x[i] + x2[i])/2); //For October it looks like this needs to be changed.
    ydiff[i] = TMath::Abs( y[i] - y2[i] );
    //if(ydiff[i] > 0.005) cout << "x value is " << xave[i] << endl;
  }


  TH1F *new_interpolation = new TH1F("interpolated mean waveform - C3 - 50 GeV - scaled","interpolated mean waveform - C3 - 50 GeV - scaled",nTempBins,-0.125 - tmax,937.375 - tmax);
  TH1F *new_interpolation2 = new TH1F("interpolated mean waveform - C3 - 100 GeV - scaled","interpolated mean waveform - C3 - 100 GeV - scaled",nTempBins,-0.125 - tmax2,937.375 - tmax2);
  TH1F *new_int_diff = new TH1F("interpolated mean waveform - Oct 50/100 Difference","interpolated mean waveform - Oct 50/100 Difference",nTempBins,-0.125 - tmax,937.375 - tmax);

//  for (int i = 0;i < nTempBins;i++) new_interpolation->SetBinContent(i+1,interpolated_mean_waveform->GetBinContent(i+1) / amp_max);
//  for (int i = 0;i < nTempBins;i++) new_interpolation2->SetBinContent(i+1,interpolated_mean_waveform2->GetBinContent(i+1) / amp_max2);


  TAxis *xaxis = new_interpolation->GetXaxis();
  TAxis *xaxis2 = new_interpolation2->GetXaxis();

  int t_1, t_3;
  double t_2, t_4;

  ROOT::Math::Interpolator inter(nsamples+1, ROOT::Math::Interpolation::kAKIMA);
  inter.SetData(nsamples+1,x,y);
  for (int i = 0;i < nTempBins;i++) new_interpolation->SetBinContent(i+1,inter.Eval(xaxis->GetBinCenter(i+1)));
  t_1 = new_interpolation->GetMaximumBin();
  t_2 = new_interpolation->GetXaxis()->GetBinCenter(t_1);
  cout << "Maximum of the interpolation(C3) is " << t_2 << " The bin number is this " << t_1 << endl;

  ROOT::Math::Interpolator inter2(nsamples+1, ROOT::Math::Interpolation::kAKIMA);
  inter2.SetData(nsamples+1,x2,y2);
  for (int i = 0;i < nTempBins;i++) new_interpolation2->SetBinContent(i+1,inter2.Eval(xaxis2->GetBinCenter(i+1)));

  t_3 = new_interpolation2->GetMaximumBin();
  t_4 = new_interpolation2->GetXaxis()->GetBinCenter(t_3);
  cout << "Maximum of the interpolation(C3) is " << t_4 << " The bin number is this " << t_3 << endl;

  int offset = t_3 - t_1;

//  double new_int_diff[nTempBins];
//  for (int i = 0;i < nTempBins;i++) new_int_diff[i] = TMath::Abs( new_interpolation->GetBinContent(i) - new_interpolation2->GetBinContent(i) ) ;

  for (int i = offset;i < nTempBins - TMath::Abs(offset) ;i++) new_int_diff->SetBinContent(i - offset, TMath::Abs( new_interpolation->GetBinContent(i - offset) - new_interpolation2->GetBinContent(i) ) ); //Note, 456 = Abs(t_3 - t_1) There was a bin offset, and I had to account for that. 

TGraph *graph = new TGraph(nsamples, x, y);
TGraph *graph2 = new TGraph(nsamples, x2, y2);
TGraph *grdiff = new TGraph(nsamples, xave, ydiff);

graph->GetXaxis()->SetRangeUser(-25,55);
graph2->GetXaxis()->SetRangeUser(-25,55);
grdiff->GetXaxis()->SetRangeUser(-25,55);
new_interpolation->GetXaxis()->SetRangeUser(-25,60);
new_interpolation2->GetXaxis()->SetRangeUser(-25,60);
new_int_diff->GetXaxis()->SetRangeUser(-25,60);

//graph->SetLineColor(4);
//graph2->SetLineColor(3);


graph->SetMarkerColor(4);
graph->SetMarkerStyle(1);

graph2->SetMarkerColor(2);
graph2->SetMarkerStyle(1);

graph->Draw("AP");
graph2->Draw("same P");
//grdiff->Draw("same");

new_interpolation->SetLineColor(38);
new_interpolation->Draw("same");

new_int_diff->Draw("same");

new_interpolation2->SetLineColor(46);
new_interpolation2->Draw("same");
//In the future, I want to figure out the spline for the scaled plots or scale the interpolation plots the same way (which might be easier)

}
