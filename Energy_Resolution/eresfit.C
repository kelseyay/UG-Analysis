int eresfit(void)
{

 c1 = new TCanvas("c1","Canvas",200,10,700,500);

  const Int_t n = 5;

  Double_t x[n] = {10, 20, 50, 100, 150};
  Double_t y[n];

  Double_t sig[n] = {41.7, 43.9, 50.38, 60.74, 69.17};
  Double_t peak[n] = {444, 878, 2290, 4672, 6987};

  Double_t ex[n] = {0, 0, 0, 0, 0};
  Double_t ey[n] = {0.007, 0.003, .001,.0005,.0002};

  for(int i = 0; i < n; i++) y[i] = sqrt( pow(sig[i],2) - pow(39.25,2) ) / peak[i] ; 

//  for(int i = 0; i < n; i++)  cout << x[i]*4626/100 << endl;
  for(int i = 0; i < n; i++)  cout << y[i] << endl;

//  TF1 *myfit = new TF1("myfit","sqrt(pow([0],2) + pow([1],2)/x)", 0, 160);
//RMS 13.38

  gr = new TGraphErrors(n,x,y,ex,ey);
  gr->GetXaxis()->SetTitle( "Beam  Energy GeV ");
  gr->GetYaxis()->SetTitle(" Sigma/Mean  ");
  gr->GetYaxis()->SetRangeUser(0, 0.04);

  gr->SetTitle("C3 Sigma/Mean");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");

  TF1 *myfit = new TF1("myfit","sqrt(pow([0],2) + pow([1],2)/(x))", 0, 160);
  myfit->SetParName(0,"a");
  myfit->SetParName(1,"b");
  myfit->SetParLimits(0, 0, 100);
  myfit->SetParLimits(1, 0, 100);
  myfit->SetParameter(0, 0.001);
  myfit->SetParameter(1, 30);
  myfit->SetLineColor(4);
  gr->Fit("myfit","Q");

  cout << "a = " << myfit->GetParameter(0) << " +/- " << myfit->GetParError(0) << endl;
  cout << "b = " << myfit->GetParameter(1) << " +/- " << myfit->GetParError(1) << endl;

  return 0;

}
