int tresfit(void)
{

 c1 = new TCanvas("c1","Canvas",200,10,700,500);

  const Int_t n = 5;

  Double_t A[n] = {238, 518, 1271, 2690, 4025}; //Using A_eff instead of A
  Double_t sig_n = 19.04; // = sqrt(13.22^2 + 13.75^2) = 19.04
  Double_t x[n];

  for(int i = 0; i < n; i++) x[i] = A[i]/sig_n;
  for(int i = 0; i < n; i++) cout << x[i] << endl;
 
  Double_t y[n] = {0.35, 0.19, 0.078, 0.048, 0.039}; 
  Double_t ey[n] = {0.04, 0.02, .003,.002,.002};

  for(int i = 0; i < n; i++) y[i] = y[i]*1000/sqrt(2);
  for(int i = 0; i < n; i++) ey[i] = ey[i]*1000/sqrt(2);

  Double_t ex[n] = {0, 0, 0, 0, 0};

  gr = new TGraphErrors(n,x,y,ex,ey);
  gr->GetXaxis()->SetTitle( "A_eff/sigma_noise ");
  gr->GetYaxis()->SetTitle(" Sigma of t_max[C3] - t_max[C4] in ps");

  gr->SetTitle("Timing Resolution");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");

  TF1 *myfit = new TF1("myfit","sqrt( pow([0],2) + pow(([1]/x) , 2) )", 0, 350);
  myfit->SetParName(0,"a");
  myfit->SetParName(1,"b");
  myfit->SetParameter(0, 30);
  myfit->SetParameter(1, 10000);
  myfit->SetLineColor(4);
  gr->Fit("myfit","Q");

  cout << "a = " << myfit->GetParameter(0) << " +/- " << myfit->GetParError(0) << endl;
  cout << "b = " << myfit->GetParameter(1) << " +/- " << myfit->GetParError(1) << endl;

  return 0;

}
