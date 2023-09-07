#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void fit_delta(){

  //change to 1 if running the full thing
  int debug_mode = 0;
  Long64_t nentries= 2000000;

  //////////////////////////////
  //setup input files
  //////////////////////////////
  const int nSettings = 6;
  //int nRun[nSettings]={16962};//,20862,20867};//run numbers for delta = 10,13,15%, last 2 from deep
  int nRun[nSettings]={16036,16028,16026,16962,20862,20867};//run numbers for delta = 10,13,15%, last 2 from deep
  double centAng[nSettings]={8.295,7.495,6.8,8.3,9.125,7.70};//deg
  //double centAng[nSettings]={8.3};//,9.125,7.704};//deg
  TString inputroot;
  TString outputhist="output2/delta_fit_plots.root";
  double Mp = 0.93827231;
  double pCent = 8.55*0.9979;
  double Ep[nSettings] = {1.821, 1.609, 1.431, 1.82, 2.040, 1.663};
  //double Ep[nSettings] = {1.82};//, 2.040, 1.663,};
  double Me = 0.000511;
  double eBeam = 10.549355;
  double pBeam = sqrt(pow(eBeam,2.) - Me*Me);
  int nfit=0, npar, nfit_max=2979;//1500;//7962;//nfit_max=nSettings x maxKin 
  int maxKin = 500;//500;//1327;
  int imaxKin=0;
  double Mp_rad = 0.944;//mass of proton in SIMC
  double xcentral[nSettings]={1.026,1.03,1.035,1.028,1.01,1.012};
  double xsigma[nSettings]={0.012,0.014,0.018,0.010,0.01,0.014};
  //double xcentral[nSettings]={1.028};//,1.01,1.012};10 MeV / 10 GeV
  //double xsigma[nSettings]={0.010};//,0.011,0.014};
  //make it even in yptar, 10 bins from -0.025 to 0.025 in increments of 0.005
  int nypbins = 10;
  int counteryp[nSettings][nypbins];
  int nypctmax = 50;
 
  //////////////////////////////
  //setup plots
  //////////////////////////////
  TObjArray HList(0);
  TH1F *h_delta[nSettings];
  TH1F *h_xB[nSettings];
  TH1F *h_xB_precut[nSettings];
  TH1F *h_xBcalc[nSettings];
  TH1F *h_Q2[nSettings];
  TH1F *h_angle[nSettings];
  TH1F *h_W[nSettings];
  TH1F *h_P[nSettings];
  TH1F *hnew_P[nSettings];
  TH1F *h_ytar[nSettings];
  TH1F *h_yptar[nSettings];
  TH2F *h2_delta_vs_W[nSettings];
  TH2F *h2_xptar_vs_W[nSettings];
  TH2F *h2_ypfp_vs_yfp[nSettings];
  TH2F *h2_xpfp_vs_xfp[nSettings];
  TH2F *h2_ypfp_vs_yptar[nSettings];
  TH2F *h2_xpfp_vs_yptar[nSettings];
  TH2F *h2_xfp_vs_yptar[nSettings];
  TH2F *h2_yptar_vs_W[nSettings];
  TH2F *h2_ytar_vs_W[nSettings];
  TH2F *h2_ypfp_vs_W[nSettings];
  TH2F *h2_xfp_vs_W[nSettings];
  TH2F *h2_yfp_vs_W[nSettings];
  TH2F *h2_xpfp_vs_W[nSettings];
  TH1F *h_deltaDiff[nSettings];
  TH1F *h_xpDiff[nSettings];
  TH1F *h_ypDiff[nSettings];
  TH1F *h_xDiff[nSettings];
  TH1F *h_Q2Diff[nSettings];
  TH1F *h_nuDiff[nSettings];
  TH1F *hnew_delta[nSettings];
  TH1F *hnew_xB[nSettings];
  TH1F *hnew_Q2[nSettings];
  TH1F *hnew_W[nSettings];
  TH1F *h_Wdiff[nSettings];
  TH1F *h_deltaDiffSVD[nSettings];
  TH2F *h2new_ypfp_vs_W[nSettings];
  TH2F *h2new_xfp_vs_W[nSettings];
  TH2F *h2new_yfp_vs_W[nSettings];
  TH2F *h2new_xpfp_vs_W[nSettings];
  TH1F *h_deltaDiff_pre[nSettings];
  TH1F *h_deltaDiff_post[nSettings];
  TH2F *h2new_delta_vs_W[nSettings];
  TH2F *h2new_xptar_vs_W[nSettings];
  TH2F *h2new_yptar_vs_W[nSettings];
  TH2F *h2new_ytar_vs_W[nSettings];
  TH2F *h2_yptar_vs_W_xcut[nSettings];
  TH2F *h2new_yptar_vs_W_xcut[nSettings];
  TH2F *h2_xpfp_ypfp[nSettings];
  TH2F *h2_xptar_yptar[nSettings];

  //////////////////////////////
  //read in the matrix elements
  /////////////////////////////
  string coeffsfilenameD="delta_terms.dat";
  ifstream coeffsfileD(coeffsfilenameD.c_str());
  TString currentlineD;
  int num_recon_termsD=0;

  //this file must have elements in same order as the big matrix elements file
  while( currentlineD.ReadLine(coeffsfileD,kFALSE) && !currentlineD.BeginsWith(" ----") ){
    num_recon_termsD++;
  }

  string coeffsfilename="newfit_cafe4.dat";
  ifstream coeffsfile(coeffsfilename.c_str());
  TString currentline;
  int num_recon_terms=0;

  vector<double> xptarcoeffs;
  vector<double> yptarcoeffs;
  vector<double> ytarcoeffs;
  vector<double> deltacoeffs;
  vector<int> xfpexpon;
  vector<int> xpfpexpon;
  vector<int> yfpexpon;
  vector<int> ypfpexpon;
  vector<int> xtarexpon;
  
  while( currentline.ReadLine(coeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    
    TString sc1(currentline(1,16));
    TString sc2(currentline(17,16));
    TString sc3(currentline(33,16));
    TString sc4(currentline(49,16));
    
    xptarcoeffs.push_back(sc1.Atof());
    ytarcoeffs.push_back(sc2.Atof());
    yptarcoeffs.push_back(sc3.Atof());
    deltacoeffs.push_back(sc4.Atof());
    
    int expontemp[5];
    
    for(int expon=0; expon<5; expon++){
      TString stemp(currentline(66+expon,1));
      expontemp[expon] = stemp.Atoi();
    }
    
    xfpexpon.push_back(expontemp[0]);
    xpfpexpon.push_back(expontemp[1]);
    yfpexpon.push_back(expontemp[2]);
    ypfpexpon.push_back(expontemp[3]);
    xtarexpon.push_back(expontemp[4]);
    
    cout << num_recon_terms << " " <<  xptarcoeffs[num_recon_terms] << " " << ytarcoeffs[num_recon_terms] << " " <<  yptarcoeffs[num_recon_terms] << " " << deltacoeffs[num_recon_terms] << " " << xfpexpon[num_recon_terms] << " " << xpfpexpon[num_recon_terms] << " " << yfpexpon[num_recon_terms] << " " << ypfpexpon[num_recon_terms] << " " << xtarexpon[num_recon_terms] << " " << endl;
    
    num_recon_terms++;   
  }
  cout<<"number of terms to optimize: "<<num_recon_termsD<<endl;
  npar = num_recon_termsD;
  TVectorD b_delta(npar);
  TMatrixD lambda(npar,nfit_max);
  TMatrixD Ay(npar,npar);
  TVectorD b_yptar(npar);
  TVectorD b_xptar(npar);
  
  //////////////////////////////
  //loop the input files
  //////////////////////////////
  for (int ii=0; ii<nSettings; ii++){
  //for (int ii=0; ii<1; ii++){
    //if(ii==0 || ii==1 || ii==3){continue;}
    //for (int ii=0; ii<1; ii++){
    if (ii>3){
      eBeam = 10.549355-0.005;
      pBeam = sqrt(pow(eBeam,2.) - Me*Me);}

    for (int jj=0;jj<nypbins;jj++){
      counteryp[ii][jj]=0;
    }
  //for (int ii=0; ii<1; ii++){
    if(int(nRun[ii])<20000){
      inputroot=Form("Rootfiles_may/cafe_replay_optics_%d_-1.root",nRun[ii]);}
    else{
      inputroot=Form("Rootfiles_may/cafe_replay_optics_%d_500000.root",nRun[ii]);}

    cout << " input root = " << inputroot << endl;
    imaxKin=0;

    //read in variables
    TFile *fsimc = new TFile(inputroot); 
    TTree *tsimc = (TTree*) fsimc->Get("T");
    Double_t  sumnpe,sumhgnpe,etracknorm,ytar,xtar,reactx,reacty,reactz,delta,yptar,xptar,yfp,ypfp,xfp,xpfp,xbpm_tar,ybpm_tar,frx,fry,xbcalc,Q2calc,nucalc,xapos,xbpmtar,ybpmtar,yapos,xbpos,ybpos,thetacalc;
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&xbpmtar);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&ybpmtar);
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpmA",&xapos);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpmA",&yapos); 
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpmB",&xbpos);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpmB",&ybpos);
    tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
    tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
    tsimc->SetBranchAddress("P.cal.etottracknorm",&etracknorm);
    tsimc->SetBranchAddress("P.gtr.y",&ytar);
    tsimc->SetBranchAddress("P.gtr.x",&xtar);
    tsimc->SetBranchAddress("P.react.x",&reactx);
    tsimc->SetBranchAddress("P.react.y",&reacty);
    tsimc->SetBranchAddress("P.react.z",&reactz);
    tsimc->SetBranchAddress("P.gtr.dp",&delta);
    tsimc->SetBranchAddress("P.gtr.ph",&yptar);
    tsimc->SetBranchAddress("P.gtr.th",&xptar);
    tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
    tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
    tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
    tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
    tsimc->SetBranchAddress("P.kin.primary.x_bj",&xbcalc);
    tsimc->SetBranchAddress("P.kin.primary.Q2",&Q2calc);
    tsimc->SetBranchAddress("P.kin.primary.nu",&nucalc);
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&xbpm_tar);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&ybpm_tar);
    tsimc->SetBranchAddress("P.rb.raster.fr_xa",&frx);
    tsimc->SetBranchAddress("P.rb.raster.fr_ya",&fry);
    tsimc->SetBranchAddress("P.kin.primary.scat_ang_deg",&thetacalc);


    //////////////////////////////
    //make the plot
    //////////////////////////////
    h_delta[ii] = new TH1F(Form("h_delta_%d",ii),Form("Run %d; Delta; Counts",nRun[ii]),100,-5,20);
    HList.Add(h_delta[ii]);
    h_xB[ii] = new TH1F(Form("h_xB_%d",ii),Form("Run %d; xB; Counts",nRun[ii]),100,0.9,1.15);
    HList.Add(h_xB[ii]);
    h_xB_precut[ii] = new TH1F(Form("h_xB_precut_%d",ii),Form("Run %d; xB-precut; Counts",nRun[ii]),100,0.9,1.15);
    HList.Add(h_xB_precut[ii]);
    h_xBcalc[ii] = new TH1F(Form("h_xBcalc_%d",ii),Form("Run %d; xB-replay,cut; Counts",nRun[ii]),100,0.9,1.15);
    HList.Add(h_xBcalc[ii]);
    h_Q2[ii] = new TH1F(Form("h_Q2_%d",ii),Form("Run %d; Q2; Counts",nRun[ii]),100,0,6);
    HList.Add(h_Q2[ii]);
    h_angle[ii] = new TH1F(Form("h_angle_%d",ii),Form("Run %d; scattered e- angle; Counts",nRun[ii]),100,0,20);
    HList.Add(h_angle[ii]);
    h_W[ii] = new TH1F(Form("h_W_%d",ii),Form("Run %d; W; Counts",nRun[ii]),100,0.9,1.);
    HList.Add(h_W[ii]);
    h_P[ii] = new TH1F(Form("h_P_%d",ii),Form("Run %d; P [GeV]; Counts",nRun[ii]),100,0,12.0);
    HList.Add(h_P[ii]);
    h_ytar[ii] = new TH1F(Form("h_ytar_%d",ii),Form("Run %d; yTarget; Counts",nRun[ii]),100,-1.5,1.5);
    HList.Add(h_ytar[ii]);
    h_yptar[ii] = new TH1F(Form("h_yptar_%d",ii),Form("Run %d; ypTar; Counts",nRun[ii]),100,-0.1,0.1);
    HList.Add(h_yptar[ii]);
    h2_delta_vs_W[ii] = new TH2F(Form("h2_delta_vs_W_%d",ii),Form("Run %d; Delta; W",nRun[ii]),100,-5,20,100,0.8,1);
    HList.Add(h2_delta_vs_W[ii]);
    h2_xptar_vs_W[ii] = new TH2F(Form("h2_xptar_vs_W_%d",ii),Form("Run %d; xptar; W",nRun[ii]),100,-0.06,0.06,100,0.8,1);
    HList.Add(h2_xptar_vs_W[ii]);
    h2_ypfp_vs_yfp[ii] = new TH2F(Form("h2_ypfp_vs_yfp_%d",ii),Form("Run %d; ypfp; yfp",nRun[ii]),100,-0.04,0.04,100,-40.0,40.0);
    HList.Add(h2_ypfp_vs_yfp[ii]);
    h2_xpfp_vs_xfp[ii] = new TH2F(Form("h2_xpfp_vs_xfp_%d",ii),Form("Run %d; xpfp; xfp",nRun[ii]),100,-0.1,0.1,100,-40.0,40.0);
    HList.Add(h2_xpfp_vs_xfp[ii]);
    h2_ypfp_vs_yptar[ii] = new TH2F(Form("h2_ypfp_vs_yptar_%d",ii),Form("Run %d; ypfp; ypTar",nRun[ii]),100,-0.04,0.04,100,-0.05,0.05);
    HList.Add(h2_ypfp_vs_yptar[ii]);
    h2_xpfp_vs_yptar[ii] = new TH2F(Form("h2_xpfp_vs_yptar_%d",ii),Form("Run %d; xpfp; ypTar",nRun[ii]),100,-0.1,0.1,100,-0.05,0.05);
    HList.Add(h2_xpfp_vs_yptar[ii]);
    h2_xfp_vs_yptar[ii] = new TH2F(Form("h2_xfp_vs_yptar_%d",ii),Form("Run %d; xfp; ypTar",nRun[ii]),100,-40.0,40.0,100,-0.05,0.05);
    HList.Add(h2_xfp_vs_yptar[ii]);
    h2_yptar_vs_W[ii] = new TH2F(Form("h2_yptar_vs_W_%d",ii),Form("Run %d; ypTar; W",nRun[ii]),100,-0.05,0.05,100,0.8,1);
    HList.Add(h2_yptar_vs_W[ii]);
    h2_yptar_vs_W_xcut[ii] = new TH2F(Form("h2_yptar_vs_W_xcut_%d",ii),Form("Run %d (xB cut); ypTar; W",nRun[ii]),100,-0.05,0.05,100,0.8,1);
    HList.Add(h2_yptar_vs_W_xcut[ii]);
    h2_ytar_vs_W[ii] = new TH2F(Form("h2_ytar_vs_W_%d",ii),Form("Run %d; yTar; W",nRun[ii]),100,-1.5,1.5,100,0.8,1);
    HList.Add(h2_ytar_vs_W[ii]);
    h2_ypfp_vs_W[ii] = new TH2F(Form("h2_ypfp_vs_W_%d",ii),Form("Run %d; ypfp; W",nRun[ii]),100,-0.04,0.04,100,0.8,1);
    HList.Add(h2_ypfp_vs_W[ii]);
    h2_xpfp_vs_W[ii] = new TH2F(Form("h2_xpfp_vs_W_%d",ii),Form("Run %d; xpfp; W",nRun[ii]),100,-0.1,0.1,100,0.8,1);
    HList.Add(h2_xpfp_vs_W[ii]);
    h2_xfp_vs_W[ii] = new TH2F(Form("h2_xfp_vs_W_%d",ii),Form("Run %d; xfp; W",nRun[ii]),100,-40,40,100,0.8,1);
    HList.Add(h2_xfp_vs_W[ii]);
    h2_yfp_vs_W[ii] = new TH2F(Form("h2_yfp_vs_W_%d",ii),Form("Run %d; yfp; W",nRun[ii]),100,-40,40,100,0.8,1);
    HList.Add(h2_yfp_vs_W[ii]);
    h_deltaDiff[ii] = new TH1F(Form("h_deltaDiff_%d",ii),Form("Run %d; Delta difference; Counts",nRun[ii]),100,-3,3);
    HList.Add(h_deltaDiff[ii]);
    h_xpDiff[ii] = new TH1F(Form("h_xpDiff_%d",ii),Form("Run %d; xptar difference; Counts",nRun[ii]),100,-0.01,0.01);
    HList.Add(h_xpDiff[ii]);
    h_ypDiff[ii] = new TH1F(Form("h_ypDiff_%d",ii),Form("Run %d; yptar difference; Counts",nRun[ii]),100,-0.006,0.006);
    HList.Add(h_ypDiff[ii]);
    h_xDiff[ii] = new TH1F(Form("h_xDiff_%d",ii),Form("Run %d; xB difference; Counts",nRun[ii]),100,-0.1,0.1);
    HList.Add(h_xDiff[ii]);
    h_Q2Diff[ii] = new TH1F(Form("h_Q2Diff_%d",ii),Form("Run %d; Q2 difference; Counts",nRun[ii]),100,-1,1);
    HList.Add(h_Q2Diff[ii]);
    h_nuDiff[ii] = new TH1F(Form("h_nuDiff_%d",ii),Form("Run %d; nu difference; Counts",nRun[ii]),100,-0.05,0.05);
    HList.Add(h_nuDiff[ii]);
    h_deltaDiffSVD[ii] = new TH1F(Form("h_deltaDiffSVD_%d",ii),Form("Run %d; Delta true - delta; Counts",nRun[ii]),100,-0.01,0.01);
    HList.Add(h_deltaDiffSVD[ii]);
    h_deltaDiff_pre[ii] = new TH1F(Form("h_deltaDiff_pre_%d",ii),Form("Run %d; delta true - delta (orig); Counts",nRun[ii]),100,-3,3);
    HList.Add(h_deltaDiff_pre[ii]);
    h2_xpfp_ypfp[ii] = new TH2F(Form("h2_xpfp_ypfp_%d",ii),Form("Run %d; xpfp; ypfp",nRun[ii]),100,0,0.08,100,-0.01,0.01);
    HList.Add(h2_xpfp_ypfp[ii]);
    h2_xptar_yptar[ii] = new TH2F(Form("h2_xptar_yptar_%d",ii),Form("Run %d; xptar; yptar",nRun[ii]),100,-0.05,0.05,100,-0.05,0.05);
    HList.Add(h2_xptar_yptar[ii]);
    //////////////////////////////////////////////////////
    //choose events (cuts), make initial plots of events
    //////////////////////////////////////////////////////
    if (debug_mode !=0 ){nentries = tsimc->GetEntries();}
    for (int i = 0; i < nentries; i++) {
      tsimc->GetEntry(i);
      if (i%50000==0) cout << " Entry = " << i << endl;
      //if (imaxKin >= maxKin){ continue;}
      if (sumnpe>6.0 && etracknorm>0.8){

	//thetacalc -= 4.*0.05;//1mrad=0.057
	//////////////////////////////////////////////////////
	//calculate the reconstructed quantities
	//////////////////////////////////////////////////////
	Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
	Double_t etemp;
	for( int icoeff=0; icoeff<num_recon_terms; icoeff++ ){
	  etemp= 
	    pow( xfp / 100.0, xfpexpon[icoeff] ) * 
	    pow( yfp / 100.0, yfpexpon[icoeff] ) * 
	    pow( xpfp, xpfpexpon[icoeff] ) * 
	    pow( ypfp, ypfpexpon[icoeff] ) * 
	    pow( xtar/100., xtarexpon[icoeff] );
	  deltatemp += deltacoeffs[icoeff] * etemp;
	  ytartemp += ytarcoeffs[icoeff] * etemp;
	  yptartemp += yptarcoeffs[icoeff] * etemp;
	  xptartemp += xptarcoeffs[icoeff] *etemp; 
	} // for icoeffold loop
	
	Double_t delta_per = deltatemp*100.0;
	ytartemp *=100;

	//target angle offsets
	double tt = -(xapos-xbpos)/(320.17-224.81);
	double tp = (yapos-ybpos)/(320.17-224.81);
	double scaling =1./sqrt(1.0+tt*tt+tp*tp);
	double x_tar_angle = scaling*tt;
	double y_tar_angle = scaling*tp;
	//cout<<"xp beam angle: "<<x_tar_angle<<endl;

	double emom = pCent*(1.0+delta_per/100.0);//momentum
	double ep = sqrt(pow(emom,2.) + Me*Me);//energy
	double nu = eBeam - ep;
	
	//double esp_z = emom/sqrt(1.0 + xptartemp*xptartemp + yptartemp*yptartemp);
	//double p_ez = esp_z*(-yptartemp*sin(centAng[ii]*(22./7.)/180.0) + cos(centAng[ii]*(22./7.)/180.0));
	double estheta = thetacalc*(22./7.)/180.;//acos(p_ez/emom)+x_tar_angle;
	double Q2 = 4*eBeam*ep*pow(sin(estheta/2.),2.);
	
	//double Q2 = 4*eBeam*ep*pow(sin(centAng[ii]*(22./7.)/180.0/2.),2.);
	double W = sqrt(Mp*Mp + 2*Mp*nu - Q2);
	double xB = Q2/(2.0*Mp*nu);

	h_xB_precut[ii]->Fill(xbcalc);
	h2_delta_vs_W[ii]->Fill(delta_per,W);
	h2_yptar_vs_W[ii]->Fill(yptartemp,W);
	h2_ytar_vs_W[ii]->Fill(ytartemp,W);
	h2_xptar_vs_W[ii]->Fill(xptartemp,W);
	//cut on xB
	if ((xbcalc>xcentral[ii]-xsigma[ii] && xbcalc<xcentral[ii]+xsigma[ii]) && (imaxKin<maxKin) ){//&& abs(yptartemp)<0.025){
	  //if ((xB>1 && xB<1.01) && (imaxKin<maxKin)){
	  //int flag = 1;//0;
	  //if (flag != 0){
	    
	    h_delta[ii]->Fill(delta_per);
	    h_xB[ii]->Fill(xB);
	    h_xBcalc[ii]->Fill(xbcalc);
	    h_W[ii]->Fill(W);
	    h_Q2[ii]->Fill(Q2);
	    h_angle[ii]->Fill(estheta*180/(22./7));
	    h_P[ii]->Fill(emom); 
	    h_ytar[ii]->Fill(ytartemp); 
	    h_yptar[ii]->Fill(yptartemp); 
	    h2_ypfp_vs_yfp[ii]->Fill(ypfp,yfp);
	    h2_xpfp_vs_xfp[ii]->Fill(xpfp,xfp);
	    h2_ypfp_vs_yptar[ii]->Fill(ypfp,yptartemp);
	    h2_xpfp_vs_yptar[ii]->Fill(xpfp,yptartemp);
	    h2_xfp_vs_yptar[ii]->Fill(xfp,yptartemp);
	    h2_ypfp_vs_W[ii]->Fill(ypfp,W);
	    h2_yfp_vs_W[ii]->Fill(yfp,W);
	    h2_xfp_vs_W[ii]->Fill(xfp,W);
	    h2_xpfp_vs_W[ii]->Fill(xpfp,W);
	    h2_xpfp_ypfp[ii]->Fill(xpfp,ypfp);
	    h2_xptar_yptar[ii]->Fill(xptartemp,yptartemp);
	    
	    //h_deltaDiff[ii]->Fill(delta_per - delta);
	    //h_xpDiff[ii]->Fill(xptartemp - xptar);
	    //h_ypDiff[ii]->Fill(yptartemp - yptar);
	    //h_xDiff[ii]->Fill(xB - xbcalc);
	    //h_Q2Diff[ii]->Fill(Q2 - Q2calc);
	    //h_nuDiff[ii]->Fill(nu - nucalc);
	    
	    //////////////////////////////
	    //calculate the true variables
	    //////////////////////////////
	    double eptrue = (pow(Mp_rad,2.0)-pow(Mp,2.0)-2.0*Mp*eBeam)/(-2.0*Mp-4.0*eBeam*pow(sin(estheta/2.),2.0));
	    double ptrue = sqrt(eptrue*eptrue-Me*Me);
	    double delta_true = ptrue/pCent - 1;
	    h_deltaDiffSVD[ii]->Fill(delta_true - deltatemp);
	    h_deltaDiff_pre[ii]->Fill(delta_true - deltatemp);
	    
	    //assume delta, Eb, xptar, and central angle ok
	    //double theta_check = acos((Mp+eBeam)/eBeam - Mp/ep);//acos(1-(eBeam-ep-1)*Mp/eBeam);
	    //double yptrue = sqrt(pow(theta_check,2.)-pow(xptartemp,2))-centAng[ii]*(22./7.)/180.0;
	    //double xptrue = sqrt(abs(pow(theta_check,2.)-pow(centAng[ii]*(22./7.)/180.0+yptartemp,2)));
	    //double eta = 2*asin(sqrt((pow(Mp_rad,2.)-pow(Mp,2.0)-2*Mp*nu)/(-4*eBeam*ep)));
	    //double xptrue = sqrt(pow((-yptartemp*centAng[ii]+cos(centAng[ii]*(22./7.)/180.)/cos(eta)),2)-pow(yptartemp,2)-1);
	    //if (xptartemp<0){
	    //xptrue *= -1.0;
	    //}
	    
	    //reconstruct it and fill the matrices
	    int flag = 0;
	    
	    for (int jj=0; jj<nypbins; jj++){
	      //double ypmin = -0.01+0.002*jj;
	      //double ypmax = -0.01+0.002*(jj+1.0);
	      //double ypmin = -0.025+0.005*jj;
	      //double ypmax = -0.025+0.005*(jj+1.0);
	      double ypmin = -0.03+0.006*jj;
	      double ypmax = -0.03+0.006*(jj+1.0);
	      if (yptartemp>ypmin && yptartemp<ypmax){
		if(counteryp[ii][jj]<nypctmax){
		  counteryp[ii][jj]++; flag = 1;}}
	    }
	    
	    if (flag!=0){
	      etemp = 0.0;
	      for( int icoeff_fit=0; icoeff_fit<num_recon_termsD; icoeff_fit++ ){
		etemp= 
		  pow( xfp / 100.0, xfpexpon[icoeff_fit] ) * 
		  pow( yfp / 100.0, yfpexpon[icoeff_fit] ) * 
		  pow( xpfp, xpfpexpon[icoeff_fit] ) * 
		  pow( ypfp, ypfpexpon[icoeff_fit] ) * 
		  pow( xtar/100., xtarexpon[icoeff_fit] );
		if (nfit < nfit_max) {
		  lambda[icoeff_fit][nfit] = etemp;
		  b_delta[icoeff_fit] += (delta_true - deltatemp) * etemp; 
		  //cout<<"\t nfit, nfit_max: "<<nfit<<","<<nfit_max<<endl;
		  //nfit++;
		}
		h2_yptar_vs_W_xcut[ii]->Fill(yptartemp,W);
		//b_yptar[icoeff_fit] += (yptrue - yptartemp) * etemp;
		//b_xptar[icoeff_fit] += (xptrue - xptartemp) * etemp;
	      } // for icoeff_fit loop
	      nfit++;
	      imaxKin++;

	    }//end if flag!=0
	}// cut on xB
      }//cut on npesum
    }//loop entries
    cout<<"Number passed: "<<imaxKin<<endl;
  }//loop run kinematics
  
  //////////////////////////////
  //setup the SVD
  //////////////////////////////
  
  if (nfit < nfit_max) {
    cout << " nfit < nfit_max, set nfit_max = " << nfit << endl;
    return;
  }
  //
  cout << " number to fit = " << nfit << " max = " << nfit_max << endl;
  for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
  }
  for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 5000 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
      	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
      }
    }
  }
  
  TDecompSVD Ay_svd(Ay);
  bool ok;
  
  ok = Ay_svd.Solve( b_delta );
  cout << "delta solution ok = " << ok << endl;
  b_delta.Print();
  Ay.Print();
  /*   
  for (int ii=0; ii<nSettings; ii++){
    cout<<"setting: "<<ii<<endl;
    for (int jj=0; jj<nypbins; jj++){
      cout<<"   "<<(-0.03+jj*0.002)<<" "<<counteryp[ii][jj]<<endl;
    }
    }
  */
  /*
  ok = Ay_svd.Solve( b_yptar );
  cout << "yptar solution ok = " << ok << endl;
  b_yptar.Print();
  */
  /*
  ok = Ay_svd.Solve( b_xptar );
  cout << "xptar solution ok = " << ok << endl;
  b_xptar.Print();
  */
  
  //////////////////////////////
  //reconstruct results from SVD
  //////////////////////////////
  nfit=0;
  //for (int ii=0; ii<1; ii++){

    for (int ii=0; ii<nSettings; ii++){

    //for (int ii=0; ii<1; ii++){
    //if(ii!=0){continue;}

      if (ii>3){
      eBeam = 10.549355-0.005;
      pBeam = sqrt(pow(eBeam,2.) - Me*Me);}

    //  for (int ii=0; ii<1; ii++){
    if(int(nRun[ii])<20000){
      inputroot=Form("Rootfiles_may/cafe_replay_optics_%d_-1.root",nRun[ii]);}
    else{
      inputroot=Form("Rootfiles_may/cafe_replay_optics_%d_500000.root",nRun[ii]);}
    //inputroot=Form("Rootfiles_may/cafe_replay_optics_%d_-1.root",nRun[ii]);
    cout << " input root = " << inputroot << endl;
    imaxKin=0;
    
    //read in variables
    TFile *fsimc = new TFile(inputroot); 
    TTree *tsimc = (TTree*) fsimc->Get("T");
    Double_t  sumnpe,sumhgnpe,etracknorm,ytar,xtar,reactx,reacty,reactz,delta,yptar,xptar,yfp,ypfp,xfp,xpfp,xbpm_tar,ybpm_tar,frx,fry,xbcalc,Q2calc,nucalc,xapos,xbpmtar,ybpmtar,yapos,xbpos,ybpos,thetacalc;
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&xbpmtar);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&ybpmtar);
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpmA",&xapos);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpmA",&yapos); 
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpmB",&xbpos);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpmB",&ybpos);
    tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
    tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
    tsimc->SetBranchAddress("P.cal.etottracknorm",&etracknorm);
    tsimc->SetBranchAddress("P.gtr.y",&ytar);
    tsimc->SetBranchAddress("P.gtr.x",&xtar);
    tsimc->SetBranchAddress("P.react.x",&reactx);
    tsimc->SetBranchAddress("P.react.y",&reacty);
    tsimc->SetBranchAddress("P.react.z",&reactz);
    tsimc->SetBranchAddress("P.gtr.dp",&delta);
    tsimc->SetBranchAddress("P.gtr.ph",&yptar);
    tsimc->SetBranchAddress("P.gtr.th",&xptar);
    tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
    tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
    tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
    tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
    tsimc->SetBranchAddress("P.kin.primary.x_bj",&xbcalc);
    tsimc->SetBranchAddress("P.kin.primary.Q2",&Q2calc);
    tsimc->SetBranchAddress("P.kin.primary.nu",&nucalc);
    tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&xbpm_tar);
    tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&ybpm_tar);
    tsimc->SetBranchAddress("P.rb.raster.fr_xa",&frx);
    tsimc->SetBranchAddress("P.rb.raster.fr_ya",&fry);
    tsimc->SetBranchAddress("P.kin.primary.scat_ang_deg",&thetacalc);

    
    hnew_delta[ii] = new TH1F(Form("hnew_delta_%d",ii),Form("Run %d; new Delta; Counts",nRun[ii]),100,-5,20);
    HList.Add(hnew_delta[ii]);
    hnew_xB[ii] = new TH1F(Form("hnew_xB_%d",ii),Form("Run %d; new xB; Counts",nRun[ii]),100,0.9,1.15);
    HList.Add(hnew_xB[ii]);
    hnew_Q2[ii] = new TH1F(Form("hnew_Q2_%d",ii),Form("Run %d; new Q2; Counts",nRun[ii]),100,0,6);
    HList.Add(hnew_Q2[ii]);
    hnew_W[ii] = new TH1F(Form("hnew_W_%d",ii),Form("Run %d; new W; Counts",nRun[ii]),100,0.9,1.);
    HList.Add(hnew_W[ii]);
    h2new_ypfp_vs_W[ii] = new TH2F(Form("h2new_ypfp_vs_W_%d",ii),Form("Run %d; ypfp; W",nRun[ii]),100,-0.04,0.04,100,0.8,1);
    HList.Add(h2new_ypfp_vs_W[ii]);
    h2new_xpfp_vs_W[ii] = new TH2F(Form("h2new_xpfp_vs_W_%d",ii),Form("Run %d; xpfp; W",nRun[ii]),100,-0.1,0.1,100,0.8,1);
    HList.Add(h2new_xpfp_vs_W[ii]);
    h2new_xfp_vs_W[ii] = new TH2F(Form("h2new_xfp_vs_W_%d",ii),Form("Run %d; xfp; W",nRun[ii]),100,-40,40,100,0.8,1);
    HList.Add(h2new_xfp_vs_W[ii]);
    h2new_yfp_vs_W[ii] = new TH2F(Form("h2new_yfp_vs_W_%d",ii),Form("Run %d; yfp; W",nRun[ii]),100,-40,40,100,0.8,1);
    HList.Add(h2new_yfp_vs_W[ii]);
    hnew_P[ii] = new TH1F(Form("hnew_P_%d",ii),Form("Run %d; new P [GeV]; Counts",nRun[ii]),100,0,12.0);
    HList.Add(hnew_P[ii]);
    h_deltaDiff_post[ii] = new TH1F(Form("h_deltaDiff_post_%d",ii),Form("Run %d; delta true - delta (new); Counts",nRun[ii]),100,-3,3);
    HList.Add(h_deltaDiff_post[ii]);
    h2new_delta_vs_W[ii] = new TH2F(Form("h2new_delta_vs_W_%d",ii),Form("Run %d; Delta-new; W-new",nRun[ii]),100,-5,20,100,0.8,1);
    HList.Add(h2new_delta_vs_W[ii]);
    h2new_xptar_vs_W[ii] = new TH2F(Form("h2new_xptar_vs_W_%d",ii),Form("Run %d; xptar; W-new",nRun[ii]),100,-0.06,0.06,100,0.8,1);
    HList.Add(h2new_xptar_vs_W[ii]);
    h2new_yptar_vs_W[ii] = new TH2F(Form("h2new_yptar_vs_W_%d",ii),Form("Run %d; yptar; W-new",nRun[ii]),100,-0.05,0.05,100,0.8,1);
    HList.Add(h2new_yptar_vs_W[ii]);
    h2new_yptar_vs_W_xcut[ii] = new TH2F(Form("h2new_yptar_vs_W_xcut_%d",ii),Form("Run %d (xB cut); ypTar; W-new",nRun[ii]),100,-0.05,0.05,100,0.8,1);
    HList.Add(h2new_yptar_vs_W_xcut[ii]);
    h2new_ytar_vs_W[ii] = new TH2F(Form("h2new_ytar_vs_W_%d",ii),Form("Run %d; ytar; W-new",nRun[ii]),100,-1.5,1.5,100,0.8,1);
    HList.Add(h2new_ytar_vs_W[ii]);
    
    if (debug_mode !=0 ){nentries = tsimc->GetEntries();}
    for (int i = 0; i < nentries; i++) {
      tsimc->GetEntry(i);
      if (i%50000==0) cout << " Entry = " << i << endl;
      if (sumnpe>6.0 && etracknorm>0.8){

	Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
	Double_t etemp;
	for( int icoeff=0; icoeff<num_recon_terms; icoeff++ ){
	  etemp= 
	    pow( xfp / 100.0, xfpexpon[icoeff] ) * 
	    pow( yfp / 100.0, yfpexpon[icoeff] ) * 
	    pow( xpfp, xpfpexpon[icoeff] ) * 
	    pow( ypfp, ypfpexpon[icoeff] ) * 
	    pow( xtar/100., xtarexpon[icoeff] );
	  deltatemp += deltacoeffs[icoeff] * etemp;
	  ytartemp += ytarcoeffs[icoeff] * etemp;
	  yptartemp += yptarcoeffs[icoeff] * etemp;
	  xptartemp += xptarcoeffs[icoeff] *etemp;
	  if (icoeff<num_recon_termsD){
	    deltatemp += b_delta[icoeff]*etemp;
	    //yptartemp += b_yptar[icoeff]*etemp;
	    //xptartemp += b_xptar[icoeff]*etemp;
	  }
	} // for icoeffold loop

	Double_t delta_per = deltatemp*100.0;
	ytartemp *=100;
	
	//target angle offsets
	double tt = -(xapos-xbpos)/(320.17-224.81);
	double tp = (yapos-ybpos)/(320.17-224.81);
	double scaling =1./sqrt(1.0+tt*tt+tp*tp);
	double x_tar_angle = scaling*tt;
	double y_tar_angle = scaling*tp;

	double emom = pCent*(1.0+delta_per/100.0);//momentum
	double ep = sqrt(pow(emom,2.) + Me*Me);//energy
	double nu = eBeam - ep;
	
	double esp_z = emom/sqrt(1.0 + xptartemp*xptartemp + yptartemp*yptartemp);
	double p_ez = esp_z*(-yptartemp*sin(centAng[ii]*(22./7.)/180.0) + cos(centAng[ii]*(22./7.)/180.0));
	double estheta = thetacalc*(22./7.)/180.;//acos(p_ez/emom)+x_tar_angle;
	double Q2 = 4*eBeam*ep*pow(sin(estheta/2.),2.);
	
	double W = sqrt(Mp*Mp + 2*Mp*nu - Q2);
	double emiss = nu + Mp - Ep[ii];
	double xB = Q2/(2.0*Mp*nu);


	//calc true to compare:
	double eptrue = (-2*Mp*eBeam)/(-2*Mp-4*eBeam*pow(sin(estheta/2.),2.0));
	double ptrue = sqrt(eptrue*eptrue-Me*Me);
	double delta_true = ptrue/pCent - 1;
	
	h2new_delta_vs_W[ii]->Fill(delta_per,W);
	h2new_xptar_vs_W[ii]->Fill(xptartemp,W);
	h2new_yptar_vs_W[ii]->Fill(yptartemp,W);
	h2new_ytar_vs_W[ii]->Fill(ytartemp,W);


	//TODO: cut on xB again
	if ((xbcalc>xcentral[ii]-xsigma[ii] && xbcalc<xcentral[ii]+xsigma[ii]) && (imaxKin<maxKin) ){//&& abs(yptartemp)<0.025){
	//if ((xB>0.985 && xB<1) && (imaxKin<maxKin)){&& abs(xptartemp)<0.01 && abs(yptartemp)<0.01
	  hnew_delta[ii]->Fill(delta_per); 
	  hnew_xB[ii]->Fill(xB);
	  hnew_Q2[ii]->Fill(Q2);
	  hnew_W[ii]->Fill(W);
	  h2new_ypfp_vs_W[ii]->Fill(ypfp,W);
	  h2new_yfp_vs_W[ii]->Fill(yfp,W);
	  h2new_xfp_vs_W[ii]->Fill(xfp,W);
	  h2new_xpfp_vs_W[ii]->Fill(xpfp,W);
	  hnew_P[ii]->Fill(emom);
	  h_deltaDiff_post[ii]->Fill(delta_true - deltatemp);
	  h2new_yptar_vs_W_xcut[ii]->Fill(yptartemp,W);

	  imaxKin++;
	}
      }//cut on npesum
    }//loop entries
  }//loop kin settings (aka files)

  ////////////////////
  // write out coeff
  ///////////////////
  string newcoeffsfilename="new_delta_cafe.dat";
  //string newcoeffsfilename="new_yptar_cafe.dat";
  ofstream newcoeffsfile(newcoeffsfilename.c_str());
  char coeffstring[100];
  Double_t tt;
  cout << "writing new coeffs file" << endl;
  //newcoeffsfile << "! new fit to "<<endl;//+fname << endl;
  //newcoeffsfile << " ---------------------------------------------" << endl;
  for( int icoeff_fit=0; icoeff_fit<num_recon_termsD; icoeff_fit++ ){
    newcoeffsfile << " ";
    sprintf( coeffstring, "%16.9g", 0.0 );
    //sprintf( coeffstring, "%16.9g", b_xptar[icoeff_fit] );
    newcoeffsfile << coeffstring; 
    sprintf( coeffstring, "%16.9g", 0.0 );
    newcoeffsfile << coeffstring;
    sprintf( coeffstring, "%16.9g", 0.0 );
    //sprintf( coeffstring, "%16.9g", b_yptar[icoeff_fit] );
    newcoeffsfile << coeffstring;
    //sprintf( coeffstring, "%16.9g", 0.0 );
    sprintf( coeffstring, "%16.9g", b_delta[icoeff_fit] );
    newcoeffsfile << coeffstring; 
    newcoeffsfile << " ";
    newcoeffsfile << setw(1) << setprecision(1) << xfpexpon[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << yfpexpon[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xtarexpon[icoeff_fit]; 
    newcoeffsfile << endl;
    
  }
  
  TFile hsimc(outputhist,"recreate");
  HList.Write();
  
}
