#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TF1.h>
#include "dc_dead_map.h"
#include "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/svx_cent_ana/source/MyEvent.h"
#include "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/svx_cent_ana/source/Reconstruction.h" 
//#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/MyEvent.h"
//#include "/direct/phenix+u/tongzhouguo/install/include/Run14AuAuLeptonCombyConstants.h"
//#include "/direct/phenix+u/tongzhouguo/install/include/embedana.h"
using namespace std;
using namespace DileptonAnalysis;

const float ncoll[] = {
  313.0,
  129.0,
  41.8,
  10.1,
  3.6   // estimated from MB 0-93 closure
};

const float mee_min = 0.01, mee_max = 0.12;

// Scale factors for sectors 0–7
const double scale[8] = {
    0.988, 0.990, 0.985, 0.981,
    0.980, 0.985, 1.025, 1.023
};

// Smear c1 values for sectors 0–7
const double smear_c1[8] = {
    0.055, 0.066, 0.055, 0.059,
    0.057, 0.056, 0.072, 0.072
};

// Smear c2 values for sectors 0–7
const double smear_c2[8] = {
    0.011, 0.017, 0.011, 0.013,
    0.012, 0.012, 0.063, 0.062
};

Double_t Hagedorn_Yield_Func(Double_t *x, Double_t *p) {
  //
  // Hagedorn function in 1/pt dN/pt from ppg088 
  //
    Double_t pt   = x[0];
    Double_t mass = p[0];
    Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.134*0.134);
    Double_t A    = p[1];
    Double_t a    = p[2];
    Double_t b    = p[3];
    Double_t p0   = p[4];
    Double_t n    = p[5];
   
    Double_t value = 2 * TMath::Pi() * pt *A* pow( exp(-a*mt-b*mt*mt)+mt/p0 , -n);
    return value;
  }
  
TF1 *funcHagedorn(const Char_t *name, Double_t lowerlim,  Double_t upperlim,  const double dNdy = 95.7, const double mass = 0.135, const double A = 504.5, const double a = 0.5169, const double b = 0.1626, const double p0 = 0.7366, const double n = 8.274)
{
  TF1 *funcHagedorn = new TF1(name, Hagedorn_Yield_Func, lowerlim, upperlim, 6);
  funcHagedorn->SetParameters(mass, A, a, b, p0, n); 
  funcHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return funcHagedorn;
}


float EMCMAP[8][48][96];
static const double Me = 0.000510998918;
static const double Me2 = Me*Me;
const int centbin = 5;
const int sectbin = 1;
const int Sectbin = 8;
const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5;
float ecore_cut = 0.5;
float chi2_cut = 3.0;

TH1F* hpT_pi0[centbin];
TH1F* hpT_pi0_unw[centbin];
TH2F* h2d_EMCdead[Sectbin];

double getDcenter(double phi1, double z1, double phi2, double z2){
  double dcenter_phi_sigma = 0.01, dcenter_phi_offset = 0.;
  double dcenter_z_sigma = 3.6, dcenter_z_offset = 0.;

  double dcenter_z = (z1-z2-dcenter_z_offset)/dcenter_z_sigma;
  double dcenter_phi = (phi1-phi2-dcenter_phi_offset)/dcenter_phi_sigma;

  return sqrt(dcenter_phi*dcenter_phi+dcenter_z*dcenter_z);
}

void read_in_emcmap(){
  ifstream readmap("/phenix/plhf/tongzhouguo/pi02gg/all_209_deadmap.txt");//"/phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/test/offline1/AnalysisTrain/DileptonAnalysis/Run14AuAuEmcalDeadMap.txt"
  for (int i = 0; i < 8; ++i){
    for (int j = 0; j < 48; ++j){
      for (int k = 0; k < 96; ++k){
        EMCMAP[i][j][k] = 0;
      }
    }
  }

  int armsect = 0, ypos = 0, zpos = 0, status = 0;
  std::cout<<"reading EMCal dead map: "<<std::endl;
  while (readmap >> armsect >> ypos >> zpos >> status){
      EMCMAP[armsect][ypos][zpos] = status;
  }
  readmap.close();
}

bool apply_emcmap(int sect, int y, int z){
  if (y==0 || z==0) return false;
  if ( sect < 6 && ( y == 35 || z == 71) ) return false;
  if ( sect > 5 && ( y == 47 || z == 95) ) return false;
  if (EMCMAP[sect][y][z]) return false;
  return true;
}


const int cutbin = 3;
const char* cut_name[cutbin] = {"default", "tight", "verytight"};
const double DCENTER_CUT[cutbin] = {10.0, 12.0, 15.0};
const double EP_LO_CUT[cutbin] = {0.6, 0.7, 0.75};
const double EP_HI_CUT[cutbin] = {1.0e9, 1.4, 1.3};
const double N0_CUT[cutbin] = {2.0, 3.0, 3.0};
const double EMCDZ_CUT[cutbin] = {20.0, 15.0, 10.0};
const double EMCDPHI_CUT[cutbin] = {0.05, 0.04, 0.03};

TH2F* eemass[cutbin][centbin][sectbin];
TH2F* eegmass[cutbin][centbin][sectbin];
TH2F* eemass_unw[cutbin][centbin][sectbin];
TH2F* eegmass_unw[cutbin][centbin][sectbin];
TH2F* eemass_pure[cutbin][centbin][sectbin];
TH2F* eegmass_pure[cutbin][centbin][sectbin];

bool single_cut(MyTrack* mytrk, int icut){
  int _side = mytrk->GetDCside();
  int _arm = mytrk->GetArm();
  float _phi = mytrk->GetPhiDC();
  float _alpha = mytrk->GetAlphaPrime(); 
  float _zed = mytrk->GetZDC();

  float _board = -999;
  if(_arm == 0) _board = (3.72402-_phi+0.008047*cos(_phi+0.87851))/0.01963496;
  else if(_arm == 1) _board = (0.573231+_phi-0.0046*cos(_phi+0.05721))/0.01963496;
  if(_board < -99) return false;

  bool dead_area = DCdeadArea(_alpha, _phi, _board, _zed);
  if(dead_area) return false;

  double _pt = mytrk->GetPtPrime();
  if (_pt<0.20) return false;

  double _theta = mytrk->GetThe0Prime();
  double _p = _pt/sin(_theta);

  double _e = mytrk->GetEcore();
  double _ep = _e/_p;
  if(_ep < EP_LO_CUT[icut]) return false;
  if(_ep > EP_HI_CUT[icut]) return false;

  double _n0 = mytrk->GetN0();
  if(_n0 < N0_CUT[icut]) return false;

  int _charge = mytrk->GetCharge();
  if(fabs(_charge) != 1) return false;

  double _emcdz = mytrk->GetEmcdz();
  if(fabs(_emcdz) > EMCDZ_CUT[icut]) return false;

  double _emcdphi = mytrk->GetEmcdphi();
  if(fabs(_emcdphi) > EMCDPHI_CUT[icut]) return false;

  return true;
}

float radius_r;
float phi_positron;
float phi_electron;
float theta_positron;  
float theta_electron;  
bool solution(MyTrack* mytrk1, MyTrack* mytrk2, MyPair* mypair, Reconstruction* reco, float zVtx){
  float DPHI=0.005, R_LO=1, R_HI=29;
  reco->findIntersection(mytrk1, mytrk2, mypair, zVtx);
  radius_r = mypair->GetRPair();
  phi_electron = mypair->GetPhiElectron();
  phi_positron = mypair->GetPhiPositron();
  theta_electron = mypair->GetThetaElectron();
  theta_positron = mypair->GetThetaPositron();
  if ( (radius_r<=R_LO) || (radius_r>=R_HI) ) return false;

  float dphi_r = mypair->GetPhiElectron()-mypair->GetPhiPositron();
  if ( dphi_r>TMath::Pi() )  dphi_r = 2*TMath::Pi()-dphi_r; 
  if ( dphi_r<-TMath::Pi() ) dphi_r = -2*TMath::Pi()-dphi_r; 
  if ( TMath::Abs(dphi_r)>=DPHI ) return false;
  
  return true;
}



void get_pi0_from_gamma_data(const char* inFile = "../tong1.root", const char* outFile = "2d_mass_pt_05t.root", const int system = 0, int ert = 0){ 
  //gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
  //gSystem->Load("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/work/ee/offline/AnalysisTrain/Run14AuAuLeptonComby/lib/libRun14AuAuLeptonEvent.so");
  //gSystem->Load("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/work/ee/offline/AnalysisTrain/Run14AuAuLeptonComby/lib/libRun14AuAuLeptonComby.so");
  //gSystem->Load("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/work/ee/offline/AnalysisTrain/Run14AuAuLeptonComby/lib/libRun14AuAuLeptonConvReco.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonEvent.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libsvxcentana.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonConvReco.so");
  Reconstruction reco("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/analysis/emb/lookup_3D_one_phi.root");


  TF1 f_rand("f","gaus",-5,5);
  f_rand.SetParameter(0,1);
  f_rand.SetParameter(1,0);
  f_rand.SetParameter(2,1);

  read_in_emcmap();
  TFile* input = new TFile(inFile,"READ");
  if(!(input)){
    cout << "no input file" << endl;
    exit(1);
  }

  cout << input << endl;
  
  TTree* T = (TTree*)input->Get("tree");
  TBranch* br = T->GetBranch("MyEvent");
  MyEvent* event = 0;
  br->SetAddress(&event);
  cout << "Input Tree read in" << endl;

  for (int isect = 0; isect < Sectbin; ++isect){
    h2d_EMCdead[isect] = new TH2F(Form("h2d_EMCdead_sect%d",isect), "IY versus IZ", 48, -0.5, 47.5, 96, -0.5, 95.5);
  }

  for(int icent = 0; icent < centbin; icent++){   
    hpT_pi0[icent] = new TH1F(Form("htrue_pi0%d",icent), "", 100, 0.0, 10.0);
    hpT_pi0_unw[icent] = new TH1F(Form("htrue_unw_pi0_wt%d",icent), "", 100, 0.0, 10.0);
    hpT_pi0[icent]->Sumw2();
    hpT_pi0_unw[icent]->Sumw2();
    for (int icut = 0; icut < cutbin; ++icut){
      for (int isect = 0; isect < sectbin; ++isect){
        eemass[icut][icent][isect] = new TH2F(Form("eemass_%s_cent%d_sect%d",cut_name[icut],icent,isect), "ee mass versus ee pt", 100, 0.0, 1.0, 100, 0.0, 10.0);
        eegmass[icut][icent][isect] = new TH2F(Form("eegmass_%s_cent%d_sect%d",cut_name[icut],icent,isect), "pi0 mass versus pi0 pt", 100, 0.0, 1.0, 100, 0.0, 10.0);
        eemass_unw[icut][icent][isect] = new TH2F(Form("eemass_unw_%s_cent%d_sect%d",cut_name[icut],icent,isect), "ee mass versus ee pt", 100, 0.0, 1.0, 100, 0.0, 10.0);
        eegmass_unw[icut][icent][isect] = new TH2F(Form("eegmass_unw_%s_cent%d_sect%d",cut_name[icut],icent,isect), "pi0 mass versus pi0 pt", 100, 0.0, 1.0, 100, 0.0, 10.0);
        eemass_pure[icut][icent][isect] = new TH2F(Form("eemass_pure_%s_cent%d_sect%d",cut_name[icut],icent,isect), "ee mass versus ee pt", 100, 0.0, 1.0, 100, 0.0, 10.0);
        eegmass_pure[icut][icent][isect] = new TH2F(Form("eegmass_pure_%s_cent%d_sect%d",cut_name[icut],icent,isect), "pi0 mass versus pi0 pt", 100, 0.0, 1.0, 100, 0.0, 10.0);
        eemass[icut][icent][isect]->Sumw2();
        eegmass[icut][icent][isect]->Sumw2();
        eemass_unw[icut][icent][isect]->Sumw2();
        eegmass_unw[icut][icent][isect]->Sumw2();
        eemass_pure[icut][icent][isect]->Sumw2();
        eegmass_pure[icut][icent][isect]->Sumw2();
      }
    }
  }

  TF1* fHagedorn[centbin];
  for(int icent = 0; icent < centbin; icent++){
    fHagedorn[icent] = new TF1("fHagedorn","[7]*6.2831853*x*((1/(1+exp((x-[5])/[6])))*[0]/(1+x/[1])**[2] + (1-(1/(1+exp((x-[5])/[6]))))*[3]/x**[4])",0.0,10.0);
    fHagedorn[icent]->SetNpx(10000);
    fHagedorn[icent]->SetParameters(1.922170e+02, 2.106180e+00, 1.284030e+01, 1.630000e+01, 8.060480e+00, 4.033000e+00, 6.395340e-02, 1.0*ncoll[icent]); ////NEED TO UPDATE PARAMETERS
    //fHagedorn[icent]->SetParameters(364.5, 2.106180e+00, 1.284030e+01, 1.630000e+01, 8.261e+00, 4.033000e+00, 6.395340e-02, 1.0*ncoll[icent]); ////NEED TO UPDATE PARAMETERS
    //fHagedorn[icent] = funcHagedorn(Form("fHagedorn_%d", icent), 0.01, 25, 99, 0.139, 364.5*ncoll[icent], 0.433, 0.1221 , 0.7385 , 8.261);
  }

  cout << "Histograms defined" << endl;
  int nevt = T->GetEntries();

  //==============================
  //         Event A Loop
  //==============================
  for (int ievent_A = 0; ievent_A < nevt; ievent_A++){
    if (ievent_A%1000==0) cout << "Event: " << ievent_A << " / " << nevt << endl;

    event->ClearEvent();
    br->GetEntry(ievent_A);

    int ntrack_A = event->GetNtrack();
    int nclust_A = event->GetNcluster();
    int cent = event->GetCentrality();
    double zVtx_A = event->GetPreciseZ();
    int n_gen_track = event->GetNgentrack();

    if ( ntrack_A == 0 ) continue;
    if ( cent < 0 || cent > 93 ) continue;

    int cent_A = cent/(100/centbin);

    float pT_pi0;
    float pT_x = 0;
    float pT_y = 0;
    for (int i_gen_track = 0; i_gen_track < n_gen_track; ++i_gen_track){
      MyGenTrack my_gen_track = *event->GetGenTrack(i_gen_track);
      pT_x += my_gen_track.GetPx();
      pT_y += my_gen_track.GetPy();
    }
    pT_pi0 = sqrt(pT_x*pT_x+pT_y*pT_y);
    float w_pt_tr = fHagedorn[cent_A]->Eval(pT_pi0);
    hpT_pi0[cent_A]->Fill(pT_pi0, w_pt_tr);
    hpT_pi0_unw[cent_A]->Fill(pT_pi0);

    for (int icut = 0; icut < cutbin; ++icut){
      //Track Selection for Event A
      int ghost_A[40], real_A[40];
      int all_A[40];
      int nghost_A = 0, nreal_A = 0, nall_A = 0;

      for (int k1 = 0; k1 < ntrack_A; ++k1){
        if ( !single_cut(event->GetEntry(k1), icut)) continue; // eid cuts
        all_A[nall_A] = k1;
        ++nall_A;        
        int ghost_flag = 0;

        for (int k2 = 0; k2 < ntrack_A; ++k2){
          if ( k1==k2 ) continue;
          if ( !single_cut(event->GetEntry(k2), icut) ) continue;
          // found a RICH ghost pair
          if ( getDcenter((event->GetEntry(k1))->GetCrkphi(), (event->GetEntry(k1))->GetCrkz(), (event->GetEntry(k2))->GetCrkphi(), (event->GetEntry(k2))->GetCrkz()) < DCENTER_CUT[icut] ) { ghost_flag = 1; break; }
        }
        if ( ghost_flag ) { ghost_A[nghost_A] = k1; ++nghost_A; }
        else { real_A[nreal_A] = k1; ++nreal_A; }
      }

      //Track Loop
      for (int ireal_A1 = 0; ireal_A1 < nreal_A; ireal_A1++){ 
        MyTrack mytrk_A1 = *event->GetEntry(real_A[ireal_A1]);
        int arm_A1 = mytrk_A1.GetArm();
        int charge_A1 = mytrk_A1.GetCharge();
        float the_A1 = mytrk_A1.GetThe0Prime();
        float phi_A1 = mytrk_A1.GetPhi0Prime();
        float phidc_A1 = mytrk_A1.GetPhiDC();
        float zed_A1 = mytrk_A1.GetZDC();
        float pt_A1 = mytrk_A1.GetPtPrime();

        float momx_A1 = pt_A1 * cos(phi_A1);
        float momy_A1 = pt_A1 * sin(phi_A1);
        float momz_A1 = pt_A1/tan(the_A1);

        float id_trk1 = mytrk_A1.GetEmcId();

        for (int ireal_A2 = ireal_A1+1; ireal_A2 < nreal_A; ireal_A2++){
          MyTrack mytrk_A2 = *event->GetEntry(real_A[ireal_A2]);

          int arm_A2 = mytrk_A2.GetArm();
          int charge_A2 = mytrk_A2.GetCharge();
          float the_A2 = mytrk_A2.GetThe0Prime();
          float phi_A2 = mytrk_A2.GetPhi0Prime();
          float phidc_A2 = mytrk_A2.GetPhiDC();
          float zed_A2 = mytrk_A2.GetZDC();
          float pt_A2 = mytrk_A2.GetPtPrime();

          float momx_A2 = pt_A2 * cos(phi_A2);
          float momy_A2 = pt_A2 * sin(phi_A2);
          float momz_A2 = pt_A2/tan(the_A2);

          float id_trk2 = mytrk_A2.GetEmcId();

          if(id_trk1==id_trk2) continue;
          if (fabs(phidc_A1-phidc_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) continue;
          if (charge_A1 == charge_A2) continue;
          if (arm_A1 != arm_A2) continue;

          if (fabs(zed_A1-zed_A2)>4.0) continue;

          MyPair mypair_AA;
          //if (!solution(&mytrk_A1, &mytrk_A2, &mypair_AA, &reco, zVtx_A)) continue;

          TLorentzVector p_A1_DC;
          p_A1_DC.Clear();
          p_A1_DC.SetPx(momx_A1);
          p_A1_DC.SetPy(momy_A1);
          p_A1_DC.SetPz(momz_A1);
          p_A1_DC.SetE(sqrt(pow(p_A1_DC.P(),2) + Me2));

          TLorentzVector p_A2_DC;
          p_A2_DC.Clear();
          p_A2_DC.SetPx(momx_A2);
          p_A2_DC.SetPy(momy_A2);
          p_A2_DC.SetPz(momz_A2);
          p_A2_DC.SetE(sqrt(pow(p_A2_DC.P(),2) + Me2));

          //convPhoton_AA: conversion photon with virtual mass; real_convPhoton_AA: conversion photon without virtual mass
          TLorentzVector convPhoton_AA;
          convPhoton_AA.Clear();
          convPhoton_AA = p_A1_DC + p_A2_DC;
          float Pt_AA = convPhoton_AA.Pt();

          TLorentzVector real_convPhoton_AA;
          real_convPhoton_AA.Clear();
          real_convPhoton_AA.SetPx(momx_A1+momx_A2);
          real_convPhoton_AA.SetPy(momy_A1+momy_A2);
          real_convPhoton_AA.SetPz(momz_A1+momz_A2);
          real_convPhoton_AA.SetE(sqrt(pow(convPhoton_AA.P(),2)));
          float Pt_AA_real = real_convPhoton_AA.Pt();

          eemass[icut][cent_A][0]->Fill(convPhoton_AA.M(), Pt_AA, w_pt_tr);
          eemass_unw[icut][cent_A][0]->Fill(convPhoton_AA.M(), Pt_AA);
          if((convPhoton_AA.M()<mee_min)||(convPhoton_AA.M()>mee_max)) continue;

          for (int iclust_A = 0; iclust_A < nclust_A; ++iclust_A){
            MyCluster myclust_A = *event->GetClusterEntry(iclust_A);
            const int is_embd = myclust_A.GetMcId()==1 ? 1 : 0;
            if (is_embd==0) continue;
            float xclust_A = myclust_A.GetX();
            float yclust_A = myclust_A.GetY();
            float zclust_A = myclust_A.GetZ(); 
            float e_A  = myclust_A.GetEcore();
            float chi2_A   = myclust_A.GetChi2();
            //int sector = 7 - myclust_A.GetSect();
            int sector = myclust_A.GetSect();
            const double e_smear = scale[sector] *  ( 1 + f_rand.GetRandom() * sqrt( pow ( smear_c1[sector] , 2 ) + pow( smear_c2[sector] / sqrt(e_A), 2) ) );          // (1.+0.025*f.GetRandom());//// 8.1/sqrt(E)+2.1
            e_A *= e_smear;
            float id_clust = myclust_A.GetEmcId();
            if(id_clust==id_trk1||id_clust==id_trk2) continue;
            //cout<<"IY = "<<myclust_A.GetIY()<<"\t"<<"IZ = "<<myclust_A.GetIZ()<<endl;
            //cluster cuts
            if(!apply_emcmap(sector,myclust_A.GetIY(),myclust_A.GetIZ())) continue; 
            if(e_A < ecore_cut) continue;
            if(chi2_A > chi2_cut) continue;

            TLorentzVector emcPhoton_A;
            emcPhoton_A.Clear();
            emcPhoton_A.SetPx( e_A*xclust_A/sqrt(xclust_A*xclust_A+yclust_A*yclust_A+zclust_A*zclust_A) );
            emcPhoton_A.SetPy( e_A*yclust_A/sqrt(xclust_A*xclust_A+yclust_A*yclust_A+zclust_A*zclust_A) );
            emcPhoton_A.SetPz( e_A*zclust_A/sqrt(xclust_A*xclust_A+yclust_A*yclust_A+zclust_A*zclust_A) );
            emcPhoton_A.SetE( e_A );

            TLorentzVector pi0_AA;
            pi0_AA.Clear();
            pi0_AA = real_convPhoton_AA + emcPhoton_A;

            float mass_pi0 = pi0_AA.M();
            float pt_pi0_reco = pi0_AA.Pt();
            if (id_clust==id_trk1||id_clust==id_trk2||id_trk1==id_trk2) continue;

            //float diff = (pt_pi0_reco-pT_pi0)/pT_pi0;
            //if(diff>0.15||diff<-0.15) continue;

            if (icut == 0) h2d_EMCdead[sector]->Fill(myclust_A.GetIY(), myclust_A.GetIZ());
            eegmass[icut][cent_A][0]->Fill(pi0_AA.M(), Pt_AA, w_pt_tr);
            eegmass_unw[icut][cent_A][0]->Fill(pi0_AA.M(), Pt_AA);
          } 
        }
      }
    }
  }

  for (int icut = 0; icut < cutbin; ++icut) {
    for (int icent = 0; icent < centbin; ++icent) {
      for (int isect = 0; isect < sectbin; ++isect) {
        int nbinsx_ee = eemass[icut][icent][isect]->GetNbinsX();
        int nbinsy_ee = eemass[icut][icent][isect]->GetNbinsY();

        for (int ibin_x = 1; ibin_x <= nbinsx_ee; ++ibin_x) {
          for (int ibin_y = 1; ibin_y <= nbinsy_ee; ++ibin_y) {
            double n_unw = eemass_unw[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            if (n_unw == 0) continue;
            double content = eemass[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            double err = content / sqrt(n_unw);
            eemass[icut][icent][isect]->SetBinError(ibin_x, ibin_y, err);
          }
        }

        int nbinsx_eeg = eegmass[icut][icent][isect]->GetNbinsX();
        int nbinsy_eeg = eegmass[icut][icent][isect]->GetNbinsY();

        for (int ibin_x = 1; ibin_x <= nbinsx_eeg; ++ibin_x) {
          for (int ibin_y = 1; ibin_y <= nbinsy_eeg; ++ibin_y) {
            double n_unw = eegmass_unw[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            if (n_unw == 0) continue;
            double content = eegmass[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            double err = content / sqrt(n_unw);
            eegmass[icut][icent][isect]->SetBinError(ibin_x, ibin_y, err);
          }
        }

        int nbinsx_ee_pure = eemass_pure[icut][icent][isect]->GetNbinsX();
        int nbinsy_ee_pure = eemass_pure[icut][icent][isect]->GetNbinsY();
        for (int ibin_x = 1; ibin_x <= nbinsx_ee_pure; ++ibin_x) {
          for (int ibin_y = 1; ibin_y <= nbinsy_ee_pure; ++ibin_y) {
            double n_unw = eemass_unw[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            if (n_unw == 0) continue;
            double content = eemass_pure[icut][icent][isect]->GetBinContent(ibin_x, ibin_y);
            double err = content / sqrt(n_unw);
            eemass_pure[icut][icent][isect]->SetBinError(ibin_x, ibin_y, err);
          }
        }
      }
    }
  }

  TFile* output = new TFile(outFile,"RECREATE");
  output->cd();
  for(int isect = 0; isect < Sectbin; isect++){
    h2d_EMCdead[isect]->Write();
  } 

  for(int icent = 0; icent < centbin; icent++){
    hpT_pi0[icent]->Write();
    hpT_pi0_unw[icent]->Write();
    for (int icut = 0; icut < cutbin; ++icut){
      for(int isect = 0; isect < sectbin; isect++){      
        eemass[icut][icent][isect]->Write();
        eegmass[icut][icent][isect]->Write();
        eemass_unw[icut][icent][isect]->Write();
        eegmass_unw[icut][icent][isect]->Write();
        eemass_pure[icut][icent][isect]->Write();
        eegmass_pure[icut][icent][isect]->Write();
      }
    }
  }
  
  output->Close();
  input->Close();
}
