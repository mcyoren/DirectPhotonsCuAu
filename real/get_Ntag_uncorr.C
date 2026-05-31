#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TRandom3.h>

#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/MyEvent.h"
#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/Reconstruction.h"
#include "dc_dead_map.h"

using namespace std;
using namespace DileptonAnalysis;

//==============================================================
// Systematic version of the original real-data macro.
//
// This keeps the original input structure:
//   - function name: get_Ntag_uncorr
//   - tree name:     T
//   - branch name:   MyEvent
//   - track vars:    GetPx/GetPy/GetPz, GetThe0/GetPhi0/GetAlpha
//   - cluster sector convention: 7 - GetSect()
//   - BG definition: event-mixed photons, as in original code
//
// It only generalizes the cuts into configs and writes many histograms.
//==============================================================

static const double PI = 3.14159;
static const double Me = 0.000510998918;
static const double Me2 = Me*Me;
static const double c = 299792458;

const double DCENTERCUT = 10; // RICH ghost cut
const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5; // PC1 ghost cut

// Keep original centrality structure from this real-data macro.
const int centbin = 4;
int cent_low[centbin] = {0, 20, 40, 60};
int cent_hi[centbin]  = {20, 40, 60, 93};
double n[centbin] = {0};

// Baseline event cuts from original macro.
double bbcz_cut_val = 10;
const double bbcq_cut_val = 1800;

// Original default values are kept for reference.
// The actual systematic cut values are in cfgs[] below.
float ecore_cut = 0.4;
float chi2_cut = 3.0;

float EMCMAP[8][48][96];
const double POOLD = 10; // pool depth for ee+gamma mixed photon BG
const double POOLD_EE = 3; // smaller pool depth for mixed-event ee BG

// Histogram binning requested for systematic output.
const int NMASS = 100;
const int NPT = 100;
const double MASS_MIN = 0.0;
const double MASS_MAX = 1.0;
const double PT_MIN = 0.0;
const double PT_MAX = 10.0;

const double MEE_MIN_FOR_EEG = 0.01;
const double MEE_MAX_FOR_EEG = 0.12;

// Main dzed cut follows original code behavior: dzed < 4.
const double DZED_MAIN_CUT = 4.0;

// "Tight" conversion cut: keep original solution(), then add tighter R and theta.
// If your screenshot uses different values, edit these three constants only.
const double TIGHT_RCONV_MIN = 9.0;
const double TIGHT_RCONV_MAX = 21.0;
const double TIGHT_DTHETA_CUT = 0.01;

// Systematic configuration.
enum ConvType
{
  CONV_SOLUTION = 0, // solution() only
  CONV_DZED     = 1, // solution() + dzed < DZED_MAIN_CUT
  CONV_TIGHT    = 2  // solution() + dzed + tight R/dtheta
};

struct SystConfig
{
  const char* name;
  int convType;
  double ecoreCut;
  double ptCut;
  int eidType;
  double chi2Cut;
  double escale;
};

// Nominal choice:
//   conv:  solution + dzed < 4
//   ecore: 0.4
//   pt_e:  0.3
//   eid2:  n0>2, e/p>0.7, disp<4
//   chi2:  3
//   escale: 1.00
const int NCFG = 14;
SystConfig cfgs[NCFG] = {
  {"nominal",       CONV_DZED,     0.4, 0.3, 2, 3.0, 1.00},
  {"conv_solution", CONV_SOLUTION, 0.4, 0.3, 2, 3.0, 1.00},
  {"conv_tight",    CONV_TIGHT,    0.4, 0.3, 2, 3.0, 1.00},
  {"ecore03",       CONV_DZED,     0.3, 0.3, 2, 3.0, 1.00},
  {"ecore05",       CONV_DZED,     0.5, 0.3, 2, 3.0, 1.00},
  {"pte02",         CONV_DZED,     0.4, 0.2, 2, 3.0, 1.00},
  {"pte04",         CONV_DZED,     0.4, 0.4, 2, 3.0, 1.00},
  {"eid0",          CONV_DZED,     0.4, 0.3, 0, 3.0, 1.00},
  {"eid1",          CONV_DZED,     0.4, 0.3, 1, 3.0, 1.00},
  {"eid3",          CONV_DZED,     0.4, 0.3, 3, 3.0, 1.00},
  {"chi2_4",        CONV_DZED,     0.4, 0.3, 2, 4.0, 1.00},
  {"chi2_5",        CONV_DZED,     0.4, 0.3, 2, 5.0, 1.00},
  {"escale099",     CONV_DZED,     0.4, 0.3, 2, 3.0, 0.99},
  {"escale101",     CONV_DZED,     0.4, 0.3, 2, 3.0, 1.01}
};

// Original-style objects.
TH1F* h_nevents = 0;

TH2D* h_cent_temp = 0;
TH2D* h_cent = 0;
TH1D* h_n_cent = 0;
TH1D* h_nev_cent = 0;

// Systematic histograms.
// FG = same-event unlike-sign conversion ee + same-event photon.
// BG = same accepted conversion ee + mixed-event photon, preserving original macro definition.
TH2F* fg2d_eemass_reco[NCFG][centbin];
TH2F* bg2d_eemass_reco[NCFG][centbin]; // mixed-event ee background
TH2F* fg2d_eegmass_reco[NCFG][centbin];
TH2F* bg2d_eegmass_reco[NCFG][centbin];

TH1F* h_cutflow[NCFG];

// Keep old nominal names too, for quick comparison with old scripts.
TH2F* fg2d_eemass_reco_dzed[centbin];
TH2F* fg2d_eemass_reco_dphi[centbin];
TH2F* bg2d_eemass_reco_dzed[centbin];
TH2F* bg2d_eemass_reco_dphi[centbin];
TH2F* fg2d_eegmass_reco_dzed[centbin];
TH2F* fg2d_eegmass_reco_dphi[centbin];
TH2F* bg2d_eegmass_reco_dzed[centbin];
TH2F* bg2d_eegmass_reco_dphi[centbin];

double getDcenter(double phi1, double z1, double phi2, double z2)
{
  double dcenter_phi_sigma = 0.01, dcenter_phi_offset = 0.;
  double dcenter_z_sigma = 3.6, dcenter_z_offset = 0.;

  double dcenter_z = (z1-z2-dcenter_z_offset)/dcenter_z_sigma;
  double dcenter_phi = (phi1-phi2-dcenter_phi_offset)/dcenter_phi_sigma;

  return sqrt(dcenter_phi*dcenter_phi+dcenter_z*dcenter_z);
}

void read_in_emcmap()
{
  ifstream readmap("/phenix/plhf/tongzhouguo/pi02gg/209_deadmap.txt");
  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 48; ++j)
    {
      for (int k = 0; k < 96; ++k)
      {
        EMCMAP[i][j][k] = 0;
      }
    }
  }

  int armsect = 0, ypos = 0, zpos = 0, status = 0;
  std::cout << "reading EMCal dead map: " << std::endl;
  while (readmap >> armsect >> ypos >> zpos >> status)
  {
    if (armsect >= 0 && armsect < 8 &&
        ypos >= 0 && ypos < 48 &&
        zpos >= 0 && zpos < 96)
    {
      EMCMAP[armsect][ypos][zpos] = status;
    }
  }
  readmap.close();
}

bool apply_emcmap(int sect, int y, int z)
{
  if (sect < 0 || sect >= 8) return false;
  if (y < 0 || y >= 48) return false;
  if (z < 0 || z >= 96) return false;

  if (y==0 || z==0) return false;
  if (sect < 6 && (y == 35 || z == 71)) return false;
  if (sect > 5 && (y == 47 || z == 95)) return false;
  if (EMCMAP[sect][y][z]) return false;
  return true;
}

// Base track quality from original single_cut, but without pt/eid.
// Config-specific pt/eid are applied in pass_eid().
bool pass_base_track(MyTrack* mytrk)
{
  if (!mytrk) return false;

  int _arm = mytrk->GetArm();
  float _phi = mytrk->GetPhiDC();
  float _alpha = mytrk->GetAlpha();
  float _zed = mytrk->GetZDC();

  float _board = -999;
  if (_arm == 0) _board = (3.72402-_phi+0.008047*cos(_phi+0.87851))/0.01963496;
  else if (_arm == 1) _board = (0.573231+_phi-0.0046*cos(_phi+0.05721))/0.01963496;
  if (_board < -99) return false;

  bool dead_area = DCdeadArea(_alpha, _phi, _board, _zed);
  if (dead_area) return false;

  int _charge = mytrk->GetCharge();
  if (fabs((double)_charge) != 1.0) return false;

  double _emcdz = mytrk->GetEmcdz();
  if (fabs(_emcdz) > 20) return false;

  double _emcdphi = mytrk->GetEmcdphi();
  if (fabs(_emcdphi) > 0.05) return false;

  return true;
}

double track_p(MyTrack* mytrk)
{
  return sqrt(mytrk->GetPx()*mytrk->GetPx() +
              mytrk->GetPy()*mytrk->GetPy() +
              mytrk->GetPz()*mytrk->GetPz());
}

double track_pt(MyTrack* mytrk)
{
  return track_p(mytrk) * sin(mytrk->GetThe0());
}

// EID definitions requested.
// eid0: n0>=2 and e/p<0.5
// eid1: n0>=2 and e/p>0.6 and disp<5
// eid2: n0>2 and e/p>0.7 and disp<4  (nominal)
// eid3: e/p>0.8 and n0>3 and disp<4 and prob>0.01
bool pass_eid(MyTrack* mytrk, int eidType, double ptCut)
{
  if (!mytrk) return false;
  if (!pass_base_track(mytrk)) return false;

  double p = track_p(mytrk);
  if (p <= 0.0) return false;

  double pt = track_pt(mytrk);
  if (pt < ptCut) return false;

  double ep = mytrk->GetEcore() / p;
  double n0 = mytrk->GetN0();
  double disp = mytrk->GetDisp();
  double prob = mytrk->GetProb();

  if (eidType == 0)
  {
    return (n0 >= 2.0 && ep < 0.5);
  }
  else if (eidType == 1)
  {
    return (n0 >= 2.0 && ep > 0.6 && disp < 5.0);
  }
  else if (eidType == 2)
  {
    return (n0 > 2.0 && ep > 0.7 && disp < 4.0);
  }
  else if (eidType == 3)
  {
    return (ep > 0.8 && n0 > 3.0 && disp < 4.0 && prob > 0.01);
  }

  return false;
}

// Original single_cut kept for reference and for possible old checks.
bool single_cut(MyTrack* mytrk)
{
  return pass_eid(mytrk, 1, 0.20); // old: pt>0.2, e/p>0.6, n0>=2, plus base quality. Adds disp<5 only in new cfg; not used for nominal logic.
}

float radius_r;
float phi_positron;
float phi_electron;
float theta_positron;
float theta_electron;
bool solution(MyTrack* mytrk1, MyTrack* mytrk2, MyPair* mypair, Reconstruction* reco, float zVtx)
{
  float DPHI=0.005, R_LO=1, R_HI=29;
  reco->findIntersection(mytrk1, mytrk2, mypair, zVtx);
  radius_r = mypair->GetRPair();
  phi_electron = mypair->GetPhiElectron();
  phi_positron = mypair->GetPhiPositron();
  theta_electron = mypair->GetThetaElectron();
  theta_positron = mypair->GetThetaPositron();
  if ((radius_r<=R_LO) || (radius_r>=R_HI)) return false;

  float dphi_r = mypair->GetPhiElectron()-mypair->GetPhiPositron();
  if (dphi_r>TMath::Pi())  dphi_r = 2*TMath::Pi()-dphi_r;
  if (dphi_r<-TMath::Pi()) dphi_r = -2*TMath::Pi()-dphi_r;
  if (TMath::Abs(dphi_r)>=DPHI) return false;

  return true;
}

bool pass_conv_config(int icfg, double dzed)
{
  if (cfgs[icfg].convType == CONV_SOLUTION)
  {
    return true;
  }

  if (cfgs[icfg].convType == CONV_DZED)
  {
    return (dzed < DZED_MAIN_CUT);
  }

  if (cfgs[icfg].convType == CONV_TIGHT)
  {
    if (dzed >= DZED_MAIN_CUT) return false;
    if (radius_r <= TIGHT_RCONV_MIN || radius_r >= TIGHT_RCONV_MAX) return false;
    if (fabs(theta_electron - theta_positron) >= TIGHT_DTHETA_CUT) return false;
    return true;
  }

  return false;
}

int get_cent_bin(float cent_A)
{
  int icent = -999;
  if (cent_A > 0 && cent_A <= 20) icent = 0;
  else if (cent_A > 20 && cent_A <= 40) icent = 1;
  else if (cent_A > 40 && cent_A <= 60) icent = 2;
  else if (cent_A > 60 && cent_A <= 93) icent = 3;
  return icent;
}

void book_histograms()
{
  h_nevents = new TH1F("h_nevents","Event counter", 20, -0.5, 19.5);

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    h_cutflow[icfg] = new TH1F(Form("h_cutflow_%s", cfgs[icfg].name),
                               "cutflow", 12, -0.5, 11.5);
  }

  for (int icent = 0; icent < centbin; icent++)
  {
    // Old-like names for quick comparison with old output.
    fg2d_eemass_reco_dzed[icent] = new TH2F(Form("fg2d_eemass_reco_dzed_icent%d",icent), "eemass fg dzed nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    fg2d_eemass_reco_dphi[icent] = new TH2F(Form("fg2d_eemass_reco_dphi_icent%d",icent), "eemass fg solution nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    bg2d_eemass_reco_dzed[icent] = new TH2F(Form("bg2d_eemass_reco_dzed_icent%d",icent), "eemass mixed-event bg dzed nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    bg2d_eemass_reco_dphi[icent] = new TH2F(Form("bg2d_eemass_reco_dphi_icent%d",icent), "eemass mixed-event bg solution nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    fg2d_eegmass_reco_dzed[icent] = new TH2F(Form("fg2d_eegmass_reco_dzed_icent%d",icent), "eegmass fg dzed nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    fg2d_eegmass_reco_dphi[icent] = new TH2F(Form("fg2d_eegmass_reco_dphi_icent%d",icent), "eegmass fg solution nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    bg2d_eegmass_reco_dzed[icent] = new TH2F(Form("bg2d_eegmass_reco_dzed_icent%d",icent), "eegmass bg dzed nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
    bg2d_eegmass_reco_dphi[icent] = new TH2F(Form("bg2d_eegmass_reco_dphi_icent%d",icent), "eegmass bg solution nominal", NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

    fg2d_eemass_reco_dzed[icent]->Sumw2();
    fg2d_eemass_reco_dphi[icent]->Sumw2();
    bg2d_eemass_reco_dzed[icent]->Sumw2();
    bg2d_eemass_reco_dphi[icent]->Sumw2();
    fg2d_eegmass_reco_dzed[icent]->Sumw2();
    fg2d_eegmass_reco_dphi[icent]->Sumw2();
    bg2d_eegmass_reco_dzed[icent]->Sumw2();
    bg2d_eegmass_reco_dphi[icent]->Sumw2();

    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      fg2d_eemass_reco[icfg][icent] = new TH2F(Form("fg2d_eemass_reco_%s_icent%d", cfgs[icfg].name, icent),
                                                "eemass fg;M_{ee};p_{T,ee}",
                                                NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
      bg2d_eemass_reco[icfg][icent] = new TH2F(Form("bg2d_eemass_reco_%s_icent%d", cfgs[icfg].name, icent),
                                                "eemass mixed-event bg;M_{ee};p_{T,ee}",
                                                NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
      fg2d_eegmass_reco[icfg][icent] = new TH2F(Form("fg2d_eegmass_reco_%s_icent%d", cfgs[icfg].name, icent),
                                                 "eegmass fg;M_{ee#gamma};p_{T,ee#gamma}",
                                                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);
      bg2d_eegmass_reco[icfg][icent] = new TH2F(Form("bg2d_eegmass_reco_%s_icent%d", cfgs[icfg].name, icent),
                                                 "eegmass bg mixed;M_{ee#gamma};p_{T,ee#gamma}",
                                                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eemass_reco[icfg][icent]->Sumw2();
      bg2d_eemass_reco[icfg][icent]->Sumw2();
      fg2d_eegmass_reco[icfg][icent]->Sumw2();
      bg2d_eegmass_reco[icfg][icent]->Sumw2();
    }
  }
}

void write_histograms(TFile* output)
{
  output->cd();

  h_nevents->Write();

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    h_cutflow[icfg]->Write();
  }

  for (int icent = 0; icent < centbin; icent++)
  {
    fg2d_eemass_reco_dphi[icent]->Write();
    fg2d_eemass_reco_dzed[icent]->Write();
    bg2d_eemass_reco_dphi[icent]->Write();
    bg2d_eemass_reco_dzed[icent]->Write();
    fg2d_eegmass_reco_dphi[icent]->Write();
    fg2d_eegmass_reco_dzed[icent]->Write();
    bg2d_eegmass_reco_dphi[icent]->Write();
    bg2d_eegmass_reco_dzed[icent]->Write();

    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      fg2d_eemass_reco[icfg][icent]->Write();
      bg2d_eemass_reco[icfg][icent]->Write();
      fg2d_eegmass_reco[icfg][icent]->Write();
      bg2d_eegmass_reco[icfg][icent]->Write();
    }
  }

  if (h_cent) h_cent->Write();
  if (h_n_cent) h_n_cent->Write();
  if (h_nev_cent) h_nev_cent->Write();
}

void get_Ntag_uncorr(const char* inFile = "/phenix/plhf/tongzhouguo/taxi/Run12CuAu200MinBias/19325/data/372403.root",
                     const char* outFile = "test.root",
                     const char* runno = "372403",
                     const int system = 2,
                     const int mode = 0)
{ // start make_TTree
  // system = 0 for pp; system = 1 for pAu
  // mode = 0 for data; mode = 1 for simulation; mode = 2 for embedding

  TH1::AddDirectory(kFALSE);

  // Keep exactly the libraries from the correct real-data macro.
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");
  Reconstruction reco("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Tongzhou/real/lookup_3D_one_phi.root");

  read_in_emcmap();

  TString inFile_A = inFile;
  TString inFile_B = inFile;

  TFile* input_A = TFile::Open(inFile_A, "READ");
  if (!input_A || input_A->IsZombie())
  {
    cout << "no input_A file: " << inFile_A << endl;
    if (input_A) { input_A->Close(); delete input_A; }
    return;
  }

  TFile* input_B = TFile::Open(inFile_B, "READ");
  if (!input_B || input_B->IsZombie())
  {
    cout << "no input_B file: " << inFile_B << endl;
    if (input_B) { input_B->Close(); delete input_B; }
    input_A->Close();
    delete input_A;
    return;
  }

  h_cent_temp = (TH2D*)input_A->Get("hcent");
  if (h_cent_temp)
  {
    h_cent = (TH2D*)h_cent_temp->Clone("h_cent");
    h_cent->Reset("ICESM");
    h_cent->Add(h_cent_temp);
    h_cent->GetYaxis()->SetRangeUser(-10.01, 10.01);
    h_n_cent = (TH1D*)h_cent->ProjectionX("h_n_cent",200,300);
    h_n_cent->SetName("h_n_cent");

    h_nev_cent = new TH1D("h_nev_cent", "Centrality events", 11, -0.5, 10.5);
    for (int icent = 0; icent < centbin; ++icent)
    {
      n[icent] = 0.0;
      for (int ibin = cent_low[icent]+2; ibin < cent_hi[icent]+2; ++ibin)
      {
        n[icent] += h_n_cent->GetBinContent(ibin);
      }
      h_nev_cent->SetBinContent(icent+1, n[icent]);
    }
  }
  else
  {
    cout << "Warning: hcent not found. Continuing without hcent normalization histograms." << endl;
  }

  TTree* T_A = (TTree*)input_A->Get("T");
  if (!T_A)
  {
    cout << "Cannot find tree T in input_A" << endl;
    input_A->Close(); delete input_A;
    input_B->Close(); delete input_B;
    return;
  }

  TBranch* br_A = T_A->GetBranch("MyEvent");
  if (!br_A)
  {
    cout << "Cannot find branch MyEvent in T_A" << endl;
    input_A->Close(); delete input_A;
    input_B->Close(); delete input_B;
    return;
  }

  MyEvent* event_A = 0;
  br_A->SetAddress(&event_A);

  TTree* T_B = (TTree*)input_B->Get("T");
  if (!T_B)
  {
    cout << "Cannot find tree T in input_B" << endl;
    input_A->Close(); delete input_A;
    input_B->Close(); delete input_B;
    return;
  }

  TBranch* br_B = T_B->GetBranch("MyEvent");
  if (!br_B)
  {
    cout << "Cannot find branch MyEvent in T_B" << endl;
    input_A->Close(); delete input_A;
    input_B->Close(); delete input_B;
    return;
  }

  MyEvent* event_B = 0;
  br_B->SetAddress(&event_B);

  cout << "Trees read!" << endl;

  TFile* output = TFile::Open(outFile, "RECREATE");
  if (!output || output->IsZombie())
  {
    cout << "Cannot open output file: " << outFile << endl;
    if (output) { output->Close(); delete output; }
    input_A->Close(); delete input_A;
    input_B->Close(); delete input_B;
    return;
  }

  book_histograms();

  int nevt = T_A->GetEntries();
  float w_pt = 1;

  //==============================
  //         Event A Loop
  //==============================
  for (int ievent_A = 0; ievent_A < nevt; ievent_A++)
  {
    if (ievent_A%1000==0) cout << "Event: " << ievent_A << " / " << nevt << endl;
    ///if (ievent_A>10000) break; // for quick checks; remove for full processing

    h_nevents->Fill(0);
    if (event_A) event_A->ClearEvent();
    br_A->GetEntry(ievent_A);
    if (!event_A) continue;
    h_nevents->Fill(1);

    int ntrack_A = event_A->GetNtrack();
    int nclust_A = event_A->GetNcluster();
    double zVtx_A = event_A->GetVtxZ();
    float cent_A = event_A->GetCentrality();
    float bbcq_A = event_A->GetBBCcharge();

    if (cent_A<0 || cent_A>93) continue;
    if (ntrack_A == 0) continue;
    h_nevents->Fill(2);
    if (TMath::Abs(zVtx_A)>bbcz_cut_val) continue;
    h_nevents->Fill(3);
    if (bbcq_A>bbcq_cut_val) continue;
    h_nevents->Fill(4);

    int icent = get_cent_bin(cent_A);
    if (icent < 0 || icent >= centbin) continue;

    //===========================================
    //       Track ghost selection for each config
    //===========================================
    std::vector< std::vector<int> > real_tracks(NCFG);

    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      for (int k1 = 0; k1 < ntrack_A; ++k1)
      {
        MyTrack trk1_tmp = event_A->GetEntry(k1);
        if (!pass_eid(&trk1_tmp, cfgs[icfg].eidType, cfgs[icfg].ptCut)) continue;

        int ghost_flag = 0;
        for (int k2 = 0; k2 < ntrack_A; ++k2)
        {
          if (k1 == k2) continue;
          MyTrack trk2_tmp = event_A->GetEntry(k2);
          if (!pass_eid(&trk2_tmp, cfgs[icfg].eidType, cfgs[icfg].ptCut)) continue;

          if (getDcenter(trk1_tmp.GetCrkphi(), trk1_tmp.GetCrkz(),
                         trk2_tmp.GetCrkphi(), trk2_tmp.GetCrkz()) < DCENTERCUT)
          {
            ghost_flag = 1;
            break;
          }
        }

        if (!ghost_flag) real_tracks[icfg].push_back(k1);
      }
    }

    // Count events with at least one nominal real track.
    if (real_tracks[0].size() == 0) continue;
    h_nevents->Fill(5);

    //==================================================
    // Mixed-event ee background.
    // This is intentionally much cheaper than the eeg photon mixing:
    // use only POOLD_EE matched B events per A event.
    // We do not apply the conversion solution to mixed-event ee pairs,
    // because the two tracks come from different events. This histogram is
    // intended as a combinatorial ee mass background diagnostic/systematic.
    //==================================================
    if (mode == 0)
    {
      int poold_ee = 0;
      for (int ievent_B_ee = ievent_A + 10; ievent_B_ee < ievent_A + 10000; ++ievent_B_ee)
      {
        if (poold_ee >= POOLD_EE) break;
        if (ievent_B_ee >= nevt) break;

        if (event_B) event_B->ClearEvent();
        br_B->GetEntry(ievent_B_ee);
        if (!event_B) continue;

        int ntrack_B_ee = event_B->GetNtrack();
        double zVtx_B_ee = event_B->GetVtxZ();
        float cent_B_ee = event_B->GetCentrality();
        float bbcq_B_ee = event_B->GetBBCcharge();

        if (cent_B_ee < 0 || cent_B_ee > 93) continue;
        if (ntrack_B_ee == 0) continue;
        if (bbcq_B_ee > bbcq_cut_val) continue;
        if (TMath::Abs(zVtx_B_ee) > bbcz_cut_val) continue;
        if (fabs(cent_A - cent_B_ee) > 5.0) continue;
        if (fabs(zVtx_A - zVtx_B_ee) > 2.5) continue;

        poold_ee++;

        for (int icfg = 0; icfg < NCFG; ++icfg)
        {
          if (real_tracks[icfg].size() == 0) continue;

          std::vector<int> real_tracks_B;
          real_tracks_B.reserve(ntrack_B_ee);

          for (int kB = 0; kB < ntrack_B_ee; ++kB)
          {
            MyTrack trkB_tmp = event_B->GetEntry(kB);
            if (!pass_eid(&trkB_tmp, cfgs[icfg].eidType, cfgs[icfg].ptCut)) continue;

            int ghost_flag_B = 0;
            for (int jB = 0; jB < ntrack_B_ee; ++jB)
            {
              if (kB == jB) continue;
              MyTrack trkB2_tmp = event_B->GetEntry(jB);
              if (!pass_eid(&trkB2_tmp, cfgs[icfg].eidType, cfgs[icfg].ptCut)) continue;

              if (getDcenter(trkB_tmp.GetCrkphi(), trkB_tmp.GetCrkz(),
                             trkB2_tmp.GetCrkphi(), trkB2_tmp.GetCrkz()) < DCENTERCUT)
              {
                ghost_flag_B = 1;
                break;
              }
            }

            if (!ghost_flag_B) real_tracks_B.push_back(kB);
          }

          if (real_tracks_B.size() == 0) continue;

          for (unsigned int ia = 0; ia < real_tracks[icfg].size(); ++ia)
          {
            MyTrack trkA = event_A->GetEntry(real_tracks[icfg][ia]);

            const int armA = trkA.GetArm();
            const int chargeA = trkA.GetCharge();
            const double momA = track_p(&trkA);
            const double theA = trkA.GetThe0();
            const double phiA = trkA.GetPhi0();

            const double pxA = momA * sin(theA) * cos(phiA);
            const double pyA = momA * sin(theA) * sin(phiA);
            const double pzA = momA * cos(theA);

            TLorentzVector pA;
            pA.SetPxPyPzE(pxA, pyA, pzA, sqrt(pxA*pxA + pyA*pyA + pzA*pzA + Me2));

            for (unsigned int ib = 0; ib < real_tracks_B.size(); ++ib)
            {
              MyTrack trkB = event_B->GetEntry(real_tracks_B[ib]);

              const int armB = trkB.GetArm();
              const int chargeB = trkB.GetCharge();

              if (armA != armB) continue;
              if (chargeA * chargeB >= 0) continue; // unlike-sign mixed ee

              const double momB = track_p(&trkB);
              const double theB = trkB.GetThe0();
              const double phiB = trkB.GetPhi0();

              const double pxB = momB * sin(theB) * cos(phiB);
              const double pyB = momB * sin(theB) * sin(phiB);
              const double pzB = momB * cos(theB);

              TLorentzVector pB;
              pB.SetPxPyPzE(pxB, pyB, pzB, sqrt(pxB*pxB + pyB*pyB + pzB*pzB + Me2));

              TLorentzVector pairAB = pA + pB;
              bg2d_eemass_reco[icfg][icent]->Fill(pairAB.M(), pairAB.Pt(), w_pt);

              if (icfg == 1)
              {
                bg2d_eemass_reco_dphi[icent]->Fill(pairAB.M(), pairAB.Pt(), w_pt);
              }
              if (icfg == 0)
              {
                bg2d_eemass_reco_dzed[icent]->Fill(pairAB.M(), pairAB.Pt(), w_pt);
              }
            }
          }
        }
      }
    }

    //======================================
    //      Config loops
    //======================================
    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      if (real_tracks[icfg].size() == 0) continue;

      for (unsigned int ireal_A1 = 0; ireal_A1 < real_tracks[icfg].size(); ++ireal_A1)
      {
        MyTrack mytrk_A1 = event_A->GetEntry(real_tracks[icfg][ireal_A1]);

        float mom_A1 = sqrt(mytrk_A1.GetPx()*mytrk_A1.GetPx() +
                            mytrk_A1.GetPy()*mytrk_A1.GetPy() +
                            mytrk_A1.GetPz()*mytrk_A1.GetPz());

        int arm_A1 = mytrk_A1.GetArm();
        int charge_A1 = mytrk_A1.GetCharge();
        int emcid_A1 = mytrk_A1.GetEmcId();
        float the_A1 = mytrk_A1.GetThe0();
        float phi_A1 = mytrk_A1.GetPhi0();
        float phidc_A1 = mytrk_A1.GetPhiDC();
        float zed_A1 = mytrk_A1.GetZDC();

        float momx_A1 = mom_A1 * sin(the_A1) * cos(phi_A1);
        float momy_A1 = mom_A1 * sin(the_A1) * sin(phi_A1);
        float momz_A1 = mom_A1 * cos(the_A1);

        for (unsigned int ireal_A2 = ireal_A1+1; ireal_A2 < real_tracks[icfg].size(); ++ireal_A2)
        {
          MyTrack mytrk_A2 = event_A->GetEntry(real_tracks[icfg][ireal_A2]);

          float mom_A2 = sqrt(mytrk_A2.GetPx()*mytrk_A2.GetPx() +
                              mytrk_A2.GetPy()*mytrk_A2.GetPy() +
                              mytrk_A2.GetPz()*mytrk_A2.GetPz());

          int arm_A2 = mytrk_A2.GetArm();
          int charge_A2 = mytrk_A2.GetCharge();
          int emcid_A2 = mytrk_A2.GetEmcId();
          float the_A2 = mytrk_A2.GetThe0();
          float phi_A2 = mytrk_A2.GetPhi0();
          float phidc_A2 = mytrk_A2.GetPhiDC();
          float zed_A2 = mytrk_A2.GetZDC();

          float momx_A2 = mom_A2 * sin(the_A2) * cos(phi_A2);
          float momy_A2 = mom_A2 * sin(the_A2) * sin(phi_A2);
          float momz_A2 = mom_A2 * cos(the_A2);

          float dzed = fabs(zed_A1-zed_A2);

          h_cutflow[icfg]->Fill(0);

          if (fabs(phidc_A1-phidc_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) continue;
          if (charge_A1 == charge_A2) continue;
          if (arm_A1 != arm_A2) continue;
          if (emcid_A1 == emcid_A2) continue;
          h_cutflow[icfg]->Fill(1);

          MyPair mypair_AA;
          if (!solution(&mytrk_A1, &mytrk_A2, &mypair_AA, &reco, zVtx_A)) continue;
          h_cutflow[icfg]->Fill(2);

          if (!pass_conv_config(icfg, dzed)) continue;
          h_cutflow[icfg]->Fill(3);

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

          TLorentzVector convPhoton_AA;
          convPhoton_AA.Clear();
          convPhoton_AA = p_A1_DC + p_A2_DC;

          TLorentzVector real_convPhoton_AA;
          real_convPhoton_AA.Clear();
          real_convPhoton_AA.SetPx(momx_A1+momx_A2);
          real_convPhoton_AA.SetPy(momy_A1+momy_A2);
          real_convPhoton_AA.SetPz(momz_A1+momz_A2);
          real_convPhoton_AA.SetE(sqrt(pow(convPhoton_AA.P(),2)));

          const float Pt_AA = real_convPhoton_AA.Pt();

          // Fill ee denominator.
          fg2d_eemass_reco[icfg][icent]->Fill(convPhoton_AA.M(), convPhoton_AA.Pt(), w_pt);

          // Old-like outputs for nominal/solution configs.
          if (icfg == 1)
          {
            fg2d_eemass_reco_dphi[icent]->Fill(convPhoton_AA.M(), convPhoton_AA.Pt(), w_pt);
          }
          if (icfg == 0)
          {
            fg2d_eemass_reco_dzed[icent]->Fill(convPhoton_AA.M(), convPhoton_AA.Pt(), w_pt);
          }

          if ((convPhoton_AA.M()<MEE_MIN_FOR_EEG) || (convPhoton_AA.M()>MEE_MAX_FOR_EEG)) continue;
          h_cutflow[icfg]->Fill(4);

          //==================================================
          // combine with photons from same event: FG
          //==================================================
          for (int iclust_A = 0; iclust_A < nclust_A; ++iclust_A)
          {
            MyCluster myclust_A = event_A->GetClusterEntry(iclust_A);
            float xclust_A = myclust_A.GetX();
            float yclust_A = myclust_A.GetY();
            float zclust_A = myclust_A.GetZ()-zVtx_A; // original real-data behavior
            float e_A = myclust_A.GetEcore() * cfgs[icfg].escale;
            float chi2_A = myclust_A.GetChi2();

            int sector_A = 7 - myclust_A.GetSect();
            if (!apply_emcmap(sector_A, myclust_A.GetIY(), myclust_A.GetIZ())) continue;
            if (e_A < cfgs[icfg].ecoreCut) continue;
            if (chi2_A > cfgs[icfg].chi2Cut) continue;
            if (myclust_A.GetID() == emcid_A1 || myclust_A.GetID() == emcid_A2) continue;
            h_cutflow[icfg]->Fill(5);

            double r_A = sqrt(xclust_A*xclust_A + yclust_A*yclust_A + zclust_A*zclust_A);
            if (r_A <= 0.0) continue;

            TLorentzVector emcPhoton_A;
            emcPhoton_A.Clear();
            emcPhoton_A.SetPx(e_A*xclust_A/r_A);
            emcPhoton_A.SetPy(e_A*yclust_A/r_A);
            emcPhoton_A.SetPz(e_A*zclust_A/r_A);
            emcPhoton_A.SetE(e_A);

            TLorentzVector pi0_AA;
            pi0_AA.Clear();
            pi0_AA = real_convPhoton_AA + emcPhoton_A;

            float Mass_AA = pi0_AA.M();
            float Pt_pi0_AA = pi0_AA.Pt();

            fg2d_eegmass_reco[icfg][icent]->Fill(Mass_AA, Pt_pi0_AA, w_pt);

            if (icfg == 1)
            {
              fg2d_eegmass_reco_dphi[icent]->Fill(Mass_AA, Pt_pi0_AA, w_pt);
            }
            if (icfg == 0)
            {
              fg2d_eegmass_reco_dzed[icent]->Fill(Mass_AA, Pt_pi0_AA, w_pt);
            }
          } // End Photon loop Event A

          //==================================================
          // combine with photons from different event: BG
          //==================================================
          if (mode > 0) continue;
          if (ievent_A > 100000) continue;

          int poold = 0;

          for (int ievent_B = ievent_A+10; ievent_B < ievent_A + 10000; ievent_B++)
          {
            if (poold >= POOLD) break;
            if (ievent_B >= nevt) break;

            if (event_B) event_B->ClearEvent();
            br_B->GetEntry(ievent_B);
            if (!event_B) continue;

            int nclust_B = event_B->GetNcluster();
            double zVtx_B = event_B->GetVtxZ();
            float cent_B = event_B->GetCentrality();
            float bbcq_B = event_B->GetBBCcharge();

            if (nclust_B==0) continue;
            if (cent_B<0 || cent_B>93) continue;
            if (bbcq_B>bbcq_cut_val) continue;
            if (TMath::Abs(zVtx_B)>bbcz_cut_val) continue;
            if (fabs(cent_A - cent_B) > 5.0) continue;
            if (fabs(zVtx_A - zVtx_B) > 2.5) continue;

            poold++;

            for (int iclust_B = 0; iclust_B < nclust_B; ++iclust_B)
            {
              MyCluster myclust_B = event_B->GetClusterEntry(iclust_B);
              float xclust_B = myclust_B.GetX();
              float yclust_B = myclust_B.GetY();
              float zclust_B = myclust_B.GetZ()-zVtx_B;
              float e_B = myclust_B.GetEcore() * cfgs[icfg].escale;
              float chi2_B = myclust_B.GetChi2();

              int sector_B = 7 - myclust_B.GetSect();
              if (!apply_emcmap(sector_B, myclust_B.GetIY(), myclust_B.GetIZ())) continue;
              if (e_B < cfgs[icfg].ecoreCut) continue;
              if (chi2_B > cfgs[icfg].chi2Cut) continue;

              double r_B = sqrt(xclust_B*xclust_B + yclust_B*yclust_B + zclust_B*zclust_B);
              if (r_B <= 0.0) continue;

              TLorentzVector emcPhoton_B;
              emcPhoton_B.Clear();
              emcPhoton_B.SetPx(e_B*xclust_B/r_B);
              emcPhoton_B.SetPy(e_B*yclust_B/r_B);
              emcPhoton_B.SetPz(e_B*zclust_B/r_B);
              emcPhoton_B.SetE(e_B);

              TLorentzVector pi0_AB;
              pi0_AB.Clear();
              pi0_AB = real_convPhoton_AA + emcPhoton_B;

              float Mass_AB = pi0_AB.M();
              float Pt_pi0_AB = pi0_AB.Pt();

              bg2d_eegmass_reco[icfg][icent]->Fill(Mass_AB, Pt_pi0_AB, w_pt);

              if (icfg == 1)
              {
                bg2d_eegmass_reco_dphi[icent]->Fill(Mass_AB, Pt_pi0_AB, w_pt);
              }
              if (icfg == 0)
              {
                bg2d_eegmass_reco_dzed[icent]->Fill(Mass_AB, Pt_pi0_AB, w_pt);
              }
            } // End Photon loop Event B
          } // End Event B loop
        } // End Track E/P 2 Loop
      } // End E/P A1 Loop
    } // End config loop
  } // End Event A loop

  write_histograms(output);

  input_A->Close();
  input_B->Close();
  output->Close();

  delete input_A;
  delete input_B;
  delete output;

  cout << "Done. Output: " << outFile << endl;
} // End make_TTree