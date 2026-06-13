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

#include "dc_dead_map.h"
#include "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/svx_cent_ana/source/MyEvent.h"
#include "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/svx_cent_ana/source/Reconstruction.h"

using namespace std;
using namespace DileptonAnalysis;

//==============================================================
// Fast-mask systematic EMBEDDING macro.
//
// This keeps the original embedding input structure/libraries:
//   - function name: get_pi0_from_gamma_data
//   - tree name:     tree
//   - branch name:   MyEvent
//   - libraries:     libRun14AuAuLeptonEvent, libsvxcentana,
//                    libRun14AuAuLeptonConvReco
//   - embedding variables: GetPtPrime/GetThe0Prime/GetPhi0Prime,
//                          GetAlphaPrime/GetPreciseZ/GetMcId
//
// It mirrors the pure-MC fast-mask output:
//   - same 14 systematic configs as data/pure MC
//   - no BG histograms
//   - histogram names exactly like pure MC, with _sim before icent:
//       fg2d_eemass_reco_<cfg>_sim_icentX
//       fg2d_eegmass_reco_<cfg>_sim_icentX
//       fg2d_eemass_reco_<cfg>_sim_unw_icentX
//       fg2d_eegmass_reco_<cfg>_sim_unw_icentX
//   - no expensive conversion solution() call
//       * all non-tight configs: |dzed| < 4 cm
//       * conv_tight:            |dzed| < 2 cm
//   - photon energy scale/smearing kept
//   - only embedded photon clusters are used: cluster.GetMcId()==1
//   - eeg histograms are M_ee_gamma vs pT_ee, not pT_ee_gamma
//   - fast bit-mask filling: no outer config loop over tracks/pairs/clusters
//==============================================================

static const double Me = 0.000510998918;
static const double Me2 = Me*Me;

const double DCENTERCUT = 10.0;
const double PC1_DPHI_CUT = 0.02;
const double PC1_DZ_CUT = 0.5;

const int centbin = 5;
const int sectbin = 1;
const int Sectbin = 8;

const double MEE_MIN_FOR_EEG = 0.01;
const double MEE_MAX_FOR_EEG = 0.12;

const int NMASS = 400;
const int NPT = 100;
const double MASS_MIN = 0.0;
const double MASS_MAX = 1.0;
const double PT_MIN = 0.0;
const double PT_MAX = 10.0;

const float ncoll[] = {
  313.0,
  129.0,
  41.8,
  10.1,
  2.28
};

// Scale/smear arrays from original embedding macro.
const double scale[8] = {
    0.988, 0.990, 0.985, 0.981,
    0.980, 0.985, 1.025, 1.023
};

const double smear_c1[8] = {
    0.055, 0.066, 0.055, 0.059,
    0.057, 0.056, 0.072, 0.072
};

const double smear_c2[8] = {
    0.011, 0.017, 0.011, 0.013,
    0.012, 0.012, 0.063, 0.062
};

float EMCMAP[8][48][96];

// Old QA/compatibility histograms kept from original embedding macro.
TH2F* h2d_EMCdead[Sectbin];

TH1D* h1d_truepiz = 0;
TH1D* h1d_truepiz_unw = 0;
TH1F* h_nevents_sim = 0;
TH1F* h_cutflow_sim[14];

TH2F* fg2d_eemass_reco_sim[14][centbin];
TH2F* fg2d_eegmass_reco_sim[14][centbin];
TH2F* fg2d_eemass_reco_sim_unw[14][centbin];
TH2F* fg2d_eegmass_reco_sim_unw[14][centbin];

// Conversion config types, same as data/pure MC.
enum ConvType
{
  CONV_SOLUTION = 0,
  CONV_DZED     = 1,
  CONV_TIGHT    = 2
};


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

unsigned int cfg_bit(const int icfg)
{
  return (1u << icfg);
}

Double_t Hagedorn_Yield_Func(Double_t *x, Double_t *p)
{
  Double_t pt   = x[0];
  Double_t mass = p[0];
  Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.134*0.134);
  Double_t A    = p[1];
  Double_t a    = p[2];
  Double_t b    = p[3];
  Double_t p0   = p[4];
  Double_t n    = p[5];
  Double_t value = 2 * TMath::Pi() * pt * A * pow(exp(-a*mt-b*mt*mt)+mt/p0, -n);
  return value;
}

TF1 *funcHagedorn(const Char_t *name,
                  Double_t lowerlim,
                  Double_t upperlim,
                  const double dNdy = 95.7,
                  const double mass = 0.135,
                  const double A = 504.5,
                  const double a = 0.5169,
                  const double b = 0.1626,
                  const double p0 = 0.7366,
                  const double n = 8.274)
{
  TF1 *f = new TF1(name, Hagedorn_Yield_Func, lowerlim, upperlim, 6);
  f->SetParameters(mass, A, a, b, p0, n);
  f->SetParNames("mass", "A", "a", "b", "p0", "n");
  return f;
}

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
  ifstream readmap("/phenix/plhf/tongzhouguo/pi02gg/all_209_deadmap.txt");
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

  if (y == 0 || z == 0) return false;
  if (sect < 6 && (y == 35 || z == 71)) return false;
  if (sect > 5 && (y == 47 || z == 95)) return false;
  if (EMCMAP[sect][y][z]) return false;
  return true;
}

bool pass_base_track(MyTrack* mytrk)
{
  if (!mytrk) return false;

  int _arm = mytrk->GetArm();
  float _phi = mytrk->GetPhiDC();
  float _alpha = mytrk->GetAlphaPrime();
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
  if (fabs(_emcdz) > 20.0) return false;

  double _emcdphi = mytrk->GetEmcdphi();
  if (fabs(_emcdphi) > 0.05) return false;

  return true;
}

double track_pt(MyTrack* mytrk)
{
  return mytrk->GetPtPrime();
}

double track_p(MyTrack* mytrk)
{
  double pt = mytrk->GetPtPrime();
  double theta = mytrk->GetThe0Prime();
  double s = sin(theta);
  if (fabs(s) < 1e-12) return -1.0;
  return pt / s;
}

bool pass_eid(MyTrack* mytrk, int eidType, double ptCut)
{
  if (!mytrk) return false;
  if (!pass_base_track(mytrk)) return false;

  double pt = track_pt(mytrk);
  if (pt < ptCut) return false;

  double p = track_p(mytrk);
  if (p <= 0.0) return false;

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

unsigned int get_track_pass_mask(MyTrack* trk)
{
  unsigned int mask = 0u;
  if (!trk) return mask;
  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    if (pass_eid(trk, cfgs[icfg].eidType, cfgs[icfg].ptCut))
    {
      mask |= cfg_bit(icfg);
    }
  }
  return mask;
}

void compute_real_track_masks(std::vector<MyTrack*>& tracks,
                              std::vector<unsigned int>& real_masks)
{
  const int ntrack = tracks.size();
  std::vector<unsigned int> pass_masks(ntrack, 0u);
  real_masks.assign(ntrack, 0u);

  for (int i = 0; i < ntrack; ++i)
  {
    pass_masks[i] = get_track_pass_mask(tracks[i]);
  }

  for (int i = 0; i < ntrack; ++i)
  {
    if (pass_masks[i] == 0u) continue;

    unsigned int ghost_mask = 0u;
    for (int j = 0; j < ntrack; ++j)
    {
      if (i == j) continue;

      const unsigned int common_mask = pass_masks[i] & pass_masks[j];
      if (common_mask == 0u) continue;

      if (!tracks[i] || !tracks[j]) continue;
      if (getDcenter(tracks[i]->GetCrkphi(), tracks[i]->GetCrkz(),
                     tracks[j]->GetCrkphi(), tracks[j]->GetCrkz()) < DCENTERCUT)
      {
        ghost_mask |= common_mask;
      }
    }

    real_masks[i] = pass_masks[i] & (~ghost_mask);
  }
}

unsigned int get_conv_mask_embed(const double dzed)
{
  unsigned int mask = 0u;

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    bool pass = false;

    if (cfgs[icfg].convType == CONV_TIGHT)
    {
      pass = (dzed < 2.0);
    }
    else
    {
      pass = (dzed < 4.0);
    }

    if (pass) mask |= cfg_bit(icfg);
  }

  return mask;
}

struct ClusterInfo
{
  bool ok_geom;
  int id;
  int mcid;
  int iy;
  int iz;
  int sector;
  double x;
  double y;
  double z;
  double e_smeared_base;
  double chi2;
};

ClusterInfo make_cluster_info(MyCluster& cl, TF1& f_rand)
{
  ClusterInfo ci;
  ci.ok_geom = false;
  ci.id = cl.GetEmcId();
  ci.mcid = cl.GetMcId();
  ci.iy = cl.GetIY();
  ci.iz = cl.GetIZ();
  ci.sector = cl.GetSect();
  ci.x = cl.GetX();
  ci.y = cl.GetY();
  ci.z = cl.GetZ();
  ci.e_smeared_base = 0.0;
  ci.chi2 = cl.GetChi2();

  // embedding: only embedded photon clusters
  if (ci.mcid != 1) return ci;

  if (ci.sector < 0 || ci.sector >= 8) return ci;
  if (!apply_emcmap(ci.sector, ci.iy, ci.iz)) return ci;

  double ecore = cl.GetEcore();
  if (ecore <= 0.0) return ci;

  const double resolution = sqrt(pow(smear_c1[ci.sector], 2) +
                                 pow(smear_c2[ci.sector] / sqrt(ecore), 2));
  const double e_smear = scale[ci.sector] * (1.0 + f_rand.GetRandom() * resolution);
  ci.e_smeared_base = ecore * e_smear;
  //if (ci.e_smeared_base <= 0.0) return ci;
  //ci.e_smeared_base = ecore;
  if (ci.e_smeared_base <= 0.3) return ci;

  ci.ok_geom = true;
  return ci;
}

unsigned int get_cluster_mask_embed(const ClusterInfo& ci,
                                    const int emcid1,
                                    const int emcid2)
{
  if (!ci.ok_geom) return 0u;
  if (ci.id == emcid1 || ci.id == emcid2) return 0u;

  unsigned int mask = 0u;

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    const double e_cfg = ci.e_smeared_base * cfgs[icfg].escale;

    if (e_cfg < cfgs[icfg].ecoreCut) continue;
    if (ci.chi2 > cfgs[icfg].chi2Cut) continue;

    mask |= cfg_bit(icfg);
  }

  return mask;
}

void book_histograms_embed()
{
  h_nevents_sim = new TH1F("h_nevents_sim", "embedding event counter", 20, -0.5, 19.5);

  for (int isect = 0; isect < Sectbin; ++isect)
  {
    h2d_EMCdead[isect] = new TH2F(Form("h2d_EMCdead_sect%d", isect),
                                  "IY versus IZ", 48, -0.5, 47.5, 96, -0.5, 95.5);
  }

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    h_cutflow_sim[icfg] = new TH1F(Form("h_cutflow_%s_sim", cfgs[icfg].name),
                                   "embedding cutflow", 12, -0.5, 11.5);
  }

  for (int icent = 0; icent < centbin; ++icent)
  {
    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      fg2d_eemass_reco_sim[icfg][icent] =
        new TH2F(Form("fg2d_eemass_reco_%s_sim_icent%d", cfgs[icfg].name, icent),
                 "embedding eemass;M_{ee};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eegmass_reco_sim[icfg][icent] =
        new TH2F(Form("fg2d_eegmass_reco_%s_sim_icent%d", cfgs[icfg].name, icent),
                 "embedding eegmass;M_{ee#gamma};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eemass_reco_sim_unw[icfg][icent] =
        new TH2F(Form("fg2d_eemass_reco_%s_sim_unw_icent%d", cfgs[icfg].name, icent),
                 "embedding eemass unweighted;M_{ee};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eegmass_reco_sim_unw[icfg][icent] =
        new TH2F(Form("fg2d_eegmass_reco_%s_sim_unw_icent%d", cfgs[icfg].name, icent),
                 "embedding eegmass unweighted;M_{ee#gamma};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eemass_reco_sim[icfg][icent]->Sumw2();
      fg2d_eegmass_reco_sim[icfg][icent]->Sumw2();
      fg2d_eemass_reco_sim_unw[icfg][icent]->Sumw2();
      fg2d_eegmass_reco_sim_unw[icfg][icent]->Sumw2();
    }
  }

  h1d_truepiz = new TH1D("h1d_truepiz", "number of pizero generated", 400, 0.0, 40.0);
  h1d_truepiz_unw = new TH1D("h1d_truepiz_unw", "number of pizero generated", 400, 0.0, 40.0);
  h1d_truepiz->Sumw2();
  h1d_truepiz_unw->Sumw2();
}

void set_weighted_errors_embed()
{
  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    for (int icent = 0; icent < centbin; ++icent)
    {
      TH2F* h_ee = fg2d_eemass_reco_sim[icfg][icent];
      TH2F* h_ee_unw = fg2d_eemass_reco_sim_unw[icfg][icent];
      TH2F* h_eeg = fg2d_eegmass_reco_sim[icfg][icent];
      TH2F* h_eeg_unw = fg2d_eegmass_reco_sim_unw[icfg][icent];

      for (int ix = 1; ix <= h_ee->GetNbinsX(); ++ix)
      {
        for (int iy = 1; iy <= h_ee->GetNbinsY(); ++iy)
        {
          double n_unw = h_ee_unw->GetBinContent(ix, iy);
          if (n_unw > 0.0)
          {
            double content = h_ee->GetBinContent(ix, iy);
            h_ee->SetBinError(ix, iy, content / sqrt(n_unw));
          }

          n_unw = h_eeg_unw->GetBinContent(ix, iy);
          if (n_unw > 0.0)
          {
            double content = h_eeg->GetBinContent(ix, iy);
            h_eeg->SetBinError(ix, iy, content / sqrt(n_unw));
          }
        }
      }
    }
  }
}

void write_histograms_embed(TFile* output)
{
  output->cd();

  if (h_nevents_sim) h_nevents_sim->Write();
  if (h1d_truepiz) h1d_truepiz->Write();
  if (h1d_truepiz_unw) h1d_truepiz_unw->Write();

  for (int isect = 0; isect < Sectbin; ++isect)
  {
    if (h2d_EMCdead[isect]) h2d_EMCdead[isect]->Write();
  }

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    if (h_cutflow_sim[icfg]) h_cutflow_sim[icfg]->Write();
  }

  for (int icent = 0; icent < centbin; ++icent)
  {
    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      fg2d_eemass_reco_sim[icfg][icent]->Write();
      fg2d_eegmass_reco_sim[icfg][icent]->Write();
      fg2d_eemass_reco_sim_unw[icfg][icent]->Write();
      fg2d_eegmass_reco_sim_unw[icfg][icent]->Write();
    }
  }
}

void get_pi0_from_gamma_data(
                             const int system = 0,
                             const char* inFile0 = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/analysis/emb/372586/",
                             const char* outFile0 = "embedding_fastmask",
                             int ert = 0,
                             const int debug = 1)
{
  TH1::AddDirectory(kFALSE);

  // Keep exactly the active libraries from the uploaded embedding file.
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonEvent.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libsvxcentana.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonConvReco.so");

  // Constructed to keep the same environment/structure as the original file,
  // but solution() is intentionally not called in the fast systematic pass.
  Reconstruction reco("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/analysis/emb/lookup_3D_one_phi.root");

  TF1 f_rand("f_rand_embed", "gaus", -5, 5);
  f_rand.SetParameter(0, 1);
  f_rand.SetParameter(1, 0);
  f_rand.SetParameter(2, 1);

  TString inFile  = Form("%s/emb_%d.root", inFile0, system);
  TString outFile = Form("%s_%d.root", outFile0, system);

  std::cout << "Input file: " << inFile.Data()
            << "  Output file: " << outFile.Data() << std::endl;

  read_in_emcmap();

  TFile* input = TFile::Open(inFile.Data(), "READ");
  if (!input || input->IsZombie())
  {
    cout << "no input file: " << inFile.Data() << endl;
    if (input) { input->Close(); delete input; }
    return;
  }

  TTree* T = (TTree*)input->Get("tree");
  if (!T)
  {
    cout << "Cannot find tree named tree" << endl;
    input->Close(); delete input;
    return;
  }

  TBranch* br = T->GetBranch("MyEvent");
  if (!br)
  {
    cout << "Cannot find branch MyEvent" << endl;
    input->Close(); delete input;
    return;
  }

  MyEvent* event = 0;
  br->SetAddress(&event);
  cout << "Input Tree read in" << endl;

  TF1* fHagedorn[centbin];
  for (int icent = 0; icent < centbin; icent++)
  {
    fHagedorn[icent] = new TF1(Form("fHagedorn_%d", icent),
                               "[7]*6.2831853*x*((1/(1+exp((x-[5])/[6])))*[0]/(1+x/[1])**[2] + (1-(1/(1+exp((x-[5])/[6]))))*[3]/x**[4])",
                               0.0, 10.0);
    fHagedorn[icent]->SetNpx(10000);
    fHagedorn[icent]->SetParameters(1.922170e+02, 2.106180e+00, 1.284030e+01,
                                    1.630000e+01, 8.060480e+00, 4.033000e+00,
                                    6.395340e-02, 1.0*ncoll[icent]);
  }

  book_histograms_embed();
  cout << "Histograms defined" << endl;

  int nevt = T->GetEntries();
  cout << "Processing nevt = " << nevt << endl;

  // Debug counters. They are printed at the end and help locate where entries disappear.
  long long dbg_evt_read = 0;
  long long dbg_evt_ntrack = 0;
  long long dbg_evt_cent = 0;
  long long dbg_evt_has_any_real_track = 0;
  long long dbg_track_input = 0;
  long long dbg_track_base_pass = 0;
  long long dbg_track_any_mask = 0;
  long long dbg_track_real_any = 0;
  long long dbg_pair_all = 0;
  long long dbg_pair_common_mask = 0;
  long long dbg_pair_emcid_pass = 0;
  long long dbg_pair_pc1_pass = 0;
  long long dbg_pair_unlike = 0;
  long long dbg_pair_samearm = 0;
  long long dbg_pair_dzed_pass = 0;
  long long dbg_ee_fill = 0;
  long long dbg_pair_mee_window = 0;
  long long dbg_cluster_input = 0;
  long long dbg_cluster_ok_geom = 0;
  long long dbg_cluster_pass_mask = 0;
  long long dbg_eeg_fill = 0;
  long long dbg_cfg_ee[NCFG] = {0};
  long long dbg_cfg_eeg[NCFG] = {0};
  long long dbg_cent_evt[centbin] = {0};
  long long dbg_cent_ee[centbin] = {0};
  long long dbg_cent_eeg[centbin] = {0};

  for (int ievent_A = 0; ievent_A < nevt; ievent_A++)
  {
    if (ievent_A % 100000 == 0) cout << "Event: " << ievent_A << " / " << nevt << endl;
    //if (ievent_A >= 1000) break;
    if (event) event->ClearEvent();
    br->GetEntry(ievent_A);
    if (!event) continue;
    h_nevents_sim->Fill(0);
    ++dbg_evt_read;

    const int ntrack_A = event->GetNtrack();
    const int nclust_A = event->GetNcluster();
    const int cent = event->GetCentrality()*0;///podgon
    const double zVtx_A = event->GetVtxZ();
    const int n_gen_track = event->GetNgentrack();

    if (ntrack_A < 2 ) continue;
    ++dbg_evt_ntrack;
    if (cent < 0 || cent > 93) continue;
    ++dbg_evt_cent;

    int cent_A = cent / (100 / centbin);
    if (cent_A < 0) cent_A = 0;
    if (cent_A >= centbin) cent_A = centbin - 1;
    ++dbg_cent_evt[cent_A];

    // True embedded pi0 pT from generated tracks, as in original embedding macro.
    float pT_x = 0.0;
    float pT_y = 0.0;
    for (int i_gen_track = 0; i_gen_track < n_gen_track; ++i_gen_track)
    {
      MyGenTrack my_gen_track = *event->GetGenTrack(i_gen_track);
      pT_x += my_gen_track.GetPx();
      pT_y += my_gen_track.GetPy();
    }
    const float pT_pi0 = sqrt(pT_x*pT_x + pT_y*pT_y);

    float w_pt_tr = fHagedorn[cent_A]->Eval(pT_pi0);
    h1d_truepiz->Fill(pT_pi0, w_pt_tr);
    h1d_truepiz_unw->Fill(pT_pi0);
    h_nevents_sim->Fill(1);

    // Load track pointers once. Do NOT copy MyTrack objects here:
    // in this embedding class some members can be lost/invalid after by-value copy.
    std::vector<MyTrack*> tracks_A;
    tracks_A.reserve(ntrack_A);
    for (int itrk = 0; itrk < ntrack_A; ++itrk)
    {
      MyTrack* trkp = event->GetEntry(itrk);
      if (trkp) tracks_A.push_back(trkp);
    }

    const int ntracks_loaded = tracks_A.size();
    if (ntracks_loaded == 0) continue;
    dbg_track_input += ntracks_loaded;

    std::vector<unsigned int> real_mask_A;
    compute_real_track_masks(tracks_A, real_mask_A);

    bool has_any_real_track = false;
    int ev_base = 0;
    int ev_any_mask = 0;
    int ev_real_any = 0;
    int ev_real_nom = 0;
    int ev_real_eid0 = 0;
    int ev_real_eid1 = 0;
    int ev_real_eid3 = 0;

    for (unsigned int itrk = 0; itrk < real_mask_A.size(); ++itrk)
    {
      if (debug)
      {
        if (pass_base_track(tracks_A[itrk])) ++ev_base;
        const unsigned int raw_mask = get_track_pass_mask(tracks_A[itrk]);
        if (raw_mask) ++ev_any_mask;
        else if (pass_base_track(tracks_A[itrk]) && ievent_A < 20)
        {
          MyTrack* trk_dbg = tracks_A[itrk];
          const double pt_dbg = trk_dbg->GetPtPrime();
          const double theta_dbg = trk_dbg->GetThe0Prime();
          const double p_dbg = (fabs(sin(theta_dbg)) > 1e-12) ? pt_dbg / sin(theta_dbg) : -999.0;
          const double ecore_dbg = trk_dbg->GetEcore();
          const double ep_dbg = (p_dbg > 0.0) ? ecore_dbg / p_dbg : -999.0;

          cout << "DBG rawmask0 track"
               << " event=" << ievent_A
               << " itrk=" << itrk
               << " pt=" << pt_dbg
               << " theta=" << theta_dbg
               << " p=" << p_dbg
               << " ecore=" << ecore_dbg
               << " ep=" << ep_dbg
               << " n0=" << trk_dbg->GetN0()
               << " disp=" << trk_dbg->GetDisp()
               << " prob=" << trk_dbg->GetProb()
               << " charge=" << trk_dbg->GetCharge()
               << " emcdz=" << trk_dbg->GetEmcdz()
               << " emcdphi=" << trk_dbg->GetEmcdphi()
               << endl;
        }
      }

      if (real_mask_A[itrk] != 0u)
      {
        has_any_real_track = true;
        ++ev_real_any;
        if (real_mask_A[itrk] & cfg_bit(0)) ++ev_real_nom;
        if (real_mask_A[itrk] & cfg_bit(7)) ++ev_real_eid0;
        if (real_mask_A[itrk] & cfg_bit(8)) ++ev_real_eid1;
        if (real_mask_A[itrk] & cfg_bit(9)) ++ev_real_eid3;
      }
    }

    if (debug)
    {
      dbg_track_base_pass += ev_base;
      dbg_track_any_mask += ev_any_mask;
      dbg_track_real_any += ev_real_any;

      if (ievent_A < 20)
      {
        cout << "DBG event " << ievent_A
             << " cent=" << cent
             << " cent_A=" << cent_A
             << " ntrack=" << ntrack_A
             << " nclust=" << nclust_A
             << " ngen=" << n_gen_track
             << " pTpi0=" << pT_pi0
             << " w=" << w_pt_tr
             << " base=" << ev_base
             << " rawmask=" << ev_any_mask
             << " real_any=" << ev_real_any
             << " real_nom=" << ev_real_nom
             << " real_eid0=" << ev_real_eid0
             << " real_eid1=" << ev_real_eid1
             << " real_eid3=" << ev_real_eid3
             << endl;
      }
    }

    if (!has_any_real_track) continue;
    ++dbg_evt_has_any_real_track;
    h_nevents_sim->Fill(2);

    // Precompute clusters once per event. Smearing is applied once per cluster.
    std::vector<ClusterInfo> clusters_A;
    clusters_A.reserve(nclust_A);
    for (int iclust_A = 0; iclust_A < nclust_A; ++iclust_A)
    {
      MyCluster* clp = event->GetClusterEntry(iclust_A);
      if (!clp) continue;
      ClusterInfo ci = make_cluster_info(*clp, f_rand);
      if (debug)
      {
        ++dbg_cluster_input;
        if (ci.ok_geom) ++dbg_cluster_ok_geom;
      }
      clusters_A.push_back(ci);
    }

    // Pair loop, once for all configs.
    for (int ireal_A1 = 0; ireal_A1 < ntracks_loaded; ++ireal_A1)
    {
      if (real_mask_A[ireal_A1] == 0u) continue;

      MyTrack* mytrk_A1p = tracks_A[ireal_A1];
      if (!mytrk_A1p) continue;
      MyTrack& mytrk_A1 = *mytrk_A1p;

      const int arm_A1 = mytrk_A1.GetArm();
      const int charge_A1 = mytrk_A1.GetCharge();
      const int emcid_A1 = mytrk_A1.GetEmcId();
      const double the_A1 = mytrk_A1.GetThe0Prime();
      const double phi_A1 = mytrk_A1.GetPhi0Prime();
      const double phidc_A1 = mytrk_A1.GetPhiDC();
      const double zed_A1 = mytrk_A1.GetZDC();
      const double pt_A1 = mytrk_A1.GetPtPrime();

      const double momx_A1 = pt_A1 * cos(phi_A1);
      const double momy_A1 = pt_A1 * sin(phi_A1);
      const double momz_A1 = pt_A1 / tan(the_A1);

      for (int ireal_A2 = ireal_A1 + 1; ireal_A2 < ntracks_loaded; ++ireal_A2)
      {
        ++dbg_pair_all;
        const unsigned int common_track_mask = real_mask_A[ireal_A1] & real_mask_A[ireal_A2];
        if (common_track_mask == 0u) continue;
        ++dbg_pair_common_mask;

        MyTrack* mytrk_A2p = tracks_A[ireal_A2];
        if (!mytrk_A2p) continue;
        MyTrack& mytrk_A2 = *mytrk_A2p;

        const int arm_A2 = mytrk_A2.GetArm();
        const int charge_A2 = mytrk_A2.GetCharge();
        const int emcid_A2 = mytrk_A2.GetEmcId();
        const double the_A2 = mytrk_A2.GetThe0Prime();
        const double phi_A2 = mytrk_A2.GetPhi0Prime();
        const double phidc_A2 = mytrk_A2.GetPhiDC();
        const double zed_A2 = mytrk_A2.GetZDC();
        const double pt_A2 = mytrk_A2.GetPtPrime();

        if (emcid_A1 == emcid_A2) continue;
        ++dbg_pair_emcid_pass;
        if (fabs(phidc_A1 - phidc_A2) < PC1_DPHI_CUT && fabs(zed_A1 - zed_A2) < PC1_DZ_CUT) continue;
        ++dbg_pair_pc1_pass;
        if (charge_A1 == charge_A2) continue;
        ++dbg_pair_unlike;
        if (arm_A1 != arm_A2) continue;
        ++dbg_pair_samearm;

        const double dzed = fabs(zed_A1 - zed_A2);
        const unsigned int conv_mask = get_conv_mask_embed(dzed);
        const unsigned int pair_mask = common_track_mask & conv_mask;
        if (pair_mask == 0u) continue;
        ++dbg_pair_dzed_pass;

        //checking solution() 
        MyPair mypair_AA;
        mypair_AA.Clear();
        //if (!solution(&mytrk_A1, &mytrk_A2, &mypair_AA, &reco, zVtx_A)) continue;

        const double momx_A2 = pt_A2 * cos(phi_A2);
        const double momy_A2 = pt_A2 * sin(phi_A2);
        const double momz_A2 = pt_A2 / tan(the_A2);

        TLorentzVector p_A1_DC;
        p_A1_DC.Clear();
        p_A1_DC.SetPxPyPzE(momx_A1, momy_A1, momz_A1,
                           sqrt(momx_A1*momx_A1 + momy_A1*momy_A1 + momz_A1*momz_A1 + Me2));

        TLorentzVector p_A2_DC;
        p_A2_DC.Clear();
        p_A2_DC.SetPxPyPzE(momx_A2, momy_A2, momz_A2,
                           sqrt(momx_A2*momx_A2 + momy_A2*momy_A2 + momz_A2*momz_A2 + Me2));

        TLorentzVector convPhoton_AA;
        convPhoton_AA.Clear();
        convPhoton_AA = p_A1_DC + p_A2_DC;

        TLorentzVector real_convPhoton_AA;
        real_convPhoton_AA.Clear();
        real_convPhoton_AA.SetPx(momx_A1 + momx_A2);
        real_convPhoton_AA.SetPy(momy_A1 + momy_A2);
        real_convPhoton_AA.SetPz(momz_A1 + momz_A2);
        real_convPhoton_AA.SetE(sqrt(pow(convPhoton_AA.P(), 2)));

        const double Pt_AA = real_convPhoton_AA.Pt();
        const double Mass_ee = convPhoton_AA.M();

        for (int icfg = 0; icfg < NCFG; ++icfg)
        {
          if (!(pair_mask & cfg_bit(icfg))) continue;
          fg2d_eemass_reco_sim[icfg][cent_A]->Fill(Mass_ee, Pt_AA, w_pt_tr);
          fg2d_eemass_reco_sim_unw[icfg][cent_A]->Fill(Mass_ee, Pt_AA);
          h_cutflow_sim[icfg]->Fill(0);
          ++dbg_ee_fill;
          ++dbg_cfg_ee[icfg];
          ++dbg_cent_ee[cent_A];
        }

        if (Mass_ee < MEE_MIN_FOR_EEG || Mass_ee > MEE_MAX_FOR_EEG) continue;
        ++dbg_pair_mee_window;

        for (unsigned int iclust_A = 0; iclust_A < clusters_A.size(); ++iclust_A)
        {
          const ClusterInfo& ci = clusters_A[iclust_A];
          const unsigned int cl_mask = get_cluster_mask_embed(ci, emcid_A1, emcid_A2);
          const unsigned int pass_mask = pair_mask & cl_mask;
          if (pass_mask == 0u) continue;
          ++dbg_cluster_pass_mask;

          const double z_A = ci.z - zVtx_A;
          const double r_A = sqrt(ci.x*ci.x + ci.y*ci.y + z_A*z_A);
          if (r_A <= 0.0) continue;

          if (ci.sector >= 0 && ci.sector < Sectbin)
          {
            h2d_EMCdead[ci.sector]->Fill(ci.iy, ci.iz);
          }

          for (int icfg = 0; icfg < NCFG; ++icfg)
          {
            if (!(pass_mask & cfg_bit(icfg))) continue;

            const double e_A = ci.e_smeared_base * cfgs[icfg].escale;

            TLorentzVector emcPhoton_A;
            emcPhoton_A.Clear();
            emcPhoton_A.SetPx(e_A*ci.x/r_A);
            emcPhoton_A.SetPy(e_A*ci.y/r_A);
            emcPhoton_A.SetPz(e_A*z_A/r_A);
            emcPhoton_A.SetE(e_A);

            TLorentzVector pi0_AA;
            pi0_AA.Clear();
            pi0_AA = real_convPhoton_AA + emcPhoton_A;

            float Mass_AA = pi0_AA.M();

            //float pt_pi0_reco = pi0_AA.Pt();
            //float diff = (pt_pi0_reco-pT_pi0)/pT_pi0;
            //if(diff>0.1) continue;

            // IMPORTANT: y-axis is pT_ee, matching pure MC/data handling.
            fg2d_eegmass_reco_sim[icfg][cent_A]->Fill(Mass_AA, Pt_AA, w_pt_tr);
            fg2d_eegmass_reco_sim_unw[icfg][cent_A]->Fill(Mass_AA, Pt_AA);
            h_cutflow_sim[icfg]->Fill(1);
            ++dbg_eeg_fill;
            ++dbg_cfg_eeg[icfg];
            ++dbg_cent_eeg[cent_A];
          }
        }
      }
    }
  }

  if (debug)
  {
    cout << "================ EMBEDDING FASTMASK DEBUG SUMMARY ================" << endl;
    cout << "events: read=" << dbg_evt_read
         << " ntrack_pass=" << dbg_evt_ntrack
         << " cent_pass=" << dbg_evt_cent
         << " has_any_real_track=" << dbg_evt_has_any_real_track
         << endl;
    cout << "tracks: input=" << dbg_track_input
         << " base_pass=" << dbg_track_base_pass
         << " any_raw_mask=" << dbg_track_any_mask
         << " real_any=" << dbg_track_real_any
         << endl;
    cout << "pairs: all=" << dbg_pair_all
         << " common_mask=" << dbg_pair_common_mask
         << " emcid_pass=" << dbg_pair_emcid_pass
         << " pc1_pass=" << dbg_pair_pc1_pass
         << " unlike=" << dbg_pair_unlike
         << " samearm=" << dbg_pair_samearm
         << " dzed_pass=" << dbg_pair_dzed_pass
         << " ee_fills=" << dbg_ee_fill
         << " mee_window=" << dbg_pair_mee_window
         << endl;
    cout << "clusters: input=" << dbg_cluster_input
         << " ok_embed_geom_map=" << dbg_cluster_ok_geom
         << " pass_pair_cluster_mask=" << dbg_cluster_pass_mask
         << " eeg_fills=" << dbg_eeg_fill
         << endl;
    cout << "centrality summary:" << endl;
    for (int icent = 0; icent < centbin; ++icent)
    {
      cout << "  icent " << icent
           << " events=" << dbg_cent_evt[icent]
           << " ee_fills=" << dbg_cent_ee[icent]
           << " eeg_fills=" << dbg_cent_eeg[icent]
           << endl;
    }
    cout << "config summary:" << endl;
    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      double int_ee = 0.0;
      double int_eeg = 0.0;
      for (int icent = 0; icent < centbin; ++icent)
      {
        if (fg2d_eemass_reco_sim[icfg][icent]) int_ee += fg2d_eemass_reco_sim[icfg][icent]->Integral();
        if (fg2d_eegmass_reco_sim[icfg][icent]) int_eeg += fg2d_eegmass_reco_sim[icfg][icent]->Integral();
      }
      cout << "  " << cfgs[icfg].name
           << " ee_fills=" << dbg_cfg_ee[icfg]
           << " eeg_fills=" << dbg_cfg_eeg[icfg]
           << " ee_integral=" << int_ee
           << " eeg_integral=" << int_eeg
           << endl;
    }
    cout << "==================================================================" << endl;
  }

  set_weighted_errors_embed();

  TFile* output = TFile::Open(outFile.Data(), "RECREATE");
  if (!output || output->IsZombie())
  {
    cout << "Cannot open output file: " << outFile.Data() << endl;
    if (output) { output->Close(); delete output; }
    input->Close(); delete input;
    for (int ic = 0; ic < centbin; ++ic) delete fHagedorn[ic];
    return;
  }

  write_histograms_embed(output);

  output->Close();
  input->Close();

  delete output;
  delete input;
  for (int ic = 0; ic < centbin; ++ic) delete fHagedorn[ic];

  cout << "Done. Output: " << outFile.Data() << endl;
}

// Optional wrapper with a more explicit name. The original function name above is kept.
void get_Ntag_uncorr_embedding_fastmask(const char* inFile = "../tong1.root",
                                        const char* outFile = "embedding_fastmask.root",
                                        const int system = 0,
                                        int ert = 0,
                                        const int debug = 1)
{
  get_pi0_from_gamma_data(system, inFile, outFile, ert, debug);
}
