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
#include <TGraph.h>

#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/MyEvent.h"
#include "/direct/phenix+u/roli/scratch/HELIOS/WriteEvent.h"
#include "/direct/phenix+u/tongzhouguo/install/include/DileptonAnalysis/Reconstruction.h"
#include "dc_dead_map.h"

using namespace std;
using namespace DileptonAnalysis;

//==============================================================
// Fast-mask systematic SIM macro.
//
// Main choices:
//   - same 14 systematic configs as data
//   - no BG histograms in sim
//   - histogram names match data pattern with _sim before icent
//   - no expensive conversion solution() call in sim
//       * nominal/conv_solution/all dzed configs: |dzed| < 4 cm
//       * conv_tight: |dzed| < 2 cm
//   - photon energy scale/smearing is kept for sim
//   - eeg histograms are M_ee_gamma vs pT_ee, not pT_ee_gamma
//   - fast bit-mask filling: no outer loop over configs around pairs/clusters
//==============================================================

static const double PI = 3.14159;
static const double Me = 0.000510998918;
static const double Me2 = Me*Me;
static const double c = 299792458;

const double DCENTERCUT = 10.0; // RICH ghost cut, same as data fast-mask macro
const double PC1_DPHI_CUT = 0.02;
const double PC1_DZ_CUT = 0.5;

// Keep five centrality bins for output compatibility, but fill only icent=0 for now.
const int centbin = 5;
const int sectbin = 1;

float EMCMAP[8][48][96];

// Histogram binning, matching the current real-data systematic production.
const int NMASS = 400;
const int NPT = 100;
const double MASS_MIN = 0.0;
const double MASS_MAX = 1.0;
const double PT_MIN = 0.0;
const double PT_MAX = 10.0;

const double MEE_MIN_FOR_EEG = 0.01;
const double MEE_MAX_FOR_EEG = 0.12;

// Sim event cuts kept from original sim macro.
double bbcz_cut_val = 10.0;
const double bbcq_cut_val = 1800.0;

// Photon scale/smear arrays from the original sim macro.
//const double scale[8] = {
//    0.988, 0.990, 0.985, 0.981,
//    0.980, 0.985, 1.025, 1.023};
//
//const double smear_c1[8] = {
//    0.055, 0.066, 0.055, 0.059,
//    0.057, 0.056, 0.072, 0.072};
//
//const double smear_c2[8] = {
//    0.011, 0.017, 0.011, 0.013,
//    0.012, 0.012, 0.063, 0.062};

const double scale[8] = {
    0.98, 0.99, 0.97, 0.97,
    0.98, 0.97, 1.03, 1.02};

const double smear_c1[8] = {
    0.044, 0.051, 0.057, 0.054,
    0.042, 0.048, 0.104, 0.112};

const double smear_c2[8] = {
    0.,0.,0.,0.,
    0.,0.,0.,0.};


// Conversion config types, same names as data.
enum ConvType
{
  CONV_SOLUTION = 0,
  CONV_DZED     = 1,
  CONV_TIGHT    = 2
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

TH1D* h1d_truepiz = 0;
TH1D* h1d_truepiz_unw = 0;
TH1F* h_nevents_sim = 0;
TH1F* h_cutflow_sim[NCFG];

TH2F* fg2d_eemass_reco_sim[NCFG][centbin];
TH2F* fg2d_eegmass_reco_sim[NCFG][centbin];
TH2F* fg2d_eemass_reco_sim_unw[NCFG][centbin];
TH2F* fg2d_eegmass_reco_sim_unw[NCFG][centbin];

TFile *ratio[10];
TGraph *gr[10];
TLorentzVector pion;

unsigned int cfg_bit(const int icfg)
{
  return (1u << icfg);
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
  ifstream readmap("/phenix/plhf/tongzhouguo/pi02gg/all_209_deadmap.txt");////"/phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/AnaTrain/test/offline1/AnalysisTrain/DileptonAnalysis/Run14AuAuEmcalDeadMap.txt"
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

// EID definitions matched to data systematics.
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

unsigned int get_track_pass_mask(MyTrack& trk)
{
  unsigned int mask = 0u;
  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    if (pass_eid(&trk, cfgs[icfg].eidType, cfgs[icfg].ptCut))
    {
      mask |= cfg_bit(icfg);
    }
  }
  return mask;
}

void compute_real_track_masks(std::vector<MyTrack>& tracks,
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

      if (getDcenter(tracks[i].GetCrkphi(), tracks[i].GetCrkz(),
                     tracks[j].GetCrkphi(), tracks[j].GetCrkz()) < DCENTERCUT)
      {
        ghost_mask |= common_mask;
      }
    }

    real_masks[i] = pass_masks[i] & (~ghost_mask);
  }
}

unsigned int get_conv_mask_sim(const double dzed)
{
  unsigned int mask = 0u;

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    bool pass = false;

    if (cfgs[icfg].convType == CONV_TIGHT)
    {
      // In sim we skip solution() but make tight slightly tighter in dzed.
      pass = (dzed < 2.0);
    }
    else
    {
      // nominal, conv_solution, and all other configs use dzed<4 in sim.
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
  int iy;
  int iz;
  int sector;
  double x;
  double y;
  double z;
  double e_smeared_base; // includes sector scale and stochastic smearing, not cfg escale
  double chi2;
};

ClusterInfo make_cluster_info(MyCluster& cl, TF1& f_rand)
{
  ClusterInfo ci;
  ci.ok_geom = false;
  ci.id = cl.GetID();
  ci.iy = cl.GetIY();
  ci.iz = cl.GetIZ();
  ci.sector = 7-cl.GetSect(); //  std to 0-7
  ci.x = cl.GetX();
  ci.y = cl.GetY();
  ci.z = cl.GetZ();
  ci.e_smeared_base = 0.0;
  ci.chi2 = cl.GetChi2();

  if (ci.sector < 0 || ci.sector >= 8) return ci;
  if (!apply_emcmap(ci.sector, ci.iy, ci.iz)) return ci;

  double ecore = cl.GetEcore();
  if (ecore <= 0.0) return ci;

  //const double resolution = sqrt(pow(smear_c1[ci.sector], 2) +
  //                               pow(smear_c2[ci.sector] / sqrt(ecore), 2));
  //const double e_smear = scale[ci.sector] * (1.0 + f_rand.GetRandom() * resolution);
  //ci.e_smeared_base = ecore * e_smear;
  //if (ci.e_smeared_base <= 0.0) return ci;
  ci.e_smeared_base = ecore;
  ci.ok_geom = true;
  return ci;
}

unsigned int get_cluster_mask_sim(const ClusterInfo& ci,
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

void book_histograms_sim()
{
  h_nevents_sim = new TH1F("h_nevents_sim", "sim event counter", 20, -0.5, 19.5);

  for (int icfg = 0; icfg < NCFG; ++icfg)
  {
    h_cutflow_sim[icfg] = new TH1F(Form("h_cutflow_%s_sim", cfgs[icfg].name),
                                   "sim cutflow", 12, -0.5, 11.5);
  }

  for (int icent = 0; icent < centbin; ++icent)
  {
    for (int icfg = 0; icfg < NCFG; ++icfg)
    {
      fg2d_eemass_reco_sim[icfg][icent] =
        new TH2F(Form("fg2d_eemass_reco_%s_sim_icent%d", cfgs[icfg].name, icent),
                 "sim eemass;M_{ee};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eegmass_reco_sim[icfg][icent] =
        new TH2F(Form("fg2d_eegmass_reco_%s_sim_icent%d", cfgs[icfg].name, icent),
                 "sim eegmass;M_{ee#gamma};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eemass_reco_sim_unw[icfg][icent] =
        new TH2F(Form("fg2d_eemass_reco_%s_sim_unw_icent%d", cfgs[icfg].name, icent),
                 "sim eemass unweighted;M_{ee};p_{T,ee}",
                 NMASS, MASS_MIN, MASS_MAX, NPT, PT_MIN, PT_MAX);

      fg2d_eegmass_reco_sim_unw[icfg][icent] =
        new TH2F(Form("fg2d_eegmass_reco_%s_sim_unw_icent%d", cfgs[icfg].name, icent),
                 "sim eegmass unweighted;M_{ee#gamma};p_{T,ee}",
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

void set_weighted_errors_sim()
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

void write_histograms_sim(TFile* output)
{
  output->cd();

  if (h_nevents_sim) h_nevents_sim->Write();
  if (h1d_truepiz) h1d_truepiz->Write();
  if (h1d_truepiz_unw) h1d_truepiz_unw->Write();

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

void get_Ntag_uncorr_sim_fastmask(int num = 100,
                         const char *inFile = "/phenix/plhf/tongzhouguo/test_tree2/1_tree",
                         const char *outFile = "sim_fastmask.root",
                         const int system = 0,
                         int iter = 0)
{
  TH1::AddDirectory(kFALSE);

  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
  gSystem->Load("/direct/phenix+u/roli/scratch/HELIOS/WriteEvent_C.so");
  gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");

  TF1 f_rand("f_rand_sim", "gaus", -5, 5);
  f_rand.SetParameter(0, 1);
  f_rand.SetParameter(1, 0);
  f_rand.SetParameter(2, 1);

  read_in_emcmap();

  TF1 *fHagedorn = new TF1("fHagedorn", "6.2831853*x*((1/(1+exp((x-[5])/[6])))*[0]/(1+x/[1])**[2] + (1-(1/(1+exp((x-[5])/[6]))))*[3]/x**[4])", 0.0, 40.0);
  fHagedorn->SetParameters(1.922170e+02, 2.106180e+00, 1.284030e+01, 1.630000e+01, 8.060480e+00, 4.033000e+00, 6.395340e-02);
  fHagedorn->SetNpx(10000);

  TFile *output = TFile::Open(outFile, "RECREATE");
  if (!output || output->IsZombie())
  {
    cout << "Cannot open output file: " << outFile << endl;
    if (output) { output->Close(); delete output; }
    delete fHagedorn;
    return;
  }
  ///direct/phenix+u/roli/scratch/HELIOS/simulation/singlePi0.root
  TFile *helios = TFile::Open("/direct/phenix+u/roli/scratch/HELIOS/simulation/singlePi0.root", "READ");///phenix/plhf/roli/HELIOS/simulation/singlePi0.root
  if (!helios || helios->IsZombie())
  {
    cout << "Cannot open HELIOS file" << endl;
    output->Close(); delete output;
    if (helios) { helios->Close(); delete helios; }
    delete fHagedorn;
    return;
  }

  TTree *helios_T = (TTree *)helios->Get("T");
  if (!helios_T)
  {
    cout << "Cannot find HELIOS tree T" << endl;
    helios->Close(); delete helios;
    output->Close(); delete output;
    delete fHagedorn;
    return;
  }

  TBranch *helios_br = helios_T->GetBranch("MyEvent");
  if (!helios_br)
  {
    cout << "Cannot find HELIOS branch MyEvent" << endl;
    helios->Close(); delete helios;
    output->Close(); delete output;
    delete fHagedorn;
    return;
  }

  WriteEvent *helios_event = 0;
  helios_br->SetAddress(&helios_event);
  cout << "HELIOS Tree read in" << endl;

  book_histograms_sim();
  cout << "Histograms defined" << endl;

  // Keep original true-pi0 normalization behavior.
  double w_pt = 1.0;
  for (int ievt = 0; ievt < 5000; ievt++)
  {
    helios_br->GetEntry(abs(0));
    if (!helios_event) continue;

    for (int i = 0; i < helios_event->GetNEntries(); i++)
    {
      WriteTrack helios_track = helios_event->GetWriteTrack(i);
      if (helios_track.GetID() == 111)
      {
        pion.SetPxPyPzE(helios_track.GetPx(), helios_track.GetPy(), helios_track.GetPz(), helios_track.GetEnergy());
      }
    }

    w_pt = fHagedorn->Eval(pion.Pt());
    for (int i = 0; i < iter; i++)
    {
      if (gr[i]) w_pt = w_pt * gr[i]->Eval(pion.Pt());
    }
    h1d_truepiz->Fill(pion.Pt(), w_pt);
    h1d_truepiz_unw->Fill(pion.Pt());
  }

  for (int iFile = 0; iFile < num; iFile++)
  {
    const TString inFile0 = Form("%s/sim_trees_%05d.root", inFile, iFile);

    TFile *input_A = TFile::Open(inFile0, "READ");
    if (!input_A || input_A->IsZombie())
    {
      cout << "no input_A file: " << inFile0 << endl;
      if (input_A) { input_A->Close(); delete input_A; }
      continue;
    }

    TTree *T_A = (TTree *)input_A->Get("T");
    if (!T_A)
    {
      cout << "Cannot find tree T in " << inFile0 << endl;
      input_A->Close(); delete input_A;
      continue;
    }

    TBranch *br_A = T_A->GetBranch("MyEvent");
    if (!br_A)
    {
      cout << "Cannot find branch MyEvent in " << inFile0 << endl;
      input_A->Close(); delete input_A;
      continue;
    }

    MyEvent *event = 0;
    br_A->SetAddress(&event);
    int nevt = T_A->GetEntries();

    cout << "Tree from " << inFile0 << " read in, nevt=" << nevt << endl;

    int nfolder = iFile;
    int nstart = nfolder * 10000;

    for (int ievent_A = 0; ievent_A < nevt; ievent_A++)
    {
      if (ievent_A % 1000 == 0) cout << "File " << iFile << " Event: " << ievent_A << " / " << nevt << endl;

      if (event) event->ClearEvent();
      br_A->GetEntry(ievent_A);
      if (!event) continue;
      h_nevents_sim->Fill(0);

      int ntrack_A = event->GetNtrack();
      int nclust_A = event->GetNcluster();
      double zVtx_A = event->GetVtxZ();
      float bbcq_A = event->GetBBCcharge();

      int id = event->GetEvtNo() - 1;
      helios_br->GetEntry(abs(nstart + id));
      if (!helios_event) continue;

      for (int i = 0; i < helios_event->GetNEntries(); i++)
      {
        WriteTrack helios_track = helios_event->GetWriteTrack(i);
        if (helios_track.GetID() == 111)
        {
          pion.SetPxPyPzE(helios_track.GetPx(), helios_track.GetPy(), helios_track.GetPz(), helios_track.GetEnergy());
        }
      }

      if (pion.Pt() > 10) continue;

      w_pt = fHagedorn->Eval(pion.Pt());
      for (int i = 0; i < iter; i++)
      {
        if (gr[i]) w_pt = w_pt * gr[i]->Eval(pion.Pt());
      }

      if (ntrack_A == 0) continue;
      if (TMath::Abs(zVtx_A) > bbcz_cut_val) continue;
      if (bbcq_A > bbcq_cut_val) continue;
      h_nevents_sim->Fill(1);

      // Keep flat sim centrality for now.
      int icent = 0;

      // Load tracks once.
      std::vector<MyTrack> tracks_A;
      tracks_A.reserve(ntrack_A);
      for (int itrk = 0; itrk < ntrack_A; ++itrk)
      {
        tracks_A.push_back(event->GetEntry(itrk));
      }

      std::vector<unsigned int> real_mask_A;
      compute_real_track_masks(tracks_A, real_mask_A);

      bool has_nominal_track = false;
      for (unsigned int itrk = 0; itrk < real_mask_A.size(); ++itrk)
      {
        if (real_mask_A[itrk] & cfg_bit(0))
        {
          has_nominal_track = true;
          break;
        }
      }
      if (!has_nominal_track) continue;
      h_nevents_sim->Fill(2);

      // Precompute clusters once per event. The random smearing is applied once per cluster.
      std::vector<ClusterInfo> clusters_A;
      clusters_A.reserve(nclust_A);
      for (int iclust_A = 0; iclust_A < nclust_A; ++iclust_A)
      {
        MyCluster myclust_A = event->GetClusterEntry(iclust_A);
        clusters_A.push_back(make_cluster_info(myclust_A, f_rand));
      }

      // Pair loop, once for all configs.
      for (int ireal_A1 = 0; ireal_A1 < ntrack_A; ++ireal_A1)
      {
        if (real_mask_A[ireal_A1] == 0u) continue;

        MyTrack& mytrk_A1 = tracks_A[ireal_A1];

        const int arm_A1 = mytrk_A1.GetArm();
        const int charge_A1 = mytrk_A1.GetCharge();
        const int emcid_A1 = mytrk_A1.GetEmcId();
        const double the_A1 = mytrk_A1.GetThe0();
        const double phi_A1 = mytrk_A1.GetPhi0();
        const double phidc_A1 = mytrk_A1.GetPhiDC();
        const double zed_A1 = mytrk_A1.GetZDC();
        const double mom_A1 = track_p(&mytrk_A1);
        const double pt_A1 = mom_A1 * sin(the_A1);

        const double momx_A1 = pt_A1 * cos(phi_A1);
        const double momy_A1 = pt_A1 * sin(phi_A1);
        const double momz_A1 = pt_A1 / tan(the_A1);

        for (int ireal_A2 = ireal_A1 + 1; ireal_A2 < ntrack_A; ++ireal_A2)
        {
          const unsigned int common_track_mask = real_mask_A[ireal_A1] & real_mask_A[ireal_A2];
          if (common_track_mask == 0u) continue;

          MyTrack& mytrk_A2 = tracks_A[ireal_A2];

          const int arm_A2 = mytrk_A2.GetArm();
          const int charge_A2 = mytrk_A2.GetCharge();
          const int emcid_A2 = mytrk_A2.GetEmcId();
          const double the_A2 = mytrk_A2.GetThe0();
          const double phi_A2 = mytrk_A2.GetPhi0();
          const double phidc_A2 = mytrk_A2.GetPhiDC();
          const double zed_A2 = mytrk_A2.GetZDC();
          const double mom_A2 = track_p(&mytrk_A2);
          const double pt_A2 = mom_A2 * sin(the_A2);

          if (emcid_A1 == emcid_A2) continue;
          if (fabs(phidc_A1 - phidc_A2) < PC1_DPHI_CUT && fabs(zed_A1 - zed_A2) < PC1_DZ_CUT) continue;
          if (charge_A1 == charge_A2) continue;
          if (arm_A1 != arm_A2) continue;

          const double dzed = fabs(zed_A1 - zed_A2);
          const unsigned int conv_mask = get_conv_mask_sim(dzed);
          const unsigned int pair_mask = common_track_mask & conv_mask;
          if (pair_mask == 0u) continue;

          const double momx_A2 = pt_A2 * cos(phi_A2);
          const double momy_A2 = pt_A2 * sin(phi_A2);
          const double momz_A2 = pt_A2 / tan(the_A2);

          TLorentzVector p_A1_DC;
          p_A1_DC.SetPxPyPzE(momx_A1, momy_A1, momz_A1,
                             sqrt(momx_A1*momx_A1 + momy_A1*momy_A1 + momz_A1*momz_A1 + Me2));

          TLorentzVector p_A2_DC;
          p_A2_DC.SetPxPyPzE(momx_A2, momy_A2, momz_A2,
                             sqrt(momx_A2*momx_A2 + momy_A2*momy_A2 + momz_A2*momz_A2 + Me2));

          TLorentzVector convPhoton_AA = p_A1_DC + p_A2_DC;

          TLorentzVector real_convPhoton_AA;
          real_convPhoton_AA.SetPx(momx_A1 + momx_A2);
          real_convPhoton_AA.SetPy(momy_A1 + momy_A2);
          real_convPhoton_AA.SetPz(momz_A1 + momz_A2);
          real_convPhoton_AA.SetE(sqrt(pow(convPhoton_AA.P(), 2)));

          const double Pt_AA = convPhoton_AA.Pt();
          const double Mass_ee = convPhoton_AA.M();

          for (int icfg = 0; icfg < NCFG; ++icfg)
          {
            if (!(pair_mask & cfg_bit(icfg))) continue;
            fg2d_eemass_reco_sim[icfg][icent]->Fill(Mass_ee, Pt_AA, w_pt);
            fg2d_eemass_reco_sim_unw[icfg][icent]->Fill(Mass_ee, Pt_AA);
            h_cutflow_sim[icfg]->Fill(0);
          }

          if (Mass_ee < MEE_MIN_FOR_EEG || Mass_ee > MEE_MAX_FOR_EEG) continue;

          for (int iclust_A = 0; iclust_A < nclust_A; ++iclust_A)
          {
            const ClusterInfo& ci = clusters_A[iclust_A];
            const unsigned int cl_mask = get_cluster_mask_sim(ci, emcid_A1, emcid_A2);
            const unsigned int pass_mask = pair_mask & cl_mask;
            if (pass_mask == 0u) continue;

            const double r_A = sqrt(ci.x*ci.x + ci.y*ci.y + (ci.z-zVtx_A)*(ci.z-zVtx_A));
            if (r_A <= 0.0) continue;

            for (int icfg = 0; icfg < NCFG; ++icfg)
            {
              if (!(pass_mask & cfg_bit(icfg))) continue;

              const double e_A = ci.e_smeared_base * cfgs[icfg].escale;

              TLorentzVector emcPhoton_A;
              emcPhoton_A.SetPxPyPzE(e_A * ci.x / r_A,
                                     e_A * ci.y / r_A,
                                     e_A * (ci.z-zVtx_A) / r_A,
                                     e_A);

              TLorentzVector pi0_AA = real_convPhoton_AA + emcPhoton_A;
              const double Mass_AA = pi0_AA.M();

              // IMPORTANT: y-axis is pT_ee, not pT_ee_gamma.
              fg2d_eegmass_reco_sim[icfg][icent]->Fill(Mass_AA, Pt_AA, w_pt);
              fg2d_eegmass_reco_sim_unw[icfg][icent]->Fill(Mass_AA, Pt_AA);
              h_cutflow_sim[icfg]->Fill(1);
            }
          }
        }
      }
    }

    input_A->Close();
    delete input_A;
  }

  set_weighted_errors_sim();
  write_histograms_sim(output);

  output->Close();
  helios->Close();

  delete output;
  delete helios;
  delete fHagedorn;

  cout << "Done. Output: " << outFile << endl;
}
