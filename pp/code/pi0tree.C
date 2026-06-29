#define pi0tree_cxx
#include "pi0tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void mass_vs_pt(const char * outName = "mass_vs_pt_sim.root", const char * inName = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/simppizero_DSF_simonly/TTree_out_pi0_30GeV_00000.root"){
   std::cout << "I GOT HERE!!!" << std::endl;
   TChain *t_pi0tree = new TChain("pi0tree");
   t_pi0tree->Add(inName); // only the  1 files
   pi0tree t(t_pi0tree);

   TFile * fout = new TFile(outName,"RECREATE");
   // Sector dependent histograms

   THnSparseD * iyiz_a;
   Int_t bins_iyiz_a[3] =    {6, 96, 48};
   Double_t xmin_iyiz_a[3] = {-0.5, -0.5, -0.5};
   Double_t xmax_iyiz_a[3] = {5.5, 95.5, 47.5};
   iyiz_a = new THnSparseD("iyiz_a","hit dist for cluster A; sec; iz; iy;", 3, bins_iyiz_a, xmin_iyiz_a, xmax_iyiz_a);
   
   THnSparseD * x11_y11;
   Int_t bins_x11_y11[3] =    {6, 1000, 1000};
   Double_t xmin_x11_y11[3] = {-0.5, -1000, -1000};
   Double_t xmax_x11_y11[3] = {5.5, 1000, 1000};
   x11_y11 = new THnSparseD("x11_y11","hit dist X11 vs Y11; sec; x11; y11;", 3, bins_x11_y11, xmin_x11_y11, xmax_x11_y11);

   THnSparseD * mass_vs_pt;
   THnSparseD * mass_vs_pt_uw;
   Int_t bins_mass_vs_pt[3] =    {6, 160, 60};
   Double_t xmin_mass_vs_pt[3] = {-0.5, 0, 0};
   Double_t xmax_mass_vs_pt[3] = {5.5, 0.8, 30};

   mass_vs_pt = new THnSparseD("mass_vs_pt","mass vs pT; sec; m_gg [GeV/c^2]; pair p_T [GeV]", 3, bins_mass_vs_pt, xmin_mass_vs_pt, xmax_mass_vs_pt);      
   mass_vs_pt_uw = new THnSparseD("mass_vs_pt_uw","mass vs pT unweighted; sec; m_gg [GeV/c^2]; pair p_T [GeV]", 3, bins_mass_vs_pt, xmin_mass_vs_pt, xmax_mass_vs_pt);      

   
   TH1D * pT_a = new TH1D("pT_a","pT distribution of A cluster",240,0.0,30.0);
   TH1D * truepT = new TH1D("truepT","true pT distribution",240,0.0,30.0);
   
   TH1D * dist_gea = new TH1D("dist_gea","distance between geant hits and reco hits; sqrt((impx-x)^2 + (impy-y)^2 + (impz-z)^2); counts",5000,0.0,1200);

   // Loading DH maps
   // loading Dead Hot maps   
   TFile * fDeadHot = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Analysis/TTree_Analysis/output_ERT_iyvsiz2/DHmaps.root");
   TH2D * dhmaps[6];
   for(int s = 0; s < 6; s++){
      dhmaps[s] = (TH2D*)fDeadHot->Get(Form("iy_iz_masked3x3_sec%d",s));
      dhmaps[s]->SetDirectory(0);
   }
   fDeadHot->Close();
   
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_1.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_2.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_3.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_4.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_altbin_5.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_altbin_6.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_altbin_7.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_altbin_8.root","READ");
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_altbin_9.root","READ");

   // This is my last iteration
   //TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_altbin_10.root","READ");
   
   // this is norbert's as a wheight
   TFile * f_weights = new TFile("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/d131417/Simulation/Pi0_2/wheight_for_simu_Norb_1.root","READ");

   TGraph * g_weights = (TGraph*) f_weights->Get("g_pi0_weight");
   g_weights->SetName("g_weights");
   f_weights->Close();
   //TF1* fHagedorn = new TF1("fHagedorn","([0]*2*3.14159*x)/pow((exp(-[1]*x -[2]*x*x) + x/[3]),[4])", 0.1, 40);
   //fHagedorn->SetParameters(3.51316e+01, -2.46986e-01, 9.76431e-03, 2.97233e+00, 9.81597e+00);
   //fHagedorn->SetParameters(1.33695e+02, 1.80695e-01, 1.21563e-01, 6.68003e-01, 8.25566e+00);

   cout << "begin event loop" << endl;
   Long64_t nentries = t.fChain->GetEntries();

   Double_t filler_iyiz_a[3];
   Double_t filler_x11_y11[3];
   Double_t filler_mass_vs_pt[3];
   Double_t filler_mass_vs_pt_uw[3];

   for (int jentry=0; jentry<nentries;jentry++) {
      if(jentry % 1000 == 0){cout << jentry <<"/" << nentries <<endl;}
      t.fChain->GetEntry(jentry);
      truepT->Fill((t.vpt)); // filling true pT before anything, at the event loop

      if( !((*t.clusters_ecore1).size() > 0) ) {continue;}; 
      
      //cout << "foo" << endl;
      double true_pt = (t.vpt);

      for(int i = 0; i < 3; i++){
         filler_iyiz_a[i] = 0;
         filler_x11_y11[i] = 0;
         filler_mass_vs_pt[i] = 0;
         filler_mass_vs_pt_uw[i] = 0;
      }
      
      for(int a = 0; a < (*t.clusters_nGeaCluster).size() -1; a++){
         int sec_a = (*t.clusters_sec1)[a];
         if(sec_a > 5) continue;
         double ecore_a = (*t.clusters_ecore1)[a];
         if(ecore_a < 0.3 ) continue;
         double chi2_a = (*t.clusters_chi21)[a];
         if(chi2_a > 3) continue;
         // ADD APPLY DH MAPS HERE
         if(dhmaps[sec_a]->GetBinContent(dhmaps[sec_a]->FindBin((*t.clusters_iz)[a],(*t.clusters_iy)[a])) > 0.1) continue; // 0.1 just to account for possible double precision

         // Fill iyiz histos
         //iyiz_a[sec_a]->Fill((*t.clusters_iz)[a],(*t.clusters_iy)[a]);
         //x11_y11[sec_a]->Fill((*t.clusters_x11)[a],(*t.clusters_y11)[a]);
         
         filler_iyiz_a[0] = sec_a;
         filler_iyiz_a[1] = (*t.clusters_iz)[a];
         filler_iyiz_a[2] = (*t.clusters_iy)[a];
         iyiz_a->Fill(filler_iyiz_a);

         filler_x11_y11[0] = sec_a;
         filler_x11_y11[1] = (*t.clusters_x11)[a];
         filler_x11_y11[2] = (*t.clusters_y11)[a];
         x11_y11->Fill(filler_x11_y11);

         double dist_gea_calc = sqrt( 
            pow((*t.clusters_impx11)[a]-(*t.clusters_x1)[a],2) +  
            pow((*t.clusters_impy11)[a]-(*t.clusters_y1)[a],2) +  
            pow((*t.clusters_impz11)[a]-(*t.clusters_z1)[a],2)
         );

         dist_gea->Fill(dist_gea_calc);

         // tower sum energy:
         double sum_tower_energy = 0;
         for(int tower_i = 0; tower_i < (*t.towers_e).size(); tower_i++){
            if((*t.towers_clusterID)[tower_i] == (*t.clusters_cid1)[a]){
               sum_tower_energy += (*t.towers_e)[tower_i];
            }
         }

         // derived quantities for cluster a
         double x_a = (*t.clusters_x1)[a];
         double y_a = (*t.clusters_y1)[a];
         double z_a = (*t.clusters_z1)[a];
        
         //cout << "ecore_a: " << ecore_a << endl;

         double trkLength_a = sqrt(x_a*x_a + y_a*y_a + (z_a)*(z_a));

         //if(trkLength_a == 0) continue;
         double x_unit_a = x_a/trkLength_a;
         double y_unit_a = y_a/trkLength_a;
         double z_unit_a = (z_a)/trkLength_a;

         double px_a = ecore_a * x_unit_a;
         double py_a = ecore_a * y_unit_a;
         double pz_a = ecore_a * z_unit_a; 
         double pt_a = sqrt(px_a*px_a + py_a*py_a);
         pT_a->Fill(pt_a);
         for(int b = a+1; b < (*t.clusters_ecore1).size(); b++){ // looping over clusters

            int sec_b = (*t.clusters_sec1)[b];
            if(sec_b > 5) continue;
            double ecore_b = (*t.clusters_ecore1)[b];
            if(ecore_b < 0.3 ) continue;
            double chi2_b = (*t.clusters_chi21)[b];
            if(chi2_b > 3) continue; 
            if(sec_a != sec_b) continue;

            // ADD APPLY DH MAPS HERE
            if(dhmaps[sec_b]->GetBinContent(dhmaps[sec_b]->FindBin((*t.clusters_iz)[b],(*t.clusters_iy)[b])) > 0.1) continue; // 0.1 just to account for possible double precision

            double dist_gea_calc_b = sqrt( 
               pow((*t.clusters_impx11)[b]-(*t.clusters_x1)[b],2) +  
               pow((*t.clusters_impy11)[b]-(*t.clusters_y1)[b],2) +  
               pow((*t.clusters_impz11)[b]-(*t.clusters_z1)[b],2)
            );

            // derived quantitites for cluster b
            double x_b = (*t.clusters_x1)[b];
            double y_b = (*t.clusters_y1)[b];
            double z_b = (*t.clusters_z1)[b];

            double trkLength_b = sqrt(x_b*x_b + y_b*y_b + (z_b)*(z_b));
            //if(trkLength_b == 0) continue;
            double x_unit_b = x_b/trkLength_b;
            double y_unit_b = y_b/trkLength_b;
            double z_unit_b = (z_b)/trkLength_b;

            double px_b = ecore_b * x_unit_b;
            double py_b = ecore_b * y_unit_b;
            double pz_b = ecore_b * z_unit_b; 

            double pt_b = sqrt(px_b*px_b + py_b*py_b);
            //pair quantities and cuts
            //if(ecore_a + ecore_b == 0) continue;
            double asym = fabs( ecore_a - ecore_b ) / ( ecore_a + ecore_b );

            double px = ecore_a * (x_unit_a) + ecore_b * (x_unit_b);
            double py = ecore_a * (y_unit_a) + ecore_b * (y_unit_b);
            double pz = ecore_a * (z_unit_a) + ecore_b * (z_unit_b);
            double pt = sqrt( (px * px) + (py * py) );

            //double pt = sqrt( pow(ecore_a*x_a/trkLength_a + ecore_b*x_b/trkLength_b,2) 
            //                + pow(ecore_a*y_a/trkLength_a + ecore_b*y_b/trkLength_b,2) );

            double deltaR = sqrt( pow((x_a - x_b),2) + pow((y_a - y_b),2) + pow((z_a - z_b),2) );
            if(deltaR < 8) continue; 

            double cosine = x_unit_a*x_unit_b + y_unit_a*y_unit_b + z_unit_a*z_unit_b;
            //if(2 * ecore_a * ecore_b * ( 1.0 - cosine ) < 0) continue;
            double inv_mass = sqrt( 2 * ecore_a * ecore_b * ( 1.0 - cosine ) );
            if(inv_mass < 0.0425) continue;
            if (asym > 0.8) continue;

            filler_mass_vs_pt[0] = sec_a;
            filler_mass_vs_pt[1] = inv_mass;
            filler_mass_vs_pt[2] = pt;
            mass_vs_pt->Fill(filler_mass_vs_pt,g_weights->Eval(true_pt));
            mass_vs_pt_uw->Fill(filler_mass_vs_pt);

            //mass_vs_pt_sec[sec_a]->Fill(inv_mass,pt,g_weights->Eval(true_pt));
            //mass_vs_pt_sec_uw[sec_a]->Fill(inv_mass,pt);
            
         } // END OF B LOOP
      } // END OF A LOOP
   } // END OF EVENT LOOP

   fout->cd();
   truepT->Write();
   pT_a->Write();
   dist_gea->Write();
   mass_vs_pt->Write();
   mass_vs_pt_uw->Write();
   iyiz_a->Write();
   x11_y11->Write();

   fout->Close();
}

void pi0tree::Loop()
{

   // for reasons beyond me, this does not work...
   // work around functions above

//   In a ROOT session, you can do:
//      Root > .L pi0tree.C
//      Root > pi0tree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      fChain->GetEntry(jentry);
      
   }
}
