#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TPad.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"

double dRMatch = 0.5;


double deltaPhi( double phi1, double phi2) {
        
    double dPhi(phi1-phi2);
    double pi(acos(-1.0));
    if     (dPhi<=-pi) dPhi+=2.0*pi;
    else if(dPhi> pi) dPhi-=2.0*pi;
        
    return dPhi;
}
    

double deltaEta(double eta1, double eta2){
    double dEta = (eta1-eta2);
    return dEta;
}

double deltaR(double eta1, double eta2, double phi1, double phi2) {
    double dEta = deltaEta(eta1, eta2);
    double dPhi = deltaPhi(phi1, phi2);
    return sqrt(dEta*dEta+dPhi*dPhi);
}



void plotBW( TChain *chain ){

    chain->Print();

    int run_;
    int event_; 
    int lumi_;
    
    vector<float> *gen_pt_ = 0;
    vector<float> *gen_eta_ = 0;
    vector<float> *gen_phi_ = 0;
    vector<float> *gen_energy_ = 0;
    vector<int> *gen_status_ = 0;
    vector<int> *gen_id_ = 0;

    vector<float> *c3d_pt_ = 0;
    vector<float> *c3d_eta_ = 0;
    vector<float> *c3d_phi_ = 0;
    vector<float> *c3d_energy_ = 0;
    
    chain->SetBranchAddress("run",&run_);
    chain->SetBranchAddress("event",&event_);
    chain->SetBranchAddress("lumi",&lumi_);
    
    chain->SetBranchAddress("gen_pt",&gen_pt_);
    chain->SetBranchAddress("gen_eta",&gen_eta_);
    chain->SetBranchAddress("gen_phi",&gen_phi_);
    chain->SetBranchAddress("gen_energy",&gen_energy_);
    chain->SetBranchAddress("gen_status",&gen_status_);
    chain->SetBranchAddress("gen_id",&gen_id_);
    
    chain->SetBranchAddress("cl3d_pt",&c3d_pt_);
    chain->SetBranchAddress("cl3d_eta",&c3d_eta_);
    chain->SetBranchAddress("cl3d_phi",&c3d_phi_);
    chain->SetBranchAddress("cl3d_energy",&c3d_energy_);

    TH1D*h_resoPt = new TH1D("resoPt","Pt response",100, 0, 2);
    TH1D*h_resoEta = new TH1D("resoEta","Eta response",100, -0.5, 0.5);
    TH1D*h_resoPhi = new TH1D("resoPhi","Phi response",100, -0.5, 0.5);
 
//    for(int entry=0; entry<10; entry++){
    for(int entry=0; entry<chain->GetEntries(); entry++){
        chain->GetEntry(entry);
        std::cout<<"event: " << event_ << " gen-part size: "<<(*gen_pt_).size() <<std::endl;
  
        for(unsigned i_gen = 0; i_gen<(*gen_pt_).size(); i_gen++){
 
            if( fabs( (*gen_eta_)[i_gen] ) > 1.47 &&  fabs( (*gen_eta_)[i_gen] ) < 3.0 && (*gen_pt_)[i_gen]>7 && abs( (*gen_id_)[i_gen] ) == 11 && abs( (*gen_status_)[i_gen] ) == 1 ){
               
                std::cout << "\n";
                double pt_cand = -1.;
                double eta_cand = -100.;
                double phi_cand = -100.;
                for(unsigned i_c3d = 0; i_c3d<(*c3d_pt_).size(); i_c3d++){
 
                    double dR = deltaR( (*gen_eta_)[i_gen], (*c3d_eta_)[i_c3d], (*gen_phi_)[i_gen], (*c3d_phi_)[i_c3d] );                
                    if(dR<dRMatch && fabs( (*c3d_eta_)[i_c3d] ) > 1.47 &&  fabs( (*c3d_eta_)[i_c3d] ) < 3.0 && (*c3d_pt_)[i_c3d]>2 ){
//                        std::cout << "gen pt = " << (*gen_pt_)[i_gen] << " cluster pt " <<(*c3d_pt_)[i_c3d] << std::endl;
                      
                        if( (*c3d_pt_)[i_c3d]>pt_cand ){
                            pt_cand = (*c3d_pt_)[i_c3d];
                            eta_cand = (*c3d_eta_)[i_c3d];
                            phi_cand = (*c3d_phi_)[i_c3d];
                        }

                    }                    
                }
                std::cout << " ==================== > gen pt,eta,phi = " << (*gen_pt_)[i_gen] << ", "<< (*gen_eta_)[i_gen] << ", "<< (*gen_phi_)[i_gen]  
                          << " L1 candidate  pt,eta,phi " << pt_cand << ", "<< eta_cand << ", "<< phi_cand <<std::endl;
                h_resoPt->Fill( pt_cand / (*gen_pt_)[i_gen] );
                h_resoEta->Fill( eta_cand - (*gen_eta_)[i_gen] );
                h_resoPhi->Fill( phi_cand - (*gen_phi_)[i_gen] );

            }

        }
    }

    TFile*output = new TFile("response_highStat.root","RECREATE");
    h_resoPt->Write();
    h_resoEta->Write();
    h_resoPhi->Write();

    return; 

}


void plotAllResponse(){
    ifstream fileList;
    fileList.open("ntupleList_ZEE_tc999.txt");
    TString file = "";

    TChain chain("hgcalTriggerNtuplizer/HGCalTriggerNtuple");
    while(1){
        fileList >> file;
        if (!fileList.good()) break;
        std::cout << "file: " << file << std::endl;
        chain.Add("root://cms-xrd-global.cern.ch/"+file);
    }

    plotBW(&chain);
//    plotBW("ntuple_ZEE.root");
    
    return;
}

