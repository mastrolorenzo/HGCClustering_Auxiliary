#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

#include "TROOT.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "testMacro.h" 
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

/* parameters definition */
bool verbose=false;
const unsigned nDim = 52;
const unsigned particleId = 11;
double C3d_minEta = 1.5;
double C3d_maxEta = 2.9;
double gen_minEta = 1.7;
double gen_maxEta = 2.7;
double dR_match = 0.3;

void ProduceMatrix(TChain *chain){

    TFile*fout = new TFile("output_file_with_matrix_histograms.root","RECREATE");

    TH2D* h_M = new TH2D("h_M","h_M",nDim,0,nDim,nDim,0,nDim);
    TH1D* h_N = new TH1D("h_N","h_N",nDim,0,nDim);
    TH1D* h_v = new TH1D("h_v","h_v",nDim,0,nDim);
    TH2D* h_N2 = new TH2D("h_N2","h_N2",nDim,0,nDim,nDim,0,nDim);
    
    std::vector<float> *gen_pt_ = 0;
    std::vector<float> *gen_eta_ = 0;
    std::vector<float> *gen_phi_ = 0;
    std::vector<float> *gen_energy_ = 0;
    std::vector<int> *gen_status_ = 0;
    std::vector<int> *gen_id_ = 0;       
    std::vector<float> *c3d_pt_ = 0;
    std::vector<float> *c3d_eta_ = 0;
    std::vector<float> *c3d_phi_ = 0;
    std::vector<std::vector<float>> *cl3d_cl_mipPt_ = 0;
    std::vector<std::vector<int>> *cl3d_cl_layer_ = 0;
 
    chain->SetBranchAddress("gen_pt",&gen_pt_); 
    chain->SetBranchAddress("gen_eta",&gen_eta_); 
    chain->SetBranchAddress("gen_phi",&gen_phi_); 
    chain->SetBranchAddress("gen_energy",&gen_energy_); 
    chain->SetBranchAddress("gen_status",&gen_status_); 
    chain->SetBranchAddress("gen_id",&gen_id_); 
    chain->SetBranchAddress("cl3d_pt",&c3d_pt_); 
    chain->SetBranchAddress("cl3d_eta",&c3d_eta_); 
    chain->SetBranchAddress("cl3d_phi",&c3d_phi_); 
    chain->SetBranchAddress("cl3d_cl_mipPt",&cl3d_cl_mipPt_); 
    chain->SetBranchAddress("cl3d_cl_layer",&cl3d_cl_layer_); 

    TVectorD V(nDim);
    TVectorD a(nDim);
    TMatrixDSym M(nDim);
    TMatrixD Minv(nDim,nDim);
    double tmpVec[nDim];
    double tmpMtx[nDim][nDim];
    int N=0;

    if(verbose) std::cout << "N events: " << chain->GetEntries() << std::endl;
    
    /* inizialize all matrix and vector entries to 0. */
    for(unsigned i(0); i<nDim; i++){
        V(i)=0.;
        tmpVec[i]=0.;
        for(unsigned j(0); j<nDim; j++){
            M(i,j)=0.;
            tmpMtx[i][j]=0.;
        }
    }

    for(int entry = 0; entry<chain->GetEntries(); entry++){
        chain->GetEntry(entry);
        
        if(verbose) std::cout << "processing event " << entry << std::endl;

        /* select the C3d that pass the selection */
        std::vector<int> goodC3d_idx;
        for(unsigned i_c3d(0); i_c3d<(*c3d_pt_).size(); i_c3d++){
            if( fabs( (*c3d_eta_)[i_c3d] ) > C3d_minEta 
                &&  fabs( (*c3d_eta_)[i_c3d] ) < C3d_maxEta 
                && (*c3d_pt_)[i_c3d] > 0 )
            {
                goodC3d_idx.push_back(i_c3d);
            }
        }

        /* loop over the gen particle */
        for(unsigned i_gen(0); i_gen<(*gen_pt_).size(); i_gen++){
            if( fabs( (*gen_eta_)[i_gen] ) > gen_minEta && fabs( (*gen_eta_)[i_gen] ) < gen_maxEta && (*gen_pt_)[i_gen] > 10 
                && fabs( (*gen_id_)[i_gen] ) == particleId && (*gen_status_)[i_gen] == 1 ){

                bool hasMatched = false;
                double pt_cand = 0.;
            
                /* loop over the good 3D-cluster */
                if( goodC3d_idx.size() <= 0 )
                    continue;
                
                /* look for all the C3d mathcing with the gen */
                std::vector<int> bestC3dIdx;
                int bestC3dIdx_MAX = -1;
                for( unsigned i(0); i< goodC3d_idx.size(); i++){
                    int i_c3d = goodC3d_idx[i];   
                    double dR = deltaR( (*gen_eta_)[i_gen], (*gen_phi_)[i_gen], (*c3d_eta_)[i_c3d], (*c3d_phi_)[i_c3d] );                
                    if(dR<dR_match){
                        hasMatched = true;
                        if( (*c3d_pt_)[i_c3d] > pt_cand  ){
                            pt_cand = (*c3d_pt_)[i_c3d];
                            bestC3dIdx_MAX = i_c3d;
                        }
                        bestC3dIdx.push_back(i_c3d);
                    }
                }
                if(bestC3dIdx.size()<=0) continue;
                if(verbose){
                    std::cout << "gen particle = " << (*gen_pt_)[i_gen] 
                              << ", ene= " << (*gen_energy_)[i_gen]  
                              << "  maximum C3d: " <<  (*c3d_pt_)[bestC3dIdx_MAX]<<  std::endl;
                }
                
                /* loop over all the C3d in the matching cone */
                std::map<int,double> supercluster;
                for(unsigned i_cone(0); i_cone<bestC3dIdx.size(); i_cone++){
                    if(verbose) std::cout << "c3d-"<< i_cone << "  pT = " << (*c3d_pt_)[bestC3dIdx.at(i_cone)] << std::endl;
                    
                    std::map<int, double> clusterMap_tmp;
                    
                    /*---- Start computation to add mipT in the same layer for multiple clusters ----*/                        
                    /* loop over the C2ds of each secondaries */
                    for(unsigned j=0; j<(*cl3d_cl_mipPt_)[bestC3dIdx.at(i_cone)].size(); j++){
                        if(verbose){ 
                            std::cout << " layer: " << (int)(*cl3d_cl_layer_)[bestC3dIdx.at(i_cone)][j] 
                                      << " energy:  " << (*cl3d_cl_mipPt_)[bestC3dIdx.at(i_cone)][j] << std::endl;
                        }
                        
                        supercluster.insert( make_pair( (int)(*cl3d_cl_layer_)[bestC3dIdx.at(i_cone)][j], 
                                                        supercluster[(int)(*cl3d_cl_layer_)[bestC3dIdx.at(i_cone)][j]] 
                                                        += (*cl3d_cl_mipPt_)[bestC3dIdx.at(i_cone)][j] ) );
                    }
                                        
                    if(verbose) std::cout << "----------------" << std::endl;

                }
                
                if(verbose){
                    std::cout << " " << std::endl;
                    for(std::map<int, double>::const_iterator it=supercluster.begin(); it!=supercluster.end(); ++it)
                        std::cout << it->first << " => " << it->second << '\n';
                }

                N++;

                /* Fill the histograms of the matrix and vector */
                for(auto it = supercluster.begin(); it != supercluster.end(); ++it){
                    h_v->Fill( it->first - 1, (it->second) * (*gen_pt_)[i_gen] );
                    tmpVec[ it->first - 1 ] += (it->second)*(*gen_pt_)[i_gen];
                    for(auto jt = supercluster.begin(); jt != supercluster.end(); ++jt){
                        h_M->Fill( it->first - 1, jt->first -1, (it->second) * (jt->second) );
                        tmpMtx[it->first - 1][jt->first -1] += ( (it->second) * (jt->second) );                        
                    }
                }
            }
        }
    }

    for(unsigned i(0); i<nDim; i++){
        h_N->Fill(i,N);
        for(unsigned j(0); j<nDim; j++){
            h_N2->Fill(i,j,N);
        }
    }

    fout->cd();
    h_M->Write();
    h_v->Write();
    h_N->Write();
    h_N2->Write();

}


void InversionMIPt_Cluster(TString fileName){
    TFile*f = new TFile(fileName,"READ");

    const int nDim = 52;

    TH2D* h_de_i_de_j = (TH2D*)f->Get("h_M");
    TH1D* h_N = (TH1D*)f->Get("h_N");
    TH1D* h_de_j_Te = (TH1D*)f->Get("h_v");
    TH2D* h_N2 = (TH2D*)f->Get("h_N2");
    TH2F*h_M = new TH2F("h_MM","h_MM",nDim,0,nDim,nDim,0,nDim);
    TH2F*h_Minv = new TH2F("h_Minv","h_Minv",nDim,0,nDim,nDim,0,nDim);
    TH2F*h_I = new TH2F("h_I","h_I",nDim,0,nDim,nDim,0,nDim);
    TH1F*h_a = new TH1F("h_a","h_a",nDim,0,nDim);
    TH1F*h_v = new TH1F("h_vv","h_vv",nDim,0,nDim);
  
    TMatrixDSym m(nDim);
    TMatrixDSym m_original(nDim);  
    TVectorD V(nDim);

    /* define the trigger layers set */
    int k[nDim]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
  
    /* From the histogram fill the matrix and the vector */
    for(int i=0; i<nDim; i++){
        V(i)= h_de_j_Te->GetBinContent(k[i]) ;
        h_v->SetBinContent(i+1, V(i) );  
       
        for(int j=0; j<nDim; j++){
            m(i,j) = h_de_i_de_j->GetBinContent(k[i],k[j]);            
            m_original(i,j) = h_de_i_de_j->GetBinContent(k[i],k[j]);            
            h_M->SetBinContent(i+1, j+1, m(i,j));
        }

    }
    if(verbose){
        for(int i(0); i< nDim; i++){
            std::cout << "vector element (i): = "<< V(i)<< std::endl;
            for(int j(0); j< nDim; j++){
                std::cout << "matrix element (i,j): "<< m(i,j) <<std::endl; 
            }
        }
        std::cout << "Matrix Determinant = "<< m.Determinant() << std::endl;   
    }
    
    /* Invert the matrix */ 
    TMatrixD M = m.Invert();

    /* Make the cross-check */
    TMatrixD I = m_original*M;    
    for(int i=0; i<nDim; i++){
        for(int j=0; j<nDim; j++){
            h_Minv->SetBinContent(i+1, j+1, M(i,j));
            if(i==j) std::cout << I(i,j)<< std::endl; 
            h_I->SetBinContent(i+1, j+1, I(i,j));
        }
    }

    /* Print the coefficients in a .txt file */
    ofstream outTextFile;
    outTextFile.open ("coefficients_file.txt");
    TVectorD a(nDim);
    a = M * V;
    for(int i=0; i<nDim; i++){
        if(i<nDim-1) outTextFile << a(i) << ",\n";
        else  outTextFile << a(i) << "\n";
        h_a->SetBinContent(i+1, a(i));       
    }
    outTextFile.close();
  
    /* Produce and store all the relevant plots */
    TCanvas*c = new TCanvas();
    h_M->Draw("COLZ");

    TCanvas*c1 = new TCanvas();
    h_Minv->Draw("COLZ");

    TCanvas*c2 = new TCanvas();
    h_a->Draw();
    TLine *l = new TLine(0, 0, nDim, 0);
    l->SetLineColor(2);
    l->Draw("same");

    TCanvas*c4 = new TCanvas();
    h_I->Draw("COLZ");

    TCanvas*c5 = new TCanvas();
    h_v->Draw();
 
    TFile*fout = new TFile("histoFile_MatrixInversion.root","RECREATE");
    h_M->Write();
    h_Minv->Write();
    h_I->Write();
    h_v->Write();
    h_a->Write();

    return;
}




void InvertMatrix(){
    ifstream fileList;
    
    /* Specify the file list */
    fileList.open("/afs/cern.ch/work/l/lmastrol/Branch/branchTest/test/HGC_L1Calib_CMSSW/ClusterDev/performancePR/HGCTPGPerformance/list_ElePt10to150_afterCalib_NewC3dPos.txt");

    TString file = "";
    TChain* chain = new TChain("hgcalTriggerNtuplizer/HGCalTriggerNtuple");
    while(1){
        fileList >> file;
        if (!fileList.good()) break;
        chain->Add("root://cms-xrd-global.cern.ch/"+file);
    }

    ProduceMatrix(chain);
    InversionMIPt_Cluster("./histoFile_MatrixInversion.root");

    return;
}
