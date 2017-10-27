/* Simple root macro that project the number of 2D-clusters and C3d into histograms - L. Mastrolorenzo */

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"

void produceBWC2d( TChain *chain ){

    TFile*output = new TFile("histo_TT_NNC2d_SE5TE2_dRC3d_withBH_highStat.root","RECREATE");
    TH2D*h = new TH2D("h","bandwidthPlot",52, 0, 53, 500, 0, 500 );
    std::cout << "Total number of events: " << chain->GetEntries() << std::endl;
    chain->Project("h","cls_NxLayer:cls_l");
    h->Draw();
    h->Write();

    return; 

}

void produceBWC3d( TChain *chain ){

    TFile*output = new TFile("histo_TT_NNC2d_SE5TE2_dRC3d_withBH_highStat.root","RECREATE");
    TH2D*h = new TH2D("h","bandwidthPlot",52, 0, 53, 500, 0, 500 );
    std::cout << "Total number of events: " << chain->GetEntries() << std::endl;
    chain->Project("h","cl3ds_NxEndcap");
    h->Draw();
    h->Write();
    
    return; 

}

void produceSeedsBW(){

    TFile*output = new TFile("histo_TT_NNC2d_SE5TE2_dRC3d_withBH_highStat.root","RECREATE");

    TH1D* h_InclusiveNseed = new TH1D("h_InclusiveNseed","Number of seeds inside C2d", 15, 0, 16);
    TH2D* h_InclusiveNseedXlayer = new TH2D("h_InclusiveNseedXlayer","Number of seeds inside C2d vs layer", 52, 0, 53, 15, 0, 16);  
    
    chain.Project("h_InclusiveNseed","cls_Nseeds");
    chain.Project("h_InclusiveNseedXlayer","clsTc_NseedsXlayer:clsTc_l");
    
    TFile* out = new TFile("outNseedsHistos.root","RECREATE");
    
    h_InclusiveNseed->Write();
    h_InclusiveNseedXlayer->Write();

return;
}


void produceHisto(){
    ifstream fileList;
    fileList.open("list_TTbarWithBH_PU200_93X_GeoV8_TDR.txt");
    TString file = "";
    TChain chain("hgcalTriggerNtuplizer/HGCalTriggerNtuple");
    
    while(1){
        fileList >> file;
        if (!fileList.good()) break;
        std::cout << "file: " << file << std::endl;
        chain.Add("root://cms-xrd-global.cern.ch/"+file);
    }

    produceBWC2d(&chain);
    produceBWC3d(&chain);
    produceSeedsBW(&chain);

    return;
}

