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

//TString dir = "";


void plotBW( TChain *chain ){

    chain->Print();

    TH2D*h = new TH2D("h","bandwidthPlot",52, 0, 53, 500, 0, 500 );
    chain->Project("h","cls_NxLayer:cls_l");
    h->Draw();
//    TFile*output = new TFile("histo_TT_SE5TE2dR3_highStat.root","RECREATE");
    TFile*output = new TFile("histo_QCD_SE5TE2dR3_highStat.root","RECREATE");
    h->Write();

    return; 

}


void plotHisto(){
    ifstream fileList;
//    fileList.open("./ntupleList_TTbar_tc999.txt");
    fileList.open("ntupleList_QCD_tc999.txt");
    TString file = "";

    TChain chain("hgcalTriggerNtuplizer/HGCalTriggerNtuple");
    while(1){
        fileList >> file;
        if (!fileList.good()) break;
        std::cout << "file: " << file << std::endl;
        chain.Add("root://cms-xrd-global.cern.ch/"+file);
    }

    plotBW(&chain);
    
    return;
}

