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

//TString fin = "../ntuple_TT_SE5TE2dR3.root";

void getHisto(){
ifstream fileList;
    fileList.open("./ntupleList_TTbar_tc999.txt");
//fileList.open("ntupleList_QCD_tc999.txt");
TString file = "";

TChain chain("hgcalTriggerNtuplizer/HGCalTriggerNtuple");
//while(1){
for(int i=0; i<2; i++){
    fileList >> file;
    if (!fileList.good()) break;
    std::cout << "file: " << file << std::endl;
    chain.Add("root://cms-xrd-global.cern.ch/"+file);
}

//TFile*f = new TFile(fin,"READ");
//TDirectoryFile* d = (TDirectoryFile*)f->Get("hgcalTriggerNtuplizer");
//TTree*t = (TTree*)d->Get("HGCalTriggerNtuple");

TH1D* h_InclusiveNseed = new TH1D("h_InclusiveNseed","Number of seeds inside C2d", 15, 0, 15);
TH2D* h_InclusiveNseedXlayer = new TH2D("h_InclusiveNseedXlayer","Number of seeds inside C2d vs layer", 52, 0, 52, 15, 0, 15);  

chain.Project("h_InclusiveNseed","cls_Nseeds");
chain.Project("h_InclusiveNseedXlayer","clsTc_NseedsXlayer:clsTc_l");
chain.Draw("clsTc_NseedsXlayer");
//TFile* out = new TFile("outNseedsHistos.root","RECREATE");
//
//h_InclusiveNseed->Write();
//h_InclusiveNseedXlayer->Write();

return;
}
