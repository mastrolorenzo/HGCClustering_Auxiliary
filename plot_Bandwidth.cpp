#include <iostream>
#include <cstdlib>
#include <numeric>
#include <vector>

#include "TString.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"

#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TApplication.h"

using namespace std;

bool test = false;

TString linkSpeed = "16.4";
TString encoding = "64b/66b";

//int nbit_x_cl = 128; /* in bits */

double Bx = 25.0; /* in ns */
double linkBw = 16.4; /* raw link bandwidth in Gb/s */
double linkEnc = 64./66; 
int nTwMaps = 1008; /* number of tower maps */
int bitXtwMap = 8; /* bits allocated per tower map */
//int tMux = 18; /* time multiplexed */
float linkMarginFactor = 0.15;

float maxYaxis = 500.;
int ncl_max = maxYaxis;

int plotBW( TString inFileName, bool* maskTrgLayers, 
             int tMux=18, int nbit_x_cl=128 ){

    /* open file */
    TFile* f = new TFile(inFileName, "READ");
    
    /* get the histo */
    TH2D* h = (TH2D*)f->Get("h");
    h->SetTitle("");

    /* set style*/
    gStyle->SetTitleFontSize(0.035);
    gStyle->SetOptStat(0);
   
    /* define detector regions */
    TLine *EE = new TLine(29, 0, 
                          29, maxYaxis);
    EE->SetLineWidth(3);
    

    TLine *FH = new TLine(41, 0, 
                          41, maxYaxis);
    FH->SetLineWidth(3);
    
    TPaveText *ptEE = new TPaveText(1, 250, 
                                    28, 260 );    
    ptEE->AddText("EE");
    ptEE->SetBorderSize(0);
    ptEE->SetFillStyle(0);


    TPaveText *ptFH = new TPaveText(28, 250, 
                                    41, 260 );    
    ptFH->AddText("FH");
    ptFH->SetBorderSize(0);
    ptFH->SetFillStyle(0);


    TPaveText *ptBH = new TPaveText(41, 250, 
                                    52, 260 );    
    ptBH->AddText("BH");
    ptBH->SetBorderSize(0);
    ptBH->SetFillStyle(0);


    /* 
       evaluation of the bandwidth scale 
       bandwidth =  [ ncl_max * (nbit/cl) ] / Bx  Gb/s 
    */
    //original float nMinBW = ( linkBw * linkEnc * tMux);
    float minBw = ( nTwMaps * bitXtwMap / Bx ) * (1+linkMarginFactor); // Gb/s
    float maxBw = ( ncl_max * nbit_x_cl / Bx ) + minBw; // Gb/s

    /* maxNlinks = max_bw / ( link_bandwidth * encoding ) + nTowerMaps * 8bit / Bx */
    float nMinLinks = minBw / ( linkBw * linkEnc ) ;
    float nMaxLinks = maxBw / ( linkBw * linkEnc ) ; 
  
    /* drawing the new Axises */
    int axisOffset=18;
    TGaxis *axisbw = new TGaxis(41+axisOffset, 0, 41+axisOffset, maxYaxis, 
                                minBw/100, maxBw/100, 510, "+L" );
    axisbw->SetLineColor(kBlue);
    axisbw->SetLabelColor(kBlue);
    axisbw->SetTextColor(kBlue);
    axisbw->CenterTitle(true);
    axisbw->SetLabelFont(42);
    axisbw->SetTextFont(42);
    axisbw->SetTitle("Bandwidth [ 100 Gb/s ]");
    axisbw->SetLabelSize(0.02);
    axisbw->SetTitleSize(0.02);
    axisbw->SetTickLength(0.01);

    //axisbw->SetTitleOffset(1.3);

     
    TGaxis *axislink = new TGaxis(47+axisOffset, 0, 47+axisOffset, maxYaxis, 
                                  nMinLinks, nMaxLinks, 510, "+L" );
    axislink->SetLineColor(kMagenta+1);
    axislink->SetLabelColor(kMagenta+1);
    axislink->SetTextColor(kMagenta+1);
    axislink->CenterTitle(true);
    axislink->SetLabelFont(42);
    axislink->SetTextFont(42);
    axislink->SetTitle("Number of "+linkSpeed+" Gb/s links "+encoding);
    axislink->SetLabelSize(0.02);
    axislink->SetTitleSize(0.02);
    axislink->SetTickLength(0.01);
    //axislink->SetTitleOffset(1.3);
    

    TGaxis *axislinkC2dOut = new TGaxis(53, 0, 53, maxYaxis,
                                        nMinLinks/tMux, nMaxLinks/tMux, 510, "+L" );
    axislinkC2dOut->SetLineColor(kRed);
    axislinkC2dOut->SetLabelColor(kRed);
    axislinkC2dOut->SetTextColor(kRed);
    axislinkC2dOut->CenterTitle(true);
    axislinkC2dOut->SetLabelFont(42);
    axislinkC2dOut->SetTextFont(42);
    axislinkC2dOut->SetTitle("Number of "+linkSpeed+" Gb/s links "+encoding+" out from C2d-layer / board");
    axislinkC2dOut->SetLabelSize(0.03);
    axislinkC2dOut->SetTitleSize(0.03);
    axislinkC2dOut->SetTitleOffset(0.7);
    axislinkC2dOut->SetTickLength(0.01);
  

    TGraph*g0 = new TGraph();
    TGraph*g1 = new TGraph();
    TGraph*g2 = new TGraph();
    TGraph*g3 = new TGraph();


    float bins[] = { 1,2,3,4,
                     5,6,7,8,
                     9,10,11,12,
                     13,14,15,16,
                     17,18,19,20,
                     21,22,23,24,
                     25,26,27,28,
                     29,30,31,32,
                     33,34,35,36,
                     37,38,39,40,
                     41,42,43,
                     45,47,
                     49,51, 
                     53
    };
    //TH1D* hMaxY = new TH1D("hMaxY","hMaxY", 52, 1, 53);
    TH1D* hMaxY = new TH1D("hMaxY","hMaxY", 47, bins );
    
 //    double binsX[12];
 //    for(int i(0); 52; i++){
 //        binsX[i] = i;
 //    }
 //    TH1D*h1 = new TH1D("h1","h1", 11, binsX);

     for(int i_binx(1); i_binx<h->GetNbinsX(); i_binx++ ){
         TH1D* h_tmp = h->ProjectionY( Form("tmp_bin%d",i_binx), i_binx, i_binx );
         int max_x = 0; 
         double bin_centre = h->GetXaxis()->GetBinCenter(i_binx);
         for(int i_tmp(1); i_tmp<h_tmp->GetNbinsX(); i_tmp++)
             if( h_tmp->GetBinContent(i_tmp)>0 )
                 max_x = i_tmp;

         g0->SetPoint(g0->GetN(), bin_centre, max_x);
         g1->SetPoint(g1->GetN(), bin_centre, max_x + 0.05*max_x);
         g2->SetPoint(g2->GetN(), bin_centre, max_x + 0.15*max_x);
         g3->SetPoint(g3->GetN(), bin_centre, max_x + 0.20*max_x);

     }

     g0->SetLineColor(kGreen+0);
     g1->SetLineColor(kGreen+1);
     g2->SetLineColor(kGreen+2);
     g3->SetLineColor(kGreen+3);

     g0->SetMarkerColor(kGreen+0);
     g1->SetMarkerColor(kGreen+1);
     g2->SetMarkerColor(kGreen+2);
     g3->SetMarkerColor(kGreen+3);

     /*g0->Draw("same pl");
     g1->Draw("same pl");
     g2->Draw("same pl");
     g3->Draw("same pl");
     */

     double* gY = g0->GetY();
     double* gX = g0->GetX();
     int i_vbin=0;
     int max = 0;

 //    for(int ig=0; ig<g0->GetN(); ig++){
     for(int ig=0; ig<41; ig++){
         //if(gX[ig]<binsX[i_vbin]){
         if(gX[ig]<hMaxY->GetBinContent(i_vbin) ){
             if(max<gY[ig]){
                 max = gY[ig];
             }
         }
         else{
             hMaxY->SetBinContent(i_vbin, max);
             max = gY[ig];
             i_vbin++;
         }
     }
     hMaxY->SetBinContent(41, hMaxY->GetBinContent(40) );
     hMaxY->SetBinContent(42, hMaxY->GetBinContent(40) );
     hMaxY->SetBinContent(43, 2*hMaxY->GetBinContent(40) );
     hMaxY->SetBinContent(44, 2*hMaxY->GetBinContent(40) );
     hMaxY->SetBinContent(45, 2*hMaxY->GetBinContent(40) );
     hMaxY->SetBinContent(46, 2*hMaxY->GetBinContent(40) );
     hMaxY->SetBinContent(47, 2*hMaxY->GetBinContent(40) );

     hMaxY->SetLineColor(kBlack);
     TH1D* hMaxYx05 = (TH1D*)hMaxY->Clone();
     TH1D* hMaxYx15 = (TH1D*)hMaxY->Clone();
     TH1D* hMaxYx20 = (TH1D*)hMaxY->Clone();  
//     TH1D* hNlinksPerLayer = new TH1D("hNlinksPerLayer","", 52, 0, 52);
     TH1D* hNlinksPerLayer = new TH1D("hNlinksPerLayer","", 47, bins);

/* defining the main histo style and Draw */
    hMaxYx15->GetYaxis()->SetTitle("Number of 2D-Clusters");
    hMaxYx15->GetYaxis()->SetTitleOffset(1.2);
    hMaxYx15->GetXaxis()->SetTitle("Number of Trigger Links");
    hMaxYx15->GetYaxis()->CenterTitle();
    hMaxYx15->GetXaxis()->CenterTitle();

    hMaxYx15->GetYaxis()->SetTickLength(0.01);
    hMaxYx15->GetXaxis()->SetTickLength(0.01);

    hMaxYx05->SetTitle("");
    hMaxYx15->SetTitle("");
    hMaxYx20->SetTitle("");
    
    hMaxYx05->Scale(1.05);
    hMaxYx15->Scale(1.15);
    hMaxYx20->Scale(1.20);
    
     int InTotLinks_SingleFPGA_layer2 = 0;
     char *label[47];
     float TMlinkBwTMux = nMinLinks / tMux;

     for(int i_layer = 1; i_layer<=47; i_layer++){

         int nLinksTmux=0;
         
         double maxNc2d = hMaxYx15->GetBinContent( hMaxYx15->FindBin(i_layer) );
         
         //std::cout << minBw << " " << TMlinkBwTMux << " " << nMinLinks << std::endl;
         nLinksTmux = std::ceil( TMlinkBwTMux + ( maxNc2d * nbit_x_cl ) / ( Bx * linkBw * linkEnc * tMux ) );
         char *n_linksLabel = Form("%d", nLinksTmux);
         /*std::cout << nLinksTmux 
                   << " " << TMlinkBwTMux + ( maxNc2d * nbit_x_cl ) / ( Bx * linkBw * linkEnc * tMux )
                   << " " << TMlinkBwTMux 
                   << " " << ( maxNc2d * nbit_x_cl ) / ( Bx * linkBw * linkEnc * tMux )
                   << " " << nMaxLinks/tMux 
                   << " " << ( (nLinksTmux - TMlinkBwTMux ) / ( (nMaxLinks-nMinLinks)/tMux ) ) * maxYaxis
                   << std::endl;
	 */
	 std::cout << "TMlinkBwTMux: " << TMlinkBwTMux
	           << " N-c2d " << maxNc2d
		   << " pureC2d: " << ( maxNc2d * nbit_x_cl ) / ( Bx * linkBw * linkEnc * tMux )
		   << " corresp: " << ( (nLinksTmux) / ( (nMaxLinks)/tMux ) ) * maxYaxis
		   << " nLinksTmux: " <<  TMlinkBwTMux + ( maxNc2d * nbit_x_cl ) / ( Bx * linkBw * linkEnc * tMux )
		   << " ceil(nLinksTmux): " << nLinksTmux
		   << " number of c2d correspondent " << ( (nLinksTmux - TMlinkBwTMux ) / ( (nMaxLinks-nMinLinks)/tMux ) ) * maxYaxis
		   << " max linkY: " << nMaxLinks/tMux
		   << " ncl_max: " << ncl_max
		   << std::endl;
	 hNlinksPerLayer->SetBinContent( i_layer, ( (nLinksTmux - TMlinkBwTMux) / ( (nMaxLinks-nMinLinks)/tMux ) ) * maxYaxis );

         label[i_layer] = n_linksLabel;
         InTotLinks_SingleFPGA_layer2 += nLinksTmux;

     }
     std::cout << "FINE" << std::endl;
 

     for(int i=1; i<=hMaxYx15->GetNbinsX(); i++){
         hMaxYx15->GetXaxis()->SetBinLabel( i, label[i] );
     }
     hMaxYx15->GetXaxis()->LabelsOption("u");   

     TGaxis *axisX = new TGaxis(1, maxYaxis, 
                                53, maxYaxis, 
                                1, 53, 
                                510, "-");
     axisX->SetLabelFont(42);
     axisX->SetTextFont(42);
     axisX->SetTitle("Layer");
     axisX->CenterTitle();
     axisX->SetLabelSize(0.03);
     axisX->SetTitleSize(0.03);
     axisX->SetTitleOffset(1.25);
     axisX->SetTickLength(0.01);


     //std::cout << "InTotLinks_SingleFPGA_layer2 = "<<  InTotLinks_SingleFPGA_layer2 << std::endl;
     hMaxYx05->SetLineColor(kBlack);
     hMaxYx15->SetLineColor(kBlack);
     hMaxYx20->SetLineColor(kBlack);
     hNlinksPerLayer->SetLineColor(kRed);

     hMaxYx05->SetLineWidth(2);
     hMaxYx15->SetLineWidth(2);
     hMaxYx20->SetLineWidth(2);
     hNlinksPerLayer->SetLineWidth(2);

     hMaxYx05->SetLineStyle(3);
     hMaxYx15->SetLineStyle(3);
     hMaxYx20->SetLineStyle(3);
     hNlinksPerLayer->SetLineStyle(2);


     TEllipse* e =  new TEllipse(0.3, -6, 0.5, 6.);
     e->SetFillStyle(0);
     e->SetLineWidth(2);
     e->SetLineColor(kGreen+1);

     int n_linksTo3Dboard = 0;

     for(int i_e=1; i_e<47; i_e++){

         if( !maskTrgLayers[i_e-1] )
             continue;

         n_linksTo3Dboard += atoi( label[i_e-1] );

     }

     TPaveText *pt = new TPaveText(54, maxYaxis+10, 
                                   64, maxYaxis+35 );

     pt->AddText( Form("N links out: %d", n_linksTo3Dboard) );
     pt->AddText( Form("Clu 2D size: %d b", nbit_x_cl) );
     pt->SetLineColor(kBlack);
     pt->SetFillStyle(0);

     /* THE LEGEND */

     TLegend*leg = new TLegend(0.25, 0.77, 
                               0.42, 0.87 );
     leg->AddEntry(hMaxY,"Max(Number of C2d)","l");
 //    leg->AddEntry(hMaxYx05,"Max(Number of C2d)+5%","l");
     leg->AddEntry(hMaxYx15,"Max(Number of C2d)+15%","l");
 //    leg->AddEntry(hMaxYx20,"Max(Number of C2d)+20%","l");
     leg->AddEntry(hNlinksPerLayer,"Max link/layer-board","l");

     /************/
     /* Box on BH*/
     TBox* BHbox =  new TBox(41, 0, 53, maxYaxis);
     BHbox->SetFillColor(kGray+1);
     BHbox->SetFillStyle(3004);
     //BHbox->SetHatchesSpacing(2);
     BHbox->SetLineWidth(0);

     TPaveText *ptExtrapolated = new TPaveText(44, 150, 
                                               50, 200 );    
     ptExtrapolated->AddText("Extrapolated");
     ptExtrapolated->SetBorderSize(0);
     ptExtrapolated->SetFillStyle(0);
     ptExtrapolated->SetTextColor(kGray+1);

     /**********************/
     /* Draw what you need */

     TCanvas* c = new TCanvas("c", "c", 1500, 900);
     c->SetRightMargin(0.225);
          
     hMaxYx15->Draw();
     BHbox->Draw("same");
     ptExtrapolated->Draw("same");
     h->Draw("box same");
     hMaxYx15->GetYaxis()->SetRangeUser(0, maxYaxis);
     c->Update();

     hMaxY->Draw("same");

     //hMaxYx05->Draw("same hist");
     //hMaxYx20->Draw("same hist");

     hNlinksPerLayer->Draw("same");

     EE->Draw("same");
     FH->Draw("same");
     ptEE->Draw("same");
     ptFH->Draw("same");
     ptBH->Draw("same");

     axisbw->Draw();
     axislink->Draw();
     axislinkC2dOut->Draw();
     axisX->Draw();

     leg->Draw("same");
     pt->Draw("same");
   

    for(int i_e=1; i_e<=47; i_e++){

        if( !maskTrgLayers[i_e-1] )
            continue;

        TEllipse* e_copy = (TEllipse*)e->Clone();
        if (i_e<43 )
            e_copy->SetX1( e_copy->GetX1() + i_e );        
        else
            e_copy->SetX1( e_copy->GetX1() + 41.5 + 2*(i_e-42) );
        e_copy->Draw("same");

    }

      /* Save the beautiful results */
    TString outFileName;
    outFileName.Append( Form("%d_%db", tMux, nbit_x_cl) );

    TString outPdf = outFileName;
    TString outPng = outFileName;
    TString outC   = outFileName;
    outPdf.Append(".pdf");
    outPng.Append(".png");
    outC.Append(".C");

    if(!test){
        c->SaveAs(outPdf);
        c->SaveAs(outPng);
        c->SaveAs(outC);
    }


    /* make the base histogram */
    TCanvas* cc = new TCanvas("cc", "cc", 1500, 900);
    cc->SetRightMargin(0.225);
    hMaxYx15->GetXaxis()->SetTitle("");
    hMaxYx15->Draw("axis");
    hMaxYx15->GetXaxis()->SetLabelOffset(999);
    hMaxYx15->GetXaxis()->SetLabelSize(0);
    
    h->Draw("colz AH same");
    axisX->Draw();
 
    EE->Draw("same");
    FH->Draw("same");
    ptEE->Draw("same");
    ptFH->Draw("same");
    ptBH->Draw("same");
      
    if(!test){
        cc->SaveAs("baseHisto.pdf");
        cc->SaveAs("baseHisto.png");
        cc->SaveAs("baseHisto.C");
    }
    
    return n_linksTo3Dboard;

}


void plotAll(){

    bool mask[52] = {
        // EE
        0,1,0,1, //  1  2  3  4
        0,1,0,1, //  5  6  7  8
        0,1,0,1, //  9 10 11 12
        0,1,0,1, // 13 14 15 16
        0,1,0,1, // 17 18 19 20
        0,1,0,1, // 21 22 23 24
        0,1,0,1, // 25 26 27 28
        // FH
        1,1,1,1, // 29 30 31 32
        1,1,1,1, // 33 34 35 36
        1,1,1,1, // 37 38 39 40
        // BH
        1,1,1,   // 41 42 43-44
        1,1,     // 45-46 47-48
        1,1      // 49-50 51-52
    };

    std::vector<int> tMux;
    tMux.push_back(12);
    if(!test) tMux.push_back(18);
    if(!test) tMux.push_back(24);
    int NtMux = tMux.size();

    std::vector<int> cluSize;
    int maxNwords = 0;
    int minNwords = 4;
    if(!test) maxNwords = minNwords+10;
    else maxNwords = minNwords+1;
  
    for(int i=minNwords; i<=maxNwords; i++)
        cluSize.push_back(i*16);
    int NcluSize = cluSize.size();

    int results[NtMux][NcluSize];

    if(!test)
      for(int i_mux=0; i_mux<NtMux; i_mux++)
        for(int i_size=0; i_size<NcluSize; i_size++)
	  results[i_mux][i_size] = plotBW("histo_mipTcuts_SE5TE2dR3.root", mask, 
					  tMux[i_mux], cluSize[i_size] );
    
    std::cout << "\t";
    for(int i_size=0; i_size<NcluSize; i_size++)
        std::cout << cluSize[i_size] << "\t";
    std::cout << std::endl;
    
    for(int i_mux=0; i_mux<NtMux; i_mux++){
        std::cout << tMux[i_mux] << "\t";
        for(int i_size=0; i_size<NcluSize; i_size++){
            std::cout << results[i_mux][i_size] << "\t";
        }
        std::cout << std::endl;
    }

}

void plotNseed(TString fin){
  TFile*f = new TFile(fin,"READ");
  TDirectoryFile* d = (TDirectoryFile*)f->Get("hgcalTriggerNtuplizer");
  TTree*t = (TTree*)d->Get("HGCalTriggerNtuple");

  TH1D* h_InclusiveNseed = new TH1D("h_InclusiveNseed","Number of seeds inside C2d", 15, 0, 15);
  TH2D* h_InclusiveNseedXlayer = new TH2D("h_InclusiveNseedXlayer","Number of seeds inside C2d vs layer", 52, 0, 52, 15, 0, 15);  
  t->Project("h_InclusiveNseed","cls_Nseeds");
  t->Project("h_InclusiveNseedXlayer","clsTc_NseedsXlayer:clsTc_l");
  
  TCanvas*c0 = new TCanvas();
  c0->SetLogy();
  c0->SetGrid();
  h_InclusiveNseed->GetXaxis()->SetTitle("N of seeds");
  h_InclusiveNseed->GetYaxis()->SetTitle("a.u");
  h_InclusiveNseed->Scale(1./h_InclusiveNseed->Integral());
  h_InclusiveNseed->Draw();
  c0->SaveAs("testNseed.pdf");

  TH1D*h_AverageDataFormat = new TH1D("h_AverageDataFormat","Weighted Cluster Size", 200, 0, 200);
  for(int i=h_InclusiveNseed->FindBin(1); i<=h_InclusiveNseed->FindBin(5); i++){
    std::cout << "content " << i << "  " << h_InclusiveNseed->GetBinContent(i)  << "  " <<  64 + 32 * (i - h_InclusiveNseed->FindBin(1))  << std::endl;
    h_AverageDataFormat->Fill(h_AverageDataFormat->FindBin( 64 + 32 * (i - h_InclusiveNseed->FindBin(1)) ),  h_InclusiveNseed->GetBinContent(i) );
    std::cout << "----> bin found "<< h_AverageDataFormat->FindBin( 64 + 24 * (i - h_InclusiveNseed->FindBin(1)) ) << std::endl;

  }
  TCanvas*c1 = new TCanvas();
  c1->SetLogy();
  c1->SetGrid();
  gStyle->SetOptStat(0);
  h_AverageDataFormat->GetXaxis()->SetTitle("2D-cluster size [b]");
  h_AverageDataFormat->GetYaxis()->SetTitle("a.u.");
  h_AverageDataFormat->SetLineColor(kBlack);
  h_AverageDataFormat->SetFillColor(38);
  
  
  TPaveText *pC2dAverageSize = new TPaveText(140, 0.2, 190, 0.6 );
  
  //  pC2dAverageSize->AddText( Form("< N > links out: %d", n_linksTo3Dboard) );
  pC2dAverageSize->AddText( Form("Clu 2D < size >: %d b", (int)std::ceil(h_AverageDataFormat->GetMean() ) ) );
  pC2dAverageSize->SetLineColor(kBlack);
  pC2dAverageSize->SetFillStyle(0);
  
  h_AverageDataFormat->Draw();
  pC2dAverageSize->Draw("same");

  c1->SaveAs("testAverageC2dDataFormat.pdf");
  
  TCanvas*c2 = new TCanvas();
  c2->SetLogz();
  c2->SetGrid();
  h_InclusiveNseedXlayer->GetXaxis()->SetTitle("layer");
  h_InclusiveNseedXlayer->GetYaxis()->SetTitle("N of seeds");  
  h_InclusiveNseedXlayer->Scale(1./h_InclusiveNseedXlayer->Integral());
  h_InclusiveNseedXlayer->Draw("colz");
  gStyle->SetNdivisions(1000);
  c2->SaveAs("testNseed_vs_C2dLayer.pdf");
 
}

int main(int argc, char** argv){
    
    TApplication app("app", &argc, argv);

    //plotAll();
    plotNseed("testNseed.root");
    
    app.Run();

    return 0;
    
}

