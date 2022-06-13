/// C++ Headers
#include <iostream>
#include <fstream>
#include <sstream>

/// ROOT Headers
#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
/// Drawer Utility
#include "./Utility/constant.h"
#include "./Utility/functions.h"
#include "./Utility/draw.C"
#include "./Utility/style.C"

void drawflow(const Char_t *inFile = "../postflowep/flow/postflowep.txt")
{ // main 
  static const Int_t Centrality_total = 4;
  static const Int_t fitting_methods = 2;
  static const Int_t pt_total_New_phi = 8;
  int datatotal = 3 * Centrality_total * fitting_methods * pt_total_New_phi;
  cout << "datatotal: " << datatotal << endl;
  Double_t data[datatotal] ;

  Double_t mpT[Centrality_total][fitting_methods][pt_total_New_phi] ; 	
  Double_t mV2[Centrality_total][fitting_methods][pt_total_New_phi] ; 	
  Double_t mV2StatErr[Centrality_total][fitting_methods][pt_total_New_phi] ; 	

  std::ifstream input(inFile);
  std::string line;
  int iline = 0 ;
  /*for (int i = 0; i < datatotal; i+=1)
  { 
      input.ignore("#")
      input >> data[i];
  }*/

  while (!input.eof()) 
  { // read line 
    getline(input,line);
  
    if (line.length() != 0 && line[0] != '#')
	{
  	  //cout << line << "\n";
	  //cout << input << "\n";
  	  //cout << line[0] << "\n";
          std::istringstream iss(line);
          float num; // The number in the line

          //while the iss is a number 
	  int icolumn = 0;
          while ((iss >> num))
          {
              //look at the number
  	      data[3*iline + icolumn]	= num; 
  	      cout << "line: " << iline << ", column: "<< icolumn<< ", order: "<< 3*iline + icolumn << ", content: "<<  data[3*iline + icolumn] << "\n";
	
	      icolumn ++ ;
          }


	  iline ++ ;
	}
  } // read line

  for (int i = 0 ; i <Centrality_total ; i++ )
  {
    for (int j = 0 ; j <fitting_methods ; j++ )
    {
      for (int k = 0 ; k <pt_total_New_phi ; k++ )
      {
	mpT[i][j][k] = data[ 3 * (16 * i + 8 * j + k ) + 0] ;
	mV2[i][j][k] = data[ 3 * (16 * i + 8 * j + k ) + 1] ;
	mV2StatErr[i][j][k] = data[ 3 * (16 * i + 8 * j + k ) + 2] ;
      }
    }
  }
   // 19.6 GeV cent 0 counting 
  double x_phi_cent0_count[pt_total_New_phi]    = {0};
  double xErr_phi_cent0_count[pt_total_New_phi] = {0};
  double y_phi_cent0_count[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent0_count[pt_total_New_phi] = {0};
  // BES-I 19.6 GeV cent0
  Double_t pt_bin_center_cent0[5] = {0.6515,0.9975,1.3735,1.7775,2.2605};
  Double_t v2_values_cent0[5] = {0.00429508,0.0416252,0.0740384,0.0834012,0.0991101};
  Double_t v2_stat_error_cent0[5] = {0.00881035,0.00663025,0.00863229,0.0135508,0.0197109};
  Double_t v2_syst_low_error_cent0[5] = {0.0128433,0.00288614,0.0043494,0.00550575,0.00324487};
  Double_t v2_syst_high_error_cent0[5] = {0.00919987,0.00337382,0.00549493,0.00441887,0.00261673};
   // 19.6 GeV cent 0 BW 
  double x_phi_cent0_BW[pt_total_New_phi]    = {0};
  double xErr_phi_cent0_BW[pt_total_New_phi] = {0};
  double y_phi_cent0_BW[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent0_BW[pt_total_New_phi] = {0};
   // 19.6 GeV cent 1 counting 
  double x_phi_cent1_count[pt_total_New_phi]    = {0};
  double xErr_phi_cent1_count[pt_total_New_phi] = {0};
  double y_phi_cent1_count[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent1_count[pt_total_New_phi] = {0};
  // BES-I 19.6 GeV cent1
  Double_t pt_bin_center_cent1[3] = {0.7725,1.2955,2.0895};
  Double_t v2_values_cent1[3] = {-0.0131604,0.0614218,0.0682975};
  Double_t v2_stat_error_cent1[3] = {0.0159064,0.0139996,0.033429};
  Double_t v2_syst_low_error_cent1[3] = {0.0121539,0.00563686,0.0154182};
  Double_t v2_syst_high_error_cent1[3] = {0.00814727,0.00523948,0.0122273};
   // 19.6 GeV cent 1 BW 
  double x_phi_cent1_BW[pt_total_New_phi]    = {0};
  double xErr_phi_cent1_BW[pt_total_New_phi] = {0};
  double y_phi_cent1_BW[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent1_BW[pt_total_New_phi] = {0};
   // 19.6 GeV cent 2 counting 
  double x_phi_cent2_count[pt_total_New_phi]    = {0};
  double xErr_phi_cent2_count[pt_total_New_phi] = {0};
  double y_phi_cent2_count[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent2_count[pt_total_New_phi] = {0};
  // BES-I 19.6 GeV cent2
  Double_t pt_bin_center_cent2[5] = {0.6515,0.9985,1.3735,1.7775,2.2645};
  Double_t v2_values_cent2[5] = {0.014229,0.0481282,0.0734469,0.0872473,0.113327};
  Double_t v2_stat_error_cent2[5] = {0.00862292,0.00630574,0.00811914,0.0124024,0.017897};
  Double_t v2_syst_low_error_cent2[5] = {0.00529471,0.00438684,0.00180629,0.00541846,0.00201318};
  Double_t v2_syst_high_error_cent2[5] = {0.00512403,0.00249747,0.00339875,0.00440116,0.00222428};
   // 19.6 GeV cent 2 BW 
  double x_phi_cent2_BW[pt_total_New_phi]    = {0};
  double xErr_phi_cent2_BW[pt_total_New_phi] = {0};
  double y_phi_cent2_BW[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent2_BW[pt_total_New_phi] = {0};
   // 19.6 GeV cent 3 counting 
  double x_phi_cent3_count[pt_total_New_phi]    = {0};
  double xErr_phi_cent3_count[pt_total_New_phi] = {0};
  double y_phi_cent3_count[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent3_count[pt_total_New_phi] = {0};
  // BES-I 19.6 GeV cent3
  Double_t pt_bin_center_cent3[3] = {0.7565,1.2355,1.8915};
  Double_t v2_values_cent3[3] = {0.0410852,0.0768995,0.0821739};
  Double_t v2_stat_error_cent3[3] = {0.0129654,0.0149608,0.02925};
  Double_t v2_syst_low_error_cent3[3] = {0.00499771,0.00244199,0.0016754};
  Double_t v2_syst_high_error_cent3[3] = {0.0070186,0.00311403,0.00217307};

   // 19.6 GeV cent 3 BW 
  double x_phi_cent3_BW[pt_total_New_phi]    = {0};
  double xErr_phi_cent3_BW[pt_total_New_phi] = {0};
  double y_phi_cent3_BW[pt_total_New_phi]      = {0};
  double yErr_stat_phi_cent3_BW[pt_total_New_phi] = {0};

  for (int k = 0 ; k <pt_total_New_phi ; k++ )
  {
     // 19.6 GeV cent 0 counting 
    x_phi_cent0_count[k]    = mpT[0][0][k];
    y_phi_cent0_count[k]      = mV2[0][0][k];
    yErr_stat_phi_cent0_count[k] = mV2StatErr[0][0][k];
     // 19.6 GeV cent 0 BW 
    x_phi_cent0_BW[k]    = mpT[0][1][k]+0.02;
    y_phi_cent0_BW[k]      = mV2[0][1][k];
    yErr_stat_phi_cent0_BW[k] = mV2StatErr[0][1][k];
     // 19.6 GeV cent 1 counting 
    x_phi_cent1_count[k]    = mpT[1][0][k];
    y_phi_cent1_count[k]      = mV2[1][0][k];
    yErr_stat_phi_cent1_count[k] = mV2StatErr[1][0][k];
     // 19.6 GeV cent 1 BW 
    x_phi_cent1_BW[k]    = mpT[1][1][k]+0.02;
    y_phi_cent1_BW[k]      = mV2[1][1][k];
    yErr_stat_phi_cent1_BW[k] = mV2StatErr[1][1][k];
     // 19.6 GeV cent 2 counting 
    x_phi_cent2_count[k]    = mpT[2][0][k];
    y_phi_cent2_count[k]      = mV2[2][0][k];
    yErr_stat_phi_cent2_count[k] = mV2StatErr[2][0][k];
     // 19.6 GeV cent 2 BW 
    x_phi_cent2_BW[k]    = mpT[2][1][k]+0.02;
    y_phi_cent2_BW[k]      = mV2[2][1][k];
    yErr_stat_phi_cent2_BW[k] = mV2StatErr[2][1][k];
     // 19.6 GeV cent 3 counting 
    x_phi_cent3_count[k]    = mpT[3][0][k];
    y_phi_cent3_count[k]      = mV2[3][0][k];
    yErr_stat_phi_cent3_count[k] = mV2StatErr[3][0][k];
     // 19.6 GeV cent 3 BW 
    x_phi_cent3_BW[k]    = mpT[3][1][k]+0.02;
    y_phi_cent3_BW[k]      = mV2[3][1][k];
    yErr_stat_phi_cent3_BW[k] = mV2StatErr[3][1][k];
  }
  
 
    double markersize = 4;

    TGraphErrors *gr_phi_19gev[Centrality_total][fitting_methods];
    TGraphErrors *gr_phi_19gev_BES_I[Centrality_total];
    gr_phi_19gev[0][0] = new TGraphErrors(pt_total_New_phi, x_phi_cent0_count, y_phi_cent0_count, 0 , yErr_stat_phi_cent0_count);
    gr_phi_19gev[0][1] = new TGraphErrors(pt_total_New_phi, x_phi_cent0_BW, y_phi_cent0_BW, 0 , yErr_stat_phi_cent0_BW);
    gr_phi_19gev[1][0] = new TGraphErrors(pt_total_New_phi, x_phi_cent1_count, y_phi_cent1_count, 0 , yErr_stat_phi_cent1_count);
    gr_phi_19gev[1][1] = new TGraphErrors(pt_total_New_phi, x_phi_cent1_BW, y_phi_cent1_BW, 0 , yErr_stat_phi_cent1_BW);
    gr_phi_19gev[2][0] = new TGraphErrors(pt_total_New_phi, x_phi_cent2_count, y_phi_cent2_count, 0 , yErr_stat_phi_cent2_count);
    gr_phi_19gev[2][1] = new TGraphErrors(pt_total_New_phi, x_phi_cent2_BW, y_phi_cent2_BW, 0 , yErr_stat_phi_cent2_BW);
    gr_phi_19gev[3][0] = new TGraphErrors(pt_total_New_phi, x_phi_cent3_count, y_phi_cent3_count, 0 , yErr_stat_phi_cent3_count);
    gr_phi_19gev[3][1] = new TGraphErrors(pt_total_New_phi, x_phi_cent3_BW, y_phi_cent3_BW, 0 , yErr_stat_phi_cent3_BW);
    gr_phi_19gev_BES_I[0] = new TGraphErrors(5, pt_bin_center_cent0, v2_values_cent0, 0 , v2_stat_error_cent0);
    gr_phi_19gev_BES_I[1] = new TGraphErrors(3, pt_bin_center_cent1, v2_values_cent1, 0 , v2_stat_error_cent1);
    gr_phi_19gev_BES_I[2] = new TGraphErrors(5, pt_bin_center_cent2, v2_values_cent2, 0 , v2_stat_error_cent2);
    gr_phi_19gev_BES_I[3] = new TGraphErrors(3, pt_bin_center_cent3, v2_values_cent3, 0 , v2_stat_error_cent3);

    double h1_x1 = -0.1; double h1_x2 = 3.1; double h1_y1 = -0.037; double h1_y2 = 0.198;

    TCanvas *c1=new TCanvas("c1","",6000, 1600);
    c1->Draw();
    TPad *p_1=new TPad("p_1","",0.04,0.1,0.27,1.0);
    //p_1->SetGrid(1,1);
    p_1->SetTicks(1,1);
    p_1->SetLeftMargin(0.14);
    p_1->SetRightMargin(0.);
    p_1->SetBottomMargin(0.1);
    p_1->Draw();
    p_1->cd();

    TH2D *h_v2_1 = new TH2D("h_v2_1","", 9, h1_x1, h1_x2, 9, h1_y1, h1_y2);
    h_v2_1->SetStats(0);
    h_v2_1->SetNdivisions(306,"X");
    h_v2_1->SetNdivisions(306,"Y");
    h_v2_1->GetXaxis()->SetLabelSize(0.045);
    h_v2_1->GetYaxis()->SetLabelSize(0.045);
    h_v2_1->GetXaxis()->SetLabelFont(42);
    h_v2_1->GetYaxis()->SetLabelFont(42);
    h_v2_1->GetXaxis()->SetLabelSize(0.055);
    h_v2_1->GetYaxis()->SetLabelSize(0.055);
    h_v2_1->GetXaxis()->SetLabelOffset(0.01);
    h_v2_1->GetYaxis()->SetLabelOffset(0.01);
    h_v2_1->Draw();

    drawHistBox(h1_x1, h1_x2, h1_y1, h1_y2, 2222);
    drawLine(h1_x1, h1_y1, h1_x2, h1_y1, 2, 1, 1);
    drawLine(h1_x1, h1_y2, h1_x2, h1_y2, 2, 1, 1);
    
    gr_phi_19gev[0][0]->SetMarkerStyle(20);
    gr_phi_19gev[0][0]->SetMarkerSize(markersize);
    gr_phi_19gev[0][0]->SetLineColor(1);
    gr_phi_19gev[0][0]->SetMarkerColor(1);
    gr_phi_19gev[0][0]->Draw("PE");
    gr_phi_19gev[0][1]->SetMarkerStyle(21);
    gr_phi_19gev[0][1]->SetMarkerSize(markersize);
    gr_phi_19gev[0][1]->SetLineColor(2);
    gr_phi_19gev[0][1]->SetMarkerColor(2);
    gr_phi_19gev[0][1]->Draw("PE");
    gr_phi_19gev_BES_I[0]->SetMarkerStyle(24);
    gr_phi_19gev_BES_I[0]->SetMarkerSize(markersize);
    gr_phi_19gev_BES_I[0]->SetLineColor(4);
    gr_phi_19gev_BES_I[0]->SetMarkerColor(4);
    gr_phi_19gev_BES_I[0]->Draw("PE");
    drawText(0.02,0.085, "19.6 GeV (0-80%)", 42, 0.05, 0);
    
    c1->cd();
    TPad *p_2=new TPad("p_2","",0.27,0.1,0.5,1.0);
    //p_2->SetGrid(1,1);
    p_2->SetTicks(1,1);
    p_2->SetLeftMargin(0.);
    p_2->SetRightMargin(0.0);
    p_2->SetBottomMargin(0.1);
    p_2->Draw();
    p_2->cd();

    TH2D *h_v2_2 = new TH2D("h_v2_2","", 9, h1_x1, h1_x2, 9, h1_y1, h1_y2);
    h_v2_2->SetStats(0);
    h_v2_2->SetNdivisions(306,"X");
    h_v2_2->SetNdivisions(306,"Y");
    h_v2_2->GetXaxis()->SetLabelSize(0.045);
    h_v2_2->GetYaxis()->SetLabelSize(0.045);
    h_v2_2->GetXaxis()->SetLabelFont(42);
    h_v2_2->GetYaxis()->SetLabelFont(42);
    h_v2_2->GetXaxis()->SetLabelSize(0.055);
    h_v2_2->GetYaxis()->SetLabelSize(0.055);
    h_v2_2->GetXaxis()->SetLabelOffset(0.01);
    h_v2_2->GetYaxis()->SetLabelOffset(0.01);
    h_v2_2->Draw();

    drawHistBox(h1_x1, h1_x2, h1_y1, h1_y2, 2222);
    drawLine(h1_x1, h1_y1, h1_x2, h1_y1, 2, 1, 1);
    drawLine(h1_x1, h1_y2, h1_x2, h1_y2, 2, 1, 1);

    gr_phi_19gev[1][0]->SetMarkerStyle(20);
    gr_phi_19gev[1][0]->SetMarkerSize(markersize);
    gr_phi_19gev[1][0]->SetLineColor(1);
    gr_phi_19gev[1][0]->SetMarkerColor(1);
    gr_phi_19gev[1][0]->Draw("PE");
    gr_phi_19gev[1][1]->SetMarkerStyle(21);
    gr_phi_19gev[1][1]->SetMarkerSize(markersize);
    gr_phi_19gev[1][1]->SetLineColor(2);
    gr_phi_19gev[1][1]->SetMarkerColor(2);
    gr_phi_19gev[1][1]->Draw("PE");
    gr_phi_19gev_BES_I[1]->SetMarkerStyle(24);
    gr_phi_19gev_BES_I[1]->SetMarkerSize(markersize);
    gr_phi_19gev_BES_I[1]->SetLineColor(4);
    gr_phi_19gev_BES_I[1]->SetMarkerColor(4);
    gr_phi_19gev_BES_I[1]->Draw("PE");
    drawText(0.02,0.085, "(0-10%)", 42, 0.05, 0);

    c1->cd();
    TPad *p_3=new TPad("p_3","",0.5,0.1,0.73,1.0);
    //p_3->SetGrid(1,1);
    p_3->SetTicks(1,1);
    p_3->SetLeftMargin(0.);
    p_3->SetRightMargin(0.0);
    p_3->SetBottomMargin(0.1);
    p_3->Draw();
    p_3->cd();

    TH2D *h_v2_3 = new TH2D("h_v2_3","", 9, h1_x1, h1_x2, 9, h1_y1, h1_y2);
    h_v2_3->SetStats(0);
    h_v2_3->SetNdivisions(306,"X");
    h_v2_3->SetNdivisions(306,"Y");
    h_v2_3->GetXaxis()->SetLabelSize(0.045);
    h_v2_3->GetYaxis()->SetLabelSize(0.045);
    h_v2_3->GetXaxis()->SetLabelFont(42);
    h_v2_3->GetYaxis()->SetLabelFont(42);
    h_v2_3->GetXaxis()->SetLabelSize(0.055);
    h_v2_3->GetYaxis()->SetLabelSize(0.055);
    h_v2_3->GetXaxis()->SetLabelOffset(0.01);
    h_v2_3->GetYaxis()->SetLabelOffset(0.01);
    h_v2_3->Draw();

    drawHistBox(h1_x1, h1_x2, h1_y1, h1_y2, 2222);
    drawLine(h1_x1, h1_y1, h1_x2, h1_y1, 2, 1, 1);
    drawLine(h1_x1, h1_y2, h1_x2, h1_y2, 2, 1, 1);

    gr_phi_19gev[2][0]->SetMarkerStyle(20);
    gr_phi_19gev[2][0]->SetMarkerSize(markersize);
    gr_phi_19gev[2][0]->SetLineColor(1);
    gr_phi_19gev[2][0]->SetMarkerColor(1);
    gr_phi_19gev[2][0]->Draw("PE");
    gr_phi_19gev[2][1]->SetMarkerStyle(21);
    gr_phi_19gev[2][1]->SetMarkerSize(markersize);
    gr_phi_19gev[2][1]->SetLineColor(2);
    gr_phi_19gev[2][1]->SetMarkerColor(2);
    gr_phi_19gev[2][1]->Draw("PE");
    gr_phi_19gev_BES_I[2]->SetMarkerStyle(24);
    gr_phi_19gev_BES_I[2]->SetMarkerSize(markersize);
    gr_phi_19gev_BES_I[2]->SetLineColor(4);
    gr_phi_19gev_BES_I[2]->SetMarkerColor(4);
    gr_phi_19gev_BES_I[2]->Draw("PE");
    drawText(0.02,0.085, "(10-40%)", 42, 0.05, 0);

    c1->cd();
    TPad *p_4=new TPad("p_4","",0.73,0.1,0.96,1.0);
    //p_4->SetGrid(1,1);
    p_4->SetTicks(1,1);
    p_4->SetLeftMargin(0.);
    p_4->SetRightMargin(0.1);
    p_4->SetBottomMargin(0.1);
    p_4->Draw();
    p_4->cd();

    TH2D *h_v2_4 = new TH2D("h_v2_4","", 9, h1_x1, h1_x2, 9, h1_y1, h1_y2);
    h_v2_4->SetStats(0);
    h_v2_4->SetNdivisions(306,"X");
    h_v2_4->SetNdivisions(306,"Y");
    h_v2_4->GetXaxis()->SetLabelSize(0.045);
    h_v2_4->GetYaxis()->SetLabelSize(0.045);
    h_v2_4->GetXaxis()->SetLabelFont(42);
    h_v2_4->GetYaxis()->SetLabelFont(42);
    h_v2_4->GetXaxis()->SetLabelSize(0.055);
    h_v2_4->GetYaxis()->SetLabelSize(0.055);
    h_v2_4->GetXaxis()->SetLabelOffset(0.01);
    h_v2_4->GetYaxis()->SetLabelOffset(0.01);
    h_v2_4->Draw();

    drawHistBox(h1_x1, h1_x2, h1_y1, h1_y2, 2222);
    drawLine(h1_x1, h1_y1, h1_x2, h1_y1, 2, 1, 1);
    drawLine(h1_x1, h1_y2, h1_x2, h1_y2, 2, 1, 1);

    gr_phi_19gev[3][0]->SetMarkerStyle(20);
    gr_phi_19gev[3][0]->SetMarkerSize(markersize);
    gr_phi_19gev[3][0]->SetLineColor(1);
    gr_phi_19gev[3][0]->SetMarkerColor(1);
    gr_phi_19gev[3][0]->Draw("PE");
    gr_phi_19gev[3][1]->SetMarkerStyle(21);
    gr_phi_19gev[3][1]->SetMarkerSize(markersize);
    gr_phi_19gev[3][1]->SetLineColor(2);
    gr_phi_19gev[3][1]->SetMarkerColor(2);
    gr_phi_19gev[3][1]->Draw("PE");
    gr_phi_19gev_BES_I[3]->SetMarkerStyle(24);
    gr_phi_19gev_BES_I[3]->SetMarkerSize(markersize);
    gr_phi_19gev_BES_I[3]->SetLineColor(4);
    gr_phi_19gev_BES_I[3]->SetMarkerColor(4);
    gr_phi_19gev_BES_I[3]->Draw("PE");
    drawText(0.02,0.085, "(40-80%)", 42, 0.05, 0);

    c1->cd();
    drawText(0.47,0.07, "pT (GeV/c)", 42, 0.065, 0);

    c1->cd();
    drawText(0.03, 0.47, "v_{2}", 42, 0.065, 90);
    c1->cd(); p_3->cd();
    //drawText(0.35,0.020, "(a) Positive Particle", 42, 0.06, 0);
    double x1 = 0.130;
    double x2 = 0.230;
    double x3 = 0.330;

    double y1 = 0.175;
    double y2 = 0.150;
    double y3 = 0.125;

    Draw_Point(x1, y1,  0, 0, 20, kBlack, 4);
    Draw_Point(x1, y2,  0, 0, 21, kRed, 4);
    Draw_Point(x1, y3,  0, 0, 24, kBlue, 4);
    
    drawText(0.2, y1, "Bin counting", 62, 0.045, 0);
    drawText(0.2, y2, "Breit Wigner", 62, 0.045, 0);
    drawText(0.2, y3, "BES-I", 62, 0.045, 0);

    c1->SaveAs("fig_phi_v2_19GeV.pdf");
} // main
