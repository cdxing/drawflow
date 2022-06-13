#include "TLine.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphAsymmErrors.h"

void drawHistBox(double x1=0., double x2=1., double y1=0., double y2=1., int lineWidth=3333)
//void drawHistBox(double x1, double x2, double y1, double y2, int lineWidth)
{
  TLine *l1 = new TLine(x1,y1,x2,y1);
  l1->SetLineWidth(lineWidth/1000);
  l1->Draw("same");
  TLine *l2 = new TLine(x1,y2,x2,y2);
  l2->SetLineWidth((lineWidth%1000)/100);
  l2->Draw("same");
  TLine *l3 = new TLine(x1,y1,x1,y2);
  l3->SetLineWidth((lineWidth%100)/10);
  l3->Draw("same");
  TLine *l4 = new TLine(x2,y1,x2,y2);
  l4->SetLineWidth(lineWidth%10);
  l4->Draw("same");
}

void drawBox(double x1=0., double y1=0., double x2=0., double y2=0., double fillcolor=19, double linecolor=0, double linewidth=1)
{
	TBox *box = new TBox(x1, y1, x2, y2);
	box->SetFillColor(fillcolor);
	box->SetLineColor(linecolor);
	box->SetLineWidth(linewidth);
	box->Draw("same");
}


void drawText(double x1=0., double y1=0., const char* text="test", int textFont=42, double textSize=0.05, double textAngle=0)
{
  TLatex *tex = new TLatex(x1, y1, text);
  tex->SetTextFont(textFont);
  tex->SetTextSize(textSize);
  tex->SetTextAngle(textAngle);
  tex->Draw("same");
}

void DrawText(double x1=0., double y1=0., const char* text="test", int textFont=42, double textSize=0.05, double textAngle=0, double textColor=1)
{
  TLatex *tex = new TLatex(x1, y1, text);
  tex->SetTextFont(textFont);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetTextAngle(textAngle);
  tex->SetNDC();
  tex->Draw("same");
}

void drawLine(double x1=0., double y1=0., double x2=0., double y2=1., int lineWidth=1, int lineStyle=1, int lineColor=1)
//void drawLine(double x1, double y1, double x2, double y2, int lineWidth, int lineStyle, int lineColor)
{
  TLine *la = new TLine(x1,y1,x2,y2);
  la->SetLineWidth(lineWidth);
  la->SetLineStyle(lineStyle);
  la->SetLineColor(lineColor);
  la->Draw("same");
}
void Draw_TGAE_new_Symbol(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size)
 {
     TString HistName;
     Float_t size_A = 1.35*size;
     Float_t size_B = size;
     Float_t size_C = size;
     Int_t alt_marker = style;
     if(style == 24)
     {
         alt_marker = 20;
         size_A = 1.35*size;
     }
     if(style == 25)
     {
         alt_marker = 21;
         size_A = 1.35*size;
     }
     if(style == 26)
     {
         alt_marker = 22;
         size_A = 1.4*size;
     }
     if(style == 23)
     {
         alt_marker = 23;
         size_A = 1.35*size;
     }
     if(style == 30 || style == 29)
     {
         alt_marker = 29;
         size_A = 1.55*size;
     }

     // black and filled outer marker
     HistName = "tgae_dummy_A";
     TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
     ge_clone_A->SetMarkerSize(size_A);
     ge_clone_A->SetMarkerStyle(alt_marker);
  ge_clone_A->SetMarkerColor(1);
     ge_clone_A->SetLineColor(10);
     ge_clone_A->SetLineWidth(10);
     ge_clone_A->Draw("same PZ0");

     // white and filled inner marker
     HistName = "tgae_dummy_B";
     TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
     ge_clone_B->SetMarkerSize(size_B);
     ge_clone_B->SetMarkerStyle(alt_marker);
     ge_clone_B->SetMarkerColor(10);
     ge_clone_B->SetLineColor(10);
     ge_clone_B->SetLineWidth(10);
     ge_clone_B->Draw("same PZ0");

     // color inner marker
     HistName = "tgae_dummy_C";
     TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
     ge_clone_C->SetMarkerSize(size_C);
     ge_clone_C->SetMarkerStyle(style);
     ge_clone_C->SetMarkerColor(color);
     ge_clone_C->SetLineColor(1);
     ge_clone_C->SetLineWidth(5);
     ge_clone_C->Draw("same PZ0");
 }


 void Draw_TGAE_Point_new_Symbol(Double_t x_val, Double_t y_val, Double_t x_min_err, Double_t x_max_err, Double_t y_min_err, Double_t y_max_err, Int_t style, Int_t color, Float_t size)
 {
     TGraphAsymmErrors* tgae = new TGraphAsymmErrors();
     tgae->SetPoint(0,x_val,y_val);
     tgae->SetPointError(0,x_min_err,x_max_err,y_min_err,y_max_err);
	 tgae->SetLineWidth(10);

     TString HistName;
     Float_t size_A = 1.35*size;
     Float_t size_B = size;
     Float_t size_C = size;
     Int_t alt_marker = style;
     Int_t style_in = style;
	 if(style == 20)
	 {
		 size_A = 1.6*size;
	 }
     if(style == 24)
     {
         alt_marker = 20;
         size_A = 1.35*size;
     }
     if(style == 25)
     {
         alt_marker = 21;
         size_A = 1.35*size;
     }
     if(style == 26)
     {
         alt_marker = 22;
         size_A = 1.4*size;
     }
     if(style == 23)
     {
         alt_marker = 23;
         size_A = 1.35*size;
     }
     if(style == 30 || style == 29)
     {
         alt_marker = 29;
         size_A = 1.55*size;
     }
     if(style == 260)
     {
         alt_marker = 26;
         size_A = 1.15*size;
         style = 26;
     }
     if(style == 300)
 {
         alt_marker = 30;
         size_A = 1.3*size;
         style = 30;
     }

     // black and filled outer marker
     HistName = "tgae_dummy_A";
     TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
     ge_clone_A->SetMarkerSize(size_A);
     ge_clone_A->SetMarkerStyle(alt_marker);
     ge_clone_A->SetMarkerColor(1);
     ge_clone_A->SetLineColor(10);
     if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
     ge_clone_A->Draw("same PZ0");

     if(!(style_in == 260 || style_in == 300))
     {
         // white and filled inner marker
         HistName = "tgae_dummy_B";
         TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
         ge_clone_B->SetMarkerSize(size_B);
         ge_clone_B->SetMarkerStyle(alt_marker);
         ge_clone_B->SetMarkerColor(10);
         ge_clone_B->SetLineColor(10);
         ge_clone_B->Draw("same PZ0");
     }

     // color inner marker
     HistName = "tgae_dummy_C";
     TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
     ge_clone_C->SetMarkerSize(size_C);
     ge_clone_C->SetMarkerStyle(style);
     ge_clone_C->SetMarkerColor(color);
     ge_clone_C->SetLineColor(1);
     ge_clone_C->Draw("same PZ0");
 }

void Draw_Point(Double_t x_val, Double_t y_val, Double_t x_err, Double_t y_err, Int_t style, Int_t color, Float_t size)
{
	TGraphErrors* tgae = new TGraphErrors();
	tgae->SetPoint(0,x_val,y_val);
	tgae->SetPointError(0,x_err,y_err);
	tgae->SetLineWidth(1);
	tgae->SetLineColor(color);
	tgae->SetMarkerStyle(style);
	tgae->SetMarkerSize(size);
	tgae->SetMarkerColor(color);
	tgae->Draw("SAME PE");
}

void Draw_Band(Double_t x1_val, Double_t y1_val, Double_t x2_val, Double_t y2_val, Double_t x_err, Double_t y_err, Int_t fillstyle, Int_t fillcolor)
{
	TGraphErrors *tgae = new TGraphErrors();
	tgae->SetPoint(0,x1_val,y1_val);
	tgae->SetPoint(1,x2_val,y2_val);
	tgae->SetPointError(0,x_err,y_err);
	tgae->SetPointError(1,x_err,y_err);
	tgae->SetFillColor(fillcolor);
	tgae->SetFillStyle(fillstyle);
	tgae->SetLineWidth(5);
	tgae->Draw("e3 SAME");
}

void Draw_Sys(Double_t x_val, Double_t y_val, Double_t x_err, Double_t y_err, Int_t style, Int_t color, Float_t size)
{
    TGraphErrors* tgae = new TGraphErrors();
    tgae->SetPoint(0,x_val,y_val);
    tgae->SetPointError(0,x_err,y_err);
    //tgae->SetLineWidth(1);
    //tgae->SetLineColor(color);
    //tgae->SetMarkerStyle(style);
    //tgae->SetMarkerSize(size);
    //tgae->SetMarkerColor(color);
    tgae->SetFillColor(color);
    tgae->Draw("2");
}


