#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TProfile.h"



//------------------------------------------------------------------------------------------------------------
static const Float_t Pi = TMath::Pi();
static TRandom ran;
static TString HistName;
static char NoP[50];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
static const Int_t N_v2_vs_pt_BW = 14;
static Double_t mt_m0_cut;
static Double_t pt_cut_BW_global;
static Int_t flag_v2_BW_use[N_v2_vs_pt_BW];
static TGraphAsymmErrors* tgae_v2_stat_BW[N_v2_vs_pt_BW];
//------------------------------------------------------------------------------------------------------------

//Error propagation
double ErrorAdd(double x, double y){
    return sqrt(x*x+y*y);
}
double ErrDiv(double x, double y, double dx, double dy){
    if(x < 0) x = -1*x;
    if(y < 0) y = -1*y;
    return x/y*ErrorAdd(dx/x,dy/y);
}
double ErrTime(double x, double y, double dx, double dy){
    if(x < 0) x = -1*x;
    if(y < 0) y = -1*y;
    return x*y*ErrorAdd(dx/x,dy/y);
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void PlotHistErrorBand(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    // Modified March 14th
    Int_t N_points = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y != 0 && x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }
    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line = new TGraph(N_total_points);
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    Hist_line -> SetFillStyle(FillStyle);
    Hist_line -> SetFillColor(FillColor);
    Int_t N_point = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y     = Histo->GetBinContent(i);
        Double_t y_err = Histo->GetBinError(i);
        if(y != 0)
        {
            Double_t x = Histo->GetBinCenter(i);
            if(x >= x_start && x <= x_stop)
            {
                //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
                Hist_line->SetPoint(N_point,x,y-y_err);
                Hist_line->SetPoint(N_total_points-2-N_point,x,y+y_err);
                if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-y_err);
                N_point++;
            }
        }
    }
    Hist_line -> Draw("f");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t One_over_x_FitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    y = 0.0;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    //if(x != 0.0 && !((par0/x) < 0.0 && par2 < 1.0)) y = pow(par0/x,par2) + par1;
    if(x != 0.0) y = par0*TMath::Power(1.0/x,par2) + par1;
    //if(x != 0.0) y =par0/x + par1 + 0.0001*par2;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t TwoGaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + par3*TMath::Gaus(x,par4,par5,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussPolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, pol0, pol1, pol2, pol3, pol4, pol5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    pol0  = par[3];
    pol1  = par[4];
    pol2  = par[5];
    pol3  = par[6];
    pol4  = par[7];
    pol5  = par[8];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t FlowFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t phi, y, v0, v1, v2, v3, v4;
    v0  = par[0];
    v1  = par[1];
    v2  = par[2];
    v3  = par[3];
    v4  = par[4];
    phi = x_val[0];
    y = v0 * (1.0 + 2.0*v1*TMath::Cos(phi) + 2.0*v2*TMath::Cos(2.0*phi)
              + 2.0*v3*TMath::Cos(3.0*phi) + 2.0*v4*TMath::Cos(4.0*phi));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod_x(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------

Double_t LevyFitFunc_pT(Double_t* x_val, Double_t* par)
{
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2 vs. pT
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT, a, b, c, d, n;
    pT = x_val[0];
    n  = par[0]; // number-of-constituent quarks
    a  = par[1];
    b  = par[2];
    c  = par[3];
    d  = par[4];

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. pT/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT_ncq, a, b, c, d, n;
    pT_ncq = x_val[0];
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT_ncq - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc_Poly(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0, par6, par7;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass
    par6   = par[6];
    par7   = par[7];

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = (a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d) + par6*mT_ncq + par7*mT_ncq*mT_ncq;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_m0, a, b, c, d, n, m0;
    mT_m0  = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_m0 + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncE(Double_t* x_val, Double_t* par)
{
    // Original
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        //for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            //r = r_start + j*delta_r;
            r = 0.01;

            //r = R; //

            rho = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_radial(Double_t* x_val, Double_t* par)
{
    // Blast wave function with radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            r = r_start + j*delta_r;

            //r = R; //

            rho  = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------

Double_t BlastWaveFitFunc_no_mass_array(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Double_t mass = par[5];
    Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}


//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncC(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 100;
    Int_t nbins_r      = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}
    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncD(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t NCQ[14]  = {2.0,3.0,3.0,3.0,3.0,3.0,3.0,2.0,3.0,3.0,2.0,2.0,2.0,2.0};
    pt /= NCQ[PID];
    Mass[PID] /= NCQ[PID];
    //pt *= NCQ[PID];
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 25;
    Int_t nbins_r      = 25;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}

    v2 *= NCQ[PID];
    //v2 /= NCQ[PID];
    //v2 *= (1.0/(1.0 + TMath::Exp((pt-4.0/NCQ[PID])/0.6)));

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void BlastWaveSimultaneous(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
    Int_t npfits      = 0;
    Double_t chi2     = 0.0;
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};

    for(Int_t i = 0; i < 14; i++) // loop over PIDs
    {
        p[5] = i; // PID
        //Double_t pt_cut = TMath::Sqrt((Mass[i]+mt_m0_cut)*(Mass[i]+mt_m0_cut) - Mass[i]*Mass[i]); // mt-m0 cut
        Double_t pt_cut = pt_cut_BW_global;
        if(flag_v2_BW_use[i] == 1)
        {
            for(Int_t ix = 0; ix < tgae_v2_stat_BW[i]->GetN(); ix++)
            {
                Double_t x[] = {tgae_v2_stat_BW[i]->GetX()[ix]};
                if(x[0] < pt_cut)
                {
                    Double_t y      = tgae_v2_stat_BW[i]->GetY()[ix];
                    Double_t ye     = tgae_v2_stat_BW[i]->GetErrorYhigh(ix);
//                    ye = 0.01;
                    //cout << "ix = " << ix << ", x = " << x[0] << ", y = " << y << ", ye = " << ye << endl;
                    Double_t bw_val = BlastWaveFitFunc(x,p);
                    //Double_t bw_val = 0.1;
//                    ye += 0.003;
                    Double_t diff   = (y - bw_val)/ye;
                    chi2 += diff*diff;
                    npfits++;
                }
                else break;
            }
        }
    }

    fval = chi2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Hist_interpolate_and_error(TH1F* hist, Double_t x, Double_t &Int_val, Double_t &Int_err)
{
    // Linear interpolation of a histogram
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[2]   = {0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetEntries() == 0) // no data poins to interpolate
    {
        return_val = -1;
        Int_val = 0;
        Int_err = 0;
        //cout << "No entries in histogram" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t binx = 1; binx < hist->GetNbinsX(); binx++)
        {
            Double_t bin_error_val   = hist->GetBinError(binx);
            Double_t bin_x_pos       = hist->GetBinCenter(binx);
            if(bin_error_val != 0)
            {
                err_counter++;
                bin_entries[1] = hist->GetBinContent(binx);
                bin_error[1]   = hist->GetBinError(binx);
                bin_x_val[1]   = hist->GetBinCenter(binx);
                if(bin_x_pos >= x)
                {
                    binx_high = binx;
                    flag_max  = 1;
                    break;
                }
                else flag_max = 0;
            }
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val = bin_entries[1];
            Int_err = bin_error[1];
            return return_val;
        }
        for(Int_t binx_low = binx_high; binx_low > 0; binx_low--)
        {
            bin_entries[0] = hist->GetBinContent(binx_low);
            bin_error[0]   = hist->GetBinError(binx_low);
            bin_x_val[0]   = hist->GetBinCenter(binx_low);
            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        if(bin_error[0] != 0 && bin_error[1] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;
                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[1];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;
                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val = bin_entries[0];
                Int_err = bin_error[0];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t TGraphAsymmErrors_interpolate_and_error(TGraphAsymmErrors* hist, Double_t x, Double_t &Int_val, Double_t &Int_err_low,Double_t &Int_err_high)
{
    // V2: 03.12.2012 -> bug fixed to calculate error bars
    // Linear interpolation of a TGraphAsymmErrors
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[4]   = {0,0,0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetN() == 0) // no data poins to interpolate
    {
        return_val   = -1;
        Int_val      = 0;
        Int_err_low  = 0;
        Int_err_high = 0;
        //cout << "No entries in TGraphAsymmErrors" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t epoint = 0; epoint < hist->GetN(); epoint++)
        {
            hist->GetPoint(epoint,bin_x_val[1],bin_entries[1]);
            bin_error[2] = hist->GetErrorYlow(epoint);
            bin_error[3] = hist->GetErrorYhigh(epoint);

            err_counter++;
            if(bin_x_val[1] >= x)
            {
                binx_high = epoint;
                flag_max  = 1;
                break;
            }
            else flag_max = 0;
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val      = bin_entries[1];
            Int_err_low  = bin_error[2];
            Int_err_high = bin_error[3];
            return return_val;
        }
        for(Int_t epoint = binx_high; epoint >= 0; epoint--)
        {
            hist->GetPoint(epoint,bin_x_val[0],bin_entries[0]);
            bin_error[0] = hist->GetErrorYlow(epoint);
            bin_error[1] = hist->GetErrorYhigh(epoint);

            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        //cout << "bin0 = " << bin_error[0] << ", bin2 = " << bin_error[2] << endl;

        if(bin_error[0] != 0 && bin_error[2] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;

                //cout << "x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 << ", y2 = " << y2 << endl;

                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[2];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_low = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);

                y1_err = bin_error[1];
                y2_err = bin_error[3];

                termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                termC  = -(x-x2)/(x2-x1);
                termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_high = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val      = bin_entries[0];
                Int_err_low  = bin_error[0];
                Int_err_high = bin_error[1];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Get_muB(Double_t sqrt_sNN, Int_t flag_input, Double_t &muB_err)
{
    // Returns baryon chemical potential as a function of the center-of-mass energy sqrt(sNN)
    // from http://arxiv.org/pdf/1111.2406v1.pdf
    // flag_input: 0 = parametrization, 1 = STAR fits, 2 = corrected parametrization for 0-80% (see below)

    // corrected parametrization: Centrality data taken from: http://arxiv.org/pdf/0808.2041.pdf
    Double_t centrality_vals[9] = {75,65,55,45,35,25,15,7.5,2.5};
    Double_t muB_Energy[2][2][9] = // [62.4,200][values,error][70-80%,...,0-5%]
    {
        {   // Au+Au @ 62.4 GeV
            {37.7,42.5,47.0,51.3,54.2,54.5,59.4,61.0,62.7}, // values
            {6.5,5.8,5.1,5.2,5.2,5.2,5.4,5.7,6.0}  // errors
        },
        {   // Au+Au @ 200 GeV
            {14.1,15.3,17.7,18.9,18.6,21.3,21.0,22.8,21.9}, // values
            {4.2,4.2,4.2,4.2,4.2,4.2,4.2,4.5,4.5}  // errors
        }
    };

    Double_t muB     = 0.0;
    muB_err          = 0.0;
    Double_t a       = 1.482;  // +/- 0.0037 GeV
    Double_t b       = 0.3517; // +/- 0.009 GeV^-1

    if(!(flag_input == 0 || flag_input == 1 || flag_input == 2))
    {
        cout << "WARNING: Get_muB function, flag_input is wrong. Either 0 or 1. Parametrization is used." << endl;
        flag_input = 0;
    }

    if(sqrt_sNN >= 0.0 && (flag_input == 0 || flag_input == 2))
    {
        muB = a/(1.0+b*sqrt_sNN);
    }
    if(flag_input == 2) // correct the muB values for 0-80%
    {
        Double_t scale_fac[2][2];  // [62.4,200][val,error]
        Double_t muB_mean[2][2]; // [62.4,200][val,error]
        for(Int_t ebeam = 0; ebeam < 2; ebeam++)
        {
            for(Int_t val_err = 0; val_err < 2; val_err++) // calculate mean muB + error
            {
                muB_mean[ebeam][val_err] = 0.0; // approximation for 0-80%
                Double_t total_weight = 0.0;
                for(Int_t cent = 0; cent < 9; cent++)
                {
                    Double_t weight = 1.0;
                    if(cent >= 7) weight = 0.5;
                    if(val_err == 0) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent];
                    if(val_err == 1) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent]*weight*muB_Energy[ebeam][val_err][cent];
                    total_weight += weight;
                }
                if(total_weight > 0.0)
                {
                    if(val_err == 0) muB_mean[ebeam][val_err] /= total_weight;
                    if(val_err == 1)
                    {
                        muB_mean[ebeam][val_err] = TMath::Sqrt(muB_mean[ebeam][val_err]);
                        muB_mean[ebeam][val_err] /= total_weight;
                    }
                }
            }
            scale_fac[ebeam][0] = muB_mean[ebeam][0]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0); // <muB>/0-10% = 0-80%/0-10%
            Double_t term1 = muB_mean[ebeam][1]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0);
            Double_t term2 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][8]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            Double_t term3 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][7]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            scale_fac[ebeam][1] = TMath::Sqrt(term1*term1+term2*term2+term3*term3); // error on scaling factor
        }
        Double_t mean_scale_fac[2] = {(scale_fac[0][0]+scale_fac[1][0])/2.0,TMath::Sqrt(scale_fac[0][1]*scale_fac[0][1]+scale_fac[1][1]*scale_fac[1][1])/2.0}; // [value,error]
        Double_t original_muB = muB;
        muB     = muB * mean_scale_fac[0]; // 0-80% estimation
        muB_err = muB * mean_scale_fac[1]; // 0-80% error estimation
        //cout << "muB = " << original_muB << ", muB(0-80%) = " << muB << ", muB_err(0-80%) = " << muB_err << ", <scale> = " << mean_scale_fac[0] << ", <scale_err> = " << mean_scale_fac[1] << endl;
        //cout << "scale_fac(62 GeV) = " << scale_fac[0][0] << ", scale_fac(200 GeV) = " << scale_fac[1][0] << endl;
        //cout << "scale_fac_err(62 GeV) = " << scale_fac[0][1] << ", scale_fac_err(200 GeV) = " << scale_fac[1][1] << endl;
    }

    if(flag_input == 1)
    {
        // Preliminary muB parameters from Lokesh with strange particles for 0-80% Au+Au collisions
        // 7.7 :  3.61377e-01   1.51360e-02
        // 11.5:  2.59978e-01   1.06312e-02
        // 39  :  8.81846e-02   4.73062e-03
        // 200 :  22 -> Thats a guess

        // Fittet in TmuB_fit.cc macro with One_over_x_FitFunc function
        Double_t par0 = 2.13495e+00;
        Double_t par1 = 4.35916e-04;
        Double_t par2 = 8.67894e-01;
        muB = par0*TMath::Power(1.0/sqrt_sNN,par2) + par1;
    }

    return muB;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
void Write_HTML_header(FILE* HTML_file)
{
    fprintf(HTML_file,"%s \n","<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">");
    fprintf(HTML_file,"%s \n","<html>");
    fprintf(HTML_file,"%s \n","<head>");


    fprintf(HTML_file,"%s \n","  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">");
    fprintf(HTML_file,"%s \n","  <title>STAR Physics Database</title>");
    fprintf(HTML_file,"%s \n","</head>");
    fprintf(HTML_file,"%s \n"," <body bgcolor=\"white\">");

    fprintf(HTML_file,"%s \n","<h2><span class=\"normaltext\"><font color=\"#ff8c00\" face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\">");
    fprintf(HTML_file,"%s \n","Observation of an energy-dependent difference in elliptic flow between particles and anti-particles in relativistic heavy ion collisions");
    fprintf(HTML_file,"%s \n","</font>");

    fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\"><br> <br> </font></span></h2>");

    //fprintf(HTML_file,"%s \n","<hr><h3><font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">Figure 1");
    //fprintf(HTML_file,"%s \n","</font></h3>");
    //fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">");
    //fprintf(HTML_file,"%s \n","<A HREF=\"fig1.png\"> png </A> | <A HREF=\"fig1.eps\"> eps </A> </br>");
    fprintf(HTML_file,"%s \n","</font>");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_title(TString Label, FILE* HTML_file)
{
    TString label_html = "<hr><h3><font face=\\\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\\\" color=\\\"#0000cd\\\">";
    label_html += Label.Data();
    label_html += " </font> </h3> </hr>";
    //label_html += "\"";
    fprintf(HTML_file,"%s \n",label_html.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_tables(TString labelX, TString labelY, TGraphAsymmErrors* data_stat, TGraphAsymmErrors* syst_errors, TGraphAsymmErrors* glob_syst_errors, TString Label, FILE* HTML_file)
{
    fprintf(HTML_file,"%s %s %s \n","<tr><td>",Label.Data(),"</td> </tr>");
    TString label_html = "<table border=1> <tr><th> ";
    label_html += labelX.Data();
    label_html += " </th> <th> ";
    label_html += labelY.Data();
    label_html += " </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>";

    //fprintf(HTML_file,"%s \n","<table border=1> <tr><th> p<sub>T</sub> (GeV/c) </th> <th> v<sub>2</sub> </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>");

    fprintf(HTML_file,"%s \n",label_html.Data());
    Int_t flag_glob_syst = 1;
    if(glob_syst_errors == syst_errors) flag_glob_syst = -1;

    for(Int_t epoint = 0; epoint < data_stat->GetN(); epoint++)
    {
        Double_t x_val, y_val, x_syst, stat_error, low_syst, high_syst, glob_syst;
        data_stat                     ->GetPoint(epoint,x_val,y_val);
        stat_error = data_stat        ->GetErrorYhigh(epoint);
        high_syst  = syst_errors      ->GetErrorYhigh(epoint);
        low_syst   = syst_errors      ->GetErrorYlow(epoint);
        glob_syst  = glob_syst_errors ->GetErrorYhigh(0);
        if(flag_glob_syst == -1) glob_syst = 0.0;

        fprintf(HTML_file,"%s %f %s %f %s %f %s %f %s %f %s %f %s \n","<tr> <td> ",x_val," </td> <td> ",y_val," </td> <td> ",stat_error," </td> <td> ",
                low_syst," </td> <td> ",high_syst," </td> <td> ",glob_syst,"</td></tr>");

    }
    fprintf(HTML_file,"%s \n","</table><hr><table border=1>");
}
//----------------------------------------------------------------------------------------



