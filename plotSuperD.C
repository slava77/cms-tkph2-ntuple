#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

void printEfficiency(const TFile* f, const std::string hN, const std::string hD){

}

void plotEffOverlay(const TString ltol, const TString ltolTxt, bool include_g2 = false, bool zoomY = true, bool pt1GeV = true){
  const float minY = zoomY ? 0.9f : 0.0f;
  const float minX = pt1GeV ? 0.51 : 0.31;
  auto fb = new TFile("outHistogramsSuperD_mm3_D16.0cm24.0cm_us1.root");
  auto fg2 = include_g2 ? new TFile("outHistogramsSuperD_mm3_D2.0cm2.0cm_us1.root") : nullptr;  
  auto fg8 = new TFile("outHistogramsSuperD_mm3_D8.0cm8.0cm_us1.root");

  hd_b = (TH1F*)fb->Get("h_num8MH_"+ltol+"_loProdXY_pt");
  hn_b = (TH1F*)fb->Get("h_numSDL_4of4_"+ltol+"_loProdXY_pt");
  he_b = new TEfficiency(*hn_b, *hd_b);
  he_b->SetTitle("Barrel-barrel "+ltolTxt+", PU=140; p_{T}, GeV; Efficiency");

  TH1F* hd_g2;
  TH1F* hn_g2;
  TEfficiency* he_g2;
  if (include_g2){
    hd_g2 = (TH1F*)fg2->Get("h_num8MH_"+ltol+"_loProdXY_pt");
    hn_g2 = (TH1F*)fg2->Get("h_numSDL_4of4_"+ltol+"_loProdXY_pt");
    he_g2 = new TEfficiency(*hn_g2, *hd_g2);
  }

  hd_g8 = (TH1F*)fg8->Get("h_num8MH_"+ltol+"_loProdXY_pt");
  hn_g8 = (TH1F*)fg8->Get("h_numSDL_4of4_"+ltol+"_loProdXY_pt");
  he_g8 = new TEfficiency(*hn_g8, *hd_g8);

  he_b->SetLineColor(kBlue);
  he_b->SetLineWidth(2);
  he_b->SetMarkerStyle(20);
  
  if (include_g2){
    he_g2->SetLineWidth(2);
    he_g2->SetLineColor(kRed);
    he_g2->SetMarkerStyle(24);
  }

  he_g8->SetLineWidth(2);
  he_g8->SetLineColor(kBlack);
  he_g8->SetMarkerStyle(28);

  TCanvas* cv = new TCanvas(ltol, ltol, 600, 600);
  
  he_b->Draw();
  gPad->PaintModified();
  pg = he_b->GetPaintedGraph();
  pg->SetMinimum(minY);
  pg->SetMaximum(1.02);
  auto ax = pg->GetXaxis();

  ax->SetLimits(minX, ax->GetXmax());
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  if (include_g2) he_g2->Draw("same");
  he_g8->Draw("same");

  std::cout<<"Relative ratios:"<<std::endl;
  for (int i = pt1GeV ? 7 : 4; i < 12; ++i){//only 1 GeV and 0.5 GeV options; binning changes
    std::cout<<"p_T > "<<hn_b->GetXaxis()->GetBinLowEdge(i)
	     <<" d(b) "<<hd_b->Integral(i,100) <<" d(8) "<< hd_g8->Integral(i,100)
	     <<" e(b)/e(g8) "<<(hn_b->Integral(i,100)/hd_b->Integral(i,100))/(hn_g8->Integral(i,100)/hd_g8->Integral(i,100));
    if (include_g2) std::cout<<" d(2) "<< hd_g2->Integral(i,100)
			     <<" e(b)/e(g2) "<<(hn_b->Integral(i,100)/hd_b->Integral(i,100))/(hn_g2->Integral(i,100)/hd_g2->Integral(i,100));
    std::cout<<std::endl;
  }
    
  auto leg = new TLegend(0.6, 0.12, 0.9, 0.37);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(he_b, "Equidistant", "LP");
  leg->AddEntry(he_g8, "Grouped 8 cm", "LP");
  if (include_g2) leg->AddEntry(he_g2, "Grouped 2 cm", "LP");
  leg->Draw();
  gPad->SaveAs("heff_"+ltol+"_SDL_8MH_b_vs_g8"+ (include_g2? "_vs_g2" : "")+(zoomY ? "min0.9" : "")+".png");
}
