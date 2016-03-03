/* Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("trkTree/tree")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create tkph2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] tkph2 cms2 

   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
*/
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <array>
#include <utility>
#include <algorithm>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include "TVector3.h"

#include "tkph2.cc"

namespace tkph2consts{
  constexpr int nLayers = 10;
}
using namespace tas;
using namespace tkph2consts;

//histograms for occupancy and non-zero occupancy vs value-X
struct Hists2DProfWNZ {
  TH2D* h2D;
  TProfile* profile;
  TProfile* profileNZ;
  Hists2DProfWNZ() : h2D(nullptr), profile(nullptr), profileNZ(nullptr) {}
  Hists2DProfWNZ(const std::string& tBase, const std::vector<double>& bX,
		 int nY, double minY, double maxY, bool makeNZ = true){
    init(tBase, bX, nY, minY, maxY, makeNZ);
  }
  void init(const std::string& tBase, const std::vector<double>& bX,
	    int nY, double minY, double maxY, bool makeNZ = true){
    std::string hN = tBase+"_2D";
    h2D = new TH2D(hN.c_str(), hN.c_str(), bX.size()-1, bX.data(), nY, minY, maxY);
    hN = tBase+"_Profile";
    profile = new TProfile(hN.c_str(), hN.c_str(), bX.size()-1, bX.data());
    hN = tBase+"_ProfileNZ";    
    if (makeNZ) profileNZ = new TProfile(hN.c_str(), hN.c_str(), bX.size()-1, bX.data());
  }
  void write(TFile* f){
    TDirectory* wDir = gDirectory;
    f->cd();
    if (h2D) h2D->Write();
    if (profile) profile->Write();
    if (profileNZ) profileNZ->Write();
    wDir->cd();
  }
  void fill(double x, double occ, bool clampOverflowY2D = true){
    const TAxis* ay = h2D->GetYaxis();
    double max = ay->GetBinUpEdge(ay->GetNbins());
    double eps = 0.001*(max - ay->GetBinLowEdge(ay->GetNbins()));
    h2D->Fill(x, clampOverflowY2D ? std::min(occ, max-eps) : occ);
    if (profile) profile->Fill(x, occ);
    if (profileNZ) profileNZ->Fill(x, occ > 0);
  }
};

//a set of variables and efficiencies as a function of 1D variable
struct HistoSet1D {
  TH1D* den;
  TH1D* denDesMDFid;//fiducial wrt design mini-doublet expectation
  TH1D* denDesMDFid1TP;//fiducial wrt design mini-doublet expectation; one matching TP
  TH1D* num_aCut;
  TEfficiency* eff_aCut;

  //implies there is a truth match in inner and outer layers
  TH1D* denActAnySH;
  TH1D* denAct;
  TEfficiency* effActAnySH;//matching efficiency
  TEfficiency* effAct;//matching efficiency
  //as above, but starting from denDesMDFid1TP
  TH1D* denActAnySH1TP;
  TH1D* denAct1TP;
  TEfficiency* effActAnySH1TP;//matching efficiency
  TEfficiency* effAct1TP;//matching efficiency

  TH1D* num2SH_aCut;
  TEfficiency* eff2SH_aCut;
  TH1D* numAct_aCut;
  TEfficiency* effAct_aCut;
  TH1D* numRec_aCut;
  TEfficiency* effRec_aCut;

  TH2D* dminiDirActVs2mm;
  TH2D* dminiDirActVs2mm_pi;
  
  //plain count of other hits on the same layer
  Hists2DProfWNZ others_inModule;
  Hists2DProfWNZ othersRec_inModule;

  //minidoublets combinations of base hit with others
  //among those passing a cut given a d=2mm or design
  Hists2DProfWNZ mdOthers2mm_aCut;
  Hists2DProfWNZ mdOthersDes_aCut;
  Hists2DProfWNZ mdOthers2mm2SH_aCut;
  Hists2DProfWNZ mdOthersDes2SH_aCut;
  Hists2DProfWNZ mdOthersAct_aCut;
  Hists2DProfWNZ mdOthersRec_aCut;

  //minidoublets combinations of base hit with others closer than correct
  //among those passing a cut given a d=2mm or design
  Hists2DProfWNZ mdOthers2mm_aCutCloser;
  Hists2DProfWNZ mdOthersDes_aCutCloser;
  Hists2DProfWNZ mdOthers2mm2SH_aCutCloser;
  Hists2DProfWNZ mdOthersDes2SH_aCutCloser;
  Hists2DProfWNZ mdOthersAct_aCutCloser;
  Hists2DProfWNZ mdOthersRec_aCutCloser;
  
  void write(TFile* f){
    TDirectory* wDir = gDirectory;
    f->cd();
    if (den) den->Write();
    if (denDesMDFid) denDesMDFid->Write();
    if (denDesMDFid1TP) denDesMDFid1TP->Write();
    if (num_aCut) num_aCut->Write();

    if (denActAnySH) denActAnySH->Write();
    if (denAct) denAct->Write();
    if (denActAnySH1TP) denActAnySH1TP->Write();
    if (denAct1TP) denAct1TP->Write();
    if (num2SH_aCut) numAct_aCut->Write();
    if (numAct_aCut) numAct_aCut->Write();
    if (numRec_aCut) numRec_aCut->Write();
    //    eff_aCut->Write(); //this write doesn't work

    if (dminiDirActVs2mm) dminiDirActVs2mm->Write();
    if (dminiDirActVs2mm_pi) dminiDirActVs2mm_pi->Write();
    
    others_inModule.write(f);
    othersRec_inModule.write(f);

    mdOthers2mm_aCut.write(f);
    mdOthersDes_aCut.write(f);
    mdOthers2mm2SH_aCut.write(f);
    mdOthersDes2SH_aCut.write(f);
    mdOthersAct_aCut.write(f);
    mdOthersRec_aCut.write(f);

    mdOthers2mm_aCutCloser.write(f);
    mdOthersDes_aCutCloser.write(f);
    mdOthers2mm2SH_aCutCloser.write(f);
    mdOthersDes2SH_aCutCloser.write(f);
    mdOthersAct_aCutCloser.write(f);
    mdOthersRec_aCutCloser.write(f);

    wDir->cd();
  }
};


void createEfficiency(const TH1* pass, const TH1* all, TEfficiency* &eff){
  eff = new TEfficiency(*pass, *all);
  eff->SetTitle(Form("Eff for %s; %s; #epsilon", pass->GetName(), pass->GetXaxis()->GetTitle() ));
}

void plotEffOverlay(std::array<TEfficiency*,3> effs, const std::string& title, const std::string& ofname,
		    std::array<std::string, 3> llabel, float min, float max = - 1, bool keepCanvas = false){
  auto efh = effs[0];
  efh->SetTitle(title.c_str());
  auto cn = efh->GetTitle();
  TCanvas* cv = new TCanvas(cn, cn, 600, 600);
  cv->cd();
  efh->Draw();
  efh->SetMarkerStyle(22);
  efh->SetLineWidth(2);
  efh->SetLineColor(1);
  efh->SetMarkerSize(0.8);
  gPad->PaintModified();
  auto pg = efh->GetPaintedGraph();
  pg->SetMinimum(min);
  if (max>0) pg->SetMaximum(max);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  
  efh = effs[1];
  efh->Draw("same");
  efh->SetMarkerStyle(23);
  efh->SetLineWidth(2);
  efh->SetLineColor(2);
  efh->SetMarkerSize(0.8);
  gPad->PaintModified();
  
  efh = effs[2];
  efh->Draw("same");
  efh->SetMarkerStyle(24);
  efh->SetLineWidth(2);
  efh->SetLineColor(3);
  efh->SetMarkerSize(0.8);
  gPad->PaintModified();
  
  TLegend* leg = new TLegend(0.4, 0.3, 0.9, 0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(effs[0], llabel[0].c_str());
  leg->AddEntry(effs[1], llabel[1].c_str());
  leg->AddEntry(effs[2], llabel[2].c_str());
  leg->Draw();
  
  gPad->SaveAs(ofname.c_str());
  if (! keepCanvas) delete cv;
}

void plot2DProf(Hists2DProfWNZ hs, const std::string& aTitle, float miny, int logXYZ = 0, bool keepCanvas = false){
  auto h2 = hs.h2D;
  h2->SetTitle(aTitle.c_str());
  auto cn = h2->GetTitle();
  TCanvas* cv = new TCanvas(cn, cn, 600, 600);
  cv->cd();
  gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
  h2->Draw("colz");
  auto ay = h2->GetYaxis();
  h2->GetYaxis()->SetLimits(miny, ay->GetXmax());
  h2->SetStats(0);
  gPad->SetGridx();
  gPad->SetLogx((logXYZ & 0x1));
  gPad->SetLogy((logXYZ & 0x2));
  gPad->SetLogz((logXYZ & 0x4));
  auto hp = hs.profile;
  if (hp){
    hp->SetLineColor(2);
    hp->SetLineWidth(2);
    hp->Draw("same");
  }
  auto hz = hs.profileNZ;
  if (hz){
    hz->SetLineColor(1);
    hz->SetLineWidth(2);
    hz->Draw("same");
  }
  
  gPad->SaveAs(Form("%s.png", h2->GetName()) );
  if (! keepCanvas) delete cv;
}


void createPtHistograms(std::array<HistoSet1D, nLayers+1>& layerMD, const std::string& ext, std::vector<double>& ptBins){
  for (int iL = 1; iL <= nLayers; ++iL){
    std::ostringstream oS; oS<< "layer"<<iL<<"MD_den_pt_"<<ext;
    auto aS = oS.str().c_str();
    layerMD[iL].den = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_denDesMDFid_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].denDesMDFid = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_denDesMDFid1TP_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].denDesMDFid1TP = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_num_aCut_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].num_aCut = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    
    oS.str(""); oS<< "layer"<<iL<<"MD_denActAnySH_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].denActAnySH = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_denAct_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].denAct = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);

    oS.str(""); oS<< "layer"<<iL<<"MD_denActAnySH1TP_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].denActAnySH1TP = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_denAct1TP_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].denAct1TP = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);

    oS.str(""); oS<< "layer"<<iL<<"MD_num2SH_aCut_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].num2SH_aCut = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_numAct_aCut_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].numAct_aCut = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);
    oS.str(""); oS<< "layer"<<iL<<"MD_numRec_aCut_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].numRec_aCut = new TH1D(aS, aS, ptBins.size()-1, &ptBins[0]);

    oS.str(""); oS<< "layer"<<iL<<"MD_dminiDirActVs2mm_pt_"<<ext; aS = oS.str().c_str();
    layerMD[iL].dminiDirActVs2mm = new TH2D(aS, aS, ptBins.size()-1, &ptBins[0], 200, -0.2, 0.2);
    oS.str(""); oS<< "layer"<<iL<<"MD_dminiDirActVs2mm_pt_pi_"<<ext; aS = oS.str().c_str();
    layerMD[iL].dminiDirActVs2mm_pi = new TH2D(aS, aS, ptBins.size()-1, &ptBins[0], 200, -0.2, 0.2);
    
    layerMD[iL].others_inModule.init(Form("layer%dMD_others_inModule_pt_%s", iL, ext.c_str()), ptBins, 100, 0, 100, false);
    layerMD[iL].othersRec_inModule.init(Form("layer%dMD_othersRec_inModule_pt_%s", iL, ext.c_str()), ptBins, 100, 0, 100, false);

    layerMD[iL].mdOthers2mm_aCut.init(Form("layer%dMD_mdOthers2mm_aCut_pt_%s", iL, ext.c_str()), ptBins, 20, 0, 20);
    layerMD[iL].mdOthersDes_aCut.init(Form("layer%dMD_mdOthersDes_aCut_pt_%s", iL, ext.c_str()), ptBins, 20, 0, 20);
    layerMD[iL].mdOthers2mm2SH_aCut.init(Form("layer%dMD_mdOthers2mm2SH_aCut_pt_%s", iL, ext.c_str()), ptBins, 20, 0, 20);
    layerMD[iL].mdOthersDes2SH_aCut.init(Form("layer%dMD_mdOthersDes2SH_aCut_pt_%s", iL, ext.c_str()), ptBins, 20, 0, 20);
    layerMD[iL].mdOthersAct_aCut.init(Form("layer%dMD_mdOthersAct_aCut_pt_%s", iL, ext.c_str()), ptBins, 20, 0, 20);
    layerMD[iL].mdOthersRec_aCut.init(Form("layer%dMD_mdOthersRec_aCut_pt_%s", iL, ext.c_str()), ptBins, 20, 0, 20);

    layerMD[iL].mdOthers2mm_aCutCloser.init(Form("layer%dMD_mdOthers2mm_aCutCloser_pt_%s", iL, ext.c_str()), ptBins, 10, 0, 10);
    layerMD[iL].mdOthersDes_aCutCloser.init(Form("layer%dMD_mdOthersDes_aCutCloser_pt_%s", iL, ext.c_str()), ptBins, 10, 0, 10);
    layerMD[iL].mdOthers2mm2SH_aCutCloser.init(Form("layer%dMD_mdOthers2mm2SH_aCutCloser_pt_%s", iL, ext.c_str()), ptBins, 10, 0, 10);
    layerMD[iL].mdOthersDes2SH_aCutCloser.init(Form("layer%dMD_mdOthersDes2SH_aCutCloser_pt_%s", iL, ext.c_str()), ptBins, 10, 0, 10);
    layerMD[iL].mdOthersAct_aCutCloser.init(Form("layer%dMD_mdOthersAct_aCutCloser_pt_%s", iL, ext.c_str()), ptBins, 10, 0, 10);
    layerMD[iL].mdOthersRec_aCutCloser.init(Form("layer%dMD_mdOthersRec_aCutCloser_pt_%s", iL, ext.c_str()), ptBins, 10, 0, 10);
  }
}

struct MiniDoublet {
  int pixL;
  int pixU;
  TVector3 r3;
  double alpha;
};

struct SuperDoublet {
  MiniDoublet mdRef;
  MiniDoublet mdOut;
  int iRef;
  int iOut;
  TVector3 r3; //may be different from plain mdRef.r3
  double alpha;
};

struct SDLink {
  SuperDoublet sdIn;
  SuperDoublet sdOut;
  int lIn;//layer
  int iIn;
  int lOut;//layer
  int iOut;
  double alpha; //point-to-point wrt radial
  double betaIn;//SD angle wrt point-to-point
  double betaOut;
  double pt;
};

struct MDStats {
  int nOthers;
  int nOthersRec;
  int mdOthers2mm_aCut;
  int mdOthers2mm_aCutCloser;
  int mdOthersDes_aCut;
  int mdOthersDes_aCutCloser;
  int mdOthersAct_aCut;
  int mdOthersAct_aCutCloser;
  int mdOthersRec_aCut;
  int mdOthersRec_aCutCloser;

  bool upperMatchToTP;
  bool upperMatchFull;
  bool lowerNTP;
  bool zFiducial;
  
  float miniDir;
  float miniDirAct;
  float miniDirRec;
  float miniCut;

  int pdgId;
};

void fillLayerMD_pt(HistoSet1D& layerMD, double pt, const MDStats& md){
  layerMD.den->Fill(pt);
  if (md.zFiducial){
    layerMD.denDesMDFid->Fill(pt);
    if (md.lowerNTP == 1) layerMD.denDesMDFid1TP->Fill(pt);
  }
  if (std::abs(md.miniDir) < md. miniCut) layerMD.num_aCut->Fill(pt);
  
  
  layerMD.others_inModule.fill(pt, md.nOthers);
  layerMD.othersRec_inModule.fill(pt, md.nOthersRec);
  
  layerMD.mdOthers2mm_aCut.fill(pt, md.mdOthers2mm_aCut);
  layerMD.mdOthers2mm_aCutCloser.fill(pt, md.mdOthers2mm_aCutCloser);
  layerMD.mdOthersDes_aCut.fill(pt, md.mdOthersDes_aCut);
  layerMD.mdOthersDes_aCutCloser.fill(pt, md.mdOthersDes_aCutCloser);
  if (md.upperMatchToTP && md.zFiducial){
    layerMD.denActAnySH->Fill(pt);
    if (md.lowerNTP == 1) layerMD.denActAnySH1TP->Fill(pt);
  }
  if (md.upperMatchFull && md.zFiducial){
    layerMD.denAct->Fill(pt);
    if (md.lowerNTP == 1 ) layerMD.denAct1TP->Fill(pt);
    if (std::abs(md.miniDir)    < md.miniCut) layerMD.num2SH_aCut->Fill(pt);
    if (std::abs(md.miniDirAct) < md.miniCut) layerMD.numAct_aCut->Fill(pt);
    if (std::abs(md.miniDirRec) < md.miniCut) layerMD.numRec_aCut->Fill(pt);
    
    layerMD.dminiDirActVs2mm->Fill(pt, std::min(0.199f, std::max(-0.199f,md.miniDir - md.miniDirAct)));//clamp at +/-1
    if (std::abs(md.pdgId) == 211 ){
      layerMD.dminiDirActVs2mm_pi->Fill(pt, std::min(0.199f, std::max(-0.199f,md.miniDir - md.miniDirAct)));//clamp at +/-1
    }

    layerMD.mdOthers2mm2SH_aCut.fill(pt, md.mdOthers2mm_aCut);
    layerMD.mdOthers2mm2SH_aCutCloser.fill(pt, md.mdOthers2mm_aCutCloser);
    layerMD.mdOthersDes2SH_aCut.fill(pt, md.mdOthersDes_aCut);
    layerMD.mdOthersDes2SH_aCutCloser.fill(pt, md.mdOthersDes_aCutCloser);
    
    layerMD.mdOthersAct_aCut.fill(pt, md.mdOthersAct_aCut);
    layerMD.mdOthersAct_aCutCloser.fill(pt, md.mdOthersAct_aCutCloser);
    layerMD.mdOthersRec_aCut.fill(pt, md.mdOthersRec_aCut);
    layerMD.mdOthersRec_aCutCloser.fill(pt, md.mdOthersRec_aCutCloser);
  }
}

TVector3 linePropagate(const TVector3& r3, const TVector3& p3, double rDest, int& status, bool useClosest = true, bool verbose = false){
  double rt = r3.Pt();
  double d = rDest - rt;

  double dotPR2D = r3.x()*p3.x() + r3.y()*p3.y();

  double pt = p3.Pt();
  double p =  p3.Mag();
  
  // r3 + p3/p*x*|d| = dest : dest.t = rt + d <=> rt^2 + 2*dotPR2D/p*x*|d| + pt^2/p^2*x^2*d^2 = rt^2 + 2*rt*d + d^2
  // 2*dotPR2D/p*|d|* x + pt^2/p^2*d^2* x^2 - ( 2*rt*d + d^2) = 0
  // 2*dotPR2D/p/|d|* x + pt^2/p^2* x^2 - ( 2*rt/d + 1) = 0
  // x^2 + 2*dotPR2D/p/|d|*(p/pt)^2* x  - ( 2*rt/d + 1)*(p/pt)^2 = 0
  // - dotPR2D/p/|d|*(p/pt)^2  +/- sqrt( (dotPR2D/p/|d|*(p/pt)^2)^2 + ( 2*rt/d + 1)*(p/pt)^2 )
  // (p/pt)*( - dotPR2D/pt/|d|  +/- sqrt( (dotPR2D/pt/|d| )^2 + ( 2*rt/d + 1) ) )
  double bb = dotPR2D/pt/std::abs(d);
  double dSign = useClosest ? 1. : -1.;
  double disc = bb*bb + (2.*rt/d + 1.);
  status = 0;
  if (disc < 0){
    status = 1;
    return r3;
  }
  double xx = (p/pt)*( dSign*sqrt(bb*bb + (2.*rt/d + 1.)) - bb);
  TVector3 dest = r3 + p3*(std::abs(d)/p)*xx;
  if (verbose || std::abs(dest.Pt() - rDest)>0.001){
    std::cout<<" "<<r3.Pt()<<" "<<r3.Phi()<<" "<<r3.z()<<" "<<pt<<" "<<p
	     <<" "<<d<<" "<<r3.x()*p3.x()<<" "<<r3.y()*p3.y()<<" "<<dotPR2D<<" "<<bb<<" "<<(2.*rt/d + 1.)<<" "<<bb*bb + (2.*rt/d + 1.)
	     <<" => "<<rDest
	     <<" => "<<dest.Pt()<<" "<<dest.Phi()<<" "<<dest.z()
	     <<std::endl;
  }
  return dest;

}

int ScanChainMiniDoublets( TChain* chain, int nEvents = -1, bool drawPlots = false) {

  std::vector<double> ptBins {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};

  constexpr int minLayer = 5;
  
  std::array<HistoSet1D, nLayers+1> layerMD_pt_all;
  createPtHistograms(layerMD_pt_all, "all", ptBins);
  
  std::array<HistoSet1D, nLayers+1> layerMD_pt_prim_all;
  createPtHistograms(layerMD_pt_prim_all, "prim_all", ptBins);

  std::array<HistoSet1D, nLayers+1> layerMD_pt_prim_tt;
  createPtHistograms(layerMD_pt_prim_tt, "prim_tt", ptBins);

  std::array<float, nLayers+1> miniDeltaBarrel {0, 0, 0, 0, 0,
      0.26, 0.16, 0.16, 0.18, 0.18, 0.18};

  //p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
  //(1 sigma based on 95.45% = 2sigma at 2 GeV)
  std::array<float, nLayers+1> miniMulsPtScale {0, 0, 0, 0, 0,
      0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};

  //mean of the horizontal layer position in y; treat this as R below
  std::array<float, nLayers+1> miniRminMean {0, 0, 0, 0, 0, //may want to fill these
      21.8, 34.6, 49.6, 67.4, 87.6, 106.8};

  std::map<int, std::array<float, 4> > moduleBoundaries;
  std::map<int, int > modulePopulation;
  std::array<float, 4> dbound {999,-999,999,-999}; //zmin, zmax, phimin, phimax
  std::array<float, 4>* cbound;
  std::array<float, 4> const* cboundC;
  std::array<float, 4> const* cboundL;
  std::array<float, 4> const* cboundH;

  std::array<int, nLayers+1> maxHitsBarrelLayer {};
  
  //geomRange loop
  {
    TObjArray *listOfFiles = chain->GetListOfFiles();
    unsigned int nEventsChain=0;
    if(nEvents==-1) 
      nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    
    // file loop
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    int dpop = 0;
    int* cpop = nullptr;
    
    while (( currentFile = (TFile*)fileIter.Next() )) {
      TFile f(currentFile->GetTitle());
      TTree *tree = (TTree*)f.Get("trkTree/tree");
      cms2.Init(tree);
      
      //Event Loop
      unsigned int nEvents = tree->GetEntries();
      for( unsigned int event = 0; event < nEvents && nEventsTotal < nEventsChain; ++event) {
	cms2.GetEntry(event);
	++nEventsTotal;

	int iidOld = -1;

	std::array<int, nLayers+1> hitsBarrelLayer {};
	auto nPix = pix_isBarrel().size();
	for (auto ipix = 0U; ipix < nPix; ++ipix){
	  if (pix_isBarrel()[ipix] == false) continue;
	  int lay = pix_lay()[ipix];

	  hitsBarrelLayer[lay]++;
	  int iid = pix_detId()[ipix];
	  if (iidOld != iid){
	    iidOld = iid;
	    if (modulePopulation.find(iid) == modulePopulation.end()){
	      modulePopulation[iid] = dpop;
	      moduleBoundaries[iid] = dbound;
	    }
	    cpop = &modulePopulation[iid];
	    cbound = &moduleBoundaries[iid];
	  }
	  
	  float z = pix_zsim()[ipix];
	  if (z==0) continue;
	  float phi = atan2(pix_ysim()[ipix], pix_xsim()[ipix]);
	  if ((*cbound)[0] > z) (*cbound)[0] = z;
	  if ((*cbound)[1] < z) (*cbound)[1] = z;
	  if (*cpop ==0){
	    (*cbound)[2] = phi;
	    (*cbound)[3] = phi;
	  } else {
	    if (sin((*cbound)[2]-phi) > 0) (*cbound)[2] = phi;
	    if (sin((*cbound)[3]-phi) < 0) (*cbound)[3] = phi;
	  }
	  (*cpop)++;
	}

	for (int i = 0; i<= nLayers; ++i){
	  maxHitsBarrelLayer[i] = std::max(maxHitsBarrelLayer[i], hitsBarrelLayer[i]);
	}
      }
    }
  }//geom range map loop
  

  cout<<__LINE__<<endl;
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while (( currentFile = (TFile*)fileIter.Next() )) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("trkTree/tree");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents && nEventsTotal < nEventsChain; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

      //keep track of sims per layer 
      std::array<std::set<int>, nLayers+1> simIdxInLayer;
      int iidStart = -1;
      int iidEnd = -1;
      int iidOld = -1;
      auto nPix = pix_isBarrel().size();
      for (auto ipix = 0U; ipix < nPix; ++ipix){
	if (pix_isBarrel()[ipix] == false) continue;
	int lay = pix_lay()[ipix];
	if (lay < 5 ) continue;

	int iid = pix_detId()[ipix];
	if (iidOld != iid){
	  iidStart = ipix;
	  iidOld = iid;
	  cboundC = &moduleBoundaries[iid];
	  if ((iid & 0x4)== 4){
	    cboundL = cboundC;
	    cboundH = &moduleBoundaries[iid+4];
	  }
	}
	if ((iid & 0x4)!= 4) continue; //only lower module

	TVector3 r3Rec(pix_x()[ipix], pix_y()[ipix], pix_z()[ipix]);

	TVector3 r3Sim(pix_xsim()[ipix], pix_ysim()[ipix], pix_zsim()[ipix]);
	float rs = r3Sim.Pt();
	TVector3 p3Sim(pix_pxsim()[ipix], pix_pysim()[ipix], pix_pzsim()[ipix]);
	float pts = p3Sim.Pt();
	float ps = p3Sim.Mag();
	if (rs == 0 || ps == 0 ) continue;
	if (pts < 0.8 )  continue;
	
	//	std::cout<<__LINE__<<" "<<ipix<<std::endl;
	int iSimIdx = pix_simTrkIdx()[ipix]; //tp index in ntuple (not a g4 trackId)
	if (simIdxInLayer[lay].find(iSimIdx) != simIdxInLayer[lay].end()) continue; //only one hit per layer per track
	else {
	  simIdxInLayer[lay].insert(iSimIdx);
	}

	MDStats md;

	int iProcess = pix_process()[ipix];
	int ibx = pix_bunchXing()[ipix];
	int iev = pix_event()[ipix];
	bool isPrimaryAny = (iProcess == 2 && ibx == 0);
	bool isPrimaryTT = (isPrimaryAny && iev == 0);
	int iParticle = pix_particle()[ipix];
	md.pdgId = iParticle;
	md.lowerNTP = pix_nSimTrk()[ipix];

	double dotPR2Ds = r3Sim.x()*p3Sim.x() + r3Sim.y()*p3Sim.y();
	double dir = dotPR2Ds > 0 ? 1. : -1.;

	//(r3Sim + p3Sim*(dir*d/ps)*x).pt2 = (rs + d)^2 = rs^2 + pts^2*d^2/ps^2*x^2 + 2*dotPR2Ds/ps*dir*d*x
	// rs^2 + 2*d*rs + d^2 = rs^2 + pts^2*d^2/ps^2*x^2 + 2*dotPR2Ds/ps*dir*d*x
	// 2*rs + d = pts^2*d/ps^2*x^2 + 2*dotPR2Ds/ps*dir*x
	// x^2 + 2*dotPR2Ds*ps*dir/pts^2/d*x - (2*rs + d)*ps^2/pts^2/d = 0
	// - dotPR2DsAbs*ps/pts^2/d +/- sqrt((dotPR2DsAbs*ps/pts^2/d)^2 + (2*rs + d)*ps^2/pts^2/d )
	// sqrt((dotPR2DsAbs*ps/pts^2/d)^2 + (2*rs + d)*ps^2/pts^2/d ) - dotPR2DsAbs*ps/pts^2/d
	// = (ps/pts)*( sqrt((dotPR2DsAbs/pts/d)^2 + (2*rs/d + 1.) ) - dotPR2DsAbs/pts/d )
	//
	double bbs = std::abs(dotPR2Ds/pts/0.2);
	double xxs = (ps/pts)*( sqrt(bbs*bbs + (2.*rs/0.2 + 1.)) - bbs);
	TVector3 nextR3Sim2mm = r3Sim + p3Sim*(dir*0.2/ps)*xxs;
	bbs = std::abs(dotPR2Ds/pts/miniDeltaBarrel[lay]);
	xxs = (ps/pts)*( sqrt(bbs*bbs + (2.*rs/miniDeltaBarrel[lay] + 1.)) - bbs);
	TVector3 nextR3SimDes = r3Sim + p3Sim*(dir*miniDeltaBarrel[lay]/ps)*xxs;
	TVector3 nextR3SimAct;// filled in the loop over the upper layer
	bool nextR3SimAct_isValid = false;
	TVector3 nextR3Rec;
	bool nextR3Rec_isValid = false;
	md.upperMatchToTP = false;
	md.upperMatchFull = false;
	md.zFiducial = false;
	if (nextR3SimDes.z() < (*cboundH)[1] - 0.1 && nextR3SimDes.z() > (*cboundH)[0] + 0.1){
	  md.zFiducial = true;
	  // if (isPrimaryTT && pts> 15){
	  //   std::cout<<"Is fiducial in lay "<< lay<<" : zLo "<<r3Sim.z()<<" zHi "<<nextR3SimDes.z()
	  // 	     <<" rLo "<<r3Sim.Pt()<<" rHi "<<nextR3SimDes.Pt()
	  // 	     <<" zLLLo "<<(*cboundL)[0]<<" zLLHi "<<(*cboundL)[1]
	  // 	     <<" zULLo "<<(*cboundH)[0]<<" zULHi "<<(*cboundH)[1]<<std::endl;
	  // }
	}
	
	std::vector<TVector3> otherR3Sim;        otherR3Sim.reserve(128);
	std::vector<TVector3> nextOtherR3Sim2mm; nextOtherR3Sim2mm.reserve(128);
	std::vector<TVector3> nextOtherR3SimDes; nextOtherR3SimDes.reserve(128);
	//these do not have any synchronization with the lower or upper layer or sim
	std::vector<TVector3> otherR3Rec;        otherR3Rec.reserve(128);
	std::vector<TVector3> nextOtherR3SimAct; nextOtherR3SimAct.reserve(128);
	std::vector<TVector3> nextOtherR3Rec;    nextOtherR3Rec.reserve(128);
	//	std::cout<<__LINE__<<" "<<ipix<<std::endl;

	for (auto jpix = iidStart; jpix < static_cast<int>(nPix); ++jpix){
	  if (jpix == static_cast<int>(ipix)) continue;
	  int jid = pix_detId()[jpix];
	  if (iid == jid){
	    TVector3 ar3(pix_x()[jpix], pix_y()[jpix], pix_z()[jpix]);
	    otherR3Rec.emplace_back(ar3);
	    
	    TVector3 ar3s(pix_xsim()[jpix], pix_ysim()[jpix], pix_zsim()[jpix]);
	    double ars = ar3s.Pt();
	    if (ars == 0) continue; //use just sim for sim-based predictions
	    TVector3 ap3s(pix_pxsim()[jpix], pix_pysim()[jpix], pix_pzsim()[jpix]);
	    double aps = ap3s.Mag();
	    double apts = ap3s.Pt();
	    
	    double dotOPR2Ds = ar3s.x()*ap3s.x() + ar3s.y()*ap3s.y();
	    double odir = dotOPR2Ds > 0 ? 1. : -1.;
	    bbs = std::abs(dotOPR2Ds/apts/0.2);
	    xxs = (aps/apts)*( sqrt(bbs*bbs + (2.*ars/0.2 + 1.)) - bbs);
	    TVector3 nr3s2mm = ar3s + ap3s*(dir*0.2/aps)*xxs;
	    bbs = std::abs(dotOPR2Ds/apts/miniDeltaBarrel[lay]);
	    xxs = (aps/apts)*( sqrt(bbs*bbs + (2.*ars/miniDeltaBarrel[lay] + 1.)) - bbs);	    
	    TVector3 nr3sDes = ar3s + ap3s*(dir*miniDeltaBarrel[lay]/aps)*xxs;
	    
	    otherR3Sim.emplace_back(ar3s);
	    nextOtherR3Sim2mm.emplace_back(nr3s2mm);
	    nextOtherR3SimDes.emplace_back(nr3sDes);
	  } else if (iid+4 == jid && lay>=5) {//upper module in the doublet pairs
	    int jSimIdx = pix_simTrkIdx()[jpix];
	    if (iSimIdx == jSimIdx){
	      int jProcess = pix_process()[jpix];
	      int jParticle = pix_particle()[jpix];
	      md.upperMatchToTP = true;
	      if (jSimIdx >= 0 && iProcess == jProcess && iParticle == jParticle){
		nextR3SimAct.SetXYZ(pix_xsim()[jpix], pix_ysim()[jpix], pix_zsim()[jpix]);
		nextR3SimAct_isValid = true;
		nextR3Rec.SetXYZ(pix_x()[jpix], pix_y()[jpix], pix_z()[jpix]);
		nextR3Rec_isValid = true;
		md.upperMatchFull = true;
		// if (pts > 10 && lay == 5 && isPrimaryTT){
		//   const float amdDir = (nextR3SimAct-r3Sim).DeltaPhi(r3Sim);
		//   TVector3 op3(pix_pxsim()[jpix], pix_pysim()[jpix], pix_pzsim()[jpix]);
		//   std::cout<<__LINE__<<":"<<nEventsTotal<<" L:"<<lay
		// 	   <<" ptI:"<<pts<<" ptJ:"<<op3.Pt()
		// 	   <<" idx:"<<iSimIdx<<" dir:"<<amdDir
		// 	   <<" procI:"<<iProcess<<" procJ:"<<pix_process()[jpix]
		// 	   <<" bxI:"<<ibx<<" bxJ:"<<pix_bunchXing()[jpix]
		// 	   <<" evI:"<<iev<<" evJ:"<<pix_event()[jpix]
		// 	   <<" typeI:"<<iParticle<<" typeJ:"<<pix_particle()[jpix]
		// 	   <<std::endl;
		// }
	      }
	    } else {
	      TVector3 ar3(pix_x()[jpix], pix_y()[jpix], pix_z()[jpix]);
	      nextOtherR3Rec.emplace_back(ar3);

	      TVector3 ar3s(pix_xsim()[jpix], pix_ysim()[jpix], pix_zsim()[jpix]);
	      //keep only valid sims
	      if (ar3s.Pt() != 0) nextOtherR3SimAct.emplace_back(ar3s);
	    }
	    
	  } else {
	    //once it's not equal to current or upper(+4), we never expect another entry in the same module
	    break;
	  }
	}//loop over the other pix hits jpix

	//	std::cout<<__LINE__<<" "<<ipix<<std::endl;
	
	//these two are the same for all layers
	const float ptCut = 1.0;
	const float miniSlope = rs/175.67/ptCut;
	
	const float miniMuls = miniMulsPtScale[lay]*3./ptCut;
	const float rLayNominal = lay >= 5? miniRminMean[lay] : 1e12;
	const float miniPVoff = 0.1/rLayNominal;
	const float miniCut = miniSlope + sqrt(miniMuls*miniMuls + miniPVoff*miniPVoff);
	md.miniCut = miniCut;
	
	md.miniDir = (nextR3Sim2mm-r3Sim).DeltaPhi(r3Sim);
	const float miniDirAbs = std::abs(md.miniDir);

	md.miniDirAct = nextR3SimAct_isValid ? (nextR3SimAct-r3Sim).DeltaPhi(r3Sim) : 999;
	const float miniDirActAbs = std::abs(md.miniDirAct);
	md.miniDirRec = nextR3Rec_isValid ? (nextR3Rec-r3Rec).DeltaPhi(r3Rec) : 999;
	const float miniDirRecAbs = std::abs(md.miniDirRec);
	
	//combine this hit with the others propagated to the next layer
	//and count what's inside the cut or have a smaller window that the correct combo

	md.mdOthers2mm_aCut = 0;
	md.mdOthers2mm_aCutCloser = 0;

	//check realL->otherU combinations
	for (auto or3s : nextOtherR3Sim2mm){
	  const float ominiDir = (or3s - r3Sim).DeltaPhi(r3Sim);
	  const float ominiDirAbs = std::abs(ominiDir);

	  if (ominiDirAbs < miniCut){
	    md.mdOthers2mm_aCut++;
	    if (ominiDirAbs < miniDirAbs){
	      md.mdOthers2mm_aCutCloser++;
	    }
	  }
	}
	//check otherL->realU combinations
	for (auto or3s : otherR3Sim){
	  const float ominiDir = (nextR3Sim2mm - or3s).DeltaPhi(or3s);
	  const float ominiDirAbs = std::abs(ominiDir);

	  if (ominiDirAbs < miniCut){
	    md.mdOthers2mm_aCut++;
	    if (ominiDirAbs < miniDirAbs){
	      md.mdOthers2mm_aCutCloser++;
	    }
	  }
	}

	md.mdOthersDes_aCut = 0;
	md.mdOthersDes_aCutCloser = 0;
	for (auto or3s : nextOtherR3SimDes){
	  const float ominiDir = (or3s - r3Sim).DeltaPhi(r3Sim);
	  const float ominiDirAbs = std::abs(ominiDir);
	  
	  if (ominiDirAbs < miniCut){
	    md.mdOthersDes_aCut++;
	    if (ominiDirAbs < miniDirAbs){
	      md.mdOthersDes_aCutCloser++;
	    }
	  }
	}
	for (auto or3s : otherR3Sim){
	  const float ominiDir = (nextR3SimDes - or3s).DeltaPhi(or3s);
	  const float ominiDirAbs = std::abs(ominiDir);
	  
	  if (ominiDirAbs < miniCut){
	    md.mdOthersDes_aCut++;
	    if (ominiDirAbs < miniDirAbs){
	      md.mdOthersDes_aCutCloser++;
	    }
	  }
	}

	//now collect using hits on an actual upper layer
	md.mdOthersAct_aCut = 0;
	md.mdOthersAct_aCutCloser = 0;
	for (auto or3s : nextOtherR3SimAct){
	  const float ominiDir = (or3s - r3Sim).DeltaPhi(r3Sim);
	  const float ominiDirAbs = std::abs(ominiDir);
	  
	  if (ominiDirAbs < miniCut){
	    md.mdOthersAct_aCut++;
	    if (nextR3SimAct_isValid && ominiDirAbs < miniDirActAbs){
	      md.mdOthersAct_aCutCloser++;
	    }
	  }
	}
	if (nextR3SimAct_isValid){
	  for (auto or3s : otherR3Sim){
	    const float ominiDir = (nextR3SimAct - or3s).DeltaPhi(or3s);
	    const float ominiDirAbs = std::abs(ominiDir);
	    
	    if (ominiDirAbs < miniCut){
	      md.mdOthersAct_aCut++;
	      if (ominiDirAbs < miniDirActAbs){
		md.mdOthersAct_aCutCloser++;
	      }
	    }
	  }
	}

	//now collect using rec hits on actual layers
	//(both of the denominator track hits are matched)
	md.mdOthersRec_aCut = 0;
	md.mdOthersRec_aCutCloser = 0;
	if (nextR3SimAct_isValid && nextR3Rec_isValid){
	  for (auto or3s : nextOtherR3Rec){
	    const float ominiDir = (or3s - r3Rec).DeltaPhi(r3Rec);
	    const float ominiDirAbs = std::abs(ominiDir);
	    
	    if (ominiDirAbs < miniCut){
	      md.mdOthersRec_aCut++;
	      if (ominiDirAbs < miniDirRecAbs){
		md.mdOthersRec_aCutCloser++;
	      }
	    }
	  }
	  for (auto or3s : otherR3Rec){
	    const float ominiDir = (nextR3Rec - or3s).DeltaPhi(or3s);
	    const float ominiDirAbs = std::abs(ominiDir);
	    
	    if (ominiDirAbs < miniCut){
	      md.mdOthersRec_aCut++;
	      if (ominiDirAbs < miniDirRecAbs){
		md.mdOthersRec_aCutCloser++;
	      }
	    }
	  }
	}

	//	std::cout<<__LINE__<<" "<<ipix<<std::endl;

	md.nOthers = otherR3Sim.size();
	md.nOthersRec = otherR3Rec.size();
	fillLayerMD_pt(layerMD_pt_all[lay], pts, md);
	if (isPrimaryAny){
	  fillLayerMD_pt(layerMD_pt_prim_all[lay], pts, md);
	}
	if (isPrimaryTT){
	  fillLayerMD_pt(layerMD_pt_prim_tt[lay], pts, md);
	}
	
	//	std::cout<<__LINE__<<" "<<ipix<<std::endl;
      }//loop over ipix


      // Progress feedback to the user
      if(nEventsTotal%1 == 0) {
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)) {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
          "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }//if(nEventsTotal%20000 == 0) {


    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  std::cout<<__LINE__<<" make efficiencies "<<std::endl;
  for (int iL = 1; iL <= nLayers; ++iL){
    createEfficiency(layerMD_pt_all[iL].num_aCut, layerMD_pt_all[iL].den, layerMD_pt_all[iL].eff_aCut);
    createEfficiency(layerMD_pt_prim_all[iL].num_aCut, layerMD_pt_prim_all[iL].den, layerMD_pt_prim_all[iL].eff_aCut);
    createEfficiency(layerMD_pt_prim_tt[iL].num_aCut, layerMD_pt_prim_tt[iL].den, layerMD_pt_prim_tt[iL].eff_aCut);

    createEfficiency(layerMD_pt_all[iL].denActAnySH, layerMD_pt_all[iL].denDesMDFid, layerMD_pt_all[iL].effActAnySH);
    createEfficiency(layerMD_pt_prim_all[iL].denActAnySH, layerMD_pt_prim_all[iL].denDesMDFid, layerMD_pt_prim_all[iL].effActAnySH);
    createEfficiency(layerMD_pt_prim_tt[iL].denActAnySH, layerMD_pt_prim_tt[iL].denDesMDFid, layerMD_pt_prim_tt[iL].effActAnySH);
    createEfficiency(layerMD_pt_all[iL].denAct, layerMD_pt_all[iL].denDesMDFid, layerMD_pt_all[iL].effAct);
    createEfficiency(layerMD_pt_prim_all[iL].denAct, layerMD_pt_prim_all[iL].denDesMDFid, layerMD_pt_prim_all[iL].effAct);
    createEfficiency(layerMD_pt_prim_tt[iL].denAct, layerMD_pt_prim_tt[iL].denDesMDFid, layerMD_pt_prim_tt[iL].effAct);

    createEfficiency(layerMD_pt_all[iL].denActAnySH1TP, layerMD_pt_all[iL].denDesMDFid1TP, layerMD_pt_all[iL].effActAnySH1TP);
    createEfficiency(layerMD_pt_prim_all[iL].denActAnySH1TP, layerMD_pt_prim_all[iL].denDesMDFid1TP, layerMD_pt_prim_all[iL].effActAnySH1TP);
    createEfficiency(layerMD_pt_prim_tt[iL].denActAnySH1TP, layerMD_pt_prim_tt[iL].denDesMDFid1TP, layerMD_pt_prim_tt[iL].effActAnySH1TP);
    createEfficiency(layerMD_pt_all[iL].denAct1TP, layerMD_pt_all[iL].denDesMDFid1TP, layerMD_pt_all[iL].effAct1TP);
    createEfficiency(layerMD_pt_prim_all[iL].denAct1TP, layerMD_pt_prim_all[iL].denDesMDFid1TP, layerMD_pt_prim_all[iL].effAct1TP);
    createEfficiency(layerMD_pt_prim_tt[iL].denAct1TP, layerMD_pt_prim_tt[iL].denDesMDFid1TP, layerMD_pt_prim_tt[iL].effAct1TP);

    createEfficiency(layerMD_pt_all[iL].num2SH_aCut, layerMD_pt_all[iL].denAct, layerMD_pt_all[iL].eff2SH_aCut);
    createEfficiency(layerMD_pt_prim_all[iL].num2SH_aCut, layerMD_pt_prim_all[iL].denAct, layerMD_pt_prim_all[iL].eff2SH_aCut);
    createEfficiency(layerMD_pt_prim_tt[iL].num2SH_aCut, layerMD_pt_prim_tt[iL].denAct, layerMD_pt_prim_tt[iL].eff2SH_aCut);

    createEfficiency(layerMD_pt_all[iL].numAct_aCut, layerMD_pt_all[iL].denAct, layerMD_pt_all[iL].effAct_aCut);
    createEfficiency(layerMD_pt_prim_all[iL].numAct_aCut, layerMD_pt_prim_all[iL].denAct, layerMD_pt_prim_all[iL].effAct_aCut);
    createEfficiency(layerMD_pt_prim_tt[iL].numAct_aCut, layerMD_pt_prim_tt[iL].denAct, layerMD_pt_prim_tt[iL].effAct_aCut);

    createEfficiency(layerMD_pt_all[iL].numRec_aCut, layerMD_pt_all[iL].denAct, layerMD_pt_all[iL].effRec_aCut);
    createEfficiency(layerMD_pt_prim_all[iL].numRec_aCut, layerMD_pt_prim_all[iL].denAct, layerMD_pt_prim_all[iL].effRec_aCut);
    createEfficiency(layerMD_pt_prim_tt[iL].numRec_aCut, layerMD_pt_prim_tt[iL].denAct, layerMD_pt_prim_tt[iL].effRec_aCut);
  }


  if (drawPlots){
    std::cout<<__LINE__<<" draw and print "<<std::endl;
    for (int iL = 5; iL <= nLayers; ++iL){
      if (iL != 5 && iL != 10) continue;

      std::array<std::string, 3> lptypes {"All sim hits", "All from primaries", "All from t#bar{t}"};
      plotEffOverlay({layerMD_pt_all[iL].eff_aCut, layerMD_pt_prim_all[iL].eff_aCut, layerMD_pt_prim_tt[iL].eff_aCut},
		     Form("#alpha(1 GeV) cut efficiency in layer %d; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_all_%d_vs_prim.png", iL), lptypes, 0.0, 1.1);

      plotEffOverlay({layerMD_pt_all[iL].effAct, layerMD_pt_prim_all[iL].effAct, layerMD_pt_prim_tt[iL].effAct},
		     Form("Efficiency to match lower-upper in layer %d; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_all_%d_vs_prim_effAct.png", iL), lptypes, 0.0, 1.1);
      plotEffOverlay({layerMD_pt_all[iL].effActAnySH, layerMD_pt_prim_all[iL].effActAnySH, layerMD_pt_prim_tt[iL].effActAnySH},
		     Form("Efficiency to match lower-upper to TP in layer %d; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_all_%d_vs_prim_effActAnySH.png", iL), lptypes, 0.0, 1.1);

      plotEffOverlay({layerMD_pt_all[iL].effAct1TP, layerMD_pt_prim_all[iL].effAct1TP, layerMD_pt_prim_tt[iL].effAct1TP},
		     Form("Efficiency to match (1TP) lower-upper in layer %d; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_all_%d_vs_prim_effAct1TP.png", iL), lptypes, 0.0, 1.1);
      plotEffOverlay({layerMD_pt_all[iL].effActAnySH1TP, layerMD_pt_prim_all[iL].effActAnySH1TP, layerMD_pt_prim_tt[iL].effActAnySH1TP},
		     Form("Efficiency to match (1TP) lower-upper to TP in layer %d; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_all_%d_vs_prim_effActAnySH1TP.png", iL), lptypes, 0.0, 1.1);

      plotEffOverlay({layerMD_pt_all[iL].eff2SH_aCut, layerMD_pt_all[iL].effAct_aCut, layerMD_pt_all[iL].effRec_aCut},
		     Form("#alpha(1 GeV) cut efficiency in layer %d for all; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_all_%d_mock_sameDenAct_vs_act_rec.png", iL),
		     {"Sim hit r3 + p3", "Sim hit r3Lower-r3Upper", "Rec hit r3Lower-r3Upper"}, 0.0, 1.1);
      plotEffOverlay({layerMD_pt_prim_all[iL].eff2SH_aCut, layerMD_pt_prim_all[iL].effAct_aCut, layerMD_pt_prim_all[iL].effRec_aCut},
		     Form("#alpha(1 GeV) cut efficiency in layer %d for all primary; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_prim_all_%d_mock_sameDenAct_vs_act_rec.png", iL),
		     {"Sim hit r3 + p3", "Sim hit r3Lower-r3Upper", "Rec hit r3Lower-r3Upper"}, 0.0, 1.1);
      plotEffOverlay({layerMD_pt_prim_tt[iL].eff2SH_aCut, layerMD_pt_prim_tt[iL].effAct_aCut, layerMD_pt_prim_tt[iL].effRec_aCut},
		     Form("#alpha(1 GeV) cut efficiency in layer %d for t#bar{t}; p_{T}^{SimH}, GeV; Efficiency", iL),
		     Form("layerMD_pt_prim_tt_%d_mock_sameDenAct_vs_act_rec.png", iL),
		     {"Sim hit r3 + p3", "Sim hit r3Lower-r3Upper", "Rec hit r3Lower-r3Upper"}, 0.0, 1.1);

      {
	auto h2 = layerMD_pt_all[iL].dminiDirActVs2mm;
	h2->SetTitle(Form("#Delta Dir (true, mock) for layer %d for all;  p_{T}^{SimH}, GeV; #Delta Dir (true, mock)", iL));
	auto cn = h2->GetTitle();
	TCanvas* cv = new TCanvas(cn, cn, 600, 600);
	cv->cd();
	gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
	h2->Draw("colz");
	h2->SetStats(0);
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetLogz();
	
	gPad->SaveAs(Form("layerMD_pt_all_%d_sameDenAct_mock_vs_act.png", iL) );
      }
      {
	auto h2 = layerMD_pt_prim_all[iL].dminiDirActVs2mm;
	h2->SetTitle(Form("#Delta Dir (true, mock) for layer %d for all primary;  p_{T}^{SimH}, GeV; #Delta Dir (true, mock)", iL));
	auto cn = h2->GetTitle();
	TCanvas* cv = new TCanvas(cn, cn, 600, 600);
	cv->cd();
	gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
	h2->Draw("colz");
	h2->SetStats(0);
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetLogz();
	
	gPad->SaveAs(Form("layerMD_pt_prim_all_%d_sameDenAct_mock_vs_act.png", iL) );
      }
      {
	auto h2 = layerMD_pt_prim_tt[iL].dminiDirActVs2mm;
	h2->SetTitle(Form("#Delta Dir (true, mock) for layer %d for t#bar{t};  p_{T}^{SimH}, GeV; #Delta Dir (true, mock)", iL));
	auto cn = h2->GetTitle();
	TCanvas* cv = new TCanvas(cn, cn, 600, 600);
	cv->cd();
	gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
	h2->Draw("colz");
	h2->SetStats(0);
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetLogz();
	
	gPad->SaveAs(Form("layerMD_pt_prim_tt_%d_sameDenAct_mock_vs_act.png", iL) );
      }
      
      {
	auto h2 = layerMD_pt_all[iL].dminiDirActVs2mm_pi;
	h2->SetTitle(Form("#Delta Dir (true, mock) for layer %d for all #pi;  p_{T}^{SimH}, GeV; #Delta Dir (true, mock)", iL));
	auto cn = h2->GetTitle();
	TCanvas* cv = new TCanvas(cn, cn, 600, 600);
	cv->cd();
	gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
	h2->Draw("colz");
	h2->SetStats(0);
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetLogz();
	
	gPad->SaveAs(Form("layerMD_pt_all_pi_%d_sameDenAct_mock_vs_act.png", iL) );
      }
      {
	auto h2 = layerMD_pt_prim_all[iL].dminiDirActVs2mm_pi;
	h2->SetTitle(Form("#Delta Dir (true, mock) for layer %d for all primary #pi;  p_{T}^{SimH}, GeV; #Delta Dir (true, mock)", iL));
	auto cn = h2->GetTitle();
	TCanvas* cv = new TCanvas(cn, cn, 600, 600);
	cv->cd();
	gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
	h2->Draw("colz");
	h2->SetStats(0);
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetLogz();
	
	gPad->SaveAs(Form("layerMD_pt_prim_all_pi_%d_sameDenAct_mock_vs_act.png", iL) );
      }
      {
	auto h2 = layerMD_pt_prim_tt[iL].dminiDirActVs2mm_pi;
	h2->SetTitle(Form("#Delta Dir (true, mock) for layer %d for t#bar{t} #pi;  p_{T}^{SimH}, GeV; #Delta Dir (true, mock)", iL));
	auto cn = h2->GetTitle();
	TCanvas* cv = new TCanvas(cn, cn, 600, 600);
	cv->cd();
	gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
	h2->Draw("colz");
	h2->SetStats(0);
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetLogz();
	
	gPad->SaveAs(Form("layerMD_pt_prim_tt_pi_%d_sameDenAct_mock_vs_act.png", iL) );
      }
      
      
      std::array<std::pair<HistoSet1D const &, std::string>,3 > oinmos = {
	std::make_pair(layerMD_pt_all[iL], "all"),  std::make_pair(layerMD_pt_prim_all[iL], "all primary"),
	std::make_pair(layerMD_pt_prim_tt[iL], "t#bar{t}") };
      for (auto hset : oinmos){
	plot2DProf(hset.first.others_inModule,
		   Form("N other hits in module layer %d for %s; p_{T}^{SimH}, GeV; N other hits", iL, hset.second.c_str()),
		   0.1, 0x7, false);
	plot2DProf(hset.first.othersRec_inModule,
		   Form("N other recHits in module layer %d for %s; p_{T}^{SimH}, GeV; N other hits", iL, hset.second.c_str()),
		   0.1, 0x7, false);
	plot2DProf(hset.first.mdOthers2mm_aCut,
		   Form("N real-other MDs(2 mm) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N other hits", iL, hset.second.c_str()),
		   0.01, 0x7, false);
	plot2DProf(hset.first.mdOthersDes_aCut,
		   Form("N real-other MDs(design) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.01, 0x7, false);
	plot2DProf(hset.first.mdOthers2mm2SH_aCut,
		   Form("N real-other MDs(2 mm 2SH) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N other hits", iL, hset.second.c_str()),
		   0.01, 0x7, false);
	plot2DProf(hset.first.mdOthersDes2SH_aCut,
		   Form("N real-other MDs(design 2SH) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.01, 0x7, false);
	plot2DProf(hset.first.mdOthersAct_aCut,
		   Form("N real-other MDs(actual) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.01, 0x7, false);
	plot2DProf(hset.first.mdOthersRec_aCut,
		   Form("N real-other MDs(reco) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.01, 0x7, false);
	plot2DProf(hset.first.mdOthers2mm_aCutCloser,
		   Form("N real-other-closer MDs(2 mm) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.001, 0x7, false);
	plot2DProf(hset.first.mdOthersDes_aCutCloser,
		   Form("N real-other-closer MDs(design) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.001, 0x7, false);
	plot2DProf(hset.first.mdOthers2mm2SH_aCutCloser,
		   Form("N real-other-closer MDs(2 mm 2SH) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.001, 0x7, false);
	plot2DProf(hset.first.mdOthersDes2SH_aCutCloser,
		   Form("N real-other-closer MDs(design 2SH) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.001, 0x7, false);
	plot2DProf(hset.first.mdOthersAct_aCutCloser,
		   Form("N real-other-closer MDs(actual) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.001, 0x7, false);
	plot2DProf(hset.first.mdOthersRec_aCutCloser,
		   Form("N real-other-closer MDs(reco) passing #alpha on layer %d for %s; p_{T}^{SimH}, GeV; N real-other MDs", iL, hset.second.c_str()),
		   0.001, 0x7, false);
      }

    }
  }
  
  std::cout<<__LINE__<<" write to file "<<std::endl;
  TFile* outHistograms = new TFile("outHistograms.root", "RECREATE");
  for (int iL = 1; iL <= nLayers; ++iL){
    layerMD_pt_all[iL].write(outHistograms);
    layerMD_pt_prim_all[iL].write(outHistograms);
    layerMD_pt_prim_tt[iL].write(outHistograms);
  }
  outHistograms->Write();
  outHistograms->Close();
  std::cout<<__LINE__<<" Done "<<std::endl;

  return 0;
}

int ScanChainMockSuperDoublets( TChain* chain, int nEvents = -1, bool drawPlots = false) {

  std::vector<double> ptBins {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};

  constexpr int minLayer = 5;
  
  std::array<HistoSet1D, nLayers+1> layerMD_pt_all;
  createPtHistograms(layerMD_pt_all, "all", ptBins);
  
  std::array<HistoSet1D, nLayers+1> layerMD_pt_prim_all;
  createPtHistograms(layerMD_pt_prim_all, "prim_all", ptBins);

  std::array<HistoSet1D, nLayers+1> layerMD_pt_prim_tt;
  createPtHistograms(layerMD_pt_prim_tt, "prim_tt", ptBins);

  std::array<float, nLayers+1> miniDeltaBarrel {0, 0, 0, 0, 0,
      0.26, 0.16, 0.16, 0.18, 0.18, 0.18};

  //p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
  //(1 sigma based on 95.45% = 2sigma at 2 GeV)
  std::array<float, nLayers+1> miniMulsPtScale {0, 0, 0, 0, 0,
      0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};

  //mean of the horizontal layer position in y; treat this as R below
  std::array<float, nLayers+1> miniRminMean {0, 0, 0, 0, 0, //may want to fill these
      21.8, 34.6, 49.6, 67.4, 87.6, 106.8};

  //  cout<<__LINE__<<endl;
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  constexpr int maxHitsInLayer = 50000;
  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while (( currentFile = (TFile*)fileIter.Next() )) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("trkTree/tree");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents && nEventsTotal < nEventsChain; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

      //this could be filled on-demand and in some regions of interest at some point
      //assume somewhat hermetic tracker and the only hits in the barrel to be from the barrel
      //reference layer minidoublet
      std::array<std::vector<std::pair<int, TVector3> >, nLayers+1 > mockLayerMDfwRefLower;
      std::array<std::vector<std::pair<int, TVector3> >, nLayers+1 > mockLayerMDfwRefUpper;
      //offset layer minidoublet (2cm hardcode for now)
      std::array<std::vector<std::pair<int, TVector3> >, nLayers+1 > mockLayerMDfwD2cmLower;
      std::array<std::vector<std::pair<int, TVector3> >, nLayers+1 > mockLayerMDfwD2cmUpper;
      for (int iL = 1; iL <=nLayers; ++iL){
	mockLayerMDfwRefLower[iL].reserve(maxHitsInLayer);
	mockLayerMDfwRefUpper[iL].reserve(maxHitsInLayer);
	mockLayerMDfwD2cmLower[iL].reserve(maxHitsInLayer);
	mockLayerMDfwD2cmUpper[iL].reserve(maxHitsInLayer);
      }

      std::array<int, nLayers+1> nHitsLayer1GeV {};
      std::array<int, nLayers+1> nHitsLayer2GeV {};
      
      std::cout<<"Load hits "<<std::endl;
      std::array<std::set<int>, nLayers+1> simIdxDeltaInLayer;
      std::array<std::set<int>, nLayers+1> simIdxInLayer;
      auto nPix = pix_isBarrel().size();
      for (auto ipix = 0U; ipix < nPix; ++ipix){
	if (pix_isBarrel()[ipix] == false) continue;
	int lay = pix_lay()[ipix];
	if (lay < 5 ) continue;

	int iid = pix_detId()[ipix];
	if ((iid & 0x4)!= 4) continue; 

	TVector3 r3Sim(pix_xsim()[ipix], pix_ysim()[ipix], pix_zsim()[ipix]);
	if (r3Sim.x() == 0 && r3Sim.y() == 0) continue; //use only hits with sim info
	
	float rs = r3Sim.Pt();
	TVector3 p3Sim(pix_pxsim()[ipix], pix_pysim()[ipix], pix_pzsim()[ipix]);
	float pts = p3Sim.Pt();
	float ps = p3Sim.Mag();

	int iParticle = pix_particle()[ipix];

	int iSimIdx = pix_simTrkIdx()[ipix]; //tp index in ntuple (not a g4 trackId)

	//FIXME: consider to improve to allow multiple deltas
	if (iParticle == -11 && ps < 0.1){
	  if (simIdxDeltaInLayer[lay].find(iSimIdx) != simIdxDeltaInLayer[lay].end()) continue; //only one hit per layer per track
	  else {
	    simIdxDeltaInLayer[lay].insert(iSimIdx);
	  }
	} else {
	  if (simIdxInLayer[lay].find(iSimIdx) != simIdxInLayer[lay].end()) continue; //only one hit per layer per track
	  else {
	    simIdxInLayer[lay].insert(iSimIdx);
	  }
	}

	if (pts >1) nHitsLayer1GeV[lay]++;
	if (pts >2) nHitsLayer2GeV[lay]++;
	
	const double sdOffset = 2.0;
	const double mdOffset = miniDeltaBarrel[lay];
	const double rNominal = miniRminMean[lay]+1;//this is going to be in between the sector radii

	const double rRefLower = rNominal;
	const double rRefUpper = rNominal+mdOffset;
	const double rSDfwLower = rRefLower + sdOffset;
	const double rSDfwUpper = rRefUpper + sdOffset;

	int pstat = 0;
	auto r3RefLower = linePropagate(r3Sim, p3Sim, rRefLower, pstat);
	if (pstat == 0) mockLayerMDfwRefLower[lay].push_back(std::make_pair(ipix,r3RefLower));
	auto r3RefUpper = linePropagate(r3Sim, p3Sim, rRefUpper, pstat);
	if (pstat == 0) mockLayerMDfwRefUpper[lay].push_back(std::make_pair(ipix, r3RefUpper));
	auto r3SDfwLower = linePropagate(r3Sim, p3Sim, rSDfwLower, pstat);
	if (pstat == 0) mockLayerMDfwD2cmLower[lay].push_back(std::make_pair(ipix, r3SDfwLower));
	auto r3SDfwUpper = linePropagate(r3Sim, p3Sim, rSDfwUpper, pstat);
	if (pstat == 0) mockLayerMDfwD2cmUpper[lay].push_back(std::make_pair(ipix, r3SDfwUpper));


      }//nPix: filling mock MDs

      
      constexpr int maxMDexpected = 100;
      std::array<std::vector<MiniDoublet>, nLayers+1> mockLayerMDfwRef;
      std::array<std::vector<MiniDoublet>, nLayers+1> mockLayerMDfwD2cm;

      std::array<std::vector<SuperDoublet>, nLayers+1> mockLayerSDfwD2cm;

      std::vector<SDLink> mockLayer5to7SDLfwD2cm;
      std::vector<SDLink> mockLayer7to9SDLfwD2cm;
      
      std::cout<<"Combine MDs "<<std::endl;
      int nHitsTried = 0;
      for (int iL = minLayer; iL <= nLayers; ++iL){
	auto const&  hitsRefLower = mockLayerMDfwRefLower[iL];
	auto const&  hitsRefUpper = mockLayerMDfwRefUpper[iL];
	auto const&  hitsOutLower = mockLayerMDfwD2cmLower[iL];
	auto const&  hitsOutUpper = mockLayerMDfwD2cmUpper[iL];
	
	auto& mockMDfwRef = mockLayerMDfwRef[iL];

	auto mdCombine = [&] (decltype(hitsRefLower) hitsL, decltype(hitsRefUpper) hitsU, decltype(mockMDfwRef) mDs){
	  for (auto hL : hitsL) {
	    nHitsTried++; //if (nHitsTried>1) exit(0);
	    float rt = hL.second.Pt();
	    const float ptCut = 1.0;
	    const float miniSlope = rt/175.67/ptCut;
	    
	    const float miniMuls = miniMulsPtScale[iL]*3./ptCut;
	    const float rLayNominal = iL >= minLayer ? miniRminMean[iL] : 1e12;
	    const float miniPVoff = 0.1/rLayNominal;
	    const float miniCut = miniSlope + sqrt(miniMuls*miniMuls + miniPVoff*miniPVoff);
	    
	    const float dzCut = 10.;//may want to adjust this: PS modules are shorter
	    
	    for (auto hU : hitsU) {
	      float dz = hL.second.z() - hU.second.z();
	      if (std::abs(dz) > dzCut) continue;
	      
	      float dPhi = hL.second.DeltaPhi(hU.second-hL.second);
	      if (std::abs(dPhi) > miniCut) continue;
	      
	      MiniDoublet md;
	      md.pixL = hL.first;
	      md.pixU = hU.first;
	      md.r3 = hL.second;
	      md.alpha = dPhi;
	      mDs.push_back(md);
	    }
	  }
	};
	mdCombine(hitsRefLower, hitsRefUpper, mockMDfwRef);

	auto& mockMDfwD2cm = mockLayerMDfwD2cm[iL];
	mdCombine(hitsOutLower, hitsOutUpper, mockMDfwD2cm);

	//now make super-doublets
	auto& mockSDfwD2cm = mockLayerSDfwD2cm[iL];

	int iRef = -1;
	for (auto mdRef : mockMDfwRef){
	  iRef++;
	  float rtRef = mdRef.r3.Pt();
	  float zRef = mdRef.r3.z();
	  //
	  int iOut = -1;
	  for (auto mdOut : mockMDfwD2cm){
	    iOut++;
	    float rtOut = mdOut.r3.Pt();
	    float zOut = mdOut.r3.z();

	    //apply some loose Z compatibility
	    float zLo = rtOut/rtRef*(zRef - 15.) - 10.; //15 for the luminous ; 10 for module size
	    float zHi = rtOut/rtRef*(zRef + 15.) + 10.;
	    if (zOut < zLo || zOut > zHi) continue;
	    
	    SuperDoublet sd;
	    sd.iRef = iRef;
	    sd.iOut = iOut;

	    float rt = 0.5*(rtRef + rtOut); //take the middle: it matches better the point-to-point
	    const float ptCut = 1.0;
	    const float sdSlope = rt/175.67/ptCut;
	    const float sdMuls = miniMulsPtScale[iL]*3./ptCut*2.;//will need a better guess than x2?
	    const float sdPVoff = 0.1/rt;
	    const float sdCut = sdSlope + sqrt(sdMuls*sdMuls + sdPVoff*sdPVoff);

	    //plain SD bend cut
	    float dPhi = mdRef.r3.DeltaPhi(mdOut.r3 - mdRef.r3);
	    if (std::abs(dPhi) > sdCut ) continue;

	    
	    sd.r3 = mdRef.r3;
	    sd.alpha = dPhi;
	    //loose angle compatibility
	    float dAlpha_Bfield = (rtOut - rtRef)/175.67/ptCut;
	    float dAlpha_res = 0.04/miniDeltaBarrel[iL];//4-strip difference
	    float dAlpha_compat = dAlpha_Bfield + dAlpha_res;
	    if (std::abs(mdRef.alpha- sd.alpha) > dAlpha_compat) continue;
	    if (std::abs(mdOut.alpha- sd.alpha) > dAlpha_compat) continue;
	    if (std::abs(mdOut.alpha- mdRef.alpha) > dAlpha_compat) continue;

	    sd.mdRef = mdRef;
	    sd.mdOut = mdOut;

	    if ( (sd.mdRef.r3 - sd.mdOut.r3).Pt() > 1.5*(rtOut - rtRef)){
	      //problem in matching
	      std::cout<<__LINE__
		       <<" "<<sd.mdRef.r3.Pt()<<" "<<sd.mdRef.r3.Phi()
		       <<" "<<sd.mdOut.r3.Pt()<<" "<<sd.mdOut.r3.Phi()
		       <<std::endl;
	    }
	    mockSDfwD2cm.push_back(sd);
	  }
	}
      }//iL

      //try to link segments: do 5-7 and 7-9
      auto sdLink = [&] (int lIn, int lOut, decltype(mockLayer5to7SDLfwD2cm)& sdlV){
	auto const& sdInV  = mockLayerSDfwD2cm[lIn];
	auto const& sdOutV = mockLayerSDfwD2cm[lOut];

	int nAll = 0;
	int nDeltaZ = 0;
	int nSlope = 0;
	int nInAlphaCompat = 0;
	int nOutAlphaCompat = 0;
	int ndBeta = 0;
	
	int iIn = -1;
	for ( auto sdIn : sdInV ) {
	  iIn++;
	  //
	  float rtIn = sdIn.r3.Pt();
	  float zIn = sdIn.r3.z();

	  int iOut = -1;
	  for ( auto sdOut : sdOutV ) {
	    iOut++;
	    //
	    nAll++;
	    
	    float rtOut = sdOut.r3.Pt();
	    float zOut = sdOut.r3.z();
	    //apply some loose Z compatibility
	    //FIXME: refine using inner layer directions (can prune later)
	    float zLo = rtOut/rtIn*(zIn - 15.) - 10.; //15 for the luminous ; 10 for module size
	    float zHi = rtOut/rtIn*(zIn + 15.) + 10.;
	    if (zOut < zLo || zOut > zHi) continue;
	    nDeltaZ++;

	    auto midR3 = 0.5*(sdIn.r3 + sdOut.r3);
	    double dPhi = midR3.DeltaPhi(sdOut.r3 - sdIn.r3);
	    double rt = 0.5*(sdIn.r3.Pt() + sdOut.r3.Pt());
	    const float ptCut = 1.0;
	    const float sdlSlope = rt/175.67/ptCut;
	    const float sdlThetaMulsF = 0.015*sqrt(0.2);
	    const float sdlMuls = sdlThetaMulsF*3./ptCut*4;//will need a better guess than x4?
	    const float sdlPVoff = 0.1/rt;
	    const float sdlCut = sdlSlope + sqrt(sdlMuls*sdlMuls + sdlPVoff*sdlPVoff);
	    
	    if (std::abs(dPhi) > sdlCut ) continue;
	    nSlope++;
	    
	    double betaIn = sdIn.alpha - sdIn.r3.DeltaPhi(sdOut.r3 - sdIn.r3);
	    double betaOut = - sdOut.alpha + sdOut.r3.DeltaPhi(sdOut.r3 - sdIn.r3); //to match sign for correct match
	    
	    //loose angle compatibility
	    float dAlpha_Bfield = (rtOut - rtIn)/175.67/ptCut;
	    float dSDIn = sdIn.mdOut.r3.Pt() - sdIn.mdRef.r3.Pt();
	    float dSDOut = sdOut.mdOut.r3.Pt() - sdOut.mdRef.r3.Pt();	
	    float dAlpha_res = 0.02/std::min(dSDIn, dSDOut);//2-strip difference; use the smallest SD separation
	    float dAlpha_compat = dAlpha_Bfield + dAlpha_res;
	    if (std::abs(sdIn.alpha- dPhi) > dAlpha_compat) continue;
	    nInAlphaCompat++;
	    if (std::abs(sdOut.alpha- dPhi) > dAlpha_compat) continue;
	    nOutAlphaCompat++;
	    
	    //now the actual segment linking magic
	    float betaAv = 0.5*(betaIn + betaOut);
	    float dr = (sdOut.r3 - sdIn.r3).Mag();
	    //pt*175.67/2. = R
	    //R*sin(betaAv) = pt*175.67/2*sin(betaAv) = dr/2 => pt = dr/175.67/sin(betaAv);
	    float pt_beta = dr/175.67/sin(betaAv);
	    float dBetaRes = dAlpha_res;
	    float dBetaMuls = sdlThetaMulsF*3./std::min(pt_beta, 7.0f);//need to confirm the range-out value of 7 GeV
	    float dBetaCut = sqrt(dBetaRes*dBetaRes*2.0 + dBetaMuls*dBetaMuls);
	    if (std::abs(betaIn - betaOut) > dBetaCut) continue;
	    ndBeta++;
	    
	    SDLink sdl;
	    sdl.sdIn = sdIn;
	    sdl.sdOut = sdOut;
	    sdl.lIn = lIn;
	    sdl.iIn = iIn;
	    sdl.lOut = lOut;
	    sdl.iOut = iOut;
	    
	    sdl.alpha = dPhi;
	    sdl.betaIn = betaIn;
	    sdl.betaOut = betaOut;
	    sdl.pt = pt_beta;
	    
	    sdlV.push_back(sdl);
	  }//sdOutV
	}//sdInV
	
	std::cout<<"SD links stat "<<lIn<<"-"<<lOut
	         <<" nAll "<<nAll<<" nDeltaZ "<<nDeltaZ<<" nSlope "<<nSlope
         	 <<" nInAlphaCompat "<< nInAlphaCompat<<" nOutAlphaCompat "<<nOutAlphaCompat
	         <<" ndBeta "<<ndBeta
	         <<" final "<<sdlV.size()
	         <<std::endl;
      };//auto sdLink

      sdLink(5, 7, mockLayer5to7SDLfwD2cm);
      sdLink(7, 9, mockLayer7to9SDLfwD2cm);

      
      std::cout<<"Print stats"<<std::endl;
      for (int iL = minLayer; iL <= nLayers; ++iL){
	int nSDs1GeV = 0;
	int nSDs1GeVMatch4 = 0;
	int nSDs2GeV = 0;
	int nSDs2GeVMatch4 = 0;
	for (auto sd : mockLayerSDfwD2cm[iL]){
	  int ipix = sd.mdRef.pixL;
	  TVector3 p3RL(pix_pxsim()[ipix], pix_pysim()[ipix], pix_pzsim()[ipix]);
	  ipix = sd.mdRef.pixU;
	  TVector3 p3RU(pix_pxsim()[ipix], pix_pysim()[ipix], pix_pzsim()[ipix]);
	  ipix = sd.mdOut.pixL;
	  TVector3 p3OL(pix_pxsim()[ipix], pix_pysim()[ipix], pix_pzsim()[ipix]);
	  ipix = sd.mdOut.pixU;
	  TVector3 p3OU(pix_pxsim()[ipix], pix_pysim()[ipix], pix_pzsim()[ipix]);

	  if (p3RL.Pt() > 1 || p3RU.Pt() > 1 || p3OL.Pt() > 1 || p3OU.Pt() > 1){
	    nSDs1GeV++;
	    if (p3RL.Pt() > 2 || p3RU.Pt() > 2 || p3OL.Pt() > 2 || p3OU.Pt() > 2){
	      nSDs2GeV++;
	    }
	  }
	  if (sd.mdRef.pixL == sd.mdRef.pixU
	      && sd.mdOut.pixL == sd.mdOut.pixU
	      && sd.mdRef.pixL == sd.mdOut.pixL){
	    if ( p3RL.Pt() > 1 ) nSDs1GeVMatch4++;
	    if ( p3RL.Pt() > 2 ) nSDs2GeVMatch4++;
	  }
	}
	std::cout<<"Summary for layer "<<iL
		 <<" h1GeV "<<nHitsLayer1GeV[iL]
		 <<" h2GeV "<<nHitsLayer2GeV[iL]
		 <<" refLos "<<mockLayerMDfwRefLower[iL].size()
		 <<" refUps "<<mockLayerMDfwRefUpper[iL].size()
		 <<" refMDs " <<mockLayerMDfwRef[iL].size()
		 <<" outLos " <<mockLayerMDfwD2cmLower[iL].size()
		 <<" outUps " <<mockLayerMDfwD2cmUpper[iL].size()
		 <<" outMDs " <<mockLayerMDfwD2cm[iL].size()
		 <<" outSDs "<< mockLayerSDfwD2cm[iL].size()
		 <<" n1Any "<<nSDs1GeV<<" n1All "<<nSDs1GeVMatch4
		 <<" n2Any "<<nSDs2GeV<<" n2All "<<nSDs2GeVMatch4
		 <<std::endl;
      }//iL
      std::cout<<"\t\tsdl5-7 "<<mockLayer5to7SDLfwD2cm.size()
	       <<" sdl7-9 "<<mockLayer7to9SDLfwD2cm.size()
	       <<std::endl;
      
      // Progress feedback to the user
      if(nEventsTotal%1 == 0) {
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)) {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
          "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }//if(nEventsTotal%20000 == 0) {


    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  std::cout<<__LINE__<<" make efficiencies "<<std::endl;

  if (drawPlots){
    std::cout<<__LINE__<<" draw and print "<<std::endl;
    for (int iL = 5; iL <= nLayers; ++iL){
      if (iL != 5 && iL != 10) continue;

    }//nLayers
  }//if drawPlots
  
  std::cout<<__LINE__<<" write to file "<<std::endl;
  TFile* outHistograms = new TFile("outHistogramsSuperD.root", "RECREATE");
  for (int iL = 1; iL <= nLayers; ++iL){
    layerMD_pt_all[iL].write(outHistograms);
    layerMD_pt_prim_all[iL].write(outHistograms);
    layerMD_pt_prim_tt[iL].write(outHistograms);
  }
  outHistograms->Write();
  outHistograms->Close();
  std::cout<<__LINE__<<" Done "<<std::endl;

  return 0;
}

