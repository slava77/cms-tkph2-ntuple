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
#include "TStopWatch.h"

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
  TVector3 p3; //makes sense mostly for seed-based
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
  double ptIn;
  double ptOut;

  bool hasMatch_byHit4of4=false;
};

struct TrackLink {
  std::vector<const SDLink*> links;
  
  double pt;
  double eta;
  double phi;
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

struct PixelHit {
  PixelHit(int i):
    r3s(pix_xsim()[i], pix_ysim()[i], pix_zsim()[i]),
    p3s(pix_pxsim()[i], pix_pysim()[i], pix_pzsim()[i]),
    ind(i),
    lay(pix_lay()[i]),
    pdgId(pix_particle()[i]),
    process(pix_process()[i]),
    bx(pix_bunchXing()[i]),
    evt(pix_event()[i]),
    isBarrel(pix_isBarrel()[i])
  {}
  TVector3 r3s;
  TVector3 p3s;
  int ind;
  int lay;
  int pdgId;
  int process;
  int bx;
  int evt;
  bool isBarrel;
  void print(const std::string& pfx){
    std::cout<<pfx<<" "<<ind<<" L"<<lay<<(isBarrel? "b ": "e ")<<evt<<":"<<bx<<" "<<pdgId<<":"<<process
	     <<" ("<<r3s.Pt()<<", "<<r3s.Eta()<<","<<r3s.Phi()<<") "
	     <<" ("<<p3s.Pt()<<", "<<p3s.Eta()<<","<<p3s.Phi()<<") "
	     <<std::endl;
  }
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
  double disc = bb*bb + (2.*rt/d + 1.);
  status = 0;
  if (disc < 0){
    status = 1;
    return r3;
  }
  double dSign = useClosest ? 1. : -1.;
  double xxP = (p/pt)*( sqrt(bb*bb + (2.*rt/d + 1.)) - bb);
  double xxM = (p/pt)*( - sqrt(bb*bb + (2.*rt/d + 1.)) - bb);
  double xx;
  if (useClosest){
    xx = std::abs(xxP) < std::abs(xxM) ? xxP : xxM;
  } else {
    xx = std::abs(xxP) < std::abs(xxM) ? xxM : xxP;
  }
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

std::pair<TVector3,TVector3> helixPropagateApprox(const TVector3& r3, const TVector3& p3, double rDest, int q, int& status, bool useClosest = true, bool verbose = false){
  double epsilon = 0.001;
  double p = p3.Mag();
  double kap = (2.99792458e-3*3.8*q/p);
  
  auto lastR3 = r3;
  auto lastT3 = p3.Unit();
  int nIts = 5;

  while (std::abs(lastR3.Perp() - rDest) > epsilon && nIts >= 0){
    auto lineEst = linePropagate(lastR3, lastT3*p, rDest, status, useClosest, verbose);
    if (status){
      if (verbose) std::cout<<" failed with status "<<status<<std::endl;
      return { lineEst, lastT3*p};
    }
    if (q==0) return {lineEst, lastT3*p};
    
    double dir = (lineEst.x() - lastR3.x())*lastT3.x() + (lineEst.y() - lastR3.y())*lastT3.y() > 0 ? 1. : -1;
    double dS = (lineEst - lastR3).Mag()*dir;
    double phi = kap*dS;
    if (std::abs(phi) > 1) {
      if (verbose) std::cout<<" return line for very large angle "<<status<<std::endl;
      return { lineEst, lastT3*p};
    }
    double alpha = 1 - sin(phi)/phi;
    double beta = (1 - cos(phi))/phi;
    
    TVector3 tau = lastT3; 
    
    TVector3 hEstR3(tau.x()*(1.-alpha) + tau.y()*beta, tau.y()*(1.-alpha) - tau.x()*beta, tau.z());
    hEstR3 *= dS;
    hEstR3 += lastR3;
    lastR3 = hEstR3;
    
    TVector3 hEstT3(tau.x()*cos(phi) + tau.y()*sin(phi), tau.y()*cos(phi) - tau.x()*sin(phi), tau.z());
    lastT3 = hEstT3;
    --nIts;
    if (verbose){
      std::cout<<"nIts "<<nIts<<" rDest "<<rDest<<" dS "<<dS<<" phi "<<phi
	       <<" r3In ("<<r3.Pt()<<", "<<r3.Eta()<<", "<<r3.Phi()<<")"
	       <<" p3In ("<<p3.Pt()<<", "<<p3.Eta()<<", "<<p3.Phi()<<")"
	       <<" r3out ("<<lastR3.Pt()<<", "<<lastR3.Eta()<<", "<<lastR3.Phi()<<")"
	       <<" p3Out ("<<lastT3.Pt()*p<<", "<<lastT3.Eta()<<", "<<lastT3.Phi()<<")"
	       <<std::endl;
    }
  }
  status = (std::abs(lastR3.Perp() - rDest) > epsilon);
  return {lastR3, lastT3*p};
  
}

int ScanChainMiniDoublets( TChain* chain, int nEvents = -1, bool drawPlots = false, double ptMinGlobal = 1.0) {

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
	if (pts < 0.8*ptMinGlobal )  continue;
	
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
      }//different sim types

    }//iL from 5 to nLayers

    {
      std::array<TProfile*, 3> avOccupancyTT_LA {
	layerMD_pt_prim_all[5].mdOthers2mm2SH_aCut.profile,
	layerMD_pt_prim_all[7].mdOthers2mm2SH_aCut.profile,
	layerMD_pt_prim_all[9].mdOthers2mm2SH_aCut.profile  };

      auto h5 = avOccupancyTT_LA[0];
      auto h7 = avOccupancyTT_LA[1];
      auto h9 = avOccupancyTT_LA[2];
      h5->SetTitle("Mini-doublet wrong combination occupancy; p_{T}^{Sim} (GeV); Wrong combination rate");
      auto cn = h5->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      gPad->SetLogx();
      h5->SetStats(0);
      h5->Draw();
      auto ax = h5->GetYaxis();
      ax->SetTitleOffset(1.5);
      h5->SetMinimum(0.01);
      h5->SetMaximum(0.5);
      h5->SetLineWidth(2);
      h5->SetLineColor(kBlack);

      h7->SetLineWidth(2);
      h7->SetLineColor(kRed);
      
      h9->SetLineWidth(2);
      h9->SetLineColor(kBlue);

      h7->Draw("same");
      h9->Draw("same");

      TLegend* leg = new TLegend(0.4, 0.8, 0.8, 0.9);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(3);
      leg->AddEntry(h5, "L5");
      leg->AddEntry(h7, "L7");
      leg->AddEntry(h9, "L9");
      leg->Draw();
      gPad->SaveAs("h_mdWrong_L5L7L9_primary.png");
      delete cv;
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

int ScanChainMockSuperDoublets( TChain* chain, int nEvents = -1, const bool drawPlots = false, const int mockMode = 0, const double sdOffset = 2.0,
				const int useSeeds = 0) {
  //mockMode:
  //0 for helix to ref and then straight line;
  //1 for helix to all ref layers;
  //3: as 1, for first hit offset to the rechit position

  //useSeeds:
  //0: not used
  //1: convert to tangent SuperDoublet at outer hit
  
  std::cout<<"Running in mockMode "<<mockMode<<std::endl;
  std::cout<<"Running with SD distance "<<sdOffset<<std::endl;
  std::cout<<"Running with useSeeds "<<useSeeds<<std::endl;
  
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

  constexpr float pixelZpitch = 0.15;
  //p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
  //(1 sigma based on 95.45% = 2sigma at 2 GeV)
  std::array<float, nLayers+1> miniMulsPtScale {0, 0, 0, 0, 0,
      0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};

  //mean of the horizontal layer position in y; treat this as R below
  std::array<float, nLayers+1> miniRminMean {0, 0, 0, 0, 0, //may want to fill these
      21.8, 34.6, 49.6, 67.4, 87.6, 106.8};

  //0-5
  TH2F* h2_sdl0to5_dBeta_betaIn_NM1dBeta_all = new TH2F("h2_sdl0to5_dBeta_betaIn_NM1dBeta_all", "h2_sdl0to5_dBeta_betaIn_NM1dBeta_all", 400, -1, 1, 400, -0.5, 0.5);
  
  TH1F* h_sdl0to5_dBeta_NM1dBeta_all = new TH1F("h_sdl0to5_dBeta_NM1dBeta_all", "h_sdl0to5_dBeta_NM1dBeta_all", 400, -0.5, 0.5);
  TH1F* h_sdl0to5_dBeta_NM1dBeta_pass = new TH1F("h_sdl0to5_dBeta_NM1dBeta_pass", "h_sdl0to5_dBeta_NM1dBeta_pass", 400, -0.5, 0.5);

  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_all = new TH1F("h_sdl0to5_dBeta_zoom_NM1dBeta_all", "h_sdl0to5_dBeta_zoom_NM1dBeta_all", 400, -0.15, 0.15);
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_pass = new TH1F("h_sdl0to5_dBeta_zoom_NM1dBeta_pass", "h_sdl0to5_dBeta_zoom_NM1dBeta_pass", 400, -0.15, 0.15);

  std::string hns = "h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_all";
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_pass";
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_all";
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_pass";
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_all";
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_pass";
  TH1F* h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  //0-7
  TH2F* h2_sdl0to7_dBeta_betaIn_NM1dBeta_all = new TH2F("h2_sdl0to7_dBeta_betaIn_NM1dBeta_all", "h2_sdl0to7_dBeta_betaIn_NM1dBeta_all", 400, -1, 1, 400, -0.5, 0.5);
  
  TH1F* h_sdl0to7_dBeta_NM1dBeta_all = new TH1F("h_sdl0to7_dBeta_NM1dBeta_all", "h_sdl0to7_dBeta_NM1dBeta_all", 400, -0.5, 0.5);
  TH1F* h_sdl0to7_dBeta_NM1dBeta_pass = new TH1F("h_sdl0to7_dBeta_NM1dBeta_pass", "h_sdl0to7_dBeta_NM1dBeta_pass", 400, -0.5, 0.5);

  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_all = new TH1F("h_sdl0to7_dBeta_zoom_NM1dBeta_all", "h_sdl0to7_dBeta_zoom_NM1dBeta_all", 400, -0.15, 0.15);
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_pass = new TH1F("h_sdl0to7_dBeta_zoom_NM1dBeta_pass", "h_sdl0to7_dBeta_zoom_NM1dBeta_pass", 400, -0.15, 0.15);

  hns = "h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_all";
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass";
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_all";
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass";
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all";
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass";
  TH1F* h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  //5-7  
  TH2F* h2_sdl5to7_dBeta_betaIn_NM1dBeta_all = new TH2F("h2_sdl5to7_dBeta_betaIn_NM1dBeta_all", "h2_sdl5to7_dBeta_betaIn_NM1dBeta_all", 400, -1, 1, 400, -0.5, 0.5);
  TH2F* h2_sdl7to9_dBeta_betaIn_NM1dBeta_all = new TH2F("h2_sdl7to9_dBeta_betaIn_NM1dBeta_all", "h2_sdl7to9_dBeta_betaIn_NM1dBeta_all", 400, -1, 1, 400, -0.5, 0.5);
  
  TH1F* h_sdl5to7_dBeta_NM1dBeta_all = new TH1F("h_sdl5to7_dBeta_NM1dBeta_all", "h_sdl5to7_dBeta_NM1dBeta_all", 400, -0.5, 0.5);
  TH1F* h_sdl5to7_dBeta_NM1dBeta_pass = new TH1F("h_sdl5to7_dBeta_NM1dBeta_pass", "h_sdl5to7_dBeta_NM1dBeta_pass", 400, -0.5, 0.5);

  TH1F* h_sdl7to9_dBeta_NM1dBeta_all = new TH1F("h_sdl7to9_dBeta_NM1dBeta_all", "h_sdl7to9_dBeta_NM1dBeta_all", 400, -0.5, 0.5);
  TH1F* h_sdl7to9_dBeta_NM1dBeta_pass = new TH1F("h_sdl7to9_dBeta_NM1dBeta_pass", "h_sdl7to9_dBeta_NM1dBeta_pass", 400, -0.5, 0.5);

  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_all = new TH1F("h_sdl5to7_dBeta_zoom_NM1dBeta_all", "h_sdl5to7_dBeta_zoom_NM1dBeta_all", 400, -0.15, 0.15);
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_pass = new TH1F("h_sdl5to7_dBeta_zoom_NM1dBeta_pass", "h_sdl5to7_dBeta_zoom_NM1dBeta_pass", 400, -0.15, 0.15);

  hns = "h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_all";
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass";
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_all";
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass";
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all";
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass";
  TH1F* h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_all = new TH1F("h_sdl7to9_dBeta_zoom_NM1dBeta_all", "h_sdl7to9_dBeta_zoom_NM1dBeta_all", 400, -0.15, 0.15);
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_pass = new TH1F("h_sdl7to9_dBeta_zoom_NM1dBeta_pass", "h_sdl7to9_dBeta_zoom_NM1dBeta_pass", 400, -0.15, 0.15);

  hns = "h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_all";
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_pass";
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_all";
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_pass";
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  hns = "h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_all";
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_all = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);
  hns = "h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_pass";
  TH1F* h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_pass = new TH1F(hns.c_str(), hns.c_str(), 400, -0.15, 0.15);

  
  enum SDLayers {SDL_L0to5=0, SDL_L0to7, SDL_L5to7, SDL_L7to9, SDL_L5to9, SDL_LMAX};
  std::array<TH1F*, SDL_LMAX> ha_denSDL_pt;
  std::array<TH1F*, SDL_LMAX> ha_num8MH_pt;
  std::array<TH1F*, SDL_LMAX> ha_num4MD_pt;
  std::array<TH1F*, SDL_LMAX> ha_num2SD_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_num2SD_3of4_any_pt;

  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w01_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w012_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0123_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w01234_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w012345_4of4_pt;//see SDLSelectFlags definitions

  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_8MH;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom2_NM1dBeta_8MH;

  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_pass;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_zoom_pass;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaOut_zoom_pass;

  std::array<TH1F*, SDL_LMAX> ha_numSDL_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_numSDL_3of4_any_pt;

  std::array<TH1F*, SDL_LMAX> ha_SDLreco_all_pt;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_no4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_all_eta;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_4of4_eta;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_no4of4_eta;

  

  
  std::array<std::array<int, 2>, SDL_LMAX> layersSDL {{ {0, 5}, {0, 7}, {5, 7}, {7, 9}, {5, 9} }};
  for (int i = 0; i< SDL_LMAX; ++i){
    auto iMin = layersSDL[i][0];
    auto iMax = layersSDL[i][1];
    std::string hn = Form("h_denSDL_%dto%d_pt", iMin, iMax);
    ha_denSDL_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());

    hn = Form("h_num8MH_%dto%d_pt", iMin, iMax);
    ha_num8MH_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    

    hn = Form("h_num4MD_%dto%d_pt", iMin, iMax);
    ha_num4MD_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    

    hn = Form("h_num2SD_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_num2SD_3of4_any_%dto%d_pt", iMin, iMax);
    ha_num2SD_3of4_any_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());

    hn = Form("h_num2SD_w0_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w0_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_num2SD_w01_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w01_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_num2SD_w012_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w012_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_num2SD_w0123_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w0123_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_num2SD_w01234_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w01234_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_num2SD_w012345_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w012345_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());

    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_8MH[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    hn = Form("h_SDL_dBeta_zoom2_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom2_NM1dBeta_8MH[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.02, 0.02);

    hn = Form("h2_SDL_dBeta_betaIn_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -1, 1, 400, -0.5, 0.5);

    hn = Form("h2_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);

    hn = Form("h2_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);

    hn = Form("h2_SDL_dBeta_betaIn_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_pass[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -1, 1, 400, -0.5, 0.5);

    hn = Form("h2_SDL_dBeta_betaIn_zoom_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_zoom_pass[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);

    hn = Form("h2_SDL_dBeta_betaOut_zoom_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaOut_zoom_pass[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);

    hn = Form("h_numSDL_4of4_%dto%d_pt", iMin, iMax);
    ha_numSDL_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    hn = Form("h_numSDL_3of4_any_%dto%d_pt", iMin, iMax);
    ha_numSDL_3of4_any_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());

    hn = Form("h_SDLreco_all_%dto%d_pt", iMin, iMax);
    ha_SDLreco_all_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    hn = Form("h_SDLreco_4of4_%dto%d_pt", iMin, iMax);
    ha_SDLreco_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    hn = Form("h_SDLreco_no4of4_%dto%d_pt", iMin, iMax);
    ha_SDLreco_no4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());

    hn = Form("h_SDLreco_all_%dto%d_eta", iMin, iMax);
    ha_SDLreco_all_eta[i] = new TH1F(hn.c_str(), hn.c_str(), 40, -2.0, 2.0);
    hn = Form("h_SDLreco_4of4_%dto%d_eta", iMin, iMax);
    ha_SDLreco_4of4_eta[i] = new TH1F(hn.c_str(), hn.c_str(), 40, -2.0, 2.0);
    hn = Form("h_SDLreco_no4of4_%dto%d_eta", iMin, iMax);
    ha_SDLreco_no4of4_eta[i] = new TH1F(hn.c_str(), hn.c_str(), 40, -2.0, 2.0);
  }

  enum TimerTypes {T_timeLayout=0, T_timeReco, T_timeValidation, T_N};
  std::array<TStopwatch, T_N> timerA {};
  std::array<std::string, T_N> timerNameA {"timeLayout", "timeReco", "timeValidation"};

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
      //offset layer minidoublet (Ncm hardcode for now)
      std::array<std::vector<std::pair<int, TVector3> >, nLayers+1 > mockLayerMDfwDNcmLower;
      std::array<std::vector<std::pair<int, TVector3> >, nLayers+1 > mockLayerMDfwDNcmUpper;
      for (int iL = 1; iL <=nLayers; ++iL){
	mockLayerMDfwRefLower[iL].reserve(maxHitsInLayer);
	mockLayerMDfwRefUpper[iL].reserve(maxHitsInLayer);
	mockLayerMDfwDNcmLower[iL].reserve(maxHitsInLayer);
	mockLayerMDfwDNcmUpper[iL].reserve(maxHitsInLayer);
      }

      std::array<int, nLayers+1> nHitsLayer1GeV {};
      std::array<int, nLayers+1> nHitsLayer2GeV {};

      timerA[T_timeLayout].Start(kFALSE);
      std::cout<<"Load hits "<<std::endl;
      std::array<std::set<int>, nLayers+1> simIdxDeltaInLayer;
      std::array<std::set<int>, nLayers+1> simIdxPositronLowPInLayer;
      std::array<std::set<int>, nLayers+1> simIdxInLayer;
      auto nPix = pix_isBarrel().size();
      for (auto ipix = 0U; ipix < nPix; ++ipix){
	if (pix_isBarrel()[ipix] == false) continue;
	int lay = pix_lay()[ipix];
	if (lay < 5 ) continue;

	int iid = pix_detId()[ipix];
	//Too restrictive for reverse TP matching??//	if ((iid & 0x4)!= 4) continue; 

	TVector3 r3Rec(pix_x()[ipix], pix_y()[ipix], pix_z()[ipix]);

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
	} else if (iParticle == 11 && ps < 0.1) {
	  if (simIdxPositronLowPInLayer[lay].find(iSimIdx) != simIdxPositronLowPInLayer[lay].end()) continue; //only one hit per layer per track
	  else {
	    simIdxPositronLowPInLayer[lay].insert(iSimIdx);
	  }	  
	} else {
	  if (simIdxInLayer[lay].find(iSimIdx) != simIdxInLayer[lay].end()) continue; //only one hit per layer per track
	  else {
	    simIdxInLayer[lay].insert(iSimIdx);
	  }
	}

	if (pts >1) nHitsLayer1GeV[lay]++;
	if (pts >2) nHitsLayer2GeV[lay]++;
	
	const double mdOffset = miniDeltaBarrel[lay];
	const double rNominal = miniRminMean[lay]+1;//this is going to be in between the sector radii

	const double rRefLower = rNominal;
	const double rRefUpper = rNominal+mdOffset;
	const double rSDfwLower = rRefLower + sdOffset;
	const double rSDfwUpper = rRefUpper + sdOffset;

	int q = 0;
	if (iParticle == -11 || iParticle == -13 || iParticle == 211 || iParticle == 321 || iParticle == 2212) q = 1;
	else if (iParticle == 11 || iParticle == 13 || iParticle == -211 || iParticle == -321 || iParticle == -2212) q = -1;
	int pstat = 0;
	TVector3 r3RefLower;
	TVector3 p3RefLower;
	if (q == 0){
	  r3RefLower = linePropagate(r3Sim, p3Sim, rRefLower, pstat);
	  p3RefLower = p3Sim;
	} else {
	  auto resProp = helixPropagateApprox(r3Sim, p3Sim, rRefLower, q, pstat);
	  r3RefLower = resProp.first;
	  p3RefLower = resProp.second;
	  
	}
	
	auto propagateMH = [&](float rDest){
	  if (mockMode == 0) return linePropagate(r3RefLower, p3RefLower, rDest, pstat);
	  else if (mockMode == 1 || mockMode == 3) return helixPropagateApprox(r3RefLower, p3RefLower, rDest, q, pstat).first;
	  else {pstat = 99; return TVector3();}
	};

	auto r3RefLowerMock = r3RefLower;
	if (mockMode == 3){
	  r3RefLowerMock += (r3Rec - r3Sim);
	  if ((iid & 0x4)!= 4 && (lay >=5 && lay <=7)){	    
	    r3RefLowerMock.SetZ(r3RefLowerMock.z() - r3Rec.z() + r3Sim.z());// undo the rec-sim shift
	    //there was no simhit on the pixel layer; make up z by roundoff of sim
	    float binnedZshift = std::round(r3Sim.z()/pixelZpitch)*pixelZpitch - r3Sim.z();
	    r3RefLowerMock.SetZ(r3RefLowerMock.z() + binnedZshift);
	  }
	}
        if (pstat == 0) mockLayerMDfwRefLower[lay].push_back(std::make_pair(ipix,r3RefLowerMock));
	//decide to go full scale helix based on config:
	auto r3RefUpper = propagateMH(rRefUpper);
	if (pstat == 0) mockLayerMDfwRefUpper[lay].push_back(std::make_pair(ipix, r3RefUpper));
	auto r3SDfwLower = propagateMH(rSDfwLower);
	if (pstat == 0) mockLayerMDfwDNcmLower[lay].push_back(std::make_pair(ipix, r3SDfwLower));
	auto r3SDfwUpper = propagateMH(rSDfwUpper);
	if (pstat == 0) mockLayerMDfwDNcmUpper[lay].push_back(std::make_pair(ipix, r3SDfwUpper));


      }//nPix: filling mock MDs
      timerA[T_timeLayout].Stop();


      timerA[T_timeReco].Start(kFALSE);
      constexpr int maxMDexpected = 100;
      std::array<std::vector<MiniDoublet>, nLayers+1> mockLayerMDfwRef;
      std::array<std::vector<MiniDoublet>, nLayers+1> mockLayerMDfwDNcm;

      std::array<std::vector<SuperDoublet>, nLayers+1> mockLayerSDfwDNcm;

      std::vector<SDLink> mockLayer0to5SDLfwDNcm;
      std::vector<SDLink> mockLayer0to7SDLfwDNcm;
      std::vector<SDLink> mockLayer5to7SDLfwDNcm;
      std::vector<SDLink> mockLayer7to9SDLfwDNcm;
      
      if (useSeeds == 1){
	std::cout<<"Convert seeds to SuperDoublets"<<std::endl;
	auto const nSeeds = see_lh_px().size();
	for (auto iSeed = 0U; iSeed < nSeeds; ++iSeed){
	  TVector3 p3LH(see_lh_px()[iSeed], see_lh_py()[iSeed], see_lh_pz()[iSeed]);
	  float ptLH = p3LH.Pt();
	  float etaLH = p3LH.Eta();
	  if (ptLH < 1.0) continue; //1 GeV cut
	  if (std::abs(etaLH) > 1.6) continue;
	  TVector3 r3LH(see_lh_x()[iSeed], see_lh_y()[iSeed], see_lh_z()[iSeed]);

	  TVector3 p3PCA(see_pca_px()[iSeed], see_pca_py()[iSeed], see_pca_pz()[iSeed]);
	  TVector3 r3PCA(see_pca_px()[iSeed], see_pca_py()[iSeed], see_pca_pz()[iSeed]);
	  SuperDoublet seedSD;
	  seedSD.mdRef.r3 = r3PCA;
	  seedSD.mdRef.alpha = r3PCA.DeltaPhi(p3PCA);
	  seedSD.mdOut.r3 = r3LH;
	  seedSD.mdOut.alpha = r3LH.DeltaPhi(p3LH);	  
	  seedSD.iRef = iSeed;
	  seedSD.iOut = iSeed;
	  seedSD.r3 = r3LH;
	  seedSD.p3 = p3LH;
	  seedSD.alpha = r3LH.DeltaPhi(p3LH);
	  mockLayerSDfwDNcm[0].emplace_back(seedSD);
	}
	std::cout<<"Loaded nSeeds (pt>1 and |eta|<1.6) "<<mockLayerSDfwDNcm[0].size()<<std::endl;
      }

      std::cout<<"Combine MDs "<<std::endl;
      int nHitsTried = 0;
      for (int iL = minLayer; iL <= nLayers; ++iL){
	auto const&  hitsRefLower = mockLayerMDfwRefLower[iL];
	auto const&  hitsRefUpper = mockLayerMDfwRefUpper[iL];
	auto const&  hitsOutLower = mockLayerMDfwDNcmLower[iL];
	auto const&  hitsOutUpper = mockLayerMDfwDNcmUpper[iL];
	
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

	auto& mockMDfwDNcm = mockLayerMDfwDNcm[iL];
	mdCombine(hitsOutLower, hitsOutUpper, mockMDfwDNcm);

	//now make super-doublets
	auto& mockSDfwDNcm = mockLayerSDfwDNcm[iL];

	int iRef = -1;
	for (auto mdRef : mockMDfwRef){
	  iRef++;
	  float rtRef = mdRef.r3.Pt();
	  float zRef = mdRef.r3.z();
	  //
	  int iOut = -1;
	  for (auto mdOut : mockMDfwDNcm){
	    iOut++;
	    float rtOut = mdOut.r3.Pt();
	    float zOut = mdOut.r3.z();

	    //apply some loose Z compatibility
	    float zGeom = iL >= 5 && iL <= 6 ? 0.3 : 10;//twice the macro-pixel or strip size
	                                                //assume that the mock layer is the n+1 layer
	    float zLo = rtOut/rtRef*(zRef - 15.) - zGeom; //15 for the luminous ; 10 for module size
	    float zHi = rtOut/rtRef*(zRef + 15.) + zGeom;
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
	    mockSDfwDNcm.push_back(sd);
	  }
	}
      }//iL

      //try to link segments: do 5-7 and 7-9; or 0-X for seeds (special cases for L=0 where needed)
      auto sdLink = [&] (int lIn, int lOut, decltype(mockLayer5to7SDLfwDNcm)& sdlV){
	auto const& sdInV  = mockLayerSDfwDNcm[lIn];
	auto const& sdOutV = mockLayerSDfwDNcm[lOut];

	int nAll = 0;
	int nDeltaZ = 0;
	int nDeltaZPointed = 0;
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
	  float dSDIn = sdIn.mdOut.r3.Pt() - sdIn.mdRef.r3.Pt();
	  float dzSDIn = sdIn.mdOut.r3.z() - sdIn.mdRef.r3.z();
	  
	  int iOut = -1;
	  for ( auto sdOut : sdOutV ) {
	    iOut++;
	    //
	    nAll++;
	    
	    float rtOut = sdOut.r3.Pt();
	    float zOut = sdOut.r3.z();
	    //apply some loose Z compatibility
	    //FIXME: refine using inner layer directions (can prune later)
	    float zGeom = lIn >= 0 && lIn <= 7 && lOut >= 5 && lOut <= 7 ? 0.3 : 10;//twice the macro-pixel or strip size
	    float zLo = rtOut/rtIn*(zIn - 15.) - zGeom; //15 for the luminous ; zGeom for z geom unit size
	    float zHi = rtOut/rtIn*(zIn + 15.) + zGeom;
	    if (zOut < zLo || zOut > zHi) continue;
	    nDeltaZ++;

	    if (lIn == 0){
	      float etaErr = see_pca_etaErr()[sdIn.iRef];
	      float eta = see_lh_eta()[sdIn.iRef];
	      float dzErr = (rtOut - rtIn)*etaErr*sinh(eta);
	      dzErr *= dzErr;
	      dzErr += 0.03*0.03; // pixel size x2. ... random for now
	      dzErr *= 9; //3 sigma
	      dzErr += zGeom*zGeom;
	      dzErr = sqrt(dzErr);
	      float dzDrIn = std::abs(sdIn.p3.Z())/sdIn.p3.Pt();
	      float zLo = zIn + (dzDrIn - dzErr/dSDIn)*(rtOut - rtIn) - zGeom;
	      float zHi = zIn + (dzDrIn + dzErr/dSDIn)*(rtOut - rtIn) + zGeom;
	      if (zOut < zLo || zOut > zHi) continue;
	      //FIXME: ADD A PLOT HERE
	    }
	    else if (lIn>=5 && lIn <=6){//can point to the z pos in lOut
	      float zLo = zIn + (dzSDIn - zGeom*sqrt(2.))/dSDIn*(rtOut - rtIn) - zGeom;
	      float zHi = zIn + (dzSDIn + zGeom*sqrt(2.))/dSDIn*(rtOut - rtIn) + zGeom;
	      if (zOut < zLo || zOut > zHi) continue;
	    }
	    nDeltaZPointed++;
	    
	    auto midR3 = 0.5*(sdIn.r3 + sdOut.r3);
	    double dPhi = midR3.DeltaPhi(sdOut.r3 - sdIn.r3);
	    double rt = 0.5*(sdIn.r3.Pt() + sdOut.r3.Pt());

	    float ptSLo = 1.0;
	    if (lIn == 0){
	      //try to use seed pt: the lower bound is good
	      ptSLo = sdIn.p3.Pt();
	      float ptErr = see_pca_ptErr()[sdIn.iRef];
	      ptSLo = std::max(1.0f, ptSLo - 10.0f*ptErr);
	      ptSLo = std::min(10.0f, ptSLo); //don't let this run away either
	    }
	    const float ptCut = ptSLo;
	    const float sdlSlope = rt/175.67/ptCut;
	    const float sdlThetaMulsF = 0.015*sqrt(0.2);
	    const float sdlMuls = sdlThetaMulsF*3./ptCut*4;//will need a better guess than x4?
	    const float sdlPVoff = 0.1/rt;
	    const float sdlCut = sdlSlope + sqrt(sdlMuls*sdlMuls + sdlPVoff*sdlPVoff);
	    
	    if (std::abs(dPhi) > sdlCut ) continue;
	    nSlope++;
	    
	    double betaIn;
	    double betaOut;
	    if (mockMode == 0){
	      betaIn = sdIn.alpha - sdIn.r3.DeltaPhi(sdOut.r3 - sdIn.r3);
	      betaOut = - sdOut.alpha + sdOut.r3.DeltaPhi(sdOut.r3 - sdIn.r3); //to match sign for correct match	      
	    }
	    else if (mockMode == 1 || mockMode == 3){
	      //need a symmetric choice of end-points to achieve partial cancelation
	      betaIn  = sdIn.alpha - sdIn.r3.DeltaPhi(sdOut.mdOut.r3 - sdIn.r3);
	      betaOut = -sdOut.mdOut.r3.DeltaPhi(sdOut.mdOut.r3 - sdOut.mdRef.r3) + sdOut.mdOut.r3.DeltaPhi(sdOut.mdOut.r3 - sdIn.r3);	      
	    }
	    else{
	      betaIn = -99;
	      betaOut = 999;
	    }
	    
	    //loose angle compatibility
	    float dAlpha_Bfield = (rtOut - rtIn)/175.67/ptCut;
	    float dSDOut = sdOut.mdOut.r3.Pt() - sdOut.mdRef.r3.Pt();	
	    float dAlpha_res = 0.02/std::min(dSDIn, dSDOut);//2-strip difference; use the smallest SD separation
	    float dAlpha_compat = dAlpha_Bfield + dAlpha_res;
	    if (std::abs(sdIn.alpha- dPhi) > dAlpha_compat) continue;
	    nInAlphaCompat++;
	    if (std::abs(sdOut.alpha- dPhi) > dAlpha_compat) continue;
	    nOutAlphaCompat++;
	    
	    //now the actual segment linking magic
	    float betaAv = 0.5*(betaIn + betaOut);
	    float dr = (sdOut.r3 - sdIn.r3).Perp();
	    //pt*175.67/2. = R
	    //R*sin(betaAv) = pt*175.67/2*sin(betaAv) = dr/2 => pt = dr/175.67/sin(betaAv);
	    float pt_beta = dr/175.67/sin(betaAv);
	    if (lIn == 0) pt_beta = sdIn.p3.Pt();
	    float pt_betaIn = dr/175.67/sin(betaIn);
	    if (lIn == 0) pt_betaIn = pt_beta;
	    float pt_betaOut = dr/175.67/sin(betaOut);
	    float dBetaRes = dAlpha_res;
	    float dBetaMuls = sdlThetaMulsF*3./std::min(pt_beta, 7.0f);//need to confirm the range-out value of 7 GeV
	    float dBetaCut = sqrt(dBetaRes*dBetaRes*2.0 + dBetaMuls*dBetaMuls);
	    float dBeta = betaIn - betaOut;

	    auto ptIn = std::abs(pt_betaIn);

	    //FIXME: need to reduce copy-paste
	    if (lIn == 0 && lOut == 5){
	      h_sdl0to5_dBeta_NM1dBeta_all->Fill(dBeta);
	      h2_sdl0to5_dBeta_betaIn_NM1dBeta_all->Fill(betaIn, dBeta);
	      h_sdl0to5_dBeta_zoom_NM1dBeta_all->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_all->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_all->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_all->Fill(dBeta);
	      }
	    } else if (lIn == 0 && lOut == 7){
	      h_sdl0to7_dBeta_NM1dBeta_all->Fill(dBeta);
	      h2_sdl0to7_dBeta_betaIn_NM1dBeta_all->Fill(betaIn, dBeta);
	      h_sdl0to7_dBeta_zoom_NM1dBeta_all->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_all->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_all->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all->Fill(dBeta);
	      }
	    } else if (lIn == 5 && lOut == 7){
	      h_sdl5to7_dBeta_NM1dBeta_all->Fill(dBeta);
	      h2_sdl5to7_dBeta_betaIn_NM1dBeta_all->Fill(betaIn, dBeta);
	      h_sdl5to7_dBeta_zoom_NM1dBeta_all->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_all->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_all->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all->Fill(dBeta);
	      }
	    } else if (lIn == 7 && lOut == 9){
	      h_sdl7to9_dBeta_NM1dBeta_all->Fill(dBeta);
	      h2_sdl7to9_dBeta_betaIn_NM1dBeta_all->Fill(betaIn, dBeta);
	      h_sdl7to9_dBeta_zoom_NM1dBeta_all->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_all->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_all->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_all->Fill(dBeta);
	      }
	    }

	    if (std::abs(dBeta) > dBetaCut) continue;
	    if (lIn == 0 && lOut == 5){
	      h_sdl0to5_dBeta_NM1dBeta_pass->Fill(dBeta);
	      h_sdl0to5_dBeta_zoom_NM1dBeta_pass->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_pass->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_pass->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_pass->Fill(dBeta);
	      }
	    } else if (lIn == 0 && lOut == 7){
	      h_sdl0to7_dBeta_NM1dBeta_pass->Fill(dBeta);
	      h_sdl0to7_dBeta_zoom_NM1dBeta_pass->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass->Fill(dBeta);
	      }
	    } else if (lIn == 5 && lOut == 7){
	      h_sdl5to7_dBeta_NM1dBeta_pass->Fill(dBeta);
	      h_sdl5to7_dBeta_zoom_NM1dBeta_pass->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass->Fill(dBeta);
	      }
	    } else if (lIn == 7 && lOut == 9){
	      h_sdl7to9_dBeta_NM1dBeta_pass->Fill(dBeta);
	      h_sdl7to9_dBeta_zoom_NM1dBeta_pass->Fill(dBeta);
	      if (ptIn < 2){
		h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_pass->Fill(dBeta);
	      }
	      if (ptIn > 3 && ptIn < 5 ){
		h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_pass->Fill(dBeta);
	      }
	      if (ptIn > 7) {
		h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_pass->Fill(dBeta);
	      }
	    }

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
	    sdl.ptIn = pt_betaIn;
	    sdl.ptOut = pt_betaOut;
	    
	    sdlV.push_back(sdl);
	  }//sdOutV
	}//sdInV
	
	std::cout<<"SD links stat "<<lIn<<"-"<<lOut
	         <<" nAll "<<nAll<<" nDZ "<<nDeltaZ<<" nDZPnt "<<nDeltaZPointed<<" nAlp "<<nSlope
         	 <<" nInAlpC "<< nInAlphaCompat<<" nOutAlpC "<<nOutAlphaCompat
	         <<" ndBeta "<<ndBeta
	         <<" final "<<sdlV.size()
	         <<std::endl;
      };//auto sdLink

      sdLink(0, 5, mockLayer0to5SDLfwDNcm);
      sdLink(0, 7, mockLayer0to7SDLfwDNcm);
      sdLink(5, 7, mockLayer5to7SDLfwDNcm);
      sdLink(7, 9, mockLayer7to9SDLfwDNcm);

      std::array<decltype(mockLayer5to7SDLfwDNcm)*, SDL_LMAX> mockLayerSDLsDNcm {};
      mockLayerSDLsDNcm[SDL_L0to5] = &mockLayer0to5SDLfwDNcm;
      mockLayerSDLsDNcm[SDL_L0to7] = &mockLayer0to7SDLfwDNcm;
      mockLayerSDLsDNcm[SDL_L5to7] = &mockLayer5to7SDLfwDNcm;
      mockLayerSDLsDNcm[SDL_L7to9] = &mockLayer7to9SDLfwDNcm;

      timerA[T_timeReco].Stop();


      timerA[T_timeValidation].Start(kFALSE);
      int nSim = sim_nPixel().size();
      for (int iSim = 0; iSim < nSim; ++iSim){
	TVector3 p3(sim_px()[iSim], sim_py()[iSim], sim_pz()[iSim]);
	bool debug = false;
	auto tpPt = p3.Pt();
	auto tpEta = p3.Eta();
	auto tpPhi = p3.Phi();

	
	//	if (tpPt > 10 && std::abs(p3.Eta())< 1) debug = true;
	if (debug) std::cout<<"TP: "<<p3.Pt()<<" "<<p3.Eta()<<" "<<p3.Phi();
	
	std::map<int, int> nHitsMap;
	std::array<std::vector<SDLink>, SDL_LMAX > matchingSDLs_byHit3of4 {};
	std::array<std::vector<SDLink>, SDL_LMAX > matchingSDLs_byHit4of4 {};
	std::array<std::vector<int>, nLayers+1> simHits {};
	
	int nPix = sim_nPixel()[iSim];
	for (int iPix = 0; iPix< nPix; ++iPix){
	  
	  int iipix = sim_pixelIdx()[iSim][iPix];
	  PixelHit pixH(iipix);
	  int lay = pixH.lay;

	  // int iProcess = pix_process()[iPix];
	  // int ibx = pix_bunchXing()[iPix];	  
	  // bool isPrimaryAny = (iProcess == 2 && ibx == 0);
	  
	  if (pixH.isBarrel){
	    if (lay >= minLayer){
	      if (debug) std::cout<<" "<<lay<<" "<<iipix;
	    }
	    if (pixH.p3s.Pt()>0.8*tpPt){
	      nHitsMap[lay]++;
	      simHits[lay].push_back(iipix);
	    }
	  }
	}
	if (debug) std::cout<<std::endl;

	for (int iSDLL = 0; iSDLL< SDL_LMAX; ++iSDLL){
	  int lIn = layersSDL[iSDLL][0];
	  int lOut = layersSDL[iSDLL][1];
	  
	  bool hasMHRefInL = false;
	  bool hasMHRefInU = false;
	  bool hasMHOutInL = false;
	  bool hasMHOutInU = false;
	  bool hasMHRefOutL = false;
	  bool hasMHRefOutU = false;
	  bool hasMHOutOutL = false;
	  bool hasMHOutOutU = false;
	  bool has8MHs = false;

	  bool hasMDRefIn = false;
	  bool hasMDOutIn = false;
	  bool hasMDRefOut = false;
	  bool hasMDOutOut = false;
	  bool has4MDs = false;
	  
	  bool hasSDIn_3of4 = false;
	  bool hasSDIn_4of4 = false;
	  bool hasSDOut_3of4 = false;
	  bool hasSDOut_4of4 = false;
	  
	  if (nHitsMap[lIn] > 0 && nHitsMap[lOut] > 0){
	    ha_denSDL_pt[iSDLL]->Fill(tpPt);
	    bool debugHitLevel = false;
	    if (debug){
	      std::cout<<"\tTP is good for denSDL in L"<<lIn<<"-L"<<lOut<<std::endl;
	      if (debugHitLevel){
		for (auto i : simHits[lIn]){ auto ph = PixelHit(i); ph.print("\t");}
		for (auto i : simHits[lOut]){ auto ph = PixelHit(i); ph.print("\t");}
	      }
	    }
	    //match the 8 layers of hits
	    auto matchMH = [&](decltype(mockLayerMDfwRefLower[lIn])const& mhs, int lay){
	      auto const& shs = simHits[lay];
	      for (auto const& mh : mhs){
		if (std::find(shs.begin(), shs.end(), mh.first) != shs.end()) return true;
	      }
	      return false;
	    };
	    hasMHRefInL = matchMH(mockLayerMDfwRefLower[lIn], lIn);
	    hasMHRefInU = matchMH(mockLayerMDfwRefUpper[lIn], lIn);
	    hasMHOutInL = matchMH(mockLayerMDfwDNcmLower[lIn], lIn);
	    hasMHOutInU = matchMH(mockLayerMDfwDNcmUpper[lIn], lIn);

	    hasMHRefOutL = matchMH(mockLayerMDfwRefLower[lOut], lOut);
	    hasMHRefOutU = matchMH(mockLayerMDfwRefUpper[lOut], lOut);
	    hasMHOutOutL = matchMH(mockLayerMDfwDNcmLower[lOut], lOut);
	    hasMHOutOutU = matchMH(mockLayerMDfwDNcmUpper[lOut], lOut);

	    has8MHs = hasMHRefInL & hasMHRefInU & hasMHOutInL & hasMHOutInU
	      & hasMHRefOutL & hasMHRefOutU & hasMHOutOutL & hasMHOutOutU;
	    if (has8MHs){
	      ha_num8MH_pt[iSDLL]->Fill(tpPt);
	    }
	    
	    //match the 4 layer mini-doublets
	    auto matchMD = [&](decltype(mockLayerMDfwRef[lIn]) const& mds, int lay){
	      auto const& shs = simHits[lay];
	      for (auto const& md : mds){
		bool hasL =  std::find(shs.begin(), shs.end(), md.pixL) != shs.end();
		bool hasU =  std::find(shs.begin(), shs.end(), md.pixU) != shs.end();
		if (hasL & hasU){
		  return true;
		}
	      }
	      return false;
	    };
	    hasMDRefIn = matchMD(mockLayerMDfwRef[lIn], lIn);
	    hasMDOutIn = matchMD(mockLayerMDfwDNcm[lIn], lIn);

	    hasMDRefOut = matchMD(mockLayerMDfwRef[lOut], lOut);
	    hasMDOutOut = matchMD(mockLayerMDfwDNcm[lOut], lOut);

	    has4MDs = hasMDRefIn & hasMDOutIn & hasMDRefOut & hasMDOutOut;
	    if (has4MDs ){
	      ha_num4MD_pt[iSDLL]->Fill(tpPt);
	    }

	    std::vector<SuperDoublet> vSDIn_4of4; vSDIn_4of4.reserve(2);
	    std::vector<SuperDoublet> vSDOut_4of4; vSDOut_4of4.reserve(2);
	    
	    //match inner and outer layer super-doublets
	    for (auto sd : mockLayerSDfwDNcm[lIn]){
	      auto const& shIn = simHits[lIn];
	      
	      bool hasIRL = std::find(shIn.begin(), shIn.end(), sd.mdRef.pixL) != shIn.end();
	      if (debug && debugHitLevel && hasIRL){
		std::cout<<"SDI: found matching hit for "<<lIn<<" IRL at "<<sd.mdRef.pixL<<std::endl;
	      }
	      bool hasIRU = std::find(shIn.begin(), shIn.end(), sd.mdRef.pixU) != shIn.end();
	      if (debug && debugHitLevel && hasIRU){
		std::cout<<"SDI: found matching hit for "<<lIn<<" IRU at "<<sd.mdRef.pixU<<std::endl;
	      }
	      bool hasIOL = std::find(shIn.begin(), shIn.end(), sd.mdOut.pixL) != shIn.end();
	      if (debug && debugHitLevel && hasIOL){
		std::cout<<"SDI: found matching hit for "<<lIn<<" IOL at "<<sd.mdOut.pixL<<std::endl;
	      }
	      bool hasIOU = std::find(shIn.begin(), shIn.end(), sd.mdOut.pixU) != shIn.end();
	      if (debug && debugHitLevel && hasIOU){
		std::cout<<"SDI: found matching hit for "<<lIn<<" IOU at "<<sd.mdOut.pixU<<std::endl;
	      }
	      int scoreIn = hasIRL + hasIRU + hasIOL + hasIOU;
	      int patternIn = hasIRL + (hasIRU<<1) + (hasIOL<<2) + (hasIOU<<3);
	      if (debug && scoreIn > 2){
		std::cout<<"SDI: match on L"<< lIn <<": Have "<<scoreIn<<" matches with pattern "<<patternIn<<std::endl;
	      }
	      if (scoreIn >= 3) hasSDIn_3of4 = true;
	      if (scoreIn >= 4){
		hasSDIn_4of4 = true;
		vSDIn_4of4.push_back(sd);
	      }
	    }//SD matching Inner
	    
	    for (auto sd : mockLayerSDfwDNcm[lOut]){
	      auto const& shOut = simHits[lOut];
	      
	      bool hasIRL = std::find(shOut.begin(), shOut.end(), sd.mdRef.pixL) != shOut.end();
	      if (debug && debugHitLevel && hasIRL){
		std::cout<<"SDO: found matching hit for "<<lOut<<" IRL at "<<sd.mdRef.pixL<<std::endl;
	      }
	      bool hasIRU = std::find(shOut.begin(), shOut.end(), sd.mdRef.pixU) != shOut.end();
	      if (debug && debugHitLevel && hasIRU){
		std::cout<<"SDO: found matching hit for "<<lOut<<" IRU at "<<sd.mdRef.pixU<<std::endl;
	      }
	      bool hasIOL = std::find(shOut.begin(), shOut.end(), sd.mdOut.pixL) != shOut.end();
	      if (debug && debugHitLevel && hasIOL){
		std::cout<<"SDO: found matching hit for "<<lOut<<" IOL at "<<sd.mdOut.pixL<<std::endl;
	      }
	      bool hasIOU = std::find(shOut.begin(), shOut.end(), sd.mdOut.pixU) != shOut.end();
	      if (debug && debugHitLevel && hasIOU){
		std::cout<<"SDO: found matching hit for "<<lOut<<" IOU at "<<sd.mdOut.pixU<<std::endl;
	      }
	      int scoreOut = hasIRL + hasIRU + hasIOL + hasIOU;
	      int patternOut = hasIRL + (hasIRU<<1) + (hasIOL<<2) + (hasIOU<<3);
	      if (debug && scoreOut > 2){
		std::cout<<"SDO: match on L"<< lOut <<": Have "<<scoreOut<<" matches with pattern "<<patternOut<<std::endl;
	      }
	      if (scoreOut >= 3) hasSDOut_3of4 = true;
	      if (scoreOut >= 4){
		hasSDOut_4of4 = true;
		vSDOut_4of4.push_back(sd);
	      }
	    }//SD matching Outer
	    if (hasSDIn_3of4 && hasSDOut_3of4){
	      ha_num2SD_3of4_any_pt[iSDLL]->Fill(tpPt);
	    }
	    if (hasSDIn_4of4 && hasSDOut_4of4){
	      ha_num2SD_4of4_pt[iSDLL]->Fill(tpPt);
	    }

	    enum SDLSelectFlags { deltaZ = 0, deltaZPointed, slope, dAlphaIn, dAlphaOut, dBeta};
	    std::vector<std::pair<SDLink, int> > vSDLwInfo_4of4;
	    for (auto const& sdIn : vSDIn_4of4){
	      float rtIn = sdIn.r3.Pt();
	      float zIn = sdIn.r3.z();
	      float dSDIn = sdIn.mdOut.r3.Pt() - sdIn.mdRef.r3.Pt();
	      float dzSDIn = sdIn.mdOut.r3.z() - sdIn.mdRef.r3.z();
	      for (auto const& sdOut : vSDOut_4of4){
		int sdlFlag = 0;
		
		float rtOut = sdOut.r3.Pt();
		float zOut = sdOut.r3.z();
		//apply some loose Z compatibility
		//FIXME: refine using inner layer directions (can prune later)
		float zGeom = lIn >= 5 && lIn <= 7 && lOut >= 5 && lOut <= 7 ? 0.3 : 10;//twice the macro-pixel or strip size
		float zLo = rtOut/rtIn*(zIn - 15.) - zGeom; //15 for the luminous ; zGeom for z geom unit size
		float zHi = rtOut/rtIn*(zIn + 15.) + zGeom;

		if (zOut > zLo && zOut < zHi) sdlFlag |= 1 << SDLSelectFlags::deltaZ;

		if (lIn>=5 && lIn <=6){//can point to the z pos in lOut
		  float zLo = zIn + (dzSDIn - zGeom*sqrt(2.))/dSDIn*(rtOut - rtIn) - zGeom;
		  float zHi = zIn + (dzSDIn + zGeom*sqrt(2.))/dSDIn*(rtOut - rtIn) + zGeom;
		  if (zOut > zLo && zOut < zHi) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		} else {
		  //the flag is set to pass here
		  sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		}

		auto midR3 = 0.5*(sdIn.r3 + sdOut.r3);
		double dPhi = midR3.DeltaPhi(sdOut.r3 - sdIn.r3);
		double rt = 0.5*(sdIn.r3.Pt() + sdOut.r3.Pt());
		const float ptCut = 1.0;
		const float sdlSlope = rt/175.67/ptCut;
		const float sdlThetaMulsF = 0.015*sqrt(0.2);
		const float sdlMuls = sdlThetaMulsF*3./ptCut*4;//will need a better guess than x4?
		const float sdlPVoff = 0.1/rt;
		const float sdlCut = sdlSlope + sqrt(sdlMuls*sdlMuls + sdlPVoff*sdlPVoff);
	    
		if (std::abs(dPhi) < sdlCut ) sdlFlag |= 1 << SDLSelectFlags::slope;

		double betaIn;
		double betaOut;
		if (mockMode == 0){
		  betaIn = sdIn.alpha - sdIn.r3.DeltaPhi(sdOut.r3 - sdIn.r3);
		  betaOut = - sdOut.alpha + sdOut.r3.DeltaPhi(sdOut.r3 - sdIn.r3); //to match sign for correct match	      
		}
		else if (mockMode == 1 || mockMode == 3){
		  //need a symmetric choice of end-points to achieve partial cancelation
		  betaIn  = sdIn.alpha - sdIn.r3.DeltaPhi(sdOut.mdOut.r3 - sdIn.r3);
		  betaOut = -sdOut.mdOut.r3.DeltaPhi(sdOut.mdOut.r3 - sdOut.mdRef.r3) + sdOut.mdOut.r3.DeltaPhi(sdOut.mdOut.r3 - sdIn.r3);	      
		}
		else{
		  betaIn = -99;
		  betaOut = 999;
		}
		
		//loose angle compatibility
		float dAlpha_Bfield = (rtOut - rtIn)/175.67/ptCut;
		float dSDOut = sdOut.mdOut.r3.Pt() - sdOut.mdRef.r3.Pt();	
		float dAlpha_res = 0.02/std::min(dSDIn, dSDOut);//2-strip difference; use the smallest SD separation
		float dAlpha_compat = dAlpha_Bfield + dAlpha_res;
		if (std::abs(sdIn.alpha- dPhi) < dAlpha_compat) sdlFlag |=  1 << SDLSelectFlags::dAlphaIn;
		if (std::abs(sdOut.alpha- dPhi) < dAlpha_compat) sdlFlag |= 1 << SDLSelectFlags::dAlphaOut;

		//now the actual segment linking magic
		float betaAv = 0.5*(betaIn + betaOut);
		float dr = (sdOut.r3 - sdIn.r3).Perp();
		//pt*175.67/2. = R
		//R*sin(betaAv) = pt*175.67/2*sin(betaAv) = dr/2 => pt = dr/175.67/sin(betaAv);
		float pt_beta = dr/175.67/sin(betaAv);
		float pt_betaIn = dr/175.67/sin(betaIn);
		float pt_betaOut = dr/175.67/sin(betaOut);
		float dBetaRes = dAlpha_res;
		float dBetaMuls = sdlThetaMulsF*3./std::min(pt_beta, 7.0f);//need to confirm the range-out value of 7 GeV
		float dBetaCut = sqrt(dBetaRes*dBetaRes*2.0 + dBetaMuls*dBetaMuls);
		float dBeta = betaIn - betaOut;

		if (std::abs(dBeta) < dBetaCut) sdlFlag |= 1 << SDLSelectFlags::dBeta;

		//		std::cout<<"TP with pt "<<tpPt<<" has matching SDL flags "<<sdlFlag<<std::endl;
		SDLink sdl;
		sdl.sdIn = sdIn;
		sdl.sdOut = sdOut;
		sdl.lIn = lIn;
		sdl.lOut = lOut;
		sdl.alpha = dPhi;
		sdl.betaIn = betaIn;
		sdl.betaOut = betaOut;
		sdl.pt = pt_beta;
		sdl.ptIn = pt_betaIn;
		sdl.ptOut = pt_betaOut;
		vSDLwInfo_4of4.push_back({sdl, sdlFlag});
	      }
	    }

	    const int m_0 = 1;
	    const int m_01 = m_0 | (1 << 1);
	    const int m_012 = m_01 | (1 << 2);
	    const int m_0123 = m_012 | (1 << 3);
	    const int m_01234 = m_0123 | (1 << 4);
	    const int m_012345 = m_01234 | (1 << 5);

	    bool h_0, h_01, h_012, h_0123, h_01234, h_012345;
	    h_0 = h_01 = h_012 = h_0123 = h_01234 = h_012345 = false;
	    for (auto const& sdlf : vSDLwInfo_4of4){
	      if ((sdlf.second & m_0) == m_0 ) h_0 = true;
	      if ((sdlf.second & m_01) == m_01 ) h_01 = true;
	      if ((sdlf.second & m_012) == m_012 ) h_012 = true;
	      if ((sdlf.second & m_0123) == m_0123 ) h_0123 = true;
	      if ((sdlf.second & m_01234) == m_01234 ) h_01234 = true;
	      
	      if ((sdlf.second & m_012345) == m_012345 ) h_012345 = true;

	      if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4 && h_01234){
		auto const& sdl = sdlf.first;
		ha_SDL_dBeta_zoom_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_zoom2_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaIn_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn, sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn, sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaOut, sdl.betaIn - sdl.betaOut);

		if (tpPt > 1.5 && iSDLL == SDL_L5to7){ //&& !h_012345 && std::abs(sdl.betaIn - sdl.betaOut)> 0.01){
		  std::cout<<" tPt "<<tpPt
			   <<" tEta " << tpEta<<" tPhi "<<tpPhi
			   <<" dB "<<sdl.betaIn - sdl.betaOut
			   <<" bI "<<sdl.betaIn
			   <<" bO "<<sdl.betaOut
			   <<" sI "<<sdl.sdIn.alpha
			   <<" sO "<<sdl.sdOut.alpha
			   <<" mI "<<sdl.sdIn.mdRef.alpha
			   <<" mO "<<sdl.sdOut.mdRef.alpha
			   <<" mRI ("<<sdl.sdIn.mdRef.r3.Pt()<<", "<<sdl.sdIn.mdRef.r3.Eta()<<", "<<sdl.sdIn.mdRef.r3.Phi()<<")"
			   <<" mRO ("<<sdl.sdOut.mdRef.r3.Pt()<<", "<<sdl.sdOut.mdRef.r3.Eta()<<", "<<sdl.sdOut.mdRef.r3.Phi()<<")"
			   <<std::endl;

		  std::vector<TVector3> r3Ins;
		  std::vector<TVector3> r3Outs;
		  std::vector<TVector3> p3Ins;
		  std::vector<TVector3> p3Outs;
		  
		  for (int iPix = 0; iPix< nPix; ++iPix){
		    
		    int iipix = sim_pixelIdx()[iSim][iPix];
		    PixelHit pixH(iipix);
		    int lay = pixH.lay;
		    
		    // int iProcess = pix_process()[iPix];
		    // int ibx = pix_bunchXing()[iPix];	  
		    // bool isPrimaryAny = (iProcess == 2 && ibx == 0);
		    
		    if (pixH.isBarrel){
		      if (lay >= minLayer){
			std::cout<<" "<<lay<<" "<<iipix<<std::endl;
			
			if (pixH.p3s.Pt()>0.8*tpPt){
			  if (lay == layersSDL[iSDLL][0]){
			    r3Ins.push_back(pixH.r3s);
			    p3Ins.push_back(pixH.p3s);
			  } else if (lay == layersSDL[iSDLL][1]){
			    r3Outs.push_back(pixH.r3s);
			    p3Outs.push_back(pixH.p3s);			    
			  }
			  pixH.print("\t");
			}
		      }
		    }		  
		  }//iPix; pixelhits

		  //get true values
		  int nIns = r3Ins.size();
		  int nOuts = r3Outs.size();
		  for (int iIn = 0; iIn < nIns; ++iIn){
		    for (int iOut = 0; iOut < nOuts; ++iOut){
		      float alphaIn = r3Ins[iIn].DeltaPhi(p3Ins[iIn]);
		      float alphaOut = r3Outs[iOut].DeltaPhi(p3Outs[iOut]);
		      
		      float betaIn = r3Ins[iIn].DeltaPhi(p3Ins[iIn]) - r3Ins[iIn].DeltaPhi(r3Outs[iOut] - r3Ins[iIn]);
		      float betaOut = - r3Outs[iOut].DeltaPhi(p3Outs[iOut]) + r3Outs[iOut].DeltaPhi(r3Outs[iOut] - r3Ins[iIn]);
		      std::cout<<" Sim betas: "<<betaIn<<" "<<betaOut<<" "<<betaIn - betaOut
			       <<" alphas "<<alphaIn<<" "<<alphaOut
			       <<" ("<<r3Ins[iIn].Pt()<<", "<<r3Ins[iIn].Eta()<<", "<<r3Ins[iIn].Phi()<<") "
			       <<"- ("<<r3Outs[iOut].Pt()<<", "<<r3Outs[iOut].Eta()<<", "<<r3Outs[iOut].Phi()<<") "
			       <<std::endl;
		    }
		  }
		  
		}//if (tpPt > 1.5 && iSDLL == SDL_L5to7 ... debug
	      }
	      if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4 && h_012345){
		auto const& sdl = sdlf.first;
		ha_SDL_dBeta_betaIn_pass[iSDLL]->Fill(sdl.betaIn, sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaIn_zoom_pass[iSDLL]->Fill(sdl.betaIn, sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaOut_zoom_pass[iSDLL]->Fill(sdl.betaOut, sdl.betaIn - sdl.betaOut);
	      }

	    }

	    if (h_0) ha_num2SD_w0_4of4_pt[iSDLL]->Fill(tpPt);
	    if (h_01) ha_num2SD_w01_4of4_pt[iSDLL]->Fill(tpPt);
	    if (h_012) ha_num2SD_w012_4of4_pt[iSDLL]->Fill(tpPt);
	    if (h_0123) ha_num2SD_w0123_4of4_pt[iSDLL]->Fill(tpPt);
	    if (h_01234) ha_num2SD_w01234_4of4_pt[iSDLL]->Fill(tpPt);
	    if (h_012345) ha_num2SD_w012345_4of4_pt[iSDLL]->Fill(tpPt);
	    
	    if (mockLayerSDLsDNcm[iSDLL]){
	      for (auto& sdl : *mockLayerSDLsDNcm[iSDLL]){
		if (! (sdl.lIn == lIn && sdl.lOut == lOut )) continue;
		auto const& shIn = simHits[sdl.lIn];
		auto const& shOut = simHits[sdl.lOut];

		bool hasIRL = std::find(shIn.begin(), shIn.end(), sdl.sdIn.mdRef.pixL) != shIn.end();
		if (debug && debugHitLevel && hasIRL){
		  std::cout<<"SDL: found matching hit for "<<sdl.lIn<<" IRL at "<<sdl.sdIn.mdRef.pixL<<std::endl;
		}
		bool hasIRU = std::find(shIn.begin(), shIn.end(), sdl.sdIn.mdRef.pixU) != shIn.end();
		if (debug && debugHitLevel && hasIRU){
		  std::cout<<"SDL: found matching hit for "<<sdl.lIn<<" IRU at "<<sdl.sdIn.mdRef.pixU<<std::endl;
		}
		bool hasIOL = std::find(shIn.begin(), shIn.end(), sdl.sdIn.mdOut.pixL) != shIn.end();
		if (debug && debugHitLevel && hasIOL){
		  std::cout<<"SDL: found matching hit for "<<sdl.lIn<<" IOL at "<<sdl.sdIn.mdOut.pixL<<std::endl;
		}
		bool hasIOU = std::find(shIn.begin(), shIn.end(), sdl.sdIn.mdOut.pixU) != shIn.end();
		if (debug && debugHitLevel && hasIOU){
		  std::cout<<"SDL: found matching hit for "<<sdl.lIn<<" IOU at "<<sdl.sdIn.mdOut.pixU<<std::endl;
		}
		int scoreIn = hasIRL + hasIRU + hasIOL + hasIOU;
		int patternIn = hasIRL + (hasIRU<<1) + (hasIOL<<2) + (hasIOU<<3);
		if (debug && scoreIn > 1){
		  std::cout<<"Inner match on L"<< sdl.lIn <<": Have "<<scoreIn<<" matches with pattern "<<patternIn<<std::endl;
		}

		bool hasORL = std::find(shOut.begin(), shOut.end(), sdl.sdOut.mdRef.pixL) != shOut.end();
		if (debug && debugHitLevel && hasORL){
		  std::cout<<"SDL: found matching hit for "<<sdl.lOut<<" IRL at "<<sdl.sdOut.mdRef.pixL<<std::endl;
		}
		bool hasORU = std::find(shOut.begin(), shOut.end(), sdl.sdOut.mdRef.pixU) != shOut.end();
		if (debug && debugHitLevel && hasORU){
		  std::cout<<"SDL: found matching hit for "<<sdl.lOut<<" IRU at "<<sdl.sdOut.mdRef.pixU<<std::endl;
		}
		bool hasOOL = std::find(shOut.begin(), shOut.end(), sdl.sdOut.mdOut.pixL) != shOut.end();
		if (debug && debugHitLevel && hasOOL){
		  std::cout<<"SDL: found matching hit for "<<sdl.lOut<<" IOL at "<<sdl.sdOut.mdOut.pixL<<std::endl;
		}
		bool hasOOU = std::find(shOut.begin(), shOut.end(), sdl.sdOut.mdOut.pixU) != shOut.end();
		if (debug && debugHitLevel && hasOOU){
		  std::cout<<"SDL: found matching hit for "<<sdl.lOut<<" IOU at "<<sdl.sdOut.mdOut.pixU<<std::endl;
		}
		int scoreOut = hasORL + hasORU + hasOOL + hasOOU;
		int patternOut = hasORL + (hasORU<<1) + (hasOOL<<2) + (hasOOU<<3);
		if (debug && scoreOut > 1){
		  std::cout<<"Outer match on L"<< sdl.lOut<<": Have "<<scoreOut<<" matches with pattern "<<patternOut<<std::endl;
		}

		if (scoreIn >=3 && scoreOut>= 3){
		  matchingSDLs_byHit3of4[iSDLL].push_back(sdl);
		  if (scoreIn >=4 && scoreOut>= 4){
		    matchingSDLs_byHit4of4[iSDLL].push_back(sdl);
		    sdl.hasMatch_byHit4of4 = true;
		  }
		}
	      }//for (auto sdl : mockLayerSDLsDNcm[iSDLL]){
	    }//	if (mockLayerSDLsDNcm[iSDLL]){

	    bool hasMatch = false;
	    //matching is done: fill numerators
	    if (! matchingSDLs_byHit3of4[iSDLL].empty()){
	      ha_numSDL_3of4_any_pt[iSDLL]->Fill(tpPt);
	      if (debug) std::cout<<"\t have 3/4 match "<<std::endl;
	      hasMatch = true;
	    }
	    if (! matchingSDLs_byHit4of4[iSDLL].empty()){
	      ha_numSDL_4of4_pt[iSDLL]->Fill(tpPt);
	      if (debug) std::cout<<"\t have 4/4 match "<<std::endl;
	      hasMatch = true;
	    }
	    
	    if (debug && ! hasMatch && iSDLL != SDL_L5to9){
	      for (auto i : simHits[lIn]){ auto ph = PixelHit(i); ph.print("\tNM for: ");}
	      for (auto i : simHits[lOut]){ auto ph = PixelHit(i); ph.print("\tNM for: ");}
	    }

	  }// TP has hits in SDL layers	 	  
	}//for (int iSDLL = 0; iSDLL< SDL_LMAX; ++iSDLL){
      }//TPs
	

      //fake rates
      for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
	auto mSDLs = mockLayerSDLsDNcm[iSDL];
	if (mSDLs == nullptr ) continue;
	for (auto sdl : *mSDLs){
	  auto pt = sdl.pt;
	  auto eta = (sdl.sdIn.r3.Eta() + sdl.sdOut.r3.Eta())*0.5;

	  ha_SDLreco_all_pt[iSDL]->Fill(pt);
	  ha_SDLreco_all_eta[iSDL]->Fill(eta);
	  if (sdl.hasMatch_byHit4of4){
	    ha_SDLreco_4of4_pt[iSDL]->Fill(pt);
	    ha_SDLreco_4of4_eta[iSDL]->Fill(eta);	  
	  } else {
	    ha_SDLreco_no4of4_pt[iSDL]->Fill(pt);
	    ha_SDLreco_no4of4_eta[iSDL]->Fill(eta);	  
	  }
	}
      }//iSDL
      
      timerA[T_timeValidation].Stop();
      
      //link the links to TrackLinks
      std::vector<TrackLink> tracks;
      int countMatchMidPoint = 0;
      for (auto sdlIn : mockLayer5to7SDLfwDNcm ){
	bool hasOuter = false;

	for (auto sdlOut : mockLayer7to9SDLfwDNcm){
	  if (sdlIn.lOut == sdlOut.lIn && sdlIn.iOut == sdlOut.iIn){
	    //shared mid-point
	    double dPt = sdlIn.pt - sdlOut.pt;
	    int iirL = sdlIn.sdIn.mdRef.pixL;
	    int iirU = sdlIn.sdIn.mdRef.pixU;
	    int iioL = sdlIn.sdIn.mdOut.pixL;
	    int iioU = sdlIn.sdIn.mdOut.pixU;

	    int iorL = sdlIn.sdOut.mdRef.pixL;
	    int iorU = sdlIn.sdOut.mdRef.pixU;
	    int iooL = sdlIn.sdOut.mdOut.pixL;
	    int iooU = sdlIn.sdOut.mdOut.pixU;
	    
	    int oorL = sdlOut.sdOut.mdRef.pixL;
	    int oorU = sdlOut.sdOut.mdRef.pixU;
	    int oooL = sdlOut.sdOut.mdOut.pixL;
	    int oooU = sdlOut.sdOut.mdOut.pixU;
	    
	    TVector3 iirLp3(pix_pxsim()[iirL], pix_pysim()[iirL], pix_pzsim()[iirL]);
	    TVector3 iirUp3(pix_pxsim()[iirU], pix_pysim()[iirU], pix_pzsim()[iirU]);
	    TVector3 iioLp3(pix_pxsim()[iioL], pix_pysim()[iioL], pix_pzsim()[iioL]);
	    TVector3 iioUp3(pix_pxsim()[iioU], pix_pysim()[iioU], pix_pzsim()[iioU]);

	    TVector3 iorLp3(pix_pxsim()[iorL], pix_pysim()[iorL], pix_pzsim()[iorL]);
	    TVector3 iorUp3(pix_pxsim()[iorU], pix_pysim()[iorU], pix_pzsim()[iorU]);
	    TVector3 iooLp3(pix_pxsim()[iooL], pix_pysim()[iooL], pix_pzsim()[iooL]);
	    TVector3 iooUp3(pix_pxsim()[iooU], pix_pysim()[iooU], pix_pzsim()[iooU]);

	    TVector3 oorLp3(pix_pxsim()[oorL], pix_pysim()[oorL], pix_pzsim()[oorL]);
	    TVector3 oorUp3(pix_pxsim()[oorU], pix_pysim()[oorU], pix_pzsim()[oorU]);
	    TVector3 oooLp3(pix_pxsim()[oooL], pix_pysim()[oooL], pix_pzsim()[oooL]);
	    TVector3 oooUp3(pix_pxsim()[oooU], pix_pysim()[oooU], pix_pzsim()[oooU]);


	    /*
	    std::cout<<countMatchMidPoint<<"\t"<<sdlIn.pt<<" "<<sdlOut.pt
		     <<" vs MC "
		     <<" "<<iirL<<" "<<iirLp3.Pt()<<";"
		     <<" "<<iirU<<" "<<iirUp3.Pt()<<";"
		     <<" "<<iioL<<" "<<iioLp3.Pt()<<";"
		     <<" "<<iioU<<" "<<iioUp3.Pt()<<";;"

		     <<" "<<iorL<<" "<<iorLp3.Pt()<<";"
		     <<" "<<iorU<<" "<<iorUp3.Pt()<<";"
		     <<" "<<iooL<<" "<<iooLp3.Pt()<<";"
		     <<" "<<iooU<<" "<<iooUp3.Pt()<<";;"
	      
		     <<" "<<oorL<<" "<<oorLp3.Pt()<<";"
		     <<" "<<oorU<<" "<<oorUp3.Pt()<<";"
		     <<" "<<oooL<<" "<<oooLp3.Pt()<<";"
		     <<" "<<oooU<<" "<<oooUp3.Pt()<<""
		     <<std::endl;
	    */
	    countMatchMidPoint++;
	  }//match in-out at mid-point
	  
	}//sdlOut
      }//sdlIn
      
      
      std::cout<<"Print stats"<<std::endl;
      for (int iL = minLayer; iL <= nLayers; ++iL){
	int nSDs1GeV = 0;
	int nSDs1GeVMatch4 = 0;
	int nSDs2GeV = 0;
	int nSDs2GeVMatch4 = 0;
	for (auto sd : mockLayerSDfwDNcm[iL]){
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
		 <<" outLos " <<mockLayerMDfwDNcmLower[iL].size()
		 <<" outUps " <<mockLayerMDfwDNcmUpper[iL].size()
		 <<" outMDs " <<mockLayerMDfwDNcm[iL].size()
		 <<" outSDs "<< mockLayerSDfwDNcm[iL].size()
		 <<" n1Any "<<nSDs1GeV<<" n1All "<<nSDs1GeVMatch4
		 <<" n2Any "<<nSDs2GeV<<" n2All "<<nSDs2GeVMatch4
		 <<std::endl;
      }//iL
      std::cout<<"\t\tsdl5-7 "<<mockLayer5to7SDLfwDNcm.size()
	       <<" sdl7-9 "<<mockLayer7to9SDLfwDNcm.size()
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

  for (int iT = 0; iT < T_N; ++iT){
    std::cout<<"Timer results for "<<timerNameA[iT]<<std::endl;
    timerA[iT].Print("m");
    std::cout<<"----------------------------------"<<std::endl;
  }
  
  std::cout<<__LINE__<<" make efficiencies "<<std::endl;
  std::array<TEfficiency*, SDL_LMAX> ha_eff8MH_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff4MD_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_3of4_any_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_3of4_any_pt;

  std::array<TEfficiency*, SDL_LMAX> ha_eff4MD_den8MH_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_den8MH_3of4_any_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w01_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w012_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0123_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w01234_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w012345_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_den8MH_3of4_any_pt;

  std::array<TEfficiency*, SDL_LMAX> ha_fakeSDL_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_fakeSDL_4of4_eta;
  for (int i = 0; i< SDL_LMAX; ++i){
    auto iMin = layersSDL[i][0];
    auto iMax = layersSDL[i][1];
    ha_eff8MH_pt[i] = new TEfficiency(*ha_num8MH_pt[i], *ha_denSDL_pt[i]);
    ha_eff8MH_pt[i]->SetTitle(Form("h_eff8MH_%dto%d_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_eff4MD_pt[i] = new TEfficiency(*ha_num4MD_pt[i], *ha_denSDL_pt[i]);
    ha_eff4MD_pt[i]->SetTitle(Form("h_eff4MD_%dto%d_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_4of4_pt[i] = new TEfficiency(*ha_num2SD_4of4_pt[i], *ha_denSDL_pt[i]);
    ha_eff2SD_4of4_pt[i]->SetTitle(Form("h_eff2SD_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_3of4_any_pt[i] = new TEfficiency(*ha_num2SD_3of4_any_pt[i], *ha_denSDL_pt[i]);
    ha_eff2SD_3of4_any_pt[i]->SetTitle(Form("h_eff2SD_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_4of4_pt[i] = new TEfficiency(*ha_numSDL_4of4_pt[i], *ha_denSDL_pt[i]);
    ha_effSDL_4of4_pt[i]->SetTitle(Form("h_effSDL_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_3of4_any_pt[i] = new TEfficiency(*ha_numSDL_3of4_any_pt[i], *ha_denSDL_pt[i]);
    ha_effSDL_3of4_any_pt[i]->SetTitle(Form("h_effSDL_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_eff4MD_den8MH_pt[i] = new TEfficiency(*ha_num4MD_pt[i], *ha_num8MH_pt[i]);
    ha_eff4MD_den8MH_pt[i]->SetTitle(Form("h_eff4MD_den8MH_%dto%d_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_den8MH_3of4_any_pt[i] = new TEfficiency(*ha_num2SD_3of4_any_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_den8MH_3of4_any_pt[i]->SetTitle(Form("h_eff2SD_den8MH_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w0_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w0_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w0_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w0_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w01_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w01_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w01_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w01_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w012_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w012_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w012_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w012_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w0123_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w0123_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w0123_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w0123_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w01234_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w01234_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w01234_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w01234_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w012345_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w012345_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w012345_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w012345_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_den8MH_4of4_pt[i] = new TEfficiency(*ha_numSDL_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_effSDL_den8MH_4of4_pt[i]->SetTitle(Form("h_effSDL_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_den8MH_3of4_any_pt[i] = new TEfficiency(*ha_numSDL_3of4_any_pt[i], *ha_num8MH_pt[i]);
    ha_effSDL_den8MH_3of4_any_pt[i]->SetTitle(Form("h_effSDL_den8MH_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_fakeSDL_4of4_pt[i] = new TEfficiency(*ha_SDLreco_no4of4_pt[i], *ha_SDLreco_all_pt[i]);
    ha_fakeSDL_4of4_pt[i]->SetTitle(Form("h_fakeSDL_%dto%d_4of4_pt", iMin, iMax));

    ha_fakeSDL_4of4_eta[i] = new TEfficiency(*ha_SDLreco_no4of4_eta[i], *ha_SDLreco_all_eta[i]);
    ha_fakeSDL_4of4_eta[i]->SetTitle(Form("h_fakeSDL_%dto%d_4of4_eta", iMin, iMax));
  }

  
  if (drawPlots){
    std::cout<<__LINE__<<" draw and print "<<std::endl;
    for (int iL = 5; iL <= nLayers; ++iL){
      if (iL != 5 && iL != 10) continue;

    }//nLayers

    {
      auto h_all = h_sdl0to5_dBeta_NM1dBeta_all;
      auto h_pass = h_sdl0to5_dBeta_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to5_dBeta_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to7_dBeta_NM1dBeta_all;
      auto h_pass = h_sdl0to7_dBeta_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to7_dBeta_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl5to7_dBeta_NM1dBeta_all;
      auto h_pass = h_sdl5to7_dBeta_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl5to7_dBeta_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl7to9_dBeta_NM1dBeta_all;
      auto h_pass = h_sdl7to9_dBeta_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl7to9_dBeta_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to5_dBeta_zoom_NM1dBeta_all;
      auto h_pass = h_sdl0to5_dBeta_zoom_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to5_dBeta_zoom_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to7_dBeta_zoom_NM1dBeta_all;
      auto h_pass = h_sdl0to7_dBeta_zoom_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to7_dBeta_zoom_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl5to7_dBeta_zoom_NM1dBeta_all;
      auto h_pass = h_sdl5to7_dBeta_zoom_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl5to7_dBeta_zoom_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl7to9_dBeta_zoom_NM1dBeta_all;
      auto h_pass = h_sdl7to9_dBeta_zoom_NM1dBeta_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl7to9_dBeta_zoom_NM1dBeta_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    // << slices in betaIn or ptIn
    {
      auto h_all = h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_all;
      auto h_pass = h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn0to2_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_all;
      auto h_pass = h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn0to2_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_all;
      auto h_pass = h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn0to2_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_all;
      auto h_pass = h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn0to2_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_all;
      auto h_pass = h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn3to5_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_all;
      auto h_pass = h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn3to5_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_all;
      auto h_pass = h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn3to5_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_all;
      auto h_pass = h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn3to5_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_all;
      auto h_pass = h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to5_dBeta_zoom_NM1dBeta_ptIn7toInf_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all;
      auto h_pass = h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl0to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all;
      auto h_pass = h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl5to7_dBeta_zoom_NM1dBeta_ptIn7toInf_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h_all = h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_all;
      auto h_pass = h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_pass;
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_all->SetMinimum(0);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_sdl7to9_dBeta_zoom_NM1dBeta_ptIn7toInf_all_vs_pass_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    // >> in pt slices 
    {
      auto h2 = h2_sdl0to5_dBeta_betaIn_NM1dBeta_all;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_sdl0to5_dBeta_betaIn_NM1dBeta_all_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h2 = h2_sdl0to7_dBeta_betaIn_NM1dBeta_all;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_sdl0to7_dBeta_betaIn_NM1dBeta_all_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h2 = h2_sdl5to7_dBeta_betaIn_NM1dBeta_all;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_sdl5to7_dBeta_betaIn_NM1dBeta_all_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    {
      auto h2 = h2_sdl7to9_dBeta_betaIn_NM1dBeta_all;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_sdl7to9_dBeta_betaIn_NM1dBeta_all_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }

    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_NM1dBeta_8MH[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom2_NM1dBeta_8MH[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom2_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_NM1dBeta_8MH[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }

    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_pass[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_pass_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_zoom_pass[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_zoom_pass_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaOut_zoom_pass[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaOut_zoom_pass_%dto%d_mm%d_D%1.1fcm_us%d.png",  layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
      
    }

    //efficiencies: num/den plots
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto hden = ha_denSDL_pt[iSDL];
      auto hnum1 = ha_numSDL_3of4_any_pt[iSDL];
      auto hnum2 = ha_numSDL_4of4_pt[iSDL];

      auto cn = hden->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      hden->SetStats(0);

      hden->SetLineWidth(2);
      hnum1->SetLineWidth(2);
      hnum2->SetLineWidth(2);
      hden->SetLineColor(kBlack);
      hnum1->SetLineColor(kRed);
      hnum2->SetLineColor(kBlue);

      hden->Draw();
      hden->SetMinimum(0.9);
      auto ax = hden->GetXaxis();
      ax->SetRangeUser(0.51, ax->GetXmax());
      gPad->SetLogx();
      gPad->SetLogy();
      hnum1->Draw("same");
      hnum2->Draw("same");

      auto leg = new TLegend(0.7, 0.65, 0.9, 0.9);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(hden, "Denominator");
      leg->AddEntry(hnum1, "Num: 3 of 4");
      leg->AddEntry(hnum2, "Num: 4 of 4");
      leg->Draw();
      gPad->SaveAs(Form("h_denVSnums_SDL_%dto%d_pt_mm%d_D%1.1fcm_us%d.png", layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
    }

    //efficiencies: eff plots 3/4 vs 4/4
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff1 = ha_effSDL_3of4_any_pt[iSDL];
      auto heff2 = ha_effSDL_4of4_pt[iSDL];

      auto cn = heff1->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff1->Draw();
      heff2->Draw("same");

      heff1->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff1->SetLineColor(kRed);
      heff2->SetLineColor(kBlue);
      heff1->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff1->SetMarkerStyle(22);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff1->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(0.51, ax->GetXmax());

      auto leg = new TLegend(0.7, 0.15, 0.9, 0.3);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff1, "3 of 4");
      leg->AddEntry(heff2, "4 of 4");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_%dto%d_pt_mm%d_D%1.1fcm_us%d.png", layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
    }
    
    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heffH =  ha_eff8MH_pt[iSDL];
      auto heff =  ha_eff4MD_pt[iSDL];
      auto heff1 = ha_eff2SD_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_4of4_pt[iSDL];

      auto cn = heffH->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heffH->Draw();
      heff->Draw("same");
      heff1->Draw("same");
      heff2->Draw("same");

      heffH->SetLineWidth(2);
      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heffH->SetLineColor(kCyan);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff2->SetLineColor(kBlue);
      heffH->SetMarkerSize(0.8);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heffH->SetMarkerStyle(24);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heffH->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(0.51, ax->GetXmax());

      auto leg = new TLegend(0.6, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heffH, "8 MH match");
      leg->AddEntry(heff, "4 MD match");
      leg->AddEntry(heff1, "2 SD match");
      leg->AddEntry(heff2, "SDLink match");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_steps_4of4_%dto%d_pt_mm%d_D%1.1fcm_us%d.png", layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
    }

    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff01234->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff01234->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff01234->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff01234->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff01234->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(0.51, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff01234, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_4of4_%dto%d_pt_mm%d_D%1.1fcm_us%d.png", layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffset, useSeeds));
    }

    //fakes vs pt
    {      
      auto h57 =  ha_fakeSDL_4of4_pt[SDL_L5to7];
      auto h79 =  ha_fakeSDL_4of4_pt[SDL_L7to9];
      
      auto cn = h57->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h57->Draw();
      h79->Draw("same");

      h57->SetLineWidth(2);
      h79->SetLineWidth(2);
      h57->SetLineColor(kBlack);
      h79->SetLineColor(kRed);
      h57->SetMarkerSize(0.8);
      h79->SetMarkerSize(0.8);
      h57->SetMarkerStyle(21);
      h79->SetMarkerStyle(22);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = h57->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(1.0, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.86, 0.22);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->AddEntry(h57, "L5-L7", "LP");
      leg->AddEntry(h79, "L7-L9", "LP");
      leg->Draw();
      gPad->SaveAs(Form("h_fake_SDL_4of4_pt_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }
    //fakes vs pt
    {      
      auto h57 =  ha_fakeSDL_4of4_eta[SDL_L5to7];
      auto h79 =  ha_fakeSDL_4of4_eta[SDL_L7to9];
      
      auto cn = h57->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h57->Draw();
      h79->Draw("same");

      h57->SetLineWidth(2);
      h79->SetLineWidth(2);
      h57->SetLineColor(kBlack);
      h79->SetLineColor(kRed);
      h57->SetMarkerSize(0.8);
      h79->SetMarkerSize(0.8);
      h57->SetMarkerStyle(21);
      h79->SetMarkerStyle(22);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->PaintModified();
      auto pg = h57->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(-1.5, 1.5);

      auto leg = new TLegend(0.5, 0.15, 0.86, 0.22);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->AddEntry(h57, "L5-L7", "LP");
      leg->AddEntry(h79, "L7-L9", "LP");
      leg->Draw();
      gPad->SaveAs(Form("h_fake_SDL_4of4_eta_mm%d_D%1.1fcm_us%d.png", mockMode, sdOffset, useSeeds));
    }

  }//if drawPlots
  
  std::cout<<__LINE__<<" write to file "<<std::endl;
  TFile* outHistograms = new TFile(Form("outHistogramsSuperD_mm%d_D%1.1fcm_us%d.root", mockMode, sdOffset, useSeeds), "RECREATE");
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

