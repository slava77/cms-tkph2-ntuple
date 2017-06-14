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
#include <unistd.h>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStopwatch.h"

#include "TVector3.h"

#include "tkph2.cc"

namespace tkph2consts{
  constexpr int nLayersA = 15; //All layers: barrel 5-10, endcap 11-15
  constexpr int nLayersB = 10; //barrel layers
  constexpr float kRinv1GeVf = (2.99792458e-3*3.8);
  constexpr float k2Rinv1GeVf = kRinv1GeVf/2.;

  constexpr int miniMask = 0x3;
  constexpr int lowerId = 1;
  constexpr int upperIdDelta = 1;
  //in SLHC these were 0x4, 4, 4
}
using namespace tas;
using namespace tkph2consts;

//copy-pasted from reco::deltaPhi
template<typename T> inline
T  reduceRange(T x) {
  constexpr T o2pi = 1./(2.*M_PI);
  if (std::abs(x) <= T(M_PI)) return x;
  T n = std::round(x*o2pi);
  return x - n*T(2.*M_PI);
}
inline double deltaPhi(double phi1, double phi2) { 
  return reduceRange(phi1 - phi2);
}
inline double deltaPhi(float phi1, double phi2) {
  return deltaPhi(static_cast<double>(phi1), phi2);
}

inline double deltaPhi(double phi1, float phi2) {
  return deltaPhi(phi1, static_cast<double>(phi2));
}
inline float deltaPhi(float phi1, float phi2) { 
  return reduceRange(phi1 - phi2);
}

//check if this is the same as in the release
enum class HitType {
  Pixel = 0,
  Strip = 1,
  Glued = 2,
  Invalid = 3,
  Phase2OT = 4,
  Unknown = 99
};

/// track algorithm; partial copy from TrackBase.h
enum class TrackAlgorithm {
  undefAlgorithm = 0,
  ctf = 1, 
  duplicateMerge = 2,
  cosmics = 3,
  initialStep = 4,
  lowPtTripletStep = 5,
  pixelPairStep = 6,
  detachedTripletStep = 7,
  mixedTripletStep = 8,
  pixelLessStep = 9,
  tobTecStep = 10,
  jetCoreRegionalStep = 11,
  conversionStep = 12,
  muonSeededStepInOut = 13,
  muonSeededStepOutIn = 14,
  outInEcalSeededConv = 15, inOutEcalSeededConv = 16,
  nuclInter = 17,
  standAloneMuon = 18, globalMuon = 19, cosmicStandAloneMuon = 20, cosmicGlobalMuon = 21,
  // Phase1
  highPtTripletStep = 22, lowPtQuadStep = 23, detachedQuadStep = 24,
  reservedForUpgrades1 = 25, reservedForUpgrades2 = 26,
  bTagGhostTracks = 27,
  beamhalo = 28,
  gsf = 29
};

struct HitIndexWithType {
  HitIndexWithType(): indexWithType(((int(HitType::Unknown) & indexMask) << typeOffset)) {}
  HitIndexWithType(const int index, HitType htype): indexWithType(((int(htype) & indexMask) << typeOffset)
								  | (index & indexMask)) {}
  int index() const {return indexWithType & indexMask;}
  int type() const {return (indexWithType >> typeOffset) & typeMask;}

  static int index(const int indexWT) {return indexWT & indexMask;}
  static int typeInt(const int indexWT) {return (indexWT >> typeOffset) & typeMask;}
  static HitType type(const int indexWT) {return HitType(typeInt(indexWT));}
  static constexpr int indexOffset = 0;
  static constexpr int indexMask = 0x7FFFFF;
  static constexpr int typeOffset = 23;
  static constexpr int typeMask = 0x7F;

  const int indexWithType;

};

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

TVector3 r3FromPCA(const TVector3& p3, const float dxy, const float dz){
  const float pt = p3.Pt();
  const float p = p3.Mag();
  const float vz = dz*pt*pt/p/p;

  const float vx = dxy*p3.x()/pt - p3.y()/p*p3.z()/p*dz;
  const float vy = -dxy*p3.y()/pt - p3.x()/p*p3.z()/p*dz;
  return TVector3(vx, vy, vz);
}

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

template <typename AT>
  void createPtHistograms(AT& layerMD, const std::string& ext, std::vector<double>& ptBins){

  const int NL = layerMD.size();
  for (int iL = 1; iL <= NL-1; ++iL){
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


struct V3WithCache {
  TVector3 r3;
  float rt;
  float r;
  float phi;
  float rtRHin;
  float rtRHout;
  float phiRHin;
  float phiRHout;
  V3WithCache(TVector3 const& o3) : r3(o3), rt(r3.Pt()), r(r3.Mag()), phi(r3.Phi()),
				    rtRHin(0), rtRHout(0), phiRHin(0), phiRHout(0) {}
};

struct MiniDoublet {
  int pixL;
  int pixU;
  TVector3 r3;
  float alpha;
  float alphaRHmin = 0.f;
  float alphaRHmax = 0.f;
  float rt;
  float z;
  float r;
  float phi;
  float rtRHin = 0.f;
  float rtRHout = 0.f;
  float phiRHin = 0.f;
  float phiRHout = 0.f;
  int itp;//tp with most hits
  int ntp;//n hits with itp
  int itpLL;//the lower layer iTP
};

struct SuperDoublet {
  MiniDoublet mdRef;
  MiniDoublet mdOut;
  int iRef;
  int iOut;
  TVector3 r3; //may be different from plain mdRef.r3
  float rt;
  float rtInv;
  float z;
  TVector3 p3; //makes sense mostly for seed-based
  float alpha;
  float alphaOut;
  float alphaRHmin;
  float alphaOutRHmin;
  float alphaRHmax;
  float alphaOutRHmax;
  float dr;
  float d;
  float zeta;
  int itp;//tp with most hits
  int ntp;//n hits with itp
  int itpLL;//lower-layer TP indexing for post-ghost cleanup tracking
  int ntpLL;
};

struct SDLink {
  SuperDoublet sdIn;
  SuperDoublet sdOut;
  int lIn;//layer
  int iIn;
  int lOut;//layer
  int iOut;
  float alpha; //point-to-point wrt radial
  float alphaRHmin;
  float alphaRHmax;
  float betaIn;//SD angle wrt point-to-point
  float betaOut;
  float betaInRHmin;
  float betaInRHmax;
  float betaOutRHmin;
  float betaOutRHmax;
  float pt;
  float ptIn;
  float ptOut;

  int itp;//tp with most hits
  int ntp;//n hits with itp
  bool hasMatch_byHit4of4=false;
  bool sdInGhost=false;
  bool sdOutGhost=false;
  int itpLL;//lower-layer (wrt MD components) accounting; for simple ghost cleaning
  int ntpLL;
};

struct TrackLink {
  std::vector<const SDLink*> links;
  
  float pt;
  float eta;
  float phi;
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

int layer(int lay, int det){
  //try to restore the SLHC indexing:
  // barrel: 5-10 OT
  // endcap: 11-15 OT disks
  // IT is not handled: "4" layers are filled from seed layer counting (no geometry ref)
  if (det == 4) return 10+lay;//OT endcap
  if (det == 5) return 4+lay;//OT barrel
  
  return lay;
}

int ph2_isBarrel(int iph2){
  return ph2_det()[iph2] == 5;
}

struct SimHit {
  SimHit() {}
  SimHit(int i):
    r3s(simhit_x()[i], simhit_y()[i], simhit_z()[i]),
    p3s(simhit_px()[i], simhit_py()[i], simhit_pz()[i]),
    ind(i),
    lay(layer(simhit_lay()[i], simhit_det()[i])),
    pdgId(simhit_particle()[i]),
    process(simhit_process()[i]),
    bx(999),
    evt(999),
    isBarrel(false)
  {
    simTkIdx = simhit_simTrkIdx()[i];
    if (simTkIdx >= 0){
      bx = sim_bunchCrossing()[simTkIdx];
      evt = sim_event()[simTkIdx];
    }
    auto det = simhit_det()[i];
    isBarrel = (det == 1 || det == 3 || det == 5);//Magically, it works for phase-0/1

    hitType = HitType::Unknown;
    for (auto rht : simhit_hitType()[i]){
      hitType = HitType(rht);
      break;
    }
  }
  TVector3 r3s;
  TVector3 p3s;
  int ind = -1;
  int lay;
  int pdgId;
  int process;
  int bx;
  int evt;
  int simTkIdx = -1;
  HitType hitType = HitType::Unknown;
  bool isBarrel;
  void print(const std::string& pfx){
    std::cout<<pfx<<" "<<ind<<" L"<<lay<<(isBarrel? "b ": "e ")<<evt<<":"<<bx<<" "<<pdgId<<":"<<process<<":"<<simTkIdx
	     <<" ("<<r3s.Pt()<<", "<<r3s.Eta()<<","<<r3s.Phi()<<","<<r3s.Z() <<") "
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

TVector3 linePropagateR(const TVector3& r3, const TVector3& p3, double rDest, int& status, bool useClosest = true, bool verbose = false){
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
    std::cout<<"linePropagateR "<<r3.Pt()<<" "<<r3.Phi()<<" "<<r3.z()<<" "<<pt<<" "<<p
	     <<" "<<d<<" "<<r3.x()*p3.x()<<" "<<r3.y()*p3.y()<<" "<<dotPR2D<<" "<<bb<<" "<<(2.*rt/d + 1.)<<" "<<bb*bb + (2.*rt/d + 1.)
	     <<" => "<<rDest
	     <<" => "<<dest.Pt()<<" "<<dest.Phi()<<" "<<dest.z()
	     <<std::endl;
  }
  return dest;

}

TVector3 linePropagateZ(const TVector3& r3, const TVector3& p3, double zDest, int& status, bool useClosest = true, bool verbose = false){
  double rt = r3.Pt();
  if (useClosest && r3.z()*zDest < 0) zDest = - zDest;
  double d = zDest - r3.z();
  
  double p =  p3.Mag();
  
  TVector3 dest = p3.z() != 0 ? r3 + p3*(d*p/p3.z()/p) : r3;
  if (verbose || std::abs(dest.z() - zDest)>0.001){
    std::cout<<"linePropagateZ "<<r3.Pt()<<" "<<r3.Phi()<<" "<<r3.z()<<" "<<p3.z()<<" "<<p
	     <<" "<<r3.x()*p3.x()<<" "<<r3.y()*p3.y()
	     <<" => "<<zDest
	     <<" => "<<dest.Pt()<<" "<<dest.Phi()<<" "<<dest.z()
	     <<std::endl;
  }
  status = std::abs(dest.z() - zDest)>0.001;
  return dest;

}

std::pair<TVector3,TVector3> helixPropagateApproxR(const TVector3& r3, const TVector3& p3, double rDest, int q, int& status, bool useClosest = true, bool verbose = false){
  double epsilon = 0.001;
  double p = p3.Mag();
  double kap = (2.99792458e-3*3.8*q/p);
  
  auto lastR3 = r3;
  auto lastT3 = p3.Unit();
  int nIts = 7;

  while (std::abs(lastR3.Perp() - rDest) > epsilon && nIts >= 0){
    auto lineEst = linePropagateR(lastR3, lastT3*p, rDest, status, useClosest, verbose);
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

std::pair<TVector3,TVector3> helixPropagateApproxZ(const TVector3& r3, const TVector3& p3, double zDest, int q, int& status, bool useClosest = true, bool verbose = false){
  if (useClosest && r3.z()*zDest < 0) zDest = - zDest;
  
  double epsilon = 0.001;
  double p = p3.Mag();
  double kap = (2.99792458e-3*3.8*q/p);

  auto lastR3 = r3;
  auto lastT3 = p3.Unit();
  int nIts = 7;

  while (std::abs(lastR3.z() - zDest) > epsilon && nIts >= 0){
    auto lineEst = linePropagateZ(lastR3, lastT3*p, zDest, status, useClosest, verbose);
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
      std::cout<<"nIts "<<nIts<<" zDest "<<zDest<<" dS "<<dS<<" phi "<<phi
	       <<" r3In ("<<r3.Pt()<<", "<<r3.Eta()<<", "<<r3.Phi()<<")"
	       <<" p3In ("<<p3.Pt()<<", "<<p3.Eta()<<", "<<p3.Phi()<<")"
	       <<" r3out ("<<lastR3.Pt()<<", "<<lastR3.Eta()<<", "<<lastR3.Phi()<<")"
	       <<" p3Out ("<<lastT3.Pt()*p<<", "<<lastT3.Eta()<<", "<<lastT3.Phi()<<")"
	       <<std::endl;
    }
  }
  status = (std::abs(lastR3.z() - zDest) > epsilon);
  return {lastR3, lastT3*p};
  
}

bool sameLowerLayerID(const SuperDoublet& a, const SuperDoublet& b){
  return (a.mdRef.pixL == b.mdRef.pixL && a.mdOut.pixL == b.mdOut.pixL );
}

int ScanChainMiniDoublets( TChain* chain, int nEvents = -1, bool drawPlots = false, double ptMinGlobal = 1.0) {

  std::vector<double> ptBins {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};

  constexpr int minLayer = 5;
  
  std::array<HistoSet1D, nLayersB+1> layerMD_pt_all;
  createPtHistograms(layerMD_pt_all, "all", ptBins);
  
  std::array<HistoSet1D, nLayersB+1> layerMD_pt_prim_all;
  createPtHistograms(layerMD_pt_prim_all, "prim_all", ptBins);

  std::array<HistoSet1D, nLayersB+1> layerMD_pt_prim_tt;
  createPtHistograms(layerMD_pt_prim_tt, "prim_tt", ptBins);

  std::array<float, nLayersB+1> miniDelta {0, 0, 0, 0, 0,
      0.26, 0.16, 0.16, 0.18, 0.18, 0.18};

  //p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
  //(1 sigma based on 95.45% = 2sigma at 2 GeV)
  std::array<float, nLayersB+1> miniMulsPtScale {0, 0, 0, 0, 0,
      0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};

  //mean of the horizontal layer position in y; treat this as R below
  std::array<float, nLayersB+1> miniRminMean {0, 0, 0, 0, 0, //may want to fill these
      21.8, 34.6, 49.6, 67.4, 87.6, 106.8};

  std::map<int, std::array<float, 4> > moduleBoundaries;
  std::map<int, int > modulePopulation;
  std::array<float, 4> dbound {999,-999,999,-999}; //zmin, zmax, phimin, phimax
  std::array<float, 4>* cbound;
  std::array<float, 4> const* cboundC;
  std::array<float, 4> const* cboundL;
  std::array<float, 4> const* cboundH;

  std::array<int, nLayersB+1> maxHitsBarrelLayer {};
  
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
      TTree *tree = (TTree*)f.Get("trackingNtuple/tree");
      cms2.Init(tree);
      
      //Event Loop
      unsigned int nEvents = tree->GetEntries();
      for( unsigned int event = 0; event < nEvents && nEventsTotal < nEventsChain; ++event) {
	cms2.GetEntry(event);
	++nEventsTotal;

	int iidOld = -1;

	std::array<int, nLayersB+1> hitsBarrelLayer {};
	auto nPh2 = ph2_isBarrel().size();
	for (auto iph2 = 0U; iph2 < nPh2; ++iph2){
	  if (ph2_isBarrel(iph2) == false) continue;
	  int lay = layer(ph2_lay()[iph2], ph2_det()[iph2]);

	  hitsBarrelLayer[lay]++;
	  int iid = ph2_detId()[iph2];
	  if (iidOld != iid){
	    iidOld = iid;
	    if (modulePopulation.find(iid) == modulePopulation.end()){
	      modulePopulation[iid] = dpop;
	      moduleBoundaries[iid] = dbound;
	    }
	    cpop = &modulePopulation[iid];
	    cbound = &moduleBoundaries[iid];
	  }

	  //look at all associated simhits, ignoring double-counting
	  auto const& ph2shV = ph2_simHitIdx()[iph2];
	  auto nPh2sh = ph2shV.size();
	  for (auto iph2sh = 0U; iph2sh < nPh2sh; ++iph2sh){
	    float z = simhit_z()[ph2shV[iph2sh]];
	    if (z==0) continue;
	    float phi = atan2(simhit_y()[ph2shV[iph2sh]], simhit_x()[ph2shV[iph2sh]]);
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
	}

	for (int i = 0; i<= nLayersB; ++i){
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
      std::array<std::set<int>, nLayersB+1> simIdxInLayer;
      int iidStart = -1;
      int iidEnd = -1;
      int iidOld = -1;
      auto nPh2 = ph2_isBarrel().size();
      for (auto iph2 = 0U; iph2 < nPh2; ++iph2){
	if (ph2_isBarrel(iph2) == false) continue;
	int lay = layer(ph2_lay()[iph2], ph2_det()[iph2]);
	if (lay < 5 ) continue;

	int iid = ph2_detId()[iph2];
	if (iidOld != iid){
	  iidStart = iph2;
	  iidOld = iid;
	  cboundC = &moduleBoundaries[iid];
	  if ((iid & miniMask)== lowerId){
	    cboundL = cboundC;
	    cboundH = &moduleBoundaries[iid+upperIdDelta];
	  }
	}
	if ((iid & miniMask)!= lowerId) continue; //only lower module

	TVector3 r3Rec(ph2_x()[iph2], ph2_y()[iph2], ph2_z()[iph2]);

	TVector3 r3Sim, p3Sim;
	int iProcess = 999; int iParticle = 0;
	int ibx = 999; int iev = 999;
	float rs = 0; float ps = 0; float pts = 0;
	int iph2sh = -1;//keep track of this rec-hit simhit
	int iSimIdx = -1;
	auto const& iph2shV = ph2_simHitIdx()[iph2];
	for (auto ph2sh : iph2shV){//look for the first good simhit from a sim-track which was not used yet
	  SimHit sh(ph2sh);
	  float ars = sh.r3s.Pt();
	  float apts = sh.p3s.Pt();
	  float aps = sh.p3s.Mag();
	  if (ars == 0 || aps == 0 ) continue;
	  if (apts < 0.8*ptMinGlobal )  continue;
	  iSimIdx = sh.simTkIdx;
	  
	  r3Sim = sh.r3s;
	  p3Sim = sh.p3s;
	  rs = ars;
	  ps = aps;
	  pts = apts;
	  iProcess = sh.process;
	  iParticle = sh.pdgId;
	  ibx = sh.bx;
	  iev = sh.evt;
	  iph2sh = ph2sh;
	  if (simIdxInLayer[lay].find(iSimIdx) != simIdxInLayer[lay].end()) continue; //only one hit per layer per track
	  else break;
	}
	if (rs == 0 || ps == 0 ) continue;
	if (pts < 0.8*ptMinGlobal )  continue;

	//	std::cout<<__LINE__<<" "<<iph2<<std::endl;
	if (simIdxInLayer[lay].find(iSimIdx) != simIdxInLayer[lay].end()) continue; //redundant check wrt the one above?
	else {
	  simIdxInLayer[lay].insert(iSimIdx);
	}

	MDStats md;

	bool isPrimaryAny = (iProcess == 2 && ibx == 0);
	bool isPrimaryTT = (isPrimaryAny && iev == 0);
	md.pdgId = iParticle;
	md.lowerNTP = iph2shV.size();//FIXME: SHOULD THIS BE RECOMPUTED WITH SOME THRESHOLD?

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
	bbs = std::abs(dotPR2Ds/pts/miniDelta[lay]);
	xxs = (ps/pts)*( sqrt(bbs*bbs + (2.*rs/miniDelta[lay] + 1.)) - bbs);
	TVector3 nextR3SimDes = r3Sim + p3Sim*(dir*miniDelta[lay]/ps)*xxs;
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
	//	std::cout<<__LINE__<<" "<<iph2<<std::endl;

	for (auto jph2 = iidStart; jph2 < static_cast<int>(nPh2); ++jph2){
	  if (jph2 == static_cast<int>(iph2)) continue;
	  int jid = ph2_detId()[jph2];
	  TVector3 ar3(ph2_x()[jph2], ph2_y()[jph2], ph2_z()[jph2]);
	  if (iid == jid){
	    otherR3Rec.emplace_back(ar3);
	    
	    TVector3 ar3s;
	    TVector3 ap3s;
	    auto const& jph2shV = ph2_simHitIdx()[jph2];
	    for (auto jph2sh : jph2shV){
	      if (jph2sh == iph2sh ) continue;
	      SimHit sh(jph2sh);
	      float ars = sh.r3s.Pt();
	      if (ars == 0) continue;

	      ar3s = sh.r3s;
	      ap3s = sh.p3s;
	    
	      double aps = ap3s.Mag();
	      double apts = ap3s.Pt();
	      
	      double dotOPR2Ds = ar3s.x()*ap3s.x() + ar3s.y()*ap3s.y();
	      double odir = dotOPR2Ds > 0 ? 1. : -1.;
	      bbs = std::abs(dotOPR2Ds/apts/0.2);
	      xxs = (aps/apts)*( sqrt(bbs*bbs + (2.*ars/0.2 + 1.)) - bbs);
	      TVector3 nr3s2mm = ar3s + ap3s*(dir*0.2/aps)*xxs;
	      bbs = std::abs(dotOPR2Ds/apts/miniDelta[lay]);
	      xxs = (aps/apts)*( sqrt(bbs*bbs + (2.*ars/miniDelta[lay] + 1.)) - bbs);	    
	      TVector3 nr3sDes = ar3s + ap3s*(dir*miniDelta[lay]/aps)*xxs;
	      
	      otherR3Sim.emplace_back(ar3s);
	      if (aps > 0.001){//line-propagate the others only with p > 1 MeV
		nextOtherR3Sim2mm.emplace_back(nr3s2mm);
		nextOtherR3SimDes.emplace_back(nr3sDes);
	      }
	    }
	  } else if (iid+upperIdDelta == jid && lay>=5) {//upper module in the doublet pairs
	    auto const& jph2shV = ph2_simHitIdx()[jph2];
	    for (auto jph2sh : jph2shV){
	      if (jph2sh == iph2sh ) continue;//can really happen .. assert instead?
	      SimHit sh(jph2sh);
	      float ars = sh.r3s.Pt();
	      if (ars == 0) continue;

	      if (sh.simTkIdx == iSimIdx){
		if (md.upperMatchFull) continue;//match only once
		//pick the first hit consistent in p with iph2sh
		auto jp3s = sh.p3s;
		if (jp3s.Pt() > 0.5*p3Sim.Pt() ){
		  int jProcess = sh.process;
		  int jParticle = sh.pdgId;
		  
		  md.upperMatchToTP = true;
		  if (iSimIdx >= 0 && iProcess == jProcess && iParticle == jParticle){
		    nextR3SimAct = sh.r3s;
		    nextR3SimAct_isValid = true;
		    nextR3Rec = ar3;
		    nextR3Rec_isValid = true;
		    md.upperMatchFull = true;
		    // if (pts > 10 && lay == 5 && isPrimaryTT){
		    //   const float amdDir = (nextR3SimAct-r3Sim).DeltaPhi(r3Sim);
		    //   TVector3 op3(sh.p3s);
		    //   std::cout<<__LINE__<<":"<<nEventsTotal<<" L:"<<lay
		    // 	   <<" ptI:"<<pts<<" ptJ:"<<op3.Pt()
		    // 	   <<" idx:"<<iSimIdx<<" dir:"<<amdDir
		    // 	   <<" procI:"<<iProcess<<" procJ:"<<jProcess
		    // 	   <<" bxI:"<<ibx<<" bxJ:"<<sh.bx
		    // 	   <<" evI:"<<iev<<" evJ:"<<sh.event
		    // 	   <<" typeI:"<<iParticle<<" typeJ:"<<jParticle
		    // 	   <<std::endl;
		    // }
		  }//same process as well
		}//momentum is consistent with p3Sim
	      } else {
		if (ars != 0) nextOtherR3SimAct.emplace_back(sh.r3s);
	      }//not matching to the lower level sim
	    }//end of or (auto jph2sh : jph2shV){

	    if (! md.upperMatchToTP){
	      nextOtherR3Rec.emplace_back(ar3);
	    }
	    
	  } else {
	    //once it's not equal to current or upper(+upperIdDelta), we never expect another entry in the same module
	    break;
	  }//end of module ids to check
	}//loop over the other ph2 hits jph2

	//	std::cout<<__LINE__<<" "<<iph2<<std::endl;
	
	//these two are the same for all layers
	const float ptCut = 1.0;
	const float miniSlope = rs*k2Rinv1GeVf/ptCut;
	
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

	//	std::cout<<__LINE__<<" "<<iph2<<std::endl;

	md.nOthers = otherR3Sim.size();
	md.nOthersRec = otherR3Rec.size();
	fillLayerMD_pt(layerMD_pt_all[lay], pts, md);
	if (isPrimaryAny){
	  fillLayerMD_pt(layerMD_pt_prim_all[lay], pts, md);
	}
	if (isPrimaryTT){
	  fillLayerMD_pt(layerMD_pt_prim_tt[lay], pts, md);
	}
	
	//	std::cout<<__LINE__<<" "<<iph2<<std::endl;
      }//loop over iph2


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
  for (int iL = 1; iL <= nLayersB; ++iL){
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
    for (int iL = 5; iL <= nLayersB; ++iL){
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

    }//iL from 5 to nLayersB

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
  for (int iL = 1; iL <= nLayersB; ++iL){
    layerMD_pt_all[iL].write(outHistograms);
    layerMD_pt_prim_all[iL].write(outHistograms);
    layerMD_pt_prim_tt[iL].write(outHistograms);
  }
  outHistograms->Write();
  outHistograms->Close();
  std::cout<<__LINE__<<" Done "<<std::endl;

  return 0;
}

int ScanChainMockSuperDoublets( TChain* chain, int nEvents = -1, const bool drawPlots = false, const int mockMode = 0,
				const double sdOffsetB = 2.0, const double sdOffsetE = 2.0,
				const int useSeeds = 0, const bool layoutOnly = false, const bool cumulativeCuts = true, const bool addEndcaps = false,
				const bool useFullR3Endcap = false, const bool effForPromptTracks = false){
  
  const float minTPdxy = 0;
  const float maxTPdxy = 1e9;
  const float maxTPdz = 1e9;
  const float ptCutAll = 1.0f;
  const float sinAlphaMax = 0.95f;
  //mockMode:
  //0 for helix to ref and then straight line;
  //1 for helix to all ref layers;
  //3: as 1, for first hit offset to the rechit position

  //useSeeds:
  //0: not used
  //1: convert to tangent SuperDoublet at outer hit

  //cumulativeCuts:
  //false   : computations or analysis of building step is done inclusively (all combinations tried; cuts applied last)
  //true (D): computations or analysis of building step stop after combination fails a cut on a computed value
  
  std::cout<<"Running in mockMode "<<mockMode<<std::endl;
  std::cout<<"Running with SD distance "<<sdOffsetB<<" "<<sdOffsetE<<std::endl;
  std::cout<<"Running with useSeeds "<<useSeeds<<std::endl;
  if (addEndcaps){
    std::cout<<"Endcap is enabled "<<std::endl;
    if (useFullR3Endcap) std::cout<<"Use full R3 in endcap"<<std::endl;
    else std::cout<<"Use dPhiPos in endcap"<<std::endl;
  }
  
  bool debugReco = false;

  std::vector<TH1*> outputHV; outputHV.reserve(1024);
  
  std::vector<double> ptBins {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};
  if (ptCutAll == 0.5) ptBins = {0, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};

  constexpr int minLayer = 5;
  
  std::array<HistoSet1D, nLayersA+1> layerMD_pt_all;
  createPtHistograms(layerMD_pt_all, "all", ptBins);
  
  std::array<HistoSet1D, nLayersA+1> layerMD_pt_prim_all;
  createPtHistograms(layerMD_pt_prim_all, "prim_all", ptBins);

  std::array<HistoSet1D, nLayersA+1> layerMD_pt_prim_tt;
  createPtHistograms(layerMD_pt_prim_tt, "prim_tt", ptBins);

  std::array<float, nLayersA+1> miniDelta {0, 0, 0, 0, 0,
      0.26, 0.16, 0.16, 0.18, 0.18, 0.18,
      0.4, 0.4, 0.4, 0.4, 0.4}; //endcap, reuse the same spacing

  constexpr float deltaZLum = 15.0;
  constexpr float pixelPSZpitch = 0.15;
  constexpr float stripPSZpitch = 2.5;
  constexpr float strip2SZpitch = 5.0;
  constexpr float disks2SMinRadius = 60.f;
  //p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
  //(1 sigma based on 95.45% = 2sigma at 2 GeV)
  std::array<float, nLayersA+1> miniMulsPtScale {0, 0, 0, 0, 0,
      0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034,
      0.006, 0.006, 0.006, 0.006, 0.006}; //inter/extra-polated from L11 and L13 both roughly 0.006 [larger R have smaller value by ~50%]

  //mean of the horizontal layer position in y; treat this as R below
  std::array<float, nLayersA+1> miniRminMean {0, 0, 0, 0, 0, //may want to fill these
      21.8, 34.6, 49.6, 67.4, 87.6, 106.8,
      131.4, 156.2, 185.6, 220.3, 261.5};// use z for endcaps


  enum SDLayers {SDL_L0to5=0, SDL_L0to7, SDL_L0to11, SDL_L5to7, SDL_L7to9, SDL_L5to9,
		 SDL_L5to11,  SDL_L5to13, SDL_L7to11, SDL_L7to13, SDL_L11to13, SDL_LMAX};
  std::array<TH1F*, SDL_LMAX> ha_denSDL_pt;
  std::array<TH1F*, SDL_LMAX> ha_num8MH_pt;
  std::array<TH1F*, SDL_LMAX> ha_num4MD_pt;
  std::array<TH1F*, SDL_LMAX> ha_numInSD_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_numInSD_3of4_any_pt;
  std::array<TH1F*, SDL_LMAX> ha_num2SD_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_num2SD_3of4_any_pt;

  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w01_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w012_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0123_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w01234_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w012345_4of4_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0123456_4of4_pt;//see SDLSelectFlags definitions

  std::array<TH1F*, SDL_LMAX> ha_denSDL_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_num8MH_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_num4MD_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_numInSD_4of4_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_numInSD_3of4_any_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_num2SD_4of4_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_num2SD_3of4_any_loProdXY_pt;

  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0_4of4_loProdXY_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w01_4of4_loProdXY_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w012_4of4_loProdXY_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0123_4of4_loProdXY_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w01234_4of4_loProdXY_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w012345_4of4_loProdXY_pt;//see SDLSelectFlags definitions
  std::array<TH1F*, SDL_LMAX> ha_num2SD_w0123456_4of4_loProdXY_pt;//see SDLSelectFlags definitions

  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_NM1dBeta_all;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_zoom_NM1dBeta_all;

  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_NM1dBeta_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_NM1dBeta_pass;

  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass;

  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_0_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_0_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_0_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_0_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_01_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_01_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_01_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_01_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_012_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_012_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_012_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_012_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_0123_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_0123_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_0123_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_0123_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_01234_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_01234_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_01234_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_01234_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_012345_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_012345_pass;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_012345_all;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_012345_pass;
  
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_NM1dBeta_8MH;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom2_NM1dBeta_8MH;
  
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt0p7to1p0;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt1p0to1p2;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt1p2to1p5;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt1p5to2p0;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt2p0to4p0;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt4p0to7p0;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dBeta_zoom_8MH_pt7p0toInf;

  
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_pass_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaIn_zoom_pass_NM1dBeta_8MH;
  std::array<TH2F*, SDL_LMAX> ha_SDL_dBeta_betaOut_zoom_pass_NM1dBeta_8MH;

  std::array<TH1F*, SDL_LMAX> ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3;
  std::array<TH1F*, SDL_LMAX> ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2;

  std::array<TH1F*, SDL_LMAX> ha_numSDL_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_numSDL_3of4_any_pt;

  std::array<TH1F*, SDL_LMAX> ha_numSDL_4of4_loProdXY_pt;
  std::array<TH1F*, SDL_LMAX> ha_numSDL_3of4_any_loProdXY_pt;

  std::array<TH1F*, SDL_LMAX> ha_SDLreco_all_pt;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_no4of4_pt;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_all_eta;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_4of4_eta;
  std::array<TH1F*, SDL_LMAX> ha_SDLreco_no4of4_eta;

  TH2F* h2_hitsXY_ITrec_OTmockLL = new TH2F("h2_hitsXY_ITrec_OTmockLL", "h2_hitsXY_ITrec_OTmockLL", 1200, 0, 120, 1200, 0, 120);
  TH2F* h2_hitsRZ_ITrec_OTmockLL = new TH2F("h2_hitsRZ_ITrec_OTmockLL", "h2_hitsRZ_ITrec_OTmockLL", 1200, 0, 120, 1200, 0, 120);

  TH2F* h2_hitsRZ_ITrec_OTmockLL_BE = new TH2F("h2_hitsRZ_ITrec_OTmockLL_BE", "h2_hitsRZ_ITrec_OTmockLL_BE", 1200, 0, 270, 1200, 0, 120);

  
  std::array<std::array<int, 2>, SDL_LMAX> layersSDL{{
      {0, 5}, {0, 7}, {0, 11}, {5, 7}, {7, 9}, {5, 9},
      {5, 11}, {5, 13}, {7, 11}, {7, 13}, {11, 13} }};
  std::array<std::array<SDLayers, nLayersA+1>, nLayersA+1> sdlFromLayers;
  for (int i = 0U; i<= nLayersA; ++i){
    for (int j = 0U; j<= nLayersA; ++j){
      sdlFromLayers[i][j] = SDL_LMAX;
      for (auto k = 0U; k < SDL_LMAX; ++k){
	if (i == layersSDL[k][0] && j == layersSDL[k][1]){
	  sdlFromLayers[i][j] = (SDLayers)k;
	  break;
	}
      }
    }
  }
  
  for (int i = 0; i< SDL_LMAX; ++i){
    auto iMin = layersSDL[i][0];
    auto iMax = layersSDL[i][1];
    std::string hn = Form("h_denSDL_%dto%d_pt", iMin, iMax);
    ha_denSDL_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_denSDL_pt[i]);
    
    hn = Form("h_num8MH_%dto%d_pt", iMin, iMax);
    ha_num8MH_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num8MH_pt[i]);
    
    hn = Form("h_num4MD_%dto%d_pt", iMin, iMax);
    ha_num4MD_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num4MD_pt[i]);
    
    hn = Form("h_numInSD_4of4_%dto%d_pt", iMin, iMax);
    ha_numInSD_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_numInSD_4of4_pt[i]);    
    hn = Form("h_numInSD_3of4_any_%dto%d_pt", iMin, iMax);
    ha_numInSD_3of4_any_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_numInSD_3of4_any_pt[i]);    

    hn = Form("h_num2SD_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_4of4_pt[i]);    
    hn = Form("h_num2SD_3of4_any_%dto%d_pt", iMin, iMax);
    ha_num2SD_3of4_any_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_3of4_any_pt[i]);    

    hn = Form("h_num2SD_w0_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w0_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w0_4of4_pt[i]);    
    hn = Form("h_num2SD_w01_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w01_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w01_4of4_pt[i]);    
    hn = Form("h_num2SD_w012_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w012_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w012_4of4_pt[i]);    
    hn = Form("h_num2SD_w0123_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w0123_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w0123_4of4_pt[i]);    
    hn = Form("h_num2SD_w01234_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w01234_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w01234_4of4_pt[i]);    
    hn = Form("h_num2SD_w012345_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w012345_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_w012345_4of4_pt[i]);    
    hn = Form("h_num2SD_w0123456_4of4_%dto%d_pt", iMin, iMax);
    ha_num2SD_w0123456_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_w0123456_4of4_pt[i]);    

    //same set with low prodXY 
    hn = Form("h_denSDL_%dto%d_loProdXY_pt", iMin, iMax);
    ha_denSDL_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_denSDL_loProdXY_pt[i]);
    
    hn = Form("h_num8MH_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num8MH_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num8MH_loProdXY_pt[i]);
    
    hn = Form("h_num4MD_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num4MD_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num4MD_loProdXY_pt[i]);
    
    hn = Form("h_numInSD_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_numInSD_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_numInSD_4of4_loProdXY_pt[i]);    
    hn = Form("h_numInSD_3of4_any_%dto%d_loProdXY_pt", iMin, iMax);
    ha_numInSD_3of4_any_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_numInSD_3of4_any_loProdXY_pt[i]);    

    hn = Form("h_num2SD_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_3of4_any_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_3of4_any_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_3of4_any_loProdXY_pt[i]);    

    hn = Form("h_num2SD_w0_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w0_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w0_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_w01_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w01_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w01_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_w012_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w012_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w012_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_w0123_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w0123_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w0123_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_w01234_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w01234_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_num2SD_w01234_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_w012345_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w012345_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_w012345_4of4_loProdXY_pt[i]);    
    hn = Form("h_num2SD_w0123456_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_num2SD_w0123456_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_num2SD_w0123456_4of4_loProdXY_pt[i]);    

    //plots without MC truth requirement
    hn = Form("h2_SDL_dBeta_betaIn_NM1dBeta_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_NM1dBeta_all[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -1, 1, 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_betaIn_NM1dBeta_all[i]);    
    hn = Form("h2_SDL_dBeta_betaIn_zoom_NM1dBeta_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_zoom_NM1dBeta_all[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_betaIn_zoom_NM1dBeta_all[i]);    

    hn = Form("h_SDL_dBeta_NM1dBeta_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_NM1dBeta_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_NM1dBeta_all[i]);    
    hn = Form("h_SDL_dBeta_NM1dBeta_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_NM1dBeta_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_NM1dBeta_pass[i]);    

    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_all_NGLL_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2[i]);    

    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass[i]);    

    //looser selections
    hn = Form("h_SDL_dBeta_0_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_0_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_0_all[i]);    
    hn = Form("h_SDL_dBeta_0_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_0_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_0_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_0_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_0_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_0_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_0_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_0_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_0_pass[i]);    

    //01
    hn = Form("h_SDL_dBeta_01_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_01_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_01_all[i]);    
    hn = Form("h_SDL_dBeta_01_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_01_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_01_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_01_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_01_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_01_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_01_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_01_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_01_pass[i]);    
    //012
    hn = Form("h_SDL_dBeta_012_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_012_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_012_all[i]);    
    hn = Form("h_SDL_dBeta_012_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_012_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_012_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_012_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_012_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_012_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_012_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_012_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_012_pass[i]);    
    //0123
    hn = Form("h_SDL_dBeta_0123_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_0123_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_0123_all[i]);    
    hn = Form("h_SDL_dBeta_0123_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_0123_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_0123_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_0123_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_0123_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_0123_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_0123_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_0123_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_0123_pass[i]);    
    //01234
    hn = Form("h_SDL_dBeta_01234_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_01234_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_01234_all[i]);    
    hn = Form("h_SDL_dBeta_01234_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_01234_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_01234_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_01234_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_01234_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_01234_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_01234_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_01234_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_01234_pass[i]);    
    //012345
    hn = Form("h_SDL_dBeta_012345_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_012345_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_012345_all[i]);    
    hn = Form("h_SDL_dBeta_012345_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_012345_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_012345_pass[i]);    
    hn = Form("h_SDL_dBeta_zoom_012345_all_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_012345_all[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_012345_all[i]);    
    hn = Form("h_SDL_dBeta_zoom_012345_pass_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_012345_pass[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dBeta_zoom_012345_pass[i]);    


    //plots with MC truth requirement
    hn = Form("h_SDL_dBeta_zoom_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom_NM1dBeta_8MH[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_NM1dBeta_8MH[i]);    
    hn = Form("h_SDL_dBeta_zoom2_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_zoom2_NM1dBeta_8MH[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.02, 0.02);
    outputHV.push_back(ha_SDL_dBeta_zoom2_NM1dBeta_8MH[i]);    

    //pt slices, mostly unbiased
    hn = Form("h_SDL_dBeta_zoom_8MH_pt0p7to1p0_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt0p7to1p0[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt0p7to1p0[i]);    
    hn = Form("h_SDL_dBeta_zoom_8MH_pt1p0to1p2_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt1p0to1p2[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt1p0to1p2[i]);    
    hn = Form("h_SDL_dBeta_zoom_8MH_pt1p2to1p5_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt1p2to1p5[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt1p2to1p5[i]);    
    hn = Form("h_SDL_dBeta_zoom_8MH_pt1p5to2p0_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt1p5to2p0[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt1p5to2p0[i]);    
    hn = Form("h_SDL_dBeta_zoom_8MH_pt2p0to4p0_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt2p0to4p0[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt2p0to4p0[i]);    
    hn = Form("h_SDL_dBeta_zoom_8MH_pt4p0to7p0_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt4p0to7p0[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt4p0to7p0[i]);    
    hn = Form("h_SDL_dBeta_zoom_8MH_pt7p0toInf_%dto%d", iMin, iMax);
    ha_SDL_dBeta_zoom_8MH_pt7p0toInf[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_zoom_8MH_pt7p0toInf[i]);    

    //2D
    hn = Form("h2_SDL_dBeta_betaIn_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -1, 1, 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_betaIn_NM1dBeta_8MH[i]);    

    hn = Form("h2_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH[i]);    

    hn = Form("h2_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH[i]);    

    hn = Form("h2_SDL_dBeta_betaIn_pass_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_pass_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -1, 1, 400, -0.5, 0.5);
    outputHV.push_back(ha_SDL_dBeta_betaIn_pass_NM1dBeta_8MH[i]);    

    hn = Form("h2_SDL_dBeta_betaIn_zoom_pass_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaIn_zoom_pass_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_betaIn_zoom_pass_NM1dBeta_8MH[i]);    

    hn = Form("h2_SDL_dBeta_betaOut_zoom_pass_NM1dBeta_8MH_%dto%d_pt", iMin, iMax);
    ha_SDL_dBeta_betaOut_zoom_pass_NM1dBeta_8MH[i] = new TH2F(hn.c_str(), hn.c_str(), 400, -0.35, 0.35, 400, -0.06, 0.06);
    outputHV.push_back(ha_SDL_dBeta_betaOut_zoom_pass_NM1dBeta_8MH[i]);    

    //dZeta : angle in r-z
    hn = Form("h_SDL_dZeta_zoom_NM1dBeta_all_NGLL_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL[i]);    
    hn = Form("h_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4[i]);    
    hn = Form("h_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3[i]);    
    hn = Form("h_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2_%dto%d_pt", iMin, iMax); //not-ghost, based on lower-layer ID
    ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2[i] = new TH1F(hn.c_str(), hn.c_str(), 400, -0.15, 0.15);
    outputHV.push_back(ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2[i]);    

    hn = Form("h_numSDL_4of4_%dto%d_pt", iMin, iMax);
    ha_numSDL_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_numSDL_4of4_pt[i]);    
    hn = Form("h_numSDL_3of4_any_%dto%d_pt", iMin, iMax);
    ha_numSDL_3of4_any_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_numSDL_3of4_any_pt[i]);    

    hn = Form("h_numSDL_4of4_%dto%d_loProdXY_pt", iMin, iMax);
    ha_numSDL_4of4_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());    
    outputHV.push_back(ha_numSDL_4of4_loProdXY_pt[i]);    
    hn = Form("h_numSDL_3of4_any_%dto%d_loProdXY_pt", iMin, iMax);
    ha_numSDL_3of4_any_loProdXY_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_numSDL_3of4_any_loProdXY_pt[i]);    

    hn = Form("h_SDLreco_all_%dto%d_pt", iMin, iMax);
    ha_SDLreco_all_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_SDLreco_all_pt[i]);    
    hn = Form("h_SDLreco_4of4_%dto%d_pt", iMin, iMax);
    ha_SDLreco_4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_SDLreco_4of4_pt[i]);    
    hn = Form("h_SDLreco_no4of4_%dto%d_pt", iMin, iMax);
    ha_SDLreco_no4of4_pt[i] = new TH1F(hn.c_str(), hn.c_str(), ptBins.size()-1, ptBins.data());
    outputHV.push_back(ha_SDLreco_no4of4_pt[i]);    

    hn = Form("h_SDLreco_all_%dto%d_eta", iMin, iMax);
    ha_SDLreco_all_eta[i] = new TH1F(hn.c_str(), hn.c_str(), 40, -2.0, 2.0);
    outputHV.push_back(ha_SDLreco_all_eta[i]);    
    hn = Form("h_SDLreco_4of4_%dto%d_eta", iMin, iMax);
    ha_SDLreco_4of4_eta[i] = new TH1F(hn.c_str(), hn.c_str(), 40, -2.0, 2.0);
    outputHV.push_back(ha_SDLreco_4of4_eta[i]);    
    hn = Form("h_SDLreco_no4of4_%dto%d_eta", iMin, iMax);
    ha_SDLreco_no4of4_eta[i] = new TH1F(hn.c_str(), hn.c_str(), 40, -2.0, 2.0);
    outputHV.push_back(ha_SDLreco_no4of4_eta[i]);    
  }

  enum TimerTypes {T_timeLayout=0, T_timeReco, T_timeValidation,
		   T_timeVal_SHLoad, T_timeVal_MHMDSDMatch, T_timeVal_SDL,  T_timeVal_SDLMatch, T_N};
  std::array<TStopwatch, T_N> timerA {};
  std::array<std::string, T_N> timerNameA {"timeLayout", "timeReco", "timeValidation",
      "timeVal_SHLoad", "timeVal_MHMDSDMatch", "timeVal_SDL", "timeVal_SDLMatch"};
  bool runDetailedTimers = false;


  std::map<int, std::array<float, 4> > moduleBoundaries;
  std::map<int, int > modulePopulation;
  std::array<float, 4> dbound {999,-999,999,-999}; //zmin, zmax, phimin, phimax
  std::array<float, 4>* cbound;
  std::array<float, 4> const* cboundC;
  std::array<float, 4> const* cboundL;
  std::array<float, 4> const* cboundH;
  
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
      TTree *tree = (TTree*)f.Get("trackingNtuple/tree");
      cms2.Init(tree);
      
      //Event Loop
      unsigned int nEvents = tree->GetEntries();
      for( unsigned int event = 0; event < nEvents*10 && nEventsTotal < nEventsChain*10; ++event) {
	cms2.GetEntry(event);
	++nEventsTotal;

	int iidOld = -1;

	std::array<int, nLayersB+1> hitsBarrelLayer {};
	auto nPh2 = ph2_isBarrel().size();
	for (auto iph2 = 0U; iph2 < nPh2; ++iph2){
	  bool isBarrel = ph2_isBarrel(iph2);
	  int lay = layer(ph2_lay()[iph2], ph2_det()[iph2]);

	  hitsBarrelLayer[lay]++;
	  int iid = ph2_detId()[iph2];
	  if (iidOld != iid){
	    iidOld = iid;
	    if (modulePopulation.find(iid) == modulePopulation.end()){
	      modulePopulation[iid] = dpop;
	      moduleBoundaries[iid] = dbound;
	    }
	    cpop = &modulePopulation[iid];
	    cbound = &moduleBoundaries[iid];
	  }

	  //look at all associated simhits, ignoring double-counting
	  auto const& ph2shV = ph2_simHitIdx()[iph2];
	  auto nPh2sh = ph2shV.size();
	  for (auto iph2sh = 0U; iph2sh < nPh2sh; ++iph2sh){
	    const float x = simhit_x()[ph2shV[iph2sh]];
	    const float y = simhit_y()[ph2shV[iph2sh]];
	    const float z = isBarrel ? simhit_z()[ph2shV[iph2sh]] : sqrt(x*x + y*y);
	    if (z==0) continue;
	    float phi = atan2(y, x);
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
	  }//simhits for a given ph2
	}//for (auto iph2 = 0U; iph2 < nPh2; ++iph2){

      }//for( unsigned int event = 0
    }//files in chain
  }//geom range map loop scope
  
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
    TTree *tree = (TTree*)f.Get("trackingNtuple/tree");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents && nEventsTotal < nEventsChain; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

      std::cout<<"Event "<<nEventsTotal<<std::endl;
      //this could be filled on-demand and in some regions of interest at some point
      //assume somewhat hermetic tracker and the only hits in the barrel to be from the barrel
      //reference layer minidoublet
      std::array<std::vector<std::pair<int, const V3WithCache> >, nLayersA+1 > mockLayerMDfwRefLower;
      std::array<std::vector<std::pair<int, const V3WithCache> >, nLayersA+1 > mockLayerMDfwRefUpper;
      //offset layer minidoublet (Ncm hardcode for now)
      std::array<std::vector<std::pair<int, const V3WithCache> >, nLayersA+1 > mockLayerMDfwDNcmLower;
      std::array<std::vector<std::pair<int, const V3WithCache> >, nLayersA+1 > mockLayerMDfwDNcmUpper;
      for (int iL = 1; iL <=nLayersA; ++iL){
	mockLayerMDfwRefLower[iL].reserve(maxHitsInLayer);
	mockLayerMDfwRefUpper[iL].reserve(maxHitsInLayer);
	mockLayerMDfwDNcmLower[iL].reserve(maxHitsInLayer);
	mockLayerMDfwDNcmUpper[iL].reserve(maxHitsInLayer);
      }

      auto const nSimhit = simhit_lay().size();
      assert(nSimhit < HitIndexWithType::indexMask);
      auto const nSim = sim_nPixel().size();

      // extract TP index per simhit .. not very useful
      std::vector<int> simsPerSimHit(nSimhit,-1);
      std::vector<int> simsPerSimHitAll(nSimhit,-1);

      int errCount = 0;
      for (auto i = 0U; i<nSim; ++i){
	TVector3 p3(sim_px()[i], sim_py()[i], sim_pz()[i]);
	float tpPt = p3.Pt();

	auto const& simHitIdxV = sim_simHitIdx()[i];
	for (auto ishIdx : simHitIdxV){	  
	  SimHit simH(ishIdx);

	  if (simH.p3s.Pt()>=0.8*tpPt){
	    int& iSim = simsPerSimHit[ishIdx];
	    
	    if (iSim == -1 ) iSim = i;
	    else if (iSim != int(i)){
	      std::cout<<"ERROR: repeated hit-tp "<<ishIdx<<" has tp "<<iSim<<" and "<<i<<std::endl;
	      errCount++;
	    }
	  }
	  int& iSimAll = simsPerSimHitAll[ishIdx];
	  
	  if (iSimAll == -1 ) iSimAll = i;
	  else if (iSimAll != int(i)){
	    std::cout<<"ERROR: repeated hit-tp "<<ishIdx<<" has tp "<<iSimAll<<" and "<<i<<std::endl;
	    errCount++;
	  }
	}//simhits for a given TP
	if (errCount > 1000) return 1;
      }//TPs
      
      auto const nPh2 = ph2_lay().size();
      std::vector<int> simsPerPh2Hit(nPh2,-1);
      std::vector<int> simsPerPh2HitAll(nPh2,-1);
      std::vector<int> simHitsPerPh2Hit(nPh2,-1);
      std::vector<int> simHitsPerPh2HitAll(nPh2,-1);
      for (auto iPh2 = 0U; iPh2 < nPh2; ++iPh2){
	auto const& r2sV = ph2_simHitIdx()[iPh2];//reco->sim vector
	if (r2sV.empty()) continue;
	SimHit bestSH;
	for (auto r2s : r2sV){
	  SimHit r2sSH(r2s);
	  //FIXME: check if it makes sense to pick also by charge fraction
	  if (r2sSH.p3s.Mag() > bestSH.p3s.Mag()) bestSH = r2sSH;
	}
	if (bestSH.ind >= 0){
	  simsPerPh2HitAll[iPh2] = bestSH.simTkIdx;
	  simHitsPerPh2HitAll[iPh2] = bestSH.ind;
	  if (simsPerSimHit[bestSH.ind] >= 0){//well-matching to TP
	    simsPerPh2Hit[iPh2] = bestSH.simTkIdx;
	    simHitsPerPh2Hit[iPh2] = bestSH.ind;
	  }
	}
      }

      auto const nPix = pix_lay().size();
      std::vector<int> simsPerPixHit(nPix,-1);
      std::vector<int> simsPerPixHitAll(nPix,-1);
      std::vector<int> simHitsPerPixHit(nPix,-1);
      std::vector<int> simHitsPerPixHitAll(nPix,-1);
      for (auto iPix = 0U; iPix < nPix; ++iPix){
	auto const& r2sV = pix_simHitIdx()[iPix];//reco->sim vector
	if (r2sV.empty()) continue;
	SimHit bestSH;
	for (auto r2s : r2sV){
	  SimHit r2sSH(r2s);
	  //FIXME: check if it makes sense to pick also by charge fraction
	  if (r2sSH.p3s.Mag() > bestSH.p3s.Mag()) bestSH = r2sSH;
	}
	if (bestSH.ind >= 0){
	  simsPerPixHitAll[iPix] = bestSH.simTkIdx;
	  simHitsPerPixHitAll[iPix] = bestSH.ind;
	  if (simsPerSimHit[bestSH.ind] >= 0){//well-matching to TP
	    simsPerPixHit[iPix] = bestSH.simTkIdx;
	    simHitsPerPixHit[iPix] = bestSH.ind;
	  }
	}
      }

      auto simsPerHit = [&](const int iwt){
	HitType htype = HitIndexWithType::type(iwt);
	int index = HitIndexWithType::index(iwt);
	switch (htype){
	case HitType::Pixel: return simsPerPixHit[index]; break;
	case HitType::Phase2OT: return simsPerPh2Hit[index]; break;
	default: assert(0);//should not happen
	}
	return -1;
      };
      auto simsPerHitAll = [&](const int iwt){
	HitType htype = HitIndexWithType::type(iwt);
	int index = HitIndexWithType::index(iwt);
	switch (htype){
	case HitType::Pixel: return simsPerPixHitAll[index]; break;
	case HitType::Phase2OT: return simsPerPh2HitAll[index]; break;
	default: assert(0);//should not happen
	}
	return -1;
      };

      auto simHitsPerHit = [&](const int iwt){
	HitType htype = HitIndexWithType::type(iwt);
	int index = HitIndexWithType::index(iwt);
	switch (htype){
	case HitType::Pixel: return simHitsPerPixHit[index]; break;
	case HitType::Phase2OT: return simHitsPerPh2Hit[index]; break;
	default: assert(0);//should not happen
	}
	return -1;
      };
      auto simHitsPerHitAll = [&](const int iwt){
	HitType htype = HitIndexWithType::type(iwt);
	int index = HitIndexWithType::index(iwt);
	switch (htype){
	case HitType::Pixel: return simHitsPerPixHitAll[index]; break;
	case HitType::Phase2OT: return simHitsPerPh2HitAll[index]; break;
	default: assert(0);//should not happen
	}
	return -1;
      };

      std::array<int, nLayersA+1> nHitsLayer1GeV {};
      std::array<int, nLayersA+1> nHitsLayer2GeV {};

      timerA[T_timeLayout].Start(kFALSE);
      std::cout<<"Load hits "<<std::endl;
      std::array<std::set<int>, nLayersA+1> simIdxDeltaInLayer;
      std::array<std::set<int>, nLayersA+1> simIdxPositronLowPInLayer;
      std::array<std::set<int>, nLayersA+1> simIdxInLayer;
      int iidStart = -1;
      int iidEnd = -1;
      int iidOld = -1;
      for (auto iph2 = 0U; iph2 < nPh2; ++iph2){
	bool isBarrel = ph2_isBarrel(iph2);
	if (!addEndcaps &&  !isBarrel) continue;
	int lay = layer(ph2_lay()[iph2], ph2_det()[iph2]);
	if (addEndcaps && !isBarrel && (lay < 11 || lay > nLayersA) ) continue;

	int iid = ph2_detId()[iph2];
	if (iidOld != iid){
	  iidStart = iph2;
	  iidOld = iid;
	  cboundC = &moduleBoundaries[iid];
	  if ((iid & miniMask)== lowerId){
	    cboundL = cboundC;
	    cboundH = &moduleBoundaries[iid+upperIdDelta];
	  }
	}

	//Too restrictive for reverse TP matching??//	if ((iid & miniMask)!= lowerId) continue; 
		
	TVector3 r3Rec(ph2_x()[iph2], ph2_y()[iph2], ph2_z()[iph2]);

	auto iSHAll = simHitsPerPh2HitAll[iph2];
	auto iSH = simHitsPerPh2Hit[iph2];
	if (iSHAll == -1) continue; //use only hits with sim info
	SimHit simHit(iSHAll);
	
	TVector3 r3Sim(simHit.r3s);
	float rs = r3Sim.Pt();

	TVector3 p3Sim(simHit.p3s);
	float pts = p3Sim.Pt();
	float ps = p3Sim.Mag();

	int iParticle = simHit.pdgId;

	int iSimIdx = simHit.simTkIdx; //tp index in ntuple (not a g4 trackId)

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
	
	const float mdOffset = miniDelta[lay];
	const float rNominal = miniRminMean[lay];

	const float rRefLower = rNominal;
	const float rRefUpper = rNominal+mdOffset;
	const float rSDfwLower = rRefLower + (isBarrel ? sdOffsetB : sdOffsetE);
	const float rSDfwUpper = rRefUpper + (isBarrel ? sdOffsetB : sdOffsetE);

	int q = 0;
	if (iParticle == -11 || iParticle == -13 || iParticle == 211 || iParticle == 321 || iParticle == 2212
	    || iParticle == -3112 || iParticle == 3222 || iParticle == -3312 ) q = 1;
	else if (iParticle == 11 || iParticle == 13 || iParticle == -211 || iParticle == -321 || iParticle == -2212
		 || iParticle == 3112 || iParticle == -3222 || iParticle == 3312 ) q = -1;
	int pstat = 0;
	TVector3 r3RefLower;
	TVector3 p3RefLower;
	if (q == 0){
	  r3RefLower = isBarrel ? linePropagateR(r3Sim, p3Sim, rRefLower, pstat)
	    : linePropagateZ(r3Sim, p3Sim, rRefLower, pstat);
	  p3RefLower = p3Sim;
	} else {
	  auto resProp = isBarrel ? helixPropagateApproxR(r3Sim, p3Sim, rRefLower, q, pstat)
	    :  helixPropagateApproxZ(r3Sim, p3Sim, rRefLower, q, pstat);
	  r3RefLower = resProp.first;
	  p3RefLower = resProp.second;
	  
	}
	
	auto propagateMH = [&](float rDest){
	  if (mockMode == 0) return isBarrel ? linePropagateR(r3RefLower, p3RefLower, rDest, pstat)
			       : linePropagateZ(r3RefLower, p3RefLower, rDest, pstat) ;
	  else if (mockMode == 1 || mockMode == 3) return isBarrel ? helixPropagateApproxR(r3RefLower, p3RefLower, rDest, q, pstat).first
						     : helixPropagateApproxZ(r3RefLower, p3RefLower, rDest, q, pstat).first;
	  else {pstat = 99; return TVector3();}
	};

	auto r3RefLowerMock = r3RefLower;

	if (mockMode == 3){
	  r3RefLowerMock += (r3Rec - r3Sim);
	  if ((iid & miniMask)!= lowerId){//there was no simhit on the pixel layer; make up z or r by roundoff of sim
	    //this mitigates the issue coming from recovery of inefficiency in simHit-rec matching done by using outer mini-layers as well
	    //if the outer mini-layer is used in place of the inner one, the rec-sim shift is too large in the coarse direction
	    //2S layers are not affected because both mini-layers have the same segmentation
	    if (lay >=5 && lay <=7){	    
	      r3RefLowerMock.SetZ(r3RefLowerMock.z() - r3Rec.z() + r3Sim.z());// undo the rec-sim shift in the coarse (bad) direction
	      float binnedZshift = std::round(r3Sim.z()/pixelPSZpitch)*pixelPSZpitch - r3Sim.z();
	      r3RefLowerMock.SetZ(r3RefLowerMock.z() + binnedZshift);
	    } else if (lay >= 11 && r3Rec.Pt() < disks2SMinRadius){
	      const float phiAv = 0.5f*((*cboundC)[2] + (*cboundC)[3]);
	      const float sinPhiAv = sin(phiAv);
	      const float cosPhiAv = cos(phiAv);
	      const float xLocRecMock = r3RefLowerMock.y()*cosPhiAv - r3RefLowerMock.x()*sinPhiAv;//keep this fixed

	      const float yLocRecMock = std::round((r3RefLower.y()*sinPhiAv + r3RefLower.x()*cosPhiAv)/pixelPSZpitch)*pixelPSZpitch;//round-off the sim
	      
	      r3RefLowerMock.SetX(yLocRecMock*cosPhiAv -  xLocRecMock*sinPhiAv);
	      r3RefLowerMock.SetY(yLocRecMock*sinPhiAv +  xLocRecMock*cosPhiAv);

	    }
	  }
	}
	V3WithCache v3RefLowerMock(r3RefLowerMock);

	//set the strip edges
	if (mockMode == 3 && lay >= 11 && r3Rec.Pt() > disks2SMinRadius){
	  const float phiAv = 0.5f*((*cboundC)[2] + (*cboundC)[3]);
	  const float sinPhiAv = sin(phiAv);
	  const float cosPhiAv = cos(phiAv);
	  const float xLocRecMock = r3RefLowerMock.y()*cosPhiAv - r3RefLowerMock.x()*sinPhiAv;//keep this fixed

	  //this is supposedly the middle of the strip (not necessarily the case in SLHC setup)
	  const float yLocRecMock = r3RefLowerMock.y()*sinPhiAv + r3RefLowerMock.x()*cosPhiAv;
	  const float yLocRecMockIn = yLocRecMock - 0.5f*strip2SZpitch;
	  const float yLocRecMockOut = yLocRecMock + 0.5f*strip2SZpitch;
	  
	  TVector2  r2In(yLocRecMockIn*cosPhiAv - xLocRecMock*sinPhiAv, yLocRecMockIn*sinPhiAv + xLocRecMock*cosPhiAv);
	  TVector2  r2Out(yLocRecMockOut*cosPhiAv - xLocRecMock*sinPhiAv, yLocRecMockOut*sinPhiAv + xLocRecMock*cosPhiAv);

	  v3RefLowerMock.rtRHin = r2In.Mod();
	  v3RefLowerMock.phiRHin = r2In.Phi();
	  v3RefLowerMock.rtRHout = r2Out.Mod();
	  v3RefLowerMock.phiRHout = r2Out.Phi();
	}
	
	if (pstat == 0) mockLayerMDfwRefLower[lay].push_back(std::make_pair(HitIndexWithType(iph2, HitType::Phase2OT).indexWithType,
									    v3RefLowerMock));

	auto r3RefUpper = propagateMH(rRefUpper);
	if (mockMode == 3 && lay >= 11 && r3Rec.Pt() > disks2SMinRadius){
	  //place r3RefUpper on the same local y as the r3RefLowerMock to emulate geometry more appropriately
	  const float phiAv = 0.5f*((*cboundC)[2] + (*cboundC)[3]);
	  const float sinPhiAv = sin(phiAv);
	  const float cosPhiAv = cos(phiAv);
	  
	  const float yLocRecLower = r3RefLowerMock.y()*sinPhiAv + r3RefLowerMock.x()*cosPhiAv;
	  const float xLocRecUpper = r3RefUpper.y()*cosPhiAv - r3RefUpper.x()*sinPhiAv;
	  
	  r3RefUpper.SetX(yLocRecLower*cosPhiAv -  xLocRecUpper*sinPhiAv);
	  r3RefUpper.SetY(yLocRecLower*sinPhiAv +  xLocRecUpper*cosPhiAv);
	  
	}
	//FIXME: put more appropriate limits
	const bool r3RefUpperIsGood = (isBarrel && std::abs(r3RefUpper.z()) < 120.f) || (!isBarrel && r3RefUpper.Pt() > 23.f && r3RefUpper.Pt() < 110.f);
	if (pstat == 0 && r3RefUpperIsGood) mockLayerMDfwRefUpper[lay].push_back(std::make_pair(HitIndexWithType(iph2, HitType::Phase2OT).indexWithType,
												V3WithCache(r3RefUpper)));

	//keep the mock outer doublet at its SIM state
	auto r3SDfwLower = propagateMH(rSDfwLower);
	const bool r3SDfwLowerIsGood = (isBarrel && std::abs(r3SDfwLower.z()) < 120.f) || (!isBarrel && r3SDfwLower.Pt() > 23.f && r3SDfwLower.Pt() < 110.f);
	if (pstat == 0 && r3SDfwLowerIsGood) mockLayerMDfwDNcmLower[lay].push_back(std::make_pair(HitIndexWithType(iph2, HitType::Phase2OT).indexWithType,
												  V3WithCache(r3SDfwLower)));
	auto r3SDfwUpper = propagateMH(rSDfwUpper);
	const bool r3SDfwUpperIsGood = (isBarrel && std::abs(r3SDfwUpper.z()) < 120.f) || (!isBarrel && r3SDfwUpper.Pt() > 23.f && r3SDfwUpper.Pt() < 110.f);
	if (pstat == 0 && r3SDfwUpperIsGood)
	  mockLayerMDfwDNcmUpper[lay].push_back(std::make_pair(HitIndexWithType(iph2, HitType::Phase2OT).indexWithType,
							       V3WithCache(r3SDfwUpper)));

	if (pstat == 0 && q != 0 && pts>0.8){
	  if (lay == 5 || lay == 7 || lay == 9){
	    h2_hitsXY_ITrec_OTmockLL->Fill(r3RefLowerMock.X(), r3RefLowerMock.Y());
	    if (r3SDfwLowerIsGood) h2_hitsXY_ITrec_OTmockLL->Fill(r3SDfwLower.X(), r3SDfwLower.Y());
	    h2_hitsRZ_ITrec_OTmockLL->Fill(std::abs(r3RefLowerMock.Z()), r3RefLowerMock.Pt());
	    if (r3SDfwLowerIsGood) h2_hitsRZ_ITrec_OTmockLL->Fill(std::abs(r3SDfwLower.Z()), r3SDfwLower.Pt());
	  }
	  if (lay == 5 || lay == 7 || lay == 9 || lay == 11 || lay == 13 || lay == 15){
	    h2_hitsRZ_ITrec_OTmockLL_BE->Fill(std::abs(r3RefLowerMock.Z()), r3RefLowerMock.Pt());
	    if (lay != 15 && r3SDfwLowerIsGood) h2_hitsRZ_ITrec_OTmockLL_BE->Fill(std::abs(r3SDfwLower.Z()), r3SDfwLower.Pt());
	  }
	}
      }//nPh2: filling mock hits
      std::cout<<"MH stat: nPh2 "<<nPh2<<std::endl;
      for (int iL = 0; iL < nLayersA+1; ++iL){
	std::cout<<" L "<<iL
		 <<" RfL "<<mockLayerMDfwRefLower[iL].size()
		 <<" RfL "<<mockLayerMDfwRefUpper[iL].size()
		 <<" RfL "<<mockLayerMDfwDNcmLower[iL].size()
		 <<" RfL "<<mockLayerMDfwDNcmUpper[iL].size()
		 <<std::endl;	  
      }

      for (auto iPix = 0U; iPix < nPix; ++iPix){
	int lay = layer(pix_lay()[iPix], pix_det()[iPix]);
	bool isBarrel = pix_isBarrel()[iPix];
	TVector3 r3Rec(pix_x()[iPix], pix_y()[iPix], pix_z()[iPix]);
	if (lay < 5 && isBarrel){//fill pixel barrel plot	  
	  h2_hitsXY_ITrec_OTmockLL->Fill(r3Rec.X(), r3Rec.Y());
	  h2_hitsRZ_ITrec_OTmockLL->Fill(std::abs(r3Rec.Z()), r3Rec.Pt());
	  h2_hitsRZ_ITrec_OTmockLL_BE->Fill(std::abs(r3Rec.Z()), r3Rec.Pt());
	}
      }

      timerA[T_timeLayout].Stop();

      if (layoutOnly ) continue;

      timerA[T_timeReco].Start(kFALSE);
      constexpr int maxMDexpected = 100;
      std::array<std::vector<MiniDoublet>, nLayersA+1> mockLayerMDfwRef;
      for (auto& m : mockLayerMDfwRef ) m.reserve(100000);
      std::array<std::vector<MiniDoublet>, nLayersA+1> mockLayerMDfwDNcm;
      for (auto& m : mockLayerMDfwDNcm ) m.reserve(100000);      

      std::array<std::vector<SuperDoublet>, nLayersA+1> mockLayerSDfwDNcm;
      for (auto& m : mockLayerSDfwDNcm ) m.reserve(1000000);      

      std::array<std::vector<bool>, nLayersA+1> mockLayerSDfwDNcm_isSecondaryGhost;
      for (auto& m : mockLayerSDfwDNcm_isSecondaryGhost ) m.reserve(1000000);  
      
      std::array<std::vector<SDLink>, SDL_LMAX> mockLayerSDLsDNcm;
      for (auto& m : mockLayerSDLsDNcm) m.reserve(1000000);

      
      if (useSeeds == 1){
	std::cout<<"Convert seeds to SuperDoublets"<<std::endl;
	auto const nSeeds = see_px().size();
	for (auto iSeed = 0U; iSeed < nSeeds; ++iSeed){
	  TVector3 p3LH(see_stateTrajGlbPx()[iSeed], see_stateTrajGlbPy()[iSeed], see_stateTrajGlbPz()[iSeed]);
	  float ptLH = p3LH.Pt();
	  float etaLH = p3LH.Eta();
	  if (ptLH < ptCutAll) continue;
	  if (std::abs(etaLH) > (addEndcaps ? 3.0f : 1.6f)) continue;
	  TVector3 r3LH(see_stateTrajGlbX()[iSeed], see_stateTrajGlbY()[iSeed], see_stateTrajGlbZ()[iSeed]);

	  TVector3 p3PCA(see_px()[iSeed], see_py()[iSeed], see_pz()[iSeed]);
	  TVector3 r3PCA(r3FromPCA(p3PCA, see_dxy()[iSeed], see_dz()[iSeed]));

	  auto const& seedHitsV = see_hitIdx()[iSeed];
	  auto const& seedHitTypesV = see_hitType()[iSeed];
	  auto seedAlgo = TrackAlgorithm(see_algo()[iSeed]);
	  if (seedAlgo != TrackAlgorithm::initialStep) continue;//FIXME: make configurable
	  int nHits = seedHitsV.size();
	  if (debugReco){
	    std::cout<<"Seed with nHits "<<nHits<<" algo "<<int(seedAlgo)<<std::endl;
	    for (int iH = 0U; iH < nHits; ++iH){
	      int simIdx = -1;
	      if (HitType(seedHitTypesV[iH]) == HitType::Pixel){
		auto ipix = seedHitsV[iH];
		simIdx = simHitsPerPixHitAll[ipix];
		std::cout<<" isPixel "<<"\t ";
	      }
	      if (simIdx >= 0 ) SimHit(simIdx).print("");
	      else std::cout<<std::endl;
	    }
	  }
	  
	  assert(nHits == 4);
	  for (int iH = 0; iH < nHits; ++iH){
	    //FIXME: make this configurable
	    assert(HitType(seedHitTypesV[iH]) == HitType::Pixel);
	  }
	  SuperDoublet seedSD;
	  seedSD.mdRef.pixL = HitIndexWithType(see_hitIdx()[iSeed][0], HitType(see_hitType()[iSeed][0])).indexWithType;
	  seedSD.mdRef.pixU = HitIndexWithType(see_hitIdx()[iSeed][1], HitType(see_hitType()[iSeed][1])).indexWithType;
	  seedSD.mdRef.r3 = r3PCA;
	  seedSD.mdRef.rt = r3PCA.Pt();	  
	  seedSD.mdRef.z = r3PCA.Z();	  
	  seedSD.mdRef.r = r3PCA.Mag();	  
	  seedSD.mdRef.phi = r3PCA.Phi();	  
	  seedSD.mdRef.alpha = r3PCA.DeltaPhi(p3PCA);
	  const int itpRL = simsPerHitAll(seedSD.mdRef.pixL);
	  const int itpRU = simsPerHitAll(seedSD.mdRef.pixU);
	  seedSD.mdRef.itp = itpRL;
	  seedSD.mdRef.ntp = 1;
	  seedSD.mdRef.itpLL = itpRL;
	  if (itpRL >= 0 && itpRU >= 0){
	    if (itpRL == itpRU ){
	      seedSD.mdRef.ntp = 2;
	    }
	  } else if (itpRL == -1 && itpRU == -1){
	    seedSD.mdRef.ntp = 0;
	  } else if (itpRU >= 0){
	    seedSD.mdRef.itp = itpRU;
	  }
	  seedSD.mdOut.pixL = HitIndexWithType(see_hitIdx()[iSeed][2], HitType(see_hitType()[iSeed][2])).indexWithType;
	  if (nPix >= 4) seedSD.mdOut.pixU = HitIndexWithType(see_hitIdx()[iSeed][3], HitType(see_hitType()[iSeed][3])).indexWithType;
	  seedSD.mdOut.r3 = r3LH;
	  seedSD.mdOut.rt = r3LH.Pt();
	  seedSD.mdOut.z = r3LH.Z();
	  seedSD.mdOut.r = r3LH.Mag();
	  seedSD.mdOut.phi = r3LH.Phi();
	  seedSD.mdOut.alpha = r3LH.DeltaPhi(p3LH);	  
	  const int itpOL = simsPerHitAll(seedSD.mdOut.pixL);
	  const int itpOU = simsPerHitAll(seedSD.mdOut.pixU);
	  seedSD.mdOut.itp = itpOL;
	  seedSD.mdOut.ntp = 1;
	  seedSD.mdOut.itpLL = itpOL;
	  if (itpOL >= 0 && itpOU >= 0){
	    if (itpOL == itpOU ){
	      seedSD.mdOut.ntp = 2;
	    }
	  } else if (itpOL == -1 && itpOU == -1){
	    seedSD.mdOut.ntp = 0;
	  } else if (itpOU >= 0){
	    seedSD.mdOut.itp = itpOU;
	  }
	  seedSD.iRef = iSeed;
	  seedSD.iOut = iSeed;
	  seedSD.r3 = r3LH;
	  seedSD.rt = r3LH.Pt();
	  seedSD.rtInv = 1.f/seedSD.rt;
	  seedSD.z = seedSD.r3.Z();
	  seedSD.p3 = p3LH;
	  seedSD.alpha = r3LH.DeltaPhi(p3LH);
	  seedSD.dr = (r3LH - r3PCA).Pt();
	  seedSD.d = seedSD.rt - r3PCA.Pt();
	  seedSD.zeta = seedSD.p3.Pt()/seedSD.p3.Z();
	  //
	  std::map<int, int> tps;
	  seedSD.itp = -1;
	  seedSD.ntp = 0;
	  tps[itpRL]++;  tps[itpRU]++;  tps[itpOL]++;  tps[itpOU]++;
	  for ( auto m : tps){
	    if (m.first >= 0 && m.second > seedSD.ntp){
	      seedSD.itp = m.first;
	      seedSD.ntp = m.second;
	    }
	  }
	  //LL indexing is not very useful for seeds; to it anyways for "full" coverage
	  seedSD.itpLL = -1;
	  seedSD.ntpLL = 0;
	  tps.clear();
	  tps[itpRL]++;  tps[itpOL]++;
	  for ( auto m : tps){
	    if (m.first >= 0 && m.second > seedSD.ntpLL){
	      seedSD.itpLL = m.first;
	      seedSD.ntpLL = m.second;
	    }
	  }
	  mockLayerSDfwDNcm[0].emplace_back(seedSD);
	}
	std::cout<<"Loaded nSeeds (pt>1 and |eta|<1.6) "<<mockLayerSDfwDNcm[0].size()<<std::endl;
      }      
      
      std::cout<<"Combine MDs "<<std::endl;
      int nHitsTried = 0;
      for (int iL = minLayer; iL <= nLayersA; ++iL){
	if (!addEndcaps && iL > 10) continue;

	auto const&  hitsRefLower = mockLayerMDfwRefLower[iL];
	auto const&  hitsRefUpper = mockLayerMDfwRefUpper[iL];
	auto const&  hitsOutLower = mockLayerMDfwDNcmLower[iL];
	auto const&  hitsOutUpper = mockLayerMDfwDNcmUpper[iL];
	
	auto& mockMDfwRef = mockLayerMDfwRef[iL];

	int n_all = 0;
	int n_dz = 0;
	int n_dr = 0;
	int n_dPhiPos = 0;
	int n_dPhi = 0;
	
	auto mdCombine = [&] (decltype(hitsRefLower) hitsL, decltype(hitsRefUpper) hitsU, decltype(mockMDfwRef) mDs){
	  n_all = 0;
	  n_dz = 0;
	  n_dr = 0;
	  n_dPhiPos = 0;
	  n_dPhi = 0;
	  for (auto const& hL : hitsL) {
	    nHitsTried++; //if (nHitsTried>1) exit(0);
	    float rt = hL.second.rt;
	    const float ptCut = ptCutAll;
	    const float miniSlope = std::asin(std::min(rt*k2Rinv1GeVf/ptCut, sinAlphaMax));
	    
	    const float miniMuls = miniMulsPtScale[iL]*3.f/ptCut;
	    const float rLayNominal = iL >= minLayer ? (iL < 11 ? miniRminMean[iL] : rt ) : 1e12f;
	    const float miniPVoff = 0.1f/rLayNominal;
	    const float miniCut = miniSlope + sqrt(miniMuls*miniMuls + miniPVoff*miniPVoff);
	    
	    const float dzCut = 10.f;//may want to adjust this: PS modules are shorter
	    
	    for (auto const& hU : hitsU) {
	      float dPhi = 0;
	      n_all++;
	      if (iL< 11){ //barrel
		auto const dz = hL.second.r3.z() - hU.second.r3.z();
		if (std::abs(dz) > dzCut) continue;
		n_dz++;
		n_dr++;
		
		const float dPhiPos = std::abs(deltaPhi(hU.second.phi, hL.second.phi));
		//FIXME: can be tighter
		if (dPhiPos > miniCut) continue;
		n_dPhiPos++;
		
		dPhi = hL.second.r3.DeltaPhi(hU.second.r3-hL.second.r3);
		if (std::abs(dPhi) > miniCut) continue;
		n_dPhi++;
	      } else { //endcap
		auto const dz = hU.second.r3.z() - hL.second.r3.z();//could enforce dz from geometry
		if (std::abs(dz) > 1.0f) continue; //max mini-layer separation is 4 mm
		n_dz++;
		
		auto const dr = hL.second.rt - hU.second.rt;
		if (std::abs(dr) > dzCut) continue;
		n_dr++;

		const float miniLum = useFullR3Endcap ? 0.f : deltaZLum/std::abs(hL.second.r3.z());
		const float miniCutE = miniSlope + sqrt(miniMuls*miniMuls + miniPVoff*miniPVoff + miniLum*miniLum);
		
		const float dPhiPos = deltaPhi(hL.second.phi, hU.second.phi);
		//FIXME: can be tighter
		if (std::abs(dPhiPos) > miniCutE) continue;
		n_dPhiPos++;

		const float dzFrac = dz/hL.second.r3.z();
		dPhi = dPhiPos/dzFrac*(1.f + dzFrac);		
		if (useFullR3Endcap){
		  dPhi = hL.second.r3.DeltaPhi(hU.second.r3-hL.second.r3);//NOTE: this changes combinatorial component. Use only for efficiency studies
		}
		if (std::abs(dPhi) > miniCutE) continue;
		n_dPhi++;
	      }
	      
	      MiniDoublet md;
	      md.pixL = hL.first;
	      md.pixU = hU.first;
	      md.r3 = hL.second.r3;
	      md.rt = hL.second.rt;
	      md.z = hL.second.r3.Z();
	      md.r = hL.second.r;
	      md.phi = hL.second.phi;
	      if (hL.second.rtRHin > 0){
		md.rtRHin = hL.second.rtRHin;
		md.rtRHout = hL.second.rtRHout;
		md.phiRHin = hL.second.phiRHin;
		md.phiRHout = hL.second.phiRHout;
	      } else {
		md.rtRHin = hL.second.rt;
		md.rtRHout = hL.second.rt;
		md.phiRHin = hL.second.phi;
		md.phiRHout = hL.second.phi;
	      }
	      md.alpha = dPhi;
	      md.alphaRHmin = md.alpha;//FIXME: make this better (when it matters)
	      md.alphaRHmax = md.alpha;
	      const int itpRL = simsPerHit(md.pixL);
	      const int itpRU = simsPerHit(md.pixU);
	      md.itp = itpRL;
	      md.ntp = 1;
	      md.itpLL = itpRL;
	      if (itpRL >= 0 && itpRU >= 0){
		if (itpRL == itpRU ){
		  md.ntp = 2;
		}
	      } else if (itpRL == -1 && itpRU == -1){
		md.ntp = 0;
	      } else if (itpRU >= 0){
		md.itp = itpRU;
	      }

	      
	      if (debugReco){
		const float miniLum = useFullR3Endcap ? 0.f : deltaZLum/std::abs(hL.second.r3.z());
		const float miniCutE = miniSlope + sqrt(miniMuls*miniMuls + miniPVoff*miniPVoff + miniLum*miniLum);
		if (iL==11&& n_dPhi<10){
		  float simPtL = itpRL >= 0 ? sqrt(sim_px()[itpRL]*sim_px()[itpRL]+sim_py()[itpRL]*sim_py()[itpRL]) : 0;
		  float simPtU = itpRU >= 0 ? sqrt(sim_px()[itpRU]*sim_px()[itpRU]+sim_py()[itpRU]*sim_py()[itpRU]) : 0;
		  float simDxyL = itpRL >= 0 ? sim_pca_dxy()[itpRL] : 99;
		  float simDxyU = itpRU >= 0 ? sim_pca_dxy()[itpRU] : 99;


		  std::cout<<"MD on "<<iL<<" i "<<mDs.size()
			   <<" "<<md.pixL<<" : "<<md.pixU
			   <<" "<<itpRL<<" : "<<itpRU
			   <<" pt "<<simPtL<<" "<<simPtU
			   <<" dxy "<<simDxyL<<" "<<simDxyU
			   <<" r "<<hL.second.rt<<" "<<hU.second.rt
			   <<" z "<<hL.second.r3.z()<<" "<<hU.second.r3.z()
			   <<" phi "<<hL.second.phi<<" "<<hU.second.phi
			   <<" phi in out "<<md.phiRHin<<" "<<md.phiRHout
			   <<" alpha "<<md.alpha<<" cutM "<<miniCutE
			   <<" alphaR3 "<<hL.second.r3.DeltaPhi(hU.second.r3-hL.second.r3)
			   <<std::endl;
		}
	      }
	      mDs.emplace_back(md);
	    }
	  }
	};
	mdCombine(hitsRefLower, hitsRefUpper, mockMDfwRef);
	std::cout<<"MD ref stat "<<iL<<" "<<mockMDfwRef.size()
		 <<" all "<<n_all
		 <<" dz "<<n_dz
		 <<" dr "<<n_dr
		 <<" dPhiPos "<<n_dPhiPos
		 <<" dPhi "<<n_dPhi
		 <<std::endl;

	auto& mockMDfwDNcm = mockLayerMDfwDNcm[iL];
	mdCombine(hitsOutLower, hitsOutUpper, mockMDfwDNcm);
	std::cout<<"MD out stat "<<iL<<" "<<mockMDfwDNcm.size()
		 <<" all "<<n_all
		 <<" dz "<<n_dz
		 <<" dr "<<n_dr
		 <<" dPhiPos "<<n_dPhiPos
 		 <<" dPhi "<<n_dPhi
		 <<std::endl;
	
	//now make super-doublets
	auto& mockSDfwDNcm = mockLayerSDfwDNcm[iL];
	
	enum SDSelectFlags {deltaZ = 0, deltaPhiPos, slope, alphaRef, alphaOut, alphaRefOut, max};
	std::array<string, SDSelectFlags::max> sdFlagString {"deltaZ", "deltaPhiPos",
	    "slope", "alphaRef", "alphaOut", "alphaRefOut"};
	std::array<int, SDSelectFlags::max> sdMasksCumulative {};
	for (int i = 0; i < SDSelectFlags::max; ++i){
	  if (i > 0 ) sdMasksCumulative[i] = sdMasksCumulative[i-1];
	  sdMasksCumulative[i] |= 1 << i;
	}
	
	std::array<int, SDSelectFlags::max> nPass {};
	int nAll = 0;
	
	int iRef = -1;
	const float zGeom = mockMode == 3 ?
	  (iL >= 5 && iL <= 6 ? 2.f*pixelPSZpitch : 2.f*strip2SZpitch)//twice the macro-pixel or strip size
	  : 0.f;//precise hits are used otherwise
	                                                //assume that the mock layer is the n+1 layer

	for (auto const& mdRef : mockMDfwRef){
	  iRef++;
	  float rtRef = mdRef.r3.Pt();
	  float zRef = mdRef.r3.z();
	  //
	  int iOut = -1;
	  auto dr3 = mdRef.r3;
	  for (auto const& mdOut : mockMDfwDNcm){
	    nAll++;
	    int sdFlag = 0;
	    iOut++;
	    float rtOut = mdOut.rt;
	    float zOut = mdOut.z;

	    bool debug_sdBuild = false;
	    if (mdRef.itp == mdOut.itp && mdRef.itp >= 0 && mdRef.ntp == 2 && mdOut.ntp == 2){
	      TVector3 p3TP(sim_px()[mdRef.itp], sim_py()[mdRef.itp], sim_pz()[mdRef.itp]);
	      if (std::abs(p3TP.Pt()-1.65171) < 1e-5 && (iL == 5 || iL == 11)){
		debug_sdBuild = true;
		std::cout<<"Debug for  "<<p3TP.Pt()<<" "<<p3TP.Eta()<<" "<<p3TP.Phi()<<" "<<iL<<std::endl;
	      }
	    }
	    
	    float dPhi;
	    unsigned int iFlag;
	    const float ptCut = ptCutAll;
	    float dPhiOut;

	    float dPhiRHmin;
	    float dPhiRHmax;
	    float dPhiOutRHmin;
	    float dPhiOutRHmax;
	    
	    if (iL < 11){//barrel
	      float rt = rtOut; //FIXME: should be different for mockMode = 0
	      const float sdSlope = std::asin(std::min(rt*k2Rinv1GeVf/ptCut, sinAlphaMax));
	      const float dzDrtScale = tan(sdSlope)/sdSlope;//FIXME: need approximate value

	      //apply some loose Z compatibility
	      float zLo = zRef + (zRef - deltaZLum)*(rtOut/rtRef - 1.f)*(zRef > 0.f ? 1.f : dzDrtScale) - zGeom; //slope-correction only on outer end
	      float zHi = zRef + (zRef + deltaZLum)*(rtOut/rtRef - 1.f)*(zRef < 0.f ? 1.f : dzDrtScale) + zGeom;
	      iFlag = SDSelectFlags::deltaZ;
	      if (!(zOut < zLo || zOut > zHi)) sdFlag |= 1 << iFlag;
	      else if (cumulativeCuts ){
		if (debug_sdBuild) std::cout<<"Failed SDSelectFlags::deltaZ "<<zOut<<" "<<zLo<<" "<<zHi<<std::endl;
		continue;
	      }
	      
	      if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;
	      
	      const float sdMuls = miniMulsPtScale[iL]*3.f/ptCut*2.f;//will need a better guess than x2?
	      const float sdPVoff = 0.1f/rt;
	      const float sdCut = sdSlope + sqrt(sdMuls*sdMuls + sdPVoff*sdPVoff);
	      
	      iFlag = SDSelectFlags::deltaPhiPos;
	      //FIXME: should be tighter than the local sdCut
	      if (!(std::abs(deltaPhi(mdRef.phi, mdOut.phi)) > sdCut )) sdFlag |= 1 << iFlag;
	      else if (cumulativeCuts ){
		if (debug_sdBuild) std::cout<<"Failed SelectFlags::deltaPhiPos "<<mdRef.phi<<" "<<mdOut.phi<<" "<<sdCut<<std::endl;
		continue;
	      }
	      if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;
	      
	      dr3 = mdOut.r3; dr3 -= mdRef.r3;
	      //plain SD bend cut
	      dPhi = mdRef.r3.DeltaPhi(dr3);
	      dPhiRHmin = dPhi;
	      dPhiRHmax = dPhi;
	      
	      iFlag = SDSelectFlags::slope;
	      if (!(std::abs(dPhi) > sdCut)) sdFlag |= 1 << iFlag;
	      else if (cumulativeCuts ){
		if (debug_sdBuild) std::cout<<"Failed SelectFlags::slope "<<dPhi<<" "<<sdCut<<std::endl;
		continue;
	      }
	      if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;

	      dPhiOut = mdOut.r3.DeltaPhi(dr3);
	      dPhiOutRHmin = dPhiOut;
	      dPhiOutRHmax = dPhiOut;
	    } else {//endcap
	      const float rtGeom = mockMode == 3 ?
		(rtRef < disks2SMinRadius && rtOut < disks2SMinRadius ? 2.f*pixelPSZpitch
		 : (rtRef < disks2SMinRadius || rtOut < disks2SMinRadius ) ? (pixelPSZpitch+strip2SZpitch)
		 : 2.f*strip2SZpitch) //FIXME: make this chosen by configuration for lay11,12 full PS
		: 0.f;//precise hits are used
	      
	      if (zRef*zOut < 0) continue; //do not even accumulate stats for wrong side combinations
	      
	      const float dz = zOut - zRef;
	      const float rt = rtOut; //FIXME: should be different for mockMode = 0
	      const float sdSlope = std::asin(std::min(rt*k2Rinv1GeVf/ptCut, sinAlphaMax));

	      //apply some loose Z compatibility
	      const float dLum = std::copysign(deltaZLum, zRef);
	      const float drtDzScale = sdSlope/tan(sdSlope);//FIXME: need approximate value
	      float rtLo = std::max(rtRef*(1.f + dz/(zRef + dLum)*drtDzScale) - rtGeom, rtRef - 0.5f*rtGeom); //rt should increase
	      float rtHi = rtRef*(zOut - dLum)/(zRef - dLum) + rtGeom; //dLum for luminous; rGeom for measurement size; no tanTheta_loc(pt) correction
	      iFlag = SDSelectFlags::deltaZ; //some unfortunate naming
	      if (!(rtOut < rtLo || rtOut > rtHi)) sdFlag |= 1 << iFlag;
	      else if (cumulativeCuts ){
		if (debug_sdBuild) std::cout<<"Failed SelectFlags::deltaZ "<<rtOut<<" "<<rtLo<<" "<<rtHi<<std::endl;
		continue;
	      }
	      
	      if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;

	      const float sdMuls = miniMulsPtScale[iL]*3.f/ptCut*2.f;//will need a better guess than x2?
	      const float sdPVoff = 0.1f/rt;
	      const float sdLum = useFullR3Endcap ? 0.f : deltaZLum/std::abs(zRef);
	      const float sdCut = sdSlope + sqrt(sdMuls*sdMuls + sdPVoff*sdPVoff + sdLum*sdLum);

	      const float dPhiPos = deltaPhi(mdRef.phi, mdOut.phi);
	      const float dPhiPosRHin = deltaPhi(mdRef.phiRHin, mdOut.phi);
	      const float dPhiPosRHout = deltaPhi(mdRef.phiRHout, mdOut.phi);
	      iFlag = SDSelectFlags::deltaPhiPos;
	      //FIXME: should be tighter than the local sdCut
	      if (!(std::abs(dPhiPos) > sdCut )) sdFlag |= 1 << iFlag;
	      else if (cumulativeCuts ){
		if (debug_sdBuild) std::cout<<"Failed SelectFlags::deltaPhiPos "<<mdRef.phi<<" "<<mdOut.phi<<" "<<sdCut<<std::endl;
		continue;
	      }
	      if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;
	      
	      //equivalent SD bend cut
	      const float dzFrac = dz/mdRef.r3.z();
	      dPhi = dPhiPos/dzFrac*(1.f + dzFrac);
	      const float dPhiRHin = dPhiPosRHin/dzFrac*(1.f + dzFrac);
	      const float dPhiRHout = dPhiPosRHout/dzFrac*(1.f + dzFrac);
	      dPhiRHmin = dPhi;
	      dPhiRHmax = dPhi;

	      if (mdRef.rtRHin > 0){assert(mockMode == 3);//complete range check is needed for full RH use
		dPhiRHmin = dPhiRHin;
		dPhiRHmax = dPhiRHout;
		if (std::abs(dPhiRHmin) > std::abs(dPhiRHmax)) std::swap(dPhiRHmax, dPhiRHmin);
	      }
	      
	      dr3 = mdOut.r3; dr3 -= mdRef.r3;//not needed for cuts but needed below
	      if (useFullR3Endcap) dPhi = mdRef.r3.DeltaPhi(dr3);

	      iFlag = SDSelectFlags::slope;
	      if (!(std::abs(dPhi) > sdCut )) sdFlag |= 1 << iFlag;
	      else if (cumulativeCuts ){
		if (debug_sdBuild) std::cout<<"Failed SelectFlags::slope "<<dPhi<<" "<<dPhiPos<<" "<<dzFrac<<" "<<sdCut<<std::endl;
		continue;
	      }
	      if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;

	      //	      const float dPhiAsBarrel = mdRef.r3.DeltaPhi(dr3);
	      //	      dPhi = dPhiAsBarrel;

	      dPhiOut = deltaPhi(dPhi, dPhiPos);
	      if (useFullR3Endcap) dPhiOut = mdOut.r3.DeltaPhi(dr3);

	      dPhiOutRHmin = dPhiOut;
	      dPhiOutRHmax = dPhiOut;

	      if (mdRef.rtRHin > 0){assert(mockMode == 3);//complete range check is needed for full RH use
		dPhiOutRHmin = deltaPhi(dPhiRHin, dPhiPosRHin);
		dPhiOutRHmax = deltaPhi(dPhiRHout, dPhiPosRHout);
		if (std::abs(dPhiOutRHmin) > std::abs(dPhiOutRHmax)) std::swap(dPhiOutRHmax, dPhiOutRHmin);
	      }
	      if (debug_sdBuild){
		std::cout<<"sdBuild mdRef "<<mdRef.phi<<" "<<mdRef.phiRHin<<" "<<mdRef.phiRHout<<" r "<<mdRef.rt<<" "<<mdRef.rtRHin<<" "<<mdRef.rtRHout
			 <<std::endl;
		std::cout<<"sdBuild range: "<<dPhiPos<<" "<<dPhi<<" "<<deltaPhi(dPhi, dPhiPos)
			 <<" in: "<<dPhiPosRHin<<" "<<dPhiRHin<<" "<<deltaPhi(dPhiRHin, dPhiPosRHin)
			 <<" out: "<<dPhiPosRHout<<" "<<dPhiRHout<<" "<<deltaPhi(dPhiRHout, dPhiPosRHout)
			 <<std::endl;
	      }
	    }
	    if (debug_sdBuild){
	      std::cout<<" sdBuild dPhis: "<<dPhi<<" "<<dPhiRHmin<<" "<<dPhiRHmax
		       <<" dPhiOut "<<dPhiOut<<" "<<dPhiOutRHmin<<" "<<dPhiOutRHmax
		       <<std::endl;
	    }

	    
	    SuperDoublet sd;
	    sd.iRef = iRef;
	    sd.iOut = iOut;

	    sd.r3 = mdRef.r3;
	    sd.rt = mdRef.rt;
	    sd.rtInv = 1.f/mdRef.rt;
	    sd.z = sd.r3.Z();
	    sd.alpha = dPhi;
	    sd.alphaRHmin = dPhiRHmin;
	    sd.alphaRHmax = dPhiRHmax;
	    sd.dr = dr3.Pt();//FIXME: get a better estimator for endcap
	    sd.d = mdOut.rt - sd.rt; //FIXME: get a better estimator for endcap
	    sd.zeta = sd.d/(mdOut.r3.Z() - sd.z);

	    //loose angle compatibility
	    float dAlpha_Bfield = std::asin(std::min(sd.dr*k2Rinv1GeVf/ptCut, sinAlphaMax));
	    float dAlpha_res = 0.04f/miniDelta[iL] * (iL < 11 ? 1.0f : std::abs(sd.z/sd.rt) )*(mockMode == 3 ? 1.0f : 0.0f);//4-strip difference
	    float dAlpha_compat = dAlpha_Bfield + dAlpha_res;
	    if (debug_sdBuild) dAlpha_compat *= 2;
	    iFlag = SDSelectFlags::alphaRef;
	    if (!((std::abs(mdRef.alpha- sd.alphaRHmax) > dAlpha_compat)
		  && (std::abs(mdRef.alpha- sd.alphaRHmin) > dAlpha_compat)) ) sdFlag |= 1 << iFlag;
	    else if (cumulativeCuts ){
	      if (debug_sdBuild) std::cout<<"Failed SelectFlags::alphaRef "<<mdRef.alpha<<" "<<sd.alphaRHmax<<" "<<sd.alphaRHmin<<" "<<dAlpha_compat<<std::endl;
	      continue;
	    }

	    if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;
	      
	    iFlag = SDSelectFlags::alphaOut;
	    if (!((std::abs(mdOut.alpha- sd.alphaRHmax) > dAlpha_compat)
		  && (std::abs(mdOut.alpha- sd.alphaRHmin) > dAlpha_compat)) ) sdFlag |= 1 << iFlag;//FIXME: this could be more restrictive: dBfiled cancels out
	    else if (cumulativeCuts ){
	      if (debug_sdBuild) std::cout<<"Failed SelectFlags::alphaOut "<<mdOut.alpha<<" "<<sd.alphaRHmax<<" "<<sd.alphaRHmin<<" "<<dAlpha_compat<<std::endl;
	      continue;
	    }

	    if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;

	    iFlag = SDSelectFlags::alphaRefOut;
	    if (!(std::abs(mdOut.alpha- mdRef.alpha) > dAlpha_compat)) sdFlag |= 1 << iFlag;
	    else if (cumulativeCuts ){
	      if (debug_sdBuild) std::cout<<"Failed SelectFlags::alphaRefOut "<<mdRef.alpha<<" "<<mdOut.alpha<<" "<<dAlpha_compat<<std::endl;
	      continue;
	    }

	    if (sdFlag == sdMasksCumulative[iFlag]) nPass[iFlag]++;

	    if (sdFlag != sdMasksCumulative[SDSelectFlags::max-1]) continue; //apply all cuts
	    if (debug_sdBuild) std::cout<<"Passed all selections! "<<std::endl;

	    sd.mdRef = mdRef;
	    sd.mdOut = mdOut;
	    sd.alphaOut = dPhiOut;
	    sd.alphaOutRHmin = dPhiOutRHmin;
	    sd.alphaOutRHmax = dPhiOutRHmax;

	    const int itpRL = simsPerHit(sd.mdRef.pixL);
	    const int itpRU = simsPerHit(sd.mdRef.pixU);
	    const int itpOL = simsPerHit(sd.mdOut.pixL);
	    const int itpOU = simsPerHit(sd.mdOut.pixU);
	    //
	    std::map<int, int> tps;
	    sd.itp = -1;
	    sd.ntp = 0;
	    tps[itpRL]++;  tps[itpRU]++;  tps[itpOL]++;  tps[itpOU]++;
	    for ( auto m : tps){
	      if (m.first >= 0 && m.second > sd.ntp){
		sd.itp = m.first;
		sd.ntp = m.second;
	      }
	    }
	    sd.itpLL = -1;
	    sd.ntpLL = 0;
	    tps.clear();
	    tps[itpRL]++;  tps[itpOL]++;
	    for ( auto m : tps){
	      if (m.first >= 0 && m.second > sd.ntpLL){
		sd.itpLL = m.first;
		sd.ntpLL = m.second;
	      }
	    }
	    
	    if (debugReco && iL == 11 && sd.ntp == 4){
	      TVector3 p3TP(sim_px()[sd.itp], sim_py()[sd.itp], sim_pz()[sd.itp]);
	      if (p3TP.Pt() > 2){
		std::cout<<" TP "<<p3TP.Pt()<<" "<<p3TP.Eta()<<" "<<p3TP.Phi()
			 <<" a "<<sd.alpha<<" aMin "<<sd.alphaRHmin<<" aMax "<<sd.alphaRHmax
			 <<" mrA "<<mdRef.alpha<<" moA "<<mdOut.alpha<<" aCompat "<<dAlpha_compat<<std::endl;
	      }
	    }
	    
	    if ( sd.dr > 4.0f*(rtOut - rtRef) && iL < 11){
	      //problem in matching
	      std::cout<<__LINE__
		       <<" "<<sd.mdRef.r3.Pt()<<" "<<sd.mdRef.r3.Phi()
		       <<" "<<sd.mdOut.r3.Pt()<<" "<<sd.mdOut.r3.Phi()
		       <<std::endl;
	    }
	    mockSDfwDNcm.emplace_back(sd);
	  }//mdOut : mockMDfwDNcm
	}//mdRef : mockMDfwRef

	int nSDs = mockSDfwDNcm.size();
	mockLayerSDfwDNcm_isSecondaryGhost[iL].resize(nSDs);
	for (int i = 0; i< nSDs; ++i) mockLayerSDfwDNcm_isSecondaryGhost[iL][i] = false;

	//fill isSecondaryGhost flags such that the first encounter is "false"
	for (int i = 0; i< nSDs; ++i){
	  bool isSecondaryGhost = mockLayerSDfwDNcm_isSecondaryGhost[iL][i];
	  if (! isSecondaryGhost){
	    auto const& iSD = mockSDfwDNcm[i];
	    for (int j = i+1; j < nSDs; ++j){
	      auto const& jSD = mockSDfwDNcm[j];
	      if (sameLowerLayerID(iSD, jSD)){
		mockLayerSDfwDNcm_isSecondaryGhost[iL][j] = true;
		break;
	      }
	    }
	  }
	}
	std::cout<<"SD stat "<<iL
		 <<" nAll "<<nAll;
	for (unsigned int i = 0; i< SDSelectFlags::max; ++i){
	  std::cout<<" "<<sdFlagString[i]<<" "<<nPass[i];
	}
	std::cout<<std::endl;
      }//iL

      std::cout<<"Link SuperDoublet pairs to SDLinks"<<std::endl;

      enum SDLSelectFlags { deltaZ = 0, deltaZPointed=1, deltaPhiPos=2,
			    slope=3, dAlphaIn=4, dAlphaOut=5, dBeta=6, max=7};
      std::array<int, SDLSelectFlags::max> sdlMasksCumulative {};
      for (int i = 0; i < SDLSelectFlags::max; ++i){
	if (i > 0 ) sdlMasksCumulative[i] = sdlMasksCumulative[i-1];
	sdlMasksCumulative[i] |= 1 << i;
      }

      
      //try to link segments: (special cases for L=0 where needed)
      auto sdLink = [&] (int lIn, int lOut, typename decltype(mockLayerSDLsDNcm)::value_type& sdlV){
	auto const& sdInV  = mockLayerSDfwDNcm[lIn];
	auto const& sdOutV = mockLayerSDfwDNcm[lOut];	
	
	int iSDL = sdlFromLayers[lIn][lOut];
	assert(iSDL != SDL_LMAX);// "Layer pair is not mapped");

	int nAll = 0;
	int nDeltaZ = 0;
	int nDeltaZPointed = 0;
	int nDeltaPhiPos = 0;
	int nSlope = 0;
	int nInAlphaCompat = 0;
	int nOutAlphaCompat = 0;
	int ndBeta = 0;
	
	const float zGeom = mockMode == 3 ?
	(lIn == 0 ? 0.05f : ( (lIn >= 5 && lIn <= 7) ? pixelPSZpitch : ( (lIn >= 8 && lIn <= 10) ? strip2SZpitch : 0.f)))
	+//add the macro-pixel or strip size
	(lOut == 0 ? 0.05f : ( (lOut >= 5 && lOut <= 7) ? pixelPSZpitch : ( (lOut >= 8 && lOut <= 10) ? strip2SZpitch : 0.f)))	      
	: 0.f;//precise hits otherwise
	
	int iIn = -1;
	for ( auto const& sdIn : sdInV ) {
	  iIn++;
	  //
	  const float rtIn = sdIn.rt;
	  const float rtInvIn = sdIn.rtInv;
	  const float ptIn = sdIn.p3.Pt();
	  const float zIn = sdIn.r3.z();
	  const float rIn = sqrt(zIn*zIn + rtIn*rtIn);
	  const float drtSDIn = sdIn.d;
	  const float dzSDIn = sdIn.mdOut.z - sdIn.mdRef.z;
	  const float dr3SDIn = sdIn.mdOut.r - sdIn.mdRef.r;
	  
	  float ptSLo = ptCutAll;
	  if (lIn == 0){
	    //try to use seed pt: the lower bound is good
	    ptSLo = ptIn;
	    float ptErr = see_ptErr()[sdIn.iRef];
	    ptSLo = std::max(ptCutAll, ptSLo - 10.0f*std::max(ptErr, 0.005f*ptSLo));//FIXME: check high-pt behavior
	    ptSLo = std::min(10.0f, ptSLo); //don't let this run away either
	  }
	  const float ptCut = ptSLo;
	  
	  int iOut = -1;
	  for ( auto const& sdOut : sdOutV ) {
	    int sdlFlag = 0;
	    iOut++;
	    //
	    
	    const float rtOut = sdOut.rt;
	    const float zOut = sdOut.z;

	    //don't even track stats for the opposite ends in z
	    if (lOut >= 11
		&& ( (lIn > 0 && zIn*zOut < 0) || (lIn == 0 && sdIn.p3.Z()*zOut < 0))) continue;
	    nAll++;

	    const float rt = rtOut;
	    const float sdlSlope = std::asin(std::min(rt*k2Rinv1GeVf/ptCut, sinAlphaMax));
	    const float dzDrtScale = tan(sdlSlope)/sdlSlope;//FIXME: need approximate value

	    //FIXME (later) more realistic accounting of material effects is needed
	    const float sdlThetaMulsF = 0.015f*sqrt(0.1f + 0.2*(rtOut-rtIn)/50.f) * (lIn < 11 ? sqrt(rIn/rtIn) : 1.f);
	    const float sdlMuls = sdlThetaMulsF*3.f/ptCut*4.f;//will need a better guess than x4?

	    if (lOut < 11){//barrel: match to Z proper
	      //apply some loose Z compatibility
	      //FIXME: refine using inner layer directions (can prune later)
	      const float rtOut_o_rtIn = rtOut*rtInvIn;
	      const float zLo = zIn + (zIn - deltaZLum)*(rtOut_o_rtIn - 1.f)*(zIn > 0.f ? 1.f : dzDrtScale) - zGeom; //slope-correction only on outer end
	      if (zOut < zLo && cumulativeCuts) continue;

	      const float zHi = zIn + (zIn + deltaZLum)*(rtOut_o_rtIn - 1.f)*(zIn < 0.f ? 1.f : dzDrtScale) + zGeom;
	      if (!(zOut < zLo || zOut > zHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZ;
	      else if (cumulativeCuts ) continue;
	    } else {//endcap uses matching to r
	      const float dLum = std::copysign(deltaZLum, zIn);
	      if (lIn < 11){//B-E
		const float rtGeom1 = mockMode == 3 ?
		  (rtOut < disks2SMinRadius ? pixelPSZpitch : strip2SZpitch)//FIXME: make this chosen by configuration for lay11,12 full PS
		  : 0.f;//precise hits otherwise
		const float zGeom1 = std::copysign(zGeom,zIn); //used in B-E region

		const float rtLo = rtIn*(1.f + (zOut - zIn - zGeom1)/(zIn + zGeom1 + dLum)/dzDrtScale) - rtGeom1;//slope correction only on the lower end
		if (rtOut < rtLo && cumulativeCuts) continue;
		float zInForHi = zIn - zGeom1 - dLum;
		if (zInForHi*zIn < 0 ) zInForHi = std::copysign(0.1f, zIn);
		const float rtHi = rtIn*(1.f + (zOut - zIn + zGeom1)/zInForHi) + rtGeom1;
		if (!(rtOut < rtLo || rtOut > rtHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZ;
		else if (cumulativeCuts ) continue;
	      } else {//E-E
		const float rtGeom = mockMode == 3 ?
		  (rtIn < disks2SMinRadius && rtOut < disks2SMinRadius ? 2.f*pixelPSZpitch
		   : (rtIn < disks2SMinRadius || rtOut < disks2SMinRadius ) ? (pixelPSZpitch + strip2SZpitch)
		   : 2.f*strip2SZpitch) //FIXME: make this chosen by configuration for lay11,12 full PS
		  : 0.f;//precise hits otherwise
		
		const float dz = zOut - zIn;

		const float rtLo = rtIn*(1.f + dz/(zIn + dLum)/dzDrtScale) - rtGeom;//slope correction only on the lower end
		if (rtOut < rtLo && cumulativeCuts) continue;
		const float rtHi = rtIn*(1.f + dz/(zIn - dLum)) + rtGeom;
		if (!(rtOut < rtLo || rtOut > rtHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZ;
		else if (cumulativeCuts ) continue;		
	      }//if (lIn < 11){
	    }//if (lOut < 11){//barrel: match to Z proper
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::deltaZ]) nDeltaZ++;

	    const float drOutIn = (rtOut - rtIn);
	    
	    if (lIn == 0){
	      const float etaErr = see_etaErr()[sdIn.iRef];
	      const float seedPtOut = std::hypot(see_stateTrajGlbPx()[sdIn.iRef], see_stateTrajGlbPy()[sdIn.iRef]);
	      const float coshEta = std::hypot(seedPtOut, see_stateTrajGlbPz()[sdIn.iRef])/seedPtOut;

	      if (lOut < 11){//barrel
		float dzErr = drOutIn*etaErr*coshEta; //FIXME: check with the calc in the endcap
		dzErr *= dzErr;
		dzErr += 0.03f*0.03f; // pixel size x2. ... random for now
		dzErr *= 9.f; //3 sigma
		dzErr += sdlMuls*sdlMuls*drOutIn*drOutIn/3.f*coshEta*coshEta;//sloppy
		dzErr += zGeom*zGeom;
		dzErr = sqrt(dzErr);
		const float dzDrIn = sdIn.p3.Z()/ptIn;
		const float zWindow = dzErr/drtSDIn*drOutIn + zGeom;
		const float dzMean = dzDrIn*drOutIn*(1.f + drOutIn*drOutIn*kRinv1GeVf*kRinv1GeVf/ptIn/ptIn/24.f);//with curved path correction
		const float zLo = zIn + dzMean - zWindow;
		const float zHi = zIn + dzMean + zWindow;
		if (!(zOut < zLo || zOut > zHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		else if (cumulativeCuts ) continue;
	      } else {//endcap; !!! THIS IS BLIND, WAS NOT TESTED !!!
		const float dzOutInAbs = std::abs(zOut - zIn);
		const float multDzDr = dzOutInAbs*coshEta/(coshEta*coshEta - 1.f);
		const float rtGeom1 = mockMode == 3 ?
		  (rtOut < disks2SMinRadius ? pixelPSZpitch : strip2SZpitch)//FIXME: make this chosen by configuration for lay11,12 full PS
		  : 0.f;//precise hits otherwise
		
		float drtErr = etaErr*multDzDr;
		drtErr *= drtErr;
		drtErr += 0.03f*0.03f; // pixel size x2. ... random for now
		drtErr *= 9.f; //3 sigma
		drtErr += sdlMuls*sdlMuls*multDzDr*multDzDr/3.f*coshEta*coshEta;//sloppy: relative muls is 1/3 of total muls
		drtErr = sqrt(drtErr);
		const float drtDzIn = std::abs(ptIn/sdIn.p3.Z());//all tracks are out-going in endcaps?

		const float rtWindow = drtErr + rtGeom1;
		const float drtMean = drtDzIn*dzOutInAbs*(1.f - drOutIn*drOutIn*kRinv1GeVf*kRinv1GeVf/ptIn/ptIn/24.f);//with curved path correction
		const float rtLo = rtIn + drtMean - rtWindow;
		const float rtHi = rtIn + drtMean + rtWindow;
		if (!(rtOut < rtLo || rtOut > rtHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		else if (cumulativeCuts ) continue;
	      }//if (lOut < 11){//barrel
	    }
	    else if (lIn>=5 && lIn <=6){//can point to the z pos in lOut
	      if (lOut<11){//barrel
		const float coshEta = dr3SDIn/drtSDIn;//direction estimate
		float dzErr = zGeom*zGeom*2.f;//both sides contribute to direction uncertainty
		dzErr += sdlMuls*sdlMuls*drOutIn*drOutIn/3.f*coshEta*coshEta;//sloppy
		dzErr = sqrt(dzErr);
		const float dzMean = dzSDIn/drtSDIn*drOutIn;
		const float zWindow = dzErr/drtSDIn*drOutIn + zGeom; //FIXME for ptCut lower than ~0.8 need to add curv path correction
		const float zLo = zIn + dzMean*(zIn > 0.f ? 1.f : dzDrtScale) - zWindow;
		const float zHi = zIn + dzMean*(zIn < 0.f ? 1.f : dzDrtScale) + zWindow;
		if (!(zOut < zLo || zOut > zHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		else if (cumulativeCuts ) continue;
	      } else {//endcap
		const float coshEta = dr3SDIn/drtSDIn;//direction estimate
		const float dzOutInAbs = std::abs(zOut - zIn);
		const float multDzDr = dzOutInAbs*coshEta/(coshEta*coshEta - 1.f);
		const float rtGeom1 = mockMode == 3 ?
		  (rtOut < disks2SMinRadius ? pixelPSZpitch : strip2SZpitch)//FIXME: make this chosen by configuration for lay11,12 full PS
		  : 0.f;//precise hits otherwise
		const float zGeom1 = mockMode == 3 ? pixelPSZpitch : 0.f;

		const float kZ = (zOut - zIn)/dzSDIn;
		float drtErr = zGeom1*zGeom1*drtSDIn*drtSDIn/dzSDIn/dzSDIn * (1.f - 2.f*kZ + 2.f*kZ*kZ);//Notes:122316
		drtErr += sdlMuls*sdlMuls*multDzDr*multDzDr/3.f*coshEta*coshEta;//sloppy: relative muls is 1/3 of total muls
		drtErr = sqrt(drtErr);
		const float drtMean = drtSDIn*dzOutInAbs/std::abs(dzSDIn); //
		const float rtWindow = drtErr + rtGeom1;
		const float rtLo = rtIn + drtMean/dzDrtScale - rtWindow;
		const float rtHi = rtIn + drtMean + rtWindow;
		if (!(kZ < 0 || rtOut < rtLo || rtOut > rtHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		else if (cumulativeCuts ) continue;
	      }
	    }
	    else if (lIn>=11 &&
		     (sdIn.mdRef.rt < disks2SMinRadius && sdIn.mdOut.rt < disks2SMinRadius)//FIXME: make configurable for full inner disks PS
		     ){//can point to the r pos in lOut;
	      assert(lOut >= 13);//only endcaps to match to
	      const float coshEta = dr3SDIn/drtSDIn;//direction estimate
	      const float dzOutInAbs = std::abs(zOut - zIn);
	      const float multDzDr = dzOutInAbs*coshEta/(coshEta*coshEta - 1.f);
	      const float rtGeom = mockMode == 3 ?
		(rtIn < disks2SMinRadius && rtOut < disks2SMinRadius ? 2.f*pixelPSZpitch
		 : (rtIn < disks2SMinRadius || rtOut < disks2SMinRadius ) ? (pixelPSZpitch+strip2SZpitch)
		 : 2.f*strip2SZpitch) //FIXME: make this chosen by configuration for lay11,12 full PS
		: 0.f;//precise hits otherwise

	      float drtErr = pixelPSZpitch*pixelPSZpitch*2.f/dzSDIn/dzSDIn*dzOutInAbs*dzOutInAbs;//both sides contribute to direction uncertainty
	      drtErr += sdlMuls*sdlMuls*multDzDr*multDzDr/3.f*coshEta*coshEta;//sloppy: relative muls is 1/3 of total muls
	      drtErr = sqrt(drtErr);
	      const float drtMean = drtSDIn*dzOutInAbs/std::abs(dzSDIn);
	      const float rtWindow = drtErr + rtGeom; //
	      const float rtLo = rtIn + drtMean/dzDrtScale - rtWindow;
	      const float rtHi = rtIn + drtMean + rtWindow;
	      if (!(rtOut < rtLo || rtOut > rtHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
	      else if (cumulativeCuts ) continue;
	    }
	    else {
	      sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
	    }
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::deltaZPointed]) nDeltaZPointed++;
	    
	    const float sdlPVoff = 0.1f/rt;
	    const float sdlCut = sdlSlope + sqrt(sdlMuls*sdlMuls + sdlPVoff*sdlPVoff);

	    if (lIn == 0){
	      sdlFlag |= 1 << SDLSelectFlags::deltaPhiPos;
	    } else {
	      //FIXME: can be tighter
	      if (! (std::abs(deltaPhi(sdIn.mdOut.phi, sdOut.mdOut.phi)) > sdlCut) ) sdlFlag |= 1 << SDLSelectFlags::deltaPhiPos;
	      else if (cumulativeCuts ) continue;
	    }
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::deltaPhiPos]) nDeltaPhiPos++;

	    auto const midR3 = 0.5f*(sdIn.r3 + sdOut.r3);
	    const float dPhi = midR3.DeltaPhi(sdOut.r3 - sdIn.r3);

	    if (! (std::abs(dPhi) > sdlCut) ) sdlFlag |= 1 << SDLSelectFlags::slope;
	    else if (cumulativeCuts ) continue;
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::slope]) nSlope++;
	    
	    float betaIn;
	    float betaOut;
	    float betaInRHmin;
	    float betaInRHmax;
	    float betaOutRHmin;
	    float betaOutRHmax;

	    auto const dr3 = mockMode == 0 ? sdOut.r3 - sdIn.r3 : sdOut.mdOut.r3 - sdIn.r3;
	    if (mockMode == 0){
	      betaIn = sdIn.alpha - sdIn.r3.DeltaPhi(dr3);
	      betaOut = - sdOut.alpha + sdOut.r3.DeltaPhi(dr3); //to match sign for correct match

	      betaInRHmin = betaIn;
	      betaInRHmax = betaIn;
	      betaOutRHmin = betaOut;
	      betaOutRHmax = betaOut;
	    }
	    else if (mockMode == 1 || mockMode == 3){
	      //plain segment-level definitions; uneven rotation corrections are applied later using a better estimate of pt
	      if (lIn == 0){
		betaIn  = -sdIn.p3.DeltaPhi(dr3);
		betaOut = -sdOut.alphaOut + sdOut.mdOut.r3.DeltaPhi(dr3);
	      } else {
		//need a symmetric choice of end-points to achieve partial cancelation
		betaIn  = sdIn.alpha - sdIn.r3.DeltaPhi(dr3);
		betaOut = -sdOut.alphaOut + sdOut.mdOut.r3.DeltaPhi(dr3);
	      }
	      if (mockMode == 1 || lOut < 11){
		betaInRHmin = betaIn;
		betaInRHmax = betaIn;
		betaOutRHmin = betaOut;
		betaOutRHmax = betaOut;
	      } else {
		//FIXME: dr3 part should be varied as well
		if (lIn < 11){
		  betaInRHmin = betaIn;
		  betaInRHmax = betaIn;		  
		} else {
		  betaInRHmin = betaIn + sdIn.alphaRHmin - sdIn.alpha;
		  betaInRHmax = betaIn + sdIn.alphaRHmax - sdIn.alpha;
		  if (std::abs(betaInRHmin) > std::abs(betaInRHmax)) std::swap(betaInRHmax, betaInRHmin);
		}
		if (lOut < 11){
		  betaOutRHmin = betaOut;
		  betaOutRHmax = betaOut;		  
		} else {
		  betaOutRHmin = betaOut - sdOut.alphaOutRHmin + sdOut.alphaOut;
		  betaOutRHmax = betaOut - sdOut.alphaOutRHmax + sdOut.alphaOut;
		  if (std::abs(betaOutRHmin) > std::abs(betaOutRHmax)) std::swap(betaOutRHmax, betaOutRHmin);
		}
	      }
	    }
	    else{
	      betaIn = -99;
	      betaOut = 999;
	      betaInRHmin = betaIn;
	      betaInRHmax = betaIn;
	      betaOutRHmin = betaOut;
	      betaOutRHmax = betaOut;
	    }
	    
	    const float dr = dr3.Perp();
	    //beta upper cuts: 2-strip difference for direction resolution
	    const float corrF = mockMode == 0 ? 0.f : 1.f;
	    bool pass_betaIn_cut = lIn == 0;//pixel seeds were already selected
	    if (lIn != 0){
	      const float betaIn_cut = (-sdIn.dr*corrF + dr)*k2Rinv1GeVf/ptCut + (mockMode == 3 ? 0.02f/sdIn.d : 0.f);
	      pass_betaIn_cut = std::abs(betaInRHmin) < betaIn_cut;
	    }
	    if (pass_betaIn_cut) sdlFlag |=  1 << SDLSelectFlags::dAlphaIn;
	    else if (cumulativeCuts ) continue;
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::dAlphaIn]) nInAlphaCompat++;

	    
	    
	    //now the actual segment linking magic
	    float betaAv = 0.5f*(betaIn + betaOut);
	    //pt/k2Rinv1GeVf/2. = R
	    //R*sin(betaAv) = pt/k2Rinv1GeVf/2*sin(betaAv) = dr/2 => pt = dr*k2Rinv1GeVf/sin(betaAv);
	    float pt_beta = dr*k2Rinv1GeVf/sin(betaAv);
	    if (lIn == 0) pt_beta = ptIn;

	    const float pt_betaMax = 7.0f;
	    
	    //apply segment (SD) bend correction
	    if (mockMode == 1 || mockMode == 3){
	      if (lIn == 0){
		betaOut += copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaOut);//FIXME: need a faster version
	      } else {
		const float diffDr = std::abs(sdIn.dr - sdOut.dr)/std::abs(sdIn.dr + sdOut.dr);
		if (true //do it for all//diffDr > 0.05 //only if segment length is different significantly
		    && betaIn*betaOut > 0.f
		    && (std::abs(pt_beta) < 4.f*pt_betaMax
			|| (lIn >= 11 && std::abs(pt_beta) < 8.f*pt_betaMax) )){ //and the pt_beta is well-defined; less strict for endcap-endcap
		  const float betaInUpd  = betaIn + copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		  const float betaOutUpd = betaOut + copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaOut);//FIXME: need a faster version
		  betaAv = 0.5f*(betaInUpd + betaOutUpd);
		  pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate
		  betaIn  += copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		  betaOut += copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaOut);//FIXME: need a faster version
		  //update the av and pt
		  betaAv = 0.5f*(betaIn + betaOut);
		  pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate		  
		}  else if (lIn < 11 && std::abs(betaOut) < 0.2* std::abs(betaIn) && std::abs(pt_beta) < 12.f*pt_betaMax){//use betaIn sign as ref
		  const float pt_betaIn = dr*k2Rinv1GeVf/sin(betaIn);
		  const float betaInUpd  = betaIn + copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_betaIn), sinAlphaMax)), betaIn);//FIXME: need a faster version
		  const float betaOutUpd = betaOut + copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_betaIn), sinAlphaMax)), betaIn);//FIXME: need a faster version
		  betaAv = std::abs(betaOut) > 0.2f*std::abs(betaIn) ? 0.5f*(betaInUpd + betaOutUpd) : betaInUpd;
		  pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate
		  betaIn  += copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		  betaOut += copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		  //update the av and pt
		  betaAv = 0.5f*(betaIn + betaOut);
		  pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate		  		      
		}
	      }
	    }
	    //rescale the ranges proportionally
	    const float betaInMMSF = 2.f*betaIn/std::abs(betaInRHmin+betaInRHmax);//mean value of min,max is the old betaIn
	    const float betaOutMMSF = 2.f*betaOut/std::abs(betaOutRHmin+betaOutRHmax);
	    betaInRHmin*= betaInMMSF;
	    betaInRHmax*= betaInMMSF;
	    betaOutRHmin*= betaOutMMSF;
	    betaOutRHmax*= betaOutMMSF;

	    const float dBetaMuls = sdlThetaMulsF*4.f/std::min(std::abs(pt_beta), pt_betaMax);//need to confirm the range-out value of 7 GeV
	    
	    //regularize to alpha of pt_betaMax .. review may want to add resolution
	    const float alphaInAbsReg = useFullR3Endcap ? 0.f : std::max(std::abs(sdIn.alpha), std::asin(std::min(sdIn.rt*k2Rinv1GeVf/3.0f, sinAlphaMax)));
	    const float alphaOutAbsReg = useFullR3Endcap ? 0.f : std::max(std::abs(sdOut.alpha), std::asin(std::min(sdOut.rt*k2Rinv1GeVf/3.0f, sinAlphaMax)));
	    const float dBetaInLum = lIn < 11 ? 0.0f : std::abs(alphaInAbsReg*deltaZLum/sdIn.z);
	    const float dBetaOutLum = lOut < 11 ? 0.0f : std::abs(alphaOutAbsReg*deltaZLum/sdOut.z);
	    const float dBetaLum2 = (dBetaInLum + dBetaOutLum)*(dBetaInLum + dBetaOutLum);

	    const float sinDPhi = std::sin(dPhi);
	    const float dBetaRIn2 = std::pow((sdIn.mdRef.rtRHout - sdIn.mdRef.rtRHin)*sinDPhi/dr, 2);
	    const float dBetaROut2 = std::pow((sdOut.mdOut.rtRHout - sdOut.mdOut.rtRHin)*sinDPhi/dr, 2);
	    
	    const float betaOut_cut = std::asin(std::min(dr*k2Rinv1GeVf/ptCut, sinAlphaMax))//FIXME: need faster version
	      + (mockMode == 3 ? 0.02f/sdOut.d : 0.f) + sqrt(dBetaLum2 + dBetaMuls*dBetaMuls);	    
	    if (std::abs(betaOut) < betaOut_cut) sdlFlag |=  1 << SDLSelectFlags::dAlphaOut;
	    else if (cumulativeCuts ) continue;
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::dAlphaOut]) nOutAlphaCompat++;

	    float pt_betaIn = dr*k2Rinv1GeVf/sin(betaIn);
	    if (lIn == 0) pt_betaIn = pt_beta;
	    const float pt_betaOut = dr*k2Rinv1GeVf/sin(betaOut);

	    const float dBetaRes = mockMode == 3 ? 0.02f/std::min(sdOut.d,sdIn.d) : 0.f;

	    const float dBetaCut2 = (dBetaRes*dBetaRes*2.0f + dBetaMuls*dBetaMuls + dBetaLum2 + dBetaRIn2 + dBetaROut2
				     + 0.25*std::pow(std::abs(betaInRHmin - betaInRHmax) + std::abs(betaOutRHmin - betaOutRHmax),2));

	    const float dBeta = betaIn - betaOut;

	    const float dZeta = sdIn.zeta - sdOut.zeta;
	    
	    if (!(dBeta*dBeta > dBetaCut2)) sdlFlag |= 1 << SDLSelectFlags::dBeta;
	    //stay on for N-1 plots	    else if (cumulativeCuts ) continue;
	    if (sdlFlag == sdlMasksCumulative[SDLSelectFlags::dBeta]) ndBeta++;

	    auto const ptInEst = std::abs(pt_betaIn);


	    //	    std::cout<<"Fill histograms for sdlFlag "<<sdlFlag<<std::endl;
	    //special case no "-1"
	    ha_SDL_dBeta_0_all[iSDL]->Fill(dBeta);
	    ha_SDL_dBeta_zoom_0_all[iSDL]->Fill(dBeta);
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::deltaZ]) == sdlMasksCumulative[SDLSelectFlags::deltaZ]){
	      ha_SDL_dBeta_0_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_0_pass[iSDL]->Fill(dBeta);
	    }

	    
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::deltaZPointed-1]) == sdlMasksCumulative[SDLSelectFlags::deltaZPointed-1] ){
	      ha_SDL_dBeta_01_all[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_01_all[iSDL]->Fill(dBeta);
	    }
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::deltaZPointed]) == sdlMasksCumulative[SDLSelectFlags::deltaZPointed] ){
	      ha_SDL_dBeta_01_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_01_pass[iSDL]->Fill(dBeta);
	    }
	    
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::deltaPhiPos-1]) == sdlMasksCumulative[SDLSelectFlags::deltaPhiPos-1]){
	      ha_SDL_dBeta_012_all[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_012_all[iSDL]->Fill(dBeta);
	    }
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::deltaPhiPos]) == sdlMasksCumulative[SDLSelectFlags::deltaPhiPos]){
	      ha_SDL_dBeta_012_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_012_pass[iSDL]->Fill(dBeta);
	    }

	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::slope-1]) == sdlMasksCumulative[SDLSelectFlags::slope-1]){
	      ha_SDL_dBeta_0123_all[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_0123_all[iSDL]->Fill(dBeta);
	    }
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::slope]) == sdlMasksCumulative[SDLSelectFlags::slope]){
	      ha_SDL_dBeta_0123_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_0123_pass[iSDL]->Fill(dBeta);
	    }
	    
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::dAlphaIn-1]) == sdlMasksCumulative[SDLSelectFlags::dAlphaIn-1]){
	      ha_SDL_dBeta_01234_all[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_01234_all[iSDL]->Fill(dBeta);
	    }
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::dAlphaIn]) == sdlMasksCumulative[SDLSelectFlags::dAlphaIn]){
	      ha_SDL_dBeta_01234_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_01234_pass[iSDL]->Fill(dBeta);
	    }
	    
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::dAlphaOut-1]) == sdlMasksCumulative[SDLSelectFlags::dAlphaOut-1]){
	      ha_SDL_dBeta_012345_all[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_012345_all[iSDL]->Fill(dBeta);
	    }
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::dAlphaOut]) == sdlMasksCumulative[SDLSelectFlags::dAlphaOut]){
	      ha_SDL_dBeta_012345_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_012345_pass[iSDL]->Fill(dBeta);
	    }
	    
	    //all inputs to make an SDL are available earlier
	    //but for performance reason, this is kept just before the place it's needed
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

	    const int itpIRL = lIn == 0 ? simsPerHitAll(sdl.sdIn.mdRef.pixL) : simsPerHit(sdl.sdIn.mdRef.pixL);
	    const int itpIRU = lIn == 0 ? simsPerHitAll(sdl.sdIn.mdRef.pixU) : simsPerHit(sdl.sdIn.mdRef.pixU);
	    const int itpIOL = lIn == 0 ? simsPerHitAll(sdl.sdIn.mdOut.pixL) : simsPerHit(sdl.sdIn.mdOut.pixL);
	    const int itpIOU = lIn == 0 ? simsPerHitAll(sdl.sdIn.mdOut.pixU) : simsPerHit(sdl.sdIn.mdOut.pixU);
	    const int itpORL = simsPerHit(sdl.sdOut.mdRef.pixL);
	    const int itpORU = simsPerHit(sdl.sdOut.mdRef.pixU);
	    const int itpOOL = simsPerHit(sdl.sdOut.mdOut.pixL);
	    const int itpOOU = simsPerHit(sdl.sdOut.mdOut.pixU);
	    //
	    std::map<int, int> tps;
	    sdl.itp = -1;
	    sdl.ntp = 0;
	    tps[itpIRL]++;  tps[itpIRU]++;  tps[itpIOL]++;  tps[itpIOU]++;
	    tps[itpORL]++;  tps[itpORU]++;  tps[itpOOL]++;  tps[itpOOU]++;
	    for ( auto m : tps){
	      if (m.first >= 0 && m.second > sdl.ntp){
		sdl.itp = m.first;
		sdl.ntp = m.second;
	      }
	    }
	    sdl.sdInGhost  = lIn != 0 && mockLayerSDfwDNcm_isSecondaryGhost[lIn][iIn];
	    sdl.sdOutGhost = mockLayerSDfwDNcm_isSecondaryGhost[lOut][iOut];
	    sdl.itpLL = -1;
	    sdl.ntpLL = 0;
	    tps.clear();
	    tps[itpIRL]++;  tps[itpIOL]++;
	    tps[itpORL]++;  tps[itpOOL]++;
	    for ( auto m : tps){
	      if (m.first >= 0 && m.second > sdl.ntpLL){
		sdl.itpLL = m.first;
		sdl.ntpLL = m.second;
	      }
	    }

	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::dBeta-1]) == sdlMasksCumulative[SDLSelectFlags::dBeta-1]){
	      ha_SDL_dBeta_NM1dBeta_all[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_betaIn_NM1dBeta_all[iSDL]->Fill(betaIn, dBeta);
	      
	      ha_SDL_dBeta_zoom_NM1dBeta_all[iSDL]->Fill(dBeta);
	      if (! ( (lIn != 0 && mockLayerSDfwDNcm_isSecondaryGhost[lIn][iIn])
		      || mockLayerSDfwDNcm_isSecondaryGhost[lOut][iOut])
		  ){
		//		if (!(lOut >= 11 && sdl.sdOut.mdRef.rt < disks2SMinRadius))
		  {
		    ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[iSDL]->Fill(dBeta);
		    //NB: no pt cut applied here on the TP
		    if (sdl.ntpLL == 4 ) ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL]->Fill(dBeta);
		    else if (sdl.ntpLL == 3) ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL]->Fill(dBeta);
		    else ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL]->Fill(dBeta);
		    
		    ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL[iSDL]->Fill(dZeta);
		    //NB: no pt cut applied here on the TP
		    if (sdl.ntpLL == 4 ) ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL]->Fill(dZeta);
		    else if (sdl.ntpLL == 3) ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL]->Fill(dZeta);
		    else ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL]->Fill(dZeta);
		  }
		
	      }
	      
	      if (ptInEst < 2){
		ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all[iSDL]->Fill(dBeta);
	      }
	      if (ptInEst > 3 && ptInEst < 5 ){
		ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all[iSDL]->Fill(dBeta);
	      }
	      if (ptInEst > 7) {
		ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all[iSDL]->Fill(dBeta);
	      }
	    }
	    
	    if ((sdlFlag & sdlMasksCumulative[SDLSelectFlags::dBeta]) == sdlMasksCumulative[SDLSelectFlags::dBeta]){
	      ha_SDL_dBeta_NM1dBeta_pass[iSDL]->Fill(dBeta);
	      ha_SDL_dBeta_zoom_NM1dBeta_pass[iSDL]->Fill(dBeta);
	      if (ptInEst < 2){
		ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass[iSDL]->Fill(dBeta);
	      }
	      if (ptInEst > 3 && ptInEst < 5 ){
		ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass[iSDL]->Fill(dBeta);
	      }
	      if (ptInEst > 7) {
		ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass[iSDL]->Fill(dBeta);
	      }
	    }
	    
	    //	    std::cout<<"Done filling histograms"<<std::endl;
	    if (sdlFlag != sdlMasksCumulative[SDLSelectFlags::dBeta]) continue; //apply all cuts up to including dBeta
	    
	    sdlV.emplace_back(sdl);
	    //	    std::cout<<"Appended a new SDLink "<<std::endl;
	  }//sdOutV
	}//sdInV

	std::cout<<"SD links stat "<<lIn<<"-"<<lOut
	         <<" nAll "<<nAll<<" nDZ "<<nDeltaZ<<" nDZPnt "<<nDeltaZPointed
	         <<" nPPos "<<nDeltaPhiPos<<" nAlp "<<nSlope
         	 <<" nInAlpC "<< nInAlphaCompat<<" nOutAlpC "<<nOutAlphaCompat
	         <<" ndBeta "<<ndBeta
	         <<" final "<<sdlV.size()
	         <<" noGhost "<<std::count_if(sdlV.begin(), sdlV.end(), [](SDLink const& sdl )->bool{return ! (sdl.sdInGhost || sdl.sdOutGhost);})
	         <<std::endl;
      };//auto sdLink

      auto isdll = SDL_L0to5; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
      isdll = SDL_L0to7; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
      isdll = SDL_L0to11; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
      if (ptCutAll >= 1.0f || (sdOffsetB <= 6.0f && sdOffsetE <= 6.0f)){//FIXME: NEED A MORE OPTIMAL SOLUTION
	isdll = SDL_L5to7; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
	isdll = SDL_L7to9; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
	
	isdll = SDL_L5to11; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
	isdll = SDL_L5to13; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
	isdll = SDL_L7to11; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
	isdll = SDL_L7to13; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
	isdll = SDL_L11to13; sdLink(layersSDL[isdll][0], layersSDL[isdll][1], mockLayerSDLsDNcm[isdll]);
      }
      
      //link the links to TrackLinks
      std::vector<TrackLink> tracks;
      int countMatchMidPoint = 0;
      int nSDLL_all = 0;
      int nSDLL_dC = 0;
      int nSDLL_pass = 0;
      int isdlIn = -1;
      if (ptCutAll >= 1.0f){//FIXME: NEED A MORE OPTIMAL SOLUTION
	for (auto const& sdlIn : mockLayerSDLsDNcm[SDL_L5to7] ){
	  bool hasOuter = false;
	  isdlIn++;
	  
	  int isdlOut = -1;
	  for (auto const& sdlOut : mockLayerSDLsDNcm[SDL_L7to9]){
	    isdlOut++;
	    if (sdlIn.lOut == sdlOut.lIn && sdlIn.iOut == sdlOut.iIn){
	      //shared mid-point
	      countMatchMidPoint++;
	      nSDLL_all++;
	      
	      auto const dPt = sdlIn.pt - sdlOut.pt;
	      auto const dC = 1.f/sdlIn.pt - 1.f/sdlOut.pt;
	      
	      if (std::abs(dC) > 0.05 ) continue;
	      nSDLL_dC++;
	      if (sdlIn.sdInGhost || sdlIn.sdOutGhost || sdlOut.sdOutGhost) continue;
	      nSDLL_pass++;
	      
	      if (debugReco){
		if (sdlIn.itp >=0 && sdlOut.itp >= 0){
		  TVector3 p3In(sim_px()[sdlIn.itp], sim_py()[sdlIn.itp], sim_pz()[sdlIn.itp]);
		  TVector3 p3Out(sim_px()[sdlOut.itp], sim_py()[sdlOut.itp], sim_pz()[sdlOut.itp]);
		  std::cout<<"TT: i "<<isdlIn<<" - "<<isdlOut<<" c "<<countMatchMidPoint
			   <<" itp "<<sdlIn.itp<<" - "<<sdlOut.itp
			   <<" ntp "<<sdlIn.ntp<<" - "<<sdlOut.ntp
			   <<" tppt "<<p3In.Pt()<<" - "<<p3Out.Pt()
			   <<" sdlpt "<<sdlIn.pt<<" - "<<sdlOut.pt
			   <<" # "<<sdlIn.sdIn.mdRef.pixL
			   <<" "<<sdlIn.sdIn.mdRef.pixU
			   <<" "<<sdlIn.sdIn.mdOut.pixL
			   <<" "<<sdlIn.sdIn.mdOut.pixU
			   <<" # "<<sdlIn.sdOut.mdRef.pixL
			   <<" "<<sdlIn.sdOut.mdRef.pixU
			   <<" "<<sdlIn.sdOut.mdOut.pixL
			   <<" "<<sdlIn.sdOut.mdOut.pixU
			   <<" # "<<sdlOut.sdOut.mdRef.pixL
			   <<" "<<sdlOut.sdOut.mdRef.pixU
			   <<" "<<sdlOut.sdOut.mdOut.pixL
			   <<" "<<sdlOut.sdOut.mdOut.pixU
			   <<std::endl;
		} else if (sdlIn.itp >=0){
		  TVector3 p3In(sim_px()[sdlIn.itp], sim_py()[sdlIn.itp], sim_pz()[sdlIn.itp]);
		  std::cout<<"TF: i "<<isdlIn<<" - "<<isdlOut<<" c "<<countMatchMidPoint
			   <<" itp "<<sdlIn.itp<<" - "<<sdlOut.itp
			   <<" ntp "<<sdlIn.ntp<<" - "<<sdlOut.ntp
			   <<" tppt "<<p3In.Pt()<<" - "<<0.f
			   <<" sdlpt "<<sdlIn.pt<<" - "<<sdlOut.pt
			   <<" # "<<sdlIn.sdIn.mdRef.pixL
			   <<" "<<sdlIn.sdIn.mdRef.pixU
			   <<" "<<sdlIn.sdIn.mdOut.pixL
			   <<" "<<sdlIn.sdIn.mdOut.pixU
			   <<" # "<<sdlIn.sdOut.mdRef.pixL
			   <<" "<<sdlIn.sdOut.mdRef.pixU
			   <<" "<<sdlIn.sdOut.mdOut.pixL
			   <<" "<<sdlIn.sdOut.mdOut.pixU
			   <<" # "<<sdlOut.sdOut.mdRef.pixL
			   <<" "<<sdlOut.sdOut.mdRef.pixU
			   <<" "<<sdlOut.sdOut.mdOut.pixL
			   <<" "<<sdlOut.sdOut.mdOut.pixU
			   <<std::endl;
		} else if (sdlOut.itp >= 0){
		  TVector3 p3Out(sim_px()[sdlOut.itp], sim_py()[sdlOut.itp], sim_pz()[sdlOut.itp]);
		  std::cout<<"FT: i "<<isdlIn<<" - "<<isdlOut<<" c "<<countMatchMidPoint
			   <<" itp "<<sdlIn.itp<<" - "<<sdlOut.itp
			   <<" ntp "<<sdlIn.ntp<<" - "<<sdlOut.ntp
			   <<" tppt "<<0.f<<" - "<<p3Out.Pt()
			   <<" sdlpt "<<sdlIn.pt<<" - "<<sdlOut.pt
			   <<" # "<<sdlIn.sdIn.mdRef.pixL
			   <<" "<<sdlIn.sdIn.mdRef.pixU
			   <<" "<<sdlIn.sdIn.mdOut.pixL
			   <<" "<<sdlIn.sdIn.mdOut.pixU
			   <<" # "<<sdlIn.sdOut.mdRef.pixL
			   <<" "<<sdlIn.sdOut.mdRef.pixU
			   <<" "<<sdlIn.sdOut.mdOut.pixL
			   <<" "<<sdlIn.sdOut.mdOut.pixU
			   <<" # "<<sdlOut.sdOut.mdRef.pixL
			   <<" "<<sdlOut.sdOut.mdRef.pixU
			   <<" "<<sdlOut.sdOut.mdOut.pixL
			   <<" "<<sdlOut.sdOut.mdOut.pixU
			   <<std::endl;	      
		} else {
		  std::cout<<"FF: i "<<isdlIn<<" - "<<isdlOut<<" c "<<countMatchMidPoint
			   <<" sdlpt "<<sdlIn.pt<<" - "<<sdlOut.pt
			   <<" # "<<sdlIn.sdIn.mdRef.pixL
			   <<" "<<sdlIn.sdIn.mdRef.pixU
			   <<" "<<sdlIn.sdIn.mdOut.pixL
			   <<" "<<sdlIn.sdIn.mdOut.pixU
			   <<" # "<<sdlIn.sdOut.mdRef.pixL
			   <<" "<<sdlIn.sdOut.mdRef.pixU
			   <<" "<<sdlIn.sdOut.mdOut.pixL
			   <<" "<<sdlIn.sdOut.mdOut.pixU
			   <<" # "<<sdlOut.sdOut.mdRef.pixL
			   <<" "<<sdlOut.sdOut.mdRef.pixU
			   <<" "<<sdlOut.sdOut.mdOut.pixL
			   <<" "<<sdlOut.sdOut.mdOut.pixU
			   <<std::endl;	      
		}
	      }//if debugReco
	      
	    }//match in-out at mid-point
	    
	  }//sdlOut
	}//sdlIn

	std::cout<<"nSDLL 5-7-9 "
		 <<" all "<<nSDLL_all
		 <<" dC "<<nSDLL_dC
		 <<" final "<<nSDLL_pass<<std::endl;
      }//ptCutAll >= 1.0f


      timerA[T_timeReco].Stop();

      
      std::map<int, int> nHitsStatSumMap;
      std::map<int, int> nHitsStatCntMap;

      timerA[T_timeValidation].Start(kFALSE);
      for (int iSim = 0; iSim < int(nSim); ++iSim){
	TVector3 p3(sim_px()[iSim], sim_py()[iSim], sim_pz()[iSim]);
	bool debug = false;
	auto tpPt = p3.Pt();
	if (tpPt < 0.5*ptCutAll ) continue;
	
	auto tpEta = p3.Eta();
	auto tpPhi = p3.Phi();

	auto tpDxy = sim_pca_dxy()[iSim];
	auto prodX = simvtx_x()[sim_parentVtxIdx()[iSim]];
	auto prodY = simvtx_y()[sim_parentVtxIdx()[iSim]];
	auto prodR2 = prodX*prodX + prodY*prodY;
	if (effForPromptTracks){
	  if (std::abs(tpDxy) > 0.05) continue;
	  if (prodR2 > 0.05*0.05) continue;
	}
	if (maxTPdxy > 0){
	  if (std::abs(tpDxy) > maxTPdxy) continue;
	  if (prodR2 > maxTPdxy*maxTPdxy) continue;
	}
	if (minTPdxy != 0){
	  if (std::abs(tpDxy) < minTPdxy) continue;
	  if (prodR2 < minTPdxy*minTPdxy) continue;
	}
	
	auto tpDz = sim_pca_dz()[iSim];
	auto prodZ = simvtx_z()[sim_parentVtxIdx()[iSim]];
	if (maxTPdz > 0){
	  if (std::abs(tpDz) > maxTPdz) continue;
	  if (std::abs(prodZ) > maxTPdz) continue;
	}

	bool isLowProdXY = prodR2 < 4.0;
	
	if (tpPt > 25 && std::abs(p3.Eta())< 1 && isLowProdXY && sim_bunchCrossing()[iSim] == 0) debug = false;
	if (debug) std::cout<<"TP "<< iSim<<" : "<<p3.Pt()<<" "<<p3.Eta()<<" "<<p3.Phi();
	
	std::map<int, int> nHitsMap;
	std::array<std::vector<SDLink>, SDL_LMAX > matchingSDLs_byHit3of4 {};
	for (auto& s : matchingSDLs_byHit3of4) s.reserve(10);
	std::array<std::vector<SDLink>, SDL_LMAX > matchingSDLs_byHit4of4 {};
	for (auto& s : matchingSDLs_byHit4of4) s.reserve(10);
	std::array<std::set<int>, nLayersA+1> simHits {};
	
	if (runDetailedTimers) timerA[T_timeVal_SHLoad].Start(kFALSE);

	auto const& simhitIdxV = sim_simHitIdx()[iSim];
	int iPix = -1;
	for (auto ish : simhitIdxV){
	  iPix++;

	  SimHit simH(ish);
	  int lay = simH.lay;

	  // int iProcess = simH.iProcess;
	  // int ibx = simH.bx;	  
	  // bool isPrimaryAny = (iProcess == 2 && ibx == 0);
	  
	  if ((simH.isBarrel && lay >= minLayer) || lay < minLayer
	      || (addEndcaps && !simH.isBarrel && lay > 10)){
	    if (debug) std::cout<<" "<<lay<<" "<<ish;
	    //require pt consistency and a presence of an associated rechit
	    if (simH.p3s.Pt()>0.8*tpPt && !simhit_hitIdx()[ish].empty()){//FIXME: opt/configure for BX=0 requirement instead
	      nHitsMap[lay]++;
	      simHits[lay].emplace(ish);
	      if (lay < minLayer){//this goes to the seeds list
		nHitsMap[0]++;
		simHits[0].emplace(ish);		
	      }
	    }
	  }
	}//iPix : 0 to sim_nPixel
	if (runDetailedTimers) timerA[T_timeVal_SHLoad].Stop();
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
	  
	  if ( ( (lIn == 0 && nHitsMap[1] + nHitsMap[2] + nHitsMap[3] + nHitsMap[4] > 0) || ( lIn != 0 && nHitsMap[lIn] > 0))
	      && nHitsMap[lOut] > 0){
	    ha_denSDL_pt[iSDLL]->Fill(tpPt);
	    if (isLowProdXY) ha_denSDL_loProdXY_pt[iSDLL]->Fill(tpPt);
	    if (lIn == 0){
	      nHitsStatCntMap[lIn]++;
	      nHitsStatSumMap[lIn] += nHitsMap[1] + nHitsMap[2] + nHitsMap[3] + nHitsMap[4];
	      for (int l = 1; l<= 4; ++l){
		nHitsStatCntMap[l]++;
		nHitsStatSumMap[l] += nHitsMap[l];
	      }
	    } else {
	      nHitsStatCntMap[lIn]++;
	      nHitsStatSumMap[lIn] += nHitsMap[lIn];
	      nHitsStatCntMap[lOut]++;
	      nHitsStatSumMap[lOut] += nHitsMap[lOut];
	    }
	    
	    bool debugHitLevel = true;
	    if (debug){
	      std::cout<<"\tTP is good for denSDL in L"<<lIn<<"-L"<<lOut<<std::endl;
	      if (debugHitLevel){
		if (lIn == 0){
		  for (int iL = 1; iL < minLayer; ++iL){
		    for (auto i : simHits[iL]){ auto ph = SimHit(i); ph.print("\t");}
		  }
		} else {
		  for (auto i : simHits[lIn]){ auto ph = SimHit(i); ph.print("\t");}
		}
		for (auto i : simHits[lOut]){ auto ph = SimHit(i); ph.print("\t");}
	      }
	    }
	    //match the 8 layers of hits
	    auto matchMH = [&](decltype(mockLayerMDfwRefLower)::const_reference mhs){
	      for (auto const& mh : mhs) if (iSim == simsPerHit(mh.first)) return true;
	      return false;
	    };
	    auto matchIPix = [&](int ipix){
	      return (iSim == simsPerHit(ipix));
	    };
	    auto matchIPixAll = [&](int ipix){
	      return (iSim == simsPerHitAll(ipix));
	    };

	    const bool matchAllCombinations = false;

	    if (runDetailedTimers) timerA[T_timeVal_MHMDSDMatch].Start(kFALSE);	    
	    if (lIn == 0){
	      for (auto const& sd: mockLayerSDfwDNcm[lIn]){
		if (!hasMHRefInL) hasMHRefInL = iSim == simsPerHit(sd.mdRef.pixL);
	      }
	      if (hasMHRefInL || matchAllCombinations){
		for (auto const& sd: mockLayerSDfwDNcm[lIn]){
		  if (!hasMHRefInU) hasMHRefInU = iSim == simsPerHit(sd.mdRef.pixU);
		}
	      }
	      if (hasMHRefInU || matchAllCombinations){
		for (auto const& sd: mockLayerSDfwDNcm[lIn]){
		  if (!hasMHOutInL) hasMHOutInL = iSim == simsPerHit(sd.mdOut.pixL);
		}
	      }
	      if (hasMHOutInL || matchAllCombinations){
		for (auto const& sd: mockLayerSDfwDNcm[lIn]){
		  if (!hasMHOutInU) hasMHOutInU = iSim == simsPerHit(sd.mdOut.pixU);
		}
	      }
	    } else {
	      hasMHRefInL = matchMH(mockLayerMDfwRefLower[lIn]);
	      if (hasMHRefInL || matchAllCombinations) hasMHRefInU = matchMH(mockLayerMDfwRefUpper[lIn]);
	      if (hasMHRefInU || matchAllCombinations) hasMHOutInL = matchMH(mockLayerMDfwDNcmLower[lIn]);
	      if (hasMHOutInL || matchAllCombinations) hasMHOutInU = matchMH(mockLayerMDfwDNcmUpper[lIn]);
	    }

	    if (hasMHOutInU || matchAllCombinations) hasMHRefOutL = matchMH(mockLayerMDfwRefLower[lOut]);
	    if (hasMHRefOutL || matchAllCombinations) hasMHRefOutU = matchMH(mockLayerMDfwRefUpper[lOut]);
	    if (hasMHRefOutU || matchAllCombinations) hasMHOutOutL = matchMH(mockLayerMDfwDNcmLower[lOut]);
	    if (hasMHOutOutL || matchAllCombinations) hasMHOutOutU = matchMH(mockLayerMDfwDNcmUpper[lOut]);

	    has8MHs = hasMHRefInL & hasMHRefInU & hasMHOutInL & hasMHOutInU
	      & hasMHRefOutL & hasMHRefOutU & hasMHOutOutL & hasMHOutOutU;
	    if (debug){
	      std::cout<<"has8MHs  "<<has8MHs <<" = hasMHRefInL  "<<hasMHRefInL <<" & hasMHRefInU  "<<hasMHRefInU
		       <<" & hasMHOutInL  "<<hasMHOutInL <<" & hasMHOutInU  "<<hasMHOutInU
		       <<" & hasMHRefOutL  "<<hasMHRefOutL <<" & hasMHRefOutU  "<<hasMHRefOutU
		       <<" & hasMHOutOutL  "<<hasMHOutOutL <<" & hasMHOutOutU "<<hasMHOutOutU <<std::endl;
	    }
	    if (has8MHs){
	      ha_num8MH_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_num8MH_loProdXY_pt[iSDLL]->Fill(tpPt);
	    }
	    
	    //match the 4 layer mini-doublets
	    auto matchMD = [&](decltype(mockLayerMDfwRef)::const_reference mds){	      
	      for (auto const& md : mds) if (md.itp == iSim && md.ntp == 2 ) return true;
	      return false;
	    };
	    if (lIn == 0){
	      for (auto const& sd: mockLayerSDfwDNcm[lIn]){
		if (!hasMDRefIn) hasMDRefIn = matchIPixAll(sd.mdRef.pixL) & matchIPixAll(sd.mdRef.pixU);
		if (!hasMDOutIn) hasMDOutIn = matchIPixAll(sd.mdOut.pixL) & matchIPixAll(sd.mdOut.pixU);
	      }
	    } else {
	      hasMDRefIn = matchMD(mockLayerMDfwRef[lIn]);
	      hasMDOutIn = matchMD(mockLayerMDfwDNcm[lIn]);
	    }

	    hasMDRefOut = matchMD(mockLayerMDfwRef[lOut]);
	    hasMDOutOut = matchMD(mockLayerMDfwDNcm[lOut]);

	    has4MDs = hasMDRefIn & hasMDOutIn & hasMDRefOut & hasMDOutOut & has8MHs;
	    if (debug){
	      std::cout <<"has4MDs "<<has4MDs <<" = hasMDRefIn  "<<hasMDRefIn <<" & hasMDOutIn  "<<hasMDOutIn
			<<" & hasMDRefOut  "<<hasMDRefOut <<" & hasMDOutOut  "<<hasMDOutOut  << std::endl;
	    }
	    if (has4MDs ){
	      ha_num4MD_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_num4MD_loProdXY_pt[iSDLL]->Fill(tpPt);
	    }

	    std::vector<SuperDoublet> vSDIn_4of4; vSDIn_4of4.reserve(2);
	    std::vector<SuperDoublet> vSDOut_4of4; vSDOut_4of4.reserve(2);

	    bool debug2SD = false;
	    //	    if (tpPt>1.5 && tpPt< 2.0 && has4MDs && lIn == 0 && lOut == 11) debug2SD = true;
	    
	    //match inner and outer layer super-doublets
	    for (auto const& sd : mockLayerSDfwDNcm[lIn]){
	      int score = sd.itp == iSim ? sd.ntp : 0;
	      if (debug && score > 2){
		std::cout<<"SDI: match on L"<< lIn <<": Have "<<score<<" matches"<<std::endl;
	      }
	      if (score >= 3) hasSDIn_3of4 = true;
	      if (score >= 4){
		hasSDIn_4of4 = true;
		vSDIn_4of4.push_back(sd);
	      }
	    }//SD matching Inner
	    if (not hasSDIn_4of4 && debug2SD){
	      std::cout<<"missing/poor SDI match on L"<< lIn<<" for tpPt "<<tpPt<<" tpEta "<<tpEta<<std::endl;
	      const auto& simhitIdxV = sim_simHitIdx()[iSim];
	      for (auto ish : simhitIdxV){
		
		SimHit simH(ish);
		int lay = simH.lay;
		
		std::cout<<" "<<lay<<" "<<ish<<std::endl;		      
		if (simH.p3s.Pt()>-1.0f){
		  simH.print("\t");
		}
	      }//simhitIdxV
	      
	      for (auto const& sd : mockLayerSDfwDNcm[lIn]){
		int score = sd.itp == iSim ? sd.ntp : 0;
		if (score >= 2){
		  std::cout<<"            Have "<<score<<" matches"
			   <<" tpPt "<<tpPt<<" tpEta "<<tpEta
			   <<" ptIn "<<sd.p3.Pt()<<" EtaIn "<<sd.p3.Eta()
			   <<std::endl;
		  SimHit pRL(sd.mdRef.pixL); pRL.print("\tRL \t");
		  SimHit pRU(sd.mdRef.pixU); pRU.print("\tRU \t");
		  SimHit pOL(sd.mdOut.pixL); pOL.print("\tOL \t");
		  SimHit pOU(sd.mdOut.pixU); pOU.print("\tOU \t");
		}
	      }//SD matching Inner
	    }
	    
	    for (auto const& sd : mockLayerSDfwDNcm[lOut]){
	      int score = sd.itp == iSim ? sd.ntp : 0;
	      if (debug && score > 2){
		std::cout<<"SDO: match on L"<< lOut <<": Have "<<score<<" matches"<<std::endl;
	      }

	      if (score >= 3) hasSDOut_3of4 = true;
	      if (score >= 4){
		hasSDOut_4of4 = true;
		vSDOut_4of4.push_back(sd);
	      }
	    }//SD matching Outer
	    if (runDetailedTimers) timerA[T_timeVal_MHMDSDMatch].Stop();

	    //inner SD numerator
	    if (hasSDIn_3of4  && has4MDs){
	      ha_numInSD_3of4_any_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_numInSD_3of4_any_loProdXY_pt[iSDLL]->Fill(tpPt);
	    }
	    if (hasSDIn_4of4 && has4MDs){
	      ha_numInSD_4of4_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_numInSD_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
	    }

	    //both SD present
	    if (hasSDIn_3of4 && hasSDOut_3of4 && has4MDs){
	      ha_num2SD_3of4_any_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_num2SD_3of4_any_loProdXY_pt[iSDLL]->Fill(tpPt);
	    }
	    if (hasSDIn_4of4 && hasSDOut_4of4 && has4MDs){
	      ha_num2SD_4of4_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_num2SD_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
	    }

	    if (runDetailedTimers) timerA[T_timeVal_SDL].Start(kFALSE);
	    //enum SDLSelectFlags { deltaZ = 0, deltaZPointed, slope, dAlphaIn, dAlphaOut, dBeta};
	    std::vector<std::pair<SDLink, int> > vSDLwInfo_4of4;

	    const float zGeom = mockMode == 3 ?
	      (lIn == 0 ? 0.05f : ( (lIn >= 5 && lIn <= 7) ? pixelPSZpitch : ( (lIn >= 8 && lIn <= 10) ? strip2SZpitch : 0.f)))
	      +//add the macro-pixel or strip size
	      (lOut == 0 ? 0.05f : ( (lOut >= 5 && lOut <= 7) ? pixelPSZpitch : ( (lOut >= 8 && lOut <= 10) ? strip2SZpitch : 0.f)))	      
	      : 0.f;//precise hits otherwise

	    bool debugSimMatching = tpPt < 3 && tpPt > 2 && hasSDIn_4of4 && hasSDOut_4of4 && has8MHs && has4MDs && lIn == 5 && lOut == 11 && debug;
	    if (tpPt < 1.65172 && tpPt > 1.65170 && tpEta > 1.51427 && tpEta < 1.51429){
	      std::cout<<"debugSimMatching "<<debugSimMatching<<" "<<hasSDIn_4of4<<" "<<hasSDOut_4of4<<" "<<has8MHs<<" "<<has4MDs<<" "<<lIn<<" "<<lOut<<std::endl;
	    }
	    bool debugAllDbeta = false;
	    for (auto const& sdIn : vSDIn_4of4){
	      const float rtIn = sdIn.rt;
	      const float rtInvIn = sdIn.rtInv;
	      const float ptIn = sdIn.p3.Pt();
	      const float zIn = sdIn.r3.z();
	      const float rIn = sqrt(zIn*zIn + rtIn*rtIn);
	      const float drtSDIn = sdIn.d;
	      const float dzSDIn = sdIn.mdOut.z - sdIn.mdRef.z;
	      const float dr3SDIn = sdIn.mdOut.r - sdIn.mdRef.r;
	      
	      float ptSLo = ptCutAll;
	      if (lIn == 0){
		//try to use seed pt: the lower bound is good
		ptSLo = ptIn;
		float ptErr = see_ptErr()[sdIn.iRef];
		ptSLo = std::max(ptCutAll, ptSLo - 10.0f*std::max(ptErr, 0.005f*ptSLo));//FIXME: check high-pt behavior
		ptSLo = std::min(10.0f, ptSLo); //don't let this run away either
	      }
	      const float ptCut = ptSLo;
	      
	      for (auto const& sdOut : vSDOut_4of4){
		int sdlFlag = 0;
		
		const float rtOut = sdOut.rt;
		const float zOut = sdOut.r3.z();	    

		//don't even track stats for the opposite ends in z
		if (lOut >= 11
		    && ( (lIn > 0 && zIn*zOut < 0) || (lIn == 0 && sdIn.p3.Z()*zOut < 0))) continue;

		const float rt = rtOut;
		const float sdlSlope = std::asin(std::min(rt*k2Rinv1GeVf/ptCut, sinAlphaMax));
		const float dzDrtScale = tan(sdlSlope)/sdlSlope;//FIXME: need approximate value

		//FIXME (later) more realistic accounting of material effects is needed
		const float sdlThetaMulsF = 0.015f*sqrt(0.1f + 0.2*(rtOut-rtIn)/50.f) * (lIn < 11 ? sqrt(rIn/rtIn) : 1.f);
		const float sdlMuls = sdlThetaMulsF*3.f/ptCut*4.f;//will need a better guess than x4?

		if (lOut < 11){//barrel: match to Z proper
		  //apply some loose Z compatibility
		  //FIXME: refine using inner layer directions (can prune later)		
		  const float rtOut_o_rtIn = rtOut*rtInvIn;
		  const float zLo = zIn + (zIn - deltaZLum)*(rtOut_o_rtIn - 1.f)*(zIn > 0.f ? 1.f : dzDrtScale) - zGeom; //slope-correction only on outer end
		  const float zHi = zIn + (zIn + deltaZLum)*(rtOut_o_rtIn - 1.f)*(zIn < 0.f ? 1.f : dzDrtScale) + zGeom;
		
		  if (zOut > zLo && zOut < zHi) sdlFlag |= 1 << SDLSelectFlags::deltaZ;
		  if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZ) ){
		    std::cout<<"Lum region failed: tpPt "<<tpPt<<" lIn "<<lIn <<" zLo "<<zLo<<" zHi "<<zHi<<" vs zOut "<<zOut<<std::endl;
		  }
		} else {//endcap uses matching to r
		  const float dLum = std::copysign(deltaZLum, zIn);
		  if (lIn < 11){//B-E
		    const float rtGeom1 = mockMode == 3 ?
		      (rtOut < disks2SMinRadius ? pixelPSZpitch : strip2SZpitch)//FIXME: make this chosen by configuration for lay11,12 full PS
		      : 0.f;//precise hits otherwise
		    const float zGeom1 = std::copysign(zGeom,zIn);

		    const float rtLo = rtIn*(1.f + (zOut - zIn - zGeom1)/(zIn + zGeom1 + dLum)/dzDrtScale) - rtGeom1;//slope correction only on the lower end
		    float zInForHi = zIn - zGeom1 - dLum;
		    if (zInForHi*zIn < 0 ) zInForHi = std::copysign(0.1f, zIn);		    
		    const float rtHi = rtIn*(1.f + (zOut - zIn + zGeom1)/zInForHi) + rtGeom1;

		    if (rtOut > rtLo && rtOut < rtHi) sdlFlag |= 1 << SDLSelectFlags::deltaZ;
		    if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZ) ){
		      std::cout<<"Lum region failed: tpPt "<<tpPt<<" tpEta "<<tpEta<<" tpZ "<<prodZ
			       <<" lIn "<<lIn <<" rtLo "<<rtLo<<" rtHi "<<rtHi<<" vs rtOut "<<rtOut
			       <<" ptIn "<<ptIn<<" ptCut "<<ptCut
			       <<" rtIn "<<rtIn<<" zO "<<zOut<<" zI "<<zIn<<" zG1 "<<zGeom1<<" rtG1 "<<rtGeom1
			       <<std::endl;
		      auto const& simhitIdxV = sim_simHitIdx()[iSim];
		      for (auto ish : simhitIdxV){

			SimHit simH(ish);
			int lay = simH.lay;
			
			std::cout<<" "<<lay<<" "<<ish<<std::endl;		      
			if (simH.p3s.Pt()>0.8*tpPt){
			  simH.print("\t");
			}
		      }//iPix; pixelhits
		    }
		  } else {//E-E
		    const float rtGeom = mockMode == 3 ?
		      (rtIn < disks2SMinRadius && rtOut < disks2SMinRadius ? 2.f*pixelPSZpitch
		       : (rtIn < disks2SMinRadius || rtOut < disks2SMinRadius ) ? (pixelPSZpitch + strip2SZpitch)
		       : 2.f*strip2SZpitch) //FIXME: make this chosen by configuration for lay11,12 full PS
		      : 0.f;//precise hits otherwise
		    
		    const float dz = zOut - zIn;
		    
		    const float rtLo = rtIn*(1.f + dz/(zIn + dLum)/dzDrtScale) - rtGeom;//slope correction only on the lower end
		    const float rtHi = rtIn*(1.f + dz/(zIn - dLum)) + rtGeom;
		    if (rtOut > rtLo && rtOut < rtHi) sdlFlag |= 1 << SDLSelectFlags::deltaZ;
		    if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZ) ){
		      std::cout<<"Lum region failed: tpPt "<<tpPt<<" lIn "<<lIn <<" rtLo "<<rtLo<<" rtHi "<<rtHi<<" vs rtOut "<<rtOut<<std::endl;
		    }
		  }//if (lIn < 11){
		}//if (lOut < 11){//barrel: match to Z proper
		
		const float drOutIn = (rtOut - rtIn);

		if (lIn == 0){
		  const float etaErr = see_etaErr()[sdIn.iRef];
		  const float seedPtOut = std::hypot(see_stateTrajGlbPx()[sdIn.iRef], see_stateTrajGlbPy()[sdIn.iRef]);
		  const float coshEta = std::hypot(seedPtOut, see_stateTrajGlbPz()[sdIn.iRef])/seedPtOut;

		  if (lOut < 11){//barrel
		    float dzErr = drOutIn*etaErr*coshEta;
		    dzErr *= dzErr;
		    dzErr += 0.03f*0.03f; // pixel size x2. ... random for now
		    dzErr *= 9.f; //3 sigma
		    dzErr += sdlMuls*sdlMuls*drOutIn*drOutIn/3.f*coshEta*coshEta;//sloppy
		    dzErr += zGeom*zGeom;
		    dzErr = sqrt(dzErr);
		    const float dzDrIn = sdIn.p3.Z()/ptIn;
		    const float zWindow = dzErr/drtSDIn*drOutIn + zGeom;
		    const float dzMean = dzDrIn*drOutIn*(1.f + drOutIn*drOutIn*kRinv1GeVf*kRinv1GeVf/ptIn/ptIn/24.f);//with curved path correction
		    const float zLo = zIn + dzMean - zWindow;
		    const float zHi = zIn + dzMean + zWindow;
		    if (!(zOut < zLo || zOut > zHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed; //continue;
		    if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZPointed)){
		      std::cout<<"ZPointing failed: tpPt "<<tpPt<<" lIn "<<lIn <<" zLo "<<zLo<<" zHi "<<zHi<<" vs zOut "<<zOut
			       <<" : dzDrIn "<<dzDrIn<<" drtSDIn "<<drtSDIn<<" dzErr "<<dzErr<<" (rtOut - rtIn) "<<(rtOut - rtIn)<<std::endl;
		      std::cout<<"\t\t RefL "<<sdIn.mdRef.pixL<<" RefU "<<sdIn.mdRef.pixU
			       <<" OutL "<<sdIn.mdOut.pixL<<" OutU "<<sdIn.mdOut.pixU<<std::endl;
		    }
		  } else {//endcap; !!! THIS IS BLIND, WAS NOT TESTED !!!
		    const float dzOutInAbs = std::abs(zOut - zIn);
		    const float multDzDr = dzOutInAbs*coshEta/(coshEta*coshEta - 1.f);
		    const float rtGeom1 = mockMode == 3 ?
		      (rtOut < disks2SMinRadius ? pixelPSZpitch : strip2SZpitch)//FIXME: make this chosen by configuration for lay11,12 full PS
		      : 0.f;//precise hits otherwise
		    
		    float drtErr = etaErr*multDzDr;
		    drtErr *= drtErr;
		    drtErr += 0.03f*0.03f; // pixel size x2. ... random for now
		    drtErr *= 9.f; //3 sigma
		    drtErr += sdlMuls*sdlMuls*multDzDr*multDzDr/3.f*coshEta*coshEta;//sloppy; mulsThetaPos ~ 1/sqrt(3.)*mulsTheta
		    drtErr = sqrt(drtErr);
		    const float drtDzIn = std::abs(ptIn/sdIn.p3.Z());//all tracks are out-going in endcaps?
		    
		    const float rtWindow = drtErr + rtGeom1;
		    const float drtMean = drtDzIn*dzOutInAbs*(1.f - drOutIn*drOutIn*kRinv1GeVf*kRinv1GeVf/ptIn/ptIn/24.f);//with curved path correction
		    const float rtLo = rtIn + drtMean - rtWindow;
		    const float rtHi = rtIn + drtMean + rtWindow;
		    if (!(rtOut < rtLo || rtOut > rtHi)) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		    if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZPointed)){
		      std::cout<<"rtPointing failed: tpPt "<<tpPt<<" lIn "<<lIn <<" rtLo "<<rtLo<<" rtHi "<<rtHi<<" vs rtOut "<<rtOut
			       <<" : drtDzIn "<<drtDzIn<<" drtErr "<<drtErr<<" (rtOut - rtIn) "<<(rtOut - rtIn)<<std::endl;
		      std::cout<<"\t\t RefL "<<sdIn.mdRef.pixL<<" RefU "<<sdIn.mdRef.pixU
			       <<" OutL "<<sdIn.mdOut.pixL<<" OutU "<<sdIn.mdOut.pixU<<std::endl;
		    }
		  }//if (lOut < 11){//barrel
		}
		else if (lIn>=5 && lIn <=6){//can point to the z pos in lOut
		  if (lOut<11){//barrel
		    const float coshEta = dr3SDIn/drtSDIn;//direction estimate
		    float dzErr = zGeom*zGeom*2.f;//both sides contribute to direction uncertainty
		    dzErr += sdlMuls*sdlMuls*drOutIn*drOutIn/3.f*coshEta*coshEta;//sloppy
		    dzErr = sqrt(dzErr);
		    const float dzMean = dzSDIn/drtSDIn*drOutIn;
		    const float zWindow = dzErr/drtSDIn*drOutIn + zGeom; //FIXME for ptCut lower than ~0.8 need to add curv path correction
		    const float zLo = zIn + dzMean*(zIn > 0.f ? 1.f : dzDrtScale) - zWindow;
		    const float zHi = zIn + dzMean*(zIn < 0.f ? 1.f : dzDrtScale) + zWindow;
		    if (zOut > zLo && zOut < zHi) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		  } else {//endcap
		    const float coshEta = dr3SDIn/drtSDIn;//direction estimate
		    const float dzOutInAbs = std::abs(zOut - zIn);
		    const float multDzDr = dzOutInAbs*coshEta/(coshEta*coshEta - 1.f);
		    const float rtGeom1 = mockMode == 3 ?
		      (rtOut < disks2SMinRadius ? pixelPSZpitch : strip2SZpitch)//FIXME: make this chosen by configuration for lay11,12 full PS
		      : 0.f;//precise hits otherwise
		    const float zGeom1 = mockMode == 3 ? pixelPSZpitch : 0.f;
		    
		    const float kZ = (zOut - zIn)/dzSDIn;
		    float drtErr = zGeom1*zGeom1*drtSDIn*drtSDIn/dzSDIn/dzSDIn * (1.f - 2.f*kZ + 2.f*kZ*kZ);//Notes:122316
		    drtErr += sdlMuls*sdlMuls*multDzDr*multDzDr/3.f*coshEta*coshEta;//sloppy: relative muls is 1/3 of total muls
		    drtErr = sqrt(drtErr);
		    const float drtMean = drtSDIn*dzOutInAbs/std::abs(dzSDIn); //
		    const float rtWindow = drtErr + rtGeom1;
		    const float rtLo = rtIn + drtMean/dzDrtScale - rtWindow;
		    const float rtHi = rtIn + drtMean + rtWindow;
		    if (kZ > 0 && rtOut > rtLo && rtOut < rtHi) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		    if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZPointed)){
		      std::cout<<"rtPointing failed: tpPt "<<tpPt<<" lIn "<<lIn <<" rtLo "<<rtLo<<" rtHi "<<rtHi<<" vs rtOut "<<rtOut
			       <<" : drtErr "<<drtErr<<" (rtOut - rtIn) "<<(rtOut - rtIn)<<std::endl;
		      std::cout<<"\t\t RefL "<<sdIn.mdRef.pixL<<" RefU "<<sdIn.mdRef.pixU
			       <<" OutL "<<sdIn.mdOut.pixL<<" OutU "<<sdIn.mdOut.pixU<<std::endl;
		    }
		  }//else {//endcap
		}
		else if (lIn>=11 &&
			 (sdIn.mdRef.rt < disks2SMinRadius && sdIn.mdOut.rt < disks2SMinRadius)//FIXME: make configurable for full inner disks PS
			 ){//can point to the r pos in lOut;
		  assert(lOut >= 13);//only endcaps to match to
		  const float coshEta = dr3SDIn/drtSDIn;//direction estimate
		  const float dzOutInAbs = std::abs(zOut - zIn);
		  const float multDzDr = dzOutInAbs*coshEta/(coshEta*coshEta - 1.f);
		  const float rtGeom = mockMode == 3 ?
		    (rtIn < disks2SMinRadius && rtOut < disks2SMinRadius ? 2.f*pixelPSZpitch
		     : (rtIn < disks2SMinRadius || rtOut < disks2SMinRadius ) ? (pixelPSZpitch+strip2SZpitch)
		     : 2.f*strip2SZpitch) //FIXME: make this chosen by configuration for lay11,12 full PS
		    : 0.f;//precise hits otherwise
		  
		  float drtErr = pixelPSZpitch*pixelPSZpitch*2.f/dzSDIn/dzSDIn*dzOutInAbs*dzOutInAbs;//both sides contribute to direction uncertainty
		  drtErr += sdlMuls*sdlMuls*multDzDr*multDzDr/3.f*coshEta*coshEta;//sloppy: relative muls is 1/3 of total muls
		  drtErr = sqrt(drtErr);
		  const float drtMean = drtSDIn*dzOutInAbs/std::abs(dzSDIn);
		  const float rtWindow = drtErr + rtGeom; //
		  const float rtLo = rtIn + drtMean/dzDrtScale - rtWindow;
		  const float rtHi = rtIn + drtMean + rtWindow;
		  if (rtOut > rtLo && rtOut < rtHi) sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		  
		  if (debugSimMatching && !(sdlFlag & 1 << SDLSelectFlags::deltaZPointed)){
		    std::cout<<"rtPointing failed: tpPt "<<tpPt<<" lIn "<<lIn <<" rtLo "<<rtLo<<" rtHi "<<rtHi<<" vs rtOut "<<rtOut
			     <<" : drtErr "<<drtErr<<" (rtOut - rtIn) "<<(rtOut - rtIn)<<std::endl;
		    std::cout<<"\t\t RefL "<<sdIn.mdRef.pixL<<" RefU "<<sdIn.mdRef.pixU
			     <<" OutL "<<sdIn.mdOut.pixL<<" OutU "<<sdIn.mdOut.pixU<<std::endl;
		  }
		}
		else {
		  //the flag is set to pass here
		  sdlFlag |= 1 << SDLSelectFlags::deltaZPointed;
		}

		const float sdlPVoff = 0.1f/rt;
		const float sdlCut = sdlSlope + sqrt(sdlMuls*sdlMuls + sdlPVoff*sdlPVoff);

		if (lIn == 0){
		  sdlFlag |= 1 << SDLSelectFlags::deltaPhiPos;
		} else {
		  //FIXME: can be tighter
		  if (! (std::abs(deltaPhi(sdIn.mdOut.phi, sdOut.mdOut.phi)) > sdlCut) ) sdlFlag |= 1 << SDLSelectFlags::deltaPhiPos;
		}
		
		auto const midR3 = 0.5f*(sdIn.r3 + sdOut.r3);
		const float dPhi = midR3.DeltaPhi(sdOut.r3 - sdIn.r3);

		if (std::abs(dPhi) < sdlCut ) sdlFlag |= 1 << SDLSelectFlags::slope;

		float betaIn;
		float betaOut;
		float betaInRHmin;
		float betaInRHmax;
		float betaOutRHmin;
		float betaOutRHmax;

		auto const dr3 = mockMode == 0 ? sdOut.r3 - sdIn.r3 : sdOut.mdOut.r3 - sdIn.r3;
		if (mockMode == 0){
		  betaIn = sdIn.alpha - sdIn.r3.DeltaPhi(dr3);
		  betaOut = - sdOut.alpha + sdOut.r3.DeltaPhi(dr3); //to match sign for correct match	      

		  betaInRHmin = betaIn;
		  betaInRHmax = betaIn;
		  betaOutRHmin = betaOut;
		  betaOutRHmax = betaOut;
		}
		else if (mockMode == 1 || mockMode == 3){
		  //plain segment-level definitions; uneven rotation corrections are applied later using a better estimate of pt
		  if (lIn == 0){
		    betaIn  = -sdIn.p3.DeltaPhi(dr3);
		    betaOut = -sdOut.alphaOut + sdOut.mdOut.r3.DeltaPhi(dr3);
		  } else {
		    //need a symmetric choice of end-points to achieve partial cancelation
		    betaIn  = sdIn.alpha - sdIn.r3.DeltaPhi(dr3);
		    betaOut = -sdOut.alphaOut + sdOut.mdOut.r3.DeltaPhi(dr3);
		  }
		  if (mockMode == 1 || lOut < 11){
		    betaInRHmin = betaIn;
		    betaInRHmax = betaIn;
		    betaOutRHmin = betaOut;
		    betaOutRHmax = betaOut;
		  } else {
		    //FIXME: dr3 part should be varied as well
		    if (lIn < 11){
		      betaInRHmin = betaIn;
		      betaInRHmax = betaIn;		  
		    } else {
		      betaInRHmin = betaIn + sdIn.alphaRHmin - sdIn.alpha;
		      betaInRHmax = betaIn + sdIn.alphaRHmax - sdIn.alpha;
		      if (std::abs(betaInRHmin) > std::abs(betaInRHmax)) std::swap(betaInRHmax, betaInRHmin);
		    }
		    if (lOut < 11){
		      betaOutRHmin = betaOut;
		      betaOutRHmax = betaOut;		  
		    } else {
		      betaOutRHmin = betaOut - sdOut.alphaOutRHmin + sdOut.alphaOut;
		      betaOutRHmax = betaOut - sdOut.alphaOutRHmax + sdOut.alphaOut;
		      if (std::abs(betaOutRHmin) > std::abs(betaOutRHmax)) std::swap(betaOutRHmax, betaOutRHmin);
		    }
		  }
		}
		else{
		  betaIn = -99;
		  betaOut = 999;
		  betaInRHmin = betaIn;
		  betaInRHmax = betaIn;
		  betaOutRHmin = betaOut;
		  betaOutRHmax = betaOut;
		}

		if (debugSimMatching){
		  std::cout<<"beta initial: tpPt "<<tpPt<<" eta "<<tpEta<<" bI "<<betaIn<<" bO "<<betaOut<<std::endl;
		}
		const float dr = dr3.Perp();
		//beta upper cuts: 2-strip difference for direction resolution
		const float corrF = mockMode == 0 ? 0.f : 1.f;
		bool pass_betaIn_cut = lIn == 0;//pixel seeds were already selected
		if (lIn != 0){
		  const float betaIn_cut = (-sdIn.dr*corrF + dr)*k2Rinv1GeVf/ptCut + (mockMode == 3 ? 0.02f/sdIn.d : 0.f);
		  pass_betaIn_cut = std::abs(betaInRHmin) < betaIn_cut;
		}
		if (pass_betaIn_cut) sdlFlag |=  1 << SDLSelectFlags::dAlphaIn;
		
		//now the actual segment linking magic
		float betaAv = 0.5f*(betaIn + betaOut);
		//pt/k2Rinv1GeVf/2. = R
		//R*sin(betaAv) = pt/k2Rinv1GeVf/2*sin(betaAv) = dr/2 => pt = dr*k2Rinv1GeVf/sin(betaAv);
		float pt_beta = dr*k2Rinv1GeVf/sin(betaAv);
		if (lIn == 0) pt_beta = ptIn;

		const float pt_betaMax = 7.0f;

		//apply segment (SD) bend correction
		if (mockMode == 1 || mockMode == 3){
		  if (lIn == 0){
		    betaOut += copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaOut);//FIXME: need a faster version
		  } else {
		    const float diffDr = std::abs(sdIn.dr - sdOut.dr)/std::abs(sdIn.dr + sdOut.dr);
		    if (true //do it for all//diffDr > 0.05 //only if segment length is different significantly
			&& betaIn*betaOut > 0.f
			&& (std::abs(pt_beta) < 4.f*pt_betaMax
			    || (lIn >= 11 && std::abs(pt_beta) < 8.f*pt_betaMax) )){ //and the pt_beta is well-defined; less strict for endcap-endcap
		      const float betaInUpd  = betaIn + copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		      const float betaOutUpd = betaOut + copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaOut);//FIXME: need a faster version
		      if (debugSimMatching){
			std::cout<<"beta update Opt1: tpPt "<<tpPt<<" eta "<<tpEta<<" bIu "<<betaInUpd<<" bOu "<<betaOutUpd<<std::endl;
		      }
		      betaAv = 0.5f*(betaInUpd + betaOutUpd);
		      pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate
		      betaIn  += copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		      betaOut += copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaOut);//FIXME: need a faster version
		      //update the av and pt
		      betaAv = 0.5f*(betaIn + betaOut);
		      pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate		  
		    } else if (lIn < 11 && std::abs(betaOut) < 0.2f* std::abs(betaIn) && std::abs(pt_beta) < 12.f*pt_betaMax){//use betaIn sign as ref
		      const float pt_betaIn = dr*k2Rinv1GeVf/sin(betaIn);
		      const float betaInUpd  = betaIn + copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_betaIn), sinAlphaMax)), betaIn);//FIXME: need a faster version
		      const float betaOutUpd = betaOut + copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_betaIn), sinAlphaMax)), betaIn);//FIXME: need a faster version
		      if (debugSimMatching){
			std::cout<<"beta update Opt2: tpPt "<<tpPt<<" eta "<<tpEta<<" bIu "<<betaInUpd<<" bOu "<<betaOutUpd<<std::endl;
		      }
		      betaAv = std::abs(betaOut) > 0.2f*std::abs(betaIn) ? 0.5f*(betaInUpd + betaOutUpd) : betaInUpd;
		      pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate
		      betaIn  += copysign(std::asin(std::min(sdIn.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		      betaOut += copysign(std::asin(std::min(sdOut.dr*k2Rinv1GeVf/std::abs(pt_beta), sinAlphaMax)), betaIn);//FIXME: need a faster version
		      //update the av and pt
		      betaAv = 0.5f*(betaIn + betaOut);
		      pt_beta = dr*k2Rinv1GeVf/sin(betaAv);//get a better pt estimate		  		      
		    }
		  }
		}
		if (debugSimMatching){
		  std::cout<<"beta corrF: tpPt "<<tpPt<<" eta "<<tpEta<<" bI "<<betaIn<<" bO "<<betaOut<<std::endl;
		}
		//rescale the ranges proportionally
		const float betaInMMSF = 2.f*betaIn/std::abs(betaInRHmin+betaInRHmax);//mean value of min,max is the old betaIn
		const float betaOutMMSF = 2.f*betaOut/std::abs(betaOutRHmin+betaOutRHmax);
		if (debugSimMatching){
		  std::cout<<"Scale betaMM. Before "<<betaInRHmin<<" "<<betaInRHmax<<" "<<betaOutRHmin<<" "<<betaOutRHmax
			   <<" SF "<<betaInMMSF<<" "<<betaOutMMSF <<std::endl;
		}
		betaInRHmin*= betaInMMSF;
		betaInRHmax*= betaInMMSF;
		betaOutRHmin*= betaOutMMSF;
		betaOutRHmax*= betaOutMMSF;
		
		const float dBetaMuls = sdlThetaMulsF*4.f/std::min(std::abs(pt_beta), pt_betaMax); //need to confirm the range-out value of 7 GeV

		//regularize to alpha of pt_betaMax .. review may want to add resolution
		const float alphaInAbsReg = useFullR3Endcap ? 0.f : std::max(std::abs(sdIn.alpha), std::asin(std::min(sdIn.rt*k2Rinv1GeVf/3.0f, sinAlphaMax)));
		const float alphaOutAbsReg = useFullR3Endcap ? 0.f : std::max(std::abs(sdOut.alpha), std::asin(std::min(sdOut.rt*k2Rinv1GeVf/3.0f, sinAlphaMax)));
		const float dBetaInLum = lIn < 11 ? 0.0f : std::abs(alphaInAbsReg*deltaZLum/sdIn.z);
		const float dBetaOutLum = lOut < 11 ? 0.0f : std::abs(alphaOutAbsReg*deltaZLum/sdOut.z);
		const float dBetaLum2 = (dBetaInLum + dBetaOutLum)*(dBetaInLum + dBetaOutLum);

		const float sinDPhi = std::sin(dPhi);
		const float dBetaRIn2 = std::pow((sdIn.mdRef.rtRHout - sdIn.mdRef.rtRHin)*sinDPhi/dr, 2);
		const float dBetaROut2 = std::pow((sdOut.mdOut.rtRHout - sdOut.mdOut.rtRHin)*sinDPhi/dr, 2);

		const float betaOut_cut = std::asin(std::min(dr*k2Rinv1GeVf/ptCut, sinAlphaMax))//FIXME: need faster version
		  + (mockMode == 3 ? 0.02f/sdOut.d : 0.f) + sqrt(dBetaLum2 + dBetaMuls*dBetaMuls);
		if (std::abs(betaOut) < betaOut_cut) sdlFlag |=  1 << SDLSelectFlags::dAlphaOut;

		float pt_betaIn = dr*k2Rinv1GeVf/sin(betaIn);
		if (lIn == 0) pt_betaIn = pt_beta;
		const float pt_betaOut = dr*k2Rinv1GeVf/sin(betaOut);
		
		const float dBetaRes = mockMode == 3 ? 0.02f/std::min(sdIn.d, sdOut.d) : 0.f;

		const float dBetaCut2 = (dBetaRes*dBetaRes*2.0f + dBetaMuls*dBetaMuls + dBetaLum2 + dBetaRIn2 + dBetaROut2
					 + 0.25*std::pow(std::abs(betaInRHmin - betaInRHmax) + std::abs(betaOutRHmin - betaOutRHmax),2));

		const float dBeta = betaIn - betaOut;
		
		if (dBeta*dBeta < dBetaCut2) sdlFlag |= 1 << SDLSelectFlags::dBeta;
		if (debugSimMatching ){
		  if (std::abs(betaOut) > betaOut_cut){
		    std::cout<<"betaOut failed pt "<<tpPt<<" eta "<<tpEta
			     <<" bo "<<betaOut<<" boCut "<<betaOut_cut<<" bOutPt "<<std::asin(std::min(dr*k2Rinv1GeVf/ptCut, sinAlphaMax))
			     <<" dr "<<dr<<" ptCut "<<ptCut<<" ptIn "<<ptIn<<" ptErr "<<see_ptErr()[sdIn.iRef]
			     <<" bPt "<<pt_beta<<" dr "<<dr<<" bAv "<<betaAv<<" bI "<<betaIn<<" bO "<<betaOut
			     <<" dRI "<<sdIn.dr<<" dRO "<<sdOut.dr
			     <<" aI "<<sdIn.alpha<<" aO "<<sdOut.alpha
			     <<std::endl;
		    auto const& simhitIdxV = sim_simHitIdx()[iSim];
		    for (auto ish : simhitIdxV){
		      SimHit simH(ish);
		      int lay = simH.lay;
		      
		      std::cout<<" "<<lay<<" "<<ish<<std::endl;		      
		      if (simH.p3s.Pt()>0.8*tpPt){
			simH.print("\t");
		      }
		    }//iPix; pixelhits

		  }
		  if (dBeta*dBeta >= dBetaCut2 || debugAllDbeta){
		    std::cout<<"dB failed pt "<<tpPt<<" "<<tpEta<<" "<<tpPhi<<" rz "<<tpDxy<<" "<<tpDz
			     <<" tpProd "<<prodX<<" "<<prodY<<" "<<prodZ
			     <<" dbCut "<<sqrt(dBetaCut2)<<" res "<<dBetaRes<<" muls "<<dBetaMuls		    
			     <<" mulPt "<<std::min(std::abs(pt_beta), pt_betaMax)<<" thetaM "<<sdlThetaMulsF
			     <<" dbLum "<<sqrt(dBetaLum2)<<" dbILum "<<dBetaInLum<<" dbOLum "<<dBetaOutLum
			     <<" dBInMM "<<std::abs(betaInRHmin - betaInRHmax)<<" dBOutMM "<<std::abs(betaOutRHmin - betaOutRHmax)
			     <<" bPt "<<pt_beta<<" dr "<<dr<<" bAv "<<betaAv<<" bI "<<betaIn<<" bO "<<betaOut
			     <<" dRI "<<sdIn.dr<<" dRO "<<sdOut.dr
			     <<" aI "<<sdIn.alpha<<" aO "<<sdOut.alpha
			     <<" aIO "<<sdIn.alphaOut<<" aOO "<<sdOut.alphaOut
			     <<std::endl;
		  }
		}
		// if (lIn == 0 && lOut == 7 && tpPt > 1.5f && tpPt < 2.0f){
		//   std::cout<<__LINE__<<" 5-7: tPt "<<tpPt <<"bO "<<betaOut<<" bOc "<<betaOut_cut<<" = "<<std::asin(std::min(dr*k2Rinv1GeVf/ptCut, sinAlphaMax)<<" + "<<0.02f/sdOut.d
		// 	   <<" dr "<<dr<<" ptCut "<<ptCut<<" ptIn "<<ptIn
		// 	   <<" bI "<<betaIn<<" dB "<<dBeta<<" dBc "<<sqrt(dBetaCut2)
		// 	   <<std::endl;
		// }
		
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
	      }//sdOut : vSDOut_4of4
	    }//sdIn : vSDIn_4of4
	    if (runDetailedTimers) timerA[T_timeVal_SDL].Stop();


	    const int m_0 = 1;
	    const int m_01 = m_0 | (1 << 1);
	    const int m_012 = m_01 | (1 << 2);
	    const int m_0123 = m_012 | (1 << 3);
	    const int m_01234 = m_0123 | (1 << 4);
	    const int m_012345 = m_01234 | (1 << 5);
	    const int m_0123456 = m_012345 | (1 << 6);

	    bool h_0, h_01, h_012, h_0123, h_01234, h_012345, h_0123456;
	    h_0 = h_01 = h_012 = h_0123 = h_01234 = h_012345 = h_0123456 = false;
	    bool debugMatchingSim = false;
	    for (auto const& sdlf : vSDLwInfo_4of4){
	      if ((sdlf.second & m_0) == m_0 ) h_0 = true;
	      if ((sdlf.second & m_01) == m_01 ) h_01 = true;
	      if ((sdlf.second & m_012) == m_012 ) h_012 = true;
	      if ((sdlf.second & m_0123) == m_0123 ) h_0123 = true;
	      if ((sdlf.second & m_01234) == m_01234 ) h_01234 = true;
	      if ((sdlf.second & m_012345) == m_012345 ) h_012345 = true;
	      if ((sdlf.second & m_0123456) == m_0123456 ) h_0123456 = true;

	      
	      auto const& sdl = sdlf.first;
	      auto dBeta = sdl.betaIn - sdl.betaOut;

	      if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4){
		if (tpPt > 0.7 && tpPt < 1.0) ha_SDL_dBeta_zoom_8MH_pt0p7to1p0[iSDLL]->Fill(dBeta);
		if (tpPt > 1.0 && tpPt < 1.2) ha_SDL_dBeta_zoom_8MH_pt1p0to1p2[iSDLL]->Fill(dBeta);
		if (tpPt > 1.2 && tpPt < 1.5) ha_SDL_dBeta_zoom_8MH_pt1p2to1p5[iSDLL]->Fill(dBeta);
		if (tpPt > 1.5 && tpPt < 2.0) ha_SDL_dBeta_zoom_8MH_pt1p5to2p0[iSDLL]->Fill(dBeta);
		if (tpPt > 2.0 && tpPt < 4.0) ha_SDL_dBeta_zoom_8MH_pt2p0to4p0[iSDLL]->Fill(dBeta);
		if (tpPt > 4.0 && tpPt < 7.0) ha_SDL_dBeta_zoom_8MH_pt4p0to7p0[iSDLL]->Fill(dBeta);
		if (tpPt > 7.0              ) ha_SDL_dBeta_zoom_8MH_pt7p0toInf[iSDLL]->Fill(dBeta);

		if (tpPt > 1.0 && tpPt < 1.2 && (lIn >= 11 || lOut >= 11) && debugMatchingSim){
		  std::cout<<lIn<<" "<<lOut<<" pt "<<tpPt<<" db "<<dBeta<<" bi "<<sdl.betaIn<<" bo "<<sdl.betaOut
			   <<" ai "<<sdl.sdIn.alpha<<" ao "<<sdl.sdOut.alpha
			   <<" aip "<<sdl.sdIn.mdRef.r3.DeltaPhi(sdl.sdIn.mdOut.r3 - sdl.sdIn.mdRef.r3)
			   <<" aop "<<sdl.sdOut.mdRef.r3.DeltaPhi(sdl.sdOut.mdOut.r3 - sdl.sdOut.mdRef.r3)
			   <<" dpi "<<deltaPhi(sdl.sdIn.mdRef.phi, sdl.sdIn.mdOut.phi)
			   <<" dpo "<<deltaPhi(sdl.sdOut.mdRef.phi, sdl.sdOut.mdOut.phi)
			   <<" rii "<<sdl.sdIn.mdRef.rt<<" rio "<<sdl.sdIn.mdOut.rt
			   <<" roi "<<sdl.sdOut.mdRef.rt<<" roo "<<sdl.sdOut.mdOut.rt
			   <<" pii "<<sdl.sdIn.mdRef.phi<<" pio "<<sdl.sdIn.mdOut.phi
			   <<" poi "<<sdl.sdOut.mdRef.phi<<" poo "<<sdl.sdOut.mdOut.phi
			   <<" zii "<<sdl.sdIn.mdRef.r3.z()<<" zio "<<sdl.sdIn.mdOut.r3.z()
			   <<" zoi "<<sdl.sdOut.mdRef.r3.z()<<" zoo "<<sdl.sdOut.mdOut.r3.z()
			   <<std::endl;
		}
	      }
	      if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4 && h_012345){
		ha_SDL_dBeta_zoom_NM1dBeta_8MH[iSDLL]->Fill(dBeta);
		ha_SDL_dBeta_zoom2_NM1dBeta_8MH[iSDLL]->Fill(dBeta);
		ha_SDL_dBeta_betaIn_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn, dBeta);
		ha_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn, dBeta);
		ha_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaOut, dBeta);
		
		//		if (std::abs(tpPt-18.0172) < 0.0001 && iSDLL == SDL_L5to11 && !h_0123456){
		if (debugSimMatching && (!h_0123456 || debugAllDbeta)){ //&& !h_0123456 && std::abs(sdl.betaIn - sdl.betaOut)> 0.01){
		  std::cout<<"dB-fail tPt "<<tpPt
			   <<" tEta " << tpEta<<" tPhi "<<tpPhi
			   <<" dB "<<sdl.betaIn - sdl.betaOut
			   <<" bI "<<sdl.betaIn
			   <<" bO "<<sdl.betaOut
			   <<" sI "<<sdl.sdIn.alpha
			   <<" sO "<<sdl.sdOut.alpha
			   <<" mI "<<sdl.sdIn.mdRef.alpha
			   <<" mO "<<sdl.sdOut.mdRef.alpha
			   <<" mRI ("<<sdl.sdIn.mdRef.r3.Pt()<<", "<<sdl.sdIn.mdRef.r3.Eta()<<", "<<sdl.sdIn.mdRef.r3.Phi()<<", "<<sdl.sdIn.mdRef.r3.Z()<<")"
			   <<" mOI ("<<sdl.sdIn.mdOut.r3.Pt()<<", "<<sdl.sdIn.mdOut.r3.Eta()<<", "<<sdl.sdIn.mdOut.r3.Phi()<<", "<<sdl.sdIn.mdOut.r3.Z()<<")"
			   <<" mRO ("<<sdl.sdOut.mdRef.r3.Pt()<<", "<<sdl.sdOut.mdRef.r3.Eta()<<", "<<sdl.sdOut.mdRef.r3.Phi()<<", "<<sdl.sdOut.mdRef.r3.Z() <<")"
			   <<" mOO ("<<sdl.sdOut.mdOut.r3.Pt()<<", "<<sdl.sdOut.mdOut.r3.Eta()<<", "<<sdl.sdOut.mdOut.r3.Phi()<<", "<<sdl.sdOut.mdOut.r3.Z() <<")"
			   <<std::endl;

		  std::vector<TVector3> r3Ins;
		  std::vector<TVector3> r3Outs;
		  std::vector<TVector3> p3Ins;
		  std::vector<TVector3> p3Outs;
		  

		  auto const& simhitIdxV = sim_simHitIdx()[iSim];
		  for (auto ish : simhitIdxV){
		    
		    SimHit simH(ish);
		    int lay = simH.lay;
		    
		    // int iProcess = simH.iProcess;
		    // int ibx = simH.bx;	  
		    // bool isPrimaryAny = (iProcess == 2 && ibx == 0);
		    
		    if (lay >= minLayer){
		      std::cout<<" "<<lay<<" "<<ish<<std::endl;
		      
		      if (simH.p3s.Pt()>0.8*tpPt){
			if (lay == layersSDL[iSDLL][0] || lay == layersSDL[iSDLL][0]+1){
			  r3Ins.push_back(simH.r3s);
			  p3Ins.push_back(simH.p3s);
			} else if (lay == layersSDL[iSDLL][1] || lay == layersSDL[iSDLL][1]+1){
			  r3Outs.push_back(simH.r3s);
			  p3Outs.push_back(simH.p3s);			    
			}
			simH.print("\t");
		      }
		    }		    
		  }//simhitIdxV

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
	      if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4 && h_0123456){
		auto const& sdl = sdlf.first;
		ha_SDL_dBeta_betaIn_pass_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn, sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaIn_zoom_pass_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaIn, sdl.betaIn - sdl.betaOut);
		ha_SDL_dBeta_betaOut_zoom_pass_NM1dBeta_8MH[iSDLL]->Fill(sdl.betaOut, sdl.betaIn - sdl.betaOut);
	      }

	    }

	    if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4){
	      if (h_0) ha_num2SD_w0_4of4_pt[iSDLL]->Fill(tpPt);
	      if (h_01) ha_num2SD_w01_4of4_pt[iSDLL]->Fill(tpPt);
	      if (h_012) ha_num2SD_w012_4of4_pt[iSDLL]->Fill(tpPt);
	      if (h_0123) ha_num2SD_w0123_4of4_pt[iSDLL]->Fill(tpPt);
	      if (h_01234) ha_num2SD_w01234_4of4_pt[iSDLL]->Fill(tpPt);
	      if (h_012345) ha_num2SD_w012345_4of4_pt[iSDLL]->Fill(tpPt);
	      if (h_0123456) ha_num2SD_w0123456_4of4_pt[iSDLL]->Fill(tpPt);
	      
	      if (isLowProdXY){
		if (h_0) ha_num2SD_w0_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
		if (h_01) ha_num2SD_w01_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
		if (h_012) ha_num2SD_w012_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
		if (h_0123) ha_num2SD_w0123_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
		if (h_01234) ha_num2SD_w01234_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
		if (h_012345) ha_num2SD_w012345_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
		if (h_0123456) ha_num2SD_w0123456_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
	      }
	    }
	    
	    if (runDetailedTimers) timerA[T_timeVal_SDLMatch].Start(kFALSE);
	    if (!mockLayerSDLsDNcm[iSDLL].empty()){
	      for (auto& sdl : mockLayerSDLsDNcm[iSDLL]){
		if (! (sdl.lIn == lIn && sdl.lOut == lOut )) continue;

		int scoreIn = sdl.sdIn.itp == iSim ? sdl.sdIn.ntp : 0;
		if (debug && scoreIn > 1){
		  std::cout<<"Inner match on L"<< sdl.lIn<<": Have "<<scoreIn<<" matches "<<std::endl;
		}
		int scoreOut = sdl.sdOut.itp == iSim ? sdl.sdOut.ntp : 0;
		if (debug && scoreOut > 1){
		  std::cout<<"Outer match on L"<< sdl.lOut<<": Have "<<scoreOut<<" matches "<<std::endl;
		}

		if (scoreIn >=3 && scoreOut>= 3){
		  matchingSDLs_byHit3of4[iSDLL].push_back(sdl);
		  if (scoreIn >=4 && scoreOut>= 4){
		    matchingSDLs_byHit4of4[iSDLL].push_back(sdl);
		    sdl.hasMatch_byHit4of4 = true;
		  }
		}
	      }//for (auto sdl : mockLayerSDLsDNcm[iSDLL]){
	    }//	if (!mockLayerSDLsDNcm[iSDLL].empty()){
	    if (runDetailedTimers) timerA[T_timeVal_SDLMatch].Stop();

	    bool hasMatch = false;
	    //matching is done: fill numerators
	    if (has8MHs && has4MDs && hasSDIn_3of4 && hasSDOut_3of4 && ! matchingSDLs_byHit3of4[iSDLL].empty()){
	      ha_numSDL_3of4_any_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_numSDL_3of4_any_loProdXY_pt[iSDLL]->Fill(tpPt);
	      if (debug) std::cout<<"\t have 3/4 match "<<std::endl;
	      hasMatch = true;
	    }
	    if (has8MHs && has4MDs && hasSDIn_4of4 && hasSDOut_4of4 && ! matchingSDLs_byHit4of4[iSDLL].empty()){
	      ha_numSDL_4of4_pt[iSDLL]->Fill(tpPt);
	      if (isLowProdXY) ha_numSDL_4of4_loProdXY_pt[iSDLL]->Fill(tpPt);
	      if (debug) std::cout<<"\t have 4/4 match "<<std::endl;
	      hasMatch = true;
	    }
	    
	    if (debug && ! hasMatch && iSDLL != SDL_L5to9){
	      for (auto i : simHits[lIn]){ auto ph = SimHit(i); ph.print("\tNM for: ");}
	      for (auto i : simHits[lOut]){ auto ph = SimHit(i); ph.print("\tNM for: ");}
	    }

	  }// TP has hits in SDL layers	 	  
	}//for (int iSDLL = 0; iSDLL< SDL_LMAX; ++iSDLL){
      }//TPs
	

      //fake rates
      for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
	auto const& mSDLs = mockLayerSDLsDNcm[iSDL];
	if (mSDLs.empty() ) continue;
	for (auto const& sdl : mSDLs){
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
      std::cout<<"SimHit stats:";
      for (int l = 0; l< nLayersA+1;++l) std::cout<<" L"<<l<<" "<<(float)nHitsStatSumMap[l]/std::max(1, nHitsStatCntMap[l]);
      std::cout<<std::endl;
      
      std::cout<<"Print stats"<<std::endl;
      for (int iL = minLayer; iL <= nLayersA; ++iL){
	int nSDs1GeV = 0;
	int nSDs1GeVMatch4 = 0;
	int nSDs2GeV = 0;
	int nSDs2GeVMatch4 = 0;
	for (auto sd : mockLayerSDfwDNcm[iL]){
	  int ish = simHitsPerHitAll(sd.mdRef.pixL);
	  TVector3 p3RL(simhit_px()[ish], simhit_py()[ish], simhit_pz()[ish]);
	  ish = simHitsPerHitAll(sd.mdRef.pixU);
	  TVector3 p3RU(simhit_px()[ish], simhit_py()[ish], simhit_pz()[ish]);
	  ish = simHitsPerHitAll(sd.mdOut.pixL);
	  TVector3 p3OL(simhit_px()[ish], simhit_py()[ish], simhit_pz()[ish]);
	  ish = simHitsPerHitAll(sd.mdOut.pixU);
	  TVector3 p3OU(simhit_px()[ish], simhit_py()[ish], simhit_pz()[ish]);

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
	int nSDNGs = 0;
	for (unsigned int i = 0; i< mockLayerSDfwDNcm_isSecondaryGhost[iL].size(); ++i) if (not mockLayerSDfwDNcm_isSecondaryGhost[iL][i]) nSDNGs++;
	std::cout<<"Summary for layer "<<iL
		 <<" h1GeV "<<nHitsLayer1GeV[iL]
		 <<" h2GeV "<<nHitsLayer2GeV[iL]
		 <<" refLos "<<mockLayerMDfwRefLower[iL].size()
		 <<" refUps "<<mockLayerMDfwRefUpper[iL].size()
		 <<" refMDs " <<mockLayerMDfwRef[iL].size()
		 <<" outLos " <<mockLayerMDfwDNcmLower[iL].size()
		 <<" outUps " <<mockLayerMDfwDNcmUpper[iL].size()
		 <<" outMDs " <<mockLayerMDfwDNcm[iL].size()
		 <<" SDs "<< mockLayerSDfwDNcm[iL].size()
		 <<" SDNGs "<< nSDNGs
		 <<" n1Any "<<nSDs1GeV<<" n1All "<<nSDs1GeVMatch4
		 <<" n2Any "<<nSDs2GeV<<" n2All "<<nSDs2GeVMatch4
		 <<std::endl;
      }//iL
      std::cout<<"\t\tsdl5-7 "<<mockLayerSDLsDNcm[SDL_L5to7].size()
	       <<" sdl7-9 "<<mockLayerSDLsDNcm[SDL_L7to9].size()
	       <<" sdl5-11 "<<mockLayerSDLsDNcm[SDL_L5to11].size()
	       <<" sdl7-11 "<<mockLayerSDLsDNcm[SDL_L7to11].size()
	       <<" sdl11-13 "<<mockLayerSDLsDNcm[SDL_L11to13].size()
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
    if (iT > T_timeValidation && ! runDetailedTimers) continue;
    std::cout<<"Timer results for "<<timerNameA[iT]<<std::endl;
    timerA[iT].Print("m");
    std::cout<<"----------------------------------"<<std::endl;
    
  }
  
  std::cout<<__LINE__<<" make efficiencies "<<std::endl;
  std::array<TEfficiency*, SDL_LMAX> ha_eff8MH_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff4MD_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_3of4_any_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_3of4_any_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_3of4_any_pt;

  std::array<TEfficiency*, SDL_LMAX> ha_eff4MD_den8MH_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_den8MH_3of4_any_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_den8MH_3of4_any_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w01_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w012_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0123_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w01234_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w012345_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0123456_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_den8MH_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_den8MH_3of4_any_pt;

  std::array<TEfficiency*, SDL_LMAX> ha_eff8MH_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff4MD_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_3of4_any_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_3of4_any_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_3of4_any_loProdXY_pt;

  std::array<TEfficiency*, SDL_LMAX> ha_eff4MD_den8MH_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effInSD_den8MH_3of4_any_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_den8MH_3of4_any_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w01_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w012_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_den8MH_4of4_loProdXY_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_effSDL_den8MH_3of4_any_loProdXY_pt;

  std::array<TEfficiency*, SDL_LMAX> ha_fakeSDL_4of4_pt;
  std::array<TEfficiency*, SDL_LMAX> ha_fakeSDL_4of4_eta;
  for (int i = 0; i< SDL_LMAX; ++i){
    auto iMin = layersSDL[i][0];
    auto iMax = layersSDL[i][1];
    ha_eff8MH_pt[i] = new TEfficiency(*ha_num8MH_pt[i], *ha_denSDL_pt[i]);
    ha_eff8MH_pt[i]->SetTitle(Form("h_eff8MH_%dto%d_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_eff4MD_pt[i] = new TEfficiency(*ha_num4MD_pt[i], *ha_denSDL_pt[i]);
    ha_eff4MD_pt[i]->SetTitle(Form("h_eff4MD_%dto%d_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_4of4_pt[i] = new TEfficiency(*ha_numInSD_4of4_pt[i], *ha_denSDL_pt[i]);
    ha_effInSD_4of4_pt[i]->SetTitle(Form("h_effInSD_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_3of4_any_pt[i] = new TEfficiency(*ha_numInSD_3of4_any_pt[i], *ha_denSDL_pt[i]);
    ha_effInSD_3of4_any_pt[i]->SetTitle(Form("h_effInSD_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
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
    ha_effInSD_den8MH_4of4_pt[i] = new TEfficiency(*ha_numInSD_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_effInSD_den8MH_4of4_pt[i]->SetTitle(Form("h_effInSD_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_den8MH_3of4_any_pt[i] = new TEfficiency(*ha_numInSD_3of4_any_pt[i], *ha_num8MH_pt[i]);
    ha_effInSD_den8MH_3of4_any_pt[i]->SetTitle(Form("h_effInSD_den8MH_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
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
    ha_eff2SD_w0123456_den8MH_4of4_pt[i] = new TEfficiency(*ha_num2SD_w0123456_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_eff2SD_w0123456_den8MH_4of4_pt[i]->SetTitle(Form("h_eff2SD_w0123456_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_den8MH_4of4_pt[i] = new TEfficiency(*ha_numSDL_4of4_pt[i], *ha_num8MH_pt[i]);
    ha_effSDL_den8MH_4of4_pt[i]->SetTitle(Form("h_effSDL_den8MH_%dto%d_4of4_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_den8MH_3of4_any_pt[i] = new TEfficiency(*ha_numSDL_3of4_any_pt[i], *ha_num8MH_pt[i]);
    ha_effSDL_den8MH_3of4_any_pt[i]->SetTitle(Form("h_effSDL_den8MH_%dto%d_3of4_any_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_eff8MH_loProdXY_pt[i] = new TEfficiency(*ha_num8MH_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_eff8MH_loProdXY_pt[i]->SetTitle(Form("h_eff8MH_%dto%d_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_eff4MD_loProdXY_pt[i] = new TEfficiency(*ha_num4MD_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_eff4MD_loProdXY_pt[i]->SetTitle(Form("h_eff4MD_%dto%d_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_4of4_loProdXY_pt[i] = new TEfficiency(*ha_numInSD_4of4_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_effInSD_4of4_loProdXY_pt[i]->SetTitle(Form("h_effInSD_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_3of4_any_loProdXY_pt[i] = new TEfficiency(*ha_numInSD_3of4_any_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_effInSD_3of4_any_loProdXY_pt[i]->SetTitle(Form("h_effInSD_%dto%d_3of4_any_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_4of4_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_eff2SD_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_3of4_any_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_3of4_any_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_eff2SD_3of4_any_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_%dto%d_3of4_any_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_4of4_loProdXY_pt[i] = new TEfficiency(*ha_numSDL_4of4_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_effSDL_4of4_loProdXY_pt[i]->SetTitle(Form("h_effSDL_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_3of4_any_loProdXY_pt[i] = new TEfficiency(*ha_numSDL_3of4_any_loProdXY_pt[i], *ha_denSDL_loProdXY_pt[i]);
    ha_effSDL_3of4_any_loProdXY_pt[i]->SetTitle(Form("h_effSDL_%dto%d_3of4_any_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_eff4MD_den8MH_loProdXY_pt[i] = new TEfficiency(*ha_num4MD_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff4MD_den8MH_loProdXY_pt[i]->SetTitle(Form("h_eff4MD_den8MH_%dto%d_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_numInSD_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_effInSD_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_effInSD_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effInSD_den8MH_3of4_any_loProdXY_pt[i] = new TEfficiency(*ha_numInSD_3of4_any_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_effInSD_den8MH_3of4_any_loProdXY_pt[i]->SetTitle(Form("h_effInSD_den8MH_%dto%d_3of4_any_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_den8MH_3of4_any_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_3of4_any_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_den8MH_3of4_any_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_den8MH_%dto%d_3of4_any_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w0_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w0_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w01_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w01_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w012_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w012_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w0123_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w0123_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w01234_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w01234_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w012345_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w012345_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_num2SD_w0123456_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_eff2SD_w0123456_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_den8MH_4of4_loProdXY_pt[i] = new TEfficiency(*ha_numSDL_4of4_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_effSDL_den8MH_4of4_loProdXY_pt[i]->SetTitle(Form("h_effSDL_den8MH_%dto%d_4of4_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );
    ha_effSDL_den8MH_3of4_any_loProdXY_pt[i] = new TEfficiency(*ha_numSDL_3of4_any_loProdXY_pt[i], *ha_num8MH_loProdXY_pt[i]);
    ha_effSDL_den8MH_3of4_any_loProdXY_pt[i]->SetTitle(Form("h_effSDL_den8MH_%dto%d_3of4_any_loProdXY_pt; p_{T} (GeV); Efficiency", iMin, iMax) );

    ha_fakeSDL_4of4_pt[i] = new TEfficiency(*ha_SDLreco_no4of4_pt[i], *ha_SDLreco_all_pt[i]);
    ha_fakeSDL_4of4_pt[i]->SetTitle(Form("h_fakeSDL_%dto%d_4of4_pt", iMin, iMax));

    ha_fakeSDL_4of4_eta[i] = new TEfficiency(*ha_SDLreco_no4of4_eta[i], *ha_SDLreco_all_eta[i]);
    ha_fakeSDL_4of4_eta[i]->SetTitle(Form("h_fakeSDL_%dto%d_4of4_eta", iMin, iMax));
  }


  const float xAxMinPt = ptCutAll >=1 ? 0.51 : 0.31;

  if (drawPlots){
    {
      auto h2 = h2_hitsXY_ITrec_OTmockLL;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      h2->SetTitle("Hits (x,y);x (cm);y (cm)");
      h2->SetStats(0);
      h2->Draw();
      auto ax = h2->GetYaxis();
      ax->SetTitleOffset(ax->GetTitleOffset()+0.25);
      gPad->SetLogx(0);
      gPad->SaveAs(Form("h2_hitsXY_ITrec_OTmockLL_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    {
      auto h2 = h2_hitsRZ_ITrec_OTmockLL;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      h2->SetTitle("Hits (r,z);z (cm);r (cm)");
      h2->SetStats(0);
      h2->Draw();
      auto ax = h2->GetYaxis();
      ax->SetTitleOffset(ax->GetTitleOffset()+0.25);
      gPad->SetLogx(0);
      gPad->SaveAs(Form("h2_hitsRZ_ITrec_OTmockLL_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    {
      auto h2 = h2_hitsRZ_ITrec_OTmockLL_BE;
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 1200, 600);
      cv->cd();
      h2->SetTitle("Hits (r,z);z (cm);r (cm)");
      h2->SetStats(0);
      h2->Draw();
      auto ax = h2->GetYaxis();
      ax->SetTitleOffset(ax->GetTitleOffset()+0.25);
      gPad->SetLogx(0);
      gPad->SaveAs(Form("h2_hitsRZ_ITrec_OTmockLL_BE_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
  }
  
  if (drawPlots && !layoutOnly){
    std::cout<<__LINE__<<" draw and print "<<std::endl;
    for (int iL = 5; iL <= nLayersA; ++iL){
      if (iL != 5 && iL != 10) continue;
      
    }//nLayersA

    
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_NM1dBeta_all[iSDL];
      auto h_pass = ha_SDL_dBeta_NM1dBeta_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_NM1dBeta_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_NM1dBeta_all[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_NM1dBeta_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_NM1dBeta_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_NM1dBeta_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto h_4 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL];
      auto h_3 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL];
      auto h_0to2 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_4->SetLineWidth(0);
      h_3->SetLineWidth(0);
      h_0to2->SetLineWidth(0);
      h_4->SetFillColor(kRed);
      h_3->SetFillColor(kGreen);
      h_0to2->SetFillColor(kBlue);
      
      THStack* h_st = new THStack("h_st", "");
      h_st->Add(h_4);
      h_st->Add(h_3);
      h_st->Add(h_0to2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_st->Draw("same");
      h_all->Draw("AXIS same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_prov_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy();
      h_all->SetMinimum(0.5);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_prov_log_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy(0);
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto h_4 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL];
      auto h_3 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL];
      auto h_0to2 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_4->SetLineWidth(0);
      h_3->SetLineWidth(0);
      h_0to2->SetLineWidth(0);
      h_4->SetFillColor(kRed);
      h_3->SetFillColor(kGreen);
      h_0to2->SetFillColor(kBlue);
      
      THStack* h_st = new THStack("h_st", "");
      h_st->Add(h_0to2);
      h_st->Add(h_3);
      h_st->Add(h_4);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_st->Draw("same");
      h_all->Draw("AXIS same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_prov2_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy();
      h_all->SetMinimum(0.5);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_prov2_log_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy(0);
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto h_4 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL];
      auto h_3 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL];
      auto h_0to2 = ha_SDL_dBeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_4->SetLineWidth(0);
      h_3->SetLineWidth(2);
      h_0to2->SetLineWidth(0);
      h_3->SetLineColor(kRed);
      h_4->SetFillColor(kWhite);
      h_3->SetFillColor(kWhite);
      h_0to2->SetFillColor(kWhite);
      
      THStack* h_st = new THStack("h_st", "");
      h_st->Add(h_0to2);
      h_st->Add(h_3);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_st->Draw("same");
      h_all->Draw("same");
      h_all->Draw("AXIS same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_prov3_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy();
      h_all->SetMinimum(0.5);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_NGLL_prov3_log_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy(0);
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto h_4 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL];
      auto h_3 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL];
      auto h_0to2 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_4->SetLineWidth(0);
      h_3->SetLineWidth(0);
      h_0to2->SetLineWidth(0);
      h_4->SetFillColor(kRed);
      h_3->SetFillColor(kGreen);
      h_0to2->SetFillColor(kBlue);
      
      THStack* h_st = new THStack("h_st", "");
      h_st->Add(h_4);
      h_st->Add(h_3);
      h_st->Add(h_0to2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_st->Draw("same");
      h_all->Draw("AXIS same");
      gPad->SaveAs(Form("h_SDL_dZeta_zoom_NM1dBeta_NGLL_prov_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy();
      h_all->SetMinimum(0.5);
      gPad->SaveAs(Form("h_SDL_dZeta_zoom_NM1dBeta_NGLL_prov_log_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy(0);
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto h_4 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL];
      auto h_3 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL];
      auto h_0to2 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_4->SetLineWidth(0);
      h_3->SetLineWidth(0);
      h_0to2->SetLineWidth(0);
      h_4->SetFillColor(kRed);
      h_3->SetFillColor(kGreen);
      h_0to2->SetFillColor(kBlue);
      
      THStack* h_st = new THStack("h_st", "");
      h_st->Add(h_0to2);
      h_st->Add(h_3);
      h_st->Add(h_4);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_st->Draw("same");
      h_all->Draw("AXIS same");
      gPad->SaveAs(Form("h_SDL_dZeta_zoom_NM1dBeta_NGLL_prov2_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy();
      h_all->SetMinimum(0.5);
      gPad->SaveAs(Form("h_SDL_dZeta_zoom_NM1dBeta_NGLL_prov2_log_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy(0);
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL[iSDL];
      auto h_4 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL4[iSDL];
      auto h_3 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL3[iSDL];
      auto h_0to2 = ha_SDL_dZeta_zoom_NM1dBeta_all_NGLL_LL0to2[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_4->SetLineWidth(0);
      h_3->SetLineWidth(2);
      h_0to2->SetLineWidth(0);
      h_3->SetLineColor(kRed);
      h_4->SetFillColor(kWhite);
      h_3->SetFillColor(kWhite);
      h_0to2->SetFillColor(kWhite);
      
      THStack* h_st = new THStack("h_st", "");
      h_st->Add(h_0to2);
      h_st->Add(h_3);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_st->Draw("same");
      h_all->Draw("same");
      h_all->Draw("AXIS same");
      gPad->SaveAs(Form("h_SDL_dZeta_zoom_NM1dBeta_NGLL_prov3_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy();
      h_all->SetMinimum(0.5);
      gPad->SaveAs(Form("h_SDL_dZeta_zoom_NM1dBeta_NGLL_prov3_log_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      gPad->SetLogy(0);
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_all[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    // << slices in betaIn or ptIn
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //all
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_all[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_all[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_all[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //no overlay, pass case with pt slices
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn0to2_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_ptIn0to2_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn3to5_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_ptIn3to5_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_NM1dBeta_ptIn7toInf_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_ptIn7toInf_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //other steps in selections
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_0_all[iSDL];
      auto h_pass = ha_SDL_dBeta_0_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_0_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_0_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_0_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_0_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_01_all[iSDL];
      auto h_pass = ha_SDL_dBeta_01_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_01_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_01_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_01_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_01_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_012_all[iSDL];
      auto h_pass = ha_SDL_dBeta_012_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_012_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_012_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_012_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_012_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_0123_all[iSDL];
      auto h_pass = ha_SDL_dBeta_0123_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_0123_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_0123_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_0123_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_0123_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_01234_all[iSDL];
      auto h_pass = ha_SDL_dBeta_01234_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_01234_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_01234_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_01234_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_01234_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_012345_all[iSDL];
      auto h_pass = ha_SDL_dBeta_012345_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_012345_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h_all = ha_SDL_dBeta_zoom_012345_all[iSDL];
      auto h_pass = ha_SDL_dBeta_zoom_012345_pass[iSDL];
      auto cn = h_all->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h_all->SetLineWidth(2);
      h_pass->SetLineWidth(2);
      h_pass->SetLineColor(2);

      h_all->SetStats(0);
      h_all->Draw();
      h_all->SetMinimum(0);
      h_pass->Draw("same");
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_012345_all_vs_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }


    // >> in pt slices
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_NM1dBeta_all[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_NM1dBeta_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
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
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
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
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom2_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    //pt slices
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt0p7to1p0[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt0p7to1p0_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt1p0to1p2[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt1p0to1p2_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt1p2to1p5[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt1p2to1p5_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt1p5to2p0[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt1p5to2p0_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt2p0to4p0[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt2p0to4p0_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt4p0to7p0[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt4p0to7p0_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h = ha_SDL_dBeta_zoom_8MH_pt7p0toInf[iSDL];
      auto cn = h->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->Draw();
      h->SetMinimum(0.5);
      gPad->SetGridx();
      gPad->SetLogy();
      gPad->SaveAs(Form("h_SDL_dBeta_zoom_8MH_pt7p0toInf_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    
    //2D plots
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
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
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
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_zoom_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
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
      gPad->SaveAs(Form("h2_SDL_dBeta_betaOut_zoom_NM1dBeta_8MH_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }

    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_pass_NM1dBeta_8MH[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaIn_zoom_pass_NM1dBeta_8MH[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaIn_zoom_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
    }
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      auto h2 = ha_SDL_dBeta_betaOut_zoom_pass_NM1dBeta_8MH[iSDL];
      auto cn = h2->GetTitle();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();
      gPad->SetRightMargin(gPad->GetRightMargin()*1.1);
      h2->SetStats(0);
      h2->Draw("colz");
      gPad->SetGridx();
      gPad->SetLogx(0);
      gPad->SetLogz();
      gPad->SaveAs(Form("h2_SDL_dBeta_betaOut_zoom_pass_%dto%d_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
      
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
      ax->SetRangeUser(xAxMinPt, ax->GetXmax());
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
      gPad->SaveAs(Form("h_denVSnums_SDL_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
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
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.7, 0.15, 0.9, 0.3);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff1, "3 of 4");
      leg->AddEntry(heff2, "4 of 4");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
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
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.6, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heffH, "8 MH match");
      leg->AddEntry(heff, "4 MD match");
      leg->AddEntry(heff1, "2 SD match");
      leg->AddEntry(heff2, "SDLink match");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_steps_4of4_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
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
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_4of4_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps; detailed0
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_pt[iSDL];
      auto heff1in = ha_effInSD_den8MH_4of4_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1in->Draw("same");
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1in->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1in->SetLineColor(kMagenta);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1in->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1in->SetMarkerStyle(20);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1in, "Inner SD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_detailed0_4of4_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
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
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff012->Draw("same");
      heff0123->Draw("same");
      heff01234->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012->SetLineWidth(2);
      heff0123->SetLineWidth(2);
      heff01234->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012->SetLineColor(kMagenta);
      heff0123->SetLineColor(kGray);
      heff01234->SetLineColor(kGreen);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012->SetMarkerSize(0.8);
      heff0123->SetMarkerSize(0.8);
      heff01234->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012->SetMarkerStyle(26);
      heff0123->SetMarkerStyle(27);
      heff01234->SetMarkerStyle(28);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.6);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012, "... and #Delta#phi_{pos} | 8 MHits");
      leg->AddEntry(heff0123, "... and slope | 8 MHits");
      leg->AddEntry(heff01234, "... and #Delta#alpha_{in} | 8 MHits");      
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_detailed_4of4_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
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
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.9);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_4of4_min0.9_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_pt[iSDL];
      auto heff1in = ha_effInSD_den8MH_4of4_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1in->Draw("same");
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1in->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1in->SetLineColor(kMagenta);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1in->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1in->SetMarkerStyle(20);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.9);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1in, "Inner SD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_detailed0_4of4_min0.9_%dto%d_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //same set with low production DXY
    //efficiencies: num/den plots
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto hden = ha_denSDL_loProdXY_pt[iSDL];
      auto hnum1 = ha_numSDL_3of4_any_loProdXY_pt[iSDL];
      auto hnum2 = ha_numSDL_4of4_loProdXY_pt[iSDL];

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
      ax->SetRangeUser(xAxMinPt, ax->GetXmax());
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
      gPad->SaveAs(Form("h_denVSnums_SDL_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 3/4 vs 4/4
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff1 = ha_effSDL_3of4_any_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_4of4_loProdXY_pt[iSDL];

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
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.7, 0.15, 0.9, 0.3);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff1, "3 of 4");
      leg->AddEntry(heff2, "4 of 4");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    
    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heffH =  ha_eff8MH_loProdXY_pt[iSDL];
      auto heff =  ha_eff4MD_loProdXY_pt[iSDL];
      auto heff1 = ha_eff2SD_4of4_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_4of4_loProdXY_pt[iSDL];

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
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.6, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heffH, "8 MH match");
      leg->AddEntry(heff, "4 MD match");
      leg->AddEntry(heff1, "2 SD match");
      leg->AddEntry(heff2, "SDLink match");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_steps_4of4_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps den8MH
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_loProdXY_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_loProdXY_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_4of4_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps den8MH
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_loProdXY_pt[iSDL];
      auto heff1in = ha_effInSD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_loProdXY_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1in->Draw("same");
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1in->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1in->SetLineColor(kMagenta);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1in->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1in->SetMarkerStyle(20);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1in, "Inner SD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_detailed0_4of4_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_loProdXY_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_loProdXY_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff012->Draw("same");
      heff0123->Draw("same");
      heff01234->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012->SetLineWidth(2);
      heff0123->SetLineWidth(2);
      heff01234->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012->SetLineColor(kMagenta);
      heff0123->SetLineColor(kGray);
      heff01234->SetLineColor(kGreen);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012->SetMarkerSize(0.8);
      heff0123->SetMarkerSize(0.8);
      heff01234->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012->SetMarkerStyle(26);
      heff0123->SetMarkerStyle(27);
      heff01234->SetMarkerStyle(28);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.6);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012, "... and #Delta#phi_{pos} | 8 MHits");
      leg->AddEntry(heff0123, "... and slope | 8 MHits");
      leg->AddEntry(heff01234, "... and #Delta#alpha_{in} | 8 MHits");      
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_detailed_4of4_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_loProdXY_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_loProdXY_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.9);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_4of4_min0.9_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //efficiencies: eff plots 4/4 in steps
    for (int iSDL = 0; iSDL < SDL_LMAX; ++iSDL){
      if (iSDL == SDL_L5to9) continue;
      
      auto heff =  ha_eff4MD_den8MH_loProdXY_pt[iSDL];
      auto heff1in = ha_effInSD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff1 = ha_eff2SD_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0      = ha_eff2SD_w0_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01     = ha_eff2SD_w01_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012    = ha_eff2SD_w012_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123   = ha_eff2SD_w0123_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff01234  = ha_eff2SD_w01234_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff012345 = ha_eff2SD_w012345_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff0123456 = ha_eff2SD_w0123456_den8MH_4of4_loProdXY_pt[iSDL];
      auto heff2 = ha_effSDL_den8MH_4of4_loProdXY_pt[iSDL];

      auto cn = heff->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      heff->Draw();
      heff1in->Draw("same");
      heff1->Draw("same");
      heff01->Draw("same");
      heff012345->Draw("same");
      heff2->Draw("same");

      heff->SetLineWidth(2);
      heff1in->SetLineWidth(2);
      heff1->SetLineWidth(2);
      heff01->SetLineWidth(2);
      heff012345->SetLineWidth(2);
      heff2->SetLineWidth(2);
      heff->SetLineColor(kBlack);
      heff1in->SetLineColor(kMagenta);
      heff1->SetLineColor(kRed);
      heff01->SetLineColor(kOrange);
      heff012345->SetLineColor(kCyan);
      heff2->SetLineColor(kBlue);
      heff->SetMarkerSize(0.8);
      heff1in->SetMarkerSize(0.8);
      heff1->SetMarkerSize(0.8);
      heff01->SetMarkerSize(0.8);
      heff012345->SetMarkerSize(0.8);
      heff2->SetMarkerSize(0.8);
      heff->SetMarkerStyle(21);
      heff1in->SetMarkerStyle(22);
      heff1->SetMarkerStyle(22);
      heff01->SetMarkerStyle(24);
      heff012345->SetMarkerStyle(25);
      heff2->SetMarkerStyle(23);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = heff->GetPaintedGraph();
      pg->SetMinimum(0.9);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(xAxMinPt, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.87, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(heff, "4 MD match | 8 MHits");
      leg->AddEntry(heff1in, "Inner SD match | 8 MHits");
      leg->AddEntry(heff1, "2 SD match | 8 MHits");
      leg->AddEntry(heff01, "2 SD w dZ | 8 MHits");
      leg->AddEntry(heff012345, "SDLink no #Delta#beta | 8 MHits");
      leg->AddEntry(heff2, "SDLink match | 8 MHits");
      leg->Draw();
      gPad->SaveAs(Form("h_effs_SDL_den8MH_steps_detailed0_4of4_min0.9_%dto%d_loProdXY_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			layersSDL[iSDL][0], layersSDL[iSDL][1], mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //fakes vs pt
    {      
      auto h05 =  ha_fakeSDL_4of4_pt[SDL_L0to5];
      auto h07 =  ha_fakeSDL_4of4_pt[SDL_L0to7];
      
      auto cn = h05->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h05->Draw();
      h07->Draw("same");

      h05->SetLineWidth(2);
      h07->SetLineWidth(2);
      h05->SetLineColor(kBlack);
      h07->SetLineColor(kRed);
      h05->SetMarkerSize(0.8);
      h07->SetMarkerSize(0.8);
      h05->SetMarkerStyle(21);
      h07->SetMarkerStyle(22);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogx();
      gPad->PaintModified();
      auto pg = h05->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(1.0, ax->GetXmax());

      auto leg = new TLegend(0.5, 0.15, 0.86, 0.22);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->AddEntry(h05, "L0-L5", "LP");
      leg->AddEntry(h07, "L0-L7", "LP");
      leg->Draw();
      gPad->SaveAs(Form("h_fake_SDL_0-5_0-7_4of4_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    //fakes vs eta
    {      
      auto h05 =  ha_fakeSDL_4of4_eta[SDL_L0to5];
      auto h07 =  ha_fakeSDL_4of4_eta[SDL_L0to7];
      
      auto cn = h05->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h05->Draw();
      h07->Draw("same");

      h05->SetLineWidth(2);
      h07->SetLineWidth(2);
      h05->SetLineColor(kBlack);
      h07->SetLineColor(kRed);
      h05->SetMarkerSize(0.8);
      h07->SetMarkerSize(0.8);
      h05->SetMarkerStyle(21);
      h07->SetMarkerStyle(22);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->PaintModified();
      auto pg = h05->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);

      auto leg = new TLegend(0.5, 0.15, 0.86, 0.22);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->AddEntry(h05, "L0-L5", "LP");
      leg->AddEntry(h07, "L0-L7", "LP");
      leg->Draw();
      gPad->SaveAs(Form("h_fake_SDL_0-5_0-7_4of4_eta2_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    {      
      auto h05 =  ha_fakeSDL_4of4_eta[SDL_L0to5];
      auto h07 =  ha_fakeSDL_4of4_eta[SDL_L0to7];
      
      auto cn = h05->GetName();
      TCanvas* cv = new TCanvas(cn, cn, 600, 600);
      cv->cd();

      h05->Draw();
      h07->Draw("same");

      h05->SetLineWidth(2);
      h07->SetLineWidth(2);
      h05->SetLineColor(kBlack);
      h07->SetLineColor(kRed);
      h05->SetMarkerSize(0.8);
      h07->SetMarkerSize(0.8);
      h05->SetMarkerStyle(21);
      h07->SetMarkerStyle(22);

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->PaintModified();
      auto pg = h05->GetPaintedGraph();
      pg->SetMinimum(0.0);
      pg->SetMaximum(1.02);
      auto ax = pg->GetXaxis();
      ax->SetLimits(-1.5, 1.5);

      auto leg = new TLegend(0.5, 0.15, 0.86, 0.22);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->AddEntry(h05, "L0-L5", "LP");
      leg->AddEntry(h07, "L0-L7", "LP");
      leg->Draw();
      gPad->SaveAs(Form("h_fake_SDL_0-5_0-7_4of4_eta_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

    //fakes vs pt 5-7 and 7-9
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
      gPad->SaveAs(Form("h_fake_SDL_5-7_7-9_4of4_pt_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
    //fakes vs eta
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

      auto leg = new TLegend(0.5, 0.15, 0.86, 0.22);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->AddEntry(h57, "L5-L7", "LP");
      leg->AddEntry(h79, "L7-L9", "LP");
      leg->Draw();
      gPad->SaveAs(Form("h_fake_SDL_5-7_7-9_4of4_eta2_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }
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
      gPad->SaveAs(Form("h_fake_SDL_5-7_7-9_4of4_eta_mm%d_D%1.1fcm%1.1fcm_us%d.png",
			mockMode, sdOffsetB, sdOffsetE,
			useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks));
    }

  }//if drawPlots
  
  std::cout<<__LINE__<<" write to file "<<std::endl;
  TFile* outHistograms = new TFile(Form("outHistogramsSuperD_mm%d_D%1.1fcm%1.1fcm_us%d.root",
					mockMode, sdOffsetB, sdOffsetE,
					useSeeds + 10*useFullR3Endcap + 100*effForPromptTracks), "RECREATE");
  for (int iL = 1; iL <= nLayersA; ++iL){
    layerMD_pt_all[iL].write(outHistograms);
    layerMD_pt_prim_all[iL].write(outHistograms);
    layerMD_pt_prim_tt[iL].write(outHistograms);

    h2_hitsXY_ITrec_OTmockLL->Write();
    h2_hitsRZ_ITrec_OTmockLL->Write();
  }
  for (auto h : outputHV){
    if (h != nullptr) h->Write();
  }
  outHistograms->Write();
  outHistograms->Close();
  std::cout<<__LINE__<<" Done "<<std::endl;

  return 0;
}

