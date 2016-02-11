// -*- C++ -*-
#ifndef tkph2_H
#define tkph2_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class tkph2 {
private: 
protected: 
	unsigned int index;
	vector<float> *trk_px_;
	TBranch *trk_px_branch;
	bool trk_px_isLoaded;
	vector<float> *trk_py_;
	TBranch *trk_py_branch;
	bool trk_py_isLoaded;
	vector<float> *trk_pz_;
	TBranch *trk_pz_branch;
	bool trk_pz_isLoaded;
	vector<float> *trk_pt_;
	TBranch *trk_pt_branch;
	bool trk_pt_isLoaded;
	vector<float> *trk_eta_;
	TBranch *trk_eta_branch;
	bool trk_eta_isLoaded;
	vector<float> *trk_phi_;
	TBranch *trk_phi_branch;
	bool trk_phi_isLoaded;
	vector<float> *trk_dxy_;
	TBranch *trk_dxy_branch;
	bool trk_dxy_isLoaded;
	vector<float> *trk_dz_;
	TBranch *trk_dz_branch;
	bool trk_dz_isLoaded;
	vector<float> *trk_ptErr_;
	TBranch *trk_ptErr_branch;
	bool trk_ptErr_isLoaded;
	vector<float> *trk_etaErr_;
	TBranch *trk_etaErr_branch;
	bool trk_etaErr_isLoaded;
	vector<float> *trk_phiErr_;
	TBranch *trk_phiErr_branch;
	bool trk_phiErr_isLoaded;
	vector<float> *trk_dxyErr_;
	TBranch *trk_dxyErr_branch;
	bool trk_dxyErr_isLoaded;
	vector<float> *trk_dzErr_;
	TBranch *trk_dzErr_branch;
	bool trk_dzErr_isLoaded;
	vector<float> *trk_nChi2_;
	TBranch *trk_nChi2_branch;
	bool trk_nChi2_isLoaded;
	vector<float> *trk_shareFrac_;
	TBranch *trk_shareFrac_branch;
	bool trk_shareFrac_isLoaded;
	vector<int> *trk_q_;
	TBranch *trk_q_branch;
	bool trk_q_isLoaded;
	vector<int> *trk_nValid_;
	TBranch *trk_nValid_branch;
	bool trk_nValid_isLoaded;
	vector<int> *trk_nInvalid_;
	TBranch *trk_nInvalid_branch;
	bool trk_nInvalid_isLoaded;
	vector<int> *trk_nPixel_;
	TBranch *trk_nPixel_branch;
	bool trk_nPixel_isLoaded;
	vector<int> *trk_nStrip_;
	TBranch *trk_nStrip_branch;
	bool trk_nStrip_isLoaded;
	vector<int> *trk_n3DLay_;
	TBranch *trk_n3DLay_branch;
	bool trk_n3DLay_isLoaded;
	vector<int> *trk_algo_;
	TBranch *trk_algo_branch;
	bool trk_algo_isLoaded;
	vector<int> *trk_isHP_;
	TBranch *trk_isHP_branch;
	bool trk_isHP_isLoaded;
	vector<int> *trk_seedIdx_;
	TBranch *trk_seedIdx_branch;
	bool trk_seedIdx_isLoaded;
	vector<int> *trk_simIdx_;
	TBranch *trk_simIdx_branch;
	bool trk_simIdx_isLoaded;
	vector<vector<int> > *trk_pixelIdx_;
	TBranch *trk_pixelIdx_branch;
	bool trk_pixelIdx_isLoaded;
	vector<vector<int> > *trk_stripIdx_;
	TBranch *trk_stripIdx_branch;
	bool trk_stripIdx_isLoaded;
	vector<float> *sim_px_;
	TBranch *sim_px_branch;
	bool sim_px_isLoaded;
	vector<float> *sim_py_;
	TBranch *sim_py_branch;
	bool sim_py_isLoaded;
	vector<float> *sim_pz_;
	TBranch *sim_pz_branch;
	bool sim_pz_isLoaded;
	vector<float> *sim_pt_;
	TBranch *sim_pt_branch;
	bool sim_pt_isLoaded;
	vector<float> *sim_eta_;
	TBranch *sim_eta_branch;
	bool sim_eta_isLoaded;
	vector<float> *sim_phi_;
	TBranch *sim_phi_branch;
	bool sim_phi_isLoaded;
	vector<float> *sim_dxy_;
	TBranch *sim_dxy_branch;
	bool sim_dxy_isLoaded;
	vector<float> *sim_dz_;
	TBranch *sim_dz_branch;
	bool sim_dz_isLoaded;
	vector<float> *sim_prodx_;
	TBranch *sim_prodx_branch;
	bool sim_prodx_isLoaded;
	vector<float> *sim_prody_;
	TBranch *sim_prody_branch;
	bool sim_prody_isLoaded;
	vector<float> *sim_prodz_;
	TBranch *sim_prodz_branch;
	bool sim_prodz_isLoaded;
	vector<float> *sim_shareFrac_;
	TBranch *sim_shareFrac_branch;
	bool sim_shareFrac_isLoaded;
	vector<int> *sim_q_;
	TBranch *sim_q_branch;
	bool sim_q_isLoaded;
	vector<int> *sim_nValid_;
	TBranch *sim_nValid_branch;
	bool sim_nValid_isLoaded;
	vector<int> *sim_nPixel_;
	TBranch *sim_nPixel_branch;
	bool sim_nPixel_isLoaded;
	vector<int> *sim_nStrip_;
	TBranch *sim_nStrip_branch;
	bool sim_nStrip_isLoaded;
	vector<int> *sim_n3DLay_;
	TBranch *sim_n3DLay_branch;
	bool sim_n3DLay_isLoaded;
	vector<int> *sim_trkIdx_;
	TBranch *sim_trkIdx_branch;
	bool sim_trkIdx_isLoaded;
	vector<vector<int> > *sim_pixelIdx_;
	TBranch *sim_pixelIdx_branch;
	bool sim_pixelIdx_isLoaded;
	vector<vector<int> > *sim_stripIdx_;
	TBranch *sim_stripIdx_branch;
	bool sim_stripIdx_isLoaded;
	vector<int> *pix_isBarrel_;
	TBranch *pix_isBarrel_branch;
	bool pix_isBarrel_isLoaded;
	vector<int> *pix_lay_;
	TBranch *pix_lay_branch;
	bool pix_lay_isLoaded;
	vector<int> *pix_detId_;
	TBranch *pix_detId_branch;
	bool pix_detId_isLoaded;
	vector<int> *pix_nSimTrk_;
	TBranch *pix_nSimTrk_branch;
	bool pix_nSimTrk_isLoaded;
	vector<int> *pix_simTrkIdx_;
	TBranch *pix_simTrkIdx_branch;
	bool pix_simTrkIdx_isLoaded;
	vector<int> *pix_particle_;
	TBranch *pix_particle_branch;
	bool pix_particle_isLoaded;
	vector<int> *pix_process_;
	TBranch *pix_process_branch;
	bool pix_process_isLoaded;
	vector<int> *pix_bunchXing_;
	TBranch *pix_bunchXing_branch;
	bool pix_bunchXing_isLoaded;
	vector<int> *pix_event_;
	TBranch *pix_event_branch;
	bool pix_event_isLoaded;
	vector<float> *pix_x_;
	TBranch *pix_x_branch;
	bool pix_x_isLoaded;
	vector<float> *pix_y_;
	TBranch *pix_y_branch;
	bool pix_y_isLoaded;
	vector<float> *pix_z_;
	TBranch *pix_z_branch;
	bool pix_z_isLoaded;
	vector<float> *pix_xx_;
	TBranch *pix_xx_branch;
	bool pix_xx_isLoaded;
	vector<float> *pix_xy_;
	TBranch *pix_xy_branch;
	bool pix_xy_isLoaded;
	vector<float> *pix_yy_;
	TBranch *pix_yy_branch;
	bool pix_yy_isLoaded;
	vector<float> *pix_yz_;
	TBranch *pix_yz_branch;
	bool pix_yz_isLoaded;
	vector<float> *pix_zz_;
	TBranch *pix_zz_branch;
	bool pix_zz_isLoaded;
	vector<float> *pix_zx_;
	TBranch *pix_zx_branch;
	bool pix_zx_isLoaded;
	vector<float> *pix_xsim_;
	TBranch *pix_xsim_branch;
	bool pix_xsim_isLoaded;
	vector<float> *pix_ysim_;
	TBranch *pix_ysim_branch;
	bool pix_ysim_isLoaded;
	vector<float> *pix_zsim_;
	TBranch *pix_zsim_branch;
	bool pix_zsim_isLoaded;
	vector<float> *pix_pxsim_;
	TBranch *pix_pxsim_branch;
	bool pix_pxsim_isLoaded;
	vector<float> *pix_pysim_;
	TBranch *pix_pysim_branch;
	bool pix_pysim_isLoaded;
	vector<float> *pix_pzsim_;
	TBranch *pix_pzsim_branch;
	bool pix_pzsim_isLoaded;
	vector<float> *pix_pathprop_;
	TBranch *pix_pathprop_branch;
	bool pix_pathprop_isLoaded;
	vector<float> *pix_xsimprop_;
	TBranch *pix_xsimprop_branch;
	bool pix_xsimprop_isLoaded;
	vector<float> *pix_ysimprop_;
	TBranch *pix_ysimprop_branch;
	bool pix_ysimprop_isLoaded;
	vector<float> *pix_zsimprop_;
	TBranch *pix_zsimprop_branch;
	bool pix_zsimprop_isLoaded;
	vector<float> *pix_pxsimprop_;
	TBranch *pix_pxsimprop_branch;
	bool pix_pxsimprop_isLoaded;
	vector<float> *pix_pysimprop_;
	TBranch *pix_pysimprop_branch;
	bool pix_pysimprop_isLoaded;
	vector<float> *pix_pzsimprop_;
	TBranch *pix_pzsimprop_branch;
	bool pix_pzsimprop_isLoaded;
	vector<float> *pix_eloss_;
	TBranch *pix_eloss_branch;
	bool pix_eloss_isLoaded;
	vector<float> *pix_radL_;
	TBranch *pix_radL_branch;
	bool pix_radL_isLoaded;
	vector<float> *pix_bbxi_;
	TBranch *pix_bbxi_branch;
	bool pix_bbxi_isLoaded;
	vector<int> *str_isBarrel_;
	TBranch *str_isBarrel_branch;
	bool str_isBarrel_isLoaded;
	vector<int> *str_isStereo_;
	TBranch *str_isStereo_branch;
	bool str_isStereo_isLoaded;
	vector<int> *str_det_;
	TBranch *str_det_branch;
	bool str_det_isLoaded;
	vector<int> *str_lay_;
	TBranch *str_lay_branch;
	bool str_lay_isLoaded;
	vector<int> *str_detId_;
	TBranch *str_detId_branch;
	bool str_detId_isLoaded;
	vector<int> *str_nSimTrk_;
	TBranch *str_nSimTrk_branch;
	bool str_nSimTrk_isLoaded;
	vector<int> *str_simTrkIdx_;
	TBranch *str_simTrkIdx_branch;
	bool str_simTrkIdx_isLoaded;
	vector<int> *str_particle_;
	TBranch *str_particle_branch;
	bool str_particle_isLoaded;
	vector<int> *str_process_;
	TBranch *str_process_branch;
	bool str_process_isLoaded;
	vector<int> *str_bunchXing_;
	TBranch *str_bunchXing_branch;
	bool str_bunchXing_isLoaded;
	vector<int> *str_event_;
	TBranch *str_event_branch;
	bool str_event_isLoaded;
	vector<float> *str_x_;
	TBranch *str_x_branch;
	bool str_x_isLoaded;
	vector<float> *str_y_;
	TBranch *str_y_branch;
	bool str_y_isLoaded;
	vector<float> *str_z_;
	TBranch *str_z_branch;
	bool str_z_isLoaded;
	vector<float> *str_xx_;
	TBranch *str_xx_branch;
	bool str_xx_isLoaded;
	vector<float> *str_xy_;
	TBranch *str_xy_branch;
	bool str_xy_isLoaded;
	vector<float> *str_yy_;
	TBranch *str_yy_branch;
	bool str_yy_isLoaded;
	vector<float> *str_yz_;
	TBranch *str_yz_branch;
	bool str_yz_isLoaded;
	vector<float> *str_zz_;
	TBranch *str_zz_branch;
	bool str_zz_isLoaded;
	vector<float> *str_zx_;
	TBranch *str_zx_branch;
	bool str_zx_isLoaded;
	vector<float> *str_xsim_;
	TBranch *str_xsim_branch;
	bool str_xsim_isLoaded;
	vector<float> *str_ysim_;
	TBranch *str_ysim_branch;
	bool str_ysim_isLoaded;
	vector<float> *str_zsim_;
	TBranch *str_zsim_branch;
	bool str_zsim_isLoaded;
	vector<float> *str_pxsim_;
	TBranch *str_pxsim_branch;
	bool str_pxsim_isLoaded;
	vector<float> *str_pysim_;
	TBranch *str_pysim_branch;
	bool str_pysim_isLoaded;
	vector<float> *str_pzsim_;
	TBranch *str_pzsim_branch;
	bool str_pzsim_isLoaded;
	vector<float> *str_eloss_;
	TBranch *str_eloss_branch;
	bool str_eloss_isLoaded;
	vector<float> *str_radL_;
	TBranch *str_radL_branch;
	bool str_radL_isLoaded;
	vector<float> *str_bbxi_;
	TBranch *str_bbxi_branch;
	bool str_bbxi_isLoaded;
	vector<int> *glu_isBarrel_;
	TBranch *glu_isBarrel_branch;
	bool glu_isBarrel_isLoaded;
	vector<int> *glu_det_;
	TBranch *glu_det_branch;
	bool glu_det_isLoaded;
	vector<int> *glu_lay_;
	TBranch *glu_lay_branch;
	bool glu_lay_isLoaded;
	vector<int> *glu_detId_;
	TBranch *glu_detId_branch;
	bool glu_detId_isLoaded;
	vector<int> *glu_monoIdx_;
	TBranch *glu_monoIdx_branch;
	bool glu_monoIdx_isLoaded;
	vector<int> *glu_stereoIdx_;
	TBranch *glu_stereoIdx_branch;
	bool glu_stereoIdx_isLoaded;
	vector<float> *glu_x_;
	TBranch *glu_x_branch;
	bool glu_x_isLoaded;
	vector<float> *glu_y_;
	TBranch *glu_y_branch;
	bool glu_y_isLoaded;
	vector<float> *glu_z_;
	TBranch *glu_z_branch;
	bool glu_z_isLoaded;
	vector<float> *glu_xx_;
	TBranch *glu_xx_branch;
	bool glu_xx_isLoaded;
	vector<float> *glu_xy_;
	TBranch *glu_xy_branch;
	bool glu_xy_isLoaded;
	vector<float> *glu_yy_;
	TBranch *glu_yy_branch;
	bool glu_yy_isLoaded;
	vector<float> *glu_yz_;
	TBranch *glu_yz_branch;
	bool glu_yz_isLoaded;
	vector<float> *glu_zz_;
	TBranch *glu_zz_branch;
	bool glu_zz_isLoaded;
	vector<float> *glu_zx_;
	TBranch *glu_zx_branch;
	bool glu_zx_isLoaded;
	vector<float> *glu_radL_;
	TBranch *glu_radL_branch;
	bool glu_radL_isLoaded;
	vector<float> *glu_bbxi_;
	TBranch *glu_bbxi_branch;
	bool glu_bbxi_isLoaded;
	float	bsp_x_;
	TBranch *bsp_x_branch;
	bool bsp_x_isLoaded;
	float	bsp_y_;
	TBranch *bsp_y_branch;
	bool bsp_y_isLoaded;
	float	bsp_z_;
	TBranch *bsp_z_branch;
	bool bsp_z_isLoaded;
	float	bsp_sigmax_;
	TBranch *bsp_sigmax_branch;
	bool bsp_sigmax_isLoaded;
	float	bsp_sigmay_;
	TBranch *bsp_sigmay_branch;
	bool bsp_sigmay_isLoaded;
	float	bsp_sigmaz_;
	TBranch *bsp_sigmaz_branch;
	bool bsp_sigmaz_isLoaded;
	vector<float> *see_px_;
	TBranch *see_px_branch;
	bool see_px_isLoaded;
	vector<float> *see_py_;
	TBranch *see_py_branch;
	bool see_py_isLoaded;
	vector<float> *see_pz_;
	TBranch *see_pz_branch;
	bool see_pz_isLoaded;
	vector<float> *see_pt_;
	TBranch *see_pt_branch;
	bool see_pt_isLoaded;
	vector<float> *see_eta_;
	TBranch *see_eta_branch;
	bool see_eta_isLoaded;
	vector<float> *see_phi_;
	TBranch *see_phi_branch;
	bool see_phi_isLoaded;
	vector<float> *see_dxy_;
	TBranch *see_dxy_branch;
	bool see_dxy_isLoaded;
	vector<float> *see_dz_;
	TBranch *see_dz_branch;
	bool see_dz_isLoaded;
	vector<float> *see_ptErr_;
	TBranch *see_ptErr_branch;
	bool see_ptErr_isLoaded;
	vector<float> *see_etaErr_;
	TBranch *see_etaErr_branch;
	bool see_etaErr_isLoaded;
	vector<float> *see_phiErr_;
	TBranch *see_phiErr_branch;
	bool see_phiErr_isLoaded;
	vector<float> *see_dxyErr_;
	TBranch *see_dxyErr_branch;
	bool see_dxyErr_isLoaded;
	vector<float> *see_dzErr_;
	TBranch *see_dzErr_branch;
	bool see_dzErr_isLoaded;
	vector<float> *see_chi2_;
	TBranch *see_chi2_branch;
	bool see_chi2_isLoaded;
	vector<int> *see_q_;
	TBranch *see_q_branch;
	bool see_q_isLoaded;
	vector<int> *see_nValid_;
	TBranch *see_nValid_branch;
	bool see_nValid_isLoaded;
	vector<int> *see_nPixel_;
	TBranch *see_nPixel_branch;
	bool see_nPixel_isLoaded;
	vector<int> *see_nGlued_;
	TBranch *see_nGlued_branch;
	bool see_nGlued_isLoaded;
	vector<int> *see_nStrip_;
	TBranch *see_nStrip_branch;
	bool see_nStrip_isLoaded;
	vector<int> *see_algo_;
	TBranch *see_algo_branch;
	bool see_algo_isLoaded;
	vector<vector<int> > *see_pixelIdx_;
	TBranch *see_pixelIdx_branch;
	bool see_pixelIdx_isLoaded;
	vector<vector<int> > *see_gluedIdx_;
	TBranch *see_gluedIdx_branch;
	bool see_gluedIdx_isLoaded;
	vector<vector<int> > *see_stripIdx_;
	TBranch *see_stripIdx_branch;
	bool see_stripIdx_isLoaded;
	vector<int> *algo_offset_;
	TBranch *algo_offset_branch;
	bool algo_offset_isLoaded;
public: 
void Init(TTree *tree) {
  tree->SetMakeClass(1);
	trk_px_branch = 0;
	if (tree->GetBranch("trk_px") != 0) {
		trk_px_branch = tree->GetBranch("trk_px");
		trk_px_branch->SetAddress(&trk_px_);
	}
	if(trk_px_branch == 0 ) {
	cout << "Branch trk_px does not exist." << endl;
	}
	trk_py_branch = 0;
	if (tree->GetBranch("trk_py") != 0) {
		trk_py_branch = tree->GetBranch("trk_py");
		trk_py_branch->SetAddress(&trk_py_);
	}
	if(trk_py_branch == 0 ) {
	cout << "Branch trk_py does not exist." << endl;
	}
	trk_pz_branch = 0;
	if (tree->GetBranch("trk_pz") != 0) {
		trk_pz_branch = tree->GetBranch("trk_pz");
		trk_pz_branch->SetAddress(&trk_pz_);
	}
	if(trk_pz_branch == 0 ) {
	cout << "Branch trk_pz does not exist." << endl;
	}
	trk_pt_branch = 0;
	if (tree->GetBranch("trk_pt") != 0) {
		trk_pt_branch = tree->GetBranch("trk_pt");
		trk_pt_branch->SetAddress(&trk_pt_);
	}
	if(trk_pt_branch == 0 ) {
	cout << "Branch trk_pt does not exist." << endl;
	}
	trk_eta_branch = 0;
	if (tree->GetBranch("trk_eta") != 0) {
		trk_eta_branch = tree->GetBranch("trk_eta");
		trk_eta_branch->SetAddress(&trk_eta_);
	}
	if(trk_eta_branch == 0 ) {
	cout << "Branch trk_eta does not exist." << endl;
	}
	trk_phi_branch = 0;
	if (tree->GetBranch("trk_phi") != 0) {
		trk_phi_branch = tree->GetBranch("trk_phi");
		trk_phi_branch->SetAddress(&trk_phi_);
	}
	if(trk_phi_branch == 0 ) {
	cout << "Branch trk_phi does not exist." << endl;
	}
	trk_dxy_branch = 0;
	if (tree->GetBranch("trk_dxy") != 0) {
		trk_dxy_branch = tree->GetBranch("trk_dxy");
		trk_dxy_branch->SetAddress(&trk_dxy_);
	}
	if(trk_dxy_branch == 0 ) {
	cout << "Branch trk_dxy does not exist." << endl;
	}
	trk_dz_branch = 0;
	if (tree->GetBranch("trk_dz") != 0) {
		trk_dz_branch = tree->GetBranch("trk_dz");
		trk_dz_branch->SetAddress(&trk_dz_);
	}
	if(trk_dz_branch == 0 ) {
	cout << "Branch trk_dz does not exist." << endl;
	}
	trk_ptErr_branch = 0;
	if (tree->GetBranch("trk_ptErr") != 0) {
		trk_ptErr_branch = tree->GetBranch("trk_ptErr");
		trk_ptErr_branch->SetAddress(&trk_ptErr_);
	}
	if(trk_ptErr_branch == 0 ) {
	cout << "Branch trk_ptErr does not exist." << endl;
	}
	trk_etaErr_branch = 0;
	if (tree->GetBranch("trk_etaErr") != 0) {
		trk_etaErr_branch = tree->GetBranch("trk_etaErr");
		trk_etaErr_branch->SetAddress(&trk_etaErr_);
	}
	if(trk_etaErr_branch == 0 ) {
	cout << "Branch trk_etaErr does not exist." << endl;
	}
	trk_phiErr_branch = 0;
	if (tree->GetBranch("trk_phiErr") != 0) {
		trk_phiErr_branch = tree->GetBranch("trk_phiErr");
		trk_phiErr_branch->SetAddress(&trk_phiErr_);
	}
	if(trk_phiErr_branch == 0 ) {
	cout << "Branch trk_phiErr does not exist." << endl;
	}
	trk_dxyErr_branch = 0;
	if (tree->GetBranch("trk_dxyErr") != 0) {
		trk_dxyErr_branch = tree->GetBranch("trk_dxyErr");
		trk_dxyErr_branch->SetAddress(&trk_dxyErr_);
	}
	if(trk_dxyErr_branch == 0 ) {
	cout << "Branch trk_dxyErr does not exist." << endl;
	}
	trk_dzErr_branch = 0;
	if (tree->GetBranch("trk_dzErr") != 0) {
		trk_dzErr_branch = tree->GetBranch("trk_dzErr");
		trk_dzErr_branch->SetAddress(&trk_dzErr_);
	}
	if(trk_dzErr_branch == 0 ) {
	cout << "Branch trk_dzErr does not exist." << endl;
	}
	trk_nChi2_branch = 0;
	if (tree->GetBranch("trk_nChi2") != 0) {
		trk_nChi2_branch = tree->GetBranch("trk_nChi2");
		trk_nChi2_branch->SetAddress(&trk_nChi2_);
	}
	if(trk_nChi2_branch == 0 ) {
	cout << "Branch trk_nChi2 does not exist." << endl;
	}
	trk_shareFrac_branch = 0;
	if (tree->GetBranch("trk_shareFrac") != 0) {
		trk_shareFrac_branch = tree->GetBranch("trk_shareFrac");
		trk_shareFrac_branch->SetAddress(&trk_shareFrac_);
	}
	if(trk_shareFrac_branch == 0 ) {
	cout << "Branch trk_shareFrac does not exist." << endl;
	}
	trk_q_branch = 0;
	if (tree->GetBranch("trk_q") != 0) {
		trk_q_branch = tree->GetBranch("trk_q");
		trk_q_branch->SetAddress(&trk_q_);
	}
	if(trk_q_branch == 0 ) {
	cout << "Branch trk_q does not exist." << endl;
	}
	trk_nValid_branch = 0;
	if (tree->GetBranch("trk_nValid") != 0) {
		trk_nValid_branch = tree->GetBranch("trk_nValid");
		trk_nValid_branch->SetAddress(&trk_nValid_);
	}
	if(trk_nValid_branch == 0 ) {
	cout << "Branch trk_nValid does not exist." << endl;
	}
	trk_nInvalid_branch = 0;
	if (tree->GetBranch("trk_nInvalid") != 0) {
		trk_nInvalid_branch = tree->GetBranch("trk_nInvalid");
		trk_nInvalid_branch->SetAddress(&trk_nInvalid_);
	}
	if(trk_nInvalid_branch == 0 ) {
	cout << "Branch trk_nInvalid does not exist." << endl;
	}
	trk_nPixel_branch = 0;
	if (tree->GetBranch("trk_nPixel") != 0) {
		trk_nPixel_branch = tree->GetBranch("trk_nPixel");
		trk_nPixel_branch->SetAddress(&trk_nPixel_);
	}
	if(trk_nPixel_branch == 0 ) {
	cout << "Branch trk_nPixel does not exist." << endl;
	}
	trk_nStrip_branch = 0;
	if (tree->GetBranch("trk_nStrip") != 0) {
		trk_nStrip_branch = tree->GetBranch("trk_nStrip");
		trk_nStrip_branch->SetAddress(&trk_nStrip_);
	}
	if(trk_nStrip_branch == 0 ) {
	cout << "Branch trk_nStrip does not exist." << endl;
	}
	trk_n3DLay_branch = 0;
	if (tree->GetBranch("trk_n3DLay") != 0) {
		trk_n3DLay_branch = tree->GetBranch("trk_n3DLay");
		trk_n3DLay_branch->SetAddress(&trk_n3DLay_);
	}
	if(trk_n3DLay_branch == 0 ) {
	cout << "Branch trk_n3DLay does not exist." << endl;
	}
	trk_algo_branch = 0;
	if (tree->GetBranch("trk_algo") != 0) {
		trk_algo_branch = tree->GetBranch("trk_algo");
		trk_algo_branch->SetAddress(&trk_algo_);
	}
	if(trk_algo_branch == 0 ) {
	cout << "Branch trk_algo does not exist." << endl;
	}
	trk_isHP_branch = 0;
	if (tree->GetBranch("trk_isHP") != 0) {
		trk_isHP_branch = tree->GetBranch("trk_isHP");
		trk_isHP_branch->SetAddress(&trk_isHP_);
	}
	if(trk_isHP_branch == 0 ) {
	cout << "Branch trk_isHP does not exist." << endl;
	}
	trk_seedIdx_branch = 0;
	if (tree->GetBranch("trk_seedIdx") != 0) {
		trk_seedIdx_branch = tree->GetBranch("trk_seedIdx");
		trk_seedIdx_branch->SetAddress(&trk_seedIdx_);
	}
	if(trk_seedIdx_branch == 0 ) {
	cout << "Branch trk_seedIdx does not exist." << endl;
	}
	trk_simIdx_branch = 0;
	if (tree->GetBranch("trk_simIdx") != 0) {
		trk_simIdx_branch = tree->GetBranch("trk_simIdx");
		trk_simIdx_branch->SetAddress(&trk_simIdx_);
	}
	if(trk_simIdx_branch == 0 ) {
	cout << "Branch trk_simIdx does not exist." << endl;
	}
	trk_pixelIdx_branch = 0;
	if (tree->GetBranch("trk_pixelIdx") != 0) {
		trk_pixelIdx_branch = tree->GetBranch("trk_pixelIdx");
		trk_pixelIdx_branch->SetAddress(&trk_pixelIdx_);
	}
	if(trk_pixelIdx_branch == 0 ) {
	cout << "Branch trk_pixelIdx does not exist." << endl;
	}
	trk_stripIdx_branch = 0;
	if (tree->GetBranch("trk_stripIdx") != 0) {
		trk_stripIdx_branch = tree->GetBranch("trk_stripIdx");
		trk_stripIdx_branch->SetAddress(&trk_stripIdx_);
	}
	if(trk_stripIdx_branch == 0 ) {
	cout << "Branch trk_stripIdx does not exist." << endl;
	}
	sim_px_branch = 0;
	if (tree->GetBranch("sim_px") != 0) {
		sim_px_branch = tree->GetBranch("sim_px");
		sim_px_branch->SetAddress(&sim_px_);
	}
	if(sim_px_branch == 0 ) {
	cout << "Branch sim_px does not exist." << endl;
	}
	sim_py_branch = 0;
	if (tree->GetBranch("sim_py") != 0) {
		sim_py_branch = tree->GetBranch("sim_py");
		sim_py_branch->SetAddress(&sim_py_);
	}
	if(sim_py_branch == 0 ) {
	cout << "Branch sim_py does not exist." << endl;
	}
	sim_pz_branch = 0;
	if (tree->GetBranch("sim_pz") != 0) {
		sim_pz_branch = tree->GetBranch("sim_pz");
		sim_pz_branch->SetAddress(&sim_pz_);
	}
	if(sim_pz_branch == 0 ) {
	cout << "Branch sim_pz does not exist." << endl;
	}
	sim_pt_branch = 0;
	if (tree->GetBranch("sim_pt") != 0) {
		sim_pt_branch = tree->GetBranch("sim_pt");
		sim_pt_branch->SetAddress(&sim_pt_);
	}
	if(sim_pt_branch == 0 ) {
	cout << "Branch sim_pt does not exist." << endl;
	}
	sim_eta_branch = 0;
	if (tree->GetBranch("sim_eta") != 0) {
		sim_eta_branch = tree->GetBranch("sim_eta");
		sim_eta_branch->SetAddress(&sim_eta_);
	}
	if(sim_eta_branch == 0 ) {
	cout << "Branch sim_eta does not exist." << endl;
	}
	sim_phi_branch = 0;
	if (tree->GetBranch("sim_phi") != 0) {
		sim_phi_branch = tree->GetBranch("sim_phi");
		sim_phi_branch->SetAddress(&sim_phi_);
	}
	if(sim_phi_branch == 0 ) {
	cout << "Branch sim_phi does not exist." << endl;
	}
	sim_dxy_branch = 0;
	if (tree->GetBranch("sim_dxy") != 0) {
		sim_dxy_branch = tree->GetBranch("sim_dxy");
		sim_dxy_branch->SetAddress(&sim_dxy_);
	}
	if(sim_dxy_branch == 0 ) {
	cout << "Branch sim_dxy does not exist." << endl;
	}
	sim_dz_branch = 0;
	if (tree->GetBranch("sim_dz") != 0) {
		sim_dz_branch = tree->GetBranch("sim_dz");
		sim_dz_branch->SetAddress(&sim_dz_);
	}
	if(sim_dz_branch == 0 ) {
	cout << "Branch sim_dz does not exist." << endl;
	}
	sim_prodx_branch = 0;
	if (tree->GetBranch("sim_prodx") != 0) {
		sim_prodx_branch = tree->GetBranch("sim_prodx");
		sim_prodx_branch->SetAddress(&sim_prodx_);
	}
	if(sim_prodx_branch == 0 ) {
	cout << "Branch sim_prodx does not exist." << endl;
	}
	sim_prody_branch = 0;
	if (tree->GetBranch("sim_prody") != 0) {
		sim_prody_branch = tree->GetBranch("sim_prody");
		sim_prody_branch->SetAddress(&sim_prody_);
	}
	if(sim_prody_branch == 0 ) {
	cout << "Branch sim_prody does not exist." << endl;
	}
	sim_prodz_branch = 0;
	if (tree->GetBranch("sim_prodz") != 0) {
		sim_prodz_branch = tree->GetBranch("sim_prodz");
		sim_prodz_branch->SetAddress(&sim_prodz_);
	}
	if(sim_prodz_branch == 0 ) {
	cout << "Branch sim_prodz does not exist." << endl;
	}
	sim_shareFrac_branch = 0;
	if (tree->GetBranch("sim_shareFrac") != 0) {
		sim_shareFrac_branch = tree->GetBranch("sim_shareFrac");
		sim_shareFrac_branch->SetAddress(&sim_shareFrac_);
	}
	if(sim_shareFrac_branch == 0 ) {
	cout << "Branch sim_shareFrac does not exist." << endl;
	}
	sim_q_branch = 0;
	if (tree->GetBranch("sim_q") != 0) {
		sim_q_branch = tree->GetBranch("sim_q");
		sim_q_branch->SetAddress(&sim_q_);
	}
	if(sim_q_branch == 0 ) {
	cout << "Branch sim_q does not exist." << endl;
	}
	sim_nValid_branch = 0;
	if (tree->GetBranch("sim_nValid") != 0) {
		sim_nValid_branch = tree->GetBranch("sim_nValid");
		sim_nValid_branch->SetAddress(&sim_nValid_);
	}
	if(sim_nValid_branch == 0 ) {
	cout << "Branch sim_nValid does not exist." << endl;
	}
	sim_nPixel_branch = 0;
	if (tree->GetBranch("sim_nPixel") != 0) {
		sim_nPixel_branch = tree->GetBranch("sim_nPixel");
		sim_nPixel_branch->SetAddress(&sim_nPixel_);
	}
	if(sim_nPixel_branch == 0 ) {
	cout << "Branch sim_nPixel does not exist." << endl;
	}
	sim_nStrip_branch = 0;
	if (tree->GetBranch("sim_nStrip") != 0) {
		sim_nStrip_branch = tree->GetBranch("sim_nStrip");
		sim_nStrip_branch->SetAddress(&sim_nStrip_);
	}
	if(sim_nStrip_branch == 0 ) {
	cout << "Branch sim_nStrip does not exist." << endl;
	}
	sim_n3DLay_branch = 0;
	if (tree->GetBranch("sim_n3DLay") != 0) {
		sim_n3DLay_branch = tree->GetBranch("sim_n3DLay");
		sim_n3DLay_branch->SetAddress(&sim_n3DLay_);
	}
	if(sim_n3DLay_branch == 0 ) {
	cout << "Branch sim_n3DLay does not exist." << endl;
	}
	sim_trkIdx_branch = 0;
	if (tree->GetBranch("sim_trkIdx") != 0) {
		sim_trkIdx_branch = tree->GetBranch("sim_trkIdx");
		sim_trkIdx_branch->SetAddress(&sim_trkIdx_);
	}
	if(sim_trkIdx_branch == 0 ) {
	cout << "Branch sim_trkIdx does not exist." << endl;
	}
	sim_pixelIdx_branch = 0;
	if (tree->GetBranch("sim_pixelIdx") != 0) {
		sim_pixelIdx_branch = tree->GetBranch("sim_pixelIdx");
		sim_pixelIdx_branch->SetAddress(&sim_pixelIdx_);
	}
	if(sim_pixelIdx_branch == 0 ) {
	cout << "Branch sim_pixelIdx does not exist." << endl;
	}
	sim_stripIdx_branch = 0;
	if (tree->GetBranch("sim_stripIdx") != 0) {
		sim_stripIdx_branch = tree->GetBranch("sim_stripIdx");
		sim_stripIdx_branch->SetAddress(&sim_stripIdx_);
	}
	if(sim_stripIdx_branch == 0 ) {
	cout << "Branch sim_stripIdx does not exist." << endl;
	}
	pix_isBarrel_branch = 0;
	if (tree->GetBranch("pix_isBarrel") != 0) {
		pix_isBarrel_branch = tree->GetBranch("pix_isBarrel");
		pix_isBarrel_branch->SetAddress(&pix_isBarrel_);
	}
	if(pix_isBarrel_branch == 0 ) {
	cout << "Branch pix_isBarrel does not exist." << endl;
	}
	pix_lay_branch = 0;
	if (tree->GetBranch("pix_lay") != 0) {
		pix_lay_branch = tree->GetBranch("pix_lay");
		pix_lay_branch->SetAddress(&pix_lay_);
	}
	if(pix_lay_branch == 0 ) {
	cout << "Branch pix_lay does not exist." << endl;
	}
	pix_detId_branch = 0;
	if (tree->GetBranch("pix_detId") != 0) {
		pix_detId_branch = tree->GetBranch("pix_detId");
		pix_detId_branch->SetAddress(&pix_detId_);
	}
	if(pix_detId_branch == 0 ) {
	cout << "Branch pix_detId does not exist." << endl;
	}
	pix_nSimTrk_branch = 0;
	if (tree->GetBranch("pix_nSimTrk") != 0) {
		pix_nSimTrk_branch = tree->GetBranch("pix_nSimTrk");
		pix_nSimTrk_branch->SetAddress(&pix_nSimTrk_);
	}
	if(pix_nSimTrk_branch == 0 ) {
	cout << "Branch pix_nSimTrk does not exist." << endl;
	}
	pix_simTrkIdx_branch = 0;
	if (tree->GetBranch("pix_simTrkIdx") != 0) {
		pix_simTrkIdx_branch = tree->GetBranch("pix_simTrkIdx");
		pix_simTrkIdx_branch->SetAddress(&pix_simTrkIdx_);
	}
	if(pix_simTrkIdx_branch == 0 ) {
	cout << "Branch pix_simTrkIdx does not exist." << endl;
	}
	pix_particle_branch = 0;
	if (tree->GetBranch("pix_particle") != 0) {
		pix_particle_branch = tree->GetBranch("pix_particle");
		pix_particle_branch->SetAddress(&pix_particle_);
	}
	if(pix_particle_branch == 0 ) {
	cout << "Branch pix_particle does not exist." << endl;
	}
	pix_process_branch = 0;
	if (tree->GetBranch("pix_process") != 0) {
		pix_process_branch = tree->GetBranch("pix_process");
		pix_process_branch->SetAddress(&pix_process_);
	}
	if(pix_process_branch == 0 ) {
	cout << "Branch pix_process does not exist." << endl;
	}
	pix_bunchXing_branch = 0;
	if (tree->GetBranch("pix_bunchXing") != 0) {
		pix_bunchXing_branch = tree->GetBranch("pix_bunchXing");
		pix_bunchXing_branch->SetAddress(&pix_bunchXing_);
	}
	if(pix_bunchXing_branch == 0 ) {
	cout << "Branch pix_bunchXing does not exist." << endl;
	}
	pix_event_branch = 0;
	if (tree->GetBranch("pix_event") != 0) {
		pix_event_branch = tree->GetBranch("pix_event");
		pix_event_branch->SetAddress(&pix_event_);
	}
	if(pix_event_branch == 0 ) {
	cout << "Branch pix_event does not exist." << endl;
	}
	pix_x_branch = 0;
	if (tree->GetBranch("pix_x") != 0) {
		pix_x_branch = tree->GetBranch("pix_x");
		pix_x_branch->SetAddress(&pix_x_);
	}
	if(pix_x_branch == 0 ) {
	cout << "Branch pix_x does not exist." << endl;
	}
	pix_y_branch = 0;
	if (tree->GetBranch("pix_y") != 0) {
		pix_y_branch = tree->GetBranch("pix_y");
		pix_y_branch->SetAddress(&pix_y_);
	}
	if(pix_y_branch == 0 ) {
	cout << "Branch pix_y does not exist." << endl;
	}
	pix_z_branch = 0;
	if (tree->GetBranch("pix_z") != 0) {
		pix_z_branch = tree->GetBranch("pix_z");
		pix_z_branch->SetAddress(&pix_z_);
	}
	if(pix_z_branch == 0 ) {
	cout << "Branch pix_z does not exist." << endl;
	}
	pix_xx_branch = 0;
	if (tree->GetBranch("pix_xx") != 0) {
		pix_xx_branch = tree->GetBranch("pix_xx");
		pix_xx_branch->SetAddress(&pix_xx_);
	}
	if(pix_xx_branch == 0 ) {
	cout << "Branch pix_xx does not exist." << endl;
	}
	pix_xy_branch = 0;
	if (tree->GetBranch("pix_xy") != 0) {
		pix_xy_branch = tree->GetBranch("pix_xy");
		pix_xy_branch->SetAddress(&pix_xy_);
	}
	if(pix_xy_branch == 0 ) {
	cout << "Branch pix_xy does not exist." << endl;
	}
	pix_yy_branch = 0;
	if (tree->GetBranch("pix_yy") != 0) {
		pix_yy_branch = tree->GetBranch("pix_yy");
		pix_yy_branch->SetAddress(&pix_yy_);
	}
	if(pix_yy_branch == 0 ) {
	cout << "Branch pix_yy does not exist." << endl;
	}
	pix_yz_branch = 0;
	if (tree->GetBranch("pix_yz") != 0) {
		pix_yz_branch = tree->GetBranch("pix_yz");
		pix_yz_branch->SetAddress(&pix_yz_);
	}
	if(pix_yz_branch == 0 ) {
	cout << "Branch pix_yz does not exist." << endl;
	}
	pix_zz_branch = 0;
	if (tree->GetBranch("pix_zz") != 0) {
		pix_zz_branch = tree->GetBranch("pix_zz");
		pix_zz_branch->SetAddress(&pix_zz_);
	}
	if(pix_zz_branch == 0 ) {
	cout << "Branch pix_zz does not exist." << endl;
	}
	pix_zx_branch = 0;
	if (tree->GetBranch("pix_zx") != 0) {
		pix_zx_branch = tree->GetBranch("pix_zx");
		pix_zx_branch->SetAddress(&pix_zx_);
	}
	if(pix_zx_branch == 0 ) {
	cout << "Branch pix_zx does not exist." << endl;
	}
	pix_xsim_branch = 0;
	if (tree->GetBranch("pix_xsim") != 0) {
		pix_xsim_branch = tree->GetBranch("pix_xsim");
		pix_xsim_branch->SetAddress(&pix_xsim_);
	}
	if(pix_xsim_branch == 0 ) {
	cout << "Branch pix_xsim does not exist." << endl;
	}
	pix_ysim_branch = 0;
	if (tree->GetBranch("pix_ysim") != 0) {
		pix_ysim_branch = tree->GetBranch("pix_ysim");
		pix_ysim_branch->SetAddress(&pix_ysim_);
	}
	if(pix_ysim_branch == 0 ) {
	cout << "Branch pix_ysim does not exist." << endl;
	}
	pix_zsim_branch = 0;
	if (tree->GetBranch("pix_zsim") != 0) {
		pix_zsim_branch = tree->GetBranch("pix_zsim");
		pix_zsim_branch->SetAddress(&pix_zsim_);
	}
	if(pix_zsim_branch == 0 ) {
	cout << "Branch pix_zsim does not exist." << endl;
	}
	pix_pxsim_branch = 0;
	if (tree->GetBranch("pix_pxsim") != 0) {
		pix_pxsim_branch = tree->GetBranch("pix_pxsim");
		pix_pxsim_branch->SetAddress(&pix_pxsim_);
	}
	if(pix_pxsim_branch == 0 ) {
	cout << "Branch pix_pxsim does not exist." << endl;
	}
	pix_pysim_branch = 0;
	if (tree->GetBranch("pix_pysim") != 0) {
		pix_pysim_branch = tree->GetBranch("pix_pysim");
		pix_pysim_branch->SetAddress(&pix_pysim_);
	}
	if(pix_pysim_branch == 0 ) {
	cout << "Branch pix_pysim does not exist." << endl;
	}
	pix_pzsim_branch = 0;
	if (tree->GetBranch("pix_pzsim") != 0) {
		pix_pzsim_branch = tree->GetBranch("pix_pzsim");
		pix_pzsim_branch->SetAddress(&pix_pzsim_);
	}
	if(pix_pzsim_branch == 0 ) {
	cout << "Branch pix_pzsim does not exist." << endl;
	}
	pix_pathprop_branch = 0;
	if (tree->GetBranch("pix_pathprop") != 0) {
		pix_pathprop_branch = tree->GetBranch("pix_pathprop");
		pix_pathprop_branch->SetAddress(&pix_pathprop_);
	}
	if(pix_pathprop_branch == 0 ) {
	cout << "Branch pix_pathprop does not exist." << endl;
	}
	pix_xsimprop_branch = 0;
	if (tree->GetBranch("pix_xsimprop") != 0) {
		pix_xsimprop_branch = tree->GetBranch("pix_xsimprop");
		pix_xsimprop_branch->SetAddress(&pix_xsimprop_);
	}
	if(pix_xsimprop_branch == 0 ) {
	cout << "Branch pix_xsimprop does not exist." << endl;
	}
	pix_ysimprop_branch = 0;
	if (tree->GetBranch("pix_ysimprop") != 0) {
		pix_ysimprop_branch = tree->GetBranch("pix_ysimprop");
		pix_ysimprop_branch->SetAddress(&pix_ysimprop_);
	}
	if(pix_ysimprop_branch == 0 ) {
	cout << "Branch pix_ysimprop does not exist." << endl;
	}
	pix_zsimprop_branch = 0;
	if (tree->GetBranch("pix_zsimprop") != 0) {
		pix_zsimprop_branch = tree->GetBranch("pix_zsimprop");
		pix_zsimprop_branch->SetAddress(&pix_zsimprop_);
	}
	if(pix_zsimprop_branch == 0 ) {
	cout << "Branch pix_zsimprop does not exist." << endl;
	}
	pix_pxsimprop_branch = 0;
	if (tree->GetBranch("pix_pxsimprop") != 0) {
		pix_pxsimprop_branch = tree->GetBranch("pix_pxsimprop");
		pix_pxsimprop_branch->SetAddress(&pix_pxsimprop_);
	}
	if(pix_pxsimprop_branch == 0 ) {
	cout << "Branch pix_pxsimprop does not exist." << endl;
	}
	pix_pysimprop_branch = 0;
	if (tree->GetBranch("pix_pysimprop") != 0) {
		pix_pysimprop_branch = tree->GetBranch("pix_pysimprop");
		pix_pysimprop_branch->SetAddress(&pix_pysimprop_);
	}
	if(pix_pysimprop_branch == 0 ) {
	cout << "Branch pix_pysimprop does not exist." << endl;
	}
	pix_pzsimprop_branch = 0;
	if (tree->GetBranch("pix_pzsimprop") != 0) {
		pix_pzsimprop_branch = tree->GetBranch("pix_pzsimprop");
		pix_pzsimprop_branch->SetAddress(&pix_pzsimprop_);
	}
	if(pix_pzsimprop_branch == 0 ) {
	cout << "Branch pix_pzsimprop does not exist." << endl;
	}
	pix_eloss_branch = 0;
	if (tree->GetBranch("pix_eloss") != 0) {
		pix_eloss_branch = tree->GetBranch("pix_eloss");
		pix_eloss_branch->SetAddress(&pix_eloss_);
	}
	if(pix_eloss_branch == 0 ) {
	cout << "Branch pix_eloss does not exist." << endl;
	}
	pix_radL_branch = 0;
	if (tree->GetBranch("pix_radL") != 0) {
		pix_radL_branch = tree->GetBranch("pix_radL");
		pix_radL_branch->SetAddress(&pix_radL_);
	}
	if(pix_radL_branch == 0 ) {
	cout << "Branch pix_radL does not exist." << endl;
	}
	pix_bbxi_branch = 0;
	if (tree->GetBranch("pix_bbxi") != 0) {
		pix_bbxi_branch = tree->GetBranch("pix_bbxi");
		pix_bbxi_branch->SetAddress(&pix_bbxi_);
	}
	if(pix_bbxi_branch == 0 ) {
	cout << "Branch pix_bbxi does not exist." << endl;
	}
	str_isBarrel_branch = 0;
	if (tree->GetBranch("str_isBarrel") != 0) {
		str_isBarrel_branch = tree->GetBranch("str_isBarrel");
		str_isBarrel_branch->SetAddress(&str_isBarrel_);
	}
	if(str_isBarrel_branch == 0 ) {
	cout << "Branch str_isBarrel does not exist." << endl;
	}
	str_isStereo_branch = 0;
	if (tree->GetBranch("str_isStereo") != 0) {
		str_isStereo_branch = tree->GetBranch("str_isStereo");
		str_isStereo_branch->SetAddress(&str_isStereo_);
	}
	if(str_isStereo_branch == 0 ) {
	cout << "Branch str_isStereo does not exist." << endl;
	}
	str_det_branch = 0;
	if (tree->GetBranch("str_det") != 0) {
		str_det_branch = tree->GetBranch("str_det");
		str_det_branch->SetAddress(&str_det_);
	}
	if(str_det_branch == 0 ) {
	cout << "Branch str_det does not exist." << endl;
	}
	str_lay_branch = 0;
	if (tree->GetBranch("str_lay") != 0) {
		str_lay_branch = tree->GetBranch("str_lay");
		str_lay_branch->SetAddress(&str_lay_);
	}
	if(str_lay_branch == 0 ) {
	cout << "Branch str_lay does not exist." << endl;
	}
	str_detId_branch = 0;
	if (tree->GetBranch("str_detId") != 0) {
		str_detId_branch = tree->GetBranch("str_detId");
		str_detId_branch->SetAddress(&str_detId_);
	}
	if(str_detId_branch == 0 ) {
	cout << "Branch str_detId does not exist." << endl;
	}
	str_nSimTrk_branch = 0;
	if (tree->GetBranch("str_nSimTrk") != 0) {
		str_nSimTrk_branch = tree->GetBranch("str_nSimTrk");
		str_nSimTrk_branch->SetAddress(&str_nSimTrk_);
	}
	if(str_nSimTrk_branch == 0 ) {
	cout << "Branch str_nSimTrk does not exist." << endl;
	}
	str_simTrkIdx_branch = 0;
	if (tree->GetBranch("str_simTrkIdx") != 0) {
		str_simTrkIdx_branch = tree->GetBranch("str_simTrkIdx");
		str_simTrkIdx_branch->SetAddress(&str_simTrkIdx_);
	}
	if(str_simTrkIdx_branch == 0 ) {
	cout << "Branch str_simTrkIdx does not exist." << endl;
	}
	str_particle_branch = 0;
	if (tree->GetBranch("str_particle") != 0) {
		str_particle_branch = tree->GetBranch("str_particle");
		str_particle_branch->SetAddress(&str_particle_);
	}
	if(str_particle_branch == 0 ) {
	cout << "Branch str_particle does not exist." << endl;
	}
	str_process_branch = 0;
	if (tree->GetBranch("str_process") != 0) {
		str_process_branch = tree->GetBranch("str_process");
		str_process_branch->SetAddress(&str_process_);
	}
	if(str_process_branch == 0 ) {
	cout << "Branch str_process does not exist." << endl;
	}
	str_bunchXing_branch = 0;
	if (tree->GetBranch("str_bunchXing") != 0) {
		str_bunchXing_branch = tree->GetBranch("str_bunchXing");
		str_bunchXing_branch->SetAddress(&str_bunchXing_);
	}
	if(str_bunchXing_branch == 0 ) {
	cout << "Branch str_bunchXing does not exist." << endl;
	}
	str_event_branch = 0;
	if (tree->GetBranch("str_event") != 0) {
		str_event_branch = tree->GetBranch("str_event");
		str_event_branch->SetAddress(&str_event_);
	}
	if(str_event_branch == 0 ) {
	cout << "Branch str_event does not exist." << endl;
	}
	str_x_branch = 0;
	if (tree->GetBranch("str_x") != 0) {
		str_x_branch = tree->GetBranch("str_x");
		str_x_branch->SetAddress(&str_x_);
	}
	if(str_x_branch == 0 ) {
	cout << "Branch str_x does not exist." << endl;
	}
	str_y_branch = 0;
	if (tree->GetBranch("str_y") != 0) {
		str_y_branch = tree->GetBranch("str_y");
		str_y_branch->SetAddress(&str_y_);
	}
	if(str_y_branch == 0 ) {
	cout << "Branch str_y does not exist." << endl;
	}
	str_z_branch = 0;
	if (tree->GetBranch("str_z") != 0) {
		str_z_branch = tree->GetBranch("str_z");
		str_z_branch->SetAddress(&str_z_);
	}
	if(str_z_branch == 0 ) {
	cout << "Branch str_z does not exist." << endl;
	}
	str_xx_branch = 0;
	if (tree->GetBranch("str_xx") != 0) {
		str_xx_branch = tree->GetBranch("str_xx");
		str_xx_branch->SetAddress(&str_xx_);
	}
	if(str_xx_branch == 0 ) {
	cout << "Branch str_xx does not exist." << endl;
	}
	str_xy_branch = 0;
	if (tree->GetBranch("str_xy") != 0) {
		str_xy_branch = tree->GetBranch("str_xy");
		str_xy_branch->SetAddress(&str_xy_);
	}
	if(str_xy_branch == 0 ) {
	cout << "Branch str_xy does not exist." << endl;
	}
	str_yy_branch = 0;
	if (tree->GetBranch("str_yy") != 0) {
		str_yy_branch = tree->GetBranch("str_yy");
		str_yy_branch->SetAddress(&str_yy_);
	}
	if(str_yy_branch == 0 ) {
	cout << "Branch str_yy does not exist." << endl;
	}
	str_yz_branch = 0;
	if (tree->GetBranch("str_yz") != 0) {
		str_yz_branch = tree->GetBranch("str_yz");
		str_yz_branch->SetAddress(&str_yz_);
	}
	if(str_yz_branch == 0 ) {
	cout << "Branch str_yz does not exist." << endl;
	}
	str_zz_branch = 0;
	if (tree->GetBranch("str_zz") != 0) {
		str_zz_branch = tree->GetBranch("str_zz");
		str_zz_branch->SetAddress(&str_zz_);
	}
	if(str_zz_branch == 0 ) {
	cout << "Branch str_zz does not exist." << endl;
	}
	str_zx_branch = 0;
	if (tree->GetBranch("str_zx") != 0) {
		str_zx_branch = tree->GetBranch("str_zx");
		str_zx_branch->SetAddress(&str_zx_);
	}
	if(str_zx_branch == 0 ) {
	cout << "Branch str_zx does not exist." << endl;
	}
	str_xsim_branch = 0;
	if (tree->GetBranch("str_xsim") != 0) {
		str_xsim_branch = tree->GetBranch("str_xsim");
		str_xsim_branch->SetAddress(&str_xsim_);
	}
	if(str_xsim_branch == 0 ) {
	cout << "Branch str_xsim does not exist." << endl;
	}
	str_ysim_branch = 0;
	if (tree->GetBranch("str_ysim") != 0) {
		str_ysim_branch = tree->GetBranch("str_ysim");
		str_ysim_branch->SetAddress(&str_ysim_);
	}
	if(str_ysim_branch == 0 ) {
	cout << "Branch str_ysim does not exist." << endl;
	}
	str_zsim_branch = 0;
	if (tree->GetBranch("str_zsim") != 0) {
		str_zsim_branch = tree->GetBranch("str_zsim");
		str_zsim_branch->SetAddress(&str_zsim_);
	}
	if(str_zsim_branch == 0 ) {
	cout << "Branch str_zsim does not exist." << endl;
	}
	str_pxsim_branch = 0;
	if (tree->GetBranch("str_pxsim") != 0) {
		str_pxsim_branch = tree->GetBranch("str_pxsim");
		str_pxsim_branch->SetAddress(&str_pxsim_);
	}
	if(str_pxsim_branch == 0 ) {
	cout << "Branch str_pxsim does not exist." << endl;
	}
	str_pysim_branch = 0;
	if (tree->GetBranch("str_pysim") != 0) {
		str_pysim_branch = tree->GetBranch("str_pysim");
		str_pysim_branch->SetAddress(&str_pysim_);
	}
	if(str_pysim_branch == 0 ) {
	cout << "Branch str_pysim does not exist." << endl;
	}
	str_pzsim_branch = 0;
	if (tree->GetBranch("str_pzsim") != 0) {
		str_pzsim_branch = tree->GetBranch("str_pzsim");
		str_pzsim_branch->SetAddress(&str_pzsim_);
	}
	if(str_pzsim_branch == 0 ) {
	cout << "Branch str_pzsim does not exist." << endl;
	}
	str_eloss_branch = 0;
	if (tree->GetBranch("str_eloss") != 0) {
		str_eloss_branch = tree->GetBranch("str_eloss");
		str_eloss_branch->SetAddress(&str_eloss_);
	}
	if(str_eloss_branch == 0 ) {
	cout << "Branch str_eloss does not exist." << endl;
	}
	str_radL_branch = 0;
	if (tree->GetBranch("str_radL") != 0) {
		str_radL_branch = tree->GetBranch("str_radL");
		str_radL_branch->SetAddress(&str_radL_);
	}
	if(str_radL_branch == 0 ) {
	cout << "Branch str_radL does not exist." << endl;
	}
	str_bbxi_branch = 0;
	if (tree->GetBranch("str_bbxi") != 0) {
		str_bbxi_branch = tree->GetBranch("str_bbxi");
		str_bbxi_branch->SetAddress(&str_bbxi_);
	}
	if(str_bbxi_branch == 0 ) {
	cout << "Branch str_bbxi does not exist." << endl;
	}
	glu_isBarrel_branch = 0;
	if (tree->GetBranch("glu_isBarrel") != 0) {
		glu_isBarrel_branch = tree->GetBranch("glu_isBarrel");
		glu_isBarrel_branch->SetAddress(&glu_isBarrel_);
	}
	if(glu_isBarrel_branch == 0 ) {
	cout << "Branch glu_isBarrel does not exist." << endl;
	}
	glu_det_branch = 0;
	if (tree->GetBranch("glu_det") != 0) {
		glu_det_branch = tree->GetBranch("glu_det");
		glu_det_branch->SetAddress(&glu_det_);
	}
	if(glu_det_branch == 0 ) {
	cout << "Branch glu_det does not exist." << endl;
	}
	glu_lay_branch = 0;
	if (tree->GetBranch("glu_lay") != 0) {
		glu_lay_branch = tree->GetBranch("glu_lay");
		glu_lay_branch->SetAddress(&glu_lay_);
	}
	if(glu_lay_branch == 0 ) {
	cout << "Branch glu_lay does not exist." << endl;
	}
	glu_detId_branch = 0;
	if (tree->GetBranch("glu_detId") != 0) {
		glu_detId_branch = tree->GetBranch("glu_detId");
		glu_detId_branch->SetAddress(&glu_detId_);
	}
	if(glu_detId_branch == 0 ) {
	cout << "Branch glu_detId does not exist." << endl;
	}
	glu_monoIdx_branch = 0;
	if (tree->GetBranch("glu_monoIdx") != 0) {
		glu_monoIdx_branch = tree->GetBranch("glu_monoIdx");
		glu_monoIdx_branch->SetAddress(&glu_monoIdx_);
	}
	if(glu_monoIdx_branch == 0 ) {
	cout << "Branch glu_monoIdx does not exist." << endl;
	}
	glu_stereoIdx_branch = 0;
	if (tree->GetBranch("glu_stereoIdx") != 0) {
		glu_stereoIdx_branch = tree->GetBranch("glu_stereoIdx");
		glu_stereoIdx_branch->SetAddress(&glu_stereoIdx_);
	}
	if(glu_stereoIdx_branch == 0 ) {
	cout << "Branch glu_stereoIdx does not exist." << endl;
	}
	glu_x_branch = 0;
	if (tree->GetBranch("glu_x") != 0) {
		glu_x_branch = tree->GetBranch("glu_x");
		glu_x_branch->SetAddress(&glu_x_);
	}
	if(glu_x_branch == 0 ) {
	cout << "Branch glu_x does not exist." << endl;
	}
	glu_y_branch = 0;
	if (tree->GetBranch("glu_y") != 0) {
		glu_y_branch = tree->GetBranch("glu_y");
		glu_y_branch->SetAddress(&glu_y_);
	}
	if(glu_y_branch == 0 ) {
	cout << "Branch glu_y does not exist." << endl;
	}
	glu_z_branch = 0;
	if (tree->GetBranch("glu_z") != 0) {
		glu_z_branch = tree->GetBranch("glu_z");
		glu_z_branch->SetAddress(&glu_z_);
	}
	if(glu_z_branch == 0 ) {
	cout << "Branch glu_z does not exist." << endl;
	}
	glu_xx_branch = 0;
	if (tree->GetBranch("glu_xx") != 0) {
		glu_xx_branch = tree->GetBranch("glu_xx");
		glu_xx_branch->SetAddress(&glu_xx_);
	}
	if(glu_xx_branch == 0 ) {
	cout << "Branch glu_xx does not exist." << endl;
	}
	glu_xy_branch = 0;
	if (tree->GetBranch("glu_xy") != 0) {
		glu_xy_branch = tree->GetBranch("glu_xy");
		glu_xy_branch->SetAddress(&glu_xy_);
	}
	if(glu_xy_branch == 0 ) {
	cout << "Branch glu_xy does not exist." << endl;
	}
	glu_yy_branch = 0;
	if (tree->GetBranch("glu_yy") != 0) {
		glu_yy_branch = tree->GetBranch("glu_yy");
		glu_yy_branch->SetAddress(&glu_yy_);
	}
	if(glu_yy_branch == 0 ) {
	cout << "Branch glu_yy does not exist." << endl;
	}
	glu_yz_branch = 0;
	if (tree->GetBranch("glu_yz") != 0) {
		glu_yz_branch = tree->GetBranch("glu_yz");
		glu_yz_branch->SetAddress(&glu_yz_);
	}
	if(glu_yz_branch == 0 ) {
	cout << "Branch glu_yz does not exist." << endl;
	}
	glu_zz_branch = 0;
	if (tree->GetBranch("glu_zz") != 0) {
		glu_zz_branch = tree->GetBranch("glu_zz");
		glu_zz_branch->SetAddress(&glu_zz_);
	}
	if(glu_zz_branch == 0 ) {
	cout << "Branch glu_zz does not exist." << endl;
	}
	glu_zx_branch = 0;
	if (tree->GetBranch("glu_zx") != 0) {
		glu_zx_branch = tree->GetBranch("glu_zx");
		glu_zx_branch->SetAddress(&glu_zx_);
	}
	if(glu_zx_branch == 0 ) {
	cout << "Branch glu_zx does not exist." << endl;
	}
	glu_radL_branch = 0;
	if (tree->GetBranch("glu_radL") != 0) {
		glu_radL_branch = tree->GetBranch("glu_radL");
		glu_radL_branch->SetAddress(&glu_radL_);
	}
	if(glu_radL_branch == 0 ) {
	cout << "Branch glu_radL does not exist." << endl;
	}
	glu_bbxi_branch = 0;
	if (tree->GetBranch("glu_bbxi") != 0) {
		glu_bbxi_branch = tree->GetBranch("glu_bbxi");
		glu_bbxi_branch->SetAddress(&glu_bbxi_);
	}
	if(glu_bbxi_branch == 0 ) {
	cout << "Branch glu_bbxi does not exist." << endl;
	}
	bsp_x_branch = 0;
	if (tree->GetBranch("bsp_x") != 0) {
		bsp_x_branch = tree->GetBranch("bsp_x");
		bsp_x_branch->SetAddress(&bsp_x_);
	}
	if(bsp_x_branch == 0 ) {
	cout << "Branch bsp_x does not exist." << endl;
	}
	bsp_y_branch = 0;
	if (tree->GetBranch("bsp_y") != 0) {
		bsp_y_branch = tree->GetBranch("bsp_y");
		bsp_y_branch->SetAddress(&bsp_y_);
	}
	if(bsp_y_branch == 0 ) {
	cout << "Branch bsp_y does not exist." << endl;
	}
	bsp_z_branch = 0;
	if (tree->GetBranch("bsp_z") != 0) {
		bsp_z_branch = tree->GetBranch("bsp_z");
		bsp_z_branch->SetAddress(&bsp_z_);
	}
	if(bsp_z_branch == 0 ) {
	cout << "Branch bsp_z does not exist." << endl;
	}
	bsp_sigmax_branch = 0;
	if (tree->GetBranch("bsp_sigmax") != 0) {
		bsp_sigmax_branch = tree->GetBranch("bsp_sigmax");
		bsp_sigmax_branch->SetAddress(&bsp_sigmax_);
	}
	if(bsp_sigmax_branch == 0 ) {
	cout << "Branch bsp_sigmax does not exist." << endl;
	}
	bsp_sigmay_branch = 0;
	if (tree->GetBranch("bsp_sigmay") != 0) {
		bsp_sigmay_branch = tree->GetBranch("bsp_sigmay");
		bsp_sigmay_branch->SetAddress(&bsp_sigmay_);
	}
	if(bsp_sigmay_branch == 0 ) {
	cout << "Branch bsp_sigmay does not exist." << endl;
	}
	bsp_sigmaz_branch = 0;
	if (tree->GetBranch("bsp_sigmaz") != 0) {
		bsp_sigmaz_branch = tree->GetBranch("bsp_sigmaz");
		bsp_sigmaz_branch->SetAddress(&bsp_sigmaz_);
	}
	if(bsp_sigmaz_branch == 0 ) {
	cout << "Branch bsp_sigmaz does not exist." << endl;
	}
	see_px_branch = 0;
	if (tree->GetBranch("see_px") != 0) {
		see_px_branch = tree->GetBranch("see_px");
		see_px_branch->SetAddress(&see_px_);
	}
	if(see_px_branch == 0 ) {
	cout << "Branch see_px does not exist." << endl;
	}
	see_py_branch = 0;
	if (tree->GetBranch("see_py") != 0) {
		see_py_branch = tree->GetBranch("see_py");
		see_py_branch->SetAddress(&see_py_);
	}
	if(see_py_branch == 0 ) {
	cout << "Branch see_py does not exist." << endl;
	}
	see_pz_branch = 0;
	if (tree->GetBranch("see_pz") != 0) {
		see_pz_branch = tree->GetBranch("see_pz");
		see_pz_branch->SetAddress(&see_pz_);
	}
	if(see_pz_branch == 0 ) {
	cout << "Branch see_pz does not exist." << endl;
	}
	see_pt_branch = 0;
	if (tree->GetBranch("see_pt") != 0) {
		see_pt_branch = tree->GetBranch("see_pt");
		see_pt_branch->SetAddress(&see_pt_);
	}
	if(see_pt_branch == 0 ) {
	cout << "Branch see_pt does not exist." << endl;
	}
	see_eta_branch = 0;
	if (tree->GetBranch("see_eta") != 0) {
		see_eta_branch = tree->GetBranch("see_eta");
		see_eta_branch->SetAddress(&see_eta_);
	}
	if(see_eta_branch == 0 ) {
	cout << "Branch see_eta does not exist." << endl;
	}
	see_phi_branch = 0;
	if (tree->GetBranch("see_phi") != 0) {
		see_phi_branch = tree->GetBranch("see_phi");
		see_phi_branch->SetAddress(&see_phi_);
	}
	if(see_phi_branch == 0 ) {
	cout << "Branch see_phi does not exist." << endl;
	}
	see_dxy_branch = 0;
	if (tree->GetBranch("see_dxy") != 0) {
		see_dxy_branch = tree->GetBranch("see_dxy");
		see_dxy_branch->SetAddress(&see_dxy_);
	}
	if(see_dxy_branch == 0 ) {
	cout << "Branch see_dxy does not exist." << endl;
	}
	see_dz_branch = 0;
	if (tree->GetBranch("see_dz") != 0) {
		see_dz_branch = tree->GetBranch("see_dz");
		see_dz_branch->SetAddress(&see_dz_);
	}
	if(see_dz_branch == 0 ) {
	cout << "Branch see_dz does not exist." << endl;
	}
	see_ptErr_branch = 0;
	if (tree->GetBranch("see_ptErr") != 0) {
		see_ptErr_branch = tree->GetBranch("see_ptErr");
		see_ptErr_branch->SetAddress(&see_ptErr_);
	}
	if(see_ptErr_branch == 0 ) {
	cout << "Branch see_ptErr does not exist." << endl;
	}
	see_etaErr_branch = 0;
	if (tree->GetBranch("see_etaErr") != 0) {
		see_etaErr_branch = tree->GetBranch("see_etaErr");
		see_etaErr_branch->SetAddress(&see_etaErr_);
	}
	if(see_etaErr_branch == 0 ) {
	cout << "Branch see_etaErr does not exist." << endl;
	}
	see_phiErr_branch = 0;
	if (tree->GetBranch("see_phiErr") != 0) {
		see_phiErr_branch = tree->GetBranch("see_phiErr");
		see_phiErr_branch->SetAddress(&see_phiErr_);
	}
	if(see_phiErr_branch == 0 ) {
	cout << "Branch see_phiErr does not exist." << endl;
	}
	see_dxyErr_branch = 0;
	if (tree->GetBranch("see_dxyErr") != 0) {
		see_dxyErr_branch = tree->GetBranch("see_dxyErr");
		see_dxyErr_branch->SetAddress(&see_dxyErr_);
	}
	if(see_dxyErr_branch == 0 ) {
	cout << "Branch see_dxyErr does not exist." << endl;
	}
	see_dzErr_branch = 0;
	if (tree->GetBranch("see_dzErr") != 0) {
		see_dzErr_branch = tree->GetBranch("see_dzErr");
		see_dzErr_branch->SetAddress(&see_dzErr_);
	}
	if(see_dzErr_branch == 0 ) {
	cout << "Branch see_dzErr does not exist." << endl;
	}
	see_chi2_branch = 0;
	if (tree->GetBranch("see_chi2") != 0) {
		see_chi2_branch = tree->GetBranch("see_chi2");
		see_chi2_branch->SetAddress(&see_chi2_);
	}
	if(see_chi2_branch == 0 ) {
	cout << "Branch see_chi2 does not exist." << endl;
	}
	see_q_branch = 0;
	if (tree->GetBranch("see_q") != 0) {
		see_q_branch = tree->GetBranch("see_q");
		see_q_branch->SetAddress(&see_q_);
	}
	if(see_q_branch == 0 ) {
	cout << "Branch see_q does not exist." << endl;
	}
	see_nValid_branch = 0;
	if (tree->GetBranch("see_nValid") != 0) {
		see_nValid_branch = tree->GetBranch("see_nValid");
		see_nValid_branch->SetAddress(&see_nValid_);
	}
	if(see_nValid_branch == 0 ) {
	cout << "Branch see_nValid does not exist." << endl;
	}
	see_nPixel_branch = 0;
	if (tree->GetBranch("see_nPixel") != 0) {
		see_nPixel_branch = tree->GetBranch("see_nPixel");
		see_nPixel_branch->SetAddress(&see_nPixel_);
	}
	if(see_nPixel_branch == 0 ) {
	cout << "Branch see_nPixel does not exist." << endl;
	}
	see_nGlued_branch = 0;
	if (tree->GetBranch("see_nGlued") != 0) {
		see_nGlued_branch = tree->GetBranch("see_nGlued");
		see_nGlued_branch->SetAddress(&see_nGlued_);
	}
	if(see_nGlued_branch == 0 ) {
	cout << "Branch see_nGlued does not exist." << endl;
	}
	see_nStrip_branch = 0;
	if (tree->GetBranch("see_nStrip") != 0) {
		see_nStrip_branch = tree->GetBranch("see_nStrip");
		see_nStrip_branch->SetAddress(&see_nStrip_);
	}
	if(see_nStrip_branch == 0 ) {
	cout << "Branch see_nStrip does not exist." << endl;
	}
	see_algo_branch = 0;
	if (tree->GetBranch("see_algo") != 0) {
		see_algo_branch = tree->GetBranch("see_algo");
		see_algo_branch->SetAddress(&see_algo_);
	}
	if(see_algo_branch == 0 ) {
	cout << "Branch see_algo does not exist." << endl;
	}
	see_pixelIdx_branch = 0;
	if (tree->GetBranch("see_pixelIdx") != 0) {
		see_pixelIdx_branch = tree->GetBranch("see_pixelIdx");
		see_pixelIdx_branch->SetAddress(&see_pixelIdx_);
	}
	if(see_pixelIdx_branch == 0 ) {
	cout << "Branch see_pixelIdx does not exist." << endl;
	}
	see_gluedIdx_branch = 0;
	if (tree->GetBranch("see_gluedIdx") != 0) {
		see_gluedIdx_branch = tree->GetBranch("see_gluedIdx");
		see_gluedIdx_branch->SetAddress(&see_gluedIdx_);
	}
	if(see_gluedIdx_branch == 0 ) {
	cout << "Branch see_gluedIdx does not exist." << endl;
	}
	see_stripIdx_branch = 0;
	if (tree->GetBranch("see_stripIdx") != 0) {
		see_stripIdx_branch = tree->GetBranch("see_stripIdx");
		see_stripIdx_branch->SetAddress(&see_stripIdx_);
	}
	if(see_stripIdx_branch == 0 ) {
	cout << "Branch see_stripIdx does not exist." << endl;
	}
	algo_offset_branch = 0;
	if (tree->GetBranch("algo_offset") != 0) {
		algo_offset_branch = tree->GetBranch("algo_offset");
		algo_offset_branch->SetAddress(&algo_offset_);
	}
	if(algo_offset_branch == 0 ) {
	cout << "Branch algo_offset does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		trk_px_isLoaded = false;
		trk_py_isLoaded = false;
		trk_pz_isLoaded = false;
		trk_pt_isLoaded = false;
		trk_eta_isLoaded = false;
		trk_phi_isLoaded = false;
		trk_dxy_isLoaded = false;
		trk_dz_isLoaded = false;
		trk_ptErr_isLoaded = false;
		trk_etaErr_isLoaded = false;
		trk_phiErr_isLoaded = false;
		trk_dxyErr_isLoaded = false;
		trk_dzErr_isLoaded = false;
		trk_nChi2_isLoaded = false;
		trk_shareFrac_isLoaded = false;
		trk_q_isLoaded = false;
		trk_nValid_isLoaded = false;
		trk_nInvalid_isLoaded = false;
		trk_nPixel_isLoaded = false;
		trk_nStrip_isLoaded = false;
		trk_n3DLay_isLoaded = false;
		trk_algo_isLoaded = false;
		trk_isHP_isLoaded = false;
		trk_seedIdx_isLoaded = false;
		trk_simIdx_isLoaded = false;
		trk_pixelIdx_isLoaded = false;
		trk_stripIdx_isLoaded = false;
		sim_px_isLoaded = false;
		sim_py_isLoaded = false;
		sim_pz_isLoaded = false;
		sim_pt_isLoaded = false;
		sim_eta_isLoaded = false;
		sim_phi_isLoaded = false;
		sim_dxy_isLoaded = false;
		sim_dz_isLoaded = false;
		sim_prodx_isLoaded = false;
		sim_prody_isLoaded = false;
		sim_prodz_isLoaded = false;
		sim_shareFrac_isLoaded = false;
		sim_q_isLoaded = false;
		sim_nValid_isLoaded = false;
		sim_nPixel_isLoaded = false;
		sim_nStrip_isLoaded = false;
		sim_n3DLay_isLoaded = false;
		sim_trkIdx_isLoaded = false;
		sim_pixelIdx_isLoaded = false;
		sim_stripIdx_isLoaded = false;
		pix_isBarrel_isLoaded = false;
		pix_lay_isLoaded = false;
		pix_detId_isLoaded = false;
		pix_nSimTrk_isLoaded = false;
		pix_simTrkIdx_isLoaded = false;
		pix_particle_isLoaded = false;
		pix_process_isLoaded = false;
		pix_bunchXing_isLoaded = false;
		pix_event_isLoaded = false;
		pix_x_isLoaded = false;
		pix_y_isLoaded = false;
		pix_z_isLoaded = false;
		pix_xx_isLoaded = false;
		pix_xy_isLoaded = false;
		pix_yy_isLoaded = false;
		pix_yz_isLoaded = false;
		pix_zz_isLoaded = false;
		pix_zx_isLoaded = false;
		pix_xsim_isLoaded = false;
		pix_ysim_isLoaded = false;
		pix_zsim_isLoaded = false;
		pix_pxsim_isLoaded = false;
		pix_pysim_isLoaded = false;
		pix_pzsim_isLoaded = false;
		pix_pathprop_isLoaded = false;
		pix_xsimprop_isLoaded = false;
		pix_ysimprop_isLoaded = false;
		pix_zsimprop_isLoaded = false;
		pix_pxsimprop_isLoaded = false;
		pix_pysimprop_isLoaded = false;
		pix_pzsimprop_isLoaded = false;
		pix_eloss_isLoaded = false;
		pix_radL_isLoaded = false;
		pix_bbxi_isLoaded = false;
		str_isBarrel_isLoaded = false;
		str_isStereo_isLoaded = false;
		str_det_isLoaded = false;
		str_lay_isLoaded = false;
		str_detId_isLoaded = false;
		str_nSimTrk_isLoaded = false;
		str_simTrkIdx_isLoaded = false;
		str_particle_isLoaded = false;
		str_process_isLoaded = false;
		str_bunchXing_isLoaded = false;
		str_event_isLoaded = false;
		str_x_isLoaded = false;
		str_y_isLoaded = false;
		str_z_isLoaded = false;
		str_xx_isLoaded = false;
		str_xy_isLoaded = false;
		str_yy_isLoaded = false;
		str_yz_isLoaded = false;
		str_zz_isLoaded = false;
		str_zx_isLoaded = false;
		str_xsim_isLoaded = false;
		str_ysim_isLoaded = false;
		str_zsim_isLoaded = false;
		str_pxsim_isLoaded = false;
		str_pysim_isLoaded = false;
		str_pzsim_isLoaded = false;
		str_eloss_isLoaded = false;
		str_radL_isLoaded = false;
		str_bbxi_isLoaded = false;
		glu_isBarrel_isLoaded = false;
		glu_det_isLoaded = false;
		glu_lay_isLoaded = false;
		glu_detId_isLoaded = false;
		glu_monoIdx_isLoaded = false;
		glu_stereoIdx_isLoaded = false;
		glu_x_isLoaded = false;
		glu_y_isLoaded = false;
		glu_z_isLoaded = false;
		glu_xx_isLoaded = false;
		glu_xy_isLoaded = false;
		glu_yy_isLoaded = false;
		glu_yz_isLoaded = false;
		glu_zz_isLoaded = false;
		glu_zx_isLoaded = false;
		glu_radL_isLoaded = false;
		glu_bbxi_isLoaded = false;
		bsp_x_isLoaded = false;
		bsp_y_isLoaded = false;
		bsp_z_isLoaded = false;
		bsp_sigmax_isLoaded = false;
		bsp_sigmay_isLoaded = false;
		bsp_sigmaz_isLoaded = false;
		see_px_isLoaded = false;
		see_py_isLoaded = false;
		see_pz_isLoaded = false;
		see_pt_isLoaded = false;
		see_eta_isLoaded = false;
		see_phi_isLoaded = false;
		see_dxy_isLoaded = false;
		see_dz_isLoaded = false;
		see_ptErr_isLoaded = false;
		see_etaErr_isLoaded = false;
		see_phiErr_isLoaded = false;
		see_dxyErr_isLoaded = false;
		see_dzErr_isLoaded = false;
		see_chi2_isLoaded = false;
		see_q_isLoaded = false;
		see_nValid_isLoaded = false;
		see_nPixel_isLoaded = false;
		see_nGlued_isLoaded = false;
		see_nStrip_isLoaded = false;
		see_algo_isLoaded = false;
		see_pixelIdx_isLoaded = false;
		see_gluedIdx_isLoaded = false;
		see_stripIdx_isLoaded = false;
		algo_offset_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (trk_px_branch != 0) trk_px();
	if (trk_py_branch != 0) trk_py();
	if (trk_pz_branch != 0) trk_pz();
	if (trk_pt_branch != 0) trk_pt();
	if (trk_eta_branch != 0) trk_eta();
	if (trk_phi_branch != 0) trk_phi();
	if (trk_dxy_branch != 0) trk_dxy();
	if (trk_dz_branch != 0) trk_dz();
	if (trk_ptErr_branch != 0) trk_ptErr();
	if (trk_etaErr_branch != 0) trk_etaErr();
	if (trk_phiErr_branch != 0) trk_phiErr();
	if (trk_dxyErr_branch != 0) trk_dxyErr();
	if (trk_dzErr_branch != 0) trk_dzErr();
	if (trk_nChi2_branch != 0) trk_nChi2();
	if (trk_shareFrac_branch != 0) trk_shareFrac();
	if (trk_q_branch != 0) trk_q();
	if (trk_nValid_branch != 0) trk_nValid();
	if (trk_nInvalid_branch != 0) trk_nInvalid();
	if (trk_nPixel_branch != 0) trk_nPixel();
	if (trk_nStrip_branch != 0) trk_nStrip();
	if (trk_n3DLay_branch != 0) trk_n3DLay();
	if (trk_algo_branch != 0) trk_algo();
	if (trk_isHP_branch != 0) trk_isHP();
	if (trk_seedIdx_branch != 0) trk_seedIdx();
	if (trk_simIdx_branch != 0) trk_simIdx();
	if (trk_pixelIdx_branch != 0) trk_pixelIdx();
	if (trk_stripIdx_branch != 0) trk_stripIdx();
	if (sim_px_branch != 0) sim_px();
	if (sim_py_branch != 0) sim_py();
	if (sim_pz_branch != 0) sim_pz();
	if (sim_pt_branch != 0) sim_pt();
	if (sim_eta_branch != 0) sim_eta();
	if (sim_phi_branch != 0) sim_phi();
	if (sim_dxy_branch != 0) sim_dxy();
	if (sim_dz_branch != 0) sim_dz();
	if (sim_prodx_branch != 0) sim_prodx();
	if (sim_prody_branch != 0) sim_prody();
	if (sim_prodz_branch != 0) sim_prodz();
	if (sim_shareFrac_branch != 0) sim_shareFrac();
	if (sim_q_branch != 0) sim_q();
	if (sim_nValid_branch != 0) sim_nValid();
	if (sim_nPixel_branch != 0) sim_nPixel();
	if (sim_nStrip_branch != 0) sim_nStrip();
	if (sim_n3DLay_branch != 0) sim_n3DLay();
	if (sim_trkIdx_branch != 0) sim_trkIdx();
	if (sim_pixelIdx_branch != 0) sim_pixelIdx();
	if (sim_stripIdx_branch != 0) sim_stripIdx();
	if (pix_isBarrel_branch != 0) pix_isBarrel();
	if (pix_lay_branch != 0) pix_lay();
	if (pix_detId_branch != 0) pix_detId();
	if (pix_nSimTrk_branch != 0) pix_nSimTrk();
	if (pix_simTrkIdx_branch != 0) pix_simTrkIdx();
	if (pix_particle_branch != 0) pix_particle();
	if (pix_process_branch != 0) pix_process();
	if (pix_bunchXing_branch != 0) pix_bunchXing();
	if (pix_event_branch != 0) pix_event();
	if (pix_x_branch != 0) pix_x();
	if (pix_y_branch != 0) pix_y();
	if (pix_z_branch != 0) pix_z();
	if (pix_xx_branch != 0) pix_xx();
	if (pix_xy_branch != 0) pix_xy();
	if (pix_yy_branch != 0) pix_yy();
	if (pix_yz_branch != 0) pix_yz();
	if (pix_zz_branch != 0) pix_zz();
	if (pix_zx_branch != 0) pix_zx();
	if (pix_xsim_branch != 0) pix_xsim();
	if (pix_ysim_branch != 0) pix_ysim();
	if (pix_zsim_branch != 0) pix_zsim();
	if (pix_pxsim_branch != 0) pix_pxsim();
	if (pix_pysim_branch != 0) pix_pysim();
	if (pix_pzsim_branch != 0) pix_pzsim();
	if (pix_pathprop_branch != 0) pix_pathprop();
	if (pix_xsimprop_branch != 0) pix_xsimprop();
	if (pix_ysimprop_branch != 0) pix_ysimprop();
	if (pix_zsimprop_branch != 0) pix_zsimprop();
	if (pix_pxsimprop_branch != 0) pix_pxsimprop();
	if (pix_pysimprop_branch != 0) pix_pysimprop();
	if (pix_pzsimprop_branch != 0) pix_pzsimprop();
	if (pix_eloss_branch != 0) pix_eloss();
	if (pix_radL_branch != 0) pix_radL();
	if (pix_bbxi_branch != 0) pix_bbxi();
	if (str_isBarrel_branch != 0) str_isBarrel();
	if (str_isStereo_branch != 0) str_isStereo();
	if (str_det_branch != 0) str_det();
	if (str_lay_branch != 0) str_lay();
	if (str_detId_branch != 0) str_detId();
	if (str_nSimTrk_branch != 0) str_nSimTrk();
	if (str_simTrkIdx_branch != 0) str_simTrkIdx();
	if (str_particle_branch != 0) str_particle();
	if (str_process_branch != 0) str_process();
	if (str_bunchXing_branch != 0) str_bunchXing();
	if (str_event_branch != 0) str_event();
	if (str_x_branch != 0) str_x();
	if (str_y_branch != 0) str_y();
	if (str_z_branch != 0) str_z();
	if (str_xx_branch != 0) str_xx();
	if (str_xy_branch != 0) str_xy();
	if (str_yy_branch != 0) str_yy();
	if (str_yz_branch != 0) str_yz();
	if (str_zz_branch != 0) str_zz();
	if (str_zx_branch != 0) str_zx();
	if (str_xsim_branch != 0) str_xsim();
	if (str_ysim_branch != 0) str_ysim();
	if (str_zsim_branch != 0) str_zsim();
	if (str_pxsim_branch != 0) str_pxsim();
	if (str_pysim_branch != 0) str_pysim();
	if (str_pzsim_branch != 0) str_pzsim();
	if (str_eloss_branch != 0) str_eloss();
	if (str_radL_branch != 0) str_radL();
	if (str_bbxi_branch != 0) str_bbxi();
	if (glu_isBarrel_branch != 0) glu_isBarrel();
	if (glu_det_branch != 0) glu_det();
	if (glu_lay_branch != 0) glu_lay();
	if (glu_detId_branch != 0) glu_detId();
	if (glu_monoIdx_branch != 0) glu_monoIdx();
	if (glu_stereoIdx_branch != 0) glu_stereoIdx();
	if (glu_x_branch != 0) glu_x();
	if (glu_y_branch != 0) glu_y();
	if (glu_z_branch != 0) glu_z();
	if (glu_xx_branch != 0) glu_xx();
	if (glu_xy_branch != 0) glu_xy();
	if (glu_yy_branch != 0) glu_yy();
	if (glu_yz_branch != 0) glu_yz();
	if (glu_zz_branch != 0) glu_zz();
	if (glu_zx_branch != 0) glu_zx();
	if (glu_radL_branch != 0) glu_radL();
	if (glu_bbxi_branch != 0) glu_bbxi();
	if (bsp_x_branch != 0) bsp_x();
	if (bsp_y_branch != 0) bsp_y();
	if (bsp_z_branch != 0) bsp_z();
	if (bsp_sigmax_branch != 0) bsp_sigmax();
	if (bsp_sigmay_branch != 0) bsp_sigmay();
	if (bsp_sigmaz_branch != 0) bsp_sigmaz();
	if (see_px_branch != 0) see_px();
	if (see_py_branch != 0) see_py();
	if (see_pz_branch != 0) see_pz();
	if (see_pt_branch != 0) see_pt();
	if (see_eta_branch != 0) see_eta();
	if (see_phi_branch != 0) see_phi();
	if (see_dxy_branch != 0) see_dxy();
	if (see_dz_branch != 0) see_dz();
	if (see_ptErr_branch != 0) see_ptErr();
	if (see_etaErr_branch != 0) see_etaErr();
	if (see_phiErr_branch != 0) see_phiErr();
	if (see_dxyErr_branch != 0) see_dxyErr();
	if (see_dzErr_branch != 0) see_dzErr();
	if (see_chi2_branch != 0) see_chi2();
	if (see_q_branch != 0) see_q();
	if (see_nValid_branch != 0) see_nValid();
	if (see_nPixel_branch != 0) see_nPixel();
	if (see_nGlued_branch != 0) see_nGlued();
	if (see_nStrip_branch != 0) see_nStrip();
	if (see_algo_branch != 0) see_algo();
	if (see_pixelIdx_branch != 0) see_pixelIdx();
	if (see_gluedIdx_branch != 0) see_gluedIdx();
	if (see_stripIdx_branch != 0) see_stripIdx();
	if (algo_offset_branch != 0) algo_offset();
}

	vector<float> &trk_px()
	{
		if (not trk_px_isLoaded) {
			if (trk_px_branch != 0) {
				trk_px_branch->GetEntry(index);
			} else { 
				printf("branch trk_px_branch does not exist!\n");
				exit(1);
			}
			trk_px_isLoaded = true;
		}
		return *trk_px_;
	}
	vector<float> &trk_py()
	{
		if (not trk_py_isLoaded) {
			if (trk_py_branch != 0) {
				trk_py_branch->GetEntry(index);
			} else { 
				printf("branch trk_py_branch does not exist!\n");
				exit(1);
			}
			trk_py_isLoaded = true;
		}
		return *trk_py_;
	}
	vector<float> &trk_pz()
	{
		if (not trk_pz_isLoaded) {
			if (trk_pz_branch != 0) {
				trk_pz_branch->GetEntry(index);
			} else { 
				printf("branch trk_pz_branch does not exist!\n");
				exit(1);
			}
			trk_pz_isLoaded = true;
		}
		return *trk_pz_;
	}
	vector<float> &trk_pt()
	{
		if (not trk_pt_isLoaded) {
			if (trk_pt_branch != 0) {
				trk_pt_branch->GetEntry(index);
			} else { 
				printf("branch trk_pt_branch does not exist!\n");
				exit(1);
			}
			trk_pt_isLoaded = true;
		}
		return *trk_pt_;
	}
	vector<float> &trk_eta()
	{
		if (not trk_eta_isLoaded) {
			if (trk_eta_branch != 0) {
				trk_eta_branch->GetEntry(index);
			} else { 
				printf("branch trk_eta_branch does not exist!\n");
				exit(1);
			}
			trk_eta_isLoaded = true;
		}
		return *trk_eta_;
	}
	vector<float> &trk_phi()
	{
		if (not trk_phi_isLoaded) {
			if (trk_phi_branch != 0) {
				trk_phi_branch->GetEntry(index);
			} else { 
				printf("branch trk_phi_branch does not exist!\n");
				exit(1);
			}
			trk_phi_isLoaded = true;
		}
		return *trk_phi_;
	}
	vector<float> &trk_dxy()
	{
		if (not trk_dxy_isLoaded) {
			if (trk_dxy_branch != 0) {
				trk_dxy_branch->GetEntry(index);
			} else { 
				printf("branch trk_dxy_branch does not exist!\n");
				exit(1);
			}
			trk_dxy_isLoaded = true;
		}
		return *trk_dxy_;
	}
	vector<float> &trk_dz()
	{
		if (not trk_dz_isLoaded) {
			if (trk_dz_branch != 0) {
				trk_dz_branch->GetEntry(index);
			} else { 
				printf("branch trk_dz_branch does not exist!\n");
				exit(1);
			}
			trk_dz_isLoaded = true;
		}
		return *trk_dz_;
	}
	vector<float> &trk_ptErr()
	{
		if (not trk_ptErr_isLoaded) {
			if (trk_ptErr_branch != 0) {
				trk_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch trk_ptErr_branch does not exist!\n");
				exit(1);
			}
			trk_ptErr_isLoaded = true;
		}
		return *trk_ptErr_;
	}
	vector<float> &trk_etaErr()
	{
		if (not trk_etaErr_isLoaded) {
			if (trk_etaErr_branch != 0) {
				trk_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch trk_etaErr_branch does not exist!\n");
				exit(1);
			}
			trk_etaErr_isLoaded = true;
		}
		return *trk_etaErr_;
	}
	vector<float> &trk_phiErr()
	{
		if (not trk_phiErr_isLoaded) {
			if (trk_phiErr_branch != 0) {
				trk_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch trk_phiErr_branch does not exist!\n");
				exit(1);
			}
			trk_phiErr_isLoaded = true;
		}
		return *trk_phiErr_;
	}
	vector<float> &trk_dxyErr()
	{
		if (not trk_dxyErr_isLoaded) {
			if (trk_dxyErr_branch != 0) {
				trk_dxyErr_branch->GetEntry(index);
			} else { 
				printf("branch trk_dxyErr_branch does not exist!\n");
				exit(1);
			}
			trk_dxyErr_isLoaded = true;
		}
		return *trk_dxyErr_;
	}
	vector<float> &trk_dzErr()
	{
		if (not trk_dzErr_isLoaded) {
			if (trk_dzErr_branch != 0) {
				trk_dzErr_branch->GetEntry(index);
			} else { 
				printf("branch trk_dzErr_branch does not exist!\n");
				exit(1);
			}
			trk_dzErr_isLoaded = true;
		}
		return *trk_dzErr_;
	}
	vector<float> &trk_nChi2()
	{
		if (not trk_nChi2_isLoaded) {
			if (trk_nChi2_branch != 0) {
				trk_nChi2_branch->GetEntry(index);
			} else { 
				printf("branch trk_nChi2_branch does not exist!\n");
				exit(1);
			}
			trk_nChi2_isLoaded = true;
		}
		return *trk_nChi2_;
	}
	vector<float> &trk_shareFrac()
	{
		if (not trk_shareFrac_isLoaded) {
			if (trk_shareFrac_branch != 0) {
				trk_shareFrac_branch->GetEntry(index);
			} else { 
				printf("branch trk_shareFrac_branch does not exist!\n");
				exit(1);
			}
			trk_shareFrac_isLoaded = true;
		}
		return *trk_shareFrac_;
	}
	vector<int> &trk_q()
	{
		if (not trk_q_isLoaded) {
			if (trk_q_branch != 0) {
				trk_q_branch->GetEntry(index);
			} else { 
				printf("branch trk_q_branch does not exist!\n");
				exit(1);
			}
			trk_q_isLoaded = true;
		}
		return *trk_q_;
	}
	vector<int> &trk_nValid()
	{
		if (not trk_nValid_isLoaded) {
			if (trk_nValid_branch != 0) {
				trk_nValid_branch->GetEntry(index);
			} else { 
				printf("branch trk_nValid_branch does not exist!\n");
				exit(1);
			}
			trk_nValid_isLoaded = true;
		}
		return *trk_nValid_;
	}
	vector<int> &trk_nInvalid()
	{
		if (not trk_nInvalid_isLoaded) {
			if (trk_nInvalid_branch != 0) {
				trk_nInvalid_branch->GetEntry(index);
			} else { 
				printf("branch trk_nInvalid_branch does not exist!\n");
				exit(1);
			}
			trk_nInvalid_isLoaded = true;
		}
		return *trk_nInvalid_;
	}
	vector<int> &trk_nPixel()
	{
		if (not trk_nPixel_isLoaded) {
			if (trk_nPixel_branch != 0) {
				trk_nPixel_branch->GetEntry(index);
			} else { 
				printf("branch trk_nPixel_branch does not exist!\n");
				exit(1);
			}
			trk_nPixel_isLoaded = true;
		}
		return *trk_nPixel_;
	}
	vector<int> &trk_nStrip()
	{
		if (not trk_nStrip_isLoaded) {
			if (trk_nStrip_branch != 0) {
				trk_nStrip_branch->GetEntry(index);
			} else { 
				printf("branch trk_nStrip_branch does not exist!\n");
				exit(1);
			}
			trk_nStrip_isLoaded = true;
		}
		return *trk_nStrip_;
	}
	vector<int> &trk_n3DLay()
	{
		if (not trk_n3DLay_isLoaded) {
			if (trk_n3DLay_branch != 0) {
				trk_n3DLay_branch->GetEntry(index);
			} else { 
				printf("branch trk_n3DLay_branch does not exist!\n");
				exit(1);
			}
			trk_n3DLay_isLoaded = true;
		}
		return *trk_n3DLay_;
	}
	vector<int> &trk_algo()
	{
		if (not trk_algo_isLoaded) {
			if (trk_algo_branch != 0) {
				trk_algo_branch->GetEntry(index);
			} else { 
				printf("branch trk_algo_branch does not exist!\n");
				exit(1);
			}
			trk_algo_isLoaded = true;
		}
		return *trk_algo_;
	}
	vector<int> &trk_isHP()
	{
		if (not trk_isHP_isLoaded) {
			if (trk_isHP_branch != 0) {
				trk_isHP_branch->GetEntry(index);
			} else { 
				printf("branch trk_isHP_branch does not exist!\n");
				exit(1);
			}
			trk_isHP_isLoaded = true;
		}
		return *trk_isHP_;
	}
	vector<int> &trk_seedIdx()
	{
		if (not trk_seedIdx_isLoaded) {
			if (trk_seedIdx_branch != 0) {
				trk_seedIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_seedIdx_branch does not exist!\n");
				exit(1);
			}
			trk_seedIdx_isLoaded = true;
		}
		return *trk_seedIdx_;
	}
	vector<int> &trk_simIdx()
	{
		if (not trk_simIdx_isLoaded) {
			if (trk_simIdx_branch != 0) {
				trk_simIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_simIdx_branch does not exist!\n");
				exit(1);
			}
			trk_simIdx_isLoaded = true;
		}
		return *trk_simIdx_;
	}
	vector<vector<int> > &trk_pixelIdx()
	{
		if (not trk_pixelIdx_isLoaded) {
			if (trk_pixelIdx_branch != 0) {
				trk_pixelIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_pixelIdx_branch does not exist!\n");
				exit(1);
			}
			trk_pixelIdx_isLoaded = true;
		}
		return *trk_pixelIdx_;
	}
	vector<vector<int> > &trk_stripIdx()
	{
		if (not trk_stripIdx_isLoaded) {
			if (trk_stripIdx_branch != 0) {
				trk_stripIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_stripIdx_branch does not exist!\n");
				exit(1);
			}
			trk_stripIdx_isLoaded = true;
		}
		return *trk_stripIdx_;
	}
	vector<float> &sim_px()
	{
		if (not sim_px_isLoaded) {
			if (sim_px_branch != 0) {
				sim_px_branch->GetEntry(index);
			} else { 
				printf("branch sim_px_branch does not exist!\n");
				exit(1);
			}
			sim_px_isLoaded = true;
		}
		return *sim_px_;
	}
	vector<float> &sim_py()
	{
		if (not sim_py_isLoaded) {
			if (sim_py_branch != 0) {
				sim_py_branch->GetEntry(index);
			} else { 
				printf("branch sim_py_branch does not exist!\n");
				exit(1);
			}
			sim_py_isLoaded = true;
		}
		return *sim_py_;
	}
	vector<float> &sim_pz()
	{
		if (not sim_pz_isLoaded) {
			if (sim_pz_branch != 0) {
				sim_pz_branch->GetEntry(index);
			} else { 
				printf("branch sim_pz_branch does not exist!\n");
				exit(1);
			}
			sim_pz_isLoaded = true;
		}
		return *sim_pz_;
	}
	vector<float> &sim_pt()
	{
		if (not sim_pt_isLoaded) {
			if (sim_pt_branch != 0) {
				sim_pt_branch->GetEntry(index);
			} else { 
				printf("branch sim_pt_branch does not exist!\n");
				exit(1);
			}
			sim_pt_isLoaded = true;
		}
		return *sim_pt_;
	}
	vector<float> &sim_eta()
	{
		if (not sim_eta_isLoaded) {
			if (sim_eta_branch != 0) {
				sim_eta_branch->GetEntry(index);
			} else { 
				printf("branch sim_eta_branch does not exist!\n");
				exit(1);
			}
			sim_eta_isLoaded = true;
		}
		return *sim_eta_;
	}
	vector<float> &sim_phi()
	{
		if (not sim_phi_isLoaded) {
			if (sim_phi_branch != 0) {
				sim_phi_branch->GetEntry(index);
			} else { 
				printf("branch sim_phi_branch does not exist!\n");
				exit(1);
			}
			sim_phi_isLoaded = true;
		}
		return *sim_phi_;
	}
	vector<float> &sim_dxy()
	{
		if (not sim_dxy_isLoaded) {
			if (sim_dxy_branch != 0) {
				sim_dxy_branch->GetEntry(index);
			} else { 
				printf("branch sim_dxy_branch does not exist!\n");
				exit(1);
			}
			sim_dxy_isLoaded = true;
		}
		return *sim_dxy_;
	}
	vector<float> &sim_dz()
	{
		if (not sim_dz_isLoaded) {
			if (sim_dz_branch != 0) {
				sim_dz_branch->GetEntry(index);
			} else { 
				printf("branch sim_dz_branch does not exist!\n");
				exit(1);
			}
			sim_dz_isLoaded = true;
		}
		return *sim_dz_;
	}
	vector<float> &sim_prodx()
	{
		if (not sim_prodx_isLoaded) {
			if (sim_prodx_branch != 0) {
				sim_prodx_branch->GetEntry(index);
			} else { 
				printf("branch sim_prodx_branch does not exist!\n");
				exit(1);
			}
			sim_prodx_isLoaded = true;
		}
		return *sim_prodx_;
	}
	vector<float> &sim_prody()
	{
		if (not sim_prody_isLoaded) {
			if (sim_prody_branch != 0) {
				sim_prody_branch->GetEntry(index);
			} else { 
				printf("branch sim_prody_branch does not exist!\n");
				exit(1);
			}
			sim_prody_isLoaded = true;
		}
		return *sim_prody_;
	}
	vector<float> &sim_prodz()
	{
		if (not sim_prodz_isLoaded) {
			if (sim_prodz_branch != 0) {
				sim_prodz_branch->GetEntry(index);
			} else { 
				printf("branch sim_prodz_branch does not exist!\n");
				exit(1);
			}
			sim_prodz_isLoaded = true;
		}
		return *sim_prodz_;
	}
	vector<float> &sim_shareFrac()
	{
		if (not sim_shareFrac_isLoaded) {
			if (sim_shareFrac_branch != 0) {
				sim_shareFrac_branch->GetEntry(index);
			} else { 
				printf("branch sim_shareFrac_branch does not exist!\n");
				exit(1);
			}
			sim_shareFrac_isLoaded = true;
		}
		return *sim_shareFrac_;
	}
	vector<int> &sim_q()
	{
		if (not sim_q_isLoaded) {
			if (sim_q_branch != 0) {
				sim_q_branch->GetEntry(index);
			} else { 
				printf("branch sim_q_branch does not exist!\n");
				exit(1);
			}
			sim_q_isLoaded = true;
		}
		return *sim_q_;
	}
	vector<int> &sim_nValid()
	{
		if (not sim_nValid_isLoaded) {
			if (sim_nValid_branch != 0) {
				sim_nValid_branch->GetEntry(index);
			} else { 
				printf("branch sim_nValid_branch does not exist!\n");
				exit(1);
			}
			sim_nValid_isLoaded = true;
		}
		return *sim_nValid_;
	}
	vector<int> &sim_nPixel()
	{
		if (not sim_nPixel_isLoaded) {
			if (sim_nPixel_branch != 0) {
				sim_nPixel_branch->GetEntry(index);
			} else { 
				printf("branch sim_nPixel_branch does not exist!\n");
				exit(1);
			}
			sim_nPixel_isLoaded = true;
		}
		return *sim_nPixel_;
	}
	vector<int> &sim_nStrip()
	{
		if (not sim_nStrip_isLoaded) {
			if (sim_nStrip_branch != 0) {
				sim_nStrip_branch->GetEntry(index);
			} else { 
				printf("branch sim_nStrip_branch does not exist!\n");
				exit(1);
			}
			sim_nStrip_isLoaded = true;
		}
		return *sim_nStrip_;
	}
	vector<int> &sim_n3DLay()
	{
		if (not sim_n3DLay_isLoaded) {
			if (sim_n3DLay_branch != 0) {
				sim_n3DLay_branch->GetEntry(index);
			} else { 
				printf("branch sim_n3DLay_branch does not exist!\n");
				exit(1);
			}
			sim_n3DLay_isLoaded = true;
		}
		return *sim_n3DLay_;
	}
	vector<int> &sim_trkIdx()
	{
		if (not sim_trkIdx_isLoaded) {
			if (sim_trkIdx_branch != 0) {
				sim_trkIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_trkIdx_branch does not exist!\n");
				exit(1);
			}
			sim_trkIdx_isLoaded = true;
		}
		return *sim_trkIdx_;
	}
	vector<vector<int> > &sim_pixelIdx()
	{
		if (not sim_pixelIdx_isLoaded) {
			if (sim_pixelIdx_branch != 0) {
				sim_pixelIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_pixelIdx_branch does not exist!\n");
				exit(1);
			}
			sim_pixelIdx_isLoaded = true;
		}
		return *sim_pixelIdx_;
	}
	vector<vector<int> > &sim_stripIdx()
	{
		if (not sim_stripIdx_isLoaded) {
			if (sim_stripIdx_branch != 0) {
				sim_stripIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_stripIdx_branch does not exist!\n");
				exit(1);
			}
			sim_stripIdx_isLoaded = true;
		}
		return *sim_stripIdx_;
	}
	vector<int> &pix_isBarrel()
	{
		if (not pix_isBarrel_isLoaded) {
			if (pix_isBarrel_branch != 0) {
				pix_isBarrel_branch->GetEntry(index);
			} else { 
				printf("branch pix_isBarrel_branch does not exist!\n");
				exit(1);
			}
			pix_isBarrel_isLoaded = true;
		}
		return *pix_isBarrel_;
	}
	vector<int> &pix_lay()
	{
		if (not pix_lay_isLoaded) {
			if (pix_lay_branch != 0) {
				pix_lay_branch->GetEntry(index);
			} else { 
				printf("branch pix_lay_branch does not exist!\n");
				exit(1);
			}
			pix_lay_isLoaded = true;
		}
		return *pix_lay_;
	}
	vector<int> &pix_detId()
	{
		if (not pix_detId_isLoaded) {
			if (pix_detId_branch != 0) {
				pix_detId_branch->GetEntry(index);
			} else { 
				printf("branch pix_detId_branch does not exist!\n");
				exit(1);
			}
			pix_detId_isLoaded = true;
		}
		return *pix_detId_;
	}
	vector<int> &pix_nSimTrk()
	{
		if (not pix_nSimTrk_isLoaded) {
			if (pix_nSimTrk_branch != 0) {
				pix_nSimTrk_branch->GetEntry(index);
			} else { 
				printf("branch pix_nSimTrk_branch does not exist!\n");
				exit(1);
			}
			pix_nSimTrk_isLoaded = true;
		}
		return *pix_nSimTrk_;
	}
	vector<int> &pix_simTrkIdx()
	{
		if (not pix_simTrkIdx_isLoaded) {
			if (pix_simTrkIdx_branch != 0) {
				pix_simTrkIdx_branch->GetEntry(index);
			} else { 
				printf("branch pix_simTrkIdx_branch does not exist!\n");
				exit(1);
			}
			pix_simTrkIdx_isLoaded = true;
		}
		return *pix_simTrkIdx_;
	}
	vector<int> &pix_particle()
	{
		if (not pix_particle_isLoaded) {
			if (pix_particle_branch != 0) {
				pix_particle_branch->GetEntry(index);
			} else { 
				printf("branch pix_particle_branch does not exist!\n");
				exit(1);
			}
			pix_particle_isLoaded = true;
		}
		return *pix_particle_;
	}
	vector<int> &pix_process()
	{
		if (not pix_process_isLoaded) {
			if (pix_process_branch != 0) {
				pix_process_branch->GetEntry(index);
			} else { 
				printf("branch pix_process_branch does not exist!\n");
				exit(1);
			}
			pix_process_isLoaded = true;
		}
		return *pix_process_;
	}
	vector<int> &pix_bunchXing()
	{
		if (not pix_bunchXing_isLoaded) {
			if (pix_bunchXing_branch != 0) {
				pix_bunchXing_branch->GetEntry(index);
			} else { 
				printf("branch pix_bunchXing_branch does not exist!\n");
				exit(1);
			}
			pix_bunchXing_isLoaded = true;
		}
		return *pix_bunchXing_;
	}
	vector<int> &pix_event()
	{
		if (not pix_event_isLoaded) {
			if (pix_event_branch != 0) {
				pix_event_branch->GetEntry(index);
			} else { 
				printf("branch pix_event_branch does not exist!\n");
				exit(1);
			}
			pix_event_isLoaded = true;
		}
		return *pix_event_;
	}
	vector<float> &pix_x()
	{
		if (not pix_x_isLoaded) {
			if (pix_x_branch != 0) {
				pix_x_branch->GetEntry(index);
			} else { 
				printf("branch pix_x_branch does not exist!\n");
				exit(1);
			}
			pix_x_isLoaded = true;
		}
		return *pix_x_;
	}
	vector<float> &pix_y()
	{
		if (not pix_y_isLoaded) {
			if (pix_y_branch != 0) {
				pix_y_branch->GetEntry(index);
			} else { 
				printf("branch pix_y_branch does not exist!\n");
				exit(1);
			}
			pix_y_isLoaded = true;
		}
		return *pix_y_;
	}
	vector<float> &pix_z()
	{
		if (not pix_z_isLoaded) {
			if (pix_z_branch != 0) {
				pix_z_branch->GetEntry(index);
			} else { 
				printf("branch pix_z_branch does not exist!\n");
				exit(1);
			}
			pix_z_isLoaded = true;
		}
		return *pix_z_;
	}
	vector<float> &pix_xx()
	{
		if (not pix_xx_isLoaded) {
			if (pix_xx_branch != 0) {
				pix_xx_branch->GetEntry(index);
			} else { 
				printf("branch pix_xx_branch does not exist!\n");
				exit(1);
			}
			pix_xx_isLoaded = true;
		}
		return *pix_xx_;
	}
	vector<float> &pix_xy()
	{
		if (not pix_xy_isLoaded) {
			if (pix_xy_branch != 0) {
				pix_xy_branch->GetEntry(index);
			} else { 
				printf("branch pix_xy_branch does not exist!\n");
				exit(1);
			}
			pix_xy_isLoaded = true;
		}
		return *pix_xy_;
	}
	vector<float> &pix_yy()
	{
		if (not pix_yy_isLoaded) {
			if (pix_yy_branch != 0) {
				pix_yy_branch->GetEntry(index);
			} else { 
				printf("branch pix_yy_branch does not exist!\n");
				exit(1);
			}
			pix_yy_isLoaded = true;
		}
		return *pix_yy_;
	}
	vector<float> &pix_yz()
	{
		if (not pix_yz_isLoaded) {
			if (pix_yz_branch != 0) {
				pix_yz_branch->GetEntry(index);
			} else { 
				printf("branch pix_yz_branch does not exist!\n");
				exit(1);
			}
			pix_yz_isLoaded = true;
		}
		return *pix_yz_;
	}
	vector<float> &pix_zz()
	{
		if (not pix_zz_isLoaded) {
			if (pix_zz_branch != 0) {
				pix_zz_branch->GetEntry(index);
			} else { 
				printf("branch pix_zz_branch does not exist!\n");
				exit(1);
			}
			pix_zz_isLoaded = true;
		}
		return *pix_zz_;
	}
	vector<float> &pix_zx()
	{
		if (not pix_zx_isLoaded) {
			if (pix_zx_branch != 0) {
				pix_zx_branch->GetEntry(index);
			} else { 
				printf("branch pix_zx_branch does not exist!\n");
				exit(1);
			}
			pix_zx_isLoaded = true;
		}
		return *pix_zx_;
	}
	vector<float> &pix_xsim()
	{
		if (not pix_xsim_isLoaded) {
			if (pix_xsim_branch != 0) {
				pix_xsim_branch->GetEntry(index);
			} else { 
				printf("branch pix_xsim_branch does not exist!\n");
				exit(1);
			}
			pix_xsim_isLoaded = true;
		}
		return *pix_xsim_;
	}
	vector<float> &pix_ysim()
	{
		if (not pix_ysim_isLoaded) {
			if (pix_ysim_branch != 0) {
				pix_ysim_branch->GetEntry(index);
			} else { 
				printf("branch pix_ysim_branch does not exist!\n");
				exit(1);
			}
			pix_ysim_isLoaded = true;
		}
		return *pix_ysim_;
	}
	vector<float> &pix_zsim()
	{
		if (not pix_zsim_isLoaded) {
			if (pix_zsim_branch != 0) {
				pix_zsim_branch->GetEntry(index);
			} else { 
				printf("branch pix_zsim_branch does not exist!\n");
				exit(1);
			}
			pix_zsim_isLoaded = true;
		}
		return *pix_zsim_;
	}
	vector<float> &pix_pxsim()
	{
		if (not pix_pxsim_isLoaded) {
			if (pix_pxsim_branch != 0) {
				pix_pxsim_branch->GetEntry(index);
			} else { 
				printf("branch pix_pxsim_branch does not exist!\n");
				exit(1);
			}
			pix_pxsim_isLoaded = true;
		}
		return *pix_pxsim_;
	}
	vector<float> &pix_pysim()
	{
		if (not pix_pysim_isLoaded) {
			if (pix_pysim_branch != 0) {
				pix_pysim_branch->GetEntry(index);
			} else { 
				printf("branch pix_pysim_branch does not exist!\n");
				exit(1);
			}
			pix_pysim_isLoaded = true;
		}
		return *pix_pysim_;
	}
	vector<float> &pix_pzsim()
	{
		if (not pix_pzsim_isLoaded) {
			if (pix_pzsim_branch != 0) {
				pix_pzsim_branch->GetEntry(index);
			} else { 
				printf("branch pix_pzsim_branch does not exist!\n");
				exit(1);
			}
			pix_pzsim_isLoaded = true;
		}
		return *pix_pzsim_;
	}
	vector<float> &pix_pathprop()
	{
		if (not pix_pathprop_isLoaded) {
			if (pix_pathprop_branch != 0) {
				pix_pathprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_pathprop_branch does not exist!\n");
				exit(1);
			}
			pix_pathprop_isLoaded = true;
		}
		return *pix_pathprop_;
	}
	vector<float> &pix_xsimprop()
	{
		if (not pix_xsimprop_isLoaded) {
			if (pix_xsimprop_branch != 0) {
				pix_xsimprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_xsimprop_branch does not exist!\n");
				exit(1);
			}
			pix_xsimprop_isLoaded = true;
		}
		return *pix_xsimprop_;
	}
	vector<float> &pix_ysimprop()
	{
		if (not pix_ysimprop_isLoaded) {
			if (pix_ysimprop_branch != 0) {
				pix_ysimprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_ysimprop_branch does not exist!\n");
				exit(1);
			}
			pix_ysimprop_isLoaded = true;
		}
		return *pix_ysimprop_;
	}
	vector<float> &pix_zsimprop()
	{
		if (not pix_zsimprop_isLoaded) {
			if (pix_zsimprop_branch != 0) {
				pix_zsimprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_zsimprop_branch does not exist!\n");
				exit(1);
			}
			pix_zsimprop_isLoaded = true;
		}
		return *pix_zsimprop_;
	}
	vector<float> &pix_pxsimprop()
	{
		if (not pix_pxsimprop_isLoaded) {
			if (pix_pxsimprop_branch != 0) {
				pix_pxsimprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_pxsimprop_branch does not exist!\n");
				exit(1);
			}
			pix_pxsimprop_isLoaded = true;
		}
		return *pix_pxsimprop_;
	}
	vector<float> &pix_pysimprop()
	{
		if (not pix_pysimprop_isLoaded) {
			if (pix_pysimprop_branch != 0) {
				pix_pysimprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_pysimprop_branch does not exist!\n");
				exit(1);
			}
			pix_pysimprop_isLoaded = true;
		}
		return *pix_pysimprop_;
	}
	vector<float> &pix_pzsimprop()
	{
		if (not pix_pzsimprop_isLoaded) {
			if (pix_pzsimprop_branch != 0) {
				pix_pzsimprop_branch->GetEntry(index);
			} else { 
				printf("branch pix_pzsimprop_branch does not exist!\n");
				exit(1);
			}
			pix_pzsimprop_isLoaded = true;
		}
		return *pix_pzsimprop_;
	}
	vector<float> &pix_eloss()
	{
		if (not pix_eloss_isLoaded) {
			if (pix_eloss_branch != 0) {
				pix_eloss_branch->GetEntry(index);
			} else { 
				printf("branch pix_eloss_branch does not exist!\n");
				exit(1);
			}
			pix_eloss_isLoaded = true;
		}
		return *pix_eloss_;
	}
	vector<float> &pix_radL()
	{
		if (not pix_radL_isLoaded) {
			if (pix_radL_branch != 0) {
				pix_radL_branch->GetEntry(index);
			} else { 
				printf("branch pix_radL_branch does not exist!\n");
				exit(1);
			}
			pix_radL_isLoaded = true;
		}
		return *pix_radL_;
	}
	vector<float> &pix_bbxi()
	{
		if (not pix_bbxi_isLoaded) {
			if (pix_bbxi_branch != 0) {
				pix_bbxi_branch->GetEntry(index);
			} else { 
				printf("branch pix_bbxi_branch does not exist!\n");
				exit(1);
			}
			pix_bbxi_isLoaded = true;
		}
		return *pix_bbxi_;
	}
	vector<int> &str_isBarrel()
	{
		if (not str_isBarrel_isLoaded) {
			if (str_isBarrel_branch != 0) {
				str_isBarrel_branch->GetEntry(index);
			} else { 
				printf("branch str_isBarrel_branch does not exist!\n");
				exit(1);
			}
			str_isBarrel_isLoaded = true;
		}
		return *str_isBarrel_;
	}
	vector<int> &str_isStereo()
	{
		if (not str_isStereo_isLoaded) {
			if (str_isStereo_branch != 0) {
				str_isStereo_branch->GetEntry(index);
			} else { 
				printf("branch str_isStereo_branch does not exist!\n");
				exit(1);
			}
			str_isStereo_isLoaded = true;
		}
		return *str_isStereo_;
	}
	vector<int> &str_det()
	{
		if (not str_det_isLoaded) {
			if (str_det_branch != 0) {
				str_det_branch->GetEntry(index);
			} else { 
				printf("branch str_det_branch does not exist!\n");
				exit(1);
			}
			str_det_isLoaded = true;
		}
		return *str_det_;
	}
	vector<int> &str_lay()
	{
		if (not str_lay_isLoaded) {
			if (str_lay_branch != 0) {
				str_lay_branch->GetEntry(index);
			} else { 
				printf("branch str_lay_branch does not exist!\n");
				exit(1);
			}
			str_lay_isLoaded = true;
		}
		return *str_lay_;
	}
	vector<int> &str_detId()
	{
		if (not str_detId_isLoaded) {
			if (str_detId_branch != 0) {
				str_detId_branch->GetEntry(index);
			} else { 
				printf("branch str_detId_branch does not exist!\n");
				exit(1);
			}
			str_detId_isLoaded = true;
		}
		return *str_detId_;
	}
	vector<int> &str_nSimTrk()
	{
		if (not str_nSimTrk_isLoaded) {
			if (str_nSimTrk_branch != 0) {
				str_nSimTrk_branch->GetEntry(index);
			} else { 
				printf("branch str_nSimTrk_branch does not exist!\n");
				exit(1);
			}
			str_nSimTrk_isLoaded = true;
		}
		return *str_nSimTrk_;
	}
	vector<int> &str_simTrkIdx()
	{
		if (not str_simTrkIdx_isLoaded) {
			if (str_simTrkIdx_branch != 0) {
				str_simTrkIdx_branch->GetEntry(index);
			} else { 
				printf("branch str_simTrkIdx_branch does not exist!\n");
				exit(1);
			}
			str_simTrkIdx_isLoaded = true;
		}
		return *str_simTrkIdx_;
	}
	vector<int> &str_particle()
	{
		if (not str_particle_isLoaded) {
			if (str_particle_branch != 0) {
				str_particle_branch->GetEntry(index);
			} else { 
				printf("branch str_particle_branch does not exist!\n");
				exit(1);
			}
			str_particle_isLoaded = true;
		}
		return *str_particle_;
	}
	vector<int> &str_process()
	{
		if (not str_process_isLoaded) {
			if (str_process_branch != 0) {
				str_process_branch->GetEntry(index);
			} else { 
				printf("branch str_process_branch does not exist!\n");
				exit(1);
			}
			str_process_isLoaded = true;
		}
		return *str_process_;
	}
	vector<int> &str_bunchXing()
	{
		if (not str_bunchXing_isLoaded) {
			if (str_bunchXing_branch != 0) {
				str_bunchXing_branch->GetEntry(index);
			} else { 
				printf("branch str_bunchXing_branch does not exist!\n");
				exit(1);
			}
			str_bunchXing_isLoaded = true;
		}
		return *str_bunchXing_;
	}
	vector<int> &str_event()
	{
		if (not str_event_isLoaded) {
			if (str_event_branch != 0) {
				str_event_branch->GetEntry(index);
			} else { 
				printf("branch str_event_branch does not exist!\n");
				exit(1);
			}
			str_event_isLoaded = true;
		}
		return *str_event_;
	}
	vector<float> &str_x()
	{
		if (not str_x_isLoaded) {
			if (str_x_branch != 0) {
				str_x_branch->GetEntry(index);
			} else { 
				printf("branch str_x_branch does not exist!\n");
				exit(1);
			}
			str_x_isLoaded = true;
		}
		return *str_x_;
	}
	vector<float> &str_y()
	{
		if (not str_y_isLoaded) {
			if (str_y_branch != 0) {
				str_y_branch->GetEntry(index);
			} else { 
				printf("branch str_y_branch does not exist!\n");
				exit(1);
			}
			str_y_isLoaded = true;
		}
		return *str_y_;
	}
	vector<float> &str_z()
	{
		if (not str_z_isLoaded) {
			if (str_z_branch != 0) {
				str_z_branch->GetEntry(index);
			} else { 
				printf("branch str_z_branch does not exist!\n");
				exit(1);
			}
			str_z_isLoaded = true;
		}
		return *str_z_;
	}
	vector<float> &str_xx()
	{
		if (not str_xx_isLoaded) {
			if (str_xx_branch != 0) {
				str_xx_branch->GetEntry(index);
			} else { 
				printf("branch str_xx_branch does not exist!\n");
				exit(1);
			}
			str_xx_isLoaded = true;
		}
		return *str_xx_;
	}
	vector<float> &str_xy()
	{
		if (not str_xy_isLoaded) {
			if (str_xy_branch != 0) {
				str_xy_branch->GetEntry(index);
			} else { 
				printf("branch str_xy_branch does not exist!\n");
				exit(1);
			}
			str_xy_isLoaded = true;
		}
		return *str_xy_;
	}
	vector<float> &str_yy()
	{
		if (not str_yy_isLoaded) {
			if (str_yy_branch != 0) {
				str_yy_branch->GetEntry(index);
			} else { 
				printf("branch str_yy_branch does not exist!\n");
				exit(1);
			}
			str_yy_isLoaded = true;
		}
		return *str_yy_;
	}
	vector<float> &str_yz()
	{
		if (not str_yz_isLoaded) {
			if (str_yz_branch != 0) {
				str_yz_branch->GetEntry(index);
			} else { 
				printf("branch str_yz_branch does not exist!\n");
				exit(1);
			}
			str_yz_isLoaded = true;
		}
		return *str_yz_;
	}
	vector<float> &str_zz()
	{
		if (not str_zz_isLoaded) {
			if (str_zz_branch != 0) {
				str_zz_branch->GetEntry(index);
			} else { 
				printf("branch str_zz_branch does not exist!\n");
				exit(1);
			}
			str_zz_isLoaded = true;
		}
		return *str_zz_;
	}
	vector<float> &str_zx()
	{
		if (not str_zx_isLoaded) {
			if (str_zx_branch != 0) {
				str_zx_branch->GetEntry(index);
			} else { 
				printf("branch str_zx_branch does not exist!\n");
				exit(1);
			}
			str_zx_isLoaded = true;
		}
		return *str_zx_;
	}
	vector<float> &str_xsim()
	{
		if (not str_xsim_isLoaded) {
			if (str_xsim_branch != 0) {
				str_xsim_branch->GetEntry(index);
			} else { 
				printf("branch str_xsim_branch does not exist!\n");
				exit(1);
			}
			str_xsim_isLoaded = true;
		}
		return *str_xsim_;
	}
	vector<float> &str_ysim()
	{
		if (not str_ysim_isLoaded) {
			if (str_ysim_branch != 0) {
				str_ysim_branch->GetEntry(index);
			} else { 
				printf("branch str_ysim_branch does not exist!\n");
				exit(1);
			}
			str_ysim_isLoaded = true;
		}
		return *str_ysim_;
	}
	vector<float> &str_zsim()
	{
		if (not str_zsim_isLoaded) {
			if (str_zsim_branch != 0) {
				str_zsim_branch->GetEntry(index);
			} else { 
				printf("branch str_zsim_branch does not exist!\n");
				exit(1);
			}
			str_zsim_isLoaded = true;
		}
		return *str_zsim_;
	}
	vector<float> &str_pxsim()
	{
		if (not str_pxsim_isLoaded) {
			if (str_pxsim_branch != 0) {
				str_pxsim_branch->GetEntry(index);
			} else { 
				printf("branch str_pxsim_branch does not exist!\n");
				exit(1);
			}
			str_pxsim_isLoaded = true;
		}
		return *str_pxsim_;
	}
	vector<float> &str_pysim()
	{
		if (not str_pysim_isLoaded) {
			if (str_pysim_branch != 0) {
				str_pysim_branch->GetEntry(index);
			} else { 
				printf("branch str_pysim_branch does not exist!\n");
				exit(1);
			}
			str_pysim_isLoaded = true;
		}
		return *str_pysim_;
	}
	vector<float> &str_pzsim()
	{
		if (not str_pzsim_isLoaded) {
			if (str_pzsim_branch != 0) {
				str_pzsim_branch->GetEntry(index);
			} else { 
				printf("branch str_pzsim_branch does not exist!\n");
				exit(1);
			}
			str_pzsim_isLoaded = true;
		}
		return *str_pzsim_;
	}
	vector<float> &str_eloss()
	{
		if (not str_eloss_isLoaded) {
			if (str_eloss_branch != 0) {
				str_eloss_branch->GetEntry(index);
			} else { 
				printf("branch str_eloss_branch does not exist!\n");
				exit(1);
			}
			str_eloss_isLoaded = true;
		}
		return *str_eloss_;
	}
	vector<float> &str_radL()
	{
		if (not str_radL_isLoaded) {
			if (str_radL_branch != 0) {
				str_radL_branch->GetEntry(index);
			} else { 
				printf("branch str_radL_branch does not exist!\n");
				exit(1);
			}
			str_radL_isLoaded = true;
		}
		return *str_radL_;
	}
	vector<float> &str_bbxi()
	{
		if (not str_bbxi_isLoaded) {
			if (str_bbxi_branch != 0) {
				str_bbxi_branch->GetEntry(index);
			} else { 
				printf("branch str_bbxi_branch does not exist!\n");
				exit(1);
			}
			str_bbxi_isLoaded = true;
		}
		return *str_bbxi_;
	}
	vector<int> &glu_isBarrel()
	{
		if (not glu_isBarrel_isLoaded) {
			if (glu_isBarrel_branch != 0) {
				glu_isBarrel_branch->GetEntry(index);
			} else { 
				printf("branch glu_isBarrel_branch does not exist!\n");
				exit(1);
			}
			glu_isBarrel_isLoaded = true;
		}
		return *glu_isBarrel_;
	}
	vector<int> &glu_det()
	{
		if (not glu_det_isLoaded) {
			if (glu_det_branch != 0) {
				glu_det_branch->GetEntry(index);
			} else { 
				printf("branch glu_det_branch does not exist!\n");
				exit(1);
			}
			glu_det_isLoaded = true;
		}
		return *glu_det_;
	}
	vector<int> &glu_lay()
	{
		if (not glu_lay_isLoaded) {
			if (glu_lay_branch != 0) {
				glu_lay_branch->GetEntry(index);
			} else { 
				printf("branch glu_lay_branch does not exist!\n");
				exit(1);
			}
			glu_lay_isLoaded = true;
		}
		return *glu_lay_;
	}
	vector<int> &glu_detId()
	{
		if (not glu_detId_isLoaded) {
			if (glu_detId_branch != 0) {
				glu_detId_branch->GetEntry(index);
			} else { 
				printf("branch glu_detId_branch does not exist!\n");
				exit(1);
			}
			glu_detId_isLoaded = true;
		}
		return *glu_detId_;
	}
	vector<int> &glu_monoIdx()
	{
		if (not glu_monoIdx_isLoaded) {
			if (glu_monoIdx_branch != 0) {
				glu_monoIdx_branch->GetEntry(index);
			} else { 
				printf("branch glu_monoIdx_branch does not exist!\n");
				exit(1);
			}
			glu_monoIdx_isLoaded = true;
		}
		return *glu_monoIdx_;
	}
	vector<int> &glu_stereoIdx()
	{
		if (not glu_stereoIdx_isLoaded) {
			if (glu_stereoIdx_branch != 0) {
				glu_stereoIdx_branch->GetEntry(index);
			} else { 
				printf("branch glu_stereoIdx_branch does not exist!\n");
				exit(1);
			}
			glu_stereoIdx_isLoaded = true;
		}
		return *glu_stereoIdx_;
	}
	vector<float> &glu_x()
	{
		if (not glu_x_isLoaded) {
			if (glu_x_branch != 0) {
				glu_x_branch->GetEntry(index);
			} else { 
				printf("branch glu_x_branch does not exist!\n");
				exit(1);
			}
			glu_x_isLoaded = true;
		}
		return *glu_x_;
	}
	vector<float> &glu_y()
	{
		if (not glu_y_isLoaded) {
			if (glu_y_branch != 0) {
				glu_y_branch->GetEntry(index);
			} else { 
				printf("branch glu_y_branch does not exist!\n");
				exit(1);
			}
			glu_y_isLoaded = true;
		}
		return *glu_y_;
	}
	vector<float> &glu_z()
	{
		if (not glu_z_isLoaded) {
			if (glu_z_branch != 0) {
				glu_z_branch->GetEntry(index);
			} else { 
				printf("branch glu_z_branch does not exist!\n");
				exit(1);
			}
			glu_z_isLoaded = true;
		}
		return *glu_z_;
	}
	vector<float> &glu_xx()
	{
		if (not glu_xx_isLoaded) {
			if (glu_xx_branch != 0) {
				glu_xx_branch->GetEntry(index);
			} else { 
				printf("branch glu_xx_branch does not exist!\n");
				exit(1);
			}
			glu_xx_isLoaded = true;
		}
		return *glu_xx_;
	}
	vector<float> &glu_xy()
	{
		if (not glu_xy_isLoaded) {
			if (glu_xy_branch != 0) {
				glu_xy_branch->GetEntry(index);
			} else { 
				printf("branch glu_xy_branch does not exist!\n");
				exit(1);
			}
			glu_xy_isLoaded = true;
		}
		return *glu_xy_;
	}
	vector<float> &glu_yy()
	{
		if (not glu_yy_isLoaded) {
			if (glu_yy_branch != 0) {
				glu_yy_branch->GetEntry(index);
			} else { 
				printf("branch glu_yy_branch does not exist!\n");
				exit(1);
			}
			glu_yy_isLoaded = true;
		}
		return *glu_yy_;
	}
	vector<float> &glu_yz()
	{
		if (not glu_yz_isLoaded) {
			if (glu_yz_branch != 0) {
				glu_yz_branch->GetEntry(index);
			} else { 
				printf("branch glu_yz_branch does not exist!\n");
				exit(1);
			}
			glu_yz_isLoaded = true;
		}
		return *glu_yz_;
	}
	vector<float> &glu_zz()
	{
		if (not glu_zz_isLoaded) {
			if (glu_zz_branch != 0) {
				glu_zz_branch->GetEntry(index);
			} else { 
				printf("branch glu_zz_branch does not exist!\n");
				exit(1);
			}
			glu_zz_isLoaded = true;
		}
		return *glu_zz_;
	}
	vector<float> &glu_zx()
	{
		if (not glu_zx_isLoaded) {
			if (glu_zx_branch != 0) {
				glu_zx_branch->GetEntry(index);
			} else { 
				printf("branch glu_zx_branch does not exist!\n");
				exit(1);
			}
			glu_zx_isLoaded = true;
		}
		return *glu_zx_;
	}
	vector<float> &glu_radL()
	{
		if (not glu_radL_isLoaded) {
			if (glu_radL_branch != 0) {
				glu_radL_branch->GetEntry(index);
			} else { 
				printf("branch glu_radL_branch does not exist!\n");
				exit(1);
			}
			glu_radL_isLoaded = true;
		}
		return *glu_radL_;
	}
	vector<float> &glu_bbxi()
	{
		if (not glu_bbxi_isLoaded) {
			if (glu_bbxi_branch != 0) {
				glu_bbxi_branch->GetEntry(index);
			} else { 
				printf("branch glu_bbxi_branch does not exist!\n");
				exit(1);
			}
			glu_bbxi_isLoaded = true;
		}
		return *glu_bbxi_;
	}
	float &bsp_x()
	{
		if (not bsp_x_isLoaded) {
			if (bsp_x_branch != 0) {
				bsp_x_branch->GetEntry(index);
			} else { 
				printf("branch bsp_x_branch does not exist!\n");
				exit(1);
			}
			bsp_x_isLoaded = true;
		}
		return bsp_x_;
	}
	float &bsp_y()
	{
		if (not bsp_y_isLoaded) {
			if (bsp_y_branch != 0) {
				bsp_y_branch->GetEntry(index);
			} else { 
				printf("branch bsp_y_branch does not exist!\n");
				exit(1);
			}
			bsp_y_isLoaded = true;
		}
		return bsp_y_;
	}
	float &bsp_z()
	{
		if (not bsp_z_isLoaded) {
			if (bsp_z_branch != 0) {
				bsp_z_branch->GetEntry(index);
			} else { 
				printf("branch bsp_z_branch does not exist!\n");
				exit(1);
			}
			bsp_z_isLoaded = true;
		}
		return bsp_z_;
	}
	float &bsp_sigmax()
	{
		if (not bsp_sigmax_isLoaded) {
			if (bsp_sigmax_branch != 0) {
				bsp_sigmax_branch->GetEntry(index);
			} else { 
				printf("branch bsp_sigmax_branch does not exist!\n");
				exit(1);
			}
			bsp_sigmax_isLoaded = true;
		}
		return bsp_sigmax_;
	}
	float &bsp_sigmay()
	{
		if (not bsp_sigmay_isLoaded) {
			if (bsp_sigmay_branch != 0) {
				bsp_sigmay_branch->GetEntry(index);
			} else { 
				printf("branch bsp_sigmay_branch does not exist!\n");
				exit(1);
			}
			bsp_sigmay_isLoaded = true;
		}
		return bsp_sigmay_;
	}
	float &bsp_sigmaz()
	{
		if (not bsp_sigmaz_isLoaded) {
			if (bsp_sigmaz_branch != 0) {
				bsp_sigmaz_branch->GetEntry(index);
			} else { 
				printf("branch bsp_sigmaz_branch does not exist!\n");
				exit(1);
			}
			bsp_sigmaz_isLoaded = true;
		}
		return bsp_sigmaz_;
	}
	vector<float> &see_px()
	{
		if (not see_px_isLoaded) {
			if (see_px_branch != 0) {
				see_px_branch->GetEntry(index);
			} else { 
				printf("branch see_px_branch does not exist!\n");
				exit(1);
			}
			see_px_isLoaded = true;
		}
		return *see_px_;
	}
	vector<float> &see_py()
	{
		if (not see_py_isLoaded) {
			if (see_py_branch != 0) {
				see_py_branch->GetEntry(index);
			} else { 
				printf("branch see_py_branch does not exist!\n");
				exit(1);
			}
			see_py_isLoaded = true;
		}
		return *see_py_;
	}
	vector<float> &see_pz()
	{
		if (not see_pz_isLoaded) {
			if (see_pz_branch != 0) {
				see_pz_branch->GetEntry(index);
			} else { 
				printf("branch see_pz_branch does not exist!\n");
				exit(1);
			}
			see_pz_isLoaded = true;
		}
		return *see_pz_;
	}
	vector<float> &see_pt()
	{
		if (not see_pt_isLoaded) {
			if (see_pt_branch != 0) {
				see_pt_branch->GetEntry(index);
			} else { 
				printf("branch see_pt_branch does not exist!\n");
				exit(1);
			}
			see_pt_isLoaded = true;
		}
		return *see_pt_;
	}
	vector<float> &see_eta()
	{
		if (not see_eta_isLoaded) {
			if (see_eta_branch != 0) {
				see_eta_branch->GetEntry(index);
			} else { 
				printf("branch see_eta_branch does not exist!\n");
				exit(1);
			}
			see_eta_isLoaded = true;
		}
		return *see_eta_;
	}
	vector<float> &see_phi()
	{
		if (not see_phi_isLoaded) {
			if (see_phi_branch != 0) {
				see_phi_branch->GetEntry(index);
			} else { 
				printf("branch see_phi_branch does not exist!\n");
				exit(1);
			}
			see_phi_isLoaded = true;
		}
		return *see_phi_;
	}
	vector<float> &see_dxy()
	{
		if (not see_dxy_isLoaded) {
			if (see_dxy_branch != 0) {
				see_dxy_branch->GetEntry(index);
			} else { 
				printf("branch see_dxy_branch does not exist!\n");
				exit(1);
			}
			see_dxy_isLoaded = true;
		}
		return *see_dxy_;
	}
	vector<float> &see_dz()
	{
		if (not see_dz_isLoaded) {
			if (see_dz_branch != 0) {
				see_dz_branch->GetEntry(index);
			} else { 
				printf("branch see_dz_branch does not exist!\n");
				exit(1);
			}
			see_dz_isLoaded = true;
		}
		return *see_dz_;
	}
	vector<float> &see_ptErr()
	{
		if (not see_ptErr_isLoaded) {
			if (see_ptErr_branch != 0) {
				see_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch see_ptErr_branch does not exist!\n");
				exit(1);
			}
			see_ptErr_isLoaded = true;
		}
		return *see_ptErr_;
	}
	vector<float> &see_etaErr()
	{
		if (not see_etaErr_isLoaded) {
			if (see_etaErr_branch != 0) {
				see_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch see_etaErr_branch does not exist!\n");
				exit(1);
			}
			see_etaErr_isLoaded = true;
		}
		return *see_etaErr_;
	}
	vector<float> &see_phiErr()
	{
		if (not see_phiErr_isLoaded) {
			if (see_phiErr_branch != 0) {
				see_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch see_phiErr_branch does not exist!\n");
				exit(1);
			}
			see_phiErr_isLoaded = true;
		}
		return *see_phiErr_;
	}
	vector<float> &see_dxyErr()
	{
		if (not see_dxyErr_isLoaded) {
			if (see_dxyErr_branch != 0) {
				see_dxyErr_branch->GetEntry(index);
			} else { 
				printf("branch see_dxyErr_branch does not exist!\n");
				exit(1);
			}
			see_dxyErr_isLoaded = true;
		}
		return *see_dxyErr_;
	}
	vector<float> &see_dzErr()
	{
		if (not see_dzErr_isLoaded) {
			if (see_dzErr_branch != 0) {
				see_dzErr_branch->GetEntry(index);
			} else { 
				printf("branch see_dzErr_branch does not exist!\n");
				exit(1);
			}
			see_dzErr_isLoaded = true;
		}
		return *see_dzErr_;
	}
	vector<float> &see_chi2()
	{
		if (not see_chi2_isLoaded) {
			if (see_chi2_branch != 0) {
				see_chi2_branch->GetEntry(index);
			} else { 
				printf("branch see_chi2_branch does not exist!\n");
				exit(1);
			}
			see_chi2_isLoaded = true;
		}
		return *see_chi2_;
	}
	vector<int> &see_q()
	{
		if (not see_q_isLoaded) {
			if (see_q_branch != 0) {
				see_q_branch->GetEntry(index);
			} else { 
				printf("branch see_q_branch does not exist!\n");
				exit(1);
			}
			see_q_isLoaded = true;
		}
		return *see_q_;
	}
	vector<int> &see_nValid()
	{
		if (not see_nValid_isLoaded) {
			if (see_nValid_branch != 0) {
				see_nValid_branch->GetEntry(index);
			} else { 
				printf("branch see_nValid_branch does not exist!\n");
				exit(1);
			}
			see_nValid_isLoaded = true;
		}
		return *see_nValid_;
	}
	vector<int> &see_nPixel()
	{
		if (not see_nPixel_isLoaded) {
			if (see_nPixel_branch != 0) {
				see_nPixel_branch->GetEntry(index);
			} else { 
				printf("branch see_nPixel_branch does not exist!\n");
				exit(1);
			}
			see_nPixel_isLoaded = true;
		}
		return *see_nPixel_;
	}
	vector<int> &see_nGlued()
	{
		if (not see_nGlued_isLoaded) {
			if (see_nGlued_branch != 0) {
				see_nGlued_branch->GetEntry(index);
			} else { 
				printf("branch see_nGlued_branch does not exist!\n");
				exit(1);
			}
			see_nGlued_isLoaded = true;
		}
		return *see_nGlued_;
	}
	vector<int> &see_nStrip()
	{
		if (not see_nStrip_isLoaded) {
			if (see_nStrip_branch != 0) {
				see_nStrip_branch->GetEntry(index);
			} else { 
				printf("branch see_nStrip_branch does not exist!\n");
				exit(1);
			}
			see_nStrip_isLoaded = true;
		}
		return *see_nStrip_;
	}
	vector<int> &see_algo()
	{
		if (not see_algo_isLoaded) {
			if (see_algo_branch != 0) {
				see_algo_branch->GetEntry(index);
			} else { 
				printf("branch see_algo_branch does not exist!\n");
				exit(1);
			}
			see_algo_isLoaded = true;
		}
		return *see_algo_;
	}
	vector<vector<int> > &see_pixelIdx()
	{
		if (not see_pixelIdx_isLoaded) {
			if (see_pixelIdx_branch != 0) {
				see_pixelIdx_branch->GetEntry(index);
			} else { 
				printf("branch see_pixelIdx_branch does not exist!\n");
				exit(1);
			}
			see_pixelIdx_isLoaded = true;
		}
		return *see_pixelIdx_;
	}
	vector<vector<int> > &see_gluedIdx()
	{
		if (not see_gluedIdx_isLoaded) {
			if (see_gluedIdx_branch != 0) {
				see_gluedIdx_branch->GetEntry(index);
			} else { 
				printf("branch see_gluedIdx_branch does not exist!\n");
				exit(1);
			}
			see_gluedIdx_isLoaded = true;
		}
		return *see_gluedIdx_;
	}
	vector<vector<int> > &see_stripIdx()
	{
		if (not see_stripIdx_isLoaded) {
			if (see_stripIdx_branch != 0) {
				see_stripIdx_branch->GetEntry(index);
			} else { 
				printf("branch see_stripIdx_branch does not exist!\n");
				exit(1);
			}
			see_stripIdx_isLoaded = true;
		}
		return *see_stripIdx_;
	}
	vector<int> &algo_offset()
	{
		if (not algo_offset_isLoaded) {
			if (algo_offset_branch != 0) {
				algo_offset_branch->GetEntry(index);
			} else { 
				printf("branch algo_offset_branch does not exist!\n");
				exit(1);
			}
			algo_offset_isLoaded = true;
		}
		return *algo_offset_;
	}
};

#ifndef __CINT__
extern tkph2 cms2;
#endif

namespace tas {
	vector<float> &trk_px();
	vector<float> &trk_py();
	vector<float> &trk_pz();
	vector<float> &trk_pt();
	vector<float> &trk_eta();
	vector<float> &trk_phi();
	vector<float> &trk_dxy();
	vector<float> &trk_dz();
	vector<float> &trk_ptErr();
	vector<float> &trk_etaErr();
	vector<float> &trk_phiErr();
	vector<float> &trk_dxyErr();
	vector<float> &trk_dzErr();
	vector<float> &trk_nChi2();
	vector<float> &trk_shareFrac();
	vector<int> &trk_q();
	vector<int> &trk_nValid();
	vector<int> &trk_nInvalid();
	vector<int> &trk_nPixel();
	vector<int> &trk_nStrip();
	vector<int> &trk_n3DLay();
	vector<int> &trk_algo();
	vector<int> &trk_isHP();
	vector<int> &trk_seedIdx();
	vector<int> &trk_simIdx();
	vector<vector<int> > &trk_pixelIdx();
	vector<vector<int> > &trk_stripIdx();
	vector<float> &sim_px();
	vector<float> &sim_py();
	vector<float> &sim_pz();
	vector<float> &sim_pt();
	vector<float> &sim_eta();
	vector<float> &sim_phi();
	vector<float> &sim_dxy();
	vector<float> &sim_dz();
	vector<float> &sim_prodx();
	vector<float> &sim_prody();
	vector<float> &sim_prodz();
	vector<float> &sim_shareFrac();
	vector<int> &sim_q();
	vector<int> &sim_nValid();
	vector<int> &sim_nPixel();
	vector<int> &sim_nStrip();
	vector<int> &sim_n3DLay();
	vector<int> &sim_trkIdx();
	vector<vector<int> > &sim_pixelIdx();
	vector<vector<int> > &sim_stripIdx();
	vector<int> &pix_isBarrel();
	vector<int> &pix_lay();
	vector<int> &pix_detId();
	vector<int> &pix_nSimTrk();
	vector<int> &pix_simTrkIdx();
	vector<int> &pix_particle();
	vector<int> &pix_process();
	vector<int> &pix_bunchXing();
	vector<int> &pix_event();
	vector<float> &pix_x();
	vector<float> &pix_y();
	vector<float> &pix_z();
	vector<float> &pix_xx();
	vector<float> &pix_xy();
	vector<float> &pix_yy();
	vector<float> &pix_yz();
	vector<float> &pix_zz();
	vector<float> &pix_zx();
	vector<float> &pix_xsim();
	vector<float> &pix_ysim();
	vector<float> &pix_zsim();
	vector<float> &pix_pxsim();
	vector<float> &pix_pysim();
	vector<float> &pix_pzsim();
	vector<float> &pix_pathprop();
	vector<float> &pix_xsimprop();
	vector<float> &pix_ysimprop();
	vector<float> &pix_zsimprop();
	vector<float> &pix_pxsimprop();
	vector<float> &pix_pysimprop();
	vector<float> &pix_pzsimprop();
	vector<float> &pix_eloss();
	vector<float> &pix_radL();
	vector<float> &pix_bbxi();
	vector<int> &str_isBarrel();
	vector<int> &str_isStereo();
	vector<int> &str_det();
	vector<int> &str_lay();
	vector<int> &str_detId();
	vector<int> &str_nSimTrk();
	vector<int> &str_simTrkIdx();
	vector<int> &str_particle();
	vector<int> &str_process();
	vector<int> &str_bunchXing();
	vector<int> &str_event();
	vector<float> &str_x();
	vector<float> &str_y();
	vector<float> &str_z();
	vector<float> &str_xx();
	vector<float> &str_xy();
	vector<float> &str_yy();
	vector<float> &str_yz();
	vector<float> &str_zz();
	vector<float> &str_zx();
	vector<float> &str_xsim();
	vector<float> &str_ysim();
	vector<float> &str_zsim();
	vector<float> &str_pxsim();
	vector<float> &str_pysim();
	vector<float> &str_pzsim();
	vector<float> &str_eloss();
	vector<float> &str_radL();
	vector<float> &str_bbxi();
	vector<int> &glu_isBarrel();
	vector<int> &glu_det();
	vector<int> &glu_lay();
	vector<int> &glu_detId();
	vector<int> &glu_monoIdx();
	vector<int> &glu_stereoIdx();
	vector<float> &glu_x();
	vector<float> &glu_y();
	vector<float> &glu_z();
	vector<float> &glu_xx();
	vector<float> &glu_xy();
	vector<float> &glu_yy();
	vector<float> &glu_yz();
	vector<float> &glu_zz();
	vector<float> &glu_zx();
	vector<float> &glu_radL();
	vector<float> &glu_bbxi();
	float &bsp_x();
	float &bsp_y();
	float &bsp_z();
	float &bsp_sigmax();
	float &bsp_sigmay();
	float &bsp_sigmaz();
	vector<float> &see_px();
	vector<float> &see_py();
	vector<float> &see_pz();
	vector<float> &see_pt();
	vector<float> &see_eta();
	vector<float> &see_phi();
	vector<float> &see_dxy();
	vector<float> &see_dz();
	vector<float> &see_ptErr();
	vector<float> &see_etaErr();
	vector<float> &see_phiErr();
	vector<float> &see_dxyErr();
	vector<float> &see_dzErr();
	vector<float> &see_chi2();
	vector<int> &see_q();
	vector<int> &see_nValid();
	vector<int> &see_nPixel();
	vector<int> &see_nGlued();
	vector<int> &see_nStrip();
	vector<int> &see_algo();
	vector<vector<int> > &see_pixelIdx();
	vector<vector<int> > &see_gluedIdx();
	vector<vector<int> > &see_stripIdx();
	vector<int> &algo_offset();
}
#endif
