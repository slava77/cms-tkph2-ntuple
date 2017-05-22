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
	unsigned long long	event_;
	TBranch *event_branch;
	bool event_isLoaded;
	unsigned int	lumi_;
	TBranch *lumi_branch;
	bool lumi_isLoaded;
	unsigned int	run_;
	TBranch *run_branch;
	bool run_isLoaded;
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
	vector<float> *trk_inner_px_;
	TBranch *trk_inner_px_branch;
	bool trk_inner_px_isLoaded;
	vector<float> *trk_inner_py_;
	TBranch *trk_inner_py_branch;
	bool trk_inner_py_isLoaded;
	vector<float> *trk_inner_pz_;
	TBranch *trk_inner_pz_branch;
	bool trk_inner_pz_isLoaded;
	vector<float> *trk_inner_pt_;
	TBranch *trk_inner_pt_branch;
	bool trk_inner_pt_isLoaded;
	vector<float> *trk_outer_px_;
	TBranch *trk_outer_px_branch;
	bool trk_outer_px_isLoaded;
	vector<float> *trk_outer_py_;
	TBranch *trk_outer_py_branch;
	bool trk_outer_py_isLoaded;
	vector<float> *trk_outer_pz_;
	TBranch *trk_outer_pz_branch;
	bool trk_outer_pz_isLoaded;
	vector<float> *trk_outer_pt_;
	TBranch *trk_outer_pt_branch;
	bool trk_outer_pt_isLoaded;
	vector<float> *trk_eta_;
	TBranch *trk_eta_branch;
	bool trk_eta_isLoaded;
	vector<float> *trk_lambda_;
	TBranch *trk_lambda_branch;
	bool trk_lambda_isLoaded;
	vector<float> *trk_cotTheta_;
	TBranch *trk_cotTheta_branch;
	bool trk_cotTheta_isLoaded;
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
	vector<float> *trk_lambdaErr_;
	TBranch *trk_lambdaErr_branch;
	bool trk_lambdaErr_isLoaded;
	vector<float> *trk_phiErr_;
	TBranch *trk_phiErr_branch;
	bool trk_phiErr_isLoaded;
	vector<float> *trk_dxyErr_;
	TBranch *trk_dxyErr_branch;
	bool trk_dxyErr_isLoaded;
	vector<float> *trk_dzErr_;
	TBranch *trk_dzErr_branch;
	bool trk_dzErr_isLoaded;
	vector<float> *trk_refpoint_x_;
	TBranch *trk_refpoint_x_branch;
	bool trk_refpoint_x_isLoaded;
	vector<float> *trk_refpoint_y_;
	TBranch *trk_refpoint_y_branch;
	bool trk_refpoint_y_isLoaded;
	vector<float> *trk_refpoint_z_;
	TBranch *trk_refpoint_z_branch;
	bool trk_refpoint_z_isLoaded;
	vector<float> *trk_nChi2_;
	TBranch *trk_nChi2_branch;
	bool trk_nChi2_isLoaded;
	vector<int> *trk_q_;
	TBranch *trk_q_branch;
	bool trk_q_isLoaded;
	vector<unsigned int> *trk_nValid_;
	TBranch *trk_nValid_branch;
	bool trk_nValid_isLoaded;
	vector<unsigned int> *trk_nInvalid_;
	TBranch *trk_nInvalid_branch;
	bool trk_nInvalid_isLoaded;
	vector<unsigned int> *trk_nPixel_;
	TBranch *trk_nPixel_branch;
	bool trk_nPixel_isLoaded;
	vector<unsigned int> *trk_nStrip_;
	TBranch *trk_nStrip_branch;
	bool trk_nStrip_isLoaded;
	vector<unsigned int> *trk_nPixelLay_;
	TBranch *trk_nPixelLay_branch;
	bool trk_nPixelLay_isLoaded;
	vector<unsigned int> *trk_nStripLay_;
	TBranch *trk_nStripLay_branch;
	bool trk_nStripLay_isLoaded;
	vector<unsigned int> *trk_n3DLay_;
	TBranch *trk_n3DLay_branch;
	bool trk_n3DLay_isLoaded;
	vector<unsigned int> *trk_nOuterLost_;
	TBranch *trk_nOuterLost_branch;
	bool trk_nOuterLost_isLoaded;
	vector<unsigned int> *trk_nInnerLost_;
	TBranch *trk_nInnerLost_branch;
	bool trk_nInnerLost_isLoaded;
	vector<unsigned int> *trk_algo_;
	TBranch *trk_algo_branch;
	bool trk_algo_isLoaded;
	vector<unsigned int> *trk_originalAlgo_;
	TBranch *trk_originalAlgo_branch;
	bool trk_originalAlgo_isLoaded;
	vector<ULong64_t> *trk_algoMask_;
	TBranch *trk_algoMask_branch;
	bool trk_algoMask_isLoaded;
	vector<unsigned short> *trk_stopReason_;
	TBranch *trk_stopReason_branch;
	bool trk_stopReason_isLoaded;
	vector<short> *trk_isHP_;
	TBranch *trk_isHP_branch;
	bool trk_isHP_isLoaded;
	vector<int> *trk_seedIdx_;
	TBranch *trk_seedIdx_branch;
	bool trk_seedIdx_isLoaded;
	vector<int> *trk_vtxIdx_;
	TBranch *trk_vtxIdx_branch;
	bool trk_vtxIdx_isLoaded;
	vector<vector<float> > *trk_shareFrac_;
	TBranch *trk_shareFrac_branch;
	bool trk_shareFrac_isLoaded;
	vector<vector<int> > *trk_simTrkIdx_;
	TBranch *trk_simTrkIdx_branch;
	bool trk_simTrkIdx_isLoaded;
	vector<vector<int> > *trk_hitIdx_;
	TBranch *trk_hitIdx_branch;
	bool trk_hitIdx_isLoaded;
	vector<vector<int> > *trk_hitType_;
	TBranch *trk_hitType_branch;
	bool trk_hitType_isLoaded;
	vector<int> *sim_event_;
	TBranch *sim_event_branch;
	bool sim_event_isLoaded;
	vector<int> *sim_bunchCrossing_;
	TBranch *sim_bunchCrossing_branch;
	bool sim_bunchCrossing_isLoaded;
	vector<int> *sim_pdgId_;
	TBranch *sim_pdgId_branch;
	bool sim_pdgId_isLoaded;
	vector<vector<int> > *sim_genPdgIds_;
	TBranch *sim_genPdgIds_branch;
	bool sim_genPdgIds_isLoaded;
	vector<int> *sim_isFromBHadron_;
	TBranch *sim_isFromBHadron_branch;
	bool sim_isFromBHadron_isLoaded;
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
	vector<float> *sim_pca_pt_;
	TBranch *sim_pca_pt_branch;
	bool sim_pca_pt_isLoaded;
	vector<float> *sim_pca_eta_;
	TBranch *sim_pca_eta_branch;
	bool sim_pca_eta_isLoaded;
	vector<float> *sim_pca_lambda_;
	TBranch *sim_pca_lambda_branch;
	bool sim_pca_lambda_isLoaded;
	vector<float> *sim_pca_cotTheta_;
	TBranch *sim_pca_cotTheta_branch;
	bool sim_pca_cotTheta_isLoaded;
	vector<float> *sim_pca_phi_;
	TBranch *sim_pca_phi_branch;
	bool sim_pca_phi_isLoaded;
	vector<float> *sim_pca_dxy_;
	TBranch *sim_pca_dxy_branch;
	bool sim_pca_dxy_isLoaded;
	vector<float> *sim_pca_dz_;
	TBranch *sim_pca_dz_branch;
	bool sim_pca_dz_isLoaded;
	vector<int> *sim_q_;
	TBranch *sim_q_branch;
	bool sim_q_isLoaded;
	vector<unsigned int> *sim_nValid_;
	TBranch *sim_nValid_branch;
	bool sim_nValid_isLoaded;
	vector<unsigned int> *sim_nPixel_;
	TBranch *sim_nPixel_branch;
	bool sim_nPixel_isLoaded;
	vector<unsigned int> *sim_nStrip_;
	TBranch *sim_nStrip_branch;
	bool sim_nStrip_isLoaded;
	vector<unsigned int> *sim_nLay_;
	TBranch *sim_nLay_branch;
	bool sim_nLay_isLoaded;
	vector<unsigned int> *sim_nPixelLay_;
	TBranch *sim_nPixelLay_branch;
	bool sim_nPixelLay_isLoaded;
	vector<unsigned int> *sim_n3DLay_;
	TBranch *sim_n3DLay_branch;
	bool sim_n3DLay_isLoaded;
	vector<vector<int> > *sim_trkIdx_;
	TBranch *sim_trkIdx_branch;
	bool sim_trkIdx_isLoaded;
	vector<vector<float> > *sim_shareFrac_;
	TBranch *sim_shareFrac_branch;
	bool sim_shareFrac_isLoaded;
	vector<vector<int> > *sim_seedIdx_;
	TBranch *sim_seedIdx_branch;
	bool sim_seedIdx_isLoaded;
	vector<int> *sim_parentVtxIdx_;
	TBranch *sim_parentVtxIdx_branch;
	bool sim_parentVtxIdx_isLoaded;
	vector<vector<int> > *sim_decayVtxIdx_;
	TBranch *sim_decayVtxIdx_branch;
	bool sim_decayVtxIdx_isLoaded;
	vector<vector<int> > *sim_simHitIdx_;
	TBranch *sim_simHitIdx_branch;
	bool sim_simHitIdx_isLoaded;
	vector<short> *pix_isBarrel_;
	TBranch *pix_isBarrel_branch;
	bool pix_isBarrel_isLoaded;
	vector<unsigned short> *pix_det_;
	TBranch *pix_det_branch;
	bool pix_det_isLoaded;
	vector<unsigned short> *pix_lay_;
	TBranch *pix_lay_branch;
	bool pix_lay_isLoaded;
	vector<unsigned int> *pix_detId_;
	TBranch *pix_detId_branch;
	bool pix_detId_isLoaded;
	vector<vector<int> > *pix_trkIdx_;
	TBranch *pix_trkIdx_branch;
	bool pix_trkIdx_isLoaded;
	vector<vector<int> > *pix_seeIdx_;
	TBranch *pix_seeIdx_branch;
	bool pix_seeIdx_isLoaded;
	vector<vector<int> > *pix_simHitIdx_;
	TBranch *pix_simHitIdx_branch;
	bool pix_simHitIdx_isLoaded;
	vector<vector<float> > *pix_xySignificance_;
	TBranch *pix_xySignificance_branch;
	bool pix_xySignificance_isLoaded;
	vector<vector<float> > *pix_chargeFraction_;
	TBranch *pix_chargeFraction_branch;
	bool pix_chargeFraction_isLoaded;
	vector<unsigned short> *pix_simType_;
	TBranch *pix_simType_branch;
	bool pix_simType_isLoaded;
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
	vector<float> *pix_radL_;
	TBranch *pix_radL_branch;
	bool pix_radL_isLoaded;
	vector<float> *pix_bbxi_;
	TBranch *pix_bbxi_branch;
	bool pix_bbxi_isLoaded;
	vector<short> *ph2_isBarrel_;
	TBranch *ph2_isBarrel_branch;
	bool ph2_isBarrel_isLoaded;
	vector<unsigned short> *ph2_det_;
	TBranch *ph2_det_branch;
	bool ph2_det_isLoaded;
	vector<unsigned short> *ph2_lay_;
	TBranch *ph2_lay_branch;
	bool ph2_lay_isLoaded;
	vector<unsigned int> *ph2_detId_;
	TBranch *ph2_detId_branch;
	bool ph2_detId_isLoaded;
	vector<vector<int> > *ph2_trkIdx_;
	TBranch *ph2_trkIdx_branch;
	bool ph2_trkIdx_isLoaded;
	vector<vector<int> > *ph2_seeIdx_;
	TBranch *ph2_seeIdx_branch;
	bool ph2_seeIdx_isLoaded;
	vector<vector<int> > *ph2_simHitIdx_;
	TBranch *ph2_simHitIdx_branch;
	bool ph2_simHitIdx_isLoaded;
	vector<vector<float> > *ph2_xySignificance_;
	TBranch *ph2_xySignificance_branch;
	bool ph2_xySignificance_isLoaded;
	vector<unsigned short> *ph2_simType_;
	TBranch *ph2_simType_branch;
	bool ph2_simType_isLoaded;
	vector<float> *ph2_x_;
	TBranch *ph2_x_branch;
	bool ph2_x_isLoaded;
	vector<float> *ph2_y_;
	TBranch *ph2_y_branch;
	bool ph2_y_isLoaded;
	vector<float> *ph2_z_;
	TBranch *ph2_z_branch;
	bool ph2_z_isLoaded;
	vector<float> *ph2_xx_;
	TBranch *ph2_xx_branch;
	bool ph2_xx_isLoaded;
	vector<float> *ph2_xy_;
	TBranch *ph2_xy_branch;
	bool ph2_xy_isLoaded;
	vector<float> *ph2_yy_;
	TBranch *ph2_yy_branch;
	bool ph2_yy_isLoaded;
	vector<float> *ph2_yz_;
	TBranch *ph2_yz_branch;
	bool ph2_yz_isLoaded;
	vector<float> *ph2_zz_;
	TBranch *ph2_zz_branch;
	bool ph2_zz_isLoaded;
	vector<float> *ph2_zx_;
	TBranch *ph2_zx_branch;
	bool ph2_zx_isLoaded;
	vector<float> *ph2_radL_;
	TBranch *ph2_radL_branch;
	bool ph2_radL_isLoaded;
	vector<float> *ph2_bbxi_;
	TBranch *ph2_bbxi_branch;
	bool ph2_bbxi_isLoaded;
	vector<short> *inv_isBarrel_;
	TBranch *inv_isBarrel_branch;
	bool inv_isBarrel_isLoaded;
	vector<unsigned short> *inv_det_;
	TBranch *inv_det_branch;
	bool inv_det_isLoaded;
	vector<unsigned short> *inv_lay_;
	TBranch *inv_lay_branch;
	bool inv_lay_isLoaded;
	vector<unsigned int> *inv_detId_;
	TBranch *inv_detId_branch;
	bool inv_detId_isLoaded;
	vector<unsigned short> *inv_type_;
	TBranch *inv_type_branch;
	bool inv_type_isLoaded;
	vector<unsigned short> *simhit_det_;
	TBranch *simhit_det_branch;
	bool simhit_det_isLoaded;
	vector<unsigned short> *simhit_lay_;
	TBranch *simhit_lay_branch;
	bool simhit_lay_isLoaded;
	vector<unsigned int> *simhit_detId_;
	TBranch *simhit_detId_branch;
	bool simhit_detId_isLoaded;
	vector<float> *simhit_x_;
	TBranch *simhit_x_branch;
	bool simhit_x_isLoaded;
	vector<float> *simhit_y_;
	TBranch *simhit_y_branch;
	bool simhit_y_isLoaded;
	vector<float> *simhit_z_;
	TBranch *simhit_z_branch;
	bool simhit_z_isLoaded;
	vector<float> *simhit_px_;
	TBranch *simhit_px_branch;
	bool simhit_px_isLoaded;
	vector<float> *simhit_py_;
	TBranch *simhit_py_branch;
	bool simhit_py_isLoaded;
	vector<float> *simhit_pz_;
	TBranch *simhit_pz_branch;
	bool simhit_pz_isLoaded;
	vector<int> *simhit_particle_;
	TBranch *simhit_particle_branch;
	bool simhit_particle_isLoaded;
	vector<short> *simhit_process_;
	TBranch *simhit_process_branch;
	bool simhit_process_isLoaded;
	vector<float> *simhit_eloss_;
	TBranch *simhit_eloss_branch;
	bool simhit_eloss_isLoaded;
	vector<float> *simhit_tof_;
	TBranch *simhit_tof_branch;
	bool simhit_tof_isLoaded;
	vector<int> *simhit_simTrkIdx_;
	TBranch *simhit_simTrkIdx_branch;
	bool simhit_simTrkIdx_isLoaded;
	vector<vector<int> > *simhit_hitIdx_;
	TBranch *simhit_hitIdx_branch;
	bool simhit_hitIdx_isLoaded;
	vector<vector<int> > *simhit_hitType_;
	TBranch *simhit_hitType_branch;
	bool simhit_hitType_isLoaded;
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
	vector<short> *see_fitok_;
	TBranch *see_fitok_branch;
	bool see_fitok_isLoaded;
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
	vector<float> *see_statePt_;
	TBranch *see_statePt_branch;
	bool see_statePt_isLoaded;
	vector<float> *see_stateTrajX_;
	TBranch *see_stateTrajX_branch;
	bool see_stateTrajX_isLoaded;
	vector<float> *see_stateTrajY_;
	TBranch *see_stateTrajY_branch;
	bool see_stateTrajY_isLoaded;
	vector<float> *see_stateTrajPx_;
	TBranch *see_stateTrajPx_branch;
	bool see_stateTrajPx_isLoaded;
	vector<float> *see_stateTrajPy_;
	TBranch *see_stateTrajPy_branch;
	bool see_stateTrajPy_isLoaded;
	vector<float> *see_stateTrajPz_;
	TBranch *see_stateTrajPz_branch;
	bool see_stateTrajPz_isLoaded;
	vector<float> *see_stateTrajGlbX_;
	TBranch *see_stateTrajGlbX_branch;
	bool see_stateTrajGlbX_isLoaded;
	vector<float> *see_stateTrajGlbY_;
	TBranch *see_stateTrajGlbY_branch;
	bool see_stateTrajGlbY_isLoaded;
	vector<float> *see_stateTrajGlbZ_;
	TBranch *see_stateTrajGlbZ_branch;
	bool see_stateTrajGlbZ_isLoaded;
	vector<float> *see_stateTrajGlbPx_;
	TBranch *see_stateTrajGlbPx_branch;
	bool see_stateTrajGlbPx_isLoaded;
	vector<float> *see_stateTrajGlbPy_;
	TBranch *see_stateTrajGlbPy_branch;
	bool see_stateTrajGlbPy_isLoaded;
	vector<float> *see_stateTrajGlbPz_;
	TBranch *see_stateTrajGlbPz_branch;
	bool see_stateTrajGlbPz_isLoaded;
	vector<float> *see_stateCcov00_;
	TBranch *see_stateCcov00_branch;
	bool see_stateCcov00_isLoaded;
	vector<float> *see_stateCcov01_;
	TBranch *see_stateCcov01_branch;
	bool see_stateCcov01_isLoaded;
	vector<float> *see_stateCcov02_;
	TBranch *see_stateCcov02_branch;
	bool see_stateCcov02_isLoaded;
	vector<float> *see_stateCcov03_;
	TBranch *see_stateCcov03_branch;
	bool see_stateCcov03_isLoaded;
	vector<float> *see_stateCcov04_;
	TBranch *see_stateCcov04_branch;
	bool see_stateCcov04_isLoaded;
	vector<float> *see_stateCcov05_;
	TBranch *see_stateCcov05_branch;
	bool see_stateCcov05_isLoaded;
	vector<float> *see_stateCcov11_;
	TBranch *see_stateCcov11_branch;
	bool see_stateCcov11_isLoaded;
	vector<float> *see_stateCcov12_;
	TBranch *see_stateCcov12_branch;
	bool see_stateCcov12_isLoaded;
	vector<float> *see_stateCcov13_;
	TBranch *see_stateCcov13_branch;
	bool see_stateCcov13_isLoaded;
	vector<float> *see_stateCcov14_;
	TBranch *see_stateCcov14_branch;
	bool see_stateCcov14_isLoaded;
	vector<float> *see_stateCcov15_;
	TBranch *see_stateCcov15_branch;
	bool see_stateCcov15_isLoaded;
	vector<float> *see_stateCcov22_;
	TBranch *see_stateCcov22_branch;
	bool see_stateCcov22_isLoaded;
	vector<float> *see_stateCcov23_;
	TBranch *see_stateCcov23_branch;
	bool see_stateCcov23_isLoaded;
	vector<float> *see_stateCcov24_;
	TBranch *see_stateCcov24_branch;
	bool see_stateCcov24_isLoaded;
	vector<float> *see_stateCcov25_;
	TBranch *see_stateCcov25_branch;
	bool see_stateCcov25_isLoaded;
	vector<float> *see_stateCcov33_;
	TBranch *see_stateCcov33_branch;
	bool see_stateCcov33_isLoaded;
	vector<float> *see_stateCcov34_;
	TBranch *see_stateCcov34_branch;
	bool see_stateCcov34_isLoaded;
	vector<float> *see_stateCcov35_;
	TBranch *see_stateCcov35_branch;
	bool see_stateCcov35_isLoaded;
	vector<float> *see_stateCcov44_;
	TBranch *see_stateCcov44_branch;
	bool see_stateCcov44_isLoaded;
	vector<float> *see_stateCcov45_;
	TBranch *see_stateCcov45_branch;
	bool see_stateCcov45_isLoaded;
	vector<float> *see_stateCcov55_;
	TBranch *see_stateCcov55_branch;
	bool see_stateCcov55_isLoaded;
	vector<int> *see_q_;
	TBranch *see_q_branch;
	bool see_q_isLoaded;
	vector<unsigned int> *see_nValid_;
	TBranch *see_nValid_branch;
	bool see_nValid_isLoaded;
	vector<unsigned int> *see_nPixel_;
	TBranch *see_nPixel_branch;
	bool see_nPixel_isLoaded;
	vector<unsigned int> *see_nGlued_;
	TBranch *see_nGlued_branch;
	bool see_nGlued_isLoaded;
	vector<unsigned int> *see_nStrip_;
	TBranch *see_nStrip_branch;
	bool see_nStrip_isLoaded;
	vector<unsigned int> *see_nPhase2OT_;
	TBranch *see_nPhase2OT_branch;
	bool see_nPhase2OT_isLoaded;
	vector<unsigned int> *see_algo_;
	TBranch *see_algo_branch;
	bool see_algo_isLoaded;
	vector<unsigned short> *see_stopReason_;
	TBranch *see_stopReason_branch;
	bool see_stopReason_isLoaded;
	vector<int> *see_trkIdx_;
	TBranch *see_trkIdx_branch;
	bool see_trkIdx_isLoaded;
	vector<vector<float> > *see_shareFrac_;
	TBranch *see_shareFrac_branch;
	bool see_shareFrac_isLoaded;
	vector<vector<int> > *see_simTrkIdx_;
	TBranch *see_simTrkIdx_branch;
	bool see_simTrkIdx_isLoaded;
	vector<vector<int> > *see_hitIdx_;
	TBranch *see_hitIdx_branch;
	bool see_hitIdx_isLoaded;
	vector<vector<int> > *see_hitType_;
	TBranch *see_hitType_branch;
	bool see_hitType_isLoaded;
	vector<unsigned int> *see_offset_;
	TBranch *see_offset_branch;
	bool see_offset_isLoaded;
	vector<float> *vtx_x_;
	TBranch *vtx_x_branch;
	bool vtx_x_isLoaded;
	vector<float> *vtx_y_;
	TBranch *vtx_y_branch;
	bool vtx_y_isLoaded;
	vector<float> *vtx_z_;
	TBranch *vtx_z_branch;
	bool vtx_z_isLoaded;
	vector<float> *vtx_xErr_;
	TBranch *vtx_xErr_branch;
	bool vtx_xErr_isLoaded;
	vector<float> *vtx_yErr_;
	TBranch *vtx_yErr_branch;
	bool vtx_yErr_isLoaded;
	vector<float> *vtx_zErr_;
	TBranch *vtx_zErr_branch;
	bool vtx_zErr_isLoaded;
	vector<float> *vtx_ndof_;
	TBranch *vtx_ndof_branch;
	bool vtx_ndof_isLoaded;
	vector<float> *vtx_chi2_;
	TBranch *vtx_chi2_branch;
	bool vtx_chi2_isLoaded;
	vector<short> *vtx_fake_;
	TBranch *vtx_fake_branch;
	bool vtx_fake_isLoaded;
	vector<short> *vtx_valid_;
	TBranch *vtx_valid_branch;
	bool vtx_valid_isLoaded;
	vector<vector<int> > *vtx_trkIdx_;
	TBranch *vtx_trkIdx_branch;
	bool vtx_trkIdx_isLoaded;
	vector<int> *simvtx_event_;
	TBranch *simvtx_event_branch;
	bool simvtx_event_isLoaded;
	vector<int> *simvtx_bunchCrossing_;
	TBranch *simvtx_bunchCrossing_branch;
	bool simvtx_bunchCrossing_isLoaded;
	vector<unsigned int> *simvtx_processType_;
	TBranch *simvtx_processType_branch;
	bool simvtx_processType_isLoaded;
	vector<float> *simvtx_x_;
	TBranch *simvtx_x_branch;
	bool simvtx_x_isLoaded;
	vector<float> *simvtx_y_;
	TBranch *simvtx_y_branch;
	bool simvtx_y_isLoaded;
	vector<float> *simvtx_z_;
	TBranch *simvtx_z_branch;
	bool simvtx_z_isLoaded;
	vector<vector<int> > *simvtx_sourceSimIdx_;
	TBranch *simvtx_sourceSimIdx_branch;
	bool simvtx_sourceSimIdx_isLoaded;
	vector<vector<int> > *simvtx_daughterSimIdx_;
	TBranch *simvtx_daughterSimIdx_branch;
	bool simvtx_daughterSimIdx_isLoaded;
	vector<int> *simpv_idx_;
	TBranch *simpv_idx_branch;
	bool simpv_idx_isLoaded;
public: 
void Init(TTree *tree) {
  tree->SetMakeClass(1);
	event_branch = 0;
	if (tree->GetBranch("event") != 0) {
		event_branch = tree->GetBranch("event");
		event_branch->SetAddress(&event_);
	}
	if(event_branch == 0 ) {
	cout << "Branch event does not exist." << endl;
	}
	lumi_branch = 0;
	if (tree->GetBranch("lumi") != 0) {
		lumi_branch = tree->GetBranch("lumi");
		lumi_branch->SetAddress(&lumi_);
	}
	if(lumi_branch == 0 ) {
	cout << "Branch lumi does not exist." << endl;
	}
	run_branch = 0;
	if (tree->GetBranch("run") != 0) {
		run_branch = tree->GetBranch("run");
		run_branch->SetAddress(&run_);
	}
	if(run_branch == 0 ) {
	cout << "Branch run does not exist." << endl;
	}
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
	trk_inner_px_branch = 0;
	if (tree->GetBranch("trk_inner_px") != 0) {
		trk_inner_px_branch = tree->GetBranch("trk_inner_px");
		trk_inner_px_branch->SetAddress(&trk_inner_px_);
	}
	if(trk_inner_px_branch == 0 ) {
	cout << "Branch trk_inner_px does not exist." << endl;
	}
	trk_inner_py_branch = 0;
	if (tree->GetBranch("trk_inner_py") != 0) {
		trk_inner_py_branch = tree->GetBranch("trk_inner_py");
		trk_inner_py_branch->SetAddress(&trk_inner_py_);
	}
	if(trk_inner_py_branch == 0 ) {
	cout << "Branch trk_inner_py does not exist." << endl;
	}
	trk_inner_pz_branch = 0;
	if (tree->GetBranch("trk_inner_pz") != 0) {
		trk_inner_pz_branch = tree->GetBranch("trk_inner_pz");
		trk_inner_pz_branch->SetAddress(&trk_inner_pz_);
	}
	if(trk_inner_pz_branch == 0 ) {
	cout << "Branch trk_inner_pz does not exist." << endl;
	}
	trk_inner_pt_branch = 0;
	if (tree->GetBranch("trk_inner_pt") != 0) {
		trk_inner_pt_branch = tree->GetBranch("trk_inner_pt");
		trk_inner_pt_branch->SetAddress(&trk_inner_pt_);
	}
	if(trk_inner_pt_branch == 0 ) {
	cout << "Branch trk_inner_pt does not exist." << endl;
	}
	trk_outer_px_branch = 0;
	if (tree->GetBranch("trk_outer_px") != 0) {
		trk_outer_px_branch = tree->GetBranch("trk_outer_px");
		trk_outer_px_branch->SetAddress(&trk_outer_px_);
	}
	if(trk_outer_px_branch == 0 ) {
	cout << "Branch trk_outer_px does not exist." << endl;
	}
	trk_outer_py_branch = 0;
	if (tree->GetBranch("trk_outer_py") != 0) {
		trk_outer_py_branch = tree->GetBranch("trk_outer_py");
		trk_outer_py_branch->SetAddress(&trk_outer_py_);
	}
	if(trk_outer_py_branch == 0 ) {
	cout << "Branch trk_outer_py does not exist." << endl;
	}
	trk_outer_pz_branch = 0;
	if (tree->GetBranch("trk_outer_pz") != 0) {
		trk_outer_pz_branch = tree->GetBranch("trk_outer_pz");
		trk_outer_pz_branch->SetAddress(&trk_outer_pz_);
	}
	if(trk_outer_pz_branch == 0 ) {
	cout << "Branch trk_outer_pz does not exist." << endl;
	}
	trk_outer_pt_branch = 0;
	if (tree->GetBranch("trk_outer_pt") != 0) {
		trk_outer_pt_branch = tree->GetBranch("trk_outer_pt");
		trk_outer_pt_branch->SetAddress(&trk_outer_pt_);
	}
	if(trk_outer_pt_branch == 0 ) {
	cout << "Branch trk_outer_pt does not exist." << endl;
	}
	trk_eta_branch = 0;
	if (tree->GetBranch("trk_eta") != 0) {
		trk_eta_branch = tree->GetBranch("trk_eta");
		trk_eta_branch->SetAddress(&trk_eta_);
	}
	if(trk_eta_branch == 0 ) {
	cout << "Branch trk_eta does not exist." << endl;
	}
	trk_lambda_branch = 0;
	if (tree->GetBranch("trk_lambda") != 0) {
		trk_lambda_branch = tree->GetBranch("trk_lambda");
		trk_lambda_branch->SetAddress(&trk_lambda_);
	}
	if(trk_lambda_branch == 0 ) {
	cout << "Branch trk_lambda does not exist." << endl;
	}
	trk_cotTheta_branch = 0;
	if (tree->GetBranch("trk_cotTheta") != 0) {
		trk_cotTheta_branch = tree->GetBranch("trk_cotTheta");
		trk_cotTheta_branch->SetAddress(&trk_cotTheta_);
	}
	if(trk_cotTheta_branch == 0 ) {
	cout << "Branch trk_cotTheta does not exist." << endl;
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
	trk_lambdaErr_branch = 0;
	if (tree->GetBranch("trk_lambdaErr") != 0) {
		trk_lambdaErr_branch = tree->GetBranch("trk_lambdaErr");
		trk_lambdaErr_branch->SetAddress(&trk_lambdaErr_);
	}
	if(trk_lambdaErr_branch == 0 ) {
	cout << "Branch trk_lambdaErr does not exist." << endl;
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
	trk_refpoint_x_branch = 0;
	if (tree->GetBranch("trk_refpoint_x") != 0) {
		trk_refpoint_x_branch = tree->GetBranch("trk_refpoint_x");
		trk_refpoint_x_branch->SetAddress(&trk_refpoint_x_);
	}
	if(trk_refpoint_x_branch == 0 ) {
	cout << "Branch trk_refpoint_x does not exist." << endl;
	}
	trk_refpoint_y_branch = 0;
	if (tree->GetBranch("trk_refpoint_y") != 0) {
		trk_refpoint_y_branch = tree->GetBranch("trk_refpoint_y");
		trk_refpoint_y_branch->SetAddress(&trk_refpoint_y_);
	}
	if(trk_refpoint_y_branch == 0 ) {
	cout << "Branch trk_refpoint_y does not exist." << endl;
	}
	trk_refpoint_z_branch = 0;
	if (tree->GetBranch("trk_refpoint_z") != 0) {
		trk_refpoint_z_branch = tree->GetBranch("trk_refpoint_z");
		trk_refpoint_z_branch->SetAddress(&trk_refpoint_z_);
	}
	if(trk_refpoint_z_branch == 0 ) {
	cout << "Branch trk_refpoint_z does not exist." << endl;
	}
	trk_nChi2_branch = 0;
	if (tree->GetBranch("trk_nChi2") != 0) {
		trk_nChi2_branch = tree->GetBranch("trk_nChi2");
		trk_nChi2_branch->SetAddress(&trk_nChi2_);
	}
	if(trk_nChi2_branch == 0 ) {
	cout << "Branch trk_nChi2 does not exist." << endl;
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
	trk_nPixelLay_branch = 0;
	if (tree->GetBranch("trk_nPixelLay") != 0) {
		trk_nPixelLay_branch = tree->GetBranch("trk_nPixelLay");
		trk_nPixelLay_branch->SetAddress(&trk_nPixelLay_);
	}
	if(trk_nPixelLay_branch == 0 ) {
	cout << "Branch trk_nPixelLay does not exist." << endl;
	}
	trk_nStripLay_branch = 0;
	if (tree->GetBranch("trk_nStripLay") != 0) {
		trk_nStripLay_branch = tree->GetBranch("trk_nStripLay");
		trk_nStripLay_branch->SetAddress(&trk_nStripLay_);
	}
	if(trk_nStripLay_branch == 0 ) {
	cout << "Branch trk_nStripLay does not exist." << endl;
	}
	trk_n3DLay_branch = 0;
	if (tree->GetBranch("trk_n3DLay") != 0) {
		trk_n3DLay_branch = tree->GetBranch("trk_n3DLay");
		trk_n3DLay_branch->SetAddress(&trk_n3DLay_);
	}
	if(trk_n3DLay_branch == 0 ) {
	cout << "Branch trk_n3DLay does not exist." << endl;
	}
	trk_nOuterLost_branch = 0;
	if (tree->GetBranch("trk_nOuterLost") != 0) {
		trk_nOuterLost_branch = tree->GetBranch("trk_nOuterLost");
		trk_nOuterLost_branch->SetAddress(&trk_nOuterLost_);
	}
	if(trk_nOuterLost_branch == 0 ) {
	cout << "Branch trk_nOuterLost does not exist." << endl;
	}
	trk_nInnerLost_branch = 0;
	if (tree->GetBranch("trk_nInnerLost") != 0) {
		trk_nInnerLost_branch = tree->GetBranch("trk_nInnerLost");
		trk_nInnerLost_branch->SetAddress(&trk_nInnerLost_);
	}
	if(trk_nInnerLost_branch == 0 ) {
	cout << "Branch trk_nInnerLost does not exist." << endl;
	}
	trk_algo_branch = 0;
	if (tree->GetBranch("trk_algo") != 0) {
		trk_algo_branch = tree->GetBranch("trk_algo");
		trk_algo_branch->SetAddress(&trk_algo_);
	}
	if(trk_algo_branch == 0 ) {
	cout << "Branch trk_algo does not exist." << endl;
	}
	trk_originalAlgo_branch = 0;
	if (tree->GetBranch("trk_originalAlgo") != 0) {
		trk_originalAlgo_branch = tree->GetBranch("trk_originalAlgo");
		trk_originalAlgo_branch->SetAddress(&trk_originalAlgo_);
	}
	if(trk_originalAlgo_branch == 0 ) {
	cout << "Branch trk_originalAlgo does not exist." << endl;
	}
	trk_algoMask_branch = 0;
	if (tree->GetBranch("trk_algoMask") != 0) {
		trk_algoMask_branch = tree->GetBranch("trk_algoMask");
		trk_algoMask_branch->SetAddress(&trk_algoMask_);
	}
	if(trk_algoMask_branch == 0 ) {
	cout << "Branch trk_algoMask does not exist." << endl;
	}
	trk_stopReason_branch = 0;
	if (tree->GetBranch("trk_stopReason") != 0) {
		trk_stopReason_branch = tree->GetBranch("trk_stopReason");
		trk_stopReason_branch->SetAddress(&trk_stopReason_);
	}
	if(trk_stopReason_branch == 0 ) {
	cout << "Branch trk_stopReason does not exist." << endl;
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
	trk_vtxIdx_branch = 0;
	if (tree->GetBranch("trk_vtxIdx") != 0) {
		trk_vtxIdx_branch = tree->GetBranch("trk_vtxIdx");
		trk_vtxIdx_branch->SetAddress(&trk_vtxIdx_);
	}
	if(trk_vtxIdx_branch == 0 ) {
	cout << "Branch trk_vtxIdx does not exist." << endl;
	}
	trk_shareFrac_branch = 0;
	if (tree->GetBranch("trk_shareFrac") != 0) {
		trk_shareFrac_branch = tree->GetBranch("trk_shareFrac");
		trk_shareFrac_branch->SetAddress(&trk_shareFrac_);
	}
	if(trk_shareFrac_branch == 0 ) {
	cout << "Branch trk_shareFrac does not exist." << endl;
	}
	trk_simTrkIdx_branch = 0;
	if (tree->GetBranch("trk_simTrkIdx") != 0) {
		trk_simTrkIdx_branch = tree->GetBranch("trk_simTrkIdx");
		trk_simTrkIdx_branch->SetAddress(&trk_simTrkIdx_);
	}
	if(trk_simTrkIdx_branch == 0 ) {
	cout << "Branch trk_simTrkIdx does not exist." << endl;
	}
	trk_hitIdx_branch = 0;
	if (tree->GetBranch("trk_hitIdx") != 0) {
		trk_hitIdx_branch = tree->GetBranch("trk_hitIdx");
		trk_hitIdx_branch->SetAddress(&trk_hitIdx_);
	}
	if(trk_hitIdx_branch == 0 ) {
	cout << "Branch trk_hitIdx does not exist." << endl;
	}
	trk_hitType_branch = 0;
	if (tree->GetBranch("trk_hitType") != 0) {
		trk_hitType_branch = tree->GetBranch("trk_hitType");
		trk_hitType_branch->SetAddress(&trk_hitType_);
	}
	if(trk_hitType_branch == 0 ) {
	cout << "Branch trk_hitType does not exist." << endl;
	}
	sim_event_branch = 0;
	if (tree->GetBranch("sim_event") != 0) {
		sim_event_branch = tree->GetBranch("sim_event");
		sim_event_branch->SetAddress(&sim_event_);
	}
	if(sim_event_branch == 0 ) {
	cout << "Branch sim_event does not exist." << endl;
	}
	sim_bunchCrossing_branch = 0;
	if (tree->GetBranch("sim_bunchCrossing") != 0) {
		sim_bunchCrossing_branch = tree->GetBranch("sim_bunchCrossing");
		sim_bunchCrossing_branch->SetAddress(&sim_bunchCrossing_);
	}
	if(sim_bunchCrossing_branch == 0 ) {
	cout << "Branch sim_bunchCrossing does not exist." << endl;
	}
	sim_pdgId_branch = 0;
	if (tree->GetBranch("sim_pdgId") != 0) {
		sim_pdgId_branch = tree->GetBranch("sim_pdgId");
		sim_pdgId_branch->SetAddress(&sim_pdgId_);
	}
	if(sim_pdgId_branch == 0 ) {
	cout << "Branch sim_pdgId does not exist." << endl;
	}
	sim_genPdgIds_branch = 0;
	if (tree->GetBranch("sim_genPdgIds") != 0) {
		sim_genPdgIds_branch = tree->GetBranch("sim_genPdgIds");
		sim_genPdgIds_branch->SetAddress(&sim_genPdgIds_);
	}
	if(sim_genPdgIds_branch == 0 ) {
	cout << "Branch sim_genPdgIds does not exist." << endl;
	}
	sim_isFromBHadron_branch = 0;
	if (tree->GetBranch("sim_isFromBHadron") != 0) {
		sim_isFromBHadron_branch = tree->GetBranch("sim_isFromBHadron");
		sim_isFromBHadron_branch->SetAddress(&sim_isFromBHadron_);
	}
	if(sim_isFromBHadron_branch == 0 ) {
	cout << "Branch sim_isFromBHadron does not exist." << endl;
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
	sim_pca_pt_branch = 0;
	if (tree->GetBranch("sim_pca_pt") != 0) {
		sim_pca_pt_branch = tree->GetBranch("sim_pca_pt");
		sim_pca_pt_branch->SetAddress(&sim_pca_pt_);
	}
	if(sim_pca_pt_branch == 0 ) {
	cout << "Branch sim_pca_pt does not exist." << endl;
	}
	sim_pca_eta_branch = 0;
	if (tree->GetBranch("sim_pca_eta") != 0) {
		sim_pca_eta_branch = tree->GetBranch("sim_pca_eta");
		sim_pca_eta_branch->SetAddress(&sim_pca_eta_);
	}
	if(sim_pca_eta_branch == 0 ) {
	cout << "Branch sim_pca_eta does not exist." << endl;
	}
	sim_pca_lambda_branch = 0;
	if (tree->GetBranch("sim_pca_lambda") != 0) {
		sim_pca_lambda_branch = tree->GetBranch("sim_pca_lambda");
		sim_pca_lambda_branch->SetAddress(&sim_pca_lambda_);
	}
	if(sim_pca_lambda_branch == 0 ) {
	cout << "Branch sim_pca_lambda does not exist." << endl;
	}
	sim_pca_cotTheta_branch = 0;
	if (tree->GetBranch("sim_pca_cotTheta") != 0) {
		sim_pca_cotTheta_branch = tree->GetBranch("sim_pca_cotTheta");
		sim_pca_cotTheta_branch->SetAddress(&sim_pca_cotTheta_);
	}
	if(sim_pca_cotTheta_branch == 0 ) {
	cout << "Branch sim_pca_cotTheta does not exist." << endl;
	}
	sim_pca_phi_branch = 0;
	if (tree->GetBranch("sim_pca_phi") != 0) {
		sim_pca_phi_branch = tree->GetBranch("sim_pca_phi");
		sim_pca_phi_branch->SetAddress(&sim_pca_phi_);
	}
	if(sim_pca_phi_branch == 0 ) {
	cout << "Branch sim_pca_phi does not exist." << endl;
	}
	sim_pca_dxy_branch = 0;
	if (tree->GetBranch("sim_pca_dxy") != 0) {
		sim_pca_dxy_branch = tree->GetBranch("sim_pca_dxy");
		sim_pca_dxy_branch->SetAddress(&sim_pca_dxy_);
	}
	if(sim_pca_dxy_branch == 0 ) {
	cout << "Branch sim_pca_dxy does not exist." << endl;
	}
	sim_pca_dz_branch = 0;
	if (tree->GetBranch("sim_pca_dz") != 0) {
		sim_pca_dz_branch = tree->GetBranch("sim_pca_dz");
		sim_pca_dz_branch->SetAddress(&sim_pca_dz_);
	}
	if(sim_pca_dz_branch == 0 ) {
	cout << "Branch sim_pca_dz does not exist." << endl;
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
	sim_nLay_branch = 0;
	if (tree->GetBranch("sim_nLay") != 0) {
		sim_nLay_branch = tree->GetBranch("sim_nLay");
		sim_nLay_branch->SetAddress(&sim_nLay_);
	}
	if(sim_nLay_branch == 0 ) {
	cout << "Branch sim_nLay does not exist." << endl;
	}
	sim_nPixelLay_branch = 0;
	if (tree->GetBranch("sim_nPixelLay") != 0) {
		sim_nPixelLay_branch = tree->GetBranch("sim_nPixelLay");
		sim_nPixelLay_branch->SetAddress(&sim_nPixelLay_);
	}
	if(sim_nPixelLay_branch == 0 ) {
	cout << "Branch sim_nPixelLay does not exist." << endl;
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
	sim_shareFrac_branch = 0;
	if (tree->GetBranch("sim_shareFrac") != 0) {
		sim_shareFrac_branch = tree->GetBranch("sim_shareFrac");
		sim_shareFrac_branch->SetAddress(&sim_shareFrac_);
	}
	if(sim_shareFrac_branch == 0 ) {
	cout << "Branch sim_shareFrac does not exist." << endl;
	}
	sim_seedIdx_branch = 0;
	if (tree->GetBranch("sim_seedIdx") != 0) {
		sim_seedIdx_branch = tree->GetBranch("sim_seedIdx");
		sim_seedIdx_branch->SetAddress(&sim_seedIdx_);
	}
	if(sim_seedIdx_branch == 0 ) {
	cout << "Branch sim_seedIdx does not exist." << endl;
	}
	sim_parentVtxIdx_branch = 0;
	if (tree->GetBranch("sim_parentVtxIdx") != 0) {
		sim_parentVtxIdx_branch = tree->GetBranch("sim_parentVtxIdx");
		sim_parentVtxIdx_branch->SetAddress(&sim_parentVtxIdx_);
	}
	if(sim_parentVtxIdx_branch == 0 ) {
	cout << "Branch sim_parentVtxIdx does not exist." << endl;
	}
	sim_decayVtxIdx_branch = 0;
	if (tree->GetBranch("sim_decayVtxIdx") != 0) {
		sim_decayVtxIdx_branch = tree->GetBranch("sim_decayVtxIdx");
		sim_decayVtxIdx_branch->SetAddress(&sim_decayVtxIdx_);
	}
	if(sim_decayVtxIdx_branch == 0 ) {
	cout << "Branch sim_decayVtxIdx does not exist." << endl;
	}
	sim_simHitIdx_branch = 0;
	if (tree->GetBranch("sim_simHitIdx") != 0) {
		sim_simHitIdx_branch = tree->GetBranch("sim_simHitIdx");
		sim_simHitIdx_branch->SetAddress(&sim_simHitIdx_);
	}
	if(sim_simHitIdx_branch == 0 ) {
	cout << "Branch sim_simHitIdx does not exist." << endl;
	}
	pix_isBarrel_branch = 0;
	if (tree->GetBranch("pix_isBarrel") != 0) {
		pix_isBarrel_branch = tree->GetBranch("pix_isBarrel");
		pix_isBarrel_branch->SetAddress(&pix_isBarrel_);
	}
	if(pix_isBarrel_branch == 0 ) {
	cout << "Branch pix_isBarrel does not exist." << endl;
	}
	pix_det_branch = 0;
	if (tree->GetBranch("pix_det") != 0) {
		pix_det_branch = tree->GetBranch("pix_det");
		pix_det_branch->SetAddress(&pix_det_);
	}
	if(pix_det_branch == 0 ) {
	cout << "Branch pix_det does not exist." << endl;
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
	pix_trkIdx_branch = 0;
	if (tree->GetBranch("pix_trkIdx") != 0) {
		pix_trkIdx_branch = tree->GetBranch("pix_trkIdx");
		pix_trkIdx_branch->SetAddress(&pix_trkIdx_);
	}
	if(pix_trkIdx_branch == 0 ) {
	cout << "Branch pix_trkIdx does not exist." << endl;
	}
	pix_seeIdx_branch = 0;
	if (tree->GetBranch("pix_seeIdx") != 0) {
		pix_seeIdx_branch = tree->GetBranch("pix_seeIdx");
		pix_seeIdx_branch->SetAddress(&pix_seeIdx_);
	}
	if(pix_seeIdx_branch == 0 ) {
	cout << "Branch pix_seeIdx does not exist." << endl;
	}
	pix_simHitIdx_branch = 0;
	if (tree->GetBranch("pix_simHitIdx") != 0) {
		pix_simHitIdx_branch = tree->GetBranch("pix_simHitIdx");
		pix_simHitIdx_branch->SetAddress(&pix_simHitIdx_);
	}
	if(pix_simHitIdx_branch == 0 ) {
	cout << "Branch pix_simHitIdx does not exist." << endl;
	}
	pix_xySignificance_branch = 0;
	if (tree->GetBranch("pix_xySignificance") != 0) {
		pix_xySignificance_branch = tree->GetBranch("pix_xySignificance");
		pix_xySignificance_branch->SetAddress(&pix_xySignificance_);
	}
	if(pix_xySignificance_branch == 0 ) {
	cout << "Branch pix_xySignificance does not exist." << endl;
	}
	pix_chargeFraction_branch = 0;
	if (tree->GetBranch("pix_chargeFraction") != 0) {
		pix_chargeFraction_branch = tree->GetBranch("pix_chargeFraction");
		pix_chargeFraction_branch->SetAddress(&pix_chargeFraction_);
	}
	if(pix_chargeFraction_branch == 0 ) {
	cout << "Branch pix_chargeFraction does not exist." << endl;
	}
	pix_simType_branch = 0;
	if (tree->GetBranch("pix_simType") != 0) {
		pix_simType_branch = tree->GetBranch("pix_simType");
		pix_simType_branch->SetAddress(&pix_simType_);
	}
	if(pix_simType_branch == 0 ) {
	cout << "Branch pix_simType does not exist." << endl;
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
	ph2_isBarrel_branch = 0;
	if (tree->GetBranch("ph2_isBarrel") != 0) {
		ph2_isBarrel_branch = tree->GetBranch("ph2_isBarrel");
		ph2_isBarrel_branch->SetAddress(&ph2_isBarrel_);
	}
	if(ph2_isBarrel_branch == 0 ) {
	cout << "Branch ph2_isBarrel does not exist." << endl;
	}
	ph2_det_branch = 0;
	if (tree->GetBranch("ph2_det") != 0) {
		ph2_det_branch = tree->GetBranch("ph2_det");
		ph2_det_branch->SetAddress(&ph2_det_);
	}
	if(ph2_det_branch == 0 ) {
	cout << "Branch ph2_det does not exist." << endl;
	}
	ph2_lay_branch = 0;
	if (tree->GetBranch("ph2_lay") != 0) {
		ph2_lay_branch = tree->GetBranch("ph2_lay");
		ph2_lay_branch->SetAddress(&ph2_lay_);
	}
	if(ph2_lay_branch == 0 ) {
	cout << "Branch ph2_lay does not exist." << endl;
	}
	ph2_detId_branch = 0;
	if (tree->GetBranch("ph2_detId") != 0) {
		ph2_detId_branch = tree->GetBranch("ph2_detId");
		ph2_detId_branch->SetAddress(&ph2_detId_);
	}
	if(ph2_detId_branch == 0 ) {
	cout << "Branch ph2_detId does not exist." << endl;
	}
	ph2_trkIdx_branch = 0;
	if (tree->GetBranch("ph2_trkIdx") != 0) {
		ph2_trkIdx_branch = tree->GetBranch("ph2_trkIdx");
		ph2_trkIdx_branch->SetAddress(&ph2_trkIdx_);
	}
	if(ph2_trkIdx_branch == 0 ) {
	cout << "Branch ph2_trkIdx does not exist." << endl;
	}
	ph2_seeIdx_branch = 0;
	if (tree->GetBranch("ph2_seeIdx") != 0) {
		ph2_seeIdx_branch = tree->GetBranch("ph2_seeIdx");
		ph2_seeIdx_branch->SetAddress(&ph2_seeIdx_);
	}
	if(ph2_seeIdx_branch == 0 ) {
	cout << "Branch ph2_seeIdx does not exist." << endl;
	}
	ph2_simHitIdx_branch = 0;
	if (tree->GetBranch("ph2_simHitIdx") != 0) {
		ph2_simHitIdx_branch = tree->GetBranch("ph2_simHitIdx");
		ph2_simHitIdx_branch->SetAddress(&ph2_simHitIdx_);
	}
	if(ph2_simHitIdx_branch == 0 ) {
	cout << "Branch ph2_simHitIdx does not exist." << endl;
	}
	ph2_xySignificance_branch = 0;
	if (tree->GetBranch("ph2_xySignificance") != 0) {
		ph2_xySignificance_branch = tree->GetBranch("ph2_xySignificance");
		ph2_xySignificance_branch->SetAddress(&ph2_xySignificance_);
	}
	if(ph2_xySignificance_branch == 0 ) {
	cout << "Branch ph2_xySignificance does not exist." << endl;
	}
	ph2_simType_branch = 0;
	if (tree->GetBranch("ph2_simType") != 0) {
		ph2_simType_branch = tree->GetBranch("ph2_simType");
		ph2_simType_branch->SetAddress(&ph2_simType_);
	}
	if(ph2_simType_branch == 0 ) {
	cout << "Branch ph2_simType does not exist." << endl;
	}
	ph2_x_branch = 0;
	if (tree->GetBranch("ph2_x") != 0) {
		ph2_x_branch = tree->GetBranch("ph2_x");
		ph2_x_branch->SetAddress(&ph2_x_);
	}
	if(ph2_x_branch == 0 ) {
	cout << "Branch ph2_x does not exist." << endl;
	}
	ph2_y_branch = 0;
	if (tree->GetBranch("ph2_y") != 0) {
		ph2_y_branch = tree->GetBranch("ph2_y");
		ph2_y_branch->SetAddress(&ph2_y_);
	}
	if(ph2_y_branch == 0 ) {
	cout << "Branch ph2_y does not exist." << endl;
	}
	ph2_z_branch = 0;
	if (tree->GetBranch("ph2_z") != 0) {
		ph2_z_branch = tree->GetBranch("ph2_z");
		ph2_z_branch->SetAddress(&ph2_z_);
	}
	if(ph2_z_branch == 0 ) {
	cout << "Branch ph2_z does not exist." << endl;
	}
	ph2_xx_branch = 0;
	if (tree->GetBranch("ph2_xx") != 0) {
		ph2_xx_branch = tree->GetBranch("ph2_xx");
		ph2_xx_branch->SetAddress(&ph2_xx_);
	}
	if(ph2_xx_branch == 0 ) {
	cout << "Branch ph2_xx does not exist." << endl;
	}
	ph2_xy_branch = 0;
	if (tree->GetBranch("ph2_xy") != 0) {
		ph2_xy_branch = tree->GetBranch("ph2_xy");
		ph2_xy_branch->SetAddress(&ph2_xy_);
	}
	if(ph2_xy_branch == 0 ) {
	cout << "Branch ph2_xy does not exist." << endl;
	}
	ph2_yy_branch = 0;
	if (tree->GetBranch("ph2_yy") != 0) {
		ph2_yy_branch = tree->GetBranch("ph2_yy");
		ph2_yy_branch->SetAddress(&ph2_yy_);
	}
	if(ph2_yy_branch == 0 ) {
	cout << "Branch ph2_yy does not exist." << endl;
	}
	ph2_yz_branch = 0;
	if (tree->GetBranch("ph2_yz") != 0) {
		ph2_yz_branch = tree->GetBranch("ph2_yz");
		ph2_yz_branch->SetAddress(&ph2_yz_);
	}
	if(ph2_yz_branch == 0 ) {
	cout << "Branch ph2_yz does not exist." << endl;
	}
	ph2_zz_branch = 0;
	if (tree->GetBranch("ph2_zz") != 0) {
		ph2_zz_branch = tree->GetBranch("ph2_zz");
		ph2_zz_branch->SetAddress(&ph2_zz_);
	}
	if(ph2_zz_branch == 0 ) {
	cout << "Branch ph2_zz does not exist." << endl;
	}
	ph2_zx_branch = 0;
	if (tree->GetBranch("ph2_zx") != 0) {
		ph2_zx_branch = tree->GetBranch("ph2_zx");
		ph2_zx_branch->SetAddress(&ph2_zx_);
	}
	if(ph2_zx_branch == 0 ) {
	cout << "Branch ph2_zx does not exist." << endl;
	}
	ph2_radL_branch = 0;
	if (tree->GetBranch("ph2_radL") != 0) {
		ph2_radL_branch = tree->GetBranch("ph2_radL");
		ph2_radL_branch->SetAddress(&ph2_radL_);
	}
	if(ph2_radL_branch == 0 ) {
	cout << "Branch ph2_radL does not exist." << endl;
	}
	ph2_bbxi_branch = 0;
	if (tree->GetBranch("ph2_bbxi") != 0) {
		ph2_bbxi_branch = tree->GetBranch("ph2_bbxi");
		ph2_bbxi_branch->SetAddress(&ph2_bbxi_);
	}
	if(ph2_bbxi_branch == 0 ) {
	cout << "Branch ph2_bbxi does not exist." << endl;
	}
	inv_isBarrel_branch = 0;
	if (tree->GetBranch("inv_isBarrel") != 0) {
		inv_isBarrel_branch = tree->GetBranch("inv_isBarrel");
		inv_isBarrel_branch->SetAddress(&inv_isBarrel_);
	}
	if(inv_isBarrel_branch == 0 ) {
	cout << "Branch inv_isBarrel does not exist." << endl;
	}
	inv_det_branch = 0;
	if (tree->GetBranch("inv_det") != 0) {
		inv_det_branch = tree->GetBranch("inv_det");
		inv_det_branch->SetAddress(&inv_det_);
	}
	if(inv_det_branch == 0 ) {
	cout << "Branch inv_det does not exist." << endl;
	}
	inv_lay_branch = 0;
	if (tree->GetBranch("inv_lay") != 0) {
		inv_lay_branch = tree->GetBranch("inv_lay");
		inv_lay_branch->SetAddress(&inv_lay_);
	}
	if(inv_lay_branch == 0 ) {
	cout << "Branch inv_lay does not exist." << endl;
	}
	inv_detId_branch = 0;
	if (tree->GetBranch("inv_detId") != 0) {
		inv_detId_branch = tree->GetBranch("inv_detId");
		inv_detId_branch->SetAddress(&inv_detId_);
	}
	if(inv_detId_branch == 0 ) {
	cout << "Branch inv_detId does not exist." << endl;
	}
	inv_type_branch = 0;
	if (tree->GetBranch("inv_type") != 0) {
		inv_type_branch = tree->GetBranch("inv_type");
		inv_type_branch->SetAddress(&inv_type_);
	}
	if(inv_type_branch == 0 ) {
	cout << "Branch inv_type does not exist." << endl;
	}
	simhit_det_branch = 0;
	if (tree->GetBranch("simhit_det") != 0) {
		simhit_det_branch = tree->GetBranch("simhit_det");
		simhit_det_branch->SetAddress(&simhit_det_);
	}
	if(simhit_det_branch == 0 ) {
	cout << "Branch simhit_det does not exist." << endl;
	}
	simhit_lay_branch = 0;
	if (tree->GetBranch("simhit_lay") != 0) {
		simhit_lay_branch = tree->GetBranch("simhit_lay");
		simhit_lay_branch->SetAddress(&simhit_lay_);
	}
	if(simhit_lay_branch == 0 ) {
	cout << "Branch simhit_lay does not exist." << endl;
	}
	simhit_detId_branch = 0;
	if (tree->GetBranch("simhit_detId") != 0) {
		simhit_detId_branch = tree->GetBranch("simhit_detId");
		simhit_detId_branch->SetAddress(&simhit_detId_);
	}
	if(simhit_detId_branch == 0 ) {
	cout << "Branch simhit_detId does not exist." << endl;
	}
	simhit_x_branch = 0;
	if (tree->GetBranch("simhit_x") != 0) {
		simhit_x_branch = tree->GetBranch("simhit_x");
		simhit_x_branch->SetAddress(&simhit_x_);
	}
	if(simhit_x_branch == 0 ) {
	cout << "Branch simhit_x does not exist." << endl;
	}
	simhit_y_branch = 0;
	if (tree->GetBranch("simhit_y") != 0) {
		simhit_y_branch = tree->GetBranch("simhit_y");
		simhit_y_branch->SetAddress(&simhit_y_);
	}
	if(simhit_y_branch == 0 ) {
	cout << "Branch simhit_y does not exist." << endl;
	}
	simhit_z_branch = 0;
	if (tree->GetBranch("simhit_z") != 0) {
		simhit_z_branch = tree->GetBranch("simhit_z");
		simhit_z_branch->SetAddress(&simhit_z_);
	}
	if(simhit_z_branch == 0 ) {
	cout << "Branch simhit_z does not exist." << endl;
	}
	simhit_px_branch = 0;
	if (tree->GetBranch("simhit_px") != 0) {
		simhit_px_branch = tree->GetBranch("simhit_px");
		simhit_px_branch->SetAddress(&simhit_px_);
	}
	if(simhit_px_branch == 0 ) {
	cout << "Branch simhit_px does not exist." << endl;
	}
	simhit_py_branch = 0;
	if (tree->GetBranch("simhit_py") != 0) {
		simhit_py_branch = tree->GetBranch("simhit_py");
		simhit_py_branch->SetAddress(&simhit_py_);
	}
	if(simhit_py_branch == 0 ) {
	cout << "Branch simhit_py does not exist." << endl;
	}
	simhit_pz_branch = 0;
	if (tree->GetBranch("simhit_pz") != 0) {
		simhit_pz_branch = tree->GetBranch("simhit_pz");
		simhit_pz_branch->SetAddress(&simhit_pz_);
	}
	if(simhit_pz_branch == 0 ) {
	cout << "Branch simhit_pz does not exist." << endl;
	}
	simhit_particle_branch = 0;
	if (tree->GetBranch("simhit_particle") != 0) {
		simhit_particle_branch = tree->GetBranch("simhit_particle");
		simhit_particle_branch->SetAddress(&simhit_particle_);
	}
	if(simhit_particle_branch == 0 ) {
	cout << "Branch simhit_particle does not exist." << endl;
	}
	simhit_process_branch = 0;
	if (tree->GetBranch("simhit_process") != 0) {
		simhit_process_branch = tree->GetBranch("simhit_process");
		simhit_process_branch->SetAddress(&simhit_process_);
	}
	if(simhit_process_branch == 0 ) {
	cout << "Branch simhit_process does not exist." << endl;
	}
	simhit_eloss_branch = 0;
	if (tree->GetBranch("simhit_eloss") != 0) {
		simhit_eloss_branch = tree->GetBranch("simhit_eloss");
		simhit_eloss_branch->SetAddress(&simhit_eloss_);
	}
	if(simhit_eloss_branch == 0 ) {
	cout << "Branch simhit_eloss does not exist." << endl;
	}
	simhit_tof_branch = 0;
	if (tree->GetBranch("simhit_tof") != 0) {
		simhit_tof_branch = tree->GetBranch("simhit_tof");
		simhit_tof_branch->SetAddress(&simhit_tof_);
	}
	if(simhit_tof_branch == 0 ) {
	cout << "Branch simhit_tof does not exist." << endl;
	}
	simhit_simTrkIdx_branch = 0;
	if (tree->GetBranch("simhit_simTrkIdx") != 0) {
		simhit_simTrkIdx_branch = tree->GetBranch("simhit_simTrkIdx");
		simhit_simTrkIdx_branch->SetAddress(&simhit_simTrkIdx_);
	}
	if(simhit_simTrkIdx_branch == 0 ) {
	cout << "Branch simhit_simTrkIdx does not exist." << endl;
	}
	simhit_hitIdx_branch = 0;
	if (tree->GetBranch("simhit_hitIdx") != 0) {
		simhit_hitIdx_branch = tree->GetBranch("simhit_hitIdx");
		simhit_hitIdx_branch->SetAddress(&simhit_hitIdx_);
	}
	if(simhit_hitIdx_branch == 0 ) {
	cout << "Branch simhit_hitIdx does not exist." << endl;
	}
	simhit_hitType_branch = 0;
	if (tree->GetBranch("simhit_hitType") != 0) {
		simhit_hitType_branch = tree->GetBranch("simhit_hitType");
		simhit_hitType_branch->SetAddress(&simhit_hitType_);
	}
	if(simhit_hitType_branch == 0 ) {
	cout << "Branch simhit_hitType does not exist." << endl;
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
	see_fitok_branch = 0;
	if (tree->GetBranch("see_fitok") != 0) {
		see_fitok_branch = tree->GetBranch("see_fitok");
		see_fitok_branch->SetAddress(&see_fitok_);
	}
	if(see_fitok_branch == 0 ) {
	cout << "Branch see_fitok does not exist." << endl;
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
	see_statePt_branch = 0;
	if (tree->GetBranch("see_statePt") != 0) {
		see_statePt_branch = tree->GetBranch("see_statePt");
		see_statePt_branch->SetAddress(&see_statePt_);
	}
	if(see_statePt_branch == 0 ) {
	cout << "Branch see_statePt does not exist." << endl;
	}
	see_stateTrajX_branch = 0;
	if (tree->GetBranch("see_stateTrajX") != 0) {
		see_stateTrajX_branch = tree->GetBranch("see_stateTrajX");
		see_stateTrajX_branch->SetAddress(&see_stateTrajX_);
	}
	if(see_stateTrajX_branch == 0 ) {
	cout << "Branch see_stateTrajX does not exist." << endl;
	}
	see_stateTrajY_branch = 0;
	if (tree->GetBranch("see_stateTrajY") != 0) {
		see_stateTrajY_branch = tree->GetBranch("see_stateTrajY");
		see_stateTrajY_branch->SetAddress(&see_stateTrajY_);
	}
	if(see_stateTrajY_branch == 0 ) {
	cout << "Branch see_stateTrajY does not exist." << endl;
	}
	see_stateTrajPx_branch = 0;
	if (tree->GetBranch("see_stateTrajPx") != 0) {
		see_stateTrajPx_branch = tree->GetBranch("see_stateTrajPx");
		see_stateTrajPx_branch->SetAddress(&see_stateTrajPx_);
	}
	if(see_stateTrajPx_branch == 0 ) {
	cout << "Branch see_stateTrajPx does not exist." << endl;
	}
	see_stateTrajPy_branch = 0;
	if (tree->GetBranch("see_stateTrajPy") != 0) {
		see_stateTrajPy_branch = tree->GetBranch("see_stateTrajPy");
		see_stateTrajPy_branch->SetAddress(&see_stateTrajPy_);
	}
	if(see_stateTrajPy_branch == 0 ) {
	cout << "Branch see_stateTrajPy does not exist." << endl;
	}
	see_stateTrajPz_branch = 0;
	if (tree->GetBranch("see_stateTrajPz") != 0) {
		see_stateTrajPz_branch = tree->GetBranch("see_stateTrajPz");
		see_stateTrajPz_branch->SetAddress(&see_stateTrajPz_);
	}
	if(see_stateTrajPz_branch == 0 ) {
	cout << "Branch see_stateTrajPz does not exist." << endl;
	}
	see_stateTrajGlbX_branch = 0;
	if (tree->GetBranch("see_stateTrajGlbX") != 0) {
		see_stateTrajGlbX_branch = tree->GetBranch("see_stateTrajGlbX");
		see_stateTrajGlbX_branch->SetAddress(&see_stateTrajGlbX_);
	}
	if(see_stateTrajGlbX_branch == 0 ) {
	cout << "Branch see_stateTrajGlbX does not exist." << endl;
	}
	see_stateTrajGlbY_branch = 0;
	if (tree->GetBranch("see_stateTrajGlbY") != 0) {
		see_stateTrajGlbY_branch = tree->GetBranch("see_stateTrajGlbY");
		see_stateTrajGlbY_branch->SetAddress(&see_stateTrajGlbY_);
	}
	if(see_stateTrajGlbY_branch == 0 ) {
	cout << "Branch see_stateTrajGlbY does not exist." << endl;
	}
	see_stateTrajGlbZ_branch = 0;
	if (tree->GetBranch("see_stateTrajGlbZ") != 0) {
		see_stateTrajGlbZ_branch = tree->GetBranch("see_stateTrajGlbZ");
		see_stateTrajGlbZ_branch->SetAddress(&see_stateTrajGlbZ_);
	}
	if(see_stateTrajGlbZ_branch == 0 ) {
	cout << "Branch see_stateTrajGlbZ does not exist." << endl;
	}
	see_stateTrajGlbPx_branch = 0;
	if (tree->GetBranch("see_stateTrajGlbPx") != 0) {
		see_stateTrajGlbPx_branch = tree->GetBranch("see_stateTrajGlbPx");
		see_stateTrajGlbPx_branch->SetAddress(&see_stateTrajGlbPx_);
	}
	if(see_stateTrajGlbPx_branch == 0 ) {
	cout << "Branch see_stateTrajGlbPx does not exist." << endl;
	}
	see_stateTrajGlbPy_branch = 0;
	if (tree->GetBranch("see_stateTrajGlbPy") != 0) {
		see_stateTrajGlbPy_branch = tree->GetBranch("see_stateTrajGlbPy");
		see_stateTrajGlbPy_branch->SetAddress(&see_stateTrajGlbPy_);
	}
	if(see_stateTrajGlbPy_branch == 0 ) {
	cout << "Branch see_stateTrajGlbPy does not exist." << endl;
	}
	see_stateTrajGlbPz_branch = 0;
	if (tree->GetBranch("see_stateTrajGlbPz") != 0) {
		see_stateTrajGlbPz_branch = tree->GetBranch("see_stateTrajGlbPz");
		see_stateTrajGlbPz_branch->SetAddress(&see_stateTrajGlbPz_);
	}
	if(see_stateTrajGlbPz_branch == 0 ) {
	cout << "Branch see_stateTrajGlbPz does not exist." << endl;
	}
	see_stateCcov00_branch = 0;
	if (tree->GetBranch("see_stateCcov00") != 0) {
		see_stateCcov00_branch = tree->GetBranch("see_stateCcov00");
		see_stateCcov00_branch->SetAddress(&see_stateCcov00_);
	}
	if(see_stateCcov00_branch == 0 ) {
	cout << "Branch see_stateCcov00 does not exist." << endl;
	}
	see_stateCcov01_branch = 0;
	if (tree->GetBranch("see_stateCcov01") != 0) {
		see_stateCcov01_branch = tree->GetBranch("see_stateCcov01");
		see_stateCcov01_branch->SetAddress(&see_stateCcov01_);
	}
	if(see_stateCcov01_branch == 0 ) {
	cout << "Branch see_stateCcov01 does not exist." << endl;
	}
	see_stateCcov02_branch = 0;
	if (tree->GetBranch("see_stateCcov02") != 0) {
		see_stateCcov02_branch = tree->GetBranch("see_stateCcov02");
		see_stateCcov02_branch->SetAddress(&see_stateCcov02_);
	}
	if(see_stateCcov02_branch == 0 ) {
	cout << "Branch see_stateCcov02 does not exist." << endl;
	}
	see_stateCcov03_branch = 0;
	if (tree->GetBranch("see_stateCcov03") != 0) {
		see_stateCcov03_branch = tree->GetBranch("see_stateCcov03");
		see_stateCcov03_branch->SetAddress(&see_stateCcov03_);
	}
	if(see_stateCcov03_branch == 0 ) {
	cout << "Branch see_stateCcov03 does not exist." << endl;
	}
	see_stateCcov04_branch = 0;
	if (tree->GetBranch("see_stateCcov04") != 0) {
		see_stateCcov04_branch = tree->GetBranch("see_stateCcov04");
		see_stateCcov04_branch->SetAddress(&see_stateCcov04_);
	}
	if(see_stateCcov04_branch == 0 ) {
	cout << "Branch see_stateCcov04 does not exist." << endl;
	}
	see_stateCcov05_branch = 0;
	if (tree->GetBranch("see_stateCcov05") != 0) {
		see_stateCcov05_branch = tree->GetBranch("see_stateCcov05");
		see_stateCcov05_branch->SetAddress(&see_stateCcov05_);
	}
	if(see_stateCcov05_branch == 0 ) {
	cout << "Branch see_stateCcov05 does not exist." << endl;
	}
	see_stateCcov11_branch = 0;
	if (tree->GetBranch("see_stateCcov11") != 0) {
		see_stateCcov11_branch = tree->GetBranch("see_stateCcov11");
		see_stateCcov11_branch->SetAddress(&see_stateCcov11_);
	}
	if(see_stateCcov11_branch == 0 ) {
	cout << "Branch see_stateCcov11 does not exist." << endl;
	}
	see_stateCcov12_branch = 0;
	if (tree->GetBranch("see_stateCcov12") != 0) {
		see_stateCcov12_branch = tree->GetBranch("see_stateCcov12");
		see_stateCcov12_branch->SetAddress(&see_stateCcov12_);
	}
	if(see_stateCcov12_branch == 0 ) {
	cout << "Branch see_stateCcov12 does not exist." << endl;
	}
	see_stateCcov13_branch = 0;
	if (tree->GetBranch("see_stateCcov13") != 0) {
		see_stateCcov13_branch = tree->GetBranch("see_stateCcov13");
		see_stateCcov13_branch->SetAddress(&see_stateCcov13_);
	}
	if(see_stateCcov13_branch == 0 ) {
	cout << "Branch see_stateCcov13 does not exist." << endl;
	}
	see_stateCcov14_branch = 0;
	if (tree->GetBranch("see_stateCcov14") != 0) {
		see_stateCcov14_branch = tree->GetBranch("see_stateCcov14");
		see_stateCcov14_branch->SetAddress(&see_stateCcov14_);
	}
	if(see_stateCcov14_branch == 0 ) {
	cout << "Branch see_stateCcov14 does not exist." << endl;
	}
	see_stateCcov15_branch = 0;
	if (tree->GetBranch("see_stateCcov15") != 0) {
		see_stateCcov15_branch = tree->GetBranch("see_stateCcov15");
		see_stateCcov15_branch->SetAddress(&see_stateCcov15_);
	}
	if(see_stateCcov15_branch == 0 ) {
	cout << "Branch see_stateCcov15 does not exist." << endl;
	}
	see_stateCcov22_branch = 0;
	if (tree->GetBranch("see_stateCcov22") != 0) {
		see_stateCcov22_branch = tree->GetBranch("see_stateCcov22");
		see_stateCcov22_branch->SetAddress(&see_stateCcov22_);
	}
	if(see_stateCcov22_branch == 0 ) {
	cout << "Branch see_stateCcov22 does not exist." << endl;
	}
	see_stateCcov23_branch = 0;
	if (tree->GetBranch("see_stateCcov23") != 0) {
		see_stateCcov23_branch = tree->GetBranch("see_stateCcov23");
		see_stateCcov23_branch->SetAddress(&see_stateCcov23_);
	}
	if(see_stateCcov23_branch == 0 ) {
	cout << "Branch see_stateCcov23 does not exist." << endl;
	}
	see_stateCcov24_branch = 0;
	if (tree->GetBranch("see_stateCcov24") != 0) {
		see_stateCcov24_branch = tree->GetBranch("see_stateCcov24");
		see_stateCcov24_branch->SetAddress(&see_stateCcov24_);
	}
	if(see_stateCcov24_branch == 0 ) {
	cout << "Branch see_stateCcov24 does not exist." << endl;
	}
	see_stateCcov25_branch = 0;
	if (tree->GetBranch("see_stateCcov25") != 0) {
		see_stateCcov25_branch = tree->GetBranch("see_stateCcov25");
		see_stateCcov25_branch->SetAddress(&see_stateCcov25_);
	}
	if(see_stateCcov25_branch == 0 ) {
	cout << "Branch see_stateCcov25 does not exist." << endl;
	}
	see_stateCcov33_branch = 0;
	if (tree->GetBranch("see_stateCcov33") != 0) {
		see_stateCcov33_branch = tree->GetBranch("see_stateCcov33");
		see_stateCcov33_branch->SetAddress(&see_stateCcov33_);
	}
	if(see_stateCcov33_branch == 0 ) {
	cout << "Branch see_stateCcov33 does not exist." << endl;
	}
	see_stateCcov34_branch = 0;
	if (tree->GetBranch("see_stateCcov34") != 0) {
		see_stateCcov34_branch = tree->GetBranch("see_stateCcov34");
		see_stateCcov34_branch->SetAddress(&see_stateCcov34_);
	}
	if(see_stateCcov34_branch == 0 ) {
	cout << "Branch see_stateCcov34 does not exist." << endl;
	}
	see_stateCcov35_branch = 0;
	if (tree->GetBranch("see_stateCcov35") != 0) {
		see_stateCcov35_branch = tree->GetBranch("see_stateCcov35");
		see_stateCcov35_branch->SetAddress(&see_stateCcov35_);
	}
	if(see_stateCcov35_branch == 0 ) {
	cout << "Branch see_stateCcov35 does not exist." << endl;
	}
	see_stateCcov44_branch = 0;
	if (tree->GetBranch("see_stateCcov44") != 0) {
		see_stateCcov44_branch = tree->GetBranch("see_stateCcov44");
		see_stateCcov44_branch->SetAddress(&see_stateCcov44_);
	}
	if(see_stateCcov44_branch == 0 ) {
	cout << "Branch see_stateCcov44 does not exist." << endl;
	}
	see_stateCcov45_branch = 0;
	if (tree->GetBranch("see_stateCcov45") != 0) {
		see_stateCcov45_branch = tree->GetBranch("see_stateCcov45");
		see_stateCcov45_branch->SetAddress(&see_stateCcov45_);
	}
	if(see_stateCcov45_branch == 0 ) {
	cout << "Branch see_stateCcov45 does not exist." << endl;
	}
	see_stateCcov55_branch = 0;
	if (tree->GetBranch("see_stateCcov55") != 0) {
		see_stateCcov55_branch = tree->GetBranch("see_stateCcov55");
		see_stateCcov55_branch->SetAddress(&see_stateCcov55_);
	}
	if(see_stateCcov55_branch == 0 ) {
	cout << "Branch see_stateCcov55 does not exist." << endl;
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
	see_nPhase2OT_branch = 0;
	if (tree->GetBranch("see_nPhase2OT") != 0) {
		see_nPhase2OT_branch = tree->GetBranch("see_nPhase2OT");
		see_nPhase2OT_branch->SetAddress(&see_nPhase2OT_);
	}
	if(see_nPhase2OT_branch == 0 ) {
	cout << "Branch see_nPhase2OT does not exist." << endl;
	}
	see_algo_branch = 0;
	if (tree->GetBranch("see_algo") != 0) {
		see_algo_branch = tree->GetBranch("see_algo");
		see_algo_branch->SetAddress(&see_algo_);
	}
	if(see_algo_branch == 0 ) {
	cout << "Branch see_algo does not exist." << endl;
	}
	see_stopReason_branch = 0;
	if (tree->GetBranch("see_stopReason") != 0) {
		see_stopReason_branch = tree->GetBranch("see_stopReason");
		see_stopReason_branch->SetAddress(&see_stopReason_);
	}
	if(see_stopReason_branch == 0 ) {
	cout << "Branch see_stopReason does not exist." << endl;
	}
	see_trkIdx_branch = 0;
	if (tree->GetBranch("see_trkIdx") != 0) {
		see_trkIdx_branch = tree->GetBranch("see_trkIdx");
		see_trkIdx_branch->SetAddress(&see_trkIdx_);
	}
	if(see_trkIdx_branch == 0 ) {
	cout << "Branch see_trkIdx does not exist." << endl;
	}
	see_shareFrac_branch = 0;
	if (tree->GetBranch("see_shareFrac") != 0) {
		see_shareFrac_branch = tree->GetBranch("see_shareFrac");
		see_shareFrac_branch->SetAddress(&see_shareFrac_);
	}
	if(see_shareFrac_branch == 0 ) {
	cout << "Branch see_shareFrac does not exist." << endl;
	}
	see_simTrkIdx_branch = 0;
	if (tree->GetBranch("see_simTrkIdx") != 0) {
		see_simTrkIdx_branch = tree->GetBranch("see_simTrkIdx");
		see_simTrkIdx_branch->SetAddress(&see_simTrkIdx_);
	}
	if(see_simTrkIdx_branch == 0 ) {
	cout << "Branch see_simTrkIdx does not exist." << endl;
	}
	see_hitIdx_branch = 0;
	if (tree->GetBranch("see_hitIdx") != 0) {
		see_hitIdx_branch = tree->GetBranch("see_hitIdx");
		see_hitIdx_branch->SetAddress(&see_hitIdx_);
	}
	if(see_hitIdx_branch == 0 ) {
	cout << "Branch see_hitIdx does not exist." << endl;
	}
	see_hitType_branch = 0;
	if (tree->GetBranch("see_hitType") != 0) {
		see_hitType_branch = tree->GetBranch("see_hitType");
		see_hitType_branch->SetAddress(&see_hitType_);
	}
	if(see_hitType_branch == 0 ) {
	cout << "Branch see_hitType does not exist." << endl;
	}
	see_offset_branch = 0;
	if (tree->GetBranch("see_offset") != 0) {
		see_offset_branch = tree->GetBranch("see_offset");
		see_offset_branch->SetAddress(&see_offset_);
	}
	if(see_offset_branch == 0 ) {
	cout << "Branch see_offset does not exist." << endl;
	}
	vtx_x_branch = 0;
	if (tree->GetBranch("vtx_x") != 0) {
		vtx_x_branch = tree->GetBranch("vtx_x");
		vtx_x_branch->SetAddress(&vtx_x_);
	}
	if(vtx_x_branch == 0 ) {
	cout << "Branch vtx_x does not exist." << endl;
	}
	vtx_y_branch = 0;
	if (tree->GetBranch("vtx_y") != 0) {
		vtx_y_branch = tree->GetBranch("vtx_y");
		vtx_y_branch->SetAddress(&vtx_y_);
	}
	if(vtx_y_branch == 0 ) {
	cout << "Branch vtx_y does not exist." << endl;
	}
	vtx_z_branch = 0;
	if (tree->GetBranch("vtx_z") != 0) {
		vtx_z_branch = tree->GetBranch("vtx_z");
		vtx_z_branch->SetAddress(&vtx_z_);
	}
	if(vtx_z_branch == 0 ) {
	cout << "Branch vtx_z does not exist." << endl;
	}
	vtx_xErr_branch = 0;
	if (tree->GetBranch("vtx_xErr") != 0) {
		vtx_xErr_branch = tree->GetBranch("vtx_xErr");
		vtx_xErr_branch->SetAddress(&vtx_xErr_);
	}
	if(vtx_xErr_branch == 0 ) {
	cout << "Branch vtx_xErr does not exist." << endl;
	}
	vtx_yErr_branch = 0;
	if (tree->GetBranch("vtx_yErr") != 0) {
		vtx_yErr_branch = tree->GetBranch("vtx_yErr");
		vtx_yErr_branch->SetAddress(&vtx_yErr_);
	}
	if(vtx_yErr_branch == 0 ) {
	cout << "Branch vtx_yErr does not exist." << endl;
	}
	vtx_zErr_branch = 0;
	if (tree->GetBranch("vtx_zErr") != 0) {
		vtx_zErr_branch = tree->GetBranch("vtx_zErr");
		vtx_zErr_branch->SetAddress(&vtx_zErr_);
	}
	if(vtx_zErr_branch == 0 ) {
	cout << "Branch vtx_zErr does not exist." << endl;
	}
	vtx_ndof_branch = 0;
	if (tree->GetBranch("vtx_ndof") != 0) {
		vtx_ndof_branch = tree->GetBranch("vtx_ndof");
		vtx_ndof_branch->SetAddress(&vtx_ndof_);
	}
	if(vtx_ndof_branch == 0 ) {
	cout << "Branch vtx_ndof does not exist." << endl;
	}
	vtx_chi2_branch = 0;
	if (tree->GetBranch("vtx_chi2") != 0) {
		vtx_chi2_branch = tree->GetBranch("vtx_chi2");
		vtx_chi2_branch->SetAddress(&vtx_chi2_);
	}
	if(vtx_chi2_branch == 0 ) {
	cout << "Branch vtx_chi2 does not exist." << endl;
	}
	vtx_fake_branch = 0;
	if (tree->GetBranch("vtx_fake") != 0) {
		vtx_fake_branch = tree->GetBranch("vtx_fake");
		vtx_fake_branch->SetAddress(&vtx_fake_);
	}
	if(vtx_fake_branch == 0 ) {
	cout << "Branch vtx_fake does not exist." << endl;
	}
	vtx_valid_branch = 0;
	if (tree->GetBranch("vtx_valid") != 0) {
		vtx_valid_branch = tree->GetBranch("vtx_valid");
		vtx_valid_branch->SetAddress(&vtx_valid_);
	}
	if(vtx_valid_branch == 0 ) {
	cout << "Branch vtx_valid does not exist." << endl;
	}
	vtx_trkIdx_branch = 0;
	if (tree->GetBranch("vtx_trkIdx") != 0) {
		vtx_trkIdx_branch = tree->GetBranch("vtx_trkIdx");
		vtx_trkIdx_branch->SetAddress(&vtx_trkIdx_);
	}
	if(vtx_trkIdx_branch == 0 ) {
	cout << "Branch vtx_trkIdx does not exist." << endl;
	}
	simvtx_event_branch = 0;
	if (tree->GetBranch("simvtx_event") != 0) {
		simvtx_event_branch = tree->GetBranch("simvtx_event");
		simvtx_event_branch->SetAddress(&simvtx_event_);
	}
	if(simvtx_event_branch == 0 ) {
	cout << "Branch simvtx_event does not exist." << endl;
	}
	simvtx_bunchCrossing_branch = 0;
	if (tree->GetBranch("simvtx_bunchCrossing") != 0) {
		simvtx_bunchCrossing_branch = tree->GetBranch("simvtx_bunchCrossing");
		simvtx_bunchCrossing_branch->SetAddress(&simvtx_bunchCrossing_);
	}
	if(simvtx_bunchCrossing_branch == 0 ) {
	cout << "Branch simvtx_bunchCrossing does not exist." << endl;
	}
	simvtx_processType_branch = 0;
	if (tree->GetBranch("simvtx_processType") != 0) {
		simvtx_processType_branch = tree->GetBranch("simvtx_processType");
		simvtx_processType_branch->SetAddress(&simvtx_processType_);
	}
	if(simvtx_processType_branch == 0 ) {
	cout << "Branch simvtx_processType does not exist." << endl;
	}
	simvtx_x_branch = 0;
	if (tree->GetBranch("simvtx_x") != 0) {
		simvtx_x_branch = tree->GetBranch("simvtx_x");
		simvtx_x_branch->SetAddress(&simvtx_x_);
	}
	if(simvtx_x_branch == 0 ) {
	cout << "Branch simvtx_x does not exist." << endl;
	}
	simvtx_y_branch = 0;
	if (tree->GetBranch("simvtx_y") != 0) {
		simvtx_y_branch = tree->GetBranch("simvtx_y");
		simvtx_y_branch->SetAddress(&simvtx_y_);
	}
	if(simvtx_y_branch == 0 ) {
	cout << "Branch simvtx_y does not exist." << endl;
	}
	simvtx_z_branch = 0;
	if (tree->GetBranch("simvtx_z") != 0) {
		simvtx_z_branch = tree->GetBranch("simvtx_z");
		simvtx_z_branch->SetAddress(&simvtx_z_);
	}
	if(simvtx_z_branch == 0 ) {
	cout << "Branch simvtx_z does not exist." << endl;
	}
	simvtx_sourceSimIdx_branch = 0;
	if (tree->GetBranch("simvtx_sourceSimIdx") != 0) {
		simvtx_sourceSimIdx_branch = tree->GetBranch("simvtx_sourceSimIdx");
		simvtx_sourceSimIdx_branch->SetAddress(&simvtx_sourceSimIdx_);
	}
	if(simvtx_sourceSimIdx_branch == 0 ) {
	cout << "Branch simvtx_sourceSimIdx does not exist." << endl;
	}
	simvtx_daughterSimIdx_branch = 0;
	if (tree->GetBranch("simvtx_daughterSimIdx") != 0) {
		simvtx_daughterSimIdx_branch = tree->GetBranch("simvtx_daughterSimIdx");
		simvtx_daughterSimIdx_branch->SetAddress(&simvtx_daughterSimIdx_);
	}
	if(simvtx_daughterSimIdx_branch == 0 ) {
	cout << "Branch simvtx_daughterSimIdx does not exist." << endl;
	}
	simpv_idx_branch = 0;
	if (tree->GetBranch("simpv_idx") != 0) {
		simpv_idx_branch = tree->GetBranch("simpv_idx");
		simpv_idx_branch->SetAddress(&simpv_idx_);
	}
	if(simpv_idx_branch == 0 ) {
	cout << "Branch simpv_idx does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		event_isLoaded = false;
		lumi_isLoaded = false;
		run_isLoaded = false;
		trk_px_isLoaded = false;
		trk_py_isLoaded = false;
		trk_pz_isLoaded = false;
		trk_pt_isLoaded = false;
		trk_inner_px_isLoaded = false;
		trk_inner_py_isLoaded = false;
		trk_inner_pz_isLoaded = false;
		trk_inner_pt_isLoaded = false;
		trk_outer_px_isLoaded = false;
		trk_outer_py_isLoaded = false;
		trk_outer_pz_isLoaded = false;
		trk_outer_pt_isLoaded = false;
		trk_eta_isLoaded = false;
		trk_lambda_isLoaded = false;
		trk_cotTheta_isLoaded = false;
		trk_phi_isLoaded = false;
		trk_dxy_isLoaded = false;
		trk_dz_isLoaded = false;
		trk_ptErr_isLoaded = false;
		trk_etaErr_isLoaded = false;
		trk_lambdaErr_isLoaded = false;
		trk_phiErr_isLoaded = false;
		trk_dxyErr_isLoaded = false;
		trk_dzErr_isLoaded = false;
		trk_refpoint_x_isLoaded = false;
		trk_refpoint_y_isLoaded = false;
		trk_refpoint_z_isLoaded = false;
		trk_nChi2_isLoaded = false;
		trk_q_isLoaded = false;
		trk_nValid_isLoaded = false;
		trk_nInvalid_isLoaded = false;
		trk_nPixel_isLoaded = false;
		trk_nStrip_isLoaded = false;
		trk_nPixelLay_isLoaded = false;
		trk_nStripLay_isLoaded = false;
		trk_n3DLay_isLoaded = false;
		trk_nOuterLost_isLoaded = false;
		trk_nInnerLost_isLoaded = false;
		trk_algo_isLoaded = false;
		trk_originalAlgo_isLoaded = false;
		trk_algoMask_isLoaded = false;
		trk_stopReason_isLoaded = false;
		trk_isHP_isLoaded = false;
		trk_seedIdx_isLoaded = false;
		trk_vtxIdx_isLoaded = false;
		trk_shareFrac_isLoaded = false;
		trk_simTrkIdx_isLoaded = false;
		trk_hitIdx_isLoaded = false;
		trk_hitType_isLoaded = false;
		sim_event_isLoaded = false;
		sim_bunchCrossing_isLoaded = false;
		sim_pdgId_isLoaded = false;
		sim_genPdgIds_isLoaded = false;
		sim_isFromBHadron_isLoaded = false;
		sim_px_isLoaded = false;
		sim_py_isLoaded = false;
		sim_pz_isLoaded = false;
		sim_pt_isLoaded = false;
		sim_eta_isLoaded = false;
		sim_phi_isLoaded = false;
		sim_pca_pt_isLoaded = false;
		sim_pca_eta_isLoaded = false;
		sim_pca_lambda_isLoaded = false;
		sim_pca_cotTheta_isLoaded = false;
		sim_pca_phi_isLoaded = false;
		sim_pca_dxy_isLoaded = false;
		sim_pca_dz_isLoaded = false;
		sim_q_isLoaded = false;
		sim_nValid_isLoaded = false;
		sim_nPixel_isLoaded = false;
		sim_nStrip_isLoaded = false;
		sim_nLay_isLoaded = false;
		sim_nPixelLay_isLoaded = false;
		sim_n3DLay_isLoaded = false;
		sim_trkIdx_isLoaded = false;
		sim_shareFrac_isLoaded = false;
		sim_seedIdx_isLoaded = false;
		sim_parentVtxIdx_isLoaded = false;
		sim_decayVtxIdx_isLoaded = false;
		sim_simHitIdx_isLoaded = false;
		pix_isBarrel_isLoaded = false;
		pix_det_isLoaded = false;
		pix_lay_isLoaded = false;
		pix_detId_isLoaded = false;
		pix_trkIdx_isLoaded = false;
		pix_seeIdx_isLoaded = false;
		pix_simHitIdx_isLoaded = false;
		pix_xySignificance_isLoaded = false;
		pix_chargeFraction_isLoaded = false;
		pix_simType_isLoaded = false;
		pix_x_isLoaded = false;
		pix_y_isLoaded = false;
		pix_z_isLoaded = false;
		pix_xx_isLoaded = false;
		pix_xy_isLoaded = false;
		pix_yy_isLoaded = false;
		pix_yz_isLoaded = false;
		pix_zz_isLoaded = false;
		pix_zx_isLoaded = false;
		pix_radL_isLoaded = false;
		pix_bbxi_isLoaded = false;
		ph2_isBarrel_isLoaded = false;
		ph2_det_isLoaded = false;
		ph2_lay_isLoaded = false;
		ph2_detId_isLoaded = false;
		ph2_trkIdx_isLoaded = false;
		ph2_seeIdx_isLoaded = false;
		ph2_simHitIdx_isLoaded = false;
		ph2_xySignificance_isLoaded = false;
		ph2_simType_isLoaded = false;
		ph2_x_isLoaded = false;
		ph2_y_isLoaded = false;
		ph2_z_isLoaded = false;
		ph2_xx_isLoaded = false;
		ph2_xy_isLoaded = false;
		ph2_yy_isLoaded = false;
		ph2_yz_isLoaded = false;
		ph2_zz_isLoaded = false;
		ph2_zx_isLoaded = false;
		ph2_radL_isLoaded = false;
		ph2_bbxi_isLoaded = false;
		inv_isBarrel_isLoaded = false;
		inv_det_isLoaded = false;
		inv_lay_isLoaded = false;
		inv_detId_isLoaded = false;
		inv_type_isLoaded = false;
		simhit_det_isLoaded = false;
		simhit_lay_isLoaded = false;
		simhit_detId_isLoaded = false;
		simhit_x_isLoaded = false;
		simhit_y_isLoaded = false;
		simhit_z_isLoaded = false;
		simhit_px_isLoaded = false;
		simhit_py_isLoaded = false;
		simhit_pz_isLoaded = false;
		simhit_particle_isLoaded = false;
		simhit_process_isLoaded = false;
		simhit_eloss_isLoaded = false;
		simhit_tof_isLoaded = false;
		simhit_simTrkIdx_isLoaded = false;
		simhit_hitIdx_isLoaded = false;
		simhit_hitType_isLoaded = false;
		bsp_x_isLoaded = false;
		bsp_y_isLoaded = false;
		bsp_z_isLoaded = false;
		bsp_sigmax_isLoaded = false;
		bsp_sigmay_isLoaded = false;
		bsp_sigmaz_isLoaded = false;
		see_fitok_isLoaded = false;
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
		see_statePt_isLoaded = false;
		see_stateTrajX_isLoaded = false;
		see_stateTrajY_isLoaded = false;
		see_stateTrajPx_isLoaded = false;
		see_stateTrajPy_isLoaded = false;
		see_stateTrajPz_isLoaded = false;
		see_stateTrajGlbX_isLoaded = false;
		see_stateTrajGlbY_isLoaded = false;
		see_stateTrajGlbZ_isLoaded = false;
		see_stateTrajGlbPx_isLoaded = false;
		see_stateTrajGlbPy_isLoaded = false;
		see_stateTrajGlbPz_isLoaded = false;
		see_stateCcov00_isLoaded = false;
		see_stateCcov01_isLoaded = false;
		see_stateCcov02_isLoaded = false;
		see_stateCcov03_isLoaded = false;
		see_stateCcov04_isLoaded = false;
		see_stateCcov05_isLoaded = false;
		see_stateCcov11_isLoaded = false;
		see_stateCcov12_isLoaded = false;
		see_stateCcov13_isLoaded = false;
		see_stateCcov14_isLoaded = false;
		see_stateCcov15_isLoaded = false;
		see_stateCcov22_isLoaded = false;
		see_stateCcov23_isLoaded = false;
		see_stateCcov24_isLoaded = false;
		see_stateCcov25_isLoaded = false;
		see_stateCcov33_isLoaded = false;
		see_stateCcov34_isLoaded = false;
		see_stateCcov35_isLoaded = false;
		see_stateCcov44_isLoaded = false;
		see_stateCcov45_isLoaded = false;
		see_stateCcov55_isLoaded = false;
		see_q_isLoaded = false;
		see_nValid_isLoaded = false;
		see_nPixel_isLoaded = false;
		see_nGlued_isLoaded = false;
		see_nStrip_isLoaded = false;
		see_nPhase2OT_isLoaded = false;
		see_algo_isLoaded = false;
		see_stopReason_isLoaded = false;
		see_trkIdx_isLoaded = false;
		see_shareFrac_isLoaded = false;
		see_simTrkIdx_isLoaded = false;
		see_hitIdx_isLoaded = false;
		see_hitType_isLoaded = false;
		see_offset_isLoaded = false;
		vtx_x_isLoaded = false;
		vtx_y_isLoaded = false;
		vtx_z_isLoaded = false;
		vtx_xErr_isLoaded = false;
		vtx_yErr_isLoaded = false;
		vtx_zErr_isLoaded = false;
		vtx_ndof_isLoaded = false;
		vtx_chi2_isLoaded = false;
		vtx_fake_isLoaded = false;
		vtx_valid_isLoaded = false;
		vtx_trkIdx_isLoaded = false;
		simvtx_event_isLoaded = false;
		simvtx_bunchCrossing_isLoaded = false;
		simvtx_processType_isLoaded = false;
		simvtx_x_isLoaded = false;
		simvtx_y_isLoaded = false;
		simvtx_z_isLoaded = false;
		simvtx_sourceSimIdx_isLoaded = false;
		simvtx_daughterSimIdx_isLoaded = false;
		simpv_idx_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (event_branch != 0) event();
	if (lumi_branch != 0) lumi();
	if (run_branch != 0) run();
	if (trk_px_branch != 0) trk_px();
	if (trk_py_branch != 0) trk_py();
	if (trk_pz_branch != 0) trk_pz();
	if (trk_pt_branch != 0) trk_pt();
	if (trk_inner_px_branch != 0) trk_inner_px();
	if (trk_inner_py_branch != 0) trk_inner_py();
	if (trk_inner_pz_branch != 0) trk_inner_pz();
	if (trk_inner_pt_branch != 0) trk_inner_pt();
	if (trk_outer_px_branch != 0) trk_outer_px();
	if (trk_outer_py_branch != 0) trk_outer_py();
	if (trk_outer_pz_branch != 0) trk_outer_pz();
	if (trk_outer_pt_branch != 0) trk_outer_pt();
	if (trk_eta_branch != 0) trk_eta();
	if (trk_lambda_branch != 0) trk_lambda();
	if (trk_cotTheta_branch != 0) trk_cotTheta();
	if (trk_phi_branch != 0) trk_phi();
	if (trk_dxy_branch != 0) trk_dxy();
	if (trk_dz_branch != 0) trk_dz();
	if (trk_ptErr_branch != 0) trk_ptErr();
	if (trk_etaErr_branch != 0) trk_etaErr();
	if (trk_lambdaErr_branch != 0) trk_lambdaErr();
	if (trk_phiErr_branch != 0) trk_phiErr();
	if (trk_dxyErr_branch != 0) trk_dxyErr();
	if (trk_dzErr_branch != 0) trk_dzErr();
	if (trk_refpoint_x_branch != 0) trk_refpoint_x();
	if (trk_refpoint_y_branch != 0) trk_refpoint_y();
	if (trk_refpoint_z_branch != 0) trk_refpoint_z();
	if (trk_nChi2_branch != 0) trk_nChi2();
	if (trk_q_branch != 0) trk_q();
	if (trk_nValid_branch != 0) trk_nValid();
	if (trk_nInvalid_branch != 0) trk_nInvalid();
	if (trk_nPixel_branch != 0) trk_nPixel();
	if (trk_nStrip_branch != 0) trk_nStrip();
	if (trk_nPixelLay_branch != 0) trk_nPixelLay();
	if (trk_nStripLay_branch != 0) trk_nStripLay();
	if (trk_n3DLay_branch != 0) trk_n3DLay();
	if (trk_nOuterLost_branch != 0) trk_nOuterLost();
	if (trk_nInnerLost_branch != 0) trk_nInnerLost();
	if (trk_algo_branch != 0) trk_algo();
	if (trk_originalAlgo_branch != 0) trk_originalAlgo();
	if (trk_algoMask_branch != 0) trk_algoMask();
	if (trk_stopReason_branch != 0) trk_stopReason();
	if (trk_isHP_branch != 0) trk_isHP();
	if (trk_seedIdx_branch != 0) trk_seedIdx();
	if (trk_vtxIdx_branch != 0) trk_vtxIdx();
	if (trk_shareFrac_branch != 0) trk_shareFrac();
	if (trk_simTrkIdx_branch != 0) trk_simTrkIdx();
	if (trk_hitIdx_branch != 0) trk_hitIdx();
	if (trk_hitType_branch != 0) trk_hitType();
	if (sim_event_branch != 0) sim_event();
	if (sim_bunchCrossing_branch != 0) sim_bunchCrossing();
	if (sim_pdgId_branch != 0) sim_pdgId();
	if (sim_genPdgIds_branch != 0) sim_genPdgIds();
	if (sim_isFromBHadron_branch != 0) sim_isFromBHadron();
	if (sim_px_branch != 0) sim_px();
	if (sim_py_branch != 0) sim_py();
	if (sim_pz_branch != 0) sim_pz();
	if (sim_pt_branch != 0) sim_pt();
	if (sim_eta_branch != 0) sim_eta();
	if (sim_phi_branch != 0) sim_phi();
	if (sim_pca_pt_branch != 0) sim_pca_pt();
	if (sim_pca_eta_branch != 0) sim_pca_eta();
	if (sim_pca_lambda_branch != 0) sim_pca_lambda();
	if (sim_pca_cotTheta_branch != 0) sim_pca_cotTheta();
	if (sim_pca_phi_branch != 0) sim_pca_phi();
	if (sim_pca_dxy_branch != 0) sim_pca_dxy();
	if (sim_pca_dz_branch != 0) sim_pca_dz();
	if (sim_q_branch != 0) sim_q();
	if (sim_nValid_branch != 0) sim_nValid();
	if (sim_nPixel_branch != 0) sim_nPixel();
	if (sim_nStrip_branch != 0) sim_nStrip();
	if (sim_nLay_branch != 0) sim_nLay();
	if (sim_nPixelLay_branch != 0) sim_nPixelLay();
	if (sim_n3DLay_branch != 0) sim_n3DLay();
	if (sim_trkIdx_branch != 0) sim_trkIdx();
	if (sim_shareFrac_branch != 0) sim_shareFrac();
	if (sim_seedIdx_branch != 0) sim_seedIdx();
	if (sim_parentVtxIdx_branch != 0) sim_parentVtxIdx();
	if (sim_decayVtxIdx_branch != 0) sim_decayVtxIdx();
	if (sim_simHitIdx_branch != 0) sim_simHitIdx();
	if (pix_isBarrel_branch != 0) pix_isBarrel();
	if (pix_det_branch != 0) pix_det();
	if (pix_lay_branch != 0) pix_lay();
	if (pix_detId_branch != 0) pix_detId();
	if (pix_trkIdx_branch != 0) pix_trkIdx();
	if (pix_seeIdx_branch != 0) pix_seeIdx();
	if (pix_simHitIdx_branch != 0) pix_simHitIdx();
	if (pix_xySignificance_branch != 0) pix_xySignificance();
	if (pix_chargeFraction_branch != 0) pix_chargeFraction();
	if (pix_simType_branch != 0) pix_simType();
	if (pix_x_branch != 0) pix_x();
	if (pix_y_branch != 0) pix_y();
	if (pix_z_branch != 0) pix_z();
	if (pix_xx_branch != 0) pix_xx();
	if (pix_xy_branch != 0) pix_xy();
	if (pix_yy_branch != 0) pix_yy();
	if (pix_yz_branch != 0) pix_yz();
	if (pix_zz_branch != 0) pix_zz();
	if (pix_zx_branch != 0) pix_zx();
	if (pix_radL_branch != 0) pix_radL();
	if (pix_bbxi_branch != 0) pix_bbxi();
	if (ph2_isBarrel_branch != 0) ph2_isBarrel();
	if (ph2_det_branch != 0) ph2_det();
	if (ph2_lay_branch != 0) ph2_lay();
	if (ph2_detId_branch != 0) ph2_detId();
	if (ph2_trkIdx_branch != 0) ph2_trkIdx();
	if (ph2_seeIdx_branch != 0) ph2_seeIdx();
	if (ph2_simHitIdx_branch != 0) ph2_simHitIdx();
	if (ph2_xySignificance_branch != 0) ph2_xySignificance();
	if (ph2_simType_branch != 0) ph2_simType();
	if (ph2_x_branch != 0) ph2_x();
	if (ph2_y_branch != 0) ph2_y();
	if (ph2_z_branch != 0) ph2_z();
	if (ph2_xx_branch != 0) ph2_xx();
	if (ph2_xy_branch != 0) ph2_xy();
	if (ph2_yy_branch != 0) ph2_yy();
	if (ph2_yz_branch != 0) ph2_yz();
	if (ph2_zz_branch != 0) ph2_zz();
	if (ph2_zx_branch != 0) ph2_zx();
	if (ph2_radL_branch != 0) ph2_radL();
	if (ph2_bbxi_branch != 0) ph2_bbxi();
	if (inv_isBarrel_branch != 0) inv_isBarrel();
	if (inv_det_branch != 0) inv_det();
	if (inv_lay_branch != 0) inv_lay();
	if (inv_detId_branch != 0) inv_detId();
	if (inv_type_branch != 0) inv_type();
	if (simhit_det_branch != 0) simhit_det();
	if (simhit_lay_branch != 0) simhit_lay();
	if (simhit_detId_branch != 0) simhit_detId();
	if (simhit_x_branch != 0) simhit_x();
	if (simhit_y_branch != 0) simhit_y();
	if (simhit_z_branch != 0) simhit_z();
	if (simhit_px_branch != 0) simhit_px();
	if (simhit_py_branch != 0) simhit_py();
	if (simhit_pz_branch != 0) simhit_pz();
	if (simhit_particle_branch != 0) simhit_particle();
	if (simhit_process_branch != 0) simhit_process();
	if (simhit_eloss_branch != 0) simhit_eloss();
	if (simhit_tof_branch != 0) simhit_tof();
	if (simhit_simTrkIdx_branch != 0) simhit_simTrkIdx();
	if (simhit_hitIdx_branch != 0) simhit_hitIdx();
	if (simhit_hitType_branch != 0) simhit_hitType();
	if (bsp_x_branch != 0) bsp_x();
	if (bsp_y_branch != 0) bsp_y();
	if (bsp_z_branch != 0) bsp_z();
	if (bsp_sigmax_branch != 0) bsp_sigmax();
	if (bsp_sigmay_branch != 0) bsp_sigmay();
	if (bsp_sigmaz_branch != 0) bsp_sigmaz();
	if (see_fitok_branch != 0) see_fitok();
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
	if (see_statePt_branch != 0) see_statePt();
	if (see_stateTrajX_branch != 0) see_stateTrajX();
	if (see_stateTrajY_branch != 0) see_stateTrajY();
	if (see_stateTrajPx_branch != 0) see_stateTrajPx();
	if (see_stateTrajPy_branch != 0) see_stateTrajPy();
	if (see_stateTrajPz_branch != 0) see_stateTrajPz();
	if (see_stateTrajGlbX_branch != 0) see_stateTrajGlbX();
	if (see_stateTrajGlbY_branch != 0) see_stateTrajGlbY();
	if (see_stateTrajGlbZ_branch != 0) see_stateTrajGlbZ();
	if (see_stateTrajGlbPx_branch != 0) see_stateTrajGlbPx();
	if (see_stateTrajGlbPy_branch != 0) see_stateTrajGlbPy();
	if (see_stateTrajGlbPz_branch != 0) see_stateTrajGlbPz();
	if (see_stateCcov00_branch != 0) see_stateCcov00();
	if (see_stateCcov01_branch != 0) see_stateCcov01();
	if (see_stateCcov02_branch != 0) see_stateCcov02();
	if (see_stateCcov03_branch != 0) see_stateCcov03();
	if (see_stateCcov04_branch != 0) see_stateCcov04();
	if (see_stateCcov05_branch != 0) see_stateCcov05();
	if (see_stateCcov11_branch != 0) see_stateCcov11();
	if (see_stateCcov12_branch != 0) see_stateCcov12();
	if (see_stateCcov13_branch != 0) see_stateCcov13();
	if (see_stateCcov14_branch != 0) see_stateCcov14();
	if (see_stateCcov15_branch != 0) see_stateCcov15();
	if (see_stateCcov22_branch != 0) see_stateCcov22();
	if (see_stateCcov23_branch != 0) see_stateCcov23();
	if (see_stateCcov24_branch != 0) see_stateCcov24();
	if (see_stateCcov25_branch != 0) see_stateCcov25();
	if (see_stateCcov33_branch != 0) see_stateCcov33();
	if (see_stateCcov34_branch != 0) see_stateCcov34();
	if (see_stateCcov35_branch != 0) see_stateCcov35();
	if (see_stateCcov44_branch != 0) see_stateCcov44();
	if (see_stateCcov45_branch != 0) see_stateCcov45();
	if (see_stateCcov55_branch != 0) see_stateCcov55();
	if (see_q_branch != 0) see_q();
	if (see_nValid_branch != 0) see_nValid();
	if (see_nPixel_branch != 0) see_nPixel();
	if (see_nGlued_branch != 0) see_nGlued();
	if (see_nStrip_branch != 0) see_nStrip();
	if (see_nPhase2OT_branch != 0) see_nPhase2OT();
	if (see_algo_branch != 0) see_algo();
	if (see_stopReason_branch != 0) see_stopReason();
	if (see_trkIdx_branch != 0) see_trkIdx();
	if (see_shareFrac_branch != 0) see_shareFrac();
	if (see_simTrkIdx_branch != 0) see_simTrkIdx();
	if (see_hitIdx_branch != 0) see_hitIdx();
	if (see_hitType_branch != 0) see_hitType();
	if (see_offset_branch != 0) see_offset();
	if (vtx_x_branch != 0) vtx_x();
	if (vtx_y_branch != 0) vtx_y();
	if (vtx_z_branch != 0) vtx_z();
	if (vtx_xErr_branch != 0) vtx_xErr();
	if (vtx_yErr_branch != 0) vtx_yErr();
	if (vtx_zErr_branch != 0) vtx_zErr();
	if (vtx_ndof_branch != 0) vtx_ndof();
	if (vtx_chi2_branch != 0) vtx_chi2();
	if (vtx_fake_branch != 0) vtx_fake();
	if (vtx_valid_branch != 0) vtx_valid();
	if (vtx_trkIdx_branch != 0) vtx_trkIdx();
	if (simvtx_event_branch != 0) simvtx_event();
	if (simvtx_bunchCrossing_branch != 0) simvtx_bunchCrossing();
	if (simvtx_processType_branch != 0) simvtx_processType();
	if (simvtx_x_branch != 0) simvtx_x();
	if (simvtx_y_branch != 0) simvtx_y();
	if (simvtx_z_branch != 0) simvtx_z();
	if (simvtx_sourceSimIdx_branch != 0) simvtx_sourceSimIdx();
	if (simvtx_daughterSimIdx_branch != 0) simvtx_daughterSimIdx();
	if (simpv_idx_branch != 0) simpv_idx();
}

	unsigned long long &event()
	{
		if (not event_isLoaded) {
			if (event_branch != 0) {
				event_branch->GetEntry(index);
			} else { 
				printf("branch event_branch does not exist!\n");
				exit(1);
			}
			event_isLoaded = true;
		}
		return event_;
	}
	unsigned int &lumi()
	{
		if (not lumi_isLoaded) {
			if (lumi_branch != 0) {
				lumi_branch->GetEntry(index);
			} else { 
				printf("branch lumi_branch does not exist!\n");
				exit(1);
			}
			lumi_isLoaded = true;
		}
		return lumi_;
	}
	unsigned int &run()
	{
		if (not run_isLoaded) {
			if (run_branch != 0) {
				run_branch->GetEntry(index);
			} else { 
				printf("branch run_branch does not exist!\n");
				exit(1);
			}
			run_isLoaded = true;
		}
		return run_;
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
	vector<float> &trk_inner_px()
	{
		if (not trk_inner_px_isLoaded) {
			if (trk_inner_px_branch != 0) {
				trk_inner_px_branch->GetEntry(index);
			} else { 
				printf("branch trk_inner_px_branch does not exist!\n");
				exit(1);
			}
			trk_inner_px_isLoaded = true;
		}
		return *trk_inner_px_;
	}
	vector<float> &trk_inner_py()
	{
		if (not trk_inner_py_isLoaded) {
			if (trk_inner_py_branch != 0) {
				trk_inner_py_branch->GetEntry(index);
			} else { 
				printf("branch trk_inner_py_branch does not exist!\n");
				exit(1);
			}
			trk_inner_py_isLoaded = true;
		}
		return *trk_inner_py_;
	}
	vector<float> &trk_inner_pz()
	{
		if (not trk_inner_pz_isLoaded) {
			if (trk_inner_pz_branch != 0) {
				trk_inner_pz_branch->GetEntry(index);
			} else { 
				printf("branch trk_inner_pz_branch does not exist!\n");
				exit(1);
			}
			trk_inner_pz_isLoaded = true;
		}
		return *trk_inner_pz_;
	}
	vector<float> &trk_inner_pt()
	{
		if (not trk_inner_pt_isLoaded) {
			if (trk_inner_pt_branch != 0) {
				trk_inner_pt_branch->GetEntry(index);
			} else { 
				printf("branch trk_inner_pt_branch does not exist!\n");
				exit(1);
			}
			trk_inner_pt_isLoaded = true;
		}
		return *trk_inner_pt_;
	}
	vector<float> &trk_outer_px()
	{
		if (not trk_outer_px_isLoaded) {
			if (trk_outer_px_branch != 0) {
				trk_outer_px_branch->GetEntry(index);
			} else { 
				printf("branch trk_outer_px_branch does not exist!\n");
				exit(1);
			}
			trk_outer_px_isLoaded = true;
		}
		return *trk_outer_px_;
	}
	vector<float> &trk_outer_py()
	{
		if (not trk_outer_py_isLoaded) {
			if (trk_outer_py_branch != 0) {
				trk_outer_py_branch->GetEntry(index);
			} else { 
				printf("branch trk_outer_py_branch does not exist!\n");
				exit(1);
			}
			trk_outer_py_isLoaded = true;
		}
		return *trk_outer_py_;
	}
	vector<float> &trk_outer_pz()
	{
		if (not trk_outer_pz_isLoaded) {
			if (trk_outer_pz_branch != 0) {
				trk_outer_pz_branch->GetEntry(index);
			} else { 
				printf("branch trk_outer_pz_branch does not exist!\n");
				exit(1);
			}
			trk_outer_pz_isLoaded = true;
		}
		return *trk_outer_pz_;
	}
	vector<float> &trk_outer_pt()
	{
		if (not trk_outer_pt_isLoaded) {
			if (trk_outer_pt_branch != 0) {
				trk_outer_pt_branch->GetEntry(index);
			} else { 
				printf("branch trk_outer_pt_branch does not exist!\n");
				exit(1);
			}
			trk_outer_pt_isLoaded = true;
		}
		return *trk_outer_pt_;
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
	vector<float> &trk_lambda()
	{
		if (not trk_lambda_isLoaded) {
			if (trk_lambda_branch != 0) {
				trk_lambda_branch->GetEntry(index);
			} else { 
				printf("branch trk_lambda_branch does not exist!\n");
				exit(1);
			}
			trk_lambda_isLoaded = true;
		}
		return *trk_lambda_;
	}
	vector<float> &trk_cotTheta()
	{
		if (not trk_cotTheta_isLoaded) {
			if (trk_cotTheta_branch != 0) {
				trk_cotTheta_branch->GetEntry(index);
			} else { 
				printf("branch trk_cotTheta_branch does not exist!\n");
				exit(1);
			}
			trk_cotTheta_isLoaded = true;
		}
		return *trk_cotTheta_;
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
	vector<float> &trk_lambdaErr()
	{
		if (not trk_lambdaErr_isLoaded) {
			if (trk_lambdaErr_branch != 0) {
				trk_lambdaErr_branch->GetEntry(index);
			} else { 
				printf("branch trk_lambdaErr_branch does not exist!\n");
				exit(1);
			}
			trk_lambdaErr_isLoaded = true;
		}
		return *trk_lambdaErr_;
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
	vector<float> &trk_refpoint_x()
	{
		if (not trk_refpoint_x_isLoaded) {
			if (trk_refpoint_x_branch != 0) {
				trk_refpoint_x_branch->GetEntry(index);
			} else { 
				printf("branch trk_refpoint_x_branch does not exist!\n");
				exit(1);
			}
			trk_refpoint_x_isLoaded = true;
		}
		return *trk_refpoint_x_;
	}
	vector<float> &trk_refpoint_y()
	{
		if (not trk_refpoint_y_isLoaded) {
			if (trk_refpoint_y_branch != 0) {
				trk_refpoint_y_branch->GetEntry(index);
			} else { 
				printf("branch trk_refpoint_y_branch does not exist!\n");
				exit(1);
			}
			trk_refpoint_y_isLoaded = true;
		}
		return *trk_refpoint_y_;
	}
	vector<float> &trk_refpoint_z()
	{
		if (not trk_refpoint_z_isLoaded) {
			if (trk_refpoint_z_branch != 0) {
				trk_refpoint_z_branch->GetEntry(index);
			} else { 
				printf("branch trk_refpoint_z_branch does not exist!\n");
				exit(1);
			}
			trk_refpoint_z_isLoaded = true;
		}
		return *trk_refpoint_z_;
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
	vector<unsigned int> &trk_nValid()
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
	vector<unsigned int> &trk_nInvalid()
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
	vector<unsigned int> &trk_nPixel()
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
	vector<unsigned int> &trk_nStrip()
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
	vector<unsigned int> &trk_nPixelLay()
	{
		if (not trk_nPixelLay_isLoaded) {
			if (trk_nPixelLay_branch != 0) {
				trk_nPixelLay_branch->GetEntry(index);
			} else { 
				printf("branch trk_nPixelLay_branch does not exist!\n");
				exit(1);
			}
			trk_nPixelLay_isLoaded = true;
		}
		return *trk_nPixelLay_;
	}
	vector<unsigned int> &trk_nStripLay()
	{
		if (not trk_nStripLay_isLoaded) {
			if (trk_nStripLay_branch != 0) {
				trk_nStripLay_branch->GetEntry(index);
			} else { 
				printf("branch trk_nStripLay_branch does not exist!\n");
				exit(1);
			}
			trk_nStripLay_isLoaded = true;
		}
		return *trk_nStripLay_;
	}
	vector<unsigned int> &trk_n3DLay()
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
	vector<unsigned int> &trk_nOuterLost()
	{
		if (not trk_nOuterLost_isLoaded) {
			if (trk_nOuterLost_branch != 0) {
				trk_nOuterLost_branch->GetEntry(index);
			} else { 
				printf("branch trk_nOuterLost_branch does not exist!\n");
				exit(1);
			}
			trk_nOuterLost_isLoaded = true;
		}
		return *trk_nOuterLost_;
	}
	vector<unsigned int> &trk_nInnerLost()
	{
		if (not trk_nInnerLost_isLoaded) {
			if (trk_nInnerLost_branch != 0) {
				trk_nInnerLost_branch->GetEntry(index);
			} else { 
				printf("branch trk_nInnerLost_branch does not exist!\n");
				exit(1);
			}
			trk_nInnerLost_isLoaded = true;
		}
		return *trk_nInnerLost_;
	}
	vector<unsigned int> &trk_algo()
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
	vector<unsigned int> &trk_originalAlgo()
	{
		if (not trk_originalAlgo_isLoaded) {
			if (trk_originalAlgo_branch != 0) {
				trk_originalAlgo_branch->GetEntry(index);
			} else { 
				printf("branch trk_originalAlgo_branch does not exist!\n");
				exit(1);
			}
			trk_originalAlgo_isLoaded = true;
		}
		return *trk_originalAlgo_;
	}
	vector<ULong64_t> &trk_algoMask()
	{
		if (not trk_algoMask_isLoaded) {
			if (trk_algoMask_branch != 0) {
				trk_algoMask_branch->GetEntry(index);
			} else { 
				printf("branch trk_algoMask_branch does not exist!\n");
				exit(1);
			}
			trk_algoMask_isLoaded = true;
		}
		return *trk_algoMask_;
	}
	vector<unsigned short> &trk_stopReason()
	{
		if (not trk_stopReason_isLoaded) {
			if (trk_stopReason_branch != 0) {
				trk_stopReason_branch->GetEntry(index);
			} else { 
				printf("branch trk_stopReason_branch does not exist!\n");
				exit(1);
			}
			trk_stopReason_isLoaded = true;
		}
		return *trk_stopReason_;
	}
	vector<short> &trk_isHP()
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
	vector<int> &trk_vtxIdx()
	{
		if (not trk_vtxIdx_isLoaded) {
			if (trk_vtxIdx_branch != 0) {
				trk_vtxIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_vtxIdx_branch does not exist!\n");
				exit(1);
			}
			trk_vtxIdx_isLoaded = true;
		}
		return *trk_vtxIdx_;
	}
	vector<vector<float> > &trk_shareFrac()
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
	vector<vector<int> > &trk_simTrkIdx()
	{
		if (not trk_simTrkIdx_isLoaded) {
			if (trk_simTrkIdx_branch != 0) {
				trk_simTrkIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_simTrkIdx_branch does not exist!\n");
				exit(1);
			}
			trk_simTrkIdx_isLoaded = true;
		}
		return *trk_simTrkIdx_;
	}
	vector<vector<int> > &trk_hitIdx()
	{
		if (not trk_hitIdx_isLoaded) {
			if (trk_hitIdx_branch != 0) {
				trk_hitIdx_branch->GetEntry(index);
			} else { 
				printf("branch trk_hitIdx_branch does not exist!\n");
				exit(1);
			}
			trk_hitIdx_isLoaded = true;
		}
		return *trk_hitIdx_;
	}
	vector<vector<int> > &trk_hitType()
	{
		if (not trk_hitType_isLoaded) {
			if (trk_hitType_branch != 0) {
				trk_hitType_branch->GetEntry(index);
			} else { 
				printf("branch trk_hitType_branch does not exist!\n");
				exit(1);
			}
			trk_hitType_isLoaded = true;
		}
		return *trk_hitType_;
	}
	vector<int> &sim_event()
	{
		if (not sim_event_isLoaded) {
			if (sim_event_branch != 0) {
				sim_event_branch->GetEntry(index);
			} else { 
				printf("branch sim_event_branch does not exist!\n");
				exit(1);
			}
			sim_event_isLoaded = true;
		}
		return *sim_event_;
	}
	vector<int> &sim_bunchCrossing()
	{
		if (not sim_bunchCrossing_isLoaded) {
			if (sim_bunchCrossing_branch != 0) {
				sim_bunchCrossing_branch->GetEntry(index);
			} else { 
				printf("branch sim_bunchCrossing_branch does not exist!\n");
				exit(1);
			}
			sim_bunchCrossing_isLoaded = true;
		}
		return *sim_bunchCrossing_;
	}
	vector<int> &sim_pdgId()
	{
		if (not sim_pdgId_isLoaded) {
			if (sim_pdgId_branch != 0) {
				sim_pdgId_branch->GetEntry(index);
			} else { 
				printf("branch sim_pdgId_branch does not exist!\n");
				exit(1);
			}
			sim_pdgId_isLoaded = true;
		}
		return *sim_pdgId_;
	}
	vector<vector<int> > &sim_genPdgIds()
	{
		if (not sim_genPdgIds_isLoaded) {
			if (sim_genPdgIds_branch != 0) {
				sim_genPdgIds_branch->GetEntry(index);
			} else { 
				printf("branch sim_genPdgIds_branch does not exist!\n");
				exit(1);
			}
			sim_genPdgIds_isLoaded = true;
		}
		return *sim_genPdgIds_;
	}
	vector<int> &sim_isFromBHadron()
	{
		if (not sim_isFromBHadron_isLoaded) {
			if (sim_isFromBHadron_branch != 0) {
				sim_isFromBHadron_branch->GetEntry(index);
			} else { 
				printf("branch sim_isFromBHadron_branch does not exist!\n");
				exit(1);
			}
			sim_isFromBHadron_isLoaded = true;
		}
		return *sim_isFromBHadron_;
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
	vector<float> &sim_pca_pt()
	{
		if (not sim_pca_pt_isLoaded) {
			if (sim_pca_pt_branch != 0) {
				sim_pca_pt_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_pt_branch does not exist!\n");
				exit(1);
			}
			sim_pca_pt_isLoaded = true;
		}
		return *sim_pca_pt_;
	}
	vector<float> &sim_pca_eta()
	{
		if (not sim_pca_eta_isLoaded) {
			if (sim_pca_eta_branch != 0) {
				sim_pca_eta_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_eta_branch does not exist!\n");
				exit(1);
			}
			sim_pca_eta_isLoaded = true;
		}
		return *sim_pca_eta_;
	}
	vector<float> &sim_pca_lambda()
	{
		if (not sim_pca_lambda_isLoaded) {
			if (sim_pca_lambda_branch != 0) {
				sim_pca_lambda_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_lambda_branch does not exist!\n");
				exit(1);
			}
			sim_pca_lambda_isLoaded = true;
		}
		return *sim_pca_lambda_;
	}
	vector<float> &sim_pca_cotTheta()
	{
		if (not sim_pca_cotTheta_isLoaded) {
			if (sim_pca_cotTheta_branch != 0) {
				sim_pca_cotTheta_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_cotTheta_branch does not exist!\n");
				exit(1);
			}
			sim_pca_cotTheta_isLoaded = true;
		}
		return *sim_pca_cotTheta_;
	}
	vector<float> &sim_pca_phi()
	{
		if (not sim_pca_phi_isLoaded) {
			if (sim_pca_phi_branch != 0) {
				sim_pca_phi_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_phi_branch does not exist!\n");
				exit(1);
			}
			sim_pca_phi_isLoaded = true;
		}
		return *sim_pca_phi_;
	}
	vector<float> &sim_pca_dxy()
	{
		if (not sim_pca_dxy_isLoaded) {
			if (sim_pca_dxy_branch != 0) {
				sim_pca_dxy_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_dxy_branch does not exist!\n");
				exit(1);
			}
			sim_pca_dxy_isLoaded = true;
		}
		return *sim_pca_dxy_;
	}
	vector<float> &sim_pca_dz()
	{
		if (not sim_pca_dz_isLoaded) {
			if (sim_pca_dz_branch != 0) {
				sim_pca_dz_branch->GetEntry(index);
			} else { 
				printf("branch sim_pca_dz_branch does not exist!\n");
				exit(1);
			}
			sim_pca_dz_isLoaded = true;
		}
		return *sim_pca_dz_;
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
	vector<unsigned int> &sim_nValid()
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
	vector<unsigned int> &sim_nPixel()
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
	vector<unsigned int> &sim_nStrip()
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
	vector<unsigned int> &sim_nLay()
	{
		if (not sim_nLay_isLoaded) {
			if (sim_nLay_branch != 0) {
				sim_nLay_branch->GetEntry(index);
			} else { 
				printf("branch sim_nLay_branch does not exist!\n");
				exit(1);
			}
			sim_nLay_isLoaded = true;
		}
		return *sim_nLay_;
	}
	vector<unsigned int> &sim_nPixelLay()
	{
		if (not sim_nPixelLay_isLoaded) {
			if (sim_nPixelLay_branch != 0) {
				sim_nPixelLay_branch->GetEntry(index);
			} else { 
				printf("branch sim_nPixelLay_branch does not exist!\n");
				exit(1);
			}
			sim_nPixelLay_isLoaded = true;
		}
		return *sim_nPixelLay_;
	}
	vector<unsigned int> &sim_n3DLay()
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
	vector<vector<int> > &sim_trkIdx()
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
	vector<vector<float> > &sim_shareFrac()
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
	vector<vector<int> > &sim_seedIdx()
	{
		if (not sim_seedIdx_isLoaded) {
			if (sim_seedIdx_branch != 0) {
				sim_seedIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_seedIdx_branch does not exist!\n");
				exit(1);
			}
			sim_seedIdx_isLoaded = true;
		}
		return *sim_seedIdx_;
	}
	vector<int> &sim_parentVtxIdx()
	{
		if (not sim_parentVtxIdx_isLoaded) {
			if (sim_parentVtxIdx_branch != 0) {
				sim_parentVtxIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_parentVtxIdx_branch does not exist!\n");
				exit(1);
			}
			sim_parentVtxIdx_isLoaded = true;
		}
		return *sim_parentVtxIdx_;
	}
	vector<vector<int> > &sim_decayVtxIdx()
	{
		if (not sim_decayVtxIdx_isLoaded) {
			if (sim_decayVtxIdx_branch != 0) {
				sim_decayVtxIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_decayVtxIdx_branch does not exist!\n");
				exit(1);
			}
			sim_decayVtxIdx_isLoaded = true;
		}
		return *sim_decayVtxIdx_;
	}
	vector<vector<int> > &sim_simHitIdx()
	{
		if (not sim_simHitIdx_isLoaded) {
			if (sim_simHitIdx_branch != 0) {
				sim_simHitIdx_branch->GetEntry(index);
			} else { 
				printf("branch sim_simHitIdx_branch does not exist!\n");
				exit(1);
			}
			sim_simHitIdx_isLoaded = true;
		}
		return *sim_simHitIdx_;
	}
	vector<short> &pix_isBarrel()
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
	vector<unsigned short> &pix_det()
	{
		if (not pix_det_isLoaded) {
			if (pix_det_branch != 0) {
				pix_det_branch->GetEntry(index);
			} else { 
				printf("branch pix_det_branch does not exist!\n");
				exit(1);
			}
			pix_det_isLoaded = true;
		}
		return *pix_det_;
	}
	vector<unsigned short> &pix_lay()
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
	vector<unsigned int> &pix_detId()
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
	vector<vector<int> > &pix_trkIdx()
	{
		if (not pix_trkIdx_isLoaded) {
			if (pix_trkIdx_branch != 0) {
				pix_trkIdx_branch->GetEntry(index);
			} else { 
				printf("branch pix_trkIdx_branch does not exist!\n");
				exit(1);
			}
			pix_trkIdx_isLoaded = true;
		}
		return *pix_trkIdx_;
	}
	vector<vector<int> > &pix_seeIdx()
	{
		if (not pix_seeIdx_isLoaded) {
			if (pix_seeIdx_branch != 0) {
				pix_seeIdx_branch->GetEntry(index);
			} else { 
				printf("branch pix_seeIdx_branch does not exist!\n");
				exit(1);
			}
			pix_seeIdx_isLoaded = true;
		}
		return *pix_seeIdx_;
	}
	vector<vector<int> > &pix_simHitIdx()
	{
		if (not pix_simHitIdx_isLoaded) {
			if (pix_simHitIdx_branch != 0) {
				pix_simHitIdx_branch->GetEntry(index);
			} else { 
				printf("branch pix_simHitIdx_branch does not exist!\n");
				exit(1);
			}
			pix_simHitIdx_isLoaded = true;
		}
		return *pix_simHitIdx_;
	}
	vector<vector<float> > &pix_xySignificance()
	{
		if (not pix_xySignificance_isLoaded) {
			if (pix_xySignificance_branch != 0) {
				pix_xySignificance_branch->GetEntry(index);
			} else { 
				printf("branch pix_xySignificance_branch does not exist!\n");
				exit(1);
			}
			pix_xySignificance_isLoaded = true;
		}
		return *pix_xySignificance_;
	}
	vector<vector<float> > &pix_chargeFraction()
	{
		if (not pix_chargeFraction_isLoaded) {
			if (pix_chargeFraction_branch != 0) {
				pix_chargeFraction_branch->GetEntry(index);
			} else { 
				printf("branch pix_chargeFraction_branch does not exist!\n");
				exit(1);
			}
			pix_chargeFraction_isLoaded = true;
		}
		return *pix_chargeFraction_;
	}
	vector<unsigned short> &pix_simType()
	{
		if (not pix_simType_isLoaded) {
			if (pix_simType_branch != 0) {
				pix_simType_branch->GetEntry(index);
			} else { 
				printf("branch pix_simType_branch does not exist!\n");
				exit(1);
			}
			pix_simType_isLoaded = true;
		}
		return *pix_simType_;
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
	vector<short> &ph2_isBarrel()
	{
		if (not ph2_isBarrel_isLoaded) {
			if (ph2_isBarrel_branch != 0) {
				ph2_isBarrel_branch->GetEntry(index);
			} else { 
				printf("branch ph2_isBarrel_branch does not exist!\n");
				exit(1);
			}
			ph2_isBarrel_isLoaded = true;
		}
		return *ph2_isBarrel_;
	}
	vector<unsigned short> &ph2_det()
	{
		if (not ph2_det_isLoaded) {
			if (ph2_det_branch != 0) {
				ph2_det_branch->GetEntry(index);
			} else { 
				printf("branch ph2_det_branch does not exist!\n");
				exit(1);
			}
			ph2_det_isLoaded = true;
		}
		return *ph2_det_;
	}
	vector<unsigned short> &ph2_lay()
	{
		if (not ph2_lay_isLoaded) {
			if (ph2_lay_branch != 0) {
				ph2_lay_branch->GetEntry(index);
			} else { 
				printf("branch ph2_lay_branch does not exist!\n");
				exit(1);
			}
			ph2_lay_isLoaded = true;
		}
		return *ph2_lay_;
	}
	vector<unsigned int> &ph2_detId()
	{
		if (not ph2_detId_isLoaded) {
			if (ph2_detId_branch != 0) {
				ph2_detId_branch->GetEntry(index);
			} else { 
				printf("branch ph2_detId_branch does not exist!\n");
				exit(1);
			}
			ph2_detId_isLoaded = true;
		}
		return *ph2_detId_;
	}
	vector<vector<int> > &ph2_trkIdx()
	{
		if (not ph2_trkIdx_isLoaded) {
			if (ph2_trkIdx_branch != 0) {
				ph2_trkIdx_branch->GetEntry(index);
			} else { 
				printf("branch ph2_trkIdx_branch does not exist!\n");
				exit(1);
			}
			ph2_trkIdx_isLoaded = true;
		}
		return *ph2_trkIdx_;
	}
	vector<vector<int> > &ph2_seeIdx()
	{
		if (not ph2_seeIdx_isLoaded) {
			if (ph2_seeIdx_branch != 0) {
				ph2_seeIdx_branch->GetEntry(index);
			} else { 
				printf("branch ph2_seeIdx_branch does not exist!\n");
				exit(1);
			}
			ph2_seeIdx_isLoaded = true;
		}
		return *ph2_seeIdx_;
	}
	vector<vector<int> > &ph2_simHitIdx()
	{
		if (not ph2_simHitIdx_isLoaded) {
			if (ph2_simHitIdx_branch != 0) {
				ph2_simHitIdx_branch->GetEntry(index);
			} else { 
				printf("branch ph2_simHitIdx_branch does not exist!\n");
				exit(1);
			}
			ph2_simHitIdx_isLoaded = true;
		}
		return *ph2_simHitIdx_;
	}
	vector<vector<float> > &ph2_xySignificance()
	{
		if (not ph2_xySignificance_isLoaded) {
			if (ph2_xySignificance_branch != 0) {
				ph2_xySignificance_branch->GetEntry(index);
			} else { 
				printf("branch ph2_xySignificance_branch does not exist!\n");
				exit(1);
			}
			ph2_xySignificance_isLoaded = true;
		}
		return *ph2_xySignificance_;
	}
	vector<unsigned short> &ph2_simType()
	{
		if (not ph2_simType_isLoaded) {
			if (ph2_simType_branch != 0) {
				ph2_simType_branch->GetEntry(index);
			} else { 
				printf("branch ph2_simType_branch does not exist!\n");
				exit(1);
			}
			ph2_simType_isLoaded = true;
		}
		return *ph2_simType_;
	}
	vector<float> &ph2_x()
	{
		if (not ph2_x_isLoaded) {
			if (ph2_x_branch != 0) {
				ph2_x_branch->GetEntry(index);
			} else { 
				printf("branch ph2_x_branch does not exist!\n");
				exit(1);
			}
			ph2_x_isLoaded = true;
		}
		return *ph2_x_;
	}
	vector<float> &ph2_y()
	{
		if (not ph2_y_isLoaded) {
			if (ph2_y_branch != 0) {
				ph2_y_branch->GetEntry(index);
			} else { 
				printf("branch ph2_y_branch does not exist!\n");
				exit(1);
			}
			ph2_y_isLoaded = true;
		}
		return *ph2_y_;
	}
	vector<float> &ph2_z()
	{
		if (not ph2_z_isLoaded) {
			if (ph2_z_branch != 0) {
				ph2_z_branch->GetEntry(index);
			} else { 
				printf("branch ph2_z_branch does not exist!\n");
				exit(1);
			}
			ph2_z_isLoaded = true;
		}
		return *ph2_z_;
	}
	vector<float> &ph2_xx()
	{
		if (not ph2_xx_isLoaded) {
			if (ph2_xx_branch != 0) {
				ph2_xx_branch->GetEntry(index);
			} else { 
				printf("branch ph2_xx_branch does not exist!\n");
				exit(1);
			}
			ph2_xx_isLoaded = true;
		}
		return *ph2_xx_;
	}
	vector<float> &ph2_xy()
	{
		if (not ph2_xy_isLoaded) {
			if (ph2_xy_branch != 0) {
				ph2_xy_branch->GetEntry(index);
			} else { 
				printf("branch ph2_xy_branch does not exist!\n");
				exit(1);
			}
			ph2_xy_isLoaded = true;
		}
		return *ph2_xy_;
	}
	vector<float> &ph2_yy()
	{
		if (not ph2_yy_isLoaded) {
			if (ph2_yy_branch != 0) {
				ph2_yy_branch->GetEntry(index);
			} else { 
				printf("branch ph2_yy_branch does not exist!\n");
				exit(1);
			}
			ph2_yy_isLoaded = true;
		}
		return *ph2_yy_;
	}
	vector<float> &ph2_yz()
	{
		if (not ph2_yz_isLoaded) {
			if (ph2_yz_branch != 0) {
				ph2_yz_branch->GetEntry(index);
			} else { 
				printf("branch ph2_yz_branch does not exist!\n");
				exit(1);
			}
			ph2_yz_isLoaded = true;
		}
		return *ph2_yz_;
	}
	vector<float> &ph2_zz()
	{
		if (not ph2_zz_isLoaded) {
			if (ph2_zz_branch != 0) {
				ph2_zz_branch->GetEntry(index);
			} else { 
				printf("branch ph2_zz_branch does not exist!\n");
				exit(1);
			}
			ph2_zz_isLoaded = true;
		}
		return *ph2_zz_;
	}
	vector<float> &ph2_zx()
	{
		if (not ph2_zx_isLoaded) {
			if (ph2_zx_branch != 0) {
				ph2_zx_branch->GetEntry(index);
			} else { 
				printf("branch ph2_zx_branch does not exist!\n");
				exit(1);
			}
			ph2_zx_isLoaded = true;
		}
		return *ph2_zx_;
	}
	vector<float> &ph2_radL()
	{
		if (not ph2_radL_isLoaded) {
			if (ph2_radL_branch != 0) {
				ph2_radL_branch->GetEntry(index);
			} else { 
				printf("branch ph2_radL_branch does not exist!\n");
				exit(1);
			}
			ph2_radL_isLoaded = true;
		}
		return *ph2_radL_;
	}
	vector<float> &ph2_bbxi()
	{
		if (not ph2_bbxi_isLoaded) {
			if (ph2_bbxi_branch != 0) {
				ph2_bbxi_branch->GetEntry(index);
			} else { 
				printf("branch ph2_bbxi_branch does not exist!\n");
				exit(1);
			}
			ph2_bbxi_isLoaded = true;
		}
		return *ph2_bbxi_;
	}
	vector<short> &inv_isBarrel()
	{
		if (not inv_isBarrel_isLoaded) {
			if (inv_isBarrel_branch != 0) {
				inv_isBarrel_branch->GetEntry(index);
			} else { 
				printf("branch inv_isBarrel_branch does not exist!\n");
				exit(1);
			}
			inv_isBarrel_isLoaded = true;
		}
		return *inv_isBarrel_;
	}
	vector<unsigned short> &inv_det()
	{
		if (not inv_det_isLoaded) {
			if (inv_det_branch != 0) {
				inv_det_branch->GetEntry(index);
			} else { 
				printf("branch inv_det_branch does not exist!\n");
				exit(1);
			}
			inv_det_isLoaded = true;
		}
		return *inv_det_;
	}
	vector<unsigned short> &inv_lay()
	{
		if (not inv_lay_isLoaded) {
			if (inv_lay_branch != 0) {
				inv_lay_branch->GetEntry(index);
			} else { 
				printf("branch inv_lay_branch does not exist!\n");
				exit(1);
			}
			inv_lay_isLoaded = true;
		}
		return *inv_lay_;
	}
	vector<unsigned int> &inv_detId()
	{
		if (not inv_detId_isLoaded) {
			if (inv_detId_branch != 0) {
				inv_detId_branch->GetEntry(index);
			} else { 
				printf("branch inv_detId_branch does not exist!\n");
				exit(1);
			}
			inv_detId_isLoaded = true;
		}
		return *inv_detId_;
	}
	vector<unsigned short> &inv_type()
	{
		if (not inv_type_isLoaded) {
			if (inv_type_branch != 0) {
				inv_type_branch->GetEntry(index);
			} else { 
				printf("branch inv_type_branch does not exist!\n");
				exit(1);
			}
			inv_type_isLoaded = true;
		}
		return *inv_type_;
	}
	vector<unsigned short> &simhit_det()
	{
		if (not simhit_det_isLoaded) {
			if (simhit_det_branch != 0) {
				simhit_det_branch->GetEntry(index);
			} else { 
				printf("branch simhit_det_branch does not exist!\n");
				exit(1);
			}
			simhit_det_isLoaded = true;
		}
		return *simhit_det_;
	}
	vector<unsigned short> &simhit_lay()
	{
		if (not simhit_lay_isLoaded) {
			if (simhit_lay_branch != 0) {
				simhit_lay_branch->GetEntry(index);
			} else { 
				printf("branch simhit_lay_branch does not exist!\n");
				exit(1);
			}
			simhit_lay_isLoaded = true;
		}
		return *simhit_lay_;
	}
	vector<unsigned int> &simhit_detId()
	{
		if (not simhit_detId_isLoaded) {
			if (simhit_detId_branch != 0) {
				simhit_detId_branch->GetEntry(index);
			} else { 
				printf("branch simhit_detId_branch does not exist!\n");
				exit(1);
			}
			simhit_detId_isLoaded = true;
		}
		return *simhit_detId_;
	}
	vector<float> &simhit_x()
	{
		if (not simhit_x_isLoaded) {
			if (simhit_x_branch != 0) {
				simhit_x_branch->GetEntry(index);
			} else { 
				printf("branch simhit_x_branch does not exist!\n");
				exit(1);
			}
			simhit_x_isLoaded = true;
		}
		return *simhit_x_;
	}
	vector<float> &simhit_y()
	{
		if (not simhit_y_isLoaded) {
			if (simhit_y_branch != 0) {
				simhit_y_branch->GetEntry(index);
			} else { 
				printf("branch simhit_y_branch does not exist!\n");
				exit(1);
			}
			simhit_y_isLoaded = true;
		}
		return *simhit_y_;
	}
	vector<float> &simhit_z()
	{
		if (not simhit_z_isLoaded) {
			if (simhit_z_branch != 0) {
				simhit_z_branch->GetEntry(index);
			} else { 
				printf("branch simhit_z_branch does not exist!\n");
				exit(1);
			}
			simhit_z_isLoaded = true;
		}
		return *simhit_z_;
	}
	vector<float> &simhit_px()
	{
		if (not simhit_px_isLoaded) {
			if (simhit_px_branch != 0) {
				simhit_px_branch->GetEntry(index);
			} else { 
				printf("branch simhit_px_branch does not exist!\n");
				exit(1);
			}
			simhit_px_isLoaded = true;
		}
		return *simhit_px_;
	}
	vector<float> &simhit_py()
	{
		if (not simhit_py_isLoaded) {
			if (simhit_py_branch != 0) {
				simhit_py_branch->GetEntry(index);
			} else { 
				printf("branch simhit_py_branch does not exist!\n");
				exit(1);
			}
			simhit_py_isLoaded = true;
		}
		return *simhit_py_;
	}
	vector<float> &simhit_pz()
	{
		if (not simhit_pz_isLoaded) {
			if (simhit_pz_branch != 0) {
				simhit_pz_branch->GetEntry(index);
			} else { 
				printf("branch simhit_pz_branch does not exist!\n");
				exit(1);
			}
			simhit_pz_isLoaded = true;
		}
		return *simhit_pz_;
	}
	vector<int> &simhit_particle()
	{
		if (not simhit_particle_isLoaded) {
			if (simhit_particle_branch != 0) {
				simhit_particle_branch->GetEntry(index);
			} else { 
				printf("branch simhit_particle_branch does not exist!\n");
				exit(1);
			}
			simhit_particle_isLoaded = true;
		}
		return *simhit_particle_;
	}
	vector<short> &simhit_process()
	{
		if (not simhit_process_isLoaded) {
			if (simhit_process_branch != 0) {
				simhit_process_branch->GetEntry(index);
			} else { 
				printf("branch simhit_process_branch does not exist!\n");
				exit(1);
			}
			simhit_process_isLoaded = true;
		}
		return *simhit_process_;
	}
	vector<float> &simhit_eloss()
	{
		if (not simhit_eloss_isLoaded) {
			if (simhit_eloss_branch != 0) {
				simhit_eloss_branch->GetEntry(index);
			} else { 
				printf("branch simhit_eloss_branch does not exist!\n");
				exit(1);
			}
			simhit_eloss_isLoaded = true;
		}
		return *simhit_eloss_;
	}
	vector<float> &simhit_tof()
	{
		if (not simhit_tof_isLoaded) {
			if (simhit_tof_branch != 0) {
				simhit_tof_branch->GetEntry(index);
			} else { 
				printf("branch simhit_tof_branch does not exist!\n");
				exit(1);
			}
			simhit_tof_isLoaded = true;
		}
		return *simhit_tof_;
	}
	vector<int> &simhit_simTrkIdx()
	{
		if (not simhit_simTrkIdx_isLoaded) {
			if (simhit_simTrkIdx_branch != 0) {
				simhit_simTrkIdx_branch->GetEntry(index);
			} else { 
				printf("branch simhit_simTrkIdx_branch does not exist!\n");
				exit(1);
			}
			simhit_simTrkIdx_isLoaded = true;
		}
		return *simhit_simTrkIdx_;
	}
	vector<vector<int> > &simhit_hitIdx()
	{
		if (not simhit_hitIdx_isLoaded) {
			if (simhit_hitIdx_branch != 0) {
				simhit_hitIdx_branch->GetEntry(index);
			} else { 
				printf("branch simhit_hitIdx_branch does not exist!\n");
				exit(1);
			}
			simhit_hitIdx_isLoaded = true;
		}
		return *simhit_hitIdx_;
	}
	vector<vector<int> > &simhit_hitType()
	{
		if (not simhit_hitType_isLoaded) {
			if (simhit_hitType_branch != 0) {
				simhit_hitType_branch->GetEntry(index);
			} else { 
				printf("branch simhit_hitType_branch does not exist!\n");
				exit(1);
			}
			simhit_hitType_isLoaded = true;
		}
		return *simhit_hitType_;
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
	vector<short> &see_fitok()
	{
		if (not see_fitok_isLoaded) {
			if (see_fitok_branch != 0) {
				see_fitok_branch->GetEntry(index);
			} else { 
				printf("branch see_fitok_branch does not exist!\n");
				exit(1);
			}
			see_fitok_isLoaded = true;
		}
		return *see_fitok_;
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
	vector<float> &see_statePt()
	{
		if (not see_statePt_isLoaded) {
			if (see_statePt_branch != 0) {
				see_statePt_branch->GetEntry(index);
			} else { 
				printf("branch see_statePt_branch does not exist!\n");
				exit(1);
			}
			see_statePt_isLoaded = true;
		}
		return *see_statePt_;
	}
	vector<float> &see_stateTrajX()
	{
		if (not see_stateTrajX_isLoaded) {
			if (see_stateTrajX_branch != 0) {
				see_stateTrajX_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajX_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajX_isLoaded = true;
		}
		return *see_stateTrajX_;
	}
	vector<float> &see_stateTrajY()
	{
		if (not see_stateTrajY_isLoaded) {
			if (see_stateTrajY_branch != 0) {
				see_stateTrajY_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajY_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajY_isLoaded = true;
		}
		return *see_stateTrajY_;
	}
	vector<float> &see_stateTrajPx()
	{
		if (not see_stateTrajPx_isLoaded) {
			if (see_stateTrajPx_branch != 0) {
				see_stateTrajPx_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajPx_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajPx_isLoaded = true;
		}
		return *see_stateTrajPx_;
	}
	vector<float> &see_stateTrajPy()
	{
		if (not see_stateTrajPy_isLoaded) {
			if (see_stateTrajPy_branch != 0) {
				see_stateTrajPy_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajPy_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajPy_isLoaded = true;
		}
		return *see_stateTrajPy_;
	}
	vector<float> &see_stateTrajPz()
	{
		if (not see_stateTrajPz_isLoaded) {
			if (see_stateTrajPz_branch != 0) {
				see_stateTrajPz_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajPz_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajPz_isLoaded = true;
		}
		return *see_stateTrajPz_;
	}
	vector<float> &see_stateTrajGlbX()
	{
		if (not see_stateTrajGlbX_isLoaded) {
			if (see_stateTrajGlbX_branch != 0) {
				see_stateTrajGlbX_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajGlbX_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajGlbX_isLoaded = true;
		}
		return *see_stateTrajGlbX_;
	}
	vector<float> &see_stateTrajGlbY()
	{
		if (not see_stateTrajGlbY_isLoaded) {
			if (see_stateTrajGlbY_branch != 0) {
				see_stateTrajGlbY_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajGlbY_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajGlbY_isLoaded = true;
		}
		return *see_stateTrajGlbY_;
	}
	vector<float> &see_stateTrajGlbZ()
	{
		if (not see_stateTrajGlbZ_isLoaded) {
			if (see_stateTrajGlbZ_branch != 0) {
				see_stateTrajGlbZ_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajGlbZ_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajGlbZ_isLoaded = true;
		}
		return *see_stateTrajGlbZ_;
	}
	vector<float> &see_stateTrajGlbPx()
	{
		if (not see_stateTrajGlbPx_isLoaded) {
			if (see_stateTrajGlbPx_branch != 0) {
				see_stateTrajGlbPx_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajGlbPx_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajGlbPx_isLoaded = true;
		}
		return *see_stateTrajGlbPx_;
	}
	vector<float> &see_stateTrajGlbPy()
	{
		if (not see_stateTrajGlbPy_isLoaded) {
			if (see_stateTrajGlbPy_branch != 0) {
				see_stateTrajGlbPy_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajGlbPy_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajGlbPy_isLoaded = true;
		}
		return *see_stateTrajGlbPy_;
	}
	vector<float> &see_stateTrajGlbPz()
	{
		if (not see_stateTrajGlbPz_isLoaded) {
			if (see_stateTrajGlbPz_branch != 0) {
				see_stateTrajGlbPz_branch->GetEntry(index);
			} else { 
				printf("branch see_stateTrajGlbPz_branch does not exist!\n");
				exit(1);
			}
			see_stateTrajGlbPz_isLoaded = true;
		}
		return *see_stateTrajGlbPz_;
	}
	vector<float> &see_stateCcov00()
	{
		if (not see_stateCcov00_isLoaded) {
			if (see_stateCcov00_branch != 0) {
				see_stateCcov00_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov00_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov00_isLoaded = true;
		}
		return *see_stateCcov00_;
	}
	vector<float> &see_stateCcov01()
	{
		if (not see_stateCcov01_isLoaded) {
			if (see_stateCcov01_branch != 0) {
				see_stateCcov01_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov01_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov01_isLoaded = true;
		}
		return *see_stateCcov01_;
	}
	vector<float> &see_stateCcov02()
	{
		if (not see_stateCcov02_isLoaded) {
			if (see_stateCcov02_branch != 0) {
				see_stateCcov02_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov02_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov02_isLoaded = true;
		}
		return *see_stateCcov02_;
	}
	vector<float> &see_stateCcov03()
	{
		if (not see_stateCcov03_isLoaded) {
			if (see_stateCcov03_branch != 0) {
				see_stateCcov03_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov03_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov03_isLoaded = true;
		}
		return *see_stateCcov03_;
	}
	vector<float> &see_stateCcov04()
	{
		if (not see_stateCcov04_isLoaded) {
			if (see_stateCcov04_branch != 0) {
				see_stateCcov04_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov04_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov04_isLoaded = true;
		}
		return *see_stateCcov04_;
	}
	vector<float> &see_stateCcov05()
	{
		if (not see_stateCcov05_isLoaded) {
			if (see_stateCcov05_branch != 0) {
				see_stateCcov05_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov05_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov05_isLoaded = true;
		}
		return *see_stateCcov05_;
	}
	vector<float> &see_stateCcov11()
	{
		if (not see_stateCcov11_isLoaded) {
			if (see_stateCcov11_branch != 0) {
				see_stateCcov11_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov11_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov11_isLoaded = true;
		}
		return *see_stateCcov11_;
	}
	vector<float> &see_stateCcov12()
	{
		if (not see_stateCcov12_isLoaded) {
			if (see_stateCcov12_branch != 0) {
				see_stateCcov12_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov12_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov12_isLoaded = true;
		}
		return *see_stateCcov12_;
	}
	vector<float> &see_stateCcov13()
	{
		if (not see_stateCcov13_isLoaded) {
			if (see_stateCcov13_branch != 0) {
				see_stateCcov13_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov13_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov13_isLoaded = true;
		}
		return *see_stateCcov13_;
	}
	vector<float> &see_stateCcov14()
	{
		if (not see_stateCcov14_isLoaded) {
			if (see_stateCcov14_branch != 0) {
				see_stateCcov14_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov14_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov14_isLoaded = true;
		}
		return *see_stateCcov14_;
	}
	vector<float> &see_stateCcov15()
	{
		if (not see_stateCcov15_isLoaded) {
			if (see_stateCcov15_branch != 0) {
				see_stateCcov15_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov15_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov15_isLoaded = true;
		}
		return *see_stateCcov15_;
	}
	vector<float> &see_stateCcov22()
	{
		if (not see_stateCcov22_isLoaded) {
			if (see_stateCcov22_branch != 0) {
				see_stateCcov22_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov22_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov22_isLoaded = true;
		}
		return *see_stateCcov22_;
	}
	vector<float> &see_stateCcov23()
	{
		if (not see_stateCcov23_isLoaded) {
			if (see_stateCcov23_branch != 0) {
				see_stateCcov23_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov23_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov23_isLoaded = true;
		}
		return *see_stateCcov23_;
	}
	vector<float> &see_stateCcov24()
	{
		if (not see_stateCcov24_isLoaded) {
			if (see_stateCcov24_branch != 0) {
				see_stateCcov24_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov24_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov24_isLoaded = true;
		}
		return *see_stateCcov24_;
	}
	vector<float> &see_stateCcov25()
	{
		if (not see_stateCcov25_isLoaded) {
			if (see_stateCcov25_branch != 0) {
				see_stateCcov25_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov25_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov25_isLoaded = true;
		}
		return *see_stateCcov25_;
	}
	vector<float> &see_stateCcov33()
	{
		if (not see_stateCcov33_isLoaded) {
			if (see_stateCcov33_branch != 0) {
				see_stateCcov33_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov33_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov33_isLoaded = true;
		}
		return *see_stateCcov33_;
	}
	vector<float> &see_stateCcov34()
	{
		if (not see_stateCcov34_isLoaded) {
			if (see_stateCcov34_branch != 0) {
				see_stateCcov34_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov34_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov34_isLoaded = true;
		}
		return *see_stateCcov34_;
	}
	vector<float> &see_stateCcov35()
	{
		if (not see_stateCcov35_isLoaded) {
			if (see_stateCcov35_branch != 0) {
				see_stateCcov35_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov35_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov35_isLoaded = true;
		}
		return *see_stateCcov35_;
	}
	vector<float> &see_stateCcov44()
	{
		if (not see_stateCcov44_isLoaded) {
			if (see_stateCcov44_branch != 0) {
				see_stateCcov44_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov44_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov44_isLoaded = true;
		}
		return *see_stateCcov44_;
	}
	vector<float> &see_stateCcov45()
	{
		if (not see_stateCcov45_isLoaded) {
			if (see_stateCcov45_branch != 0) {
				see_stateCcov45_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov45_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov45_isLoaded = true;
		}
		return *see_stateCcov45_;
	}
	vector<float> &see_stateCcov55()
	{
		if (not see_stateCcov55_isLoaded) {
			if (see_stateCcov55_branch != 0) {
				see_stateCcov55_branch->GetEntry(index);
			} else { 
				printf("branch see_stateCcov55_branch does not exist!\n");
				exit(1);
			}
			see_stateCcov55_isLoaded = true;
		}
		return *see_stateCcov55_;
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
	vector<unsigned int> &see_nValid()
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
	vector<unsigned int> &see_nPixel()
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
	vector<unsigned int> &see_nGlued()
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
	vector<unsigned int> &see_nStrip()
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
	vector<unsigned int> &see_nPhase2OT()
	{
		if (not see_nPhase2OT_isLoaded) {
			if (see_nPhase2OT_branch != 0) {
				see_nPhase2OT_branch->GetEntry(index);
			} else { 
				printf("branch see_nPhase2OT_branch does not exist!\n");
				exit(1);
			}
			see_nPhase2OT_isLoaded = true;
		}
		return *see_nPhase2OT_;
	}
	vector<unsigned int> &see_algo()
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
	vector<unsigned short> &see_stopReason()
	{
		if (not see_stopReason_isLoaded) {
			if (see_stopReason_branch != 0) {
				see_stopReason_branch->GetEntry(index);
			} else { 
				printf("branch see_stopReason_branch does not exist!\n");
				exit(1);
			}
			see_stopReason_isLoaded = true;
		}
		return *see_stopReason_;
	}
	vector<int> &see_trkIdx()
	{
		if (not see_trkIdx_isLoaded) {
			if (see_trkIdx_branch != 0) {
				see_trkIdx_branch->GetEntry(index);
			} else { 
				printf("branch see_trkIdx_branch does not exist!\n");
				exit(1);
			}
			see_trkIdx_isLoaded = true;
		}
		return *see_trkIdx_;
	}
	vector<vector<float> > &see_shareFrac()
	{
		if (not see_shareFrac_isLoaded) {
			if (see_shareFrac_branch != 0) {
				see_shareFrac_branch->GetEntry(index);
			} else { 
				printf("branch see_shareFrac_branch does not exist!\n");
				exit(1);
			}
			see_shareFrac_isLoaded = true;
		}
		return *see_shareFrac_;
	}
	vector<vector<int> > &see_simTrkIdx()
	{
		if (not see_simTrkIdx_isLoaded) {
			if (see_simTrkIdx_branch != 0) {
				see_simTrkIdx_branch->GetEntry(index);
			} else { 
				printf("branch see_simTrkIdx_branch does not exist!\n");
				exit(1);
			}
			see_simTrkIdx_isLoaded = true;
		}
		return *see_simTrkIdx_;
	}
	vector<vector<int> > &see_hitIdx()
	{
		if (not see_hitIdx_isLoaded) {
			if (see_hitIdx_branch != 0) {
				see_hitIdx_branch->GetEntry(index);
			} else { 
				printf("branch see_hitIdx_branch does not exist!\n");
				exit(1);
			}
			see_hitIdx_isLoaded = true;
		}
		return *see_hitIdx_;
	}
	vector<vector<int> > &see_hitType()
	{
		if (not see_hitType_isLoaded) {
			if (see_hitType_branch != 0) {
				see_hitType_branch->GetEntry(index);
			} else { 
				printf("branch see_hitType_branch does not exist!\n");
				exit(1);
			}
			see_hitType_isLoaded = true;
		}
		return *see_hitType_;
	}
	vector<unsigned int> &see_offset()
	{
		if (not see_offset_isLoaded) {
			if (see_offset_branch != 0) {
				see_offset_branch->GetEntry(index);
			} else { 
				printf("branch see_offset_branch does not exist!\n");
				exit(1);
			}
			see_offset_isLoaded = true;
		}
		return *see_offset_;
	}
	vector<float> &vtx_x()
	{
		if (not vtx_x_isLoaded) {
			if (vtx_x_branch != 0) {
				vtx_x_branch->GetEntry(index);
			} else { 
				printf("branch vtx_x_branch does not exist!\n");
				exit(1);
			}
			vtx_x_isLoaded = true;
		}
		return *vtx_x_;
	}
	vector<float> &vtx_y()
	{
		if (not vtx_y_isLoaded) {
			if (vtx_y_branch != 0) {
				vtx_y_branch->GetEntry(index);
			} else { 
				printf("branch vtx_y_branch does not exist!\n");
				exit(1);
			}
			vtx_y_isLoaded = true;
		}
		return *vtx_y_;
	}
	vector<float> &vtx_z()
	{
		if (not vtx_z_isLoaded) {
			if (vtx_z_branch != 0) {
				vtx_z_branch->GetEntry(index);
			} else { 
				printf("branch vtx_z_branch does not exist!\n");
				exit(1);
			}
			vtx_z_isLoaded = true;
		}
		return *vtx_z_;
	}
	vector<float> &vtx_xErr()
	{
		if (not vtx_xErr_isLoaded) {
			if (vtx_xErr_branch != 0) {
				vtx_xErr_branch->GetEntry(index);
			} else { 
				printf("branch vtx_xErr_branch does not exist!\n");
				exit(1);
			}
			vtx_xErr_isLoaded = true;
		}
		return *vtx_xErr_;
	}
	vector<float> &vtx_yErr()
	{
		if (not vtx_yErr_isLoaded) {
			if (vtx_yErr_branch != 0) {
				vtx_yErr_branch->GetEntry(index);
			} else { 
				printf("branch vtx_yErr_branch does not exist!\n");
				exit(1);
			}
			vtx_yErr_isLoaded = true;
		}
		return *vtx_yErr_;
	}
	vector<float> &vtx_zErr()
	{
		if (not vtx_zErr_isLoaded) {
			if (vtx_zErr_branch != 0) {
				vtx_zErr_branch->GetEntry(index);
			} else { 
				printf("branch vtx_zErr_branch does not exist!\n");
				exit(1);
			}
			vtx_zErr_isLoaded = true;
		}
		return *vtx_zErr_;
	}
	vector<float> &vtx_ndof()
	{
		if (not vtx_ndof_isLoaded) {
			if (vtx_ndof_branch != 0) {
				vtx_ndof_branch->GetEntry(index);
			} else { 
				printf("branch vtx_ndof_branch does not exist!\n");
				exit(1);
			}
			vtx_ndof_isLoaded = true;
		}
		return *vtx_ndof_;
	}
	vector<float> &vtx_chi2()
	{
		if (not vtx_chi2_isLoaded) {
			if (vtx_chi2_branch != 0) {
				vtx_chi2_branch->GetEntry(index);
			} else { 
				printf("branch vtx_chi2_branch does not exist!\n");
				exit(1);
			}
			vtx_chi2_isLoaded = true;
		}
		return *vtx_chi2_;
	}
	vector<short> &vtx_fake()
	{
		if (not vtx_fake_isLoaded) {
			if (vtx_fake_branch != 0) {
				vtx_fake_branch->GetEntry(index);
			} else { 
				printf("branch vtx_fake_branch does not exist!\n");
				exit(1);
			}
			vtx_fake_isLoaded = true;
		}
		return *vtx_fake_;
	}
	vector<short> &vtx_valid()
	{
		if (not vtx_valid_isLoaded) {
			if (vtx_valid_branch != 0) {
				vtx_valid_branch->GetEntry(index);
			} else { 
				printf("branch vtx_valid_branch does not exist!\n");
				exit(1);
			}
			vtx_valid_isLoaded = true;
		}
		return *vtx_valid_;
	}
	vector<vector<int> > &vtx_trkIdx()
	{
		if (not vtx_trkIdx_isLoaded) {
			if (vtx_trkIdx_branch != 0) {
				vtx_trkIdx_branch->GetEntry(index);
			} else { 
				printf("branch vtx_trkIdx_branch does not exist!\n");
				exit(1);
			}
			vtx_trkIdx_isLoaded = true;
		}
		return *vtx_trkIdx_;
	}
	vector<int> &simvtx_event()
	{
		if (not simvtx_event_isLoaded) {
			if (simvtx_event_branch != 0) {
				simvtx_event_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_event_branch does not exist!\n");
				exit(1);
			}
			simvtx_event_isLoaded = true;
		}
		return *simvtx_event_;
	}
	vector<int> &simvtx_bunchCrossing()
	{
		if (not simvtx_bunchCrossing_isLoaded) {
			if (simvtx_bunchCrossing_branch != 0) {
				simvtx_bunchCrossing_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_bunchCrossing_branch does not exist!\n");
				exit(1);
			}
			simvtx_bunchCrossing_isLoaded = true;
		}
		return *simvtx_bunchCrossing_;
	}
	vector<unsigned int> &simvtx_processType()
	{
		if (not simvtx_processType_isLoaded) {
			if (simvtx_processType_branch != 0) {
				simvtx_processType_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_processType_branch does not exist!\n");
				exit(1);
			}
			simvtx_processType_isLoaded = true;
		}
		return *simvtx_processType_;
	}
	vector<float> &simvtx_x()
	{
		if (not simvtx_x_isLoaded) {
			if (simvtx_x_branch != 0) {
				simvtx_x_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_x_branch does not exist!\n");
				exit(1);
			}
			simvtx_x_isLoaded = true;
		}
		return *simvtx_x_;
	}
	vector<float> &simvtx_y()
	{
		if (not simvtx_y_isLoaded) {
			if (simvtx_y_branch != 0) {
				simvtx_y_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_y_branch does not exist!\n");
				exit(1);
			}
			simvtx_y_isLoaded = true;
		}
		return *simvtx_y_;
	}
	vector<float> &simvtx_z()
	{
		if (not simvtx_z_isLoaded) {
			if (simvtx_z_branch != 0) {
				simvtx_z_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_z_branch does not exist!\n");
				exit(1);
			}
			simvtx_z_isLoaded = true;
		}
		return *simvtx_z_;
	}
	vector<vector<int> > &simvtx_sourceSimIdx()
	{
		if (not simvtx_sourceSimIdx_isLoaded) {
			if (simvtx_sourceSimIdx_branch != 0) {
				simvtx_sourceSimIdx_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_sourceSimIdx_branch does not exist!\n");
				exit(1);
			}
			simvtx_sourceSimIdx_isLoaded = true;
		}
		return *simvtx_sourceSimIdx_;
	}
	vector<vector<int> > &simvtx_daughterSimIdx()
	{
		if (not simvtx_daughterSimIdx_isLoaded) {
			if (simvtx_daughterSimIdx_branch != 0) {
				simvtx_daughterSimIdx_branch->GetEntry(index);
			} else { 
				printf("branch simvtx_daughterSimIdx_branch does not exist!\n");
				exit(1);
			}
			simvtx_daughterSimIdx_isLoaded = true;
		}
		return *simvtx_daughterSimIdx_;
	}
	vector<int> &simpv_idx()
	{
		if (not simpv_idx_isLoaded) {
			if (simpv_idx_branch != 0) {
				simpv_idx_branch->GetEntry(index);
			} else { 
				printf("branch simpv_idx_branch does not exist!\n");
				exit(1);
			}
			simpv_idx_isLoaded = true;
		}
		return *simpv_idx_;
	}
};

#ifndef __CINT__
extern tkph2 cms2;
#endif

namespace tas {
	unsigned long long &event();
	unsigned int &lumi();
	unsigned int &run();
	vector<float> &trk_px();
	vector<float> &trk_py();
	vector<float> &trk_pz();
	vector<float> &trk_pt();
	vector<float> &trk_inner_px();
	vector<float> &trk_inner_py();
	vector<float> &trk_inner_pz();
	vector<float> &trk_inner_pt();
	vector<float> &trk_outer_px();
	vector<float> &trk_outer_py();
	vector<float> &trk_outer_pz();
	vector<float> &trk_outer_pt();
	vector<float> &trk_eta();
	vector<float> &trk_lambda();
	vector<float> &trk_cotTheta();
	vector<float> &trk_phi();
	vector<float> &trk_dxy();
	vector<float> &trk_dz();
	vector<float> &trk_ptErr();
	vector<float> &trk_etaErr();
	vector<float> &trk_lambdaErr();
	vector<float> &trk_phiErr();
	vector<float> &trk_dxyErr();
	vector<float> &trk_dzErr();
	vector<float> &trk_refpoint_x();
	vector<float> &trk_refpoint_y();
	vector<float> &trk_refpoint_z();
	vector<float> &trk_nChi2();
	vector<int> &trk_q();
	vector<unsigned int> &trk_nValid();
	vector<unsigned int> &trk_nInvalid();
	vector<unsigned int> &trk_nPixel();
	vector<unsigned int> &trk_nStrip();
	vector<unsigned int> &trk_nPixelLay();
	vector<unsigned int> &trk_nStripLay();
	vector<unsigned int> &trk_n3DLay();
	vector<unsigned int> &trk_nOuterLost();
	vector<unsigned int> &trk_nInnerLost();
	vector<unsigned int> &trk_algo();
	vector<unsigned int> &trk_originalAlgo();
	vector<ULong64_t> &trk_algoMask();
	vector<unsigned short> &trk_stopReason();
	vector<short> &trk_isHP();
	vector<int> &trk_seedIdx();
	vector<int> &trk_vtxIdx();
	vector<vector<float> > &trk_shareFrac();
	vector<vector<int> > &trk_simTrkIdx();
	vector<vector<int> > &trk_hitIdx();
	vector<vector<int> > &trk_hitType();
	vector<int> &sim_event();
	vector<int> &sim_bunchCrossing();
	vector<int> &sim_pdgId();
	vector<vector<int> > &sim_genPdgIds();
	vector<int> &sim_isFromBHadron();
	vector<float> &sim_px();
	vector<float> &sim_py();
	vector<float> &sim_pz();
	vector<float> &sim_pt();
	vector<float> &sim_eta();
	vector<float> &sim_phi();
	vector<float> &sim_pca_pt();
	vector<float> &sim_pca_eta();
	vector<float> &sim_pca_lambda();
	vector<float> &sim_pca_cotTheta();
	vector<float> &sim_pca_phi();
	vector<float> &sim_pca_dxy();
	vector<float> &sim_pca_dz();
	vector<int> &sim_q();
	vector<unsigned int> &sim_nValid();
	vector<unsigned int> &sim_nPixel();
	vector<unsigned int> &sim_nStrip();
	vector<unsigned int> &sim_nLay();
	vector<unsigned int> &sim_nPixelLay();
	vector<unsigned int> &sim_n3DLay();
	vector<vector<int> > &sim_trkIdx();
	vector<vector<float> > &sim_shareFrac();
	vector<vector<int> > &sim_seedIdx();
	vector<int> &sim_parentVtxIdx();
	vector<vector<int> > &sim_decayVtxIdx();
	vector<vector<int> > &sim_simHitIdx();
	vector<short> &pix_isBarrel();
	vector<unsigned short> &pix_det();
	vector<unsigned short> &pix_lay();
	vector<unsigned int> &pix_detId();
	vector<vector<int> > &pix_trkIdx();
	vector<vector<int> > &pix_seeIdx();
	vector<vector<int> > &pix_simHitIdx();
	vector<vector<float> > &pix_xySignificance();
	vector<vector<float> > &pix_chargeFraction();
	vector<unsigned short> &pix_simType();
	vector<float> &pix_x();
	vector<float> &pix_y();
	vector<float> &pix_z();
	vector<float> &pix_xx();
	vector<float> &pix_xy();
	vector<float> &pix_yy();
	vector<float> &pix_yz();
	vector<float> &pix_zz();
	vector<float> &pix_zx();
	vector<float> &pix_radL();
	vector<float> &pix_bbxi();
	vector<short> &ph2_isBarrel();
	vector<unsigned short> &ph2_det();
	vector<unsigned short> &ph2_lay();
	vector<unsigned int> &ph2_detId();
	vector<vector<int> > &ph2_trkIdx();
	vector<vector<int> > &ph2_seeIdx();
	vector<vector<int> > &ph2_simHitIdx();
	vector<vector<float> > &ph2_xySignificance();
	vector<unsigned short> &ph2_simType();
	vector<float> &ph2_x();
	vector<float> &ph2_y();
	vector<float> &ph2_z();
	vector<float> &ph2_xx();
	vector<float> &ph2_xy();
	vector<float> &ph2_yy();
	vector<float> &ph2_yz();
	vector<float> &ph2_zz();
	vector<float> &ph2_zx();
	vector<float> &ph2_radL();
	vector<float> &ph2_bbxi();
	vector<short> &inv_isBarrel();
	vector<unsigned short> &inv_det();
	vector<unsigned short> &inv_lay();
	vector<unsigned int> &inv_detId();
	vector<unsigned short> &inv_type();
	vector<unsigned short> &simhit_det();
	vector<unsigned short> &simhit_lay();
	vector<unsigned int> &simhit_detId();
	vector<float> &simhit_x();
	vector<float> &simhit_y();
	vector<float> &simhit_z();
	vector<float> &simhit_px();
	vector<float> &simhit_py();
	vector<float> &simhit_pz();
	vector<int> &simhit_particle();
	vector<short> &simhit_process();
	vector<float> &simhit_eloss();
	vector<float> &simhit_tof();
	vector<int> &simhit_simTrkIdx();
	vector<vector<int> > &simhit_hitIdx();
	vector<vector<int> > &simhit_hitType();
	float &bsp_x();
	float &bsp_y();
	float &bsp_z();
	float &bsp_sigmax();
	float &bsp_sigmay();
	float &bsp_sigmaz();
	vector<short> &see_fitok();
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
	vector<float> &see_statePt();
	vector<float> &see_stateTrajX();
	vector<float> &see_stateTrajY();
	vector<float> &see_stateTrajPx();
	vector<float> &see_stateTrajPy();
	vector<float> &see_stateTrajPz();
	vector<float> &see_stateTrajGlbX();
	vector<float> &see_stateTrajGlbY();
	vector<float> &see_stateTrajGlbZ();
	vector<float> &see_stateTrajGlbPx();
	vector<float> &see_stateTrajGlbPy();
	vector<float> &see_stateTrajGlbPz();
	vector<float> &see_stateCcov00();
	vector<float> &see_stateCcov01();
	vector<float> &see_stateCcov02();
	vector<float> &see_stateCcov03();
	vector<float> &see_stateCcov04();
	vector<float> &see_stateCcov05();
	vector<float> &see_stateCcov11();
	vector<float> &see_stateCcov12();
	vector<float> &see_stateCcov13();
	vector<float> &see_stateCcov14();
	vector<float> &see_stateCcov15();
	vector<float> &see_stateCcov22();
	vector<float> &see_stateCcov23();
	vector<float> &see_stateCcov24();
	vector<float> &see_stateCcov25();
	vector<float> &see_stateCcov33();
	vector<float> &see_stateCcov34();
	vector<float> &see_stateCcov35();
	vector<float> &see_stateCcov44();
	vector<float> &see_stateCcov45();
	vector<float> &see_stateCcov55();
	vector<int> &see_q();
	vector<unsigned int> &see_nValid();
	vector<unsigned int> &see_nPixel();
	vector<unsigned int> &see_nGlued();
	vector<unsigned int> &see_nStrip();
	vector<unsigned int> &see_nPhase2OT();
	vector<unsigned int> &see_algo();
	vector<unsigned short> &see_stopReason();
	vector<int> &see_trkIdx();
	vector<vector<float> > &see_shareFrac();
	vector<vector<int> > &see_simTrkIdx();
	vector<vector<int> > &see_hitIdx();
	vector<vector<int> > &see_hitType();
	vector<unsigned int> &see_offset();
	vector<float> &vtx_x();
	vector<float> &vtx_y();
	vector<float> &vtx_z();
	vector<float> &vtx_xErr();
	vector<float> &vtx_yErr();
	vector<float> &vtx_zErr();
	vector<float> &vtx_ndof();
	vector<float> &vtx_chi2();
	vector<short> &vtx_fake();
	vector<short> &vtx_valid();
	vector<vector<int> > &vtx_trkIdx();
	vector<int> &simvtx_event();
	vector<int> &simvtx_bunchCrossing();
	vector<unsigned int> &simvtx_processType();
	vector<float> &simvtx_x();
	vector<float> &simvtx_y();
	vector<float> &simvtx_z();
	vector<vector<int> > &simvtx_sourceSimIdx();
	vector<vector<int> > &simvtx_daughterSimIdx();
	vector<int> &simpv_idx();
}
#endif
