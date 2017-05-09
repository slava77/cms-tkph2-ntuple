#include "tkph2.h"
tkph2 cms2;
namespace tas {
        unsigned long long &event() { return cms2.event(); }
	unsigned int &lumi() { return cms2.lumi(); }
	unsigned int &run() { return cms2.run(); }
	vector<float> &trk_px() { return cms2.trk_px(); }
	vector<float> &trk_py() { return cms2.trk_py(); }
	vector<float> &trk_pz() { return cms2.trk_pz(); }
	vector<float> &trk_pt() { return cms2.trk_pt(); }
	vector<float> &trk_inner_px() { return cms2.trk_inner_px(); }
	vector<float> &trk_inner_py() { return cms2.trk_inner_py(); }
	vector<float> &trk_inner_pz() { return cms2.trk_inner_pz(); }
	vector<float> &trk_inner_pt() { return cms2.trk_inner_pt(); }
	vector<float> &trk_outer_px() { return cms2.trk_outer_px(); }
	vector<float> &trk_outer_py() { return cms2.trk_outer_py(); }
	vector<float> &trk_outer_pz() { return cms2.trk_outer_pz(); }
	vector<float> &trk_outer_pt() { return cms2.trk_outer_pt(); }
	vector<float> &trk_eta() { return cms2.trk_eta(); }
	vector<float> &trk_lambda() { return cms2.trk_lambda(); }
	vector<float> &trk_cotTheta() { return cms2.trk_cotTheta(); }
	vector<float> &trk_phi() { return cms2.trk_phi(); }
	vector<float> &trk_dxy() { return cms2.trk_dxy(); }
	vector<float> &trk_dz() { return cms2.trk_dz(); }
	vector<float> &trk_ptErr() { return cms2.trk_ptErr(); }
	vector<float> &trk_etaErr() { return cms2.trk_etaErr(); }
	vector<float> &trk_lambdaErr() { return cms2.trk_lambdaErr(); }
	vector<float> &trk_phiErr() { return cms2.trk_phiErr(); }
	vector<float> &trk_dxyErr() { return cms2.trk_dxyErr(); }
	vector<float> &trk_dzErr() { return cms2.trk_dzErr(); }
	vector<float> &trk_refpoint_x() { return cms2.trk_refpoint_x(); }
	vector<float> &trk_refpoint_y() { return cms2.trk_refpoint_y(); }
	vector<float> &trk_refpoint_z() { return cms2.trk_refpoint_z(); }
	vector<float> &trk_nChi2() { return cms2.trk_nChi2(); }
	vector<int> &trk_q() { return cms2.trk_q(); }
	vector<unsigned int> &trk_nValid() { return cms2.trk_nValid(); }
	vector<unsigned int> &trk_nInvalid() { return cms2.trk_nInvalid(); }
	vector<unsigned int> &trk_nPixel() { return cms2.trk_nPixel(); }
	vector<unsigned int> &trk_nStrip() { return cms2.trk_nStrip(); }
	vector<unsigned int> &trk_nPixelLay() { return cms2.trk_nPixelLay(); }
	vector<unsigned int> &trk_nStripLay() { return cms2.trk_nStripLay(); }
	vector<unsigned int> &trk_n3DLay() { return cms2.trk_n3DLay(); }
	vector<unsigned int> &trk_nOuterLost() { return cms2.trk_nOuterLost(); }
	vector<unsigned int> &trk_nInnerLost() { return cms2.trk_nInnerLost(); }
	vector<unsigned int> &trk_algo() { return cms2.trk_algo(); }
	vector<unsigned int> &trk_originalAlgo() { return cms2.trk_originalAlgo(); }
	vector<ULong64_t> &trk_algoMask() { return cms2.trk_algoMask(); }
	vector<unsigned short> &trk_stopReason() { return cms2.trk_stopReason(); }
	vector<short> &trk_isHP() { return cms2.trk_isHP(); }
	vector<int> &trk_seedIdx() { return cms2.trk_seedIdx(); }
	vector<int> &trk_vtxIdx() { return cms2.trk_vtxIdx(); }
	vector<vector<float> > &trk_shareFrac() { return cms2.trk_shareFrac(); }
	vector<vector<int> > &trk_simTrkIdx() { return cms2.trk_simTrkIdx(); }
	vector<vector<int> > &trk_hitIdx() { return cms2.trk_hitIdx(); }
	vector<vector<int> > &trk_hitType() { return cms2.trk_hitType(); }
	vector<int> &sim_event() { return cms2.sim_event(); }
	vector<int> &sim_bunchCrossing() { return cms2.sim_bunchCrossing(); }
	vector<int> &sim_pdgId() { return cms2.sim_pdgId(); }
	vector<vector<int> > &sim_genPdgIds() { return cms2.sim_genPdgIds(); }
	vector<int> &sim_isFromBHadron() { return cms2.sim_isFromBHadron(); }
	vector<float> &sim_px() { return cms2.sim_px(); }
	vector<float> &sim_py() { return cms2.sim_py(); }
	vector<float> &sim_pz() { return cms2.sim_pz(); }
	vector<float> &sim_pt() { return cms2.sim_pt(); }
	vector<float> &sim_eta() { return cms2.sim_eta(); }
	vector<float> &sim_phi() { return cms2.sim_phi(); }
	vector<float> &sim_pca_pt() { return cms2.sim_pca_pt(); }
	vector<float> &sim_pca_eta() { return cms2.sim_pca_eta(); }
	vector<float> &sim_pca_lambda() { return cms2.sim_pca_lambda(); }
	vector<float> &sim_pca_cotTheta() { return cms2.sim_pca_cotTheta(); }
	vector<float> &sim_pca_phi() { return cms2.sim_pca_phi(); }
	vector<float> &sim_pca_dxy() { return cms2.sim_pca_dxy(); }
	vector<float> &sim_pca_dz() { return cms2.sim_pca_dz(); }
	vector<int> &sim_q() { return cms2.sim_q(); }
	vector<unsigned int> &sim_nValid() { return cms2.sim_nValid(); }
	vector<unsigned int> &sim_nPixel() { return cms2.sim_nPixel(); }
	vector<unsigned int> &sim_nStrip() { return cms2.sim_nStrip(); }
	vector<unsigned int> &sim_nLay() { return cms2.sim_nLay(); }
	vector<unsigned int> &sim_nPixelLay() { return cms2.sim_nPixelLay(); }
	vector<unsigned int> &sim_n3DLay() { return cms2.sim_n3DLay(); }
	vector<vector<int> > &sim_trkIdx() { return cms2.sim_trkIdx(); }
	vector<vector<float> > &sim_shareFrac() { return cms2.sim_shareFrac(); }
	vector<vector<int> > &sim_seedIdx() { return cms2.sim_seedIdx(); }
	vector<int> &sim_parentVtxIdx() { return cms2.sim_parentVtxIdx(); }
	vector<vector<int> > &sim_decayVtxIdx() { return cms2.sim_decayVtxIdx(); }
	vector<vector<int> > &sim_simHitIdx() { return cms2.sim_simHitIdx(); }
	vector<short> &pix_isBarrel() { return cms2.pix_isBarrel(); }
	vector<unsigned short> &pix_det() { return cms2.pix_det(); }
	vector<unsigned short> &pix_lay() { return cms2.pix_lay(); }
	vector<unsigned int> &pix_detId() { return cms2.pix_detId(); }
	vector<vector<int> > &pix_trkIdx() { return cms2.pix_trkIdx(); }
	vector<vector<int> > &pix_seeIdx() { return cms2.pix_seeIdx(); }
	vector<vector<int> > &pix_simHitIdx() { return cms2.pix_simHitIdx(); }
	vector<vector<float> > &pix_chargeFraction() { return cms2.pix_chargeFraction(); }
	vector<unsigned short> &pix_simType() { return cms2.pix_simType(); }
	vector<float> &pix_x() { return cms2.pix_x(); }
	vector<float> &pix_y() { return cms2.pix_y(); }
	vector<float> &pix_z() { return cms2.pix_z(); }
	vector<float> &pix_xx() { return cms2.pix_xx(); }
	vector<float> &pix_xy() { return cms2.pix_xy(); }
	vector<float> &pix_yy() { return cms2.pix_yy(); }
	vector<float> &pix_yz() { return cms2.pix_yz(); }
	vector<float> &pix_zz() { return cms2.pix_zz(); }
	vector<float> &pix_zx() { return cms2.pix_zx(); }
	vector<float> &pix_radL() { return cms2.pix_radL(); }
	vector<float> &pix_bbxi() { return cms2.pix_bbxi(); }
	vector<short> &ph2_isBarrel() { return cms2.ph2_isBarrel(); }
	vector<unsigned short> &ph2_det() { return cms2.ph2_det(); }
	vector<unsigned short> &ph2_lay() { return cms2.ph2_lay(); }
	vector<unsigned int> &ph2_detId() { return cms2.ph2_detId(); }
	vector<vector<int> > &ph2_trkIdx() { return cms2.ph2_trkIdx(); }
	vector<vector<int> > &ph2_seeIdx() { return cms2.ph2_seeIdx(); }
	vector<vector<int> > &ph2_simHitIdx() { return cms2.ph2_simHitIdx(); }
	vector<unsigned short> &ph2_simType() { return cms2.ph2_simType(); }
	vector<float> &ph2_x() { return cms2.ph2_x(); }
	vector<float> &ph2_y() { return cms2.ph2_y(); }
	vector<float> &ph2_z() { return cms2.ph2_z(); }
	vector<float> &ph2_xx() { return cms2.ph2_xx(); }
	vector<float> &ph2_xy() { return cms2.ph2_xy(); }
	vector<float> &ph2_yy() { return cms2.ph2_yy(); }
	vector<float> &ph2_yz() { return cms2.ph2_yz(); }
	vector<float> &ph2_zz() { return cms2.ph2_zz(); }
	vector<float> &ph2_zx() { return cms2.ph2_zx(); }
	vector<float> &ph2_radL() { return cms2.ph2_radL(); }
	vector<float> &ph2_bbxi() { return cms2.ph2_bbxi(); }
	vector<short> &inv_isBarrel() { return cms2.inv_isBarrel(); }
	vector<unsigned short> &inv_det() { return cms2.inv_det(); }
	vector<unsigned short> &inv_lay() { return cms2.inv_lay(); }
	vector<unsigned int> &inv_detId() { return cms2.inv_detId(); }
	vector<unsigned short> &inv_type() { return cms2.inv_type(); }
	vector<unsigned short> &simhit_det() { return cms2.simhit_det(); }
	vector<unsigned short> &simhit_lay() { return cms2.simhit_lay(); }
	vector<unsigned int> &simhit_detId() { return cms2.simhit_detId(); }
	vector<float> &simhit_x() { return cms2.simhit_x(); }
	vector<float> &simhit_y() { return cms2.simhit_y(); }
	vector<float> &simhit_z() { return cms2.simhit_z(); }
	vector<float> &simhit_px() { return cms2.simhit_px(); }
	vector<float> &simhit_py() { return cms2.simhit_py(); }
	vector<float> &simhit_pz() { return cms2.simhit_pz(); }
	vector<int> &simhit_particle() { return cms2.simhit_particle(); }
	vector<short> &simhit_process() { return cms2.simhit_process(); }
	vector<float> &simhit_eloss() { return cms2.simhit_eloss(); }
	vector<float> &simhit_tof() { return cms2.simhit_tof(); }
	vector<int> &simhit_simTrkIdx() { return cms2.simhit_simTrkIdx(); }
	vector<vector<int> > &simhit_hitIdx() { return cms2.simhit_hitIdx(); }
	vector<vector<int> > &simhit_hitType() { return cms2.simhit_hitType(); }
	float &bsp_x() { return cms2.bsp_x(); }
	float &bsp_y() { return cms2.bsp_y(); }
	float &bsp_z() { return cms2.bsp_z(); }
	float &bsp_sigmax() { return cms2.bsp_sigmax(); }
	float &bsp_sigmay() { return cms2.bsp_sigmay(); }
	float &bsp_sigmaz() { return cms2.bsp_sigmaz(); }
	vector<short> &see_fitok() { return cms2.see_fitok(); }
	vector<float> &see_px() { return cms2.see_px(); }
	vector<float> &see_py() { return cms2.see_py(); }
	vector<float> &see_pz() { return cms2.see_pz(); }
	vector<float> &see_pt() { return cms2.see_pt(); }
	vector<float> &see_eta() { return cms2.see_eta(); }
	vector<float> &see_phi() { return cms2.see_phi(); }
	vector<float> &see_dxy() { return cms2.see_dxy(); }
	vector<float> &see_dz() { return cms2.see_dz(); }
	vector<float> &see_ptErr() { return cms2.see_ptErr(); }
	vector<float> &see_etaErr() { return cms2.see_etaErr(); }
	vector<float> &see_phiErr() { return cms2.see_phiErr(); }
	vector<float> &see_dxyErr() { return cms2.see_dxyErr(); }
	vector<float> &see_dzErr() { return cms2.see_dzErr(); }
	vector<float> &see_chi2() { return cms2.see_chi2(); }
	vector<float> &see_statePt() { return cms2.see_statePt(); }
	vector<float> &see_stateTrajX() { return cms2.see_stateTrajX(); }
	vector<float> &see_stateTrajY() { return cms2.see_stateTrajY(); }
	vector<float> &see_stateTrajPx() { return cms2.see_stateTrajPx(); }
	vector<float> &see_stateTrajPy() { return cms2.see_stateTrajPy(); }
	vector<float> &see_stateTrajPz() { return cms2.see_stateTrajPz(); }
	vector<int> &see_q() { return cms2.see_q(); }
	vector<unsigned int> &see_nValid() { return cms2.see_nValid(); }
	vector<unsigned int> &see_nPixel() { return cms2.see_nPixel(); }
	vector<unsigned int> &see_nGlued() { return cms2.see_nGlued(); }
	vector<unsigned int> &see_nStrip() { return cms2.see_nStrip(); }
	vector<unsigned int> &see_nPhase2OT() { return cms2.see_nPhase2OT(); }
	vector<unsigned int> &see_algo() { return cms2.see_algo(); }
	vector<unsigned short> &see_stopReason() { return cms2.see_stopReason(); }
	vector<int> &see_trkIdx() { return cms2.see_trkIdx(); }
	vector<vector<float> > &see_shareFrac() { return cms2.see_shareFrac(); }
	vector<vector<int> > &see_simTrkIdx() { return cms2.see_simTrkIdx(); }
	vector<vector<int> > &see_hitIdx() { return cms2.see_hitIdx(); }
	vector<vector<int> > &see_hitType() { return cms2.see_hitType(); }
	vector<unsigned int> &see_offset() { return cms2.see_offset(); }
	vector<float> &vtx_x() { return cms2.vtx_x(); }
	vector<float> &vtx_y() { return cms2.vtx_y(); }
	vector<float> &vtx_z() { return cms2.vtx_z(); }
	vector<float> &vtx_xErr() { return cms2.vtx_xErr(); }
	vector<float> &vtx_yErr() { return cms2.vtx_yErr(); }
	vector<float> &vtx_zErr() { return cms2.vtx_zErr(); }
	vector<float> &vtx_ndof() { return cms2.vtx_ndof(); }
	vector<float> &vtx_chi2() { return cms2.vtx_chi2(); }
	vector<short> &vtx_fake() { return cms2.vtx_fake(); }
	vector<short> &vtx_valid() { return cms2.vtx_valid(); }
	vector<vector<int> > &vtx_trkIdx() { return cms2.vtx_trkIdx(); }
	vector<int> &simvtx_event() { return cms2.simvtx_event(); }
	vector<int> &simvtx_bunchCrossing() { return cms2.simvtx_bunchCrossing(); }
	vector<unsigned int> &simvtx_processType() { return cms2.simvtx_processType(); }
	vector<float> &simvtx_x() { return cms2.simvtx_x(); }
	vector<float> &simvtx_y() { return cms2.simvtx_y(); }
	vector<float> &simvtx_z() { return cms2.simvtx_z(); }
	vector<vector<int> > &simvtx_sourceSimIdx() { return cms2.simvtx_sourceSimIdx(); }
	vector<vector<int> > &simvtx_daughterSimIdx() { return cms2.simvtx_daughterSimIdx(); }
	vector<int> &simpv_idx() { return cms2.simpv_idx(); }
}
