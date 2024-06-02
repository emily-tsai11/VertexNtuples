#include "../interface/SecondaryVertex.h"


SecondaryVertex::SecondaryVertex(const reco::Vertex& sv,
    const reco::Vertex& primaryVertex) {

  // Track times aren't actually used, just initialized
  // and filled with dummy values to avoid errors
  trk_tval_ = new std::vector<double>;
  trk_terr_ = new std::vector<double>;
  trk_tsig_ = new std::vector<double>;
  trk_tqual_ = new std::vector<double>;

  trk_pt_ = new std::vector<double>;
  trk_pt2_ = new std::vector<double>;
  trk_eta_ = new std::vector<double>;
  trk_phi_ = new std::vector<double>;
  trk_dxy_ = new std::vector<double>;
  trk_dz_ = new std::vector<double>;
  trk_d3d_ = new std::vector<double>;
  trk_dxyerr_ = new std::vector<double>;
  trk_dzerr_ = new std::vector<double>;
  trk_d3derr_ = new std::vector<double>;
  trk_dxysig_ = new std::vector<double>;
  trk_dzsig_ = new std::vector<double>;
  trk_d3dsig_ = new std::vector<double>;
  trk_charge_ = new std::vector<double>;
  trk_chi2_ = new std::vector<double>;
  trk_ndof_ = new std::vector<double>;
  trk_chi2dof_ = new std::vector<double>;

  for (const reco::TrackBaseRef& trkRef : sv.tracks()) {
    double d3d = vertexntuples::d3d(trkRef);
    double d3derr = vertexntuples::d3dErr(trkRef);
    double d3dsig = d3d / d3derr;

    trk_tval_->push_back(-1.0); // Dummy value
    trk_terr_->push_back(-1.0); // Dummy value
    trk_tsig_->push_back(-1.0); // Dummy value
    trk_tqual_->push_back(-1.0); // Dummy value
    trk_pt_->push_back(trkRef->pt());
    trk_pt2_->push_back(trkRef->pt2());
    trk_eta_->push_back(trkRef->eta());
    trk_phi_->push_back(trkRef->phi());
    trk_dxy_->push_back(trkRef->dxy());
    trk_dz_->push_back(trkRef->dz());
    trk_d3d_->push_back(d3d);
    trk_dxyerr_->push_back(trkRef->dxyError());
    trk_dzerr_->push_back(trkRef->dzError());
    trk_d3derr_->push_back(d3derr);
    trk_dxysig_->push_back(trkRef->dxy() / trkRef->dxyError());
    trk_dzsig_->push_back(trkRef->dz() / trkRef->dzError());
    trk_d3dsig_->push_back(d3dsig);
    trk_charge_->push_back(trkRef->charge());
    trk_chi2_->push_back(trkRef->chi2());
    trk_ndof_->push_back(trkRef->ndof());
    trk_chi2dof_->push_back(trkRef->normalizedChi2());
  }
  Measurement1D dxy = vertexntuples::dxy(sv, primaryVertex);
  Measurement1D d3d = vertexntuples::d3d(sv, primaryVertex);

  x_ = sv.x();
  y_ = sv.y();
  z_ = sv.z();
  xerr_ = sv.xError();
  yerr_ = sv.yError();
  zerr_ = sv.zError();
  dxy_ = dxy.value();
  dz_ = vertexntuples::dz(sv, primaryVertex);
  d3d_ = d3d.value();
  dxyerr_ = dxy.error();
  dzerr_ = vertexntuples::dzErr(sv, primaryVertex);
  d3derr_ = d3d.error();
  dxysig_ = dxy.significance();
  dzsig_ = dz_ / dzerr_;
  d3dsig_ = d3d.significance();
  pt_ = sv.p4().Pt();
  pt2_ = sv.p4().Pt() * sv.p4().Pt();
  eta_ = sv.p4().Eta();
  phi_ = sv.p4().Phi();
  tavg_ = -1.0; // Dummy value
  trange_ = -1.0; // Dummy value
  chi2_ = sv.chi2();
  ndof_ = sv.ndof();
  chi2dof_ = sv.normalizedChi2();
  ntrk_ = sv.tracksSize();

  // For matching to a GenVertex
  trk_deltaR_ = new std::vector<double>;
  trk_ptResNorm_ = new std::vector<double>;
}


SecondaryVertex::SecondaryVertex(const reco::Vertex& sv,
    const reco::Vertex& primaryVertex,
    const edm::ValueMap<float>& trackTimeValueMap,
    const edm::ValueMap<float>& trackTimeErrorMap,
    const edm::ValueMap<float>& trackTimeQualityMap) {

  trk_tval_ = new std::vector<double>;
  trk_terr_ = new std::vector<double>;
  trk_tsig_ = new std::vector<double>;
  trk_tqual_ = new std::vector<double>;
  trk_pt_ = new std::vector<double>;
  trk_pt2_ = new std::vector<double>;
  trk_eta_ = new std::vector<double>;
  trk_phi_ = new std::vector<double>;
  trk_dxy_ = new std::vector<double>;
  trk_dz_ = new std::vector<double>;
  trk_d3d_ = new std::vector<double>;
  trk_dxyerr_ = new std::vector<double>;
  trk_dzerr_ = new std::vector<double>;
  trk_d3derr_ = new std::vector<double>;
  trk_dxysig_ = new std::vector<double>;
  trk_dzsig_ = new std::vector<double>;
  trk_d3dsig_ = new std::vector<double>;
  trk_charge_ = new std::vector<double>;
  trk_chi2_ = new std::vector<double>;
  trk_ndof_ = new std::vector<double>;
  trk_chi2dof_ = new std::vector<double>;

  float tavg = 0.0;
  float tmin = 99999.9;
  float tmax = -99999.9;
  for (const reco::TrackBaseRef& trkRef : sv.tracks()) {
    tavg += trackTimeValueMap[trkRef];
    tmin = std::min(tmin, trackTimeValueMap[trkRef]);
    tmax = std::max(tmax, trackTimeValueMap[trkRef]);

    double d3d = vertexntuples::d3d(trkRef);
    double d3derr = vertexntuples::d3dErr(trkRef);
    double d3dsig = d3d / d3derr;

    trk_tval_->push_back(trackTimeValueMap[trkRef]);
    trk_terr_->push_back(trackTimeErrorMap[trkRef]);
    trk_tsig_->push_back(trackTimeValueMap[trkRef] / trackTimeErrorMap[trkRef]);
    trk_tqual_->push_back(trackTimeQualityMap[trkRef]);
    trk_pt_->push_back(trkRef->pt());
    trk_pt2_->push_back(trkRef->pt2());
    trk_eta_->push_back(trkRef->eta());
    trk_phi_->push_back(trkRef->phi());
    trk_dxy_->push_back(trkRef->dxy());
    trk_dz_->push_back(trkRef->dz());
    trk_d3d_->push_back(d3d);
    trk_dxyerr_->push_back(trkRef->dxyError());
    trk_dzerr_->push_back(trkRef->dzError());
    trk_d3derr_->push_back(d3derr);
    trk_dxysig_->push_back(trkRef->dxy() / trkRef->dxyError());
    trk_dzsig_->push_back(trkRef->dz() / trkRef->dzError());
    trk_d3dsig_->push_back(d3dsig);
    trk_charge_->push_back(trkRef->charge());
    trk_chi2_->push_back(trkRef->chi2());
    trk_ndof_->push_back(trkRef->ndof());
    trk_chi2dof_->push_back(trkRef->normalizedChi2());
  }
  tavg /= (float) sv.tracksSize();
  float trange = tmax - tmin;
  Measurement1D dxy = vertexntuples::dxy(sv, primaryVertex);
  Measurement1D d3d = vertexntuples::d3d(sv, primaryVertex);

  x_ = sv.x();
  y_ = sv.y();
  z_ = sv.z();
  xerr_ = sv.xError();
  yerr_ = sv.yError();
  zerr_ = sv.zError();
  dxy_ = dxy.value();
  dz_ = vertexntuples::dz(sv, primaryVertex);
  d3d_ = d3d.value();
  dxyerr_ = dxy.error();
  dzerr_ = vertexntuples::dzErr(sv, primaryVertex);
  d3derr_ = d3d.error();
  dxysig_ = dxy.significance();
  dzsig_ = dz_ / dzerr_;
  d3dsig_ = d3d.significance();
  pt_ = sv.p4().Pt();
  pt2_ = sv.p4().Pt() * sv.p4().Pt();
  eta_ = sv.p4().Eta();
  phi_ = sv.p4().Phi();
  tavg_ = tavg;
  trange_ = trange;
  chi2_ = sv.chi2();
  ndof_ = sv.ndof();
  chi2dof_ = sv.normalizedChi2();
  ntrk_ = sv.tracksSize();

  // For matching to a GenVertex
  trk_deltaR_ = new std::vector<double>;
  trk_ptResNorm_ = new std::vector<double>;
}


// SecondaryVertex::~SecondaryVertex() {}


void SecondaryVertex::fill(std::map<TString, TH1F*>& histos1,
    std::map<TString, TH2F*>& histos2, TString prefix) {

  for (unsigned int iTrk = 0; iTrk < nTracks(); iTrk++) {
    histos1[prefix + "_trk_tval"]->Fill(trk_tval_->at(iTrk));
    histos1[prefix + "_trk_terr"]->Fill(trk_terr_->at(iTrk));
    histos1[prefix + "_trk_tsig"]->Fill(trk_tsig_->at(iTrk));
    histos1[prefix + "_trk_tqual"]->Fill(trk_tqual_->at(iTrk));
    histos1[prefix + "_trk_pt"]->Fill(trk_pt_->at(iTrk));
    histos1[prefix + "_trk_pt2"]->Fill(trk_pt2_->at(iTrk));
    histos1[prefix + "_trk_eta"]->Fill(trk_eta_->at(iTrk));
    histos1[prefix + "_trk_phi"]->Fill(trk_phi_->at(iTrk));
    histos1[prefix + "_trk_dxy"]->Fill(trk_dxy_->at(iTrk));
    histos1[prefix + "_trk_dz"]->Fill(trk_dz_->at(iTrk));
    histos1[prefix + "_trk_d3d"]->Fill(trk_d3d_->at(iTrk));
    histos1[prefix + "_trk_dxyerr"]->Fill(trk_dxyerr_->at(iTrk));
    histos1[prefix + "_trk_dzerr"]->Fill(trk_dzerr_->at(iTrk));
    histos1[prefix + "_trk_d3derr"]->Fill(trk_d3derr_->at(iTrk));
    histos1[prefix + "_trk_dxysig"]->Fill(trk_dxysig_->at(iTrk));
    histos1[prefix + "_trk_dzsig"]->Fill(trk_dzsig_->at(iTrk));
    histos1[prefix + "_trk_d3dsig"]->Fill(trk_d3dsig_->at(iTrk));
    histos1[prefix + "_trk_charge"]->Fill(trk_charge_->at(iTrk));
    histos1[prefix + "_trk_chi2"]->Fill(trk_chi2_->at(iTrk));
    histos1[prefix + "_trk_ndof"]->Fill(trk_ndof_->at(iTrk));
    histos1[prefix + "_trk_chi2dof"]->Fill(trk_chi2dof_->at(iTrk));

    histos2[prefix + "_trk_eta_tval"]->Fill(trk_eta_->at(iTrk), trk_tval_->at(iTrk));
    histos2[prefix + "_trk_eta_terr"]->Fill(trk_eta_->at(iTrk), trk_terr_->at(iTrk));
    histos2[prefix + "_trk_eta_tsig"]->Fill(trk_eta_->at(iTrk), trk_tsig_->at(iTrk));
    histos2[prefix + "_trk_eta_tqual"]->Fill(trk_eta_->at(iTrk), trk_tqual_->at(iTrk));
  }

  histos1[prefix + "_x"]->Fill(x_);
  histos1[prefix + "_y"]->Fill(y_);
  histos1[prefix + "_z"]->Fill(z_);
  histos1[prefix + "_xerr"]->Fill(xerr_);
  histos1[prefix + "_yerr"]->Fill(yerr_);
  histos1[prefix + "_zerr"]->Fill(zerr_);
  histos1[prefix + "_dxy"]->Fill(dxy_);
  histos1[prefix + "_dz"]->Fill(dz_);
  histos1[prefix + "_d3d"]->Fill(d3d_);
  histos1[prefix + "_dxyerr"]->Fill(dxyerr_);
  histos1[prefix + "_dzerr"]->Fill(dzerr_);
  histos1[prefix + "_d3derr"]->Fill(d3derr_);
  histos1[prefix + "_dxysig"]->Fill(dxysig_);
  histos1[prefix + "_dzsig"]->Fill(dzsig_);
  histos1[prefix + "_d3dsig"]->Fill(d3dsig_);
  histos1[prefix + "_pt"]->Fill(pt_);
  histos1[prefix + "_pt2"]->Fill(pt2_);
  histos1[prefix + "_eta"]->Fill(eta_);
  histos1[prefix + "_phi"]->Fill(phi_);
  histos1[prefix + "_tavg"]->Fill(tavg_);
  histos1[prefix + "_trange"]->Fill(trange_);
  histos1[prefix + "_chi2"]->Fill(chi2_);
  histos1[prefix + "_ndof"]->Fill(ndof_);
  histos1[prefix + "_chi2dof"]->Fill(chi2dof_);
  histos1[prefix + "_ntrk"]->Fill(ntrk_);

  histos2[prefix + "_eta_tavg"]->Fill(eta_, tavg_);
  histos2[prefix + "_eta_trange"]->Fill(eta_, trange_);
  histos2[prefix + "_trange_pt"]->Fill(trange_, pt_);
  histos2[prefix + "_trange_pt2"]->Fill(trange_, pt2_);
  histos2[prefix + "_trange_dxy"]->Fill(trange_, dxy_);
  histos2[prefix + "_trange_dxysig"]->Fill(trange_, dxysig_);
  histos2[prefix + "_trange_d3d"]->Fill(trange_, d3d_);
  histos2[prefix + "_trange_d3dsig"]->Fill(trange_, d3dsig_);
}
