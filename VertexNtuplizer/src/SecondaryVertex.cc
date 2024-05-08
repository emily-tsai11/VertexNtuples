#include "../interface/SecondaryVertex.h"


SecondaryVertex::SecondaryVertex(const reco::Vertex& sv,
    const std::vector<reco::TrackBaseRef>* tracks,
    const reco::Vertex& primaryVertex,
    const edm::ValueMap<float>& trackTimeValueMap,
    const edm::ValueMap<float>& trackTimeErrorMap,
    const edm::ValueMap<float>& trackTimeQualityMap) {

  trk_tval_ = new std::vector<float>;
  trk_terr_ = new std::vector<float>;
  trk_tsig_ = new std::vector<float>;
  trk_tqual_ = new std::vector<float>;
  trk_pt_ = new std::vector<float>;
  trk_pt2_ = new std::vector<float>;
  trk_eta_ = new std::vector<float>;
  trk_phi_ = new std::vector<float>;
  trk_dxy_ = new std::vector<float>;
  trk_dz_ = new std::vector<float>;
  trk_d3d_ = new std::vector<float>;
  trk_dxyerr_ = new std::vector<float>;
  trk_dzerr_ = new std::vector<float>;
  trk_d3derr_ = new std::vector<float>;
  trk_dxysig_ = new std::vector<float>;
  trk_dzsig_ = new std::vector<float>;
  trk_d3dsig_ = new std::vector<float>;
  trk_chi2_ = new std::vector<float>;
  trk_ndof_ = new std::vector<float>;
  trk_chi2dof_ = new std::vector<float>;

  float tavg = 0.0;
  float tmin = 99999.9;
  float tmax = -99999.9;
  for (const reco::TrackBaseRef& trkRef : sv.tracks()) {
    tavg += trackTimeValueMap[trkRef];
    tmin = std::min(tmin, trackTimeValueMap[trkRef]);
    tmax = std::max(tmax, trackTimeValueMap[trkRef]);

    float dxy2 = trkRef->dxy()*trkRef->dxy();
    float dz2 = trkRef->dz()*trkRef->dz();
    float dxy2err = 2*trkRef->dxyError()*trkRef->dxy();
    float dz2err = 2*trkRef->dzError()*trkRef->dz();
    float d3d2 = dxy2 + dz2;
    float d3d2err = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
    float d3d = TMath::Sqrt(d3d2);
    float d3derr = 0.5*d3d2err/d3d;

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
    trk_d3dsig_->push_back(d3d / d3derr);
    trk_chi2_->push_back(trkRef->chi2());
    trk_ndof_->push_back(trkRef->ndof());
    trk_chi2dof_->push_back(trkRef->normalizedChi2());
  }
  tavg /= (float) tracks->size();
  float trange = tmax - tmin;
  Measurement1D dxy = getDxy(sv, primaryVertex);
  Measurement1D d3d = getD3d(sv, primaryVertex);

  x_ = sv.x();
  y_ = sv.y();
  z_ = sv.z();
  xerr_ = sv.xError();
  yerr_ = sv.yError();
  zerr_ = sv.zError();
  dxy_ = dxy.value();
  dz_ = abs(sv.z() - primaryVertex.z());
  d3d_ = d3d.value();
  dxyerr_ = dxy.error();
  dzerr_ = TMath::Sqrt(sv.zError()*sv.zError() + primaryVertex.zError()*primaryVertex.zError());
  d3derr_ = dxy.error();
  dxysig_ = dxy.significance();
  dzsig_ = dz_ == 0.0 ? 0.0 : dz_ / dzerr_;
  d3dsig_ = d3d.significance();
  pt_ = sv.p4().Pt();
  pt2_ = sv.p4().Pt()*sv.p4().Pt();
  eta_ = sv.p4().Eta();
  phi_ = sv.p4().Phi();
  tavg_ = tavg;
  trange_ = trange;
  chi2_ = sv.chi2();
  ndof_ = sv.ndof();
  chi2dof_ = sv.normalizedChi2();
  ntrk_ = tracks->size();
}


// SecondaryVertex::~SecondaryVertex() {}


Measurement1D SecondaryVertex::getDxy(const reco::Vertex sv, const reco::Vertex pv) {

  VertexDistanceXY dist;
  return dist.distance(sv, pv);
}


Measurement1D SecondaryVertex::getD3d(const reco::Vertex sv, const reco::Vertex pv) {

  VertexDistance3D dist;
  return dist.distance(sv, pv);
}
