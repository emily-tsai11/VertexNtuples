#include "../interface/SecondaryVertex.h"


SecondaryVertex::SecondaryVertex(const reco::Vertex& sv,
    const std::vector<reco::TrackBaseRef>* tracks,
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

  trk_deltaR_ = new std::vector<double>;
  trk_ptresnorm_ = new std::vector<double>;

  float tavg = 0.0;
  float tmin = 99999.9;
  float tmax = -99999.9;
  for (const reco::TrackBaseRef& trkRef : *tracks) {
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
  tavg /= (float) tracks->size();
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
  ntrk_ = tracks->size();
}


// SecondaryVertex::~SecondaryVertex() {}


void SecondaryVertex::fill(std::map<TString, TH1F*>& histos1,
    std::map<TString, TH2F*>& histos2, TString prefix) {

  for (unsigned int iTrk = 0; iTrk < nTracks(); iTrk++) {
    histos1[prefix + "_trk_tval"]->Fill(trkTVal()->at(iTrk));
    histos1[prefix + "_trk_terr"]->Fill(trkTErr()->at(iTrk));
    histos1[prefix + "_trk_tsig"]->Fill(trkTSig()->at(iTrk));
    histos1[prefix + "_trk_tqual"]->Fill(trkTQual()->at(iTrk));
    histos1[prefix + "_trk_pt"]->Fill(trkPt()->at(iTrk));
    histos1[prefix + "_trk_pt2"]->Fill(trkPt2()->at(iTrk));
    histos1[prefix + "_trk_eta"]->Fill(trkEta()->at(iTrk));
    histos1[prefix + "_trk_phi"]->Fill(trkPhi()->at(iTrk));
    histos1[prefix + "_trk_dxy"]->Fill(trkDxy()->at(iTrk));
    histos1[prefix + "_trk_dz"]->Fill(trkDz()->at(iTrk));
    histos1[prefix + "_trk_d3d"]->Fill(trkD3d()->at(iTrk));
    histos1[prefix + "_trk_dxyerr"]->Fill(trkDxyErr()->at(iTrk));
    histos1[prefix + "_trk_dzerr"]->Fill(trkDzErr()->at(iTrk));
    histos1[prefix + "_trk_d3derr"]->Fill(trkD3dErr()->at(iTrk));
    histos1[prefix + "_trk_dxysig"]->Fill(trkDxySig()->at(iTrk));
    histos1[prefix + "_trk_dzsig"]->Fill(trkDzSig()->at(iTrk));
    histos1[prefix + "_trk_d3dsig"]->Fill(trkD3dSig()->at(iTrk));
    histos1[prefix + "_trk_charge"]->Fill(trkCharge()->at(iTrk));
    histos1[prefix + "_trk_chi2"]->Fill(trkChi2()->at(iTrk));
    histos1[prefix + "_trk_ndof"]->Fill(trkNDOF()->at(iTrk));
    histos1[prefix + "_trk_chi2dof"]->Fill(trkChi2DOF()->at(iTrk));
  }

  histos1[prefix + "_x"]->Fill(x());
  histos1[prefix + "_y"]->Fill(y());
  histos1[prefix + "_z"]->Fill(z());
  histos1[prefix + "_xerr"]->Fill(xErr());
  histos1[prefix + "_yerr"]->Fill(yErr());
  histos1[prefix + "_zerr"]->Fill(zErr());
  histos1[prefix + "_dxy"]->Fill(dxy());
  histos1[prefix + "_dz"]->Fill(dz());
  histos1[prefix + "_d3d"]->Fill(d3d());
  histos1[prefix + "_dxyerr"]->Fill(dxyErr());
  histos1[prefix + "_dzerr"]->Fill(dzErr());
  histos1[prefix + "_d3derr"]->Fill(d3dErr());
  histos1[prefix + "_dxysig"]->Fill(dxySig());
  histos1[prefix + "_dzsig"]->Fill(dzSig());
  histos1[prefix + "_d3dsig"]->Fill(d3dSig());
  histos1[prefix + "_pt"]->Fill(pt());
  histos1[prefix + "_pt2"]->Fill(pt2());
  histos1[prefix + "_eta"]->Fill(eta());
  histos1[prefix + "_phi"]->Fill(phi());
  histos1[prefix + "_tavg"]->Fill(tAvg());
  histos1[prefix + "_trange"]->Fill(tRange());
  histos1[prefix + "_chi2"]->Fill(chi2());
  histos1[prefix + "_ndof"]->Fill(nDOF());
  histos1[prefix + "_chi2dof"]->Fill(chi2DOF());
  histos1[prefix + "_ntrk"]->Fill(nTracks());

  histos2[prefix + "_trange_pt"]->Fill(tRange(), pt());
  histos2[prefix + "_trange_pt2"]->Fill(tRange(), pt2());
  histos2[prefix + "_trange_dxy"]->Fill(tRange(), dxy());
  histos2[prefix + "_trange_dxysig"]->Fill(tRange(), dxySig());
  histos2[prefix + "_trange_d3d"]->Fill(tRange(), d3d());
  histos2[prefix + "_trange_d3dsig"]->Fill(tRange(), d3dSig());
}
