#include "../interface/SecondaryVertex.h"


SecondaryVertex::SecondaryVertex(const reco::VertexCompositePtrCandidate& sv,
    const reco::Vertex& primaryVertex) {

  initialize(sv, primaryVertex, false);
}


SecondaryVertex::SecondaryVertex(const reco::VertexCompositePtrCandidate& sv,
    const reco::Vertex& primaryVertex,
    edm::Handle<reco::TrackCollection> generalTracksHandle,
    const edm::ValueMap<float>& trackT0,
    const edm::ValueMap<float>& trackSigmaT0
    // , const edm::ValueMap<float>& trackQuality
    ) {

  initialize(sv, primaryVertex);
  initializeTime(sv, generalTracksHandle, trackT0, trackSigmaT0);
  // initializeTime(sv, trackT0, trackSigmaT0, trackQuality);
}


void SecondaryVertex::initialize(const reco::VertexCompositePtrCandidate& sv,
    const reco::Vertex& primaryVertex, bool hasTime) {

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
  trk_charge_ = new std::vector<float>;
  trk_chi2_ = new std::vector<float>;
  trk_ndof_ = new std::vector<float>;
  trk_chi2dof_ = new std::vector<float>;

  trk_tval_ = new std::vector<float>;
  trk_terr_ = new std::vector<float>;
  trk_tsig_ = new std::vector<float>;
  // trk_tqual_ = new std::vector<float>;

  for (unsigned int iTrk = 0; iTrk < sv.daughterPtrVector().size(); iTrk++) {
    const reco::CandidatePtr& trkPtr = sv.daughterPtr(iTrk);

    float dxy = vertexntuples::dxy(trkPtr, primaryVertex);
    float dz = vertexntuples::dz(trkPtr, primaryVertex);
    float d3d = vertexntuples::d3d(trkPtr, primaryVertex);
    float d3derr = vertexntuples::d3dErr(trkPtr, primaryVertex);
    trk_pt_->push_back(trkPtr->pt());
    trk_pt2_->push_back(trkPtr->pt() * trkPtr->pt());
    trk_eta_->push_back(trkPtr->eta());
    trk_phi_->push_back(trkPtr->phi());
    trk_dxy_->push_back(dxy);
    trk_dz_->push_back(dz);
    trk_d3d_->push_back(d3d);
    trk_dxyerr_->push_back(trkPtr->dxyError());
    trk_dzerr_->push_back(trkPtr->dzError());
    trk_d3derr_->push_back(d3derr);
    trk_dxysig_->push_back(dxy / trkPtr->dxyError());
    trk_dzsig_->push_back(dz / trkPtr->dzError());
    trk_d3dsig_->push_back(d3d / d3derr);
    trk_charge_->push_back(trkPtr->charge());
    // Track quality parameters are protected in PackedCandidate.h...
    trk_chi2_->push_back(-1.0); // Dummy value
    trk_ndof_->push_back(-1.0); // Dummy value
    trk_chi2dof_->push_back(-1.0); // Dummy value

    if (!hasTime) {
      trk_tval_->push_back(-1.0); // Dummy value
      trk_terr_->push_back(-1.0); // Dummy value
      trk_tsig_->push_back(-1.0); // Dummy value
      // trk_tqual_->push_back(-1.0); // Dummy value
    }
  }

  Measurement1D dxy = vertexntuples::dxy(sv, primaryVertex);
  Measurement1D d3d = vertexntuples::d3d(sv, primaryVertex);
  x_ = sv.vx();
  y_ = sv.vy();
  z_ = sv.vz();
  // Don't know how to get xyz errors from a VertexCompositePtrCandidate object
  xerr_ = -1.0; // Dummy value
  yerr_ = -1.0; // Dummy value
  zerr_ = -1.0; // Dummy value
  dxy_ = dxy.value();
  dz_ = vertexntuples::dz(sv, primaryVertex);
  d3d_ = d3d.value();
  dxyerr_ = dxy.error();
  dzerr_ = sv.dzError();
  d3derr_ = d3d.error();
  dxysig_ = dxy.significance();
  dzsig_ = dz_ / dzerr_;
  d3dsig_ = d3d.significance();
  pt_ = sv.pt();
  pt2_ = sv.pt() * sv.pt();
  eta_ = sv.eta();
  phi_ = sv.phi();
  chi2_ = sv.vertexChi2();
  ndof_ = sv.vertexNdof();
  chi2dof_ = sv.vertexNormalizedChi2();
  ntrk_ = sv.daughterPtrVector().size();

  if (!hasTime) {
    tavg_ = -1.0;
    trange_ = -1.0;
  }

  // For matching to a GenVertex
  trk_deltaR_ = new std::vector<float>;
  trk_normPtRes_ = new std::vector<float>;
}


void SecondaryVertex::initializeTime(const reco::VertexCompositePtrCandidate& sv,
    edm::Handle<reco::TrackCollection> generalTracksHandle,
    const edm::ValueMap<float>& trackT0,
    const edm::ValueMap<float>& trackSigmaT0
    // , const edm::ValueMap<float>& trackQuality
    ) {

  float tsum = 0.0;
  float tmin = 9999.9;
  float tmax = -9999.9;
  for (unsigned int iTrk = 0; iTrk < sv.daughterPtrVector().size(); iTrk++) {
    const reco::CandidatePtr& trkPtr = sv.daughterPtr(iTrk);

    const reco::Track* bestTrk = trkPtr->bestTrack();
    int matchIdx = -1;
    if (bestTrk) {
      const reco::TrackCollection generalTracks = *generalTracksHandle;
      for (unsigned int iTrk = 0; iTrk < generalTracks.size(); iTrk++) {
        const reco::Track& trk = generalTracks.at(iTrk);
        if (trk.pt() != bestTrk->pt() || trk.eta() != bestTrk->eta() || trk.eta() != bestTrk->eta()) continue;
        matchIdx = iTrk;
      }
    }

    // Found best track
    bool filled = false;
    if (matchIdx >= 0) {
      const reco::TrackRef trkRef(generalTracksHandle, (unsigned int) matchIdx);
      if (trackT0[trkRef] != 0.0 && trackSigmaT0[trkRef] != -1.0) {
        filled = true;
        tsum += trackT0[trkRef];
        tmin = std::min(tmin, trackT0[trkRef]);
        tmax = std::max(tmax, trackT0[trkRef]);
        trk_tval_->push_back(trackT0[trkRef]*1000.0);
        trk_terr_->push_back(trackSigmaT0[trkRef]*1000.0);
        trk_tsig_->push_back(trackT0[trkRef] / trackSigmaT0[trkRef]);
        // trk_tqual_->push_back(trackQuality[trkRef]);
      }
    }

    // No best track found or no time value
    if (!filled) {
      trk_tval_->push_back(-1000.0);
      trk_terr_->push_back(-1000.0);
      trk_tsig_->push_back(-100.0);
      // trk_tqual_->push_back(-1.0);
    }
  } // End loop over daughters/tracks

  tavg_ = (tsum / (float) sv.daughterPtrVector().size())*1000.0;
  trange_ = (tmax - tmin)*1000.0;
}


// SecondaryVertex::~SecondaryVertex() {}


void SecondaryVertex::fill(std::map<TString, TH1F*>& histos1,
    std::map<TString, TH2F*>& histos2, TString prefix) {

  for (unsigned int iTrk = 0; iTrk < nTracks(); iTrk++) {
    histos1[prefix + "_trk_tval"]->Fill(trk_tval_->at(iTrk));
    histos1[prefix + "_trk_terr"]->Fill(trk_terr_->at(iTrk));
    histos1[prefix + "_trk_tsig"]->Fill(trk_tsig_->at(iTrk));
    // histos1[prefix + "_trk_tqual"]->Fill(trk_tqual_->at(iTrk));
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
    // histos2[prefix + "_trk_eta_tqual"]->Fill(trk_eta_->at(iTrk), trk_tqual_->at(iTrk));
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
