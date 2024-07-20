#include "../interface/GenVertex.h"


GenVertex::GenVertex(const reco::Candidate* mother, std::vector<const reco::Candidate*>* daughters,
    const reco::Vertex& primaryVertex, const int pdgIdBin) : daughters_(daughters) {

  if (mother->numberOfDaughters() > 0) {
    x_ = mother->daughter(0)->vx();
    y_ = mother->daughter(0)->vy();
    z_ = mother->daughter(0)->vz();
    dxy_ = vertexntuples::dxy(mother->daughter(0), primaryVertex);
    dz_ = vertexntuples::dz(mother->daughter(0), primaryVertex);
    d3d_ = vertexntuples::d3d(mother->daughter(0), primaryVertex);
    dxyerr_ = vertexntuples::dxyErr(mother->daughter(0), primaryVertex);
    dzerr_ = vertexntuples::dzErr(mother->daughter(0), primaryVertex);
    d3derr_ = vertexntuples::d3dErr(mother->daughter(0), primaryVertex);
    dxysig_ = dxy_ / dxyerr_;
    dzsig_ = dz_ / dzerr_;
    d3dsig_ = d3d_ / d3derr_;
  } else { // This code should only be used during debug runs
    x_ = 0.0;
    y_ = 0.0;
    z_ = 0.0;
    dxy_ = -1.0;
    dz_ = -1.0;
    d3d_ = -1.0;
    dxyerr_ = -1.0;
    dzerr_ = -1.0;
    d3derr_ = -1.0;
    dxysig_ = -1.0;
    dzsig_ = -1.0;
    d3dsig_ = -1.0;
  }
  xerr_ = 0.0;
  yerr_ = 0.0;
  zerr_ = 0.0;
  pt_ = mother->pt();
  pt2_ = mother->pt() * mother->pt();
  eta_ = mother->eta();
  phi_ = mother->phi();
  motherPdgId_ = mother->pdgId();
  pdgIdBin_ = pdgIdBin;
  nDaughters_ = daughters->size();

  daughterPt_ = new std::vector<float>;
  daughterPt2_ = new std::vector<float>;
  daughterEta_ = new std::vector<float>;
  daughterPhi_ = new std::vector<float>;
  daughterPdgId_ = new std::vector<float>;
  daughterCharge_ = new std::vector<float>;
  for (const reco::Candidate* dau : *daughters) {
    daughterPt_->push_back(dau->pt());
    daughterPt2_->push_back(dau->pt() * dau->pt());
    daughterEta_->push_back(dau->eta());
    daughterPhi_->push_back(dau->phi());
    daughterPdgId_->push_back(dau->pdgId());
    daughterCharge_->push_back(dau->charge());
  }

  // For matching to SecondaryVertex tracks
  daughterDeltaR_ = new std::vector<float>;
  daughterPtResNorm_ = new std::vector<float>;
}


// GenVertex::~GenVertex() {}


void GenVertex::fill(std::map<TString, TH1F*>& histos, TString prefix) {

  for (unsigned int iDau = 0; iDau < nDaughters_; iDau++) {
    histos[prefix + "_trk_pt"]->Fill(daughterPt_->at(iDau));
    histos[prefix + "_trk_pt2"]->Fill(daughterPt2_->at(iDau));
    histos[prefix + "_trk_eta"]->Fill(daughterEta_->at(iDau));
    histos[prefix + "_trk_phi"]->Fill(daughterPhi_->at(iDau));
    histos[prefix + "_trk_charge"]->Fill(daughterCharge_->at(iDau));
  }

  histos[prefix + "_x"]->Fill(x_);
  histos[prefix + "_y"]->Fill(y_);
  histos[prefix + "_z"]->Fill(z_);
  histos[prefix + "_pt"]->Fill(pt_);
  histos[prefix + "_pt2"]->Fill(pt2_);
  histos[prefix + "_eta"]->Fill(eta_);
  histos[prefix + "_phi"]->Fill(phi_);
  histos[prefix + "_dxy"]->Fill(dxy_);
  histos[prefix + "_dz"]->Fill(dz_);
  histos[prefix + "_d3d"]->Fill(d3d_);
  histos[prefix + "_dxyerr"]->Fill(dxyerr_);
  histos[prefix + "_dzerr"]->Fill(dzerr_);
  histos[prefix + "_d3derr"]->Fill(d3derr_);
  histos[prefix + "_dxysig"]->Fill(dxysig_);
  histos[prefix + "_dzsig"]->Fill(dzsig_);
  histos[prefix + "_d3dsig"]->Fill(d3dsig_);
  histos[prefix + "_motherPdgId"]->Fill(motherPdgId_);
  histos[prefix + "_pdgIdBin"]->Fill(pdgIdBin_);
  histos[prefix + "_ntrk"]->Fill(nDaughters_);
}
