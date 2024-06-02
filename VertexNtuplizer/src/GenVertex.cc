#include "../interface/GenVertex.h"


GenVertex::GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
    const reco::Vertex& primaryVertex, const int pdgIdBin) : daughters_(daughters) {

  x_ = daughters->at(0)->vx();
  y_ = daughters->at(0)->vy();
  z_ = daughters->at(0)->vz();
  xerr_ = 0.0;
  yerr_ = 0.0;
  zerr_ = 0.0;
  pt_ = mother->pt();
  pt2_ = mother->pt() * mother->pt();
  eta_ = mother->eta();
  phi_ = mother->phi();
  dxy_ = vertexntuples::dxy(daughters->at(0), primaryVertex);
  dz_ = vertexntuples::dz(daughters->at(0), primaryVertex);
  d3d_ = vertexntuples::d3d(daughters->at(0), primaryVertex);
  dxyerr_ = vertexntuples::dxyErr(daughters->at(0), primaryVertex);
  dzerr_ = vertexntuples::dzErr(daughters->at(0), primaryVertex);
  d3derr_ = vertexntuples::d3dErr(daughters->at(0), primaryVertex);
  dxysig_ = dxy_ / dxyerr_;
  dzsig_ = dz_ / dzerr_;
  d3dsig_ = d3d_ / d3derr_;
  motherPdgId_ = mother->pdgId();
  pdgIdBin_ = pdgIdBin;
  nDaughters_ = daughters->size();

  daughterPt_ = new std::vector<double>;
  daughterPt2_ = new std::vector<double>;
  daughterEta_ = new std::vector<double>;
  daughterPhi_ = new std::vector<double>;
  daughterPdgId_ = new std::vector<double>;
  daughterCharge_ = new std::vector<double>;
  for (const reco::Candidate* dau : *daughters) {
    daughterPt_->push_back(dau->pt());
    daughterPt2_->push_back(dau->pt() * dau->pt());
    daughterEta_->push_back(dau->eta());
    daughterPhi_->push_back(dau->phi());
    daughterPdgId_->push_back(dau->pdgId());
    daughterCharge_->push_back(dau->charge());
  }

  // For matching to SecondaryVertex tracks
  daughterDeltaR_ = new std::vector<double>;
  daughterPtResNorm_ = new std::vector<double>;
}


GenVertex::GenVertex(const HepMC::FourVector& position, const HepMC::GenParticle* mother,
    std::vector<const HepMC::GenParticle*>* daughters, const reco::Vertex& primaryVertex,
    const int pdgIdBin) {

  double px = mother->momentum().px();
  double py = mother->momentum().py();
  double pt = TMath::Sqrt(px * px + py * py);

  x_ = position.x();
  y_ = position.y();
  z_ = position.z();
  xerr_ = 0.0;
  yerr_ = 0.0;
  zerr_ = 0.0;
  pt_ = pt;
  pt2_ = pt * pt;
  eta_ = mother->momentum().eta();
  phi_ = mother->momentum().phi();
  dxy_ = vertexntuples::dxy(position, primaryVertex);
  dz_ = vertexntuples::dz(position, primaryVertex);
  d3d_ = vertexntuples::d3d(position, primaryVertex);
  dxyerr_ = vertexntuples::dxyErr(position, primaryVertex);
  dzerr_ = vertexntuples::dzErr(position, primaryVertex);
  d3derr_ = vertexntuples::d3dErr(position, primaryVertex);
  dxysig_ = dxy_ / dxyerr_;
  dzsig_ = dz_ / dzerr_;
  d3dsig_ = d3d_ / d3derr_;
  motherPdgId_ = mother->pdg_id();
  pdgIdBin_ = pdgIdBin;
  nDaughters_ = daughters->size();

  daughterPt_ = new std::vector<double>;
  daughterPt2_ = new std::vector<double>;
  daughterEta_ = new std::vector<double>;
  daughterPhi_ = new std::vector<double>;
  daughterPdgId_ = new std::vector<double>;
  daughterCharge_ = new std::vector<double>;
  for (const HepMC::GenParticle* dau : *daughters) {
    px = dau->momentum().px();
    py = dau->momentum().py();
    pt = TMath::Sqrt(px * px + py * py);
    daughterPt_->push_back(pt);
    daughterPt2_->push_back(pt * pt);
    daughterEta_->push_back(dau->momentum().eta());
    daughterPhi_->push_back(dau->momentum().phi());
    daughterPdgId_->push_back(dau->pdg_id());
    daughterCharge_->push_back(-1000); // Seems to not have a charge accessor for HepMC::GenParticle
  }

  // For matching to SecondaryVertex tracks
  daughterDeltaR_ = new std::vector<double>;
  daughterPtResNorm_ = new std::vector<double>;
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


void GenVertex::print() {

  std::cout << "GenVertex:" << std::endl;
  std::cout << "    vertex x          = " << x_ << std::endl;
  std::cout << "    vertex y          = " << y_ << std::endl;
  std::cout << "    vertex z          = " << z_ << std::endl;
  std::cout << "    vertex pt         = " << pt_ << std::endl;
  std::cout << "    vertex pt^2       = " << pt2_ << std::endl;
  std::cout << "    vertex eta        = " << eta_ << std::endl;
  std::cout << "    vertex phi        = " << phi_ << std::endl;
  std::cout << "    nDaughters        = " << nDaughters_ << std::endl;
  std::cout << "    mother pdg id     = " << motherPdgId_ << std::endl;
  std::cout << "    mother pdg id bin = " << pdgIdBin_ << std::endl;

  std::cout << "    daughter pdgIds = ";
  for (unsigned int iDau = 0; iDau < nDaughters_; iDau++) {
    std::cout << daughterPdgId_->at(iDau);
    if (iDau < nDaughters_ - 1) std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << "    daughter pts   = ";
  for (unsigned int iDau = 0; iDau < nDaughters_; iDau++) {
    std::cout << daughterPt_->at(iDau);
    if (iDau < nDaughters_ - 1) std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << "    daughter etas   = ";
  for (unsigned int iDau = 0; iDau < nDaughters_; iDau++) {
    std::cout << daughterEta_->at(iDau);
    if (iDau < nDaughters_ - 1) std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << "    daughter phis   = ";
  for (unsigned int iDau = 0; iDau < nDaughters_; iDau++) {
    std::cout << daughterPhi_->at(iDau);
    if (iDau < nDaughters_ - 1) std::cout << ", ";
  }
  std::cout << std::endl;
}
