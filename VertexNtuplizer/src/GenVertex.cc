#include "../interface/GenVertex.h"


GenVertex::GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
    const reco::Vertex& primaryVertex, const int pdgIdBin) : mother_(mother), daughters_(daughters) {

  float dx = abs(primaryVertex.x() - daughters->at(0)->vx());
  float dy = abs(primaryVertex.y() - daughters->at(0)->vy());
  float dxerr = primaryVertex.xError();
  float dyerr = primaryVertex.yError();
  float dx2 = dx*dx;
  float dy2 = dy*dy;
  float dx2err = 2*dx*dxerr;
  float dy2err = 2*dy*dyerr;
  float dxy2 = dx2 + dy2;
  float dxy2err = TMath::Sqrt(dx2err*dx2err + dy2err*dy2err);
  float dxy = TMath::Sqrt(dxy2);
  float dxyerr = 0.5 * dxy2err / dxy;

  float dz = abs(primaryVertex.z() - daughters->at(0)->vz());
  float dzerr = primaryVertex.zError();
  float dz2 = dz*dz;
  float dz2err = 2*dz*dzerr;
  float d3d2 = dxy2 + dz2;
  float d3d2err = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
  float d3d = TMath::Sqrt(d3d2);
  float d3derr = 0.5 * d3d2err / d3d;

  dxy_ = dxy;
  dxyerr_ = dxyerr;
  dz_ = dz;
  dzerr_ = dzerr;
  d3d_ = d3d;
  d3derr_ = d3derr;
  pdgIdBin_ = pdgIdBin;
}


// GenVertex::~GenVertex() {}


void GenVertex::fill(std::map<TString, TH1F*>& histos, TString prefix) {

  for (const reco::Candidate* dau : *daughters_) {
    histos[prefix + "_trk_pt"]->Fill(dau->pt());
    histos[prefix + "_trk_pt2"]->Fill(dau->pt() * dau->pt());
    histos[prefix + "_trk_eta"]->Fill(dau->eta());
    histos[prefix + "_trk_phi"]->Fill(dau->phi());
    histos[prefix + "_trk_charge"]->Fill(dau->charge());
  }

  histos[prefix + "_x"]->Fill(x());
  histos[prefix + "_y"]->Fill(y());
  histos[prefix + "_z"]->Fill(z());
  histos[prefix + "_pt"]->Fill(pt());
  histos[prefix + "_pt2"]->Fill(pt2());
  histos[prefix + "_eta"]->Fill(eta());
  histos[prefix + "_phi"]->Fill(phi());
  histos[prefix + "_dxy"]->Fill(dxy());
  histos[prefix + "_dz"]->Fill(dz());
  histos[prefix + "_d3d"]->Fill(d3d());
  histos[prefix + "_dxyerr"]->Fill(dxyErr());
  histos[prefix + "_dzerr"]->Fill(dzErr());
  histos[prefix + "_d3derr"]->Fill(d3dErr());
  histos[prefix + "_dxysig"]->Fill(dxySig());
  histos[prefix + "_dzsig"]->Fill(dzSig());
  histos[prefix + "_d3dsig"]->Fill(d3dSig());
  histos[prefix + "_motherPdgId"]->Fill(motherPdgId());
  histos[prefix + "_pdgIdBin"]->Fill(pdgIdBin());
  histos[prefix + "_ntrk"]->Fill(nDaughters());
}


void GenVertex::print() {

  std::cout << "GenVertex:" << std::endl;
  std::cout << "    vertex x        = " << x() << std::endl;
  std::cout << "    vertex y        = " << y() << std::endl;
  std::cout << "    vertex z        = " << z() << std::endl;
  std::cout << "    vertex pt       = " << pt() << std::endl;
  std::cout << "    vertex eta      = " << eta() << std::endl;
  std::cout << "    vertex phi      = " << phi() << std::endl;
  std::cout << "    nDaughters      = " << nDaughters() << std::endl;
  std::cout << "    mother pdg id   = " << motherPdgId() << std::endl;

  std::cout << "    daughter pdgIds = ";
  for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
    std::cout << daughters_->at(iDau)->pdgId();
    if (iDau < nDaughters() - 1) std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << "    daughter pts   = ";
  for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
    std::cout << daughters_->at(iDau)->pt();
    if (iDau < nDaughters() - 1) std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << "    daughter etas   = ";
  for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
    std::cout << daughters_->at(iDau)->eta();
    if (iDau < nDaughters() - 1) std::cout << ", ";
  }
  std::cout << std::endl;
  std::cout << "    daughter phis   = ";
  for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
    std::cout << daughters_->at(iDau)->phi();
    if (iDau < nDaughters() - 1) std::cout << ", ";
  }
  std::cout << std::endl;
}
