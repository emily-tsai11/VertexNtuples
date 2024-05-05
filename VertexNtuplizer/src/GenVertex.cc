#include "../interface/GenVertex.h"


GenVertex::GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
    const reco::Vertex primaryVertex) : mother_(mother), daughters_(daughters) {

  float dx = primaryVertex.x() - daughters->at(0)->vx();
  float dy = primaryVertex.y() - daughters->at(0)->vy();
  float dxerr = primaryVertex.xError();
  float dyerr = primaryVertex.yError();
  float dx2 = dx*dx;
  float dy2 = dy*dy;
  float dx2err = 2*dx*dxerr;
  float dy2err = 2*dy*dyerr;
  float dxy2 = dx2 + dy2;
  float dxy2err = TMath::Sqrt(dx2err*dx2err + dy2err*dy2err);
  dxy_ = TMath::Sqrt(dxy2);
  dxyerr_ = 0.5 * dxy2err / dxy_;

  dz_ = primaryVertex.z() - daughters->at(0)->vz();
  dzerr_ = primaryVertex.zError();

  float dz2 = dz_*dz_;
  float dz2err = 2*dz_*dzerr_;
  float d3d2 = dxy2 + dz2;
  float d3d2err = TMath::Sqrt(dxy2err*dxy2err + dz2err*dz2err);
  d3d_ = TMath::Sqrt(d3d2);
  d3derr_ = 0.5 * d3d2err / d3d_;
}


// GenVertex::~GenVertex() {}


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
