#include "../interface/GenVertex.h"


GenVertex::GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters) :
  mother_(mother), daughters_(daughters) {}


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
