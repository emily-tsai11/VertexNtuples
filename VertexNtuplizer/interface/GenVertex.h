#ifndef GEN_VERTEX
#define GEN_VERTEX


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class GenVertex {

  public:

    GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters);
    // ~GenVertex();

    const float x()   const { return daughters_->at(0)->vx(); }
    const float y()   const { return daughters_->at(0)->vy(); }
    const float z()   const { return daughters_->at(0)->vz(); }
    const float pt()  const { return mother_->pt(); }
    const float eta() const { return mother_->eta(); }
    const float phi() const { return mother_->phi(); }
    const unsigned int nDaughters() const { return daughters_->size(); }
    const int motherPdgId()         const { return mother_->pdgId(); }

    const reco::GenParticle* mother() const { return mother_; }
    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

    void print();

  private:

    const reco::GenParticle* mother_;
    std::vector<const reco::Candidate*>* daughters_;
};


#endif
