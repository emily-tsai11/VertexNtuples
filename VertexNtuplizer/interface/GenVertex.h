#ifndef GEN_VERTEX
#define GEN_VERTEX


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TMath.h"


class GenVertex {

  public:

    GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
        const reco::Vertex& primaryVertex);
    // ~GenVertex();

    const float x()        const { return daughters_->at(0)->vx(); }
    const float y()        const { return daughters_->at(0)->vy(); }
    const float z()        const { return daughters_->at(0)->vz(); }
    const float pt()       const { return mother_->pt(); }
    const float eta()      const { return mother_->eta(); }
    const float phi()      const { return mother_->phi(); }
    const float dxy()      const { return dxy_; }
    const float dz()       const { return dz_; }
    const float d3d()      const { return d3d_; }
    const float dxyError() const { return dxyerr_; }
    const float dzError()  const { return dzerr_; }
    const float d3dError() const { return d3derr_; }
    const unsigned int nDaughters() const { return daughters_->size(); }
    const int motherPdgId()         const { return mother_->pdgId(); }

    const reco::GenParticle* mother() const { return mother_; }
    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

    void print();

  private:

    const reco::GenParticle* mother_;
    std::vector<const reco::Candidate*>* daughters_;

    float dxy_;
    float dz_;
    float d3d_;
    float dxyerr_;
    float dzerr_;
    float d3derr_;
};


#endif
