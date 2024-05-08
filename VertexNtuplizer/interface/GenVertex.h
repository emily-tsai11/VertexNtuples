#ifndef GEN_VERTEX
#define GEN_VERTEX


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TH1.h"


class GenVertex {

  public:

    GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
        const reco::Vertex& primaryVertex, const int pdgIdBin);
    // ~GenVertex();

    void fill(std::map<TString, TH1F*>& histos, TString prefix);
    void print();

    const float x() const { return daughters_->at(0)->vx(); }
    const float y() const { return daughters_->at(0)->vy(); }
    const float z() const { return daughters_->at(0)->vz(); }
    const float pt() const { return mother_->pt(); }
    const float pt2() const { return mother_->pt() * mother_->pt(); }
    const float eta() const { return mother_->eta(); }
    const float phi() const { return mother_->phi(); }
    const float dxy() const { return dxy_; }
    const float dz() const { return dz_; }
    const float d3d() const { return d3d_; }
    const float dxyErr() const { return dxyerr_; }
    const float dzErr() const { return dzerr_; }
    const float d3dErr() const { return d3derr_; }
    const float dxySig() const { return dxy_ / dxyerr_; }
    const float dzSig() const { return dz_ / dzerr_; }
    const float d3dSig() const { return d3d_ / d3derr_; }
    const int motherPdgId() const { return mother_->pdgId(); }
    const int pdgIdBin() const { return pdgIdBin_; }
    const unsigned int nDaughters() const { return daughters_->size(); }

    const reco::GenParticle* mother() const { return mother_; }
    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

  private:

    float dxy_;
    float dz_;
    float d3d_;
    float dxyerr_;
    float dzerr_;
    float d3derr_;
    float dxysig_;
    float dzsig_;
    float d3dsig_;
    int pdgIdBin_;

    const reco::GenParticle* mother_;
    std::vector<const reco::Candidate*>* daughters_;
};


#endif
