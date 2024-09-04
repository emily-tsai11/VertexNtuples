#ifndef VertexNtuples_VertexNtuplizer_GenVertex_h
#define VertexNtuples_VertexNtuplizer_GenVertex_h


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TH1.h"

#include "VertexMatcher.h"


class GenVertex {

  public:

    GenVertex(const reco::Candidate* mother, std::vector<const reco::Candidate*>* daughters,
        const reco::Vertex& primaryVertex, const int pdgIdBin);
    // ~GenVertex();

    void fill(std::map<TString, TH1F*>& histos, TString prefix);

    void addNormPtRes(float normPtRes) { dauMatchNormPtRes_->push_back(normPtRes); }
    void addDeltaR(float deltaR) { dauMatchDeltaR_->push_back(deltaR); }

    const float x() const { return x_; }
    const float y() const { return y_; }
    const float z() const { return z_; }
    const float xErr() const { return xerr_; }
    const float yErr() const { return yerr_; }
    const float zErr() const { return zerr_; }
    const float pt() const { return pt_; }
    const float pt2() const { return pt2_; }
    const float eta() const { return eta_; }
    const float phi() const { return phi_; }
    const float dxy() const { return dxy_; }
    const float dz() const { return dz_; }
    const float d3d() const { return d3d_; }
    const float dxyErr() const { return dxyerr_; }
    const float dzErr() const { return dzerr_; }
    const float d3dErr() const { return d3derr_; }
    const float dxySig() const { return dxysig_; }
    const float dzSig() const { return dzsig_; }
    const float d3dSig() const { return d3dsig_; }
    const int motherPdgId() const { return motherPdgId_; }
    const int pdgIdBin() const { return pdgIdBin_; }
    const unsigned int nDaughters() const { return nDaughters_; }

    const std::vector<float>* daughterPt() const { return daughterPt_; }
    const std::vector<float>* daughterPt2() const { return daughterPt2_; }
    const std::vector<float>* daughterEta() const { return daughterEta_; }
    const std::vector<float>* daughterPhi() const { return daughterPhi_; }
    const std::vector<float>* daughterCharge() const { return daughterCharge_; }

    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

    const std::vector<float>* dauMatchNormPtRes() const { return dauMatchNormPtRes_; }
    const std::vector<float>* dauMatchDeltaR() const { return dauMatchDeltaR_; }

  private:

    float x_;
    float y_;
    float z_;
    float xerr_;
    float yerr_;
    float zerr_;
    float pt_;
    float pt2_;
    float eta_;
    float phi_;
    float dxy_;
    float dz_;
    float d3d_;
    float dxyerr_;
    float dzerr_;
    float d3derr_;
    float dxysig_;
    float dzsig_;
    float d3dsig_;
    int motherPdgId_;
    int pdgIdBin_;
    unsigned int nDaughters_;

    // TODO: ADD MORE PARAMETERS?
    std::vector<float>* daughterPt_;
    std::vector<float>* daughterPt2_;
    std::vector<float>* daughterEta_;
    std::vector<float>* daughterPhi_;
    std::vector<float>* daughterPdgId_;
    std::vector<float>* daughterCharge_;

    std::vector<const reco::Candidate*>* daughters_;

    // Filled when matching to SecondaryVertex tracks
    std::vector<float>* dauMatchNormPtRes_;
    std::vector<float>* dauMatchDeltaR_;
};


#endif
