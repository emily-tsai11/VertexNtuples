#ifndef VertexNtuples_VertexNtuplizer_GenVertex_h
#define VertexNtuples_VertexNtuplizer_GenVertex_h


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TH1.h"

#include "VertexCalculator.h"


class GenVertex {

  public:

    GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
        const reco::Vertex& primaryVertex, const int pdgIdBin);
    // ~GenVertex();

    void fill(std::map<TString, TH1F*>& histos, TString prefix);
    void print();

    void addDeltaR(double deltaR) { dau_deltaR_->push_back(deltaR); }
    void addPtResNorm(double ptresnorm) { dau_ptresnorm_->push_back(ptresnorm); }

    const std::vector<double>* dauMatchDeltaR() const { return dau_deltaR_; }
    const std::vector<double>* dauMatchPtResNorm() const { return dau_ptresnorm_; }

    const double x() const { return x_; }
    const double y() const { return y_; }
    const double z() const { return z_; }
    const double xErr() const { return xerr_; }
    const double yErr() const { return yerr_; }
    const double zErr() const { return zerr_; }
    const double pt() const { return mother_->pt(); }
    const double pt2() const { return mother_->pt() * mother_->pt(); }
    const double eta() const { return mother_->eta(); }
    const double phi() const { return mother_->phi(); }
    const double dxy() const { return dxy_; }
    const double dz() const { return dz_; }
    const double d3d() const { return d3d_; }
    const double dxyErr() const { return dxyerr_; }
    const double dzErr() const { return dzerr_; }
    const double d3dErr() const { return d3derr_; }
    const double dxySig() const { return dxysig_; }
    const double dzSig() const { return dzsig_; }
    const double d3dSig() const { return d3dsig_; }
    const int motherPdgId() const { return mother_->pdgId(); }
    const int pdgIdBin() const { return pdgIdBin_; }
    const unsigned int nDaughters() const { return daughters_->size(); }

    const reco::GenParticle* mother() const { return mother_; }
    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

  private:

    double x_;
    double y_;
    double z_;
    double xerr_;
    double yerr_;
    double zerr_;
    double dxy_;
    double dz_;
    double d3d_;
    double dxyerr_;
    double dzerr_;
    double d3derr_;
    double dxysig_;
    double dzsig_;
    double d3dsig_;
    int pdgIdBin_;

    // Filled when matching to SecondaryVertex
    std::vector<double>* dau_deltaR_;
    std::vector<double>* dau_ptresnorm_;

    const reco::GenParticle* mother_;
    std::vector<const reco::Candidate*>* daughters_;
};


#endif
