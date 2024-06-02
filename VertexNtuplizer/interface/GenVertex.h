#ifndef VertexNtuples_VertexNtuplizer_GenVertex_h
#define VertexNtuples_VertexNtuplizer_GenVertex_h


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TH1.h"

#include "VertexCalculator.h"


class GenVertex {

  public:

    GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
        const reco::Vertex& primaryVertex, const int pdgIdBin);
    GenVertex(const HepMC::FourVector& position, const HepMC::GenParticle* mother,
        std::vector<const HepMC::GenParticle*>* daughters,
        const reco::Vertex& primaryVertex, const int pdgIdBin);
    // ~GenVertex();

    void fill(std::map<TString, TH1F*>& histos, TString prefix);
    void print();

    void addDeltaR(double deltaR) { daughterDeltaR_->push_back(deltaR); }
    void addPtResNorm(double ptresnorm) { daughterPtResNorm_->push_back(ptresnorm); }

    const double x() const { return x_; }
    const double y() const { return y_; }
    const double z() const { return z_; }
    const double xErr() const { return xerr_; }
    const double yErr() const { return yerr_; }
    const double zErr() const { return zerr_; }
    const double pt() const { return pt_; }
    const double pt2() const { return pt2_; }
    const double eta() const { return eta_; }
    const double phi() const { return phi_; }
    const double dxy() const { return dxy_; }
    const double dz() const { return dz_; }
    const double d3d() const { return d3d_; }
    const double dxyErr() const { return dxyerr_; }
    const double dzErr() const { return dzerr_; }
    const double d3dErr() const { return d3derr_; }
    const double dxySig() const { return dxysig_; }
    const double dzSig() const { return dzsig_; }
    const double d3dSig() const { return d3dsig_; }
    const int motherPdgId() const { return motherPdgId_; }
    const int pdgIdBin() const { return pdgIdBin_; }
    const unsigned int nDaughters() const { return nDaughters_; }

    const std::vector<double>* daughterPt() const { return daughterPt_; }
    const std::vector<double>* daughterPt2() const { return daughterPt2_; }
    const std::vector<double>* daughterEta() const { return daughterEta_; }
    const std::vector<double>* daughterPhi() const { return daughterPhi_; }
    const std::vector<double>* daughterCharge() const { return daughterCharge_; }

    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }

    const std::vector<double>* dauMatchDeltaR() const { return daughterDeltaR_; }
    const std::vector<double>* dauMatchPtResNorm() const { return daughterPtResNorm_; }

  private:

    double x_;
    double y_;
    double z_;
    double xerr_;
    double yerr_;
    double zerr_;
    double pt_;
    double pt2_;
    double eta_;
    double phi_;
    double dxy_;
    double dz_;
    double d3d_;
    double dxyerr_;
    double dzerr_;
    double d3derr_;
    double dxysig_;
    double dzsig_;
    double d3dsig_;
    int motherPdgId_;
    int pdgIdBin_;
    unsigned int nDaughters_;

    std::vector<double>* daughterPt_;
    std::vector<double>* daughterPt2_;
    std::vector<double>* daughterEta_;
    std::vector<double>* daughterPhi_;
    std::vector<double>* daughterPdgId_;
    std::vector<double>* daughterCharge_;

    std::vector<const reco::Candidate*>* daughters_;

    // Filled when matching to SecondaryVertex tracks
    std::vector<double>* daughterDeltaR_;
    std::vector<double>* daughterPtResNorm_;
};


#endif
