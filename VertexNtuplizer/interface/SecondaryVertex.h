#ifndef VertexNtuples_VertexNtuplizer_SecondaryVertex_h
#define VertexNtuples_VertexNtuplizer_SecondaryVertex_h


#include <algorithm>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"

#include "VertexCalculator.h"


class SecondaryVertex {

  public:

    SecondaryVertex(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& primaryVertex);
    SecondaryVertex(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& primaryVertex,
        edm::Handle<reco::TrackCollection> generalTracksHandle,
        const edm::ValueMap<float>& trackT0,
        const edm::ValueMap<float>& trackSigmaT0
        // , const edm::ValueMap<float>& trackQuality
    );
    // ~SecondaryVertex();

    void initialize(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& primaryVertex,
        bool hasTime = true);
    void initializeTime(const reco::VertexCompositePtrCandidate& sv,
        edm::Handle<reco::TrackCollection> generalTracksHandle,
        const edm::ValueMap<float>& trackT0,
        const edm::ValueMap<float>& trackSigmaT0
        // , const edm::ValueMap<float>& trackQuality
    );

    void fill(std::map<TString, TH1F*>& histos, std::map<TString, TH2F*>& histos2, TString prefix);

    void addDeltaR(float deltaR) { trk_deltaR_->push_back(deltaR); }
    void addNormPtRes(float normPtRes) { trk_normPtRes_->push_back(normPtRes); }

    const std::vector<float>* trkPt() const { return trk_pt_; }
    const std::vector<float>* trkPt2() const { return trk_pt2_; }
    const std::vector<float>* trkEta() const { return trk_eta_; }
    const std::vector<float>* trkPhi() const { return trk_phi_; }
    const std::vector<float>* trkDxy() const { return trk_dxy_; }
    const std::vector<float>* trkDz() const { return trk_dz_; }
    const std::vector<float>* trkD3d() const { return trk_d3d_; }
    const std::vector<float>* trkDxyErr() const { return trk_dxyerr_; }
    const std::vector<float>* trkDzErr() const { return trk_dzerr_; }
    const std::vector<float>* trkD3dErr() const { return trk_d3derr_; }
    const std::vector<float>* trkDxySig() const { return trk_dxysig_; }
    const std::vector<float>* trkDzSig() const { return trk_dzsig_; }
    const std::vector<float>* trkD3dSig() const { return trk_d3dsig_; }
    const std::vector<float>* trkCharge() const { return trk_charge_; }
    const std::vector<float>* trkChi2() const { return trk_chi2_; }
    const std::vector<float>* trkNDOF() const { return trk_ndof_; }
    const std::vector<float>* trkChi2DOF() const { return trk_chi2dof_; }
    const std::vector<float>* trkTVal() const { return trk_tval_; }
    const std::vector<float>* trkTErr() const { return trk_terr_; }
    const std::vector<float>* trkTSig() const { return trk_tsig_; }
    // const std::vector<float>* trkTQual() const { return trk_tqual_; }

    const std::vector<float>* trkMatchDeltaR() const { return trk_deltaR_; }
    const std::vector<float>* trkMatchNormPtRes() const { return trk_normPtRes_; }

    const float x() const { return x_; }
    const float y() const { return y_; }
    const float z() const { return z_; }
    const float xErr() const { return xerr_; }
    const float yErr() const { return yerr_; }
    const float zErr() const { return zerr_; }
    const float dxy() const { return dxy_; }
    const float dz() const { return dz_; }
    const float d3d() const { return d3d_; }
    const float dxyErr() const { return dxyerr_; }
    const float dzErr() const { return dzerr_; }
    const float d3dErr() const { return d3derr_; }
    const float dxySig() const { return dxysig_; }
    const float dzSig() const { return dzsig_; }
    const float d3dSig() const { return d3dsig_; }
    const float pt() const { return pt_; }
    const float pt2() const { return pt2_; }
    const float eta() const { return eta_; }
    const float phi() const { return phi_; }
    const float chi2() const { return chi2_; }
    const float nDOF() const { return ndof_; }
    const float chi2DOF() const { return chi2dof_; }
    const float nTracks() const { return ntrk_; }
    const float tAvg() const { return tavg_; }
    const float tRange() const { return trange_; }

  private:

    std::vector<float>* trk_pt_;
    std::vector<float>* trk_pt2_;
    std::vector<float>* trk_eta_;
    std::vector<float>* trk_phi_;
    std::vector<float>* trk_dxy_;
    std::vector<float>* trk_dz_;
    std::vector<float>* trk_d3d_;
    std::vector<float>* trk_dxyerr_;
    std::vector<float>* trk_dzerr_;
    std::vector<float>* trk_d3derr_;
    std::vector<float>* trk_dxysig_;
    std::vector<float>* trk_dzsig_;
    std::vector<float>* trk_d3dsig_;
    std::vector<float>* trk_charge_;
    std::vector<float>* trk_chi2_;
    std::vector<float>* trk_ndof_;
    std::vector<float>* trk_chi2dof_;
    std::vector<float>* trk_tval_;
    std::vector<float>* trk_terr_;
    std::vector<float>* trk_tsig_;
    // std::vector<float>* trk_tqual_;

    float x_;
    float y_;
    float z_;
    float xerr_;
    float yerr_;
    float zerr_;
    float dxy_;
    float dz_;
    float d3d_;
    float dxyerr_;
    float dzerr_;
    float d3derr_;
    float dxysig_;
    float dzsig_;
    float d3dsig_;
    float pt_;
    float pt2_;
    float eta_;
    float phi_;
    float chi2_;
    float ndof_;
    float chi2dof_;
    unsigned int ntrk_;
    float tavg_;
    float trange_;

    // Filled when matching to GenVertex
    std::vector<float>* trk_deltaR_;
    std::vector<float>* trk_normPtRes_;
};


#endif
