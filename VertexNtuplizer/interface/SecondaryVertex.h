#ifndef VertexNtuples_VertexNtuplizer_SecondaryVertex_h
#define VertexNtuples_VertexNtuplizer_SecondaryVertex_h


#include <algorithm>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TH1.h"

#include "VertexCalculator.h"


class SecondaryVertex {

  public:

    SecondaryVertex(const reco::Vertex& sv,
        const std::vector<reco::TrackBaseRef>* tracks,
        const reco::Vertex& primaryVertex,
        const edm::ValueMap<float>& trackTimeValueMap,
        const edm::ValueMap<float>& trackTimeErrorMap,
        const edm::ValueMap<float>& trackTimeQualityMap);
    // ~SecondaryVertex();

    void fill(std::map<TString, TH1F*>& histos, TString prefix);

    const std::vector<double>* trkTVal() const { return trk_tval_; }
    const std::vector<double>* trkTErr() const { return trk_terr_; }
    const std::vector<double>* trkTSig() const { return trk_tsig_; }
    const std::vector<double>* trkTQual() const { return trk_tqual_; }
    const std::vector<double>* trkPt() const { return trk_pt_; }
    const std::vector<double>* trkPt2() const { return trk_pt2_; }
    const std::vector<double>* trkEta() const { return trk_eta_; }
    const std::vector<double>* trkPhi() const { return trk_phi_; }
    const std::vector<double>* trkDxy() const { return trk_dxy_; }
    const std::vector<double>* trkDz() const { return trk_dz_; }
    const std::vector<double>* trkD3d() const { return trk_d3d_; }
    const std::vector<double>* trkDxyErr() const { return trk_dxyerr_; }
    const std::vector<double>* trkDzErr() const { return trk_dzerr_; }
    const std::vector<double>* trkD3dErr() const { return trk_d3derr_; }
    const std::vector<double>* trkDxySig() const { return trk_dxysig_; }
    const std::vector<double>* trkDzSig() const { return trk_dzsig_; }
    const std::vector<double>* trkD3dSig() const { return trk_d3dsig_; }
    const std::vector<double>* trkCharge() const { return trk_charge_; }
    const std::vector<double>* trkChi2() const { return trk_chi2_; }
    const std::vector<double>* trkNDOF() const { return trk_ndof_; }
    const std::vector<double>* trkChi2DOF() const { return trk_chi2dof_; }

    const double x() const { return x_; }
    const double y() const { return y_; }
    const double z() const { return z_; }
    const double xErr() const { return xerr_; }
    const double yErr() const { return yerr_; }
    const double zErr() const { return zerr_; }
    const double dxy() const { return dxy_; }
    const double dz() const { return dz_; }
    const double d3d() const { return d3d_; }
    const double dxyErr() const { return dxyerr_; }
    const double dzErr() const { return dzerr_; }
    const double d3dErr() const { return d3derr_; }
    const double dxySig() const { return dxy_ / dxyerr_; }
    const double dzSig() const { return dz_ / dzerr_; }
    const double d3dSig() const { return d3d_ / d3derr_; }
    const double pt() const { return pt_; }
    const double pt2() const { return pt2_; }
    const double eta() const { return eta_; }
    const double phi() const { return phi_; }
    const double tAvg() const { return tavg_; }
    const double tRange() const { return trange_; }
    const double chi2() const { return chi2_; }
    const double nDOF() const { return ndof_; }
    const double chi2DOF() const { return chi2dof_; }
    const double nTracks() const { return ntrk_; }

  private:

    std::vector<double>* trk_tval_;
    std::vector<double>* trk_terr_;
    std::vector<double>* trk_tsig_;
    std::vector<double>* trk_tqual_;
    std::vector<double>* trk_pt_;
    std::vector<double>* trk_pt2_;
    std::vector<double>* trk_eta_;
    std::vector<double>* trk_phi_;
    std::vector<double>* trk_dxy_;
    std::vector<double>* trk_dz_;
    std::vector<double>* trk_d3d_;
    std::vector<double>* trk_dxyerr_;
    std::vector<double>* trk_dzerr_;
    std::vector<double>* trk_d3derr_;
    std::vector<double>* trk_dxysig_;
    std::vector<double>* trk_dzsig_;
    std::vector<double>* trk_d3dsig_;
    std::vector<double>* trk_charge_;
    std::vector<double>* trk_chi2_;
    std::vector<double>* trk_ndof_;
    std::vector<double>* trk_chi2dof_;

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
    double pt_;
    double pt2_;
    double eta_;
    double phi_;
    double tavg_;
    double trange_;
    double chi2_;
    double ndof_;
    double chi2dof_;
    unsigned int ntrk_;
};


#endif
