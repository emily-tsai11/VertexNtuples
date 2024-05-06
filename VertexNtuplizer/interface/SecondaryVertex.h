#ifndef SECONDARY_VERTEX
#define SECONDARY_VERTEX


#include <algorithm>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TMath.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"


class SecondaryVertex {

  public:

    SecondaryVertex(const reco::Vertex& sv,
        const std::vector<reco::TrackBaseRef>* tracks,
        const reco::Vertex& primaryVertex,
        const edm::ValueMap<float>& trackTimeValueMap,
        const edm::ValueMap<float>& trackTimeErrorMap,
        const edm::ValueMap<float>& trackTimeQualityMap);
    // ~SecondaryVertex();

  private:

    Measurement1D getDxy(reco::Vertex sv, reco::Vertex pv);
    Measurement1D getD3d(reco::Vertex sv, reco::Vertex pv);

    std::vector<float>* trk_tval_;
    std::vector<float>* trk_terr_;
    std::vector<float>* trk_tsig_;
    std::vector<float>* trk_tqal_;
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
    std::vector<float>* trk_chi2_;
    std::vector<float>* trk_ndof_;
    std::vector<float>* trk_chi2dof_;

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
    float tavg_;
    float trng_;
    float chi2_;
    float ndof_;
    float chi2dof_;
    unsigned int ntrk_;
};


#endif
