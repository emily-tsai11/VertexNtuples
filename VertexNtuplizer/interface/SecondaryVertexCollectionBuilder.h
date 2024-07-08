#ifndef VertexNtuples_VertexNtuplizer_SecondaryVertexCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_SecondaryVertexCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SecondaryVertex.h"


typedef std::vector<SecondaryVertex> SecondaryVertexCollection;


class SecondaryVertexCollectionBuilder {

  public:

    SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~SecondaryVertexCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesToken,
        edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDBSToken,
        edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDBS4Token,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeBSValueMapToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeBSErrorMapToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeBSQualityMapToken,
        edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDPVToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimePVValueMapToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimePVErrorMapToken,
        // edm::EDGetTokenT<edm::ValueMap<float>>& trackTimePVQualityMapToken,
        const reco::Vertex& primaryVertex);

    SecondaryVertexCollection getSecondaryVertexCollection() { return secondaryVertices_; }
    SecondaryVertexCollection getSecondaryVertexCollectionMTDBS() { return secondaryVerticesMTDBS_; }
    SecondaryVertexCollection getSecondaryVertexCollectionMTDBS4() { return secondaryVerticesMTDBS4_; }
    SecondaryVertexCollection getSecondaryVertexCollectionMTDPV() { return secondaryVerticesMTDPV_; }

  private:

    bool goodRecoVertex(const reco::Vertex& v);
    static bool compare(const SecondaryVertex& sva, const SecondaryVertex& svb);

    double absEtaMax_;
    double trkPtMin_;
    double trkTimeQualityCut_;
    double svChi2dofMax_;

    SecondaryVertexCollection secondaryVertices_;
    SecondaryVertexCollection secondaryVerticesMTDBS_;
    SecondaryVertexCollection secondaryVerticesMTDBS4_;
    SecondaryVertexCollection secondaryVerticesMTDPV_;

    reco::VertexCollection cmsSecondaryVertices_;
    reco::VertexCollection cmsSecondaryVerticesMTDBS_;
    reco::VertexCollection cmsSecondaryVerticesMTDBS4_;
    edm::ValueMap<float> trackTimeBSValueMap_;
    edm::ValueMap<float> trackTimeBSErrorMap_;
    edm::ValueMap<float> trackTimeBSQualityMap_;
    reco::VertexCollection cmsSecondaryVerticesMTDPV_;
    edm::ValueMap<float> trackTimePVValueMap_;
    edm::ValueMap<float> trackTimePVErrorMap_;
    // edm::ValueMap<float> trackTimePVQualityMap_;
};


#endif
