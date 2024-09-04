#ifndef VertexNtuples_VertexNtuplizer_SecondaryVertexCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_SecondaryVertexCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "SecondaryVertex.h"


typedef std::vector<SecondaryVertex> SecondaryVertexCollection;


class SecondaryVertexCollectionBuilder {

  public:

    SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~SecondaryVertexCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& inclusiveSecondaryVerticesToken,
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& inclusiveSecondaryVerticesMTDPVToken,
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& mergedSecondaryVerticesToken,
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& mergedSecondaryVerticesMTDPVToken,
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& slimmedSecondaryVerticesToken,
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& slimmedSecondaryVerticesMTDPVToken,
        edm::EDGetTokenT<reco::TrackCollection>& generalTracksToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackT0FromPVToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackSigmaT0FromPVToken,
        // edm::EDGetTokenT<edm::ValueMap<float>>& trackQualityFromPVToken,
        const reco::Vertex& primaryVertex);

    SecondaryVertexCollection getSecVerticesInclusive() { return secVerticesInclusive_; }
    SecondaryVertexCollection getSecVerticesInclusiveMTDPV() { return secVerticesInclusiveMTDPV_; }
    SecondaryVertexCollection getSecVerticesMerged() { return secVerticesMerged_; }
    SecondaryVertexCollection getSecVerticesMergedMTDPV() { return secVerticesMergedMTDPV_; }
    SecondaryVertexCollection getSecVerticesSlimmed() { return secVerticesSlimmed_; }
    SecondaryVertexCollection getSecVerticesSlimmedMTDPV() { return secVerticesSlimmedMTDPV_; }

  private:

    bool goodRecoVertex(const reco::VertexCompositePtrCandidate& v);
    // static bool compare(const SecondaryVertex& sva, const SecondaryVertex& svb);

    double trkPtMin_;
    double absEtaMax_;
    double trkTimeQualityCut_;
    double svChi2dofMax_;

    SecondaryVertexCollection secVerticesInclusive_;
    SecondaryVertexCollection secVerticesInclusiveMTDPV_;
    SecondaryVertexCollection secVerticesMerged_;
    SecondaryVertexCollection secVerticesMergedMTDPV_;
    SecondaryVertexCollection secVerticesSlimmed_;
    SecondaryVertexCollection secVerticesSlimmedMTDPV_;

    edm::Handle<reco::TrackCollection> generalTracksHandle_;

    reco::VertexCompositePtrCandidateCollection inclusiveSecondaryVertices_;
    reco::VertexCompositePtrCandidateCollection inclusiveSecondaryVerticesMTDPV_;
    reco::VertexCompositePtrCandidateCollection mergedSecondaryVertices_;
    reco::VertexCompositePtrCandidateCollection mergedSecondaryVerticesMTDPV_;
    reco::VertexCompositePtrCandidateCollection slimmedSecondaryVertices_;
    reco::VertexCompositePtrCandidateCollection slimmedSecondaryVerticesMTDPV_;
    edm::ValueMap<float> trackT0FromPV_;
    edm::ValueMap<float> trackSigmaT0FromPV_;
    // edm::ValueMap<float> trackQualityFromPV_;
};


#endif
