#ifndef GEN_VERTEX_COLLECTION_BUILDER
#define GEN_VERTEX_COLLECTION_BUILDER


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "GenVertex.h"


typedef std::vector<GenVertex> GenVertexCollection;


class GenVertexCollectionBuilder {

  public:

    GenVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~GenVertexCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::GenParticleCollection>& genParticlesToken,
        edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
        const reco::Vertex& primaryVertex);

    GenVertexCollection getGenVertexCollection() { return genVertices_; }
    GenVertexCollection getGenVertexSimMatchCollection() { return genVerticesSimMatch_; }
    GenVertexCollection getGenVertexNoNuCollection() { return genVerticesNoNu_; }
    GenVertexCollection getGenVertexNoNuSimMatchCollection() { return genVerticesNoNuSimMatch_; }

  private:

    int genPartID(int pdgId);
    template <class P> bool goodGenParticle(const P* gp, double ptCut);
    template <class P> bool matchGenToSimTrack(const P* gp, const SimTrack& st);
    bool matchGenToSimVertex(const GenVertex& gv);

    double absEtaMax_;
    double genMotherPtMin_;
    double genDaughterPtMin_;
    double drCut_;
    double ptCut_;
    double matchFrac_;

    GenVertexCollection genVertices_;
    GenVertexCollection genVerticesSimMatch_;
    GenVertexCollection genVerticesNoNu_;
    GenVertexCollection genVerticesNoNuSimMatch_;

    reco::GenParticleCollection genParticles_;
    edm::SimTrackContainer simTracks_;
};


#endif
