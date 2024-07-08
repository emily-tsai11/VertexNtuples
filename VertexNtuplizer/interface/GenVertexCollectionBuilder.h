#ifndef VertexNtuples_VertexNtuplizer_GenVertexCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_GenVertexCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
// #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "GenVertex.h"


typedef std::vector<GenVertex> GenVertexCollection;


class GenVertexCollectionBuilder {

  public:

    GenVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~GenVertexCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::GenParticleCollection>& prunedGenParticlesToken,
        edm::EDGetTokenT<pat::PackedGenParticleCollection>& packedGenParticlesToken,
        edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
        // edm::EDGetTokenT<TrackingParticleCollection>& trackingParticlesToken,
        edm::EDGetTokenT<TrackingVertexCollection>& trackingVerticesToken,
        const reco::Vertex& primaryVertex);

    GenVertexCollection getGenVertexCollection() { return genVertices_; }
    GenVertexCollection getGenVertexSimMatchCollection() { return genVerticesSimMatch_; }
    GenVertexCollection getGenVertexNoNuCollection() { return genVerticesNoNu_; }
    GenVertexCollection getGenVertexNoNuSimMatchCollection() { return genVerticesNoNuSimMatch_; }
    GenVertexCollection getGenVertexFromPackedGen() { return genVerticesFromPackedGen_; }
    GenVertexCollection getGenVertexFromPackedGenNoNu() { return genVerticesFromPackedGenNoNu_; }
    GenVertexCollection getGenVertexFromPackedGenNoNuSimMatch() { return genVerticesFromPackedGenNoNuSimMatch_; }
    GenVertexCollection getGenVertexFromTV() { return genVerticesFromTV_; }
    GenVertexCollection getGenVertexFromTVNoNu() { return genVerticesFromTVNoNu_; }

    typedef edm::Ref<edm::HepMCProduct, HepMC::GenVertex> GenVertexRef;

  private:

    typedef HepMC::GenVertex::particles_in_const_iterator it_in;
    typedef HepMC::GenVertex::particles_out_const_iterator it_out;

    template <class P> bool goodGenParticle(const P* gp, double ptCut);
    bool goodHepMCGenParticle(const HepMC::GenParticle* gp, double ptCut);
    template <class P> bool matchGenToSimTrack(const P* gp, const SimTrack& st);
    bool matchGenToSimVertex(const GenVertex& gv);
    int genPartID(int pdgId);
    // static bool compare(const GenVertex& gva, const GenVertex& gvb);

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
    GenVertexCollection genVerticesFromPackedGen_;
    GenVertexCollection genVerticesFromPackedGenNoNu_;
    GenVertexCollection genVerticesFromPackedGenNoNuSimMatch_;
    GenVertexCollection genVerticesFromTV_;
    GenVertexCollection genVerticesFromTVNoNu_;

    reco::GenParticleCollection prunedGenParticles_;
    pat::PackedGenParticleCollection packedGenParticles_;
    edm::SimTrackContainer simTracks_;
    // TrackingParticleCollection trackingParticles_;
    TrackingVertexCollection trackingVertices_;
};


#endif
