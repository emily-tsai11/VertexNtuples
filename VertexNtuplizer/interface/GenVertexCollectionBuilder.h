#ifndef VertexNtuples_VertexNtuplizer_GenVertexCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_GenVertexCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "GenVertex.h"
#include "VertexMatcher.h"


typedef std::vector<GenVertex> GenVertexCollection;


enum HadronType {
  B_MESON,
  B_BARYON,
  C_MESON,
  C_BARYON,
  // S_BARYON
};


class GenVertexCollectionBuilder {

  public:

    GenVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~GenVertexCollectionBuilder();

    unsigned int build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::GenParticleCollection>& genParticlesToken,
        edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
        const reco::Vertex& primaryVertex, VertexMatcher* matcher);

    GenVertexCollection getGenVertices() { return genVertices_; }
    GenVertexCollection getGenVerticesSimMatch() { return genVerticesSimMatch_; }
    GenVertexCollection getGenVerticesB() { return genVerticesB_; }
    GenVertexCollection getGenVerticesD() { return genVerticesD_; }

  private:

    bool goodGenPart(const reco::Candidate* gp, double ptCut);
    int genPartID(int pdgId);
    bool isAncestor(const reco::Candidate* mother, const reco::Candidate* possibleMother);
    // static bool compare(const GenVertex& gva, const GenVertex& gvb);

    double genMotherPtMin_;
    double genDaughterPtMin_;
    double absEtaMax_;

    GenVertexCollection genVertices_;
    GenVertexCollection genVerticesSimMatch_;
    GenVertexCollection genVerticesB_;
    GenVertexCollection genVerticesD_;

    reco::GenParticleCollection genParticles_;
    edm::SimTrackContainer simTracks_;
};


#endif
