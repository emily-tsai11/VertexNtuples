#ifndef VertexNtuples_VertexNtuplizer_GenVertexCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_GenVertexCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "GenVertex.h"


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

    std::vector<unsigned int> build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::GenParticleCollection>& prunedGenParticlesToken,
        edm::EDGetTokenT<reco::GenParticleCollection>& mergedGenParticlesToken,
        edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
        const reco::Vertex& primaryVertex);

    GenVertexCollection getGenVertexFromPruned() { return genVerticesFromPruned_; }
    GenVertexCollection getGenVertexFromPrunedB() { return genVerticesFromPrunedB_; }
    GenVertexCollection getGenVertexFromPrunedD() { return genVerticesFromPrunedD_; }
    GenVertexCollection getGenVertexFromPrunedSimMatch() { return genVerticesFromPrunedSimMatch_; }
    GenVertexCollection getGenVertexFromMerged() { return genVerticesFromMerged_; }
    GenVertexCollection getGenVertexFromMergedB() { return genVerticesFromMergedB_; }
    GenVertexCollection getGenVertexFromMergedD() { return genVerticesFromMergedD_; }
    GenVertexCollection getGenVertexFromMergedSimMatch() { return genVerticesFromMergedSimMatch_; }

  private:

    bool goodGenParticle(const reco::Candidate* gp, double ptCut);
    bool isAncestor(const reco::Candidate* mother, const reco::Candidate* possibleMother);
    bool matchGenToSimTrack(const reco::Candidate* gp, const SimTrack& st);
    bool matchGenToSimVertex(const GenVertex& gv);
    int genPartID(int pdgId);
    // static bool compare(const GenVertex& gva, const GenVertex& gvb);

    double absEtaMax_;
    double genMotherPtMin_;
    double genDaughterPtMin_;
    double drCut_;
    double ptCut_;
    double matchFrac_;

    GenVertexCollection genVerticesFromPruned_;
    GenVertexCollection genVerticesFromPrunedB_;
    GenVertexCollection genVerticesFromPrunedD_;
    GenVertexCollection genVerticesFromPrunedSimMatch_;
    GenVertexCollection genVerticesFromMerged_;
    GenVertexCollection genVerticesFromMergedB_;
    GenVertexCollection genVerticesFromMergedD_;
    GenVertexCollection genVerticesFromMergedSimMatch_;

    reco::GenParticleCollection prunedGenParticles_;
    reco::GenParticleCollection mergedGenParticles_;
    edm::SimTrackContainer simTracks_;
};


#endif
