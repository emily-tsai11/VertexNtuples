#ifndef VertexNtuples_VertexNtuplizer_GenJetCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_GenJetCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "GenJet.h"


typedef std::vector<GenJet> GenJetCollection;


class GenJetCollectionBuilder {

  public:

    GenJetCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~GenJetCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::GenJetCollection>& genJetsToken,
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection>& genJetsFlavourInfoToken,
        edm::EDGetTokenT<pat::JetCollection>& jetsToken);

    GenJetCollection getGenJetCollection() { return genJets_; }
    GenJetCollection getGenJetRecoMatchCollection() { return genJetsRecoMatch_; }

  private:

    template <class J> bool goodJet(const J& j);

    double jetPtMin_;
    double jetPtMax_;
    double absEtaMax_;
    double drCut_;

    GenJetCollection genJets_;
    GenJetCollection genJetsRecoMatch_;

    reco::GenJetCollection cmsGenJets_;
    reco::JetFlavourInfoMatchingCollection genJetsFlavourInfo_;
    pat::JetCollection jets_;
};


#endif
