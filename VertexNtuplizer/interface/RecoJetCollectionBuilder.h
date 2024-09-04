#ifndef VertexNtuples_VertexNtuplizer_RecoJetCollectionBuilder_h
#define VertexNtuples_VertexNtuplizer_RecoJetCollectionBuilder_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "RecoJet.h"


typedef std::vector<RecoJet> RecoJetCollection;


class RecoJetCollectionBuilder {

  public:

    RecoJetCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~RecoJetCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<pat::JetCollection>& jetsToken,
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection>& genJetsFlavourInfoToken);

    RecoJetCollection getRecoJetCollection() { return recoJets_; }
    RecoJetCollection getRecoJetGenMatchCollection() { return recoJetsGenMatch_; }

  private:

    template <class J> bool goodJet(const J& j);

    double jetPtMin_;
    double jetPtMax_;
    double absEtaMax_;
    double drCut_;

    RecoJetCollection recoJets_;
    RecoJetCollection recoJetsGenMatch_;

    pat::JetCollection jets_;
    reco::JetFlavourInfoMatchingCollection genJetsFlavourInfo_;
};


#endif
