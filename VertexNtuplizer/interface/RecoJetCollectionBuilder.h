#ifndef RECO_JET_COLLECTION_BUILDER
#define RECO_JET_COLLECTION_BUILDER


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoJet.h"


typedef std::vector<RecoJet> RecoJetCollection;


class RecoJetCollectionBuilder {

  public:

    RecoJetCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~RecoJetCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<pat::JetCollection> jetsToken,
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken);

    RecoJetCollection recoJetCollection() { return recoJets_; }
    RecoJetCollection recoJetGenMatchCollection() { return recoJetsGenMatch_; }

  private:

    bool goodRecoJet(const pat::Jet& gj);
    unsigned int getHadronFlavour(const reco::GenJet& gj);

    double absEtaMax_;
    double jetPtMin_;
    double jetPtMax_;
    double drCut_;

    RecoJetCollection recoJets_;
    RecoJetCollection recoJetsGenMatch_;

    pat::JetCollection jets_;
    reco::JetFlavourInfoMatchingCollection genJetsFlavourInfo_;
};


#endif
