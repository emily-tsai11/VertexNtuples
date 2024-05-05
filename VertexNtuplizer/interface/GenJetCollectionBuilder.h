#ifndef GEN_JET_COLLECTION_BUILDER
#define GEN_JET_COLLECTION_BUILDER


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "GenJet.h"


typedef std::vector<GenJet> GenJetCollection;


class GenJetCollectionBuilder {

  public:

    GenJetCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~GenJetCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::GenJetCollection> genJetsToken,
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken);

    GenJetCollection getGenJetCollection() { return genJets_; }
    GenJetCollection getGenJetRecoMatchCollection() { return genJetsRecoMatch_; }

  private:

    bool goodGenJet(const reco::GenJet& gj);
    unsigned int getHadronFlavour(const reco::GenJet& gj);

    double absEtaMax_;
    double jetPtMin_;
    double jetPtMax_;
    double drCut_;

    GenJetCollection genJets_;
    GenJetCollection genJetsRecoMatch_;

    reco::GenJetCollection cmsGenJets_;
    reco::JetFlavourInfoMatchingCollection genJetsFlavourInfo_;
};


#endif
