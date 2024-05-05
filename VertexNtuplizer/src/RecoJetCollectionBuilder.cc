#include "../interface/RecoJetCollectionBuilder.h"


RecoJetCollectionBuilder::RecoJetCollectionBuilder(const edm::ParameterSet& iConfig) {

  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  jetPtMin_ = iConfig.getUntrackedParameter<double>("jetPtMin");
  jetPtMax_ = iConfig.getUntrackedParameter<double>("jetPtMax");
  drCut_ = iConfig.getUntrackedParameter<double>("jetMatchDrCut");
}


// RecoJetCollectionBuilder::~RecoJetCollectionBuilder() {}


void RecoJetCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<pat::JetCollection> jetsToken,
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken) {

  jets_ = iEvent.get(jetsToken);
  genJetsFlavourInfo_ = iEvent.get(genJetsFlavourInfoToken);

  recoJets_.clear();
  recoJetsGenMatch_.clear();

  //
}


bool RecoJetCollectionBuilder::goodRecoJet(const pat::Jet& j) {

  bool pass = true;
  //
  return pass;
}


unsigned int RecoJetCollectionBuilder::getHadronFlavour(const reco::GenJet& gj) {

  return 0;
}
