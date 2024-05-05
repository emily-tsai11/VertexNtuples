#include "../interface/GenJetCollectionBuilder.h"


GenJetCollectionBuilder::GenJetCollectionBuilder(const edm::ParameterSet& iConfig) {

  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  jetPtMin_ = iConfig.getUntrackedParameter<double>("jetPtMin");
  jetPtMax_ = iConfig.getUntrackedParameter<double>("jetPtMax");
  drCut_ = iConfig.getUntrackedParameter<double>("jetMatchDrCut");
}


// GenJetCollectionBuilder::~GenJetCollectionBuilder() {}


void GenJetCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::GenJetCollection> genJetsToken,
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken) {

  cmsGenJets_ = iEvent.get(genJetsToken);
  genJetsFlavourInfo_ = iEvent.get(genJetsFlavourInfoToken);

  genJets_.clear();
  genJetsRecoMatch_.clear();

  //
}


bool GenJetCollectionBuilder::goodGenJet(const reco::GenJet& gj) {

  bool pass = true;
  //
  return pass;
}


unsigned int GenJetCollectionBuilder::getHadronFlavour(const reco::GenJet& gj) {

  return 0;
}
