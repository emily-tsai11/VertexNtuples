#include "../interface/RecoJetCollectionBuilder.h"


RecoJetCollectionBuilder::RecoJetCollectionBuilder(const edm::ParameterSet& iConfig) {

  jetPtMin_ = iConfig.getUntrackedParameter<double>("jetPtMin");
  jetPtMax_ = iConfig.getUntrackedParameter<double>("jetPtMax");
  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  drCut_ = iConfig.getUntrackedParameter<double>("jetMatchDrCut");
}


// RecoJetCollectionBuilder::~RecoJetCollectionBuilder() {}


void RecoJetCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<pat::JetCollection>& jetsToken,
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection>& genJetsFlavourInfoToken) {

  jets_ = iEvent.get(jetsToken);
  genJetsFlavourInfo_ = iEvent.get(genJetsFlavourInfoToken);

  recoJets_.clear();
  recoJetsGenMatch_.clear();

  for (const pat::Jet& jet : jets_) {
    if (!goodJet(jet)) continue;
    RecoJet newRJ(jet, genJetsFlavourInfo_, drCut_);
    recoJets_.push_back(newRJ);
    if (jet.genJet() && goodJet(*(jet.genJet()))) {
      recoJetsGenMatch_.push_back(newRJ);
    }
  }
}


template <class J>
bool RecoJetCollectionBuilder::goodJet(const J& j) {

  bool pass = true;
  if (j.pt() < jetPtMin_) pass = false;
  if (j.pt() > jetPtMax_) pass = false;
  if (abs(j.eta()) > absEtaMax_) pass = false;
  return pass;
}
