#include "../interface/GenJetCollectionBuilder.h"


GenJetCollectionBuilder::GenJetCollectionBuilder(const edm::ParameterSet& iConfig) {

  jetPtMin_ = iConfig.getUntrackedParameter<double>("jetPtMin");
  jetPtMax_ = iConfig.getUntrackedParameter<double>("jetPtMax");
  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  drCut_ = iConfig.getUntrackedParameter<double>("jetMatchDrCut");
}


// GenJetCollectionBuilder::~GenJetCollectionBuilder() {}


void GenJetCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::GenJetCollection>& genJetsToken,
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection>& genJetsFlavourInfoToken,
    edm::EDGetTokenT<pat::JetCollection>& jetsToken) {

  cmsGenJets_ = iEvent.get(genJetsToken);
  genJetsFlavourInfo_ = iEvent.get(genJetsFlavourInfoToken);
  jets_ = iEvent.get(jetsToken);

  genJets_.clear();
  genJetsRecoMatch_.clear();

  for (const reco::GenJet& genJet : cmsGenJets_) {
    if (!goodJet(genJet)) continue;
    GenJet newGJ(genJet, genJetsFlavourInfo_, drCut_);
    genJets_.push_back(newGJ);

    bool matchRecoJet = false;
    for (const pat::Jet& jet : jets_) {
      if (!goodJet(jet)) continue;
      const reco::GenJet* matchedGenJet = jet.genJet();
      if (!matchedGenJet) continue;
      if (matchedGenJet->pt() != genJet.pt() ||
          matchedGenJet->eta() != genJet.eta() ||
          matchedGenJet->phi() != genJet.phi()) continue; // No matching reco jet found
      matchRecoJet = true;
      break;
    }
    if (matchRecoJet) {
      genJetsRecoMatch_.push_back(newGJ);
    }
  }
}


template <class J>
bool GenJetCollectionBuilder::goodJet(const J& j) {

  bool pass = true;
  if (j.pt() < jetPtMin_) pass = false;
  if (j.pt() > jetPtMax_) pass = false;
  if (abs(j.eta()) > absEtaMax_) pass = false;
  return pass;
}
