#include "../interface/RecoJet.h"


RecoJet::RecoJet(const pat::Jet& jet, reco::JetFlavourInfoMatchingCollection& genJetsFlavourInfo,
    double jetMatchDrCut) {

  int hadFlav = -1000;
  const reco::GenJet* genJet = jet.genJet();
  if (genJet) {
    for (const reco::JetFlavourInfoMatching& genJetFlavInfo : genJetsFlavourInfo) {
      if (reco::deltaR(genJet->p4(), genJetFlavInfo.first->p4()) < jetMatchDrCut) {
        hadFlav = genJetFlavInfo.second.getHadronFlavour();
        break;
      }
    }
  }

  pt_ = jet.pt();
  eta_ = jet.eta();
  phi_ = jet.phi();
  hadFlav_ = hadFlav;
}


// RecoJet::~RecoJet() {}
