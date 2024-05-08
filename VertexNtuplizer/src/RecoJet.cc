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


void RecoJet::fill(std::map<TString, TH1F*>& histos, TString prefix) {

  histos[prefix + "_pt"]->Fill(pt());
  histos[prefix + "_eta"]->Fill(eta());
  histos[prefix + "_phi"]->Fill(phi());
  histos[prefix + "_hadFlav"]->Fill(hadFlav());
}
