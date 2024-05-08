#include "../interface/GenJet.h"


GenJet::GenJet(const reco::GenJet& genJet, reco::JetFlavourInfoMatchingCollection& genJetsFlavourInfo,
    double jetMatchDrCut) {

  int hadFlav = -1000;
  for (const reco::JetFlavourInfoMatching& genJetFlavInfo : genJetsFlavourInfo) {
    if (reco::deltaR(genJet.p4(), genJetFlavInfo.first->p4()) < jetMatchDrCut) {
      hadFlav = genJetFlavInfo.second.getHadronFlavour();
      break;
    }
  }

  pt_ = genJet.pt();
  eta_ = genJet.eta();
  phi_ = genJet.phi();
  hadFlav_ = hadFlav;
}


// GenJet::~GenJet() {}


void GenJet::fill(std::map<TString, TH1F*>& histos, TString prefix) {

  histos[prefix + "_pt"]->Fill(pt());
  histos[prefix + "_eta"]->Fill(eta());
  histos[prefix + "_phi"]->Fill(phi());
  histos[prefix + "_hadFlav"]->Fill(hadFlav());
}
