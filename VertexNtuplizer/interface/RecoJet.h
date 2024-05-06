#ifndef RECO_JET
#define RECO_JET


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"


class RecoJet {

  public:

    RecoJet(const pat::Jet& jet, reco::JetFlavourInfoMatchingCollection& genJetsFlavourInfo,
        double jetMatchDrCut);
    // ~RecoJet();

  private:

    double pt_;
    double eta_;
    double phi_;
    double hadFlav_;
};


#endif
