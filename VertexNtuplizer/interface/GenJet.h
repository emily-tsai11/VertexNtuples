#ifndef GEN_JET
#define GEN_JET


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"


class GenJet {

  public:

    GenJet(const reco::GenJet& genJet, reco::JetFlavourInfoMatchingCollection& genJetsFlavourInfo,
        double jetMatchDrCut);
    // ~GenJet();

    const float pt() const { return pt_; }
    const float eta() const { return eta_; }
    const float phi() const { return phi_; }
    const float hadFlav() const { return hadFlav_; }

  private:

    double pt_;
    double eta_;
    double phi_;
    double hadFlav_;
};


#endif
