#ifndef GEN_JET
#define GEN_JET


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/JetReco/interface/GenJet.h"


class GenJet {

  public:

    GenJet(const reco::GenJet& genJet);
    // ~GenJet();

  private:

    double pt_;
    double eta_;
    double phi_;
    double radius_;
    double hadFlav_;
};


#endif
