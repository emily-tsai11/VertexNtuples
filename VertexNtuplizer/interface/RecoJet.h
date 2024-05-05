#ifndef RECO_JET
#define RECO_JET


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/PatCandidates/interface/Jet.h"


class RecoJet {

  public:

    RecoJet(const pat::Jet& jet);
    // ~RecoJet();

  private:

    double pt_;
    double eta_;
    double phi_;
    double radius_;
    double hadFlav_;
};


#endif
