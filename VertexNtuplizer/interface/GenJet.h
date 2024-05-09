#ifndef VertexNtuples_VertexNtuplizer_GenJet_h
#define VertexNtuples_VertexNtuplizer_GenJet_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TH1.h"


class GenJet {

  public:

    GenJet(const reco::GenJet& genJet, reco::JetFlavourInfoMatchingCollection& genJetsFlavourInfo,
        double jetMatchDrCut);
    // ~GenJet();

    void fill(std::map<TString, TH1F*>& histos, TString prefix);

    const double pt() const { return pt_; }
    const double eta() const { return eta_; }
    const double phi() const { return phi_; }
    const double hadFlav() const { return hadFlav_; }

  private:

    double pt_;
    double eta_;
    double phi_;
    double hadFlav_;
};


#endif
