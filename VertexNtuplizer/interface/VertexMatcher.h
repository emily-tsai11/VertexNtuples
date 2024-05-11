#ifndef VertexNtuples_VertexNtuplizer_VertexMatcher_h
#define VertexNtuples_VertexNtuplizer_VertexMatcher_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1.h"

#include "VertexCalculatorExtended.h"
#include "GenVertex.h"
#include "SecondaryVertex.h"
#include "RecoJet.h"


enum MatchAlgo {
  TRACK,
  MATRIX
};


class VertexMatcher {

  public:

    VertexMatcher(const edm::ParameterSet& iConfig);
    // ~VertexMatcher();

    bool match(GenVertex& gv, SecondaryVertex& sv, MatchAlgo algo);
    bool match(SecondaryVertex& sv, RecoJet& rj);

    void fill(std::map<TString, TH1F*>& histos1, std::map<TString, TH2F*>& histos2,
        TString gvPrefix, TString svPrefix, GenVertex& gv, SecondaryVertex& sv);

  private:

    bool vtxTrackMatch(GenVertex& gv, SecondaryVertex& sv);
    // Not yet implemented, probably have to change input and return values
    bool vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv);

    double trkMatchDrCut_;
    double trkMatchPtCut_;
    double vtxMatchFrac_;
    double jetRadius_;
};


#endif
