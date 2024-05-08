#ifndef VERTEX_MATCHER
#define VERTEX_MATCHER


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

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

    bool match(GenVertex& gv, SecondaryVertex& sv, MatchAlgo algo=TRACK);
    bool match(SecondaryVertex& sv, RecoJet& rj);

    bool vtxTrackMatch(GenVertex& gv, SecondaryVertex& sv);
    bool vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv);

  private:

    double trkMatchDrCut_;
    double trkMatchPtCut_;
    double vtxMatchFrac_;
    double jetRadius_;
};


#endif
