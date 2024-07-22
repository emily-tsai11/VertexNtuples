#ifndef VertexNtuples_VertexNtuplizer_VertexMatcher_h
#define VertexNtuples_VertexNtuplizer_VertexMatcher_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TString.h"
#include "TH1.h"
#include "TH2.h"

#include "VertexCalculator.h"


class GenVertex;
class SecondaryVertex;
class RecoJet;


class VertexMatcher {

  public:

    VertexMatcher(const edm::ParameterSet& iConfig);
    // ~VertexMatcher();

    enum MatchAlgo {
      TRACK,
      // MATRIX
    };

    bool match(const reco::Candidate* daughter, const reco::Track& gt);
    bool match(const reco::Candidate* daughter, const reco::PFCandidate& pfc);

    bool match(const GenVertex& gv, const edm::SimTrackContainer& simTracks, std::vector<bool>& matches);
    bool match(GenVertex& gv, SecondaryVertex& sv, MatchAlgo algo);
    bool match(const SecondaryVertex& sv, const RecoJet& rj);

    void fill(std::map<TString, TH1F*>& histos1, std::map<TString, TH2F*>& histos2,
        TString gvPrefix, TString svPrefix, const GenVertex& gv, const SecondaryVertex& sv);

  private:

    bool match(const reco::Candidate* daughter, const SimTrack& st);
    bool match(const float pt1, const float eta1, const float phi1,
        const float pt2, const float eta2, const float phi2);
    float ptResNorm(const float pt1, const float pt2);

    bool vtxTrackMatch(GenVertex& gv, SecondaryVertex& sv);
    // Not yet implemented, probably have to change input and return values
    // bool vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv);

    double trkMatchPtCut_;
    double trkMatchDrCut_;
    double vtxMatchFrac_;
    double jetRadius_;
};


#endif
