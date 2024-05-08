#include "../interface/VertexMatcher.h"


VertexMatcher::VertexMatcher(const edm::ParameterSet& iConfig) {

  trkMatchDrCut_ = iConfig.getUntrackedParameter<double>("trkMatchDrCut");
  trkMatchPtCut_ = iConfig.getUntrackedParameter<double>("trkMatchPtCut");
  vtxMatchFrac_ = iConfig.getUntrackedParameter<double>("vtxMatchFrac");
  jetRadius_ = iConfig.getUntrackedParameter<double>("jetRadius");
}


// VertexMatcher::~VertexMatcher() {}


bool VertexMatcher::match(GenVertex& gv, SecondaryVertex& sv, MatchAlgo algo) {

  if (algo == TRACK) {
    return vtxTrackMatch(gv, sv);
  }
  else if (algo == MATRIX) {
    return vtxMatrixMatch(gv, sv);
  }
  std::cout << "WARNING in VertexMatcher::match(gv, sv): No matching algorithm specified. Returning false." << std::endl;
  return false;
}


bool VertexMatcher::match(SecondaryVertex& sv, RecoJet& rj) {

  return reco::deltaR(sv.eta(), sv.phi(), rj.eta(), rj.phi()) < jetRadius_;
}


bool VertexMatcher::vtxTrackMatch(GenVertex& gv, SecondaryVertex& sv) {

  float nmatch = 0;
  for (const reco::Candidate* dau : *(gv.daughters())) {
    for (unsigned int iTrk = 0; iTrk < sv.nTracks(); iTrk++) {
      bool match = true;
      if (reco::deltaR(sv.trkEta()->at(iTrk), sv.trkPhi()->at(iTrk), dau->eta(), dau->phi()) > trkMatchDrCut_) match = false;
      if (abs(dau->pt() - sv.trkPt()->at(iTrk)) / (dau->pt() + sv.trkPt()->at(iTrk)) > trkMatchPtCut_) match = false;
      if (match) {
        nmatch++;
        break;
      }
    }
  }
  return (nmatch/(float)gv.nDaughters()) >= vtxMatchFrac_;
}


bool VertexMatcher::vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv) {

  // Not yet implemented
  return false;
}
