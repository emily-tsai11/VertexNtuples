#include "GenVertex.h"

GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
  edm::SimTrackContainer simTracks, float matchFrac, float drCut, float ptCut) :
  mother_(mother), daughters_(daughters), isSimMatched_(false) {

  // Mother and daughter pt should already be pre-selected, and there should be >=2 tracks

  daughtersNoNu_ = new std::vector<const reco::Candidate*>;
  for (const reco::Candidate* dau : *(daughters_)) {
    if (abs(dau->pdgId()) != 12 && abs(dau->pdgId()) != 14 && abs(dau->pdgId()) != 16)
    daughtersNoNu_->push_back(dau->clone());
  }

  float nmatch = 0;
  for (const reco::Candidate* dau : *(daughters_)) {
  for (const SimTrack& st : simTracks) {
    if (reco::deltaR(dau->eta(), dau->phi(), st.momentum().Eta(), st.momentum().Phi()) > drCut) continue;
    if (abs(dau->pt() - st.momentum().Pt()) / (dau->pt() + st.momentum().Pt()) > ptCut) continue;
      nmatch++;
      break;
    }
  }
  isSimMatched_ = (nmatch/(float)daughters_->size()) >= matchFrac;

  float nmatchNoNu = 0;
  for (const reco::Candidate* dau : *(daughtersNoNu_)) {
    for (const SimTrack& st : simTracks) {
    if (reco::deltaR(dau->eta(), dau->phi(), st.momentum().Eta(), st.momentum().Phi()) > drCut) continue;
    if (abs(dau->pt() - st.momentum().Pt()) / (dau->pt() + st.momentum().Pt()) > ptCut) continue;
      nmatchNoNu++;
      break;
    }
  }
  isSimMatchedNoNu_ = (nmatchNoNu/(float)daughtersNoNu_->size()) >= matchFrac;
}

GenVertex::~GenVertex() {
  delete daughters_;
  delete daughtersNoNu_;
}
