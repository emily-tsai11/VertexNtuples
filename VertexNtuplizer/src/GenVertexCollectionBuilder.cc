#include "../interface/GenVertexCollectionBuilder.h"

#include <queue>


GenVertexCollectionBuilder::GenVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  genMotherPtMin_ = iConfig.getUntrackedParameter<double>("genMotherPtMin");
  genDaughterPtMin_ = iConfig.getUntrackedParameter<double>("genDaughterPtMin");
}


// GenVertexCollectionBuilder::~GenVertexCollectionBuilder() {}


unsigned int GenVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::GenParticleCollection>& genParticlesToken,
    edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
    const reco::Vertex& primaryVertex, VertexMatcher* matcher) {

  genParticles_ = iEvent.get(genParticlesToken);
  simTracks_ = iEvent.get(simTracksToken);

  genVertices_.clear();
  genVerticesB_.clear();
  genVerticesD_.clear();
  genVerticesSimMatch_.clear();

  unsigned int nPassingGPs = 0;
  std::vector<bool> simTrackMatches(simTracks_.size(), false);

  for (unsigned int iGP = 0; iGP < genParticles_.size(); iGP++) {
    const reco::Candidate* mother = &genParticles_[iGP];

    if (!goodGenPart(mother, genMotherPtMin_)) continue; // Kinematic cut on mother
    int motherPartID = genPartID(mother->pdgId());
    if (motherPartID < 0) continue; // Mother is not interesting hadron

    bool lastInstance = true; // Check for last instance of interesting hadron
    for (unsigned int iDau = 0; iDau < mother->numberOfDaughters(); iDau++) {
      if (genPartID((mother->daughter(iDau))->pdgId()) == motherPartID) {
        lastInstance = false; // Not last instance of interesting hadron
        break;
      }
    }
    if (!lastInstance) continue;
    nPassingGPs++;

    std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
    for (unsigned int iDau = 0; iDau < mother->numberOfDaughters(); iDau++) {
      const reco::Candidate* dau = mother->daughter(iDau)->clone();
      std::queue<const reco::Candidate*> queue;
      queue.push(dau);
      while(!queue.empty()) {
        if (queue.front()->status() == 1) { // Stable outgoing particle
          if (goodGenPart(queue.front(), genDaughterPtMin_) && queue.front()->charge() != 0) {
            goodDaughters->push_back(queue.front());
          }
        } else {
          for (unsigned int iDau = 0; iDau < queue.front()->numberOfDaughters(); iDau++) {
            queue.push(queue.front()->daughter(iDau));
          }
        }
        queue.pop();
      }
    }
    // Loop through GPs and add all stable particles with this mother (for broken connections)
    for (unsigned int iGP = 0; iGP < genParticles_.size(); iGP++) {
      const reco::Candidate* possibleDau = &genParticles_[iGP];
      if (possibleDau->status() == 1) {
        const reco::Candidate* possibleMother = possibleDau->mother(0);
        if (isAncestor(mother, possibleMother)) {
          goodDaughters->push_back(possibleDau->clone());
        }
      }
    }

    if (goodDaughters->size() >= 2) {
      GenVertex newGV(mother, goodDaughters, primaryVertex, motherPartID);
      genVertices_.push_back(newGV);
      if (motherPartID == B_MESON || motherPartID == B_BARYON) genVerticesB_.push_back(newGV);
      if (motherPartID == C_MESON || motherPartID == C_BARYON) genVerticesD_.push_back(newGV);
      if (matcher->match(newGV, simTracks_, simTrackMatches)) {
        genVerticesSimMatch_.push_back(newGV);
      }
    }
  }

  // Sort collections
  // std::sort(genVertices_.begin(), genVertices_.end(), compare);
  // std::sort(genVerticesB_.begin(), genVerticesB_.end(), compare);
  // std::sort(genVerticesD_.begin(), genVerticesD_.end(), compare);
  // std::sort(genVerticesSimMatch_.begin(), genVerticesSimMatch_.end(), compare);

  return nPassingGPs;
}


bool GenVertexCollectionBuilder::goodGenPart(const reco::Candidate* gp, double ptCut) {

  bool pass = true;
  if (gp->pt() < ptCut) pass = false;
  if (abs(gp->eta()) > absEtaMax_) pass = false;
  return pass;
}


int GenVertexCollectionBuilder::genPartID(int pdgId) {

  int checkPdgId = abs(pdgId);
  if (
    (checkPdgId == 521) ||
    (checkPdgId == 511) ||
    (checkPdgId == 531) ||
    (checkPdgId == 541)) return B_MESON;
  else if (
    (checkPdgId == 5122) ||
    (checkPdgId == 5112) ||
    (checkPdgId == 5212) ||
    (checkPdgId == 5222) ||
    (checkPdgId == 5132) ||
    (checkPdgId == 5232) ||
    (checkPdgId == 5332) ||
    (checkPdgId == 5142) ||
    (checkPdgId == 5242) ||
    (checkPdgId == 5342) ||
    (checkPdgId == 5512) ||
    (checkPdgId == 5532) ||
    (checkPdgId == 5542) ||
    (checkPdgId == 5554)) return B_BARYON;
  else if (
    (checkPdgId == 411) ||
    (checkPdgId == 421) ||
    (checkPdgId == 431)) return C_MESON;
  else if (
    (checkPdgId == 4122) ||
    (checkPdgId == 4222) ||
    (checkPdgId == 4212) ||
    (checkPdgId == 4112) ||
    (checkPdgId == 4232) ||
    (checkPdgId == 4132) ||
    (checkPdgId == 4332) ||
    (checkPdgId == 4412) ||
    (checkPdgId == 4422) ||
    (checkPdgId == 4432) ||
    (checkPdgId == 4444)) return C_BARYON;
  // else if (
  //   (checkPdgId == 3122) ||
  //   (checkPdgId == 3222) ||
  //   (checkPdgId == 3212) ||
  //   (checkPdgId == 3112) ||
  //   (checkPdgId == 3322) ||
  //   (checkPdgId == 3312) ||
  //   (checkPdgId == 3334)) return S_BARYON;
  return -1;
}


bool GenVertexCollectionBuilder::isAncestor(const reco::Candidate* mother,
    const reco::Candidate* possibleMother) {

  if (possibleMother == nullptr) return false;
  if (mother == possibleMother) return true;
  for (unsigned int iMother = 0; iMother < possibleMother->numberOfMothers(); iMother++) {
    if (isAncestor(mother, possibleMother->mother(iMother))) return true;
  }
  return false;
}


// bool GenVertexCollectionBuilder::compare(const GenVertex& gva, const GenVertex& gvb) {

//   return gvb.d3d() < gva.d3d();
// }
