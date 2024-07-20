#include "../interface/GenVertexCollectionBuilder.h"

#include <queue>


GenVertexCollectionBuilder::GenVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  genMotherPtMin_ = iConfig.getUntrackedParameter<double>("genMotherPtMin");
  genDaughterPtMin_ = iConfig.getUntrackedParameter<double>("genDaughterPtMin");
  drCut_ = iConfig.getUntrackedParameter<double>("simTrkMatchDrCut");
  ptCut_ = iConfig.getUntrackedParameter<double>("simTrkMatchPtCut");
  matchFrac_ = iConfig.getUntrackedParameter<double>("vtxMatchFrac");
}


// GenVertexCollectionBuilder::~GenVertexCollectionBuilder() {}


std::vector<unsigned int> GenVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::GenParticleCollection>& prunedGenParticlesToken,
    edm::EDGetTokenT<reco::GenParticleCollection>& mergedGenParticlesToken,
    edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
    const reco::Vertex& primaryVertex) {

  prunedGenParticles_ = iEvent.get(prunedGenParticlesToken);
  mergedGenParticles_ = iEvent.get(mergedGenParticlesToken);
  simTracks_ = iEvent.get(simTracksToken);

  genVerticesFromPruned_.clear();
  genVerticesFromPrunedB_.clear();
  genVerticesFromPrunedD_.clear();
  genVerticesFromPrunedSimMatch_.clear();
  genVerticesFromMerged_.clear();
  genVerticesFromMergedB_.clear();
  genVerticesFromMergedD_.clear();
  genVerticesFromMergedSimMatch_.clear();

  std::vector<unsigned int> nPassingGPs(2, 0);
  // 0: number of passing pruned GenParticles
  // 1: number of passing merged GenParticles

  // Construct GenVertexCollections from pruned GenParticles
  for (unsigned int iGP = 0; iGP < prunedGenParticles_.size(); iGP++) {
    const reco::Candidate* mother = &prunedGenParticles_[iGP];

    if (!goodGenParticle(mother, genMotherPtMin_)) continue; // Kinematic cut on mother
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
    nPassingGPs.at(0)++;

    std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
    for (unsigned int iDau = 0; iDau < mother->numberOfDaughters(); iDau++) {
      const reco::Candidate* dau = mother->daughter(iDau)->clone();
      std::queue<const reco::Candidate*> queue;
      queue.push(dau);
      while(!queue.empty()) {
        if (queue.front()->status() == 1) { // Stable outgoing particle
          if (goodGenParticle(queue.front(), genDaughterPtMin_) && queue.front()->charge() != 0) {
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

    if (goodDaughters->size() >= 2) {
      GenVertex newGV(mother, goodDaughters, primaryVertex, motherPartID);
      genVerticesFromPruned_.push_back(newGV);
      if (motherPartID == B_MESON || motherPartID == B_BARYON) genVerticesFromPrunedB_.push_back(newGV);
      if (motherPartID == C_MESON || motherPartID == C_BARYON) genVerticesFromPrunedD_.push_back(newGV);
      if (matchGenToSimVertex(newGV)) {
        genVerticesFromPrunedSimMatch_.push_back(newGV);
      }
    }
  }

  // Construct GenVertexCollections from merged GenParticles
  for (unsigned int iGP = 0; iGP < mergedGenParticles_.size(); iGP++) {
    const reco::Candidate* mother = &mergedGenParticles_[iGP];

    if (!goodGenParticle(mother, genMotherPtMin_)) continue; // Kinematic cut on mother
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
    nPassingGPs.at(1)++;

    std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
    for (unsigned int iDau = 0; iDau < mother->numberOfDaughters(); iDau++) {
      const reco::Candidate* dau = mother->daughter(iDau)->clone();
      std::queue<const reco::Candidate*> queue;
      queue.push(dau);
      while(!queue.empty()) {
        if (queue.front()->status() == 1) { // Stable outgoing particle
          if (goodGenParticle(queue.front(), genDaughterPtMin_) && queue.front()->charge() != 0) {
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
    for (unsigned int iGP = 0; iGP < mergedGenParticles_.size(); iGP++) {
      const reco::Candidate* possibleDau = &mergedGenParticles_[iGP];
      if (possibleDau->status() == 1) {
        const reco::Candidate* possibleMother = possibleDau->mother(0);
        if (possibleMother != nullptr && isAncestor(mother, possibleMother)) {
          goodDaughters->push_back(possibleDau->clone());
        }
      }
    }

    if (goodDaughters->size() >= 2) {
      GenVertex newGV(mother, goodDaughters, primaryVertex, motherPartID);
      genVerticesFromMerged_.push_back(newGV);
      if (motherPartID == B_MESON || motherPartID == B_BARYON) genVerticesFromMergedB_.push_back(newGV);
      if (motherPartID == C_MESON || motherPartID == C_BARYON) genVerticesFromMergedD_.push_back(newGV);
      if (matchGenToSimVertex(newGV)) {
        genVerticesFromMergedSimMatch_.push_back(newGV);
      }
    }
  }

  // Sort collections
  // std::sort(genVerticesFromPruned_.begin(), genVerticesFromPruned_.end(), compare);
  // std::sort(genVerticesFromPrunedB_.begin(), genVerticesFromPrunedB_.end(), compare);
  // std::sort(genVerticesFromPrunedD_.begin(), genVerticesFromPrunedD_.end(), compare);
  // std::sort(genVerticesFromPrunedSimMatch_.begin(), genVerticesFromPrunedSimMatch_.end(), compare);
  // std::sort(genVerticesFromMerged_.begin(), genVerticesFromMerged_.end(), compare);
  // std::sort(genVerticesFromMergedB_.begin(), genVerticesFromMergedB_.end(), compare);
  // std::sort(genVerticesFromMergedD_.begin(), genVerticesFromMergedD_.end(), compare);
  // std::sort(genVerticesFromMergedSimMatch_.begin(), genVerticesFromMergedSimMatch_.end(), compare);

  return nPassingGPs;
}


bool GenVertexCollectionBuilder::goodGenParticle(const reco::Candidate* gp, double ptCut) {

  bool pass = true;
  if (gp->pt() < ptCut) pass = false;
  if (abs(gp->eta()) > absEtaMax_) pass = false;
  return pass;
}


bool GenVertexCollectionBuilder::isAncestor(const reco::Candidate* mother,
    const reco::Candidate* possibleMother) {

  if (mother == possibleMother) return true;
  for (unsigned int iMother = 0; iMother < possibleMother->numberOfMothers(); iMother++) {
    if (isAncestor(mother, possibleMother->mother(iMother))) return true;
  }
  return false;
}


bool GenVertexCollectionBuilder::matchGenToSimTrack(const reco::Candidate* gp, const SimTrack& st) {

  bool match = true;
  if (reco::deltaR(gp->eta(), gp->phi(), st.momentum().Eta(), st.momentum().Phi()) > drCut_) match = false;
  if (abs(gp->pt() - st.momentum().Pt()) / (gp->pt() + st.momentum().Pt()) > ptCut_) match = false;
  return match;
}


bool GenVertexCollectionBuilder::matchGenToSimVertex(const GenVertex& gv) {

  float nmatch = 0;
  std::vector<bool> STmatched(simTracks_.size(), false);
  for (const reco::Candidate* dau : *(gv.daughters())) {
    for (unsigned int iST = 0; iST < simTracks_.size(); iST++) {
      const SimTrack& st = simTracks_.at(iST);
      if (matchGenToSimTrack(dau, st) && !STmatched.at(iST)) {
        STmatched.at(iST) = true;
        nmatch++;
        break;
      }
    }
  }
  return (nmatch / (float) gv.nDaughters()) >= matchFrac_;
}


int GenVertexCollectionBuilder::genPartID(int pdgId) {

  int checkPdgId = abs(pdgId);
  if ((checkPdgId == 521) ||
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


// bool GenVertexCollectionBuilder::compare(const GenVertex& gva, const GenVertex& gvb) {

//   return gvb.d3d() < gva.d3d();
// }
