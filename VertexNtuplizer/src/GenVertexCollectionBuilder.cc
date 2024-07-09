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


unsigned int GenVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::GenParticleCollection>& prunedGenParticlesToken,
    // edm::EDGetTokenT<pat::PackedGenParticleCollection>& packedGenParticlesToken,
    edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
    // edm::EDGetTokenT<TrackingParticleCollection>& trackingParticlesToken,
    edm::EDGetTokenT<TrackingVertexCollection>& trackingVerticesToken,
    const reco::Vertex& primaryVertex) {

  // std::cout << "--------------------- NEW EVENT ---------------------" << std::endl;

  prunedGenParticles_ = iEvent.get(prunedGenParticlesToken);
  // packedGenParticles_ = iEvent.get(packedGenParticlesToken);
  simTracks_ = iEvent.get(simTracksToken);
  // trackingParticles_ = iEvent.get(trackingParticlesToken);
  trackingVertices_ = iEvent.get(trackingVerticesToken);

  genVertices_.clear();
  genVerticesSimMatch_.clear();
  genVerticesNoNu_.clear();
  genVerticesNoNuSimMatch_.clear();
  genVerticesFromPrunedGenNoNu_.clear();
  genVerticesB_.clear();
  genVerticesD_.clear();
  genVerticesFromPrunedGenNoNuSimMatch_.clear();
  genVerticesFromTV_.clear();
  genVerticesFromTVNoNu_.clear();

  // Construct GenVertexCollections from pruned GenParticles
  for (const reco::GenParticle& gp : prunedGenParticles_) {
    const reco::GenParticle* mother = gp.clone();

    if (mother->numberOfDaughters() < 2) continue; // Not a vertex
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

    std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
    std::vector<const reco::Candidate*>* goodDaughtersNoNu = new std::vector<const reco::Candidate*>;
    for (unsigned int iDau = 0; iDau < mother->numberOfDaughters(); iDau++) {
      const reco::Candidate* dau = mother->daughter(iDau);
      if (goodGenParticle(dau, genDaughterPtMin_)) {
        goodDaughters->push_back(dau->clone());
        if (abs(dau->pdgId()) != 12 && abs(dau->pdgId()) != 14 && abs(dau->pdgId()) != 16) {
          goodDaughtersNoNu->push_back(dau->clone());
        }
      }
    }

    if (goodDaughters->size() >= 2) {
      GenVertex newGV(mother, goodDaughters, primaryVertex, motherPartID);
      genVertices_.push_back(newGV);
      if (matchGenToSimVertex(newGV)) {
        genVerticesSimMatch_.push_back(newGV);
      }
    }
    if (goodDaughtersNoNu->size() >= 2) {
      GenVertex newGVNoNu(mother, goodDaughtersNoNu, primaryVertex, motherPartID);
      genVerticesNoNu_.push_back(newGVNoNu);
      if (matchGenToSimVertex(newGVNoNu)) {
        genVerticesNoNuSimMatch_.push_back(newGVNoNu);
      }
    }
  } // End loop over pruned GenParticles

  // Construct GenVertexCollections from pruned GenParticles
  unsigned int nPassingPrunedGP = 0;
  for (const reco::GenParticle& gp : prunedGenParticles_) {
    const reco::GenParticle* mother = gp.clone();

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
    nPassingPrunedGP++;

    // std::cout << "--------------------- new mother ---------------------" << std::endl;
    // std::cout << "mother pdgid = " << mother->pdgId() << std::endl;
    std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
    for (unsigned int iDau = 0; iDau < mother->numberOfDaughters(); iDau++) {
      const reco::Candidate* dau = mother->daughter(iDau)->clone();
      std::queue<const reco::Candidate*> queue;
      queue.push(dau);
      // std::cout << "------------ new queue ------------" << std::endl;
      while(!queue.empty()) {
        // std::cout << "queue size = " << queue.size() << std::endl;
        // std::cout << "front pdgid = " << queue.front()->pdgId() << ", " << queue.front()->numberOfDaughters() << " daughters" << std::endl;
        if (queue.front()->status() == 1) { // Stable outgoing particle
          if (goodGenParticle(queue.front(), genDaughterPtMin_) && queue.front()->charge() != 0) {
            goodDaughters->push_back(queue.front());
            // std::cout << "stable & pass kinematic cuts pdgid = " << queue.front()->pdgId() << std::endl;
          }
        } else {
          for (unsigned int iDau = 0; iDau < queue.front()->numberOfDaughters(); iDau++) {
            queue.push(queue.front()->daughter(iDau));
          }
        }
        queue.pop();
      }
    }

    if (goodDaughters->size() >= 1) {
      GenVertex newGVNoNu(mother, goodDaughters, primaryVertex, motherPartID);
      genVerticesFromPrunedGenNoNu_.push_back(newGVNoNu);
      if (motherPartID == 0 || motherPartID == 1) genVerticesB_.push_back(newGVNoNu);
      else if (motherPartID == 2 || motherPartID == 3) genVerticesD_.push_back(newGVNoNu);
      if (matchGenToSimVertex(newGVNoNu)) {
        genVerticesFromPrunedGenNoNuSimMatch_.push_back(newGVNoNu);
      }
    }
  }

  // Construct GenVertexCollections from TrackingVertices
  // Maybe this will give pileup information?
  for (const TrackingVertex& TV : trackingVertices_) {
    for (const GenVertexRef gvRef : TV.genVertices()) {
      if (gvRef->particles_in_size() == 0) continue; // No mother particles
      if (gvRef->particles_in_size() > 1) continue; // Too many mother particles
      if (gvRef->particles_out_size() < 2) continue; // Not a vertex

      const HepMC::GenParticle* mother = *(gvRef->particles_in_const_begin());
      if (!goodHepMCGenParticle(mother, genMotherPtMin_)) continue; // Kinematic cut on mother

      const int motherPartID = genPartID(mother->pdg_id());
      if (motherPartID < 0) continue; // Not interesting mother particle

      // Iterator beginning and end for daughters
      it_out it_begin = gvRef->particles_out_const_begin();
      it_out it_end = gvRef->particles_out_const_end();

      bool lastInstance = true;
      for (it_out dau = it_begin; dau != it_end; dau++) {
        if (genPartID((*dau)->pdg_id()) == motherPartID) {
          lastInstance = false;
          break;
        }
      }
      if (!lastInstance) continue;

      std::vector<const HepMC::GenParticle*>* goodDaughters = new std::vector<const HepMC::GenParticle*>;
      std::vector<const HepMC::GenParticle*>* goodDaughtersNoNu = new std::vector<const HepMC::GenParticle*>;
      for (it_out dau = it_begin; dau != it_end; dau++) {
        if (goodHepMCGenParticle(*dau, genDaughterPtMin_)) {
          goodDaughters->push_back(*dau);
          if (abs((*dau)->pdg_id()) != 12 && abs((*dau)->pdg_id()) != 14 && abs((*dau)->pdg_id()) != 16) {
            goodDaughtersNoNu->push_back(*dau);
          }
        }
      }

      if (goodDaughters->size() >= 2) {
        GenVertex newGV(gvRef->position(), mother, goodDaughters, primaryVertex, motherPartID);
        genVerticesFromTV_.push_back(newGV);
      }
      if (goodDaughtersNoNu->size() >= 2) {
        GenVertex newGVNoNu(gvRef->position(), mother, goodDaughtersNoNu, primaryVertex, motherPartID);
        genVerticesFromTVNoNu_.push_back(newGVNoNu);
      }
    } // End loop over HepMC GenVertexs
  } // End loop over TrackingVertexs

  // Sort collections
  // std::sort(genVertices_.begin(), genVertices_.end(), compare);
  // std::sort(genVerticesSimMatch_.begin(), genVerticesSimMatch_.end(), compare);
  // std::sort(genVerticesNoNu_.begin(), genVerticesNoNu_.end(), compare);
  // std::sort(genVerticesNoNuSimMatch_.begin(), genVerticesNoNuSimMatch_.end(), compare);

  return nPassingPrunedGP;
}


template <class P>
bool GenVertexCollectionBuilder::goodGenParticle(const P* gp, double ptCut) {

  bool pass = true;
  if (gp->pt() < ptCut) pass = false;
  if (abs(gp->eta()) > absEtaMax_) pass = false;
  return pass;
}


bool GenVertexCollectionBuilder::goodHepMCGenParticle(const HepMC::GenParticle* gp, double ptCut) {

  bool pass = true;
  double px = gp->momentum().px();
  double py = gp->momentum().py();
  double pt = TMath::Sqrt(px * px + py * py);
  if (pt < ptCut) pass = false;
  if (abs(gp->momentum().eta()) > absEtaMax_) pass = false;
  return pass;
}


template <class P>
bool GenVertexCollectionBuilder::matchGenToSimTrack(const P* gp, const SimTrack& st) {

  bool match = true;
  if (reco::deltaR(gp->eta(), gp->phi(), st.momentum().Eta(), st.momentum().Phi()) > drCut_) match = false;
  if (abs(gp->pt() - st.momentum().Pt()) / (gp->pt() + st.momentum().Pt()) > ptCut_) match = false;
  return match;
}


bool GenVertexCollectionBuilder::matchGenToSimVertex(const GenVertex& gv) {

  double nmatch = 0;
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
  return (nmatch / (double) gv.nDaughters()) >= matchFrac_;
}


int GenVertexCollectionBuilder::genPartID(int pdgId) {

  int checkPdgId = abs(pdgId);
  // B meson
  if ((checkPdgId == 521) ||
    (checkPdgId == 511) ||
    (checkPdgId == 531) ||
    (checkPdgId == 541)) return 0;
  // B baryon
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
    (checkPdgId == 5554)) return 1;
  // C meson
  else if (
    (checkPdgId == 411) ||
    (checkPdgId == 421) ||
    (checkPdgId == 431)) return 2;
  // C baryon
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
    (checkPdgId == 4444)) return 3;
  // S baryon
  // else if (
  //   (checkPdgId == 3122) ||
  //   (checkPdgId == 3222) ||
  //   (checkPdgId == 3212) ||
  //   (checkPdgId == 3112) ||
  //   (checkPdgId == 3322) ||
  //   (checkPdgId == 3312) ||
  //   (checkPdgId == 3334)) return 4;
  return -1;
}


// bool GenVertexCollectionBuilder::compare(const GenVertex& gva, const GenVertex& gvb) {

//   return gvb.d3d() < gva.d3d();
// }
