#include "../interface/GenVertexCollectionBuilder.h"


GenVertexCollectionBuilder::GenVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  genMotherPtMin_ = iConfig.getUntrackedParameter<double>("genMotherPtMin");
  genDaughterPtMin_ = iConfig.getUntrackedParameter<double>("genDaughterPtMin");
  drCut_ = iConfig.getUntrackedParameter<double>("trkMatchDrCut");
  ptCut_ = iConfig.getUntrackedParameter<double>("trkMatchPtCut");
  matchFrac_ = iConfig.getUntrackedParameter<double>("vtxMatchFrac");
}


// GenVertexCollectionBuilder::~GenVertexCollectionBuilder() {}


void GenVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::GenParticleCollection>& genParticlesToken,
    edm::EDGetTokenT<edm::SimTrackContainer>& simTracksToken,
    const reco::Vertex& primaryVertex) {

  genParticles_ = iEvent.get(genParticlesToken);
  simTracks_ = iEvent.get(simTracksToken);

  genVertices_.clear();
  genVerticesSimMatch_.clear();
  genVerticesNoNu_.clear();
  genVerticesNoNuSimMatch_.clear();

  for (unsigned int iGP = 0; iGP < genParticles_.size(); iGP++) {
    const reco::GenParticle* gp = (genParticles_.at(iGP)).clone();
    int motherPartID = genPartID(gp->pdgId());

    if (motherPartID < 0) continue; // Mother is not interesting hadron
    if (!goodGenParticle(gp, genMotherPtMin_)) continue; // Cut on mother
    if (gp->numberOfDaughters() < 2) continue; // Not a vertex

    bool lastInstance = true; // Check for last instance of interesting hadron
    for (unsigned int iDau = 0; iDau < gp->numberOfDaughters(); iDau++) {
      if (genPartID((gp->daughter(iDau))->pdgId()) == motherPartID) {
        lastInstance = false; // Not last instance of interesting hadron
        break;
      }
    }
    if (!lastInstance) continue;

    std::vector<const reco::Candidate*>* goodDaughters = new std::vector<const reco::Candidate*>;
    std::vector<const reco::Candidate*>* goodDaughtersNoNu = new std::vector<const reco::Candidate*>;
    for (unsigned int iDau = 0; iDau < gp->numberOfDaughters(); iDau++) {
      const reco::Candidate* dau = gp->daughter(iDau);
      if (goodGenParticle(dau, genDaughterPtMin_)) {
        goodDaughters->push_back(dau->clone());
        if (abs(dau->pdgId()) != 12 && abs(dau->pdgId()) != 14 && abs(dau->pdgId()) != 16) {
          goodDaughtersNoNu->push_back(dau->clone());
        }
      }
    }

    if (goodDaughters->size() >= 2) {
      GenVertex newGV(gp, goodDaughters, primaryVertex);
      genVertices_.push_back(newGV);
      if (matchGenToSimVertex(newGV)) {
        genVerticesSimMatch_.push_back(newGV);
      }
    }
    if (goodDaughtersNoNu->size() >= 2) {
      GenVertex newGVNoNu(gp, goodDaughtersNoNu, primaryVertex);
      genVerticesNoNu_.push_back(newGVNoNu);
      if (matchGenToSimVertex(newGVNoNu)) {
        genVerticesNoNuSimMatch_.push_back(newGVNoNu);
      }
    }
  } // End loop over all gen particles
}


template <class P>
bool GenVertexCollectionBuilder::goodGenParticle(const P* gp, double ptCut) {

  bool pass = true;
  if (gp->pt() < ptCut) pass = false;
  if (abs(gp->eta()) > absEtaMax_) pass = false;
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
  for (const reco::Candidate* dau : *(gv.daughters())) {
    for (const SimTrack& st : simTracks_) {
      if (matchGenToSimTrack(dau, st)) {
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
