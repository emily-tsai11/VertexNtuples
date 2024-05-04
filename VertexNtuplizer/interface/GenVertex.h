#ifndef GEN_VERTEX
#define GEN_VERTEX

class GenVertex {

  public:

    GenVertex(const reco::GenParticle* mother, std::vector<const reco::Candidate*>* daughters,
        edm::SimTrackContainer simTracks, float matchFrac, float drCut, float ptCut);

    const float x()   const { return daughters_->at(0)->vx(); }
    const float y()   const { return daughters_->at(0)->vy(); }
    const float z()   const { return daughters_->at(0)->vz(); }
    const float pt()  const { return mother_->pt(); }
    const float eta() const { return mother_->eta(); }
    const float phi() const { return mother_->phi(); }
    const unsigned int nDaughters()     const { return daughters_->size(); }
    const unsigned int nDaughtersNoNu() const { return daughtersNoNu_->size(); }
    const int motherPdgId() const { return mother_->pdgId(); }

    const bool isSimMatched()     const { return isSimMatched_; }
    const bool isSimMatchedNoNu() const { return isSimMatchedNoNu_; }

    const reco::GenParticle* mother() const { return mother_; }
    const std::vector<const reco::Candidate*>* daughters() const { return daughters_; }
    const std::vector<const reco::Candidate*>* daughtersNoNu() const { return daughters_; }

    void print() {
      std::cout << "GenVertex:" << std::endl;
      std::cout << "    vertex x        = " << x() << std::endl;
      std::cout << "    vertex y        = " << y() << std::endl;
      std::cout << "    vertex z        = " << z() << std::endl;
      std::cout << "    vertex pt       = " << pt() << std::endl;
      std::cout << "    vertex eta      = " << eta() << std::endl;
      std::cout << "    vertex phi      = " << phi() << std::endl;
      std::cout << "    nDaughters      = " << nDaughters() << std::endl;
      std::cout << "    mother pdg id   = " << motherPdgId() << std::endl;

      std::cout << "    daughter pdgIds = ";
      for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
        std::cout << daughters_->at(iDau)->pdgId();
        if (iDau < nDaughters() - 1) std::cout << ", ";
      }
      std::cout << std::endl;
      std::cout << "    daughter pts   = ";
      for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
        std::cout << daughters_->at(iDau)->pt();
        if (iDau < nDaughters() - 1) std::cout << ", ";
      }
      std::cout << std::endl;
      std::cout << "    daughter etas   = ";
      for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
        std::cout << daughters_->at(iDau)->eta();
        if (iDau < nDaughters() - 1) std::cout << ", ";
      }
      std::cout << std::endl;
      std::cout << "    daughter phis   = ";
      for (unsigned int iDau = 0; iDau < nDaughters(); iDau++) {
        std::cout << daughters_->at(iDau)->phi();
        if (iDau < nDaughters() - 1) std::cout << ", ";
      }
      std::cout << std::endl;
    }

  private:

    const reco::GenParticle* mother_;
    std::vector<const reco::Candidate*>* daughters_;
    std::vector<const reco::Candidate*>* daughtersNoNu_;
    bool isSimMatched_;
    bool isSimMatchedNoNu_;
};

#endif
