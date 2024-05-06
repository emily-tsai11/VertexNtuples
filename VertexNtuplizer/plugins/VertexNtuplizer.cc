// -*- C++ -*-
//
// Package:    VertexNtuples/VertexNtuplizer
// Class:      VertexNtuplizer
//
/**\class VertexNtuplizer VertexNtuplizer.cc VertexNtuples/VertexNtuplizer/plugins/VertexNtuplizer.cc

 Description: ntuplizer for vertexing studies

 Implementation:
     EDAnalyzer in CMSSW in C++
*/
//
// Original Author:  Emily Minyun Tsai
//         Created:  Sat, 04 May 2024 12:05:41 GMT
//
//


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "../interface/GenVertex.h"
#include "../interface/SecondaryVertex.h"
#include "../interface/RecoJet.h"
#include "../interface/GenJet.h"
#include "../interface/GenVertexCollectionBuilder.h"
#include "../interface/SecondaryVertexCollectionBuilder.h"
#include "../interface/RecoJetCollectionBuilder.h"
#include "../interface/GenJetCollectionBuilder.h"


class VertexNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

  public:

    explicit VertexNtuplizer(const edm::ParameterSet&);
    ~VertexNtuplizer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:

    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
    edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDTimingToken_;
    edm::EDGetTokenT<unsigned int> IVFclustersToken_;
    edm::EDGetTokenT<unsigned int> IVFclustersMTDTimingToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeErrorMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeQualityMapToken_;
    edm::EDGetTokenT<pat::JetCollection> jetsToken_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken_;

    GenVertexCollectionBuilder* gvc_;
    SecondaryVertexCollectionBuilder* svc_;
    RecoJetCollectionBuilder* rjc_;
    GenJetCollectionBuilder* gjc_;

    std::map<TString, TH1F*> histos_;
};


static unsigned int nbins_ = 80;


VertexNtuplizer::VertexNtuplizer(const edm::ParameterSet& iConfig) :
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"))),
    simTracksToken_(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simTracks"))),
    primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
    secondaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices"))),
    secondaryVerticesMTDTimingToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDTiming"))),
    IVFclustersToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("IVFclusters"))),
    IVFclustersMTDTimingToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("IVFclustersMTDTiming"))),
    trackTimeValueMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeValueMap"))),
    trackTimeErrorMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeErrorMap"))),
    trackTimeQualityMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeQualityMap"))),
    jetsToken_(consumes<pat::JetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
    genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJets"))),
    genJetsFlavourInfoToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJetsFlavourInfo"))) {

  gvc_ = new GenVertexCollectionBuilder(iConfig);
  svc_ = new SecondaryVertexCollectionBuilder(iConfig);
  rjc_ = new RecoJetCollectionBuilder(iConfig);
  gjc_ = new GenJetCollectionBuilder(iConfig);

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  histos_["nGV"] = fs->make<TH1F>("nGV", "nGV", 28, 2, 30);
  histos_["nGVs"] = fs->make<TH1F>("nGVs", "nGVs", 28, 2, 30);
  histos_["nGVn"] = fs->make<TH1F>("nGVn", "nGVn", 28, 2, 30);
  histos_["nGVns"] = fs->make<TH1F>("nGVns", "nGVns", 28, 2, 30);

  histos_["nSV"] = fs->make<TH1F>("nSV", "nSV", 28, 2, 30);
  histos_["nSVt"] = fs->make<TH1F>("nSVt", "nSVt", 28, 2, 30);

  histos_["nC"] = fs->make<TH1F>("nC", "nC", 198, 2, 200);
  histos_["nCt"] = fs->make<TH1F>("nCt", "nCt", 198, 2, 200);

  histos_["nRJ"] = fs->make<TH1F>("nRJ", "nRJ", 20, 0, 20);
  histos_["nRJg"] = fs->make<TH1F>("nRJg", "nRJg", 20, 0, 20);

  histos_["nGJ"] = fs->make<TH1F>("nGJ", "nGJ", 20, 0, 20);
  histos_["nGJr"] = fs->make<TH1F>("nGJr", "nGJr", 20, 0, 20);
}


VertexNtuplizer::~VertexNtuplizer() {}


void VertexNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;

  const reco::VertexCollection primaryVertices = iEvent.get(primaryVerticesToken_);
  // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
  const reco::Vertex& primaryVertex = primaryVertices.at(0); // Most likely the signal vertex

  gvc_->build(iEvent, genParticlesToken_, simTracksToken_, primaryVertex);
  svc_->build(iEvent, secondaryVerticesToken_, secondaryVerticesMTDTimingToken_, primaryVertex,
      trackTimeValueMapToken_, trackTimeErrorMapToken_, trackTimeQualityMapToken_);
  rjc_->build(iEvent, jetsToken_, genJetsFlavourInfoToken_);
  gjc_->build(iEvent, genJetsToken_, genJetsFlavourInfoToken_, jetsToken_);

  GenVertexCollection genVertices = gvc_->getGenVertexCollection();
  GenVertexCollection genVerticesSimMatch = gvc_->getGenVertexSimMatchCollection();
  GenVertexCollection genVerticesNoNu = gvc_->getGenVertexNoNuCollection();
  GenVertexCollection genVerticesNoNuSimMatch = gvc_->getGenVertexNoNuSimMatchCollection();

  histos_["nGV"]->Fill(genVertices.size());
  histos_["nGVs"]->Fill(genVerticesSimMatch.size());
  histos_["nGVn"]->Fill(genVerticesNoNu.size());
  histos_["nGVns"]->Fill(genVerticesNoNuSimMatch.size());

  SecondaryVertexCollection secondaryVertices = svc_->getSecondaryVertexCollection();
  SecondaryVertexCollection secondaryVerticesMTDTiming = svc_->getSecondaryVertexCollectionMTDTiming();

  histos_["nSV"]->Fill(secondaryVertices.size());
  histos_["nSVt"]->Fill(secondaryVerticesMTDTiming.size());

  unsigned int nC = iEvent.get(IVFclustersToken_);
  unsigned int nCt = iEvent.get(IVFclustersMTDTimingToken_);

  histos_["nC"]->Fill(nC);
  histos_["nCt"]->Fill(nCt);

  RecoJetCollection recoJets = rjc_->recoJetCollection();
  RecoJetCollection recoJetsGenMatch = rjc_->recoJetGenMatchCollection();

  histos_["nRJ"]->Fill(recoJets.size());
  histos_["nRJg"]->Fill(recoJetsGenMatch.size());

  GenJetCollection genJets = gjc_->getGenJetCollection();
  GenJetCollection genJetsRecoMatch = gjc_->getGenJetRecoMatchCollection();

  histos_["nGJ"]->Fill(genJets.size());
  histos_["nGJr"]->Fill(genJetsRecoMatch.size());
}


void VertexNtuplizer::beginJob() {}


void VertexNtuplizer::endJob() {}


void VertexNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(VertexNtuplizer);
