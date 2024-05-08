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
#include "../interface/VertexMatcher.h"


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
    VertexMatcher* matcher_;

    std::vector<TString> gv_names_;
    std::vector<TString> sv_names_;
    std::vector<TString> rj_names_;
    std::vector<TString> gj_names_;

    std::map<TString, TH1F*> histos_;
};


static unsigned int nbins_ = 100;
static unsigned int nvtx_ = 30;
static unsigned int nclus_ = 200;
static unsigned int njet_ = 20;


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
  matcher_ = new VertexMatcher(iConfig);

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  gv_names_.push_back("gv"); // GenVertex
  gv_names_.push_back("gvs"); // GenVertex w/SIM match
  gv_names_.push_back("gvn"); // GenVertex without neutrino daughters
  gv_names_.push_back("gvns"); // GenVertex without neutrino daughters w/SIM match

  sv_names_.push_back("sv"); // SecondaryVertex
  sv_names_.push_back("svt"); // SecondaryVertex w/time range cut

  rj_names_.push_back("rj"); // RecoJet
  rj_names_.push_back("rjg"); // RecoJet w/GEN match

  gj_names_.push_back("gj"); // GenJet
  gj_names_.push_back("gjr"); // GenJet w/reco match

  std::vector<TString> objs_;
  for (TString obj : gv_names_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk");
  }
  for (TString obj : sv_names_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk");
  }
  for (TString obj : rj_names_) objs_.push_back(obj);
  for (TString obj : gj_names_) objs_.push_back(obj);

  std::map<TString, std::vector<float>> vars_ = {
    std::make_pair("tval", std::vector<float>{(float) nbins_, -0.8, 0.8}),
    std::make_pair("terr", std::vector<float>{(float) nbins_, 0.0, 0.1}),
    std::make_pair("tsig", std::vector<float>{(float) nbins_, -50.0, 50.0}),
    std::make_pair("tqual", std::vector<float>{(float) nbins_, 0.0, 1.0}),
    std::make_pair("tavg", std::vector<float>{(float) nbins_, -0.8, 0.8}),
    std::make_pair("trange", std::vector<float>{(float) nbins_, 0.0, 0.8}),
    std::make_pair("x", std::vector<float>{(float) nbins_, -1.0, 1.0}),
    std::make_pair("y", std::vector<float>{(float) nbins_, -1.0, 1.0}),
    std::make_pair("z", std::vector<float>{(float) nbins_, -20.0, 20.0}),
    std::make_pair("xerr", std::vector<float>{(float) nbins_, -7.0, 7.0}),
    std::make_pair("yerr", std::vector<float>{(float) nbins_, -7.0, 7.0}),
    std::make_pair("zerr", std::vector<float>{(float) nbins_, -7.0, 7.0}),
    std::make_pair("xres", std::vector<float>{(float) nbins_, -0.15, 0.15}),
    std::make_pair("yres", std::vector<float>{(float) nbins_, -0.15, 0.15}),
    std::make_pair("zres", std::vector<float>{(float) nbins_, -0.15, 0.15}),
    std::make_pair("pt", std::vector<float>{(float) nbins_, 0.0, 200.0}),
    std::make_pair("pt2", std::vector<float>{(float) nbins_, 0.0, 200.0}),
    std::make_pair("eta", std::vector<float>{(float) nbins_, -3.1, 3.1}),
    std::make_pair("phi", std::vector<float>{(float) nbins_, -3.15, 3.15}),
    std::make_pair("dxy", std::vector<float>{(float) nbins_, 0.0, 10.0}),
    std::make_pair("dz", std::vector<float>{(float) nbins_, 0.0, 20.0}),
    std::make_pair("d3d", std::vector<float>{(float) nbins_, 0.0, 20.0}),
    std::make_pair("dxyerr", std::vector<float>{(float) nbins_, 0.0, 0.1}),
    std::make_pair("dzerr", std::vector<float>{(float) nbins_, 0.0, 0.1}),
    std::make_pair("d3derr", std::vector<float>{(float) nbins_, 0.0, 0.1}),
    std::make_pair("dxysig", std::vector<float>{(float) nbins_, 0.0, 3000.0}),
    std::make_pair("dzsig", std::vector<float>{(float) nbins_, 0.0, 1000.0}),
    std::make_pair("d3dsig", std::vector<float>{(float) nbins_, 0.0, 1000.0}),
    std::make_pair("charge", std::vector<float>{(float) 5, -2.0, 3.0}),
    std::make_pair("motherPdgId", std::vector<float>{(float) nbins_, -5560.0, 5560.0}),
    std::make_pair("pdgIdBin", std::vector<float>{4, 0.0, 4.0}),
    std::make_pair("hadFlav", std::vector<float>{7, 0.0, 7.0}),
    std::make_pair("chi2", std::vector<float>{(float) nbins_, 0.0, 100.0}),
    std::make_pair("ndof", std::vector<float>{(float) nbins_, 0.0, 10.0}),
    std::make_pair("chi2dof", std::vector<float>{(float) nbins_, 0.0, 10.0}),
    std::make_pair("ntrk", std::vector<float>{(float) nbins_, 0.0, 10.0}) // Daughters for GenVertex
  };

  for (TString gv_name : gv_names_) {
    TString name = "n" + gv_name;
    histos_[name] = fs->make<TH1F>(name, name, nvtx_, 0, nvtx_);
  }
  for (TString sv_name : sv_names_) {
    TString name = "n" + sv_name;
    histos_[name] = fs->make<TH1F>(name, name, nvtx_, 0, nvtx_);
  }
  histos_["nc"] = fs->make<TH1F>("nc", "nc", nclus_, 0, nclus_);
  histos_["nct"] = fs->make<TH1F>("nct", "nct", nclus_, 0, nclus_);
  for (TString rj_name : rj_names_) {
    TString name = "n" + rj_name;
    histos_[name] = fs->make<TH1F>(name, name, njet_, 0, njet_);
  }
  for (TString gj_name : gj_names_) {
    TString name = "n" + gj_name;
    histos_[name] = fs->make<TH1F>(name, name, njet_, 0, njet_);
  }

  for (TString obj : objs_) {
    for (const auto& iter : vars_) {
      TString name = obj + "_" + iter.first;
      histos_[name] = fs->make<TH1F>(name, name, iter.second[0], iter.second[1], iter.second[2]);
    }
  }

  for (TString obj1 : objs_) {
    for (TString obj2 : objs_) {
      if (obj1.BeginsWith("gv") && obj2.BeginsWith("gv")) continue;
      if (obj1.BeginsWith("sv") && obj2.BeginsWith("sv")) continue;
      if (obj1.BeginsWith("rj") && obj2.BeginsWith("rj")) continue;
      if (obj1.BeginsWith("gj") && obj2.BeginsWith("gj")) continue;
      for (const auto& iter : vars_) {
        TString name = obj1 + "_" + obj2 + "_" + iter.first;
        histos_[name] = fs->make<TH1F>(name, name, iter.second[0], iter.second[1], iter.second[2]);
      }
    }
  }
}


VertexNtuplizer::~VertexNtuplizer() {}


void VertexNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  const reco::VertexCollection primaryVertices = iEvent.get(primaryVerticesToken_);
  // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
  const reco::Vertex& primaryVertex = primaryVertices.at(0); // Most likely the signal vertex

  gvc_->build(iEvent, genParticlesToken_, simTracksToken_, primaryVertex);
  svc_->build(iEvent, secondaryVerticesToken_, secondaryVerticesMTDTimingToken_, primaryVertex,
      trackTimeValueMapToken_, trackTimeErrorMapToken_, trackTimeQualityMapToken_);
  rjc_->build(iEvent, jetsToken_, genJetsFlavourInfoToken_);
  gjc_->build(iEvent, genJetsToken_, genJetsFlavourInfoToken_, jetsToken_);

  std::vector<GenVertexCollection> GVCollections;
  GVCollections.push_back(gvc_->getGenVertexCollection());
  GVCollections.push_back(gvc_->getGenVertexSimMatchCollection());
  GVCollections.push_back(gvc_->getGenVertexNoNuCollection());
  GVCollections.push_back(gvc_->getGenVertexNoNuSimMatchCollection());
  for (unsigned int iColl = 0; iColl < GVCollections.size(); iColl++) {
    histos_["n" + gv_names_.at(iColl)]->Fill(GVCollections.at(iColl).size());
    for (GenVertex& gv : GVCollections.at(iColl)) gv.fill(histos_, gv_names_.at(iColl));
  }

  std::vector<SecondaryVertexCollection> SVCollections;
  SVCollections.push_back(svc_->getSecondaryVertexCollection());
  SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDTiming());
  for (unsigned int iColl = 0; iColl < SVCollections.size(); iColl++) {
    histos_["n" + sv_names_.at(iColl)]->Fill(SVCollections.at(iColl).size());
    for (SecondaryVertex& sv : SVCollections.at(iColl)) sv.fill(histos_, sv_names_.at(iColl));
  }

  unsigned int nC = iEvent.get(IVFclustersToken_);
  unsigned int nCt = iEvent.get(IVFclustersMTDTimingToken_);
  histos_["nc"]->Fill(nC);
  histos_["nct"]->Fill(nCt);

  std::vector<RecoJetCollection> RJCollections;
  RJCollections.push_back(rjc_->getRecoJetCollection());
  RJCollections.push_back(rjc_->getRecoJetGenMatchCollection());
  for (unsigned int iColl = 0; iColl < RJCollections.size(); iColl++) {
    histos_["n" + rj_names_.at(iColl)]->Fill(RJCollections.at(iColl).size());
    for (RecoJet& rj : RJCollections.at(iColl)) rj.fill(histos_, rj_names_.at(iColl));
  }

  std::vector<GenJetCollection> GJCollections;
  GJCollections.push_back(gjc_->getGenJetCollection());
  GJCollections.push_back(gjc_->getGenJetRecoMatchCollection());
  for (unsigned int iColl = 0; iColl < GJCollections.size(); iColl++) {
    histos_["n" + gj_names_.at(iColl)]->Fill(GJCollections.at(iColl).size());
    for (GenJet& gj : GJCollections.at(iColl)) gj.fill(histos_, gj_names_.at(iColl));
  }

  for (unsigned int iGVs = 0; iGVs < GVCollections.size(); iGVs++) {
    for (unsigned int iSVs = 0; iSVs < SVCollections.size(); iSVs++) {
      for (GenVertex& gv : GVCollections.at(iGVs)) {
        for (SecondaryVertex& sv : SVCollections.at(iSVs)) {
          if (matcher_->match(gv, sv, TRACK)) {
            TString gv_name = gv_names_.at(iGVs) + "_" + sv_names_.at(iSVs);
            gv.fill(histos_, gv_name);
            TString sv_name = sv_names_.at(iSVs) + "_" + gv_names_.at(iGVs);
            sv.fill(histos_, sv_name);
          }
        }
      }
    }
  }
}


void VertexNtuplizer::beginJob() {}


void VertexNtuplizer::endJob() {}


void VertexNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(VertexNtuplizer);
