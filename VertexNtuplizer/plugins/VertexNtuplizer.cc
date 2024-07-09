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
// #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

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

    edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
    // edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticlesToken_;
    edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
    // edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
    edm::EDGetTokenT<TrackingVertexCollection> trackingVerticesToken_;
    edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDBSToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDBS4Token_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDPVToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDBSToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDBS4Token_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeBSValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeBSErrorMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeBSQualityMapToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDPVToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimePVValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimePVErrorMapToken_;
    // edm::EDGetTokenT<edm::ValueMap<float>> trackTimePVQualityMapToken_;
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

    std::map<TString, std::vector<float>> vars1_;
    std::map<TString, std::vector<float>> vars2_;

    std::map<TString, TH1F*> histos1_;
    std::map<TString, TH2F*> histos2_;

    bool scanCuts_;
    double trkMatchDrCut_;
    double trkMatchPtCut_;
};


VertexNtuplizer::VertexNtuplizer(const edm::ParameterSet& iConfig) :
    prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedGenParticles"))),
    // packedGenParticlesToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("packedGenParticles"))),
    simTracksToken_(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simTracks"))),
    // trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingParticles"))),
    trackingVerticesToken_(consumes<TrackingVertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingVertices"))),
    primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
    nIVFClustersToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClusters"))),
    nIVFClustersMTDBSToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDBS"))),
    nIVFClustersMTDBS4Token_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDBS4"))),
    nIVFClustersMTDPVToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDPV"))),
    secondaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices"))),
    secondaryVerticesMTDBSToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDBS"))),
    secondaryVerticesMTDBS4Token_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDBS4"))),
    trackTimeBSValueMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeBSValueMap"))),
    trackTimeBSErrorMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeBSErrorMap"))),
    trackTimeBSQualityMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeBSQualityMap"))),
    secondaryVerticesMTDPVToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDPV"))),
    trackTimePVValueMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimePVValueMap"))),
    trackTimePVErrorMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimePVErrorMap"))),
    // trackTimePVQualityMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimePVQualityMap"))),
    jetsToken_(consumes<pat::JetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
    genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJets"))),
    genJetsFlavourInfoToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJetsFlavourInfo"))) {

  scanCuts_ = iConfig.getUntrackedParameter<bool>("scanCuts");
  trkMatchDrCut_ = iConfig.getUntrackedParameter<double>("recoTrkMatchDrCut");
  trkMatchPtCut_ = iConfig.getUntrackedParameter<double>("recoTrkMatchPtCut");

  gvc_ = new GenVertexCollectionBuilder(iConfig);
  svc_ = new SecondaryVertexCollectionBuilder(iConfig);
  rjc_ = new RecoJetCollectionBuilder(iConfig);
  gjc_ = new GenJetCollectionBuilder(iConfig);
  matcher_ = new VertexMatcher(iConfig);

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  gv_names_.push_back("gv"); // GenVertex
  gv_names_.push_back("gvs"); // GenVertex w/SIM match
  gv_names_.push_back("gvn"); // GenVertex w/out neutrino daughters
  gv_names_.push_back("gvns"); // GenVertex w/out neutrino daughters w/SIM match
  gv_names_.push_back("gvpn"); // GenVertex constructed from pruned GenParticles w/out neutrino daughters
  gv_names_.push_back("gvpnB"); // GenVertex constructed from pruned GenParticles w/out neutrino daughters
  gv_names_.push_back("gvpnD"); // GenVertex constructed from pruned GenParticles w/out neutrino daughters
  gv_names_.push_back("gvpns"); // GenVertex constructed from pruned GenParticles w/out neutrino daughters w/SIM match
  gv_names_.push_back("gvt"); // GenVertex constructed from TrackingVertexs
  gv_names_.push_back("gvtn"); // GenVertex constructed from TrackingVertexs w/out neutrino daughters

  sv_names_.push_back("sv"); // SecondaryVertex
  sv_names_.push_back("svbs"); // SecondaryVertex w/track time extrapolated to the beam spot
  sv_names_.push_back("svbs4"); // SecondaryVertex w/range cut on track time extrapolated to the beam spot
  sv_names_.push_back("svpv"); // SecondaryVertex w/track time extrapolated to the primary vertex

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

  const unsigned int nbins_ = 100;
  const unsigned int ngv_ = 15;
  const unsigned int nsv_ = 200;
  const unsigned int njet_ = 20;

  vars1_["tval"] = std::vector<float>{(float) nbins_, -0.8, 0.8};
  vars1_["terr"] = std::vector<float>{(float) nbins_, 0.0, 0.1};
  vars1_["tsig"] = std::vector<float>{(float) nbins_, -50.0, 50.0};
  vars1_["tqual"] = std::vector<float>{(float) nbins_, 0.0, 1.0};
  vars1_["tavg"] = std::vector<float>{(float) nbins_, -0.8, 0.8};
  vars1_["trange"] = std::vector<float>{(float) nbins_, 0.0, 0.8};
  vars1_["pt"] = std::vector<float>{(float) nbins_, 0.0, 200.0};
  vars1_["pt2"] = std::vector<float>{(float) nbins_, 0.0, 1000.0};
  vars1_["eta"] = std::vector<float>{(float) nbins_, -3.1, 3.1};
  vars1_["phi"] = std::vector<float>{(float) nbins_, -3.15, 3.15};
  vars1_["charge"] = std::vector<float>{3, -1.0, 2.0};
  vars1_["motherPdgId"] = std::vector<float>{(float) nbins_, -5560.0, 5560.0};
  vars1_["pdgIdBin"] = std::vector<float>{5, 0.0, 5.0};
  vars1_["hadFlav"] = std::vector<float>{7, 0.0, 7.0};
  vars1_["chi2"] = std::vector<float>{(float) nbins_, 0.0, 100.0};
  vars1_["ndof"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["chi2dof"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["ntrk"] = std::vector<float>{10, 0.0, 10.0}; // Daughters for GenVertex
  // From origin
  vars1_["x"] = std::vector<float>{(float) nbins_, -1.0, 1.0};
  vars1_["y"] = std::vector<float>{(float) nbins_, -1.0, 1.0};
  vars1_["z"] = std::vector<float>{(float) nbins_, -20.0, 20.0};
  vars1_["xerr"] = std::vector<float>{(float) nbins_, 0.0, 0.5};
  vars1_["yerr"] = std::vector<float>{(float) nbins_, 0.0, 0.5};
  vars1_["zerr"] = std::vector<float>{(float) nbins_, 0.0, 0.5};
  // From primary vertex
  vars1_["dxy"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["dz"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["d3d"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["dxyerr"] = std::vector<float>{(float) nbins_, 0.0, 0.05};
  vars1_["dzerr"] = std::vector<float>{(float) nbins_, 0.0, 0.1};
  vars1_["d3derr"] = std::vector<float>{(float) nbins_, 0.0, 0.1};
  vars1_["dxysig"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["dzsig"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["d3dsig"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  // Between matched GenVertex and SecondaryVertex
  vars1_["xres"] = std::vector<float>{(float) nbins_, -0.15, 0.15};
  vars1_["yres"] = std::vector<float>{(float) nbins_, -0.15, 0.15};
  vars1_["zres"] = std::vector<float>{(float) nbins_, -0.15, 0.15};
  vars1_["xpull"] = std::vector<float>{(float) nbins_, -6.5, 6.5};
  vars1_["ypull"] = std::vector<float>{(float) nbins_, -6.5, 6.5};
  vars1_["zpull"] = std::vector<float>{(float) nbins_, -6.5, 6.5};
  vars1_["matchdxy"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["matchd3d"] = std::vector<float>{(float) nbins_, 0.0, 20.0};
  vars1_["matchdxyerr"] = std::vector<float>{(float) nbins_, 0.0, 0.1};
  vars1_["matchd3derr"] = std::vector<float>{(float) nbins_, 0.0, 0.1};
  vars1_["matchdxysig"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["matchd3dsig"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_["deltaR"] = std::vector<float>{(float) nbins_, 0, 4.0};
  vars1_["ptResNorm"] = std::vector<float>{(float) nbins_, 0.0, 1.0};

  vars2_["eta_tval"] = std::vector<float>{(float) nbins_, -3.1, 3.1, (float) nbins_, -0.8, 0.8};
  vars2_["eta_terr"] = std::vector<float>{(float) nbins_, -3.1, 3.1, (float) nbins_, 0.0, 0.1};
  vars2_["eta_tsig"] = std::vector<float>{(float) nbins_, -3.1, 3.1, (float) nbins_, -50.0, 50.0};
  vars2_["eta_tqual"] = std::vector<float>{(float) nbins_, -3.1, 3.1, (float) nbins_, 0.0, 1.0};
  vars2_["eta_tavg"] = std::vector<float>{(float) nbins_, -3.1, 3.1, (float) nbins_, -0.8, 0.8};
  vars2_["eta_trange"] = std::vector<float>{(float) nbins_, -3.1, 3.1, (float) nbins_, 0.0, 0.8};
  vars2_["trange_pt"] = std::vector<float>{(float) nbins_, 0.0, 0.8, (float) nbins_, 0.0, 200.0};
  vars2_["trange_pt2"] = std::vector<float>{(float) nbins_, 0.0, 0.8, (float) nbins_, 0.0, 200.0};
  vars2_["trange_dxy"] = std::vector<float>{(float) nbins_, 0.0, 0.8, (float) nbins_, 0.0, 10.0};
  vars2_["trange_dxysig"] = std::vector<float>{(float) nbins_, 0.0, 0.8, (float) nbins_, 0.0, 10.0};
  vars2_["trange_d3d"] = std::vector<float>{(float) nbins_, 0.0, 0.8, (float) nbins_, 0.0, 20.0};
  vars2_["trange_d3dsig"] = std::vector<float>{(float) nbins_, 0.0, 0.8, (float) nbins_, 0.0, 10.0};
  // Between matched GenVertex and SecondaryVertex
  vars2_["matchdxy_dxy"] = std::vector<float>{(float) nbins_, 0.0, 10.0, (float) nbins_, 0.0, 10.0};
  vars2_["matchdxy_d3d"] = std::vector<float>{(float) nbins_, 0.0, 10.0, (float) nbins_, 0.0, 10.0};
  vars2_["matchd3d_dxy"] = std::vector<float>{(float) nbins_, 0.0, 10.0, (float) nbins_, 0.0, 10.0};
  vars2_["matchd3d_d3d"] = std::vector<float>{(float) nbins_, 0.0, 10.0, (float) nbins_, 0.0, 10.0};
  vars2_["deltaR_ptResNorm"] = std::vector<float>{(float) nbins_, 0.0, 4.0, (float) nbins_, 0.0, 1.0};

  // Count histograms
  for (TString gv_name : gv_names_) {
    TString name = "n" + gv_name;
    histos1_[name] = fs->make<TH1F>(name, name, ngv_, 0, ngv_);
    histos1_[name]->Sumw2();
  }
  for (TString sv_name : sv_names_) {
    TString name = "n" + sv_name;
    histos1_[name] = fs->make<TH1F>(name, name, nsv_, 0, nsv_);
    histos1_[name]->Sumw2();
  }
  histos1_["nc"] = fs->make<TH1F>("nc", "nc", nsv_, 0, nsv_);
  histos1_["nc"]->Sumw2();
  histos1_["ncbs"] = fs->make<TH1F>("ncbs", "ncbs", nsv_, 0, nsv_);
  histos1_["ncbs"]->Sumw2();
  histos1_["ncbs4"] = fs->make<TH1F>("ncbs4", "ncbs4", nsv_, 0, nsv_);
  histos1_["ncbs4"]->Sumw2();
  histos1_["ncpv"] = fs->make<TH1F>("ncpv", "ncpv", nsv_, 0, nsv_);
  histos1_["ncpv"]->Sumw2();
  for (TString rj_name : rj_names_) {
    TString name = "n" + rj_name;
    histos1_[name] = fs->make<TH1F>(name, name, njet_, 0, njet_);
    histos1_[name]->Sumw2();
  }
  for (TString gj_name : gj_names_) {
    TString name = "n" + gj_name;
    histos1_[name] = fs->make<TH1F>(name, name, njet_, 0, njet_);
    histos1_[name]->Sumw2();
  }
  // Do we need this other than for comparison purposes?
  histos1_["nPrunedGPs"] = fs->make<TH1F>("nPrunedGPs", "nPrunedGPs", 13, 0, 13);

  // 1D histograms
  for (TString obj : objs_) {
    for (const auto& iter : vars1_) {
      TString name = obj + "_" + iter.first;
      histos1_[name] = fs->make<TH1F>(name, name, iter.second[0], iter.second[1], iter.second[2]);
      histos1_[name]->Sumw2();
    }
  }

  // 2D histograms
  for (TString obj : sv_names_) {
    for (const auto& iter : vars2_) {
      TString name = obj + "_" + iter.first;
      histos2_[name] = fs->make<TH2F>(name, name, iter.second[0], iter.second[1], iter.second[2], iter.second[3], iter.second[4], iter.second[5]);
      histos2_[name]->Sumw2();
      TString trkName = obj + "_trk_" + iter.first;
      histos2_[trkName] = fs->make<TH2F>(trkName, trkName, iter.second[0], iter.second[1], iter.second[2], iter.second[3], iter.second[4], iter.second[5]);
      histos2_[trkName]->Sumw2();
    }
  }

  // Matching histograms
  for (TString obj1 : objs_) {
    for (TString obj2 : objs_) {
      if (obj1.BeginsWith("gv") && obj2.BeginsWith("gv")) continue;
      if (obj1.BeginsWith("sv") && obj2.BeginsWith("sv")) continue;
      if (obj1.BeginsWith("rj") && obj2.BeginsWith("rj")) continue;
      if (obj1.BeginsWith("gj") && obj2.BeginsWith("gj")) continue;
      if (obj1.EndsWith("_trk") && obj2.EndsWith("_trk")) continue;
      if (obj1.EndsWith("_trk")) continue;

      // 1D histograms
      for (const auto& iter : vars1_) {
        TString name = obj1 + "_" + obj2 + "_" + iter.first;
        histos1_[name] = fs->make<TH1F>(name, name, iter.second[0], iter.second[1], iter.second[2]);
        histos1_[name]->Sumw2();
      }

      // 2D histograms
      if (obj1.BeginsWith("gv") || obj1.BeginsWith("sv")) {
        for (const auto& iter : vars2_) {
          TString name = obj1 + "_" + obj2 + "_" + iter.first;
          histos2_[name] = fs->make<TH2F>(name, name, iter.second[0], iter.second[1], iter.second[2], iter.second[3], iter.second[4], iter.second[5]);
          histos2_[name]->Sumw2();
        }
      }
    }
  }
}


VertexNtuplizer::~VertexNtuplizer() {}


void VertexNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  const reco::VertexCollection primaryVertices = iEvent.get(primaryVerticesToken_);
  // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
  const reco::Vertex& primaryVertex = primaryVertices.at(0); // Most likely the signal vertex

  unsigned int nPassingPrunedGP = gvc_->build(iEvent, prunedGenParticlesToken_,
    // packedGenParticlesToken_,
    simTracksToken_,
    // trackingParticlesToken_,
    trackingVerticesToken_, primaryVertex);
  svc_->build(iEvent,
      secondaryVerticesToken_,
      secondaryVerticesMTDBSToken_,
      secondaryVerticesMTDBS4Token_,
      trackTimeBSValueMapToken_, trackTimeBSErrorMapToken_, trackTimeBSQualityMapToken_,
      secondaryVerticesMTDPVToken_,
      trackTimePVValueMapToken_, trackTimePVErrorMapToken_,
      // trackTimePVQualityMapToken_,
      primaryVertex);
  rjc_->build(iEvent, jetsToken_, genJetsFlavourInfoToken_);
  gjc_->build(iEvent, genJetsToken_, genJetsFlavourInfoToken_, jetsToken_);

  std::vector<GenVertexCollection> GVCollections;
  GVCollections.push_back(gvc_->getGenVertexCollection());
  GVCollections.push_back(gvc_->getGenVertexSimMatchCollection());
  GVCollections.push_back(gvc_->getGenVertexNoNuCollection());
  GVCollections.push_back(gvc_->getGenVertexNoNuSimMatchCollection());
  GVCollections.push_back(gvc_->getGenVertexFromPrunedGenNoNu());
  GVCollections.push_back(gvc_->getGenVertexB());
  GVCollections.push_back(gvc_->getGenVertexD());
  GVCollections.push_back(gvc_->getGenVertexFromPrunedGenNoNuSimMatch());
  GVCollections.push_back(gvc_->getGenVertexFromTV());
  GVCollections.push_back(gvc_->getGenVertexFromTVNoNu());
  for (unsigned int iColl = 0; iColl < GVCollections.size(); iColl++) {
    histos1_["n" + gv_names_.at(iColl)]->Fill(GVCollections.at(iColl).size());
    for (GenVertex& gv : GVCollections.at(iColl)) gv.fill(histos1_, gv_names_.at(iColl));
  }
  histos1_["nPrunedGPs"]->Fill(nPassingPrunedGP);

  // std::cout << "GV:   " << GVCollections.at(0).size() << std::endl;
  // std::cout << "GVs:  " << GVCollections.at(1).size() << std::endl;
  // std::cout << "GVn:  " << GVCollections.at(2).size() << std::endl;
  // std::cout << "GVns: " << GVCollections.at(3).size() << std::endl;
  // std::cout << "PGn:  " << GVCollections.at(4).size() << std::endl;
  // std::cout << "PGnB: " << GVCollections.at(5).size() << std::endl;
  // std::cout << "PGnD: " << GVCollections.at(6).size() << std::endl;
  // std::cout << "PGs:  " << GVCollections.at(7).size() << std::endl;
  // std::cout << "TV:   " << GVCollections.at(8).size() << std::endl;
  // std::cout << "TVn:  " << GVCollections.at(9).size() << std::endl;

  std::vector<SecondaryVertexCollection> SVCollections;
  SVCollections.push_back(svc_->getSecondaryVertexCollection());
  SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDBS());
  SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDBS4());
  SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDPV());
  for (unsigned int iColl = 0; iColl < SVCollections.size(); iColl++) {
    histos1_["n" + sv_names_.at(iColl)]->Fill(SVCollections.at(iColl).size());
    for (SecondaryVertex& sv : SVCollections.at(iColl)) sv.fill(histos1_, histos2_, sv_names_.at(iColl));
  }

  unsigned int nC = iEvent.get(nIVFClustersToken_);
  unsigned int nCBS = iEvent.get(nIVFClustersMTDBSToken_);
  unsigned int nCBS4 = iEvent.get(nIVFClustersMTDBS4Token_);
  unsigned int nCPV = iEvent.get(nIVFClustersMTDPVToken_);
  histos1_["nc"]->Fill(nC);
  histos1_["ncbs"]->Fill(nCBS);
  histos1_["ncbs4"]->Fill(nCBS4);
  histos1_["ncpv"]->Fill(nCPV);

  std::vector<RecoJetCollection> RJCollections;
  RJCollections.push_back(rjc_->getRecoJetCollection());
  RJCollections.push_back(rjc_->getRecoJetGenMatchCollection());
  for (unsigned int iColl = 0; iColl < RJCollections.size(); iColl++) {
    histos1_["n" + rj_names_.at(iColl)]->Fill(RJCollections.at(iColl).size());
    for (RecoJet& rj : RJCollections.at(iColl)) rj.fill(histos1_, rj_names_.at(iColl));
  }

  std::vector<GenJetCollection> GJCollections;
  GJCollections.push_back(gjc_->getGenJetCollection());
  GJCollections.push_back(gjc_->getGenJetRecoMatchCollection());
  for (unsigned int iColl = 0; iColl < GJCollections.size(); iColl++) {
    histos1_["n" + gj_names_.at(iColl)]->Fill(GJCollections.at(iColl).size());
    for (GenJet& gj : GJCollections.at(iColl)) gj.fill(histos1_, gj_names_.at(iColl));
  }

  // Matching GV and SV
  for (unsigned int iGVs = 0; iGVs < GVCollections.size(); iGVs++) {
    for (unsigned int iSVs = 0; iSVs < SVCollections.size(); iSVs++) {
      std::vector<bool> SVmatched(SVCollections.at(iSVs).size(), false);
      for (GenVertex& gv : GVCollections.at(iGVs)) {
        for (unsigned int iSV = 0; iSV < SVCollections.at(iSVs).size(); iSV++) {
          SecondaryVertex& sv = SVCollections.at(iSVs).at(iSV);
          if (matcher_->match(gv, sv, TRACK) && !SVmatched.at(iSV)) {
            SVmatched.at(iSV) = true;
            TString gv_name = gv_names_.at(iGVs) + "_" + sv_names_.at(iSVs);
            TString sv_name = sv_names_.at(iSVs) + "_" + gv_names_.at(iGVs);
            gv.fill(histos1_, gv_name);
            sv.fill(histos1_, histos2_, sv_name);
            matcher_->fill(histos1_, histos2_, gv_name, sv_name, gv, sv);
            break; // Restrict to only one match
          }
        }
      }
    }
  }
}


void VertexNtuplizer::beginJob() {}


void VertexNtuplizer::endJob() {

  if (scanCuts_) {
    std::cout << "Efficiencies:" << std::endl;
    std::cout << "trkMatchDrCut = " << trkMatchDrCut_ << std::endl;
    std::cout << "trkMatchPtCut = " << trkMatchPtCut_ << std::endl;
    std::cout << "gv_sv: " << histos1_["gv_sv_pt"]->GetEntries() / histos1_["gv_pt"]->GetEntries() << " "
        << histos1_["gv_sv_matchdxy"]->GetMean() << " " << histos1_["gv_sv_matchdxy"]->GetMeanError() << " "
        << histos1_["gv_sv_matchd3d"]->GetMean() << " " << histos1_["gv_sv_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gv_svbs4: " << histos1_["gv_svbs4_pt"]->GetEntries() / histos1_["gv_pt"]->GetEntries() << " "
        << histos1_["gv_svbs4_matchdxy"]->GetMean() << " " << histos1_["gv_svbs4_matchdxy"]->GetMeanError() << " "
        << histos1_["gv_svbs4_matchd3d"]->GetMean() << " " << histos1_["gv_svbs4_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gvs_sv: " << histos1_["gvs_sv_pt"]->GetEntries() / histos1_["gvs_pt"]->GetEntries() << " "
        << histos1_["gvs_sv_matchdxy"]->GetMean() << " " << histos1_["gvs_sv_matchdxy"]->GetMeanError() << " "
        << histos1_["gvs_sv_matchd3d"]->GetMean() << " " << histos1_["gvs_sv_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gvs_svbs4: " << histos1_["gvs_svbs4_pt"]->GetEntries() / histos1_["gvs_pt"]->GetEntries() << " "
        << histos1_["gvs_svbs4_matchdxy"]->GetMean() << " " << histos1_["gvs_svbs4_matchdxy"]->GetMeanError() << " "
        << histos1_["gvs_svbs4_matchd3d"]->GetMean() << " " << histos1_["gvs_svbs4_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gvn_sv: " << histos1_["gvn_sv_pt"]->GetEntries() / histos1_["gvn_pt"]->GetEntries() << " "
        << histos1_["gvn_sv_matchdxy"]->GetMean() << " " << histos1_["gvn_sv_matchdxy"]->GetMeanError() << " "
        << histos1_["gvn_sv_matchd3d"]->GetMean() << " " << histos1_["gvn_sv_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gvn_svbs4: " << histos1_["gvn_svbs4_pt"]->GetEntries() / histos1_["gvn_pt"]->GetEntries() << " "
        << histos1_["gvn_svbs4_matchdxy"]->GetMean() << " " << histos1_["gvn_svbs4_matchdxy"]->GetMeanError() << " "
        << histos1_["gvn_svbs4_matchd3d"]->GetMean() << " " << histos1_["gvn_svbs4_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gvns_sv: " << histos1_["gvns_sv_pt"]->GetEntries() / histos1_["gvns_pt"]->GetEntries() << " "
        << histos1_["gvns_sv_matchdxy"]->GetMean() << " " << histos1_["gvns_sv_matchdxy"]->GetMeanError() << " "
        << histos1_["gvns_sv_matchd3d"]->GetMean() << " " << histos1_["gvns_sv_matchd3d"]->GetMeanError() << std::endl;
    std::cout << "gvns_svbs4: " << histos1_["gvns_svbs4_pt"]->GetEntries() / histos1_["gvns_pt"]->GetEntries() << " "
        << histos1_["gvns_svbs4_matchdxy"]->GetMean() << " " << histos1_["gvns_svbs4_matchdxy"]->GetMeanError() << " "
        << histos1_["gvns_svbs4_matchd3d"]->GetMean() << " " << histos1_["gvns_svbs4_matchd3d"]->GetMeanError() << std::endl;
  } else {
    // Print summary
    std::cout << "Efficiencies:" << std::endl;
    std::cout << "trkMatchDrCut         = " << trkMatchDrCut_ << std::endl;
    std::cout << "trkMatchPtCut         = " << trkMatchPtCut_ << std::endl;
    std::cout << "GV-SV efficiency      = " << histos1_["gv_sv_pt"]->GetEntries() / histos1_["gv_pt"]->GetEntries() << std::endl;
    std::cout << "GV-SVbs efficiency    = " << histos1_["gv_svbs_pt"]->GetEntries() / histos1_["gv_pt"]->GetEntries() << std::endl;
    std::cout << "GV-SVbs4 efficiency   = " << histos1_["gv_svbs4_pt"]->GetEntries() / histos1_["gv_pt"]->GetEntries() << std::endl;
    std::cout << "GV-SVpv efficiency    = " << histos1_["gv_svpv_pt"]->GetEntries() / histos1_["gv_pt"]->GetEntries() << std::endl;
    std::cout << "GVs-SV efficiency     = " << histos1_["gvs_sv_pt"]->GetEntries() / histos1_["gvs_pt"]->GetEntries() << std::endl;
    std::cout << "GVs-SVbs efficiency   = " << histos1_["gvs_svbs_pt"]->GetEntries() / histos1_["gvs_pt"]->GetEntries() << std::endl;
    std::cout << "GVs-SVbs4 efficiency  = " << histos1_["gvs_svbs4_pt"]->GetEntries() / histos1_["gvs_pt"]->GetEntries() << std::endl;
    std::cout << "GVs-SVpv efficiency   = " << histos1_["gvs_svpv_pt"]->GetEntries() / histos1_["gvs_pt"]->GetEntries() << std::endl;
    std::cout << "GVn-SV efficiency     = " << histos1_["gvn_sv_pt"]->GetEntries() / histos1_["gvn_pt"]->GetEntries() << std::endl;
    std::cout << "GVn-SVbs efficiency   = " << histos1_["gvn_svbs_pt"]->GetEntries() / histos1_["gvn_pt"]->GetEntries() << std::endl;
    std::cout << "GVn-SVbs4 efficiency  = " << histos1_["gvn_svbs4_pt"]->GetEntries() / histos1_["gvn_pt"]->GetEntries() << std::endl;
    std::cout << "GVn-SVpv efficiency   = " << histos1_["gvn_svpv_pt"]->GetEntries() / histos1_["gvn_pt"]->GetEntries() << std::endl;
    std::cout << "GVns-SV efficiency    = " << histos1_["gvns_sv_pt"]->GetEntries() / histos1_["gvns_pt"]->GetEntries() << std::endl;
    std::cout << "GVns-SVbs efficiency  = " << histos1_["gvns_svbs_pt"]->GetEntries() / histos1_["gvns_pt"]->GetEntries() << std::endl;
    std::cout << "GVns-SVbs4 efficiency = " << histos1_["gvns_svbs4_pt"]->GetEntries() / histos1_["gvns_pt"]->GetEntries() << std::endl;
    std::cout << "GVns-SVpv efficiency  = " << histos1_["gvns_svpv_pt"]->GetEntries() / histos1_["gvns_pt"]->GetEntries() << std::endl;
  }

  // Catch under and over flows
  // for (auto iter : histos1_) {
  //   int nBins = iter.second->GetNbinsX();
  //   float firstBinContent = iter.second->GetBinContent(0) + iter.second->GetBinContent(1);
  //   float firstBinError = TMath::Sqrt(iter.second->GetBinError(0) * iter.second->GetBinError(0) +
  //       iter.second->GetBinError(1) * iter.second->GetBinError(1));
  //   float lastBinContent = iter.second->GetBinContent(nBins) + iter.second->GetBinContent(nBins + 1);
  //   float lastBinError = TMath::Sqrt(iter.second->GetBinError(nBins) * iter.second->GetBinError(nBins) +
  //       iter.second->GetBinError(nBins + 1) * iter.second->GetBinError(nBins + 1));
  //   iter.second->SetBinContent(1, firstBinContent);
  //   iter.second->SetBinError(1, firstBinError);
  //   iter.second->SetBinContent(nBins, lastBinContent);
  //   iter.second->SetBinError(nBins, lastBinError);
  // }

  for (auto iter : histos2_) {
    int nBinsX = iter.second->GetNbinsX();
    int nBinsY = iter.second->GetNbinsY();
    for (int iBinX = 1; iBinX <= nBinsX; iBinX++) {
      float firstBinContent = iter.second->GetBinContent(iBinX, 0) + iter.second->GetBinContent(iBinX, 1);
      float firstBinError = TMath::Sqrt(iter.second->GetBinError(iBinX, 0) * iter.second->GetBinError(iBinX, 0) +
          iter.second->GetBinError(iBinX, 1) * iter.second->GetBinError(iBinX, 1));
      float lastBinContent = iter.second->GetBinContent(iBinX, nBinsY) + iter.second->GetBinContent(iBinX, nBinsY + 1);
      float lastBinError = TMath::Sqrt(iter.second->GetBinError(iBinX, nBinsY) * iter.second->GetBinError(iBinX, nBinsY) +
          iter.second->GetBinError(iBinX, nBinsY + 1) * iter.second->GetBinError(iBinX, nBinsY + 1));
      iter.second->SetBinContent(iBinX, 1, firstBinContent);
      iter.second->SetBinError(iBinX, 1, firstBinError);
      iter.second->SetBinContent(iBinX, nBinsY, lastBinContent);
      iter.second->SetBinError(iBinX, nBinsY, lastBinError);
    }
    for (int iBinY = 1; iBinY <= nBinsY; iBinY++) {
      float firstBinContent = iter.second->GetBinContent(0, iBinY) + iter.second->GetBinContent(1, iBinY);
      float firstBinError = TMath::Sqrt(iter.second->GetBinError(0, iBinY) * iter.second->GetBinError(0, iBinY) +
          iter.second->GetBinError(1, iBinY) * iter.second->GetBinError(1, iBinY));
      float lastBinContent = iter.second->GetBinContent(nBinsX, iBinY) + iter.second->GetBinContent(nBinsX + 1, iBinY);
      float lastBinError = TMath::Sqrt(iter.second->GetBinError(nBinsX, iBinY) * iter.second->GetBinError(nBinsX, iBinY) +
          iter.second->GetBinError(nBinsX + 1, iBinY) * iter.second->GetBinError(nBinsX + 1, iBinY));
      iter.second->SetBinContent(1, iBinY, firstBinContent);
      iter.second->SetBinError(1, iBinY, firstBinError);
      iter.second->SetBinContent(nBinsX, iBinY, lastBinContent);
      iter.second->SetBinError(nBinsX, iBinY, lastBinError);
    }
  }
}


void VertexNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(VertexNtuplizer);
