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
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "../interface/GenVertexCollectionBuilder.h"
#include "../interface/GenVertex.h"
#include "../interface/SecondaryVertexCollectionBuilder.h"
#include "../interface/SecondaryVertex.h"
#include "../interface/VertexMatcher.h"
#include "../interface/RecoJetCollectionBuilder.h"
#include "../interface/RecoJet.h"
#include "../interface/GenJetCollectionBuilder.h"
#include "../interface/GenJet.h"


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
    edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeBSValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeBSErrorMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimeBSQualityMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimePVValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackTimePVErrorMapToken_;
    // edm::EDGetTokenT<edm::ValueMap<float>> trackTimePVQualityMapToken_;
    edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDBSToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDBS4Token_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDPVToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDBSToken_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDBS4Token_;
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDPVToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> slimmedCandSVToken_;
    edm::EDGetTokenT<pat::JetCollection> jetsToken_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken_;

    float nGVs_;
    float nGVsB_;
    float nGVsD_;
    float nGeneralTracks_;
    float nPFCandidates_;
    float nInclusiveSVs_;
    float nSlimmedCandSVs_;

    GenVertexCollectionBuilder* gvc_;
    SecondaryVertexCollectionBuilder* svc_;
    RecoJetCollectionBuilder* rjc_;
    GenJetCollectionBuilder* gjc_;
    VertexMatcher* matcher_;

    std::vector<TString> gv_names_;
    std::vector<TString> gv_names_more_;
    std::vector<TString> sv_names_;
    std::vector<TString> rj_names_;
    std::vector<TString> gj_names_;

    std::map<TString, std::vector<float>> vars1_;
    std::map<TString, std::vector<float>> vars2_;

    std::map<TString, TH1F*> histos1_;
    std::map<TString, TH2F*> histos2_;
};


VertexNtuplizer::VertexNtuplizer(const edm::ParameterSet& iConfig) :
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"))),
    simTracksToken_(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simTracks"))),
    generalTracksToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("generalTracks"))),
    pfCandidatesToken_(consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("particleFlowCandidates"))),
    trackTimeBSValueMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeBSValueMap"))),
    trackTimeBSErrorMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeBSErrorMap"))),
    trackTimeBSQualityMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimeBSQualityMap"))),
    trackTimePVValueMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimePVValueMap"))),
    trackTimePVErrorMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimePVErrorMap"))),
    // trackTimePVQualityMapToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackTimePVQualityMap"))),
    primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
    nIVFClustersToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClusters"))),
    nIVFClustersMTDBSToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDBS"))),
    nIVFClustersMTDBS4Token_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDBS4"))),
    nIVFClustersMTDPVToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDPV"))),
    secondaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices"))),
    secondaryVerticesMTDBSToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDBS"))),
    secondaryVerticesMTDBS4Token_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDBS4"))),
    secondaryVerticesMTDPVToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVerticesMTDPV"))),
    slimmedCandSVToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("slimmedCandSVs"))),
    jetsToken_(consumes<pat::JetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
    genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJets"))),
    genJetsFlavourInfoToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJetsFlavourInfo"))) {

  nGVs_ = 0.0;
  nGVsB_ = 0.0;
  nGVsD_ = 0.0;
  nGeneralTracks_ = 0.0;
  nPFCandidates_ = 0.0;
  nInclusiveSVs_ = 0.0;
  nSlimmedCandSVs_ = 0.0;

  gvc_ = new GenVertexCollectionBuilder(iConfig);
  svc_ = new SecondaryVertexCollectionBuilder(iConfig);
  rjc_ = new RecoJetCollectionBuilder(iConfig);
  gjc_ = new GenJetCollectionBuilder(iConfig);
  matcher_ = new VertexMatcher(iConfig);

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  gv_names_.push_back("gv");
  gv_names_.push_back("gvB");
  gv_names_.push_back("gvD");
  gv_names_.push_back("gvs"); // w/SIM match

  gv_names_more_.push_back("gvMatched");
  gv_names_more_.push_back("gvBMatched");
  gv_names_more_.push_back("gvDMatched");
  gv_names_more_.push_back("gvUnmatched");
  gv_names_more_.push_back("gvBUnmatched");
  gv_names_more_.push_back("gvDUnmatched");

  sv_names_.push_back("sv"); // From inclusiveVertexFinder
  // sv_names_.push_back("svbs"); // SecondaryVertex w/track time extrapolated to the beam spot
  // sv_names_.push_back("svbs4"); // SecondaryVertex w/range cut on track time extrapolated to the beam spot
  // sv_names_.push_back("svpv"); // SecondaryVertex w/track time extrapolated to the primary vertex
  sv_names_.push_back("svSlimmedCand"); // From inclusiveCandidateVertexFinder

  rj_names_.push_back("rj"); // RecoJet
  rj_names_.push_back("rjg"); // RecoJet w/GEN match

  gj_names_.push_back("gj"); // GenJet
  gj_names_.push_back("gjr"); // GenJet w/reco match

  std::vector<TString> objs_;
  for (TString obj : gv_names_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk");
    objs_.push_back(obj + "_trk_matchgt");
    objs_.push_back(obj + "_trk_matchpfc");
  }
  for (TString obj : gv_names_more_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk");
  }
  for (TString obj : sv_names_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk");
    objs_.push_back(obj + "_trk_matchgt");
    objs_.push_back(obj + "_trk_matchpfc");
  }
  for (TString obj : rj_names_) objs_.push_back(obj);
  for (TString obj : gj_names_) objs_.push_back(obj);

  const unsigned int nbins_ = 50;
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

  // Do we need this other than for comparison purposes?
  histos1_["nGPs"] = fs->make<TH1F>("nGPs", "nGPs", 13, 0, 13);
  // Count histograms
  for (TString gv_name : gv_names_) {
    TString name = "n" + gv_name;
    histos1_[name] = fs->make<TH1F>(name, name, ngv_, 0, ngv_);
    histos1_[name]->Sumw2();
  }
  for (TString gv_name : gv_names_more_) {
    TString name = "n" + gv_name;
    histos1_[name] = fs->make<TH1F>(name, name, ngv_, 0, ngv_);
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
  for (TString sv_name : sv_names_) {
    TString name = "n" + sv_name;
    histos1_[name] = fs->make<TH1F>(name, name, nsv_, 0, nsv_);
    histos1_[name]->Sumw2();
  }
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
      if (obj1.Contains("_trk") && obj2.Contains("_trk")) continue;
      if (obj1.Contains("_trk")) continue;
      if (obj1.Contains("_match") && obj2.Contains("_match")) continue;

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

  unsigned int nPassingGPs = gvc_->build(iEvent, genParticlesToken_, simTracksToken_, primaryVertex, matcher_);
  svc_->build(iEvent,
      secondaryVerticesToken_, secondaryVerticesMTDBSToken_, secondaryVerticesMTDBS4Token_, secondaryVerticesMTDPVToken_,
      trackTimeBSValueMapToken_, trackTimeBSErrorMapToken_, trackTimeBSQualityMapToken_,
      trackTimePVValueMapToken_, trackTimePVErrorMapToken_, // trackTimePVQualityMapToken_,
      slimmedCandSVToken_, primaryVertex);
  rjc_->build(iEvent, jetsToken_, genJetsFlavourInfoToken_);
  gjc_->build(iEvent, genJetsToken_, genJetsFlavourInfoToken_, jetsToken_);

  histos1_["nGPs"]->Fill(nPassingGPs);

  const reco::TrackCollection generalTracks = iEvent.get(generalTracksToken_);
  const reco::PFCandidateCollection pfCandidates = iEvent.get(pfCandidatesToken_);
  nGeneralTracks_ += generalTracks.size();
  nPFCandidates_ += pfCandidates.size();

  std::map<TString, std::vector<bool>*> matches_;
  for (TString gv_name : gv_names_) {
    matches_[gv_name + "_gt"] = new std::vector<bool>(generalTracks.size(), false);
    matches_[gv_name + "_pfc"] = new std::vector<bool>(pfCandidates.size(), false);
    for (TString sv_name : sv_names_) {
      matches_[gv_name + "_" + sv_name + "_gt"] = new std::vector<bool>(generalTracks.size(), false);
      matches_[gv_name + "_" + sv_name + "_pfc"] = new std::vector<bool>(pfCandidates.size(), false);
    }
  }
  for (TString gv_name : gv_names_more_) {
    matches_[gv_name + "_gt"] = new std::vector<bool>(generalTracks.size(), false);
    matches_[gv_name + "_pfc"] = new std::vector<bool>(pfCandidates.size(), false);
  }

  std::vector<GenVertexCollection> GVCollections;
  GVCollections.push_back(gvc_->getGenVertices());
  GVCollections.push_back(gvc_->getGenVerticesB());
  GVCollections.push_back(gvc_->getGenVerticesD());
  GVCollections.push_back(gvc_->getGenVerticesSimMatch());
  for (unsigned int iColl = 0; iColl < GVCollections.size(); iColl++) {
    histos1_["n" + gv_names_.at(iColl)]->Fill(GVCollections.at(iColl).size());
    for (GenVertex& gv : GVCollections.at(iColl)) {
      gv.fill(histos1_, gv_names_.at(iColl), matcher_, generalTracks, pfCandidates, matches_);
    }
  }
  nGVs_ += GVCollections.at(0).size();
  nGVsB_ += GVCollections.at(1).size();
  nGVsD_ += GVCollections.at(2).size();

  unsigned int nC = iEvent.get(nIVFClustersToken_);
  unsigned int nCBS = iEvent.get(nIVFClustersMTDBSToken_);
  unsigned int nCBS4 = iEvent.get(nIVFClustersMTDBS4Token_);
  unsigned int nCPV = iEvent.get(nIVFClustersMTDPVToken_);
  histos1_["nc"]->Fill(nC);
  histos1_["ncbs"]->Fill(nCBS);
  histos1_["ncbs4"]->Fill(nCBS4);
  histos1_["ncpv"]->Fill(nCPV);

  std::vector<SecondaryVertexCollection> SVCollections;
  SVCollections.push_back(svc_->getSecondaryVertexCollection());
  // SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDBS());
  // SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDBS4());
  // SVCollections.push_back(svc_->getSecondaryVertexCollectionMTDPV());
  SVCollections.push_back(svc_->getSlimmedCandVertexCollection());
  for (unsigned int iColl = 0; iColl < SVCollections.size(); iColl++) {
    histos1_["n" + sv_names_.at(iColl)]->Fill(SVCollections.at(iColl).size());
    for (SecondaryVertex& sv : SVCollections.at(iColl)) sv.fill(histos1_, histos2_, sv_names_.at(iColl));
  }
  nInclusiveSVs_ += SVCollections.at(0).size();
  nSlimmedCandSVs_ += SVCollections.at(1).size();

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

  // Only for matching to inclusive SVs
  // TODO: expand it to include matching to all SV collections
  GenVertexCollection matchedGVs;
  GenVertexCollection matchedGVsB;
  GenVertexCollection matchedGVsD;
  GenVertexCollection unmatchedGVs;
  GenVertexCollection unmatchedGVsB;
  GenVertexCollection unmatchedGVsD;

  // Matching GV and SV
  for (unsigned int iGVs = 0; iGVs < GVCollections.size(); iGVs++) {
    for (unsigned int iSVs = 0; iSVs < SVCollections.size(); iSVs++) {
      std::vector<bool> SVmatched(SVCollections.at(iSVs).size(), false); // Prevent double matching
      for (GenVertex& gv : GVCollections.at(iGVs)) {
        bool GVmatched = false;
        for (unsigned int iSV = 0; iSV < SVCollections.at(iSVs).size(); iSV++) {
          SecondaryVertex& sv = SVCollections.at(iSVs).at(iSV);

          if (matcher_->match(gv, sv, VertexMatcher::TRACK) && !SVmatched.at(iSV)) {
            SVmatched.at(iSV) = true; // Prevent double matching
            GVmatched = true;
            TString gv_name = gv_names_.at(iGVs) + "_" + sv_names_.at(iSVs);
            TString sv_name = sv_names_.at(iSVs) + "_" + gv_names_.at(iGVs);
            gv.fill(histos1_, gv_name, matcher_, generalTracks, pfCandidates, matches_);
            sv.fill(histos1_, histos2_, sv_name);
            matcher_->fill(histos1_, histos2_, gv_name, sv_name, gv, sv);
            break; // Restrict to only one match
          } // End if matched
        } // End loop through SVs

        if (sv_names_.at(iSVs).EqualTo("sv") && gv_names_.at(iGVs).EqualTo("gv")) {
          if (GVmatched) {
            matchedGVs.push_back(gv);
            if (gv.pdgIdBin() == B_MESON || gv.pdgIdBin() == B_BARYON) matchedGVsB.push_back(gv);
            if (gv.pdgIdBin() == C_MESON || gv.pdgIdBin() == C_BARYON) matchedGVsD.push_back(gv);
          } else {
            unmatchedGVs.push_back(gv);
            if (gv.pdgIdBin() == B_MESON || gv.pdgIdBin() == B_BARYON) unmatchedGVsB.push_back(gv);
            if (gv.pdgIdBin() == C_MESON || gv.pdgIdBin() == C_BARYON) unmatchedGVsD.push_back(gv);
          } // End if GV is matched
        } // End if SV collection is inclusive SVs

      } // End loop through GVs
    } // End loop through SV collections
  } // End loop through GV collections

  std::vector<GenVertexCollection> GVCollectionsMore;
  GVCollectionsMore.push_back(matchedGVs);
  GVCollectionsMore.push_back(matchedGVsB);
  GVCollectionsMore.push_back(matchedGVsD);
  GVCollectionsMore.push_back(unmatchedGVs);
  GVCollectionsMore.push_back(unmatchedGVsB);
  GVCollectionsMore.push_back(unmatchedGVsD);
  const reco::TrackCollection emptyGTs;
  const reco::PFCandidateCollection emptyPFCs;
  for (unsigned int iColl = 0; iColl < GVCollectionsMore.size(); iColl++) {
    histos1_["n" + gv_names_more_.at(iColl)]->Fill(GVCollectionsMore.at(iColl).size());
    for (GenVertex& gv : GVCollectionsMore.at(iColl)) {
      gv.fill(histos1_, gv_names_more_.at(iColl), matcher_, emptyGTs, emptyPFCs, matches_);
    }
  }
}


void VertexNtuplizer::beginJob() {}


void VertexNtuplizer::endJob() {

  std::cout << "\nnumber of events = " << histos1_["nGPs"]->GetEntries() << std::endl;

  std::cout << "\nthere are " << nGeneralTracks_ / histos1_["nGPs"]->GetEntries() << " general tracks per event." << std::endl;
  std::cout << "there are " << nPFCandidates_ / histos1_["nGPs"]->GetEntries() << " PF candidates per event." << std::endl;

  std::cout << "mean number of GVs from               = " << histos1_["ngv"]->GetMean() << "+=" << histos1_["ngv"]->GetStdDev() << std::endl;
  std::cout << "mean number of GVs from with B mother = " << histos1_["ngvB"]->GetMean() << "+=" << histos1_["ngvB"]->GetStdDev() << std::endl;
  std::cout << "mean number of GVs from with D mother = " << histos1_["ngvD"]->GetMean() << "+=" << histos1_["ngvD"]->GetStdDev() << std::endl;

  std::cout << "\nmean number of inclusive SVs = " << histos1_["nsv"]->GetMean() << "+=" << histos1_["nsv"]->GetStdDev() << std::endl;
  std::cout << "mean number of slimmed SVs   = " << histos1_["nsvSlimmedCand"]->GetMean() << "+=" << histos1_["nsvSlimmedCand"]->GetStdDev() << std::endl;
  std::cout << "slimmed SVs / inclusive SVs  = " << nSlimmedCandSVs_ / nInclusiveSVs_ << std::endl;

  std::cout << "\nsv" << std::endl;
  std::cout << "Integrated efficiency for all GVs                = " << histos1_["gv_sv_pdgIdBin"]->GetEntries() / histos1_["gv_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs with B mother  = " << histos1_["gvB_sv_pdgIdBin"]->GetEntries() / histos1_["gvB_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs with D mother  = " << histos1_["gvD_sv_pdgIdBin"]->GetEntries() / histos1_["gvD_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs SimTrack match = " << histos1_["gvs_sv_pdgIdBin"]->GetEntries() / histos1_["gvs_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "\nsvSlimmedCand" << std::endl;
  std::cout << "Integrated efficiency for all GVs                = " << histos1_["gv_svSlimmedCand_pdgIdBin"]->GetEntries() / histos1_["gv_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs with B mother  = " << histos1_["gvB_svSlimmedCand_pdgIdBin"]->GetEntries() / histos1_["gvB_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs with D mother  = " << histos1_["gvD_svSlimmedCand_pdgIdBin"]->GetEntries() / histos1_["gvD_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs SimTrack match = " << histos1_["gvs_svSlimmedCand_pdgIdBin"]->GetEntries() / histos1_["gvs_pdgIdBin"]->GetEntries() << std::endl;

  // Catch under and over flows -- messes with mean and stddev calculations...
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
