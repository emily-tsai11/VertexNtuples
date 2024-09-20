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
#include <algorithm>

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

    void book1DHistograms(edm::Service<TFileService>& fs,
        std::map<TString, std::vector<float>>& vars,
        const TString& obj, std::map<TString, TH1F*>& histos);
    void book2DHistograms(edm::Service<TFileService>& fs,
        std::map<TString, std::vector<float>>& vars,
        const TString& obj, std::map<TString, TH2F*>& histos);

    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
    edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackT0FromBSToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackSigmaT0FromBSToken_;
    // edm::EDGetTokenT<edm::ValueMap<float>> trackQualityFromBSToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackT0FromPVToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> trackSigmaT0FromPVToken_;
    // edm::EDGetTokenT<edm::ValueMap<float>> trackQualityFromPVToken_;
    edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersToken_;
    edm::EDGetTokenT<unsigned int> nIVFClustersMTDPVToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> inclusiveSecondaryVerticesToken_;
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> inclusiveSecondaryVerticesMTDPVToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> mergedSecondaryVerticesToken_;
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> mergedSecondaryVerticesMTDPVToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> slimmedSecondaryVerticesToken_;
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> slimmedSecondaryVerticesMTDPVToken_;
    edm::EDGetTokenT<pat::JetCollection> jetsToken_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavourInfoToken_;

    float trkPtMin_;
    float jetPtMin_;
    float jetPtMax_;
    float absEtaMax_;

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

    std::map<TString, std::vector<float>> vars1_trks_;
    std::map<TString, std::vector<float>> vars1_trks_sv_;
    std::map<TString, std::vector<float>> vars1_vtxs_gv_;
    std::map<TString, std::vector<float>> vars1_vtxs_sv_;
    std::map<TString, std::vector<float>> vars1_jets_;
    std::map<TString, std::vector<float>> vars1_trkVtxMatch_;

    std::map<TString, std::vector<float>> vars2_trks_;
    std::map<TString, std::vector<float>> vars2_trks_sv_;
    std::map<TString, std::vector<float>> vars2_vtxs_;
    std::map<TString, std::vector<float>> vars2_vtxs_sv_;
    // std::map<TString, std::vector<float>> vars2_jets_;
    std::map<TString, std::vector<float>> vars2_trkVtxMatch_;

    std::map<TString, TH1F*> histos1_;
    std::map<TString, TH2F*> histos2_;
};


VertexNtuplizer::VertexNtuplizer(const edm::ParameterSet& iConfig) :
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"))),
    simTracksToken_(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simTracks"))),
    generalTracksToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("generalTracks"))),
    pfCandidatesToken_(consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidates"))),
    trackT0FromBSToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackT0FromBS"))),
    trackSigmaT0FromBSToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackSigmaT0FromBS"))),
    // trackQualityFromBSToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackQualityFromBS"))),
    trackT0FromPVToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackT0FromPV"))),
    trackSigmaT0FromPVToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackSigmaT0FromPV"))),
    // trackQualityFromPVToken_(consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("trackQualityFromPV"))),
    primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
    nIVFClustersToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClusters"))),
    nIVFClustersMTDPVToken_(consumes<unsigned int>(iConfig.getUntrackedParameter<edm::InputTag>("nIVFClustersMTDPV"))),
    inclusiveSecondaryVerticesToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("inclusiveSecondaryVertices"))),
    // inclusiveSecondaryVerticesMTDPVToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("inclusiveSecondaryVerticesMTDPV"))),
    mergedSecondaryVerticesToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("mergedSecondaryVertices"))),
    // mergedSecondaryVerticesMTDPVToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("mergedSecondaryVerticesMTDPV"))),
    slimmedSecondaryVerticesToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("slimmedSecondaryVertices"))),
    // slimmedSecondaryVerticesMTDPVToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("slimmedSecondaryVerticesMTDPV"))),
    jetsToken_(consumes<pat::JetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
    genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJets"))),
    genJetsFlavourInfoToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genJetsFlavourInfo"))) {

  trkPtMin_ = (float) iConfig.getUntrackedParameter<double>("trkPtMin");
  jetPtMin_ = (float) iConfig.getUntrackedParameter<double>("jetPtMin");
  jetPtMax_ = (float) iConfig.getUntrackedParameter<double>("jetPtMax");
  absEtaMax_ = (float) iConfig.getUntrackedParameter<double>("absEtaMax");

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
  // gv_names_.push_back("gvs"); // w/SIM match
  // gv_names_.push_back("gvsB"); // w/SIM match
  // gv_names_.push_back("gvsD"); // w/SIM match

  sv_names_.push_back("svInclusive");
  sv_names_.push_back("svInclusiveMTDBS");
  sv_names_.push_back("svInclusiveMTDPV");
  sv_names_.push_back("svSlimmed");
  sv_names_.push_back("svSlimmedMTDBS");
  sv_names_.push_back("svSlimmedMTDPV");

  rj_names_.push_back("rj"); // RecoJet
  rj_names_.push_back("rjg"); // RecoJet w/GEN match

  gj_names_.push_back("gj"); // GenJet
  gj_names_.push_back("gjr"); // GenJet w/reco match

  std::vector<TString> objs_;
  for (TString obj : gv_names_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk"); // daughters
  }
  for (TString obj : sv_names_) {
    objs_.push_back(obj);
    objs_.push_back(obj + "_trk");
  }
  for (TString obj : rj_names_) objs_.push_back(obj);
  for (TString obj : gj_names_) objs_.push_back(obj);

  const unsigned int nbins_ = 50;
  const unsigned int ngv_ = 15;
  const unsigned int nclus_ = 100;
  const unsigned int nsv_ = 60;
  const unsigned int njet_ = 20;

  const float trkPtMax = 50.0;
  const float vtxPtMin = trkPtMin_;
  const float vtxPtMax = 150.0;

  // ----------------------------------------------------------------------- //
  // 1D histograms                                                           //
  // ----------------------------------------------------------------------- //

  // -------- Tracks only -------- //
  std::map<TString, std::vector<float>> vars1_trk_;
  vars1_trk_["pt"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax};
  vars1_trk_["pt2"] = std::vector<float>{(float) nbins_, trkPtMin_*trkPtMin_, trkPtMax*trkPtMax};
  std::map<TString, std::vector<float>> vars1_trk_sv_;
  vars1_trk_sv_["tvalBTL"] = std::vector<float>{(float) nbins_, -800.0, 1000.0}; // in ps
  vars1_trk_sv_["terrBTL"] = std::vector<float>{(float) nbins_, 0.0, 80.0}; // in ps
  vars1_trk_sv_["tsigBTL"] = std::vector<float>{(float) nbins_, -20.0, 20.0};
  // vars1_trk_sv_["tqualBTL"] = std::vector<float>{(float) nbins_, 0.0, 1.0};
  vars1_trk_sv_["tvalETL"] = std::vector<float>{(float) nbins_, -800.0, 1000.0}; // in ps
  vars1_trk_sv_["terrETL"] = std::vector<float>{(float) nbins_, 0.0, 80.0}; // in ps
  vars1_trk_sv_["tsigETL"] = std::vector<float>{(float) nbins_, -20.0, 20.0};
  // vars1_trk_sv_["tqualETL"] = std::vector<float>{(float) nbins_, 0.0, 1.0};

  // -------- Vertices only -------- //
  std::map<TString, std::vector<float>> vars1_vtx_;
  vars1_vtx_["pt"] = std::vector<float>{(float) nbins_, vtxPtMin, vtxPtMax};
  vars1_vtx_["pt2"] = std::vector<float>{(float) nbins_, vtxPtMin*vtxPtMin, vtxPtMax*vtxPtMax};
  vars1_vtx_["ntrk"] = std::vector<float>{10.0, 0.0, 10.0}; // daughters for GenVertex
  std::map<TString, std::vector<float>> vars1_vtx_gv_;
  vars1_vtx_gv_["nImmediateDau"] = std::vector<float>{10.0, 0.0, 10.0};
  vars1_vtx_gv_["nFinalDau"] = std::vector<float>{10.0, 0.0, 10.0};
  std::map<TString, std::vector<float>> vars1_vtx_sv_;
  vars1_vtx_sv_["tavg"] = std::vector<float>{(float) nbins_, -800.0, 1000.0}; // in ps
  vars1_vtx_sv_["trange"] = std::vector<float>{(float) nbins_, 0.0, 1600.0}; // in ps
  vars1_vtx_sv_["mtdEff"] = std::vector<float>{(float) nbins_, 0.0, 1.0};

  // -------- Jets only -------- //
  std::map<TString, std::vector<float>> vars1_jet_;
  vars1_jet_["pt"] = std::vector<float>{(float) nbins_, jetPtMin_, jetPtMax_};
  vars1_jet_["pt2"] = std::vector<float>{(float) nbins_, jetPtMin_*jetPtMin_, jetPtMax_*jetPtMax_};
  vars1_jet_["hadFlav"] = std::vector<float>{6.0, 0.0, 6.0};

  // -------- Tracks & vertices -------- //
  std::map<TString, std::vector<float>> vars1_trkVtx_;
  vars1_trkVtx_["charge"] = std::vector<float>{3.0, -1.0, 2.0};
  vars1_trkVtx_["motherPdgId"] = std::vector<float>{(float) nbins_, -5560.0, 5560.0};
  vars1_trkVtx_["pdgIdBin"] = std::vector<float>{5.0, 0.0, 5.0};
  // From origin (in cm)
  vars1_trkVtx_["x"] = std::vector<float>{(float) nbins_, -1.0, 1.0};
  vars1_trkVtx_["y"] = std::vector<float>{(float) nbins_, -1.0, 1.0};
  vars1_trkVtx_["z"] = std::vector<float>{(float) nbins_, -20.0, 20.0};
  vars1_trkVtx_["xerr"] = std::vector<float>{(float) nbins_, 0.0, 0.5};
  vars1_trkVtx_["yerr"] = std::vector<float>{(float) nbins_, 0.0, 0.5};
  vars1_trkVtx_["zerr"] = std::vector<float>{(float) nbins_, 0.0, 0.5};
  // From primary vertex (in cm)
  vars1_trkVtx_["dxy"] = std::vector<float>{(float) nbins_, 0.0, 4.0};
  vars1_trkVtx_["dz"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_trkVtx_["d3d"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_trkVtx_["dxyerr"] = std::vector<float>{(float) nbins_, 0.0, 0.0025};
  vars1_trkVtx_["dzerr"] = std::vector<float>{(float) nbins_, 0.0, 0.005};
  vars1_trkVtx_["d3derr"] = std::vector<float>{(float) nbins_, 0.0, 0.005};
  vars1_trkVtx_["dxysig"] = std::vector<float>{(float) nbins_, 0.0, 80.0};
  vars1_trkVtx_["dzsig"] = std::vector<float>{(float) nbins_, 0.0, 200.0};
  vars1_trkVtx_["d3dsig"] = std::vector<float>{(float) nbins_, 0.0, 200.0};
  std::map<TString, std::vector<float>> vars1_trkVtx_sv_;
  vars1_trkVtx_sv_["chi2"] = std::vector<float>{(float) nbins_, 0.0, 100.0};
  vars1_trkVtx_sv_["ndof"] = std::vector<float>{(float) nbins_, 0.0, 10.0};
  vars1_trkVtx_sv_["chi2dof"] = std::vector<float>{(float) nbins_, 0.0, 10.0};

  // -------- Tracks & vertices & jets -------- //
  std::map<TString, std::vector<float>> vars1_trkVtxJet_;
  vars1_trkVtxJet_["eta"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_};
  vars1_trkVtxJet_["phi"] = std::vector<float>{(float) nbins_, -3.15, 3.15};

  // -------- Matched tracks & vertices -------- //
  // Between matched GenVertex/SecondaryVertex or daughter/pfcandidate (in cm)
  vars1_trkVtxMatch_["xres"] = std::vector<float>{(float) nbins_, -0.12, 0.12};
  vars1_trkVtxMatch_["yres"] = std::vector<float>{(float) nbins_, -0.12, 0.12};
  vars1_trkVtxMatch_["zres"] = std::vector<float>{(float) nbins_, -0.12, 0.12};
  vars1_trkVtxMatch_["xpull"] = std::vector<float>{(float) nbins_, -1, 1};
  vars1_trkVtxMatch_["ypull"] = std::vector<float>{(float) nbins_, -1, 1};
  vars1_trkVtxMatch_["zpull"] = std::vector<float>{(float) nbins_, -1, 1};
  vars1_trkVtxMatch_["matchdxy"] = std::vector<float>{(float) nbins_, 0.0, 4.0};
  vars1_trkVtxMatch_["matchd3d"] = std::vector<float>{(float) nbins_, 0.0, 20.0};
  vars1_trkVtxMatch_["matchdxyerr"] = std::vector<float>{(float) nbins_, 0.0, 0.05};
  vars1_trkVtxMatch_["matchd3derr"] = std::vector<float>{(float) nbins_, 0.0, 0.1};
  vars1_trkVtxMatch_["matchdxysig"] = std::vector<float>{(float) nbins_, 0.0, 80.0};
  vars1_trkVtxMatch_["matchd3dsig"] = std::vector<float>{(float) nbins_, 0.0, 200.0};
  vars1_trkVtxMatch_["deltaR"] = std::vector<float>{(float) nbins_, 0, 4.0};
  vars1_trkVtxMatch_["normPtRes"] = std::vector<float>{(float) nbins_, 0.0, 1.0};

  vars1_trks_.insert(vars1_trk_.begin(), vars1_trk_.end());
  vars1_trks_.insert(vars1_trkVtx_.begin(), vars1_trkVtx_.end());
  vars1_trks_.insert(vars1_trkVtxJet_.begin(), vars1_trkVtxJet_.end());

  vars1_trks_sv_.insert(vars1_trk_.begin(), vars1_trk_.end());
  vars1_trks_sv_.insert(vars1_trk_sv_.begin(), vars1_trk_sv_.end());
  vars1_trks_sv_.insert(vars1_trkVtx_.begin(), vars1_trkVtx_.end());
  vars1_trks_sv_.insert(vars1_trkVtx_sv_.begin(), vars1_trkVtx_sv_.end());
  vars1_trks_sv_.insert(vars1_trkVtxJet_.begin(), vars1_trkVtxJet_.end());

  vars1_vtxs_gv_.insert(vars1_vtx_.begin(), vars1_vtx_.end());
  vars1_vtxs_gv_.insert(vars1_vtx_gv_.begin(), vars1_vtx_gv_.end());
  vars1_vtxs_gv_.insert(vars1_trkVtx_.begin(), vars1_trkVtx_.end());
  vars1_vtxs_gv_.insert(vars1_trkVtxJet_.begin(), vars1_trkVtxJet_.end());

  vars1_vtxs_sv_.insert(vars1_vtx_.begin(), vars1_vtx_.end());
  vars1_vtxs_sv_.insert(vars1_vtx_sv_.begin(), vars1_vtx_sv_.end());
  vars1_vtxs_sv_.insert(vars1_trkVtx_.begin(), vars1_trkVtx_.end());
  vars1_vtxs_sv_.insert(vars1_trkVtx_sv_.begin(), vars1_trkVtx_sv_.end());
  vars1_vtxs_sv_.insert(vars1_trkVtxJet_.begin(), vars1_trkVtxJet_.end());

  vars1_jets_.insert(vars1_jet_.begin(), vars1_jet_.end());
  vars1_jets_.insert(vars1_trkVtxJet_.begin(), vars1_trkVtxJet_.end());

  // ----------------------------------------------------------------------- //
  // 2D histograms                                                           //
  // ----------------------------------------------------------------------- //

  // -------- Tracks only -------- //
  vars2_trks_["eta_tval"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_, (float) nbins_, -800.0, 1000.0}; // in ps
  vars2_trks_["eta_terr"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_, (float) nbins_, 0.0, 80.0}; // in ps
  vars2_trks_["eta_tsig"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_, (float) nbins_, -20.0, 20.0};
  // vars2_trks_["eta_tqual"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_, (float) nbins_, 0.0, 1.0};
  vars2_trks_["pt_tvalBTL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, -800.0, 1000.0}; // in ps
  vars2_trks_["pt_terrBTL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, 0.0, 80.0}; // in ps
  vars2_trks_["pt_tsigBTL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, -20.0, 20.0};
  // vars2_trks_["pt_tqualBTL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, 0.0, 1.0};
  vars2_trks_["pt_tvalETL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, -800.0, 1000.0}; // in ps
  vars2_trks_["pt_terrETL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, 0.0, 80.0}; // in ps
  vars2_trks_["pt_tsigETL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, -20.0, 20.0};
  // vars2_trks_["pt_tqualETL"] = std::vector<float>{(float) nbins_, trkPtMin_, trkPtMax, (float) nbins_, 0.0, 1.0};

  // -------- Vertices only -------- //
  vars2_vtxs_["eta_tavg"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_, (float) nbins_, -800.0, 1000.0}; // in ps
  vars2_vtxs_["eta_trange"] = std::vector<float>{(float) nbins_, -absEtaMax_, absEtaMax_, (float) nbins_, 0.0, 1600.0}; // in ps
  vars2_vtxs_["trangeBTL_pt"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, vtxPtMin, vtxPtMax};
  vars2_vtxs_["trangeBTL_pt2"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, vtxPtMin*vtxPtMin, vtxPtMax*vtxPtMax};
  vars2_vtxs_["trangeETL_pt"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, vtxPtMin, vtxPtMax};
  vars2_vtxs_["trangeETL_pt2"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, vtxPtMin*vtxPtMin, vtxPtMax*vtxPtMax};
  // From primary vertex (in cm)
  vars2_vtxs_["trangeBTL_dxy"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 4.0};
  vars2_vtxs_["trangeBTL_dxysig"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 80.0};
  vars2_vtxs_["trangeBTL_d3d"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 20.0};
  vars2_vtxs_["trangeBTL_d3dsig"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 200.0};
  vars2_vtxs_["trangeETL_dxy"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 4.0};
  vars2_vtxs_["trangeETL_dxysig"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 80.0};
  vars2_vtxs_["trangeETL_d3d"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 20.0};
  vars2_vtxs_["trangeETL_d3dsig"] = std::vector<float>{(float) nbins_, 0.0, 1600.0, (float) nbins_, 0.0, 200.0};

  // -------- Tracks & vertices -------- //
  // Between matched GenVertex/SecondaryVertex or daughter/pfcandidate (in cm)
  vars2_trkVtxMatch_["matchdxy_dxy"] = std::vector<float>{(float) nbins_, 0.0, 4.0, (float) nbins_, 0.0, 4.0};
  vars2_trkVtxMatch_["matchdxy_d3d"] = std::vector<float>{(float) nbins_, 0.0, 4.0, (float) nbins_, 0.0, 10.0};
  vars2_trkVtxMatch_["matchd3d_dxy"] = std::vector<float>{(float) nbins_, 0.0, 10.0, (float) nbins_, 0.0, 4.0};
  vars2_trkVtxMatch_["matchd3d_d3d"] = std::vector<float>{(float) nbins_, 0.0, 10.0, (float) nbins_, 0.0, 10.0};
  vars2_trkVtxMatch_["deltaR_normPtRes"] = std::vector<float>{(float) nbins_, 0.0, 4.0, (float) nbins_, 0.0, 1.0};

  // ----------------------------------------------------------------------- //

  // Standalone histograms
  histos1_["nGPs"] = fs->make<TH1F>("nGPs", "nGPs", 13, 0, 13);
  histos1_["nGPsB"] = fs->make<TH1F>("nGPsB", "nGPsB", 13, 0, 13);
  histos1_["nGPsD"] = fs->make<TH1F>("nGPsD", "nGPsD", 13, 0, 13);
  histos1_["nGPs"]->Sumw2();
  histos1_["nGPsB"]->Sumw2();
  histos1_["nGPsD"]->Sumw2();

  histos1_["fracReconstructableGVs"] = fs->make<TH1F>("fracReconstructableGVs", "fracReconstructableGVs", nbins_ + 1, 0.0, 1.0 + 1.0 / ((float) nbins_));
  histos1_["fracReconstructableGVsB"] = fs->make<TH1F>("fracReconstructableGVsB", "fracReconstructableGVsB", nbins_ + 1, 0.0, 1.0 + 1.0 / ((float) nbins_));
  histos1_["fracReconstructableGVsD"] = fs->make<TH1F>("fracReconstructableGVsD", "fracReconstructableGVsD", nbins_ + 1, 0.0, 1.0 + 1.0 / ((float) nbins_));
  histos1_["fracReconstructableGVs"]->Sumw2();
  histos1_["fracReconstructableGVsB"]->Sumw2();
  histos1_["fracReconstructableGVsD"]->Sumw2();

  histos1_["maxMatchEffInclusive"] = fs->make<TH1F>("maxMatchEffInclusive", "maxMatchEffInclusive", nbins_ + 1, 0.0, 1.0 + 1.0 / ((float) nbins_));
  histos1_["maxMatchEffSlimmed"] = fs->make<TH1F>("maxMatchEffSlimmed", "maxMatchEffSlimmed", nbins_ + 1, 0.0, 1.0 + 1.0 / ((float) nbins_));
  histos1_["maxMatchedWithInclusive"] = fs->make<TH1F>("maxMatchedWithInclusive", "maxMatchedWithInclusive", ngv_, 0, ngv_);
  histos1_["maxMatchedWithSlimmed"] = fs->make<TH1F>("maxMatchedWithSlimmed", "maxMatchedWithSlimmed", ngv_, 0, ngv_);
  histos1_["maxMatchEffInclusive"]->Sumw2();
  histos1_["maxMatchEffSlimmed"]->Sumw2();
  histos1_["maxMatchedWithInclusive"]->Sumw2();
  histos1_["maxMatchedWithSlimmed"]->Sumw2();

  histos1_["pfCandidate_bs_tvalBTL"] = fs->make<TH1F>("pfCandidate_bs_tvalBTL", "pfCandidate_bs_tvalBTL", nbins_, -800.0, 1000.0); // in ps
  histos1_["pfCandidate_bs_terrBTL"] = fs->make<TH1F>("pfCandidate_bs_terrBTL", "pfCandidate_bs_terrBTL", nbins_, 0.0, 80.0); // in ps
  histos1_["pfCandidate_bs_tsigBTL"] = fs->make<TH1F>("pfCandidate_bs_tsigBTL", "pfCandidate_bs_tsigBTL", nbins_, -20.0, 20.0);
  // histos1_["pfCandidate_bs_tqualBTL"] = fs->make<TH1F>("pfCandidate_bs_tqualBTL", "pfCandidate_bs_tqualBTL", nbins_, 0.0, 1.0);
  histos1_["pfCandidate_bs_tvalETL"] = fs->make<TH1F>("pfCandidate_bs_tvalETL", "pfCandidate_bs_tvalETL", nbins_, -800.0, 1000.0); // in ps
  histos1_["pfCandidate_bs_terrETL"] = fs->make<TH1F>("pfCandidate_bs_terrETL", "pfCandidate_bs_terrETL", nbins_, 0.0, 80.0); // in ps
  histos1_["pfCandidate_bs_tsigETL"] = fs->make<TH1F>("pfCandidate_bs_tsigETL", "pfCandidate_bs_tsigETL", nbins_, -20.0, 20.0);
  // histos1_["pfCandidate_bs_tqualETL"] = fs->make<TH1F>("pfCandidate_bs_tqualETL", "pfCandidate_bs_tqualETL", nbins_, 0.0, 1.0);
  histos1_["pfCandidate_pv_tvalBTL"] = fs->make<TH1F>("pfCandidate_pv_tvalBTL", "pfCandidate_pv_tvalBTL", nbins_, -800.0, 1000.0); // in ps
  histos1_["pfCandidate_pv_terrBTL"] = fs->make<TH1F>("pfCandidate_pv_terrBTL", "pfCandidate_pv_terrBTL", nbins_, 0.0, 80.0); // in ps
  histos1_["pfCandidate_pv_tsigBTL"] = fs->make<TH1F>("pfCandidate_pv_tsigBTL", "pfCandidate_pv_tsigBTL", nbins_, -20.0, 20.0);
  // histos1_["pfCandidate_pv_tqualBTL"] = fs->make<TH1F>("pfCandidate_pv_tqualBTL", "pfCandidate_pv_tqualBTL", nbins_, 0.0, 1.0);
  histos1_["pfCandidate_pv_tvalETL"] = fs->make<TH1F>("pfCandidate_pv_tvalETL", "pfCandidate_pv_tvalETL", nbins_, -800.0, 1000.0); // in ps
  histos1_["pfCandidate_pv_terrETL"] = fs->make<TH1F>("pfCandidate_pv_terrETL", "pfCandidate_pv_terrETL", nbins_, 0.0, 80.0); // in ps
  histos1_["pfCandidate_pv_tsigETL"] = fs->make<TH1F>("pfCandidate_pv_tsigETL", "pfCandidate_pv_tsigETL", nbins_, -20.0, 20.0);
  // histos1_["pfCandidate_pv_tqualETL"] = fs->make<TH1F>("pfCandidate_pv_tqualETL", "pfCandidate_pv_tqualETL", nbins_, 0.0, 1.0);
  histos1_["pfCandidate_bs_tvalBTL"]->Sumw2();
  histos1_["pfCandidate_bs_terrBTL"]->Sumw2();
  histos1_["pfCandidate_bs_tsigBTL"]->Sumw2();
  // histos1_["pfCandidate_bs_tqualBTL"]->Sumw2();
  histos1_["pfCandidate_bs_tvalETL"]->Sumw2();
  histos1_["pfCandidate_bs_terrETL"]->Sumw2();
  histos1_["pfCandidate_bs_tsigETL"]->Sumw2();
  // histos1_["pfCandidate_bs_tqualETL"]->Sumw2();
  histos1_["pfCandidate_pv_tvalBTL"]->Sumw2();
  histos1_["pfCandidate_pv_terrBTL"]->Sumw2();
  histos1_["pfCandidate_pv_tsigBTL"]->Sumw2();
  // histos1_["pfCandidate_pv_tqualBTL"]->Sumw2();
  histos1_["pfCandidate_pv_tvalETL"]->Sumw2();
  histos1_["pfCandidate_pv_terrETL"]->Sumw2();
  histos1_["pfCandidate_pv_tsigETL"]->Sumw2();
  // histos1_["pfCandidate_pv_tqualETL"]->Sumw2();

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
  histos1_["nc"] = fs->make<TH1F>("nc", "nc", nclus_, 0, nclus_);
  histos1_["nc"]->Sumw2();
  // histos1_["ncmtdpv"] = fs->make<TH1F>("ncmtdpv", "ncmtdpv", nclus_, 0, nclus_);
  // histos1_["ncmtdpv"]->Sumw2();
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
    if (obj.BeginsWith("gv")) {
      if (obj.Contains("_trk"))
        book1DHistograms(fs, vars1_trks_, obj, histos1_);
      else
        book1DHistograms(fs, vars1_vtxs_gv_, obj, histos1_);
    }
    else if (obj.BeginsWith("sv")) {
      if (obj.Contains("_trk"))
        book1DHistograms(fs, vars1_trks_sv_, obj, histos1_);
      else
        book1DHistograms(fs, vars1_vtxs_sv_, obj, histos1_);
    }
    else if (obj.BeginsWith("gj") || obj.BeginsWith("rj"))
      book1DHistograms(fs, vars1_jets_, obj, histos1_);
  }

  // 2D histograms
  for (TString obj : sv_names_) {
    book2DHistograms(fs, vars2_trks_, obj + "_trk", histos2_);
    book2DHistograms(fs, vars2_vtxs_, obj, histos2_);
  }

  // Matching histograms
  for (TString objgv : gv_names_) {
    for (TString objsv : sv_names_) {
      TString matchobjs = objgv + "_" + objsv;
      book1DHistograms(fs, vars1_trks_, matchobjs + "_trk", histos1_);
      book1DHistograms(fs, vars1_trkVtxMatch_, matchobjs + "_trk", histos1_);
      book1DHistograms(fs, vars1_vtxs_gv_, matchobjs, histos1_);
      book1DHistograms(fs, vars1_trkVtxMatch_, matchobjs, histos1_);
      book2DHistograms(fs, vars2_trks_, matchobjs + "_trk", histos2_);
      book2DHistograms(fs, vars2_trkVtxMatch_, matchobjs + "_trk", histos2_);
      book2DHistograms(fs, vars2_vtxs_, matchobjs, histos2_);
      book2DHistograms(fs, vars2_trkVtxMatch_, matchobjs, histos2_);

      matchobjs = objsv + "_" + objgv;
      book1DHistograms(fs, vars1_trks_sv_, matchobjs + "_trk", histos1_);
      book1DHistograms(fs, vars1_trkVtxMatch_, matchobjs + "_trk", histos1_);
      book1DHistograms(fs, vars1_vtxs_sv_, matchobjs, histos1_);
      book1DHistograms(fs, vars1_trkVtxMatch_, matchobjs, histos1_);
      book2DHistograms(fs, vars2_trks_, matchobjs + "_trk", histos2_);
      book2DHistograms(fs, vars2_trkVtxMatch_, matchobjs + "_trk", histos2_);
      book2DHistograms(fs, vars2_vtxs_, matchobjs, histos2_);
      book2DHistograms(fs, vars2_trkVtxMatch_, matchobjs, histos2_);
    }
  }
}


VertexNtuplizer::~VertexNtuplizer() {}


void VertexNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  const reco::VertexCollection primaryVertices = iEvent.get(primaryVerticesToken_);
  // Sorting described here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
  const reco::Vertex& primaryVertex = primaryVertices.at(0); // Most likely the signal vertex

  unsigned int nPassingGVs[3] = {0, 0, 0};
  gvc_->build(iEvent, genParticlesToken_, simTracksToken_, primaryVertex, matcher_, nPassingGVs);
  svc_->build(iEvent,
      inclusiveSecondaryVerticesToken_, // inclusiveSecondaryVerticesMTDPVToken_,
      mergedSecondaryVerticesToken_, // mergedSecondaryVerticesMTDPVToken_,
      slimmedSecondaryVerticesToken_, // slimmedSecondaryVerticesMTDPVToken_,
      generalTracksToken_,
      trackT0FromBSToken_, trackSigmaT0FromBSToken_, // trackQualityFromBSToken_,
      trackT0FromPVToken_, trackSigmaT0FromPVToken_, // trackQualityFromPVToken_,
      primaryVertex);
  rjc_->build(iEvent, jetsToken_, genJetsFlavourInfoToken_);
  gjc_->build(iEvent, genJetsToken_, genJetsFlavourInfoToken_, jetsToken_);

  std::vector<GenVertexCollection> GVCollections;
  GVCollections.push_back(gvc_->getGenVertices());
  GVCollections.push_back(gvc_->getGenVerticesB());
  GVCollections.push_back(gvc_->getGenVerticesD());
  // GVCollections.push_back(gvc_->getGenVerticesSimMatch());
  // GVCollections.push_back(gvc_->getGenVerticesSimMatchB());
  // GVCollections.push_back(gvc_->getGenVerticesSimMatchD());
  for (unsigned int iColl = 0; iColl < GVCollections.size(); iColl++) {
    histos1_["n" + gv_names_.at(iColl)]->Fill(GVCollections.at(iColl).size());
    for (GenVertex& gv : GVCollections.at(iColl)) {
      gv.fill(histos1_, gv_names_.at(iColl));
    }
  }

  unsigned int nC = iEvent.get(nIVFClustersToken_);
  // unsigned int nCMTDPV = iEvent.get(nIVFClustersMTDPVToken_);
  histos1_["nc"]->Fill(nC);
  // histos1_["ncmtdpv"]->Fill(nCMTDPV);

  std::vector<SecondaryVertexCollection> SVCollections;
  SVCollections.push_back(svc_->getSecVerticesInclusive());
  SVCollections.push_back(svc_->getSecVerticesInclusiveMTDBS());
  SVCollections.push_back(svc_->getSecVerticesInclusiveMTDPV());
  SVCollections.push_back(svc_->getSecVerticesSlimmed());
  SVCollections.push_back(svc_->getSecVerticesSlimmedMTDBS());
  SVCollections.push_back(svc_->getSecVerticesSlimmedMTDPV());
  for (unsigned int iColl = 0; iColl < SVCollections.size(); iColl++) {
    histos1_["n" + sv_names_.at(iColl)]->Fill(SVCollections.at(iColl).size());
    for (SecondaryVertex& sv : SVCollections.at(iColl)) sv.fill(histos1_, histos2_, sv_names_.at(iColl));
  }

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
    if (gv_names_.at(iGVs).EndsWith("B") || gv_names_.at(iGVs).EndsWith("D")) continue;
    for (unsigned int iSVs = 0; iSVs < SVCollections.size(); iSVs++) {
      std::vector<bool> SVmatched(SVCollections.at(iSVs).size(), false); // Prevent double matching
      for (GenVertex& gv : GVCollections.at(iGVs)) {
        for (unsigned int iSV = 0; iSV < SVCollections.at(iSVs).size(); iSV++) {
          SecondaryVertex& sv = SVCollections.at(iSVs).at(iSV);

          TString gv_name, sv_name;
          if (!SVmatched.at(iSV) && matcher_->match(gv, sv, VertexMatcher::TRACK)) {
            SVmatched.at(iSV) = true; // Prevent double matching

            gv_name = gv_names_.at(iGVs) + "_" + sv_names_.at(iSVs);
            sv_name = sv_names_.at(iSVs) + "_" + gv_names_.at(iGVs);
            gv.fill(histos1_, gv_name);
            sv.fill(histos1_, histos2_, sv_name);
            matcher_->fill(histos1_, histos2_, gv_name, sv_name, gv, sv);

            TString gv_flav = gv_names_.at(iGVs);
            if (gv.pdgIdBin() == B_MESON || gv.pdgIdBin() == B_BARYON) gv_flav += "B";
            if (gv.pdgIdBin() == C_MESON || gv.pdgIdBin() == C_BARYON) gv_flav += "D";
            gv_name = gv_flav + "_" + sv_names_.at(iSVs);
            sv_name = sv_names_.at(iSVs) + "_" + gv_flav;
            gv.fill(histos1_, gv_name);
            sv.fill(histos1_, histos2_, sv_name);
            matcher_->fill(histos1_, histos2_, gv_name, sv_name, gv, sv);

            break; // Restrict to only one match
          } // End if matched
        } // End loop through SVs
      } // End loop through GVs
    } // End loop through SV collections
  } // End loop through GV collections

  // Fill standalone histograms
  histos1_["nGPs"]->Fill(nPassingGVs[0]);
  histos1_["nGPsB"]->Fill(nPassingGVs[1]);
  histos1_["nGPsD"]->Fill(nPassingGVs[2]);
  histos1_["fracReconstructableGVs"]->Fill(((float) gvc_->getGenVertices().size()) / ((float) nPassingGVs[0])); // TODO: FIX THIS TOO
  histos1_["fracReconstructableGVsB"]->Fill(((float) gvc_->getGenVerticesB().size()) / ((float) nPassingGVs[1])); // TODO: FIX THIS TOO
  histos1_["fracReconstructableGVsD"]->Fill(((float) gvc_->getGenVerticesD().size()) / ((float) nPassingGVs[2])); // TODO: FIX THIS TOO

  float maxMatchedWithInclusive = std::min(svc_->getSecVerticesInclusive().size(), gvc_->getGenVertices().size());
  float maxMatchedWithSlimmed = std::min(svc_->getSecVerticesSlimmed().size(), gvc_->getGenVertices().size());
  histos1_["maxMatchedWithInclusive"]->Fill(maxMatchedWithInclusive);
  histos1_["maxMatchedWithSlimmed"]->Fill(maxMatchedWithSlimmed);

  float maxMatchEffInclusive = maxMatchedWithInclusive / ((float) gvc_->getGenVertices().size());
  float maxMatchEffSlimmed = maxMatchedWithSlimmed / ((float) gvc_->getGenVertices().size());
  histos1_["maxMatchEffInclusive"]->Fill(maxMatchEffInclusive);
  histos1_["maxMatchEffSlimmed"]->Fill(maxMatchEffSlimmed);

  const reco::PFCandidateCollection pfCandidates = iEvent.get(pfCandidatesToken_);
  edm::ValueMap<float> trackT0FromBS_ = iEvent.get(trackT0FromBSToken_);
  edm::ValueMap<float> trackSigmaT0FromBS_ = iEvent.get(trackSigmaT0FromBSToken_);
  // edm::ValueMap<float> trackQualityFromBS_ = iEvent.get(trackQualityFromBSToken_);
  edm::ValueMap<float> trackT0FromPV_ = iEvent.get(trackT0FromPVToken_);
  edm::ValueMap<float> trackSigmaT0FromPV_ = iEvent.get(trackSigmaT0FromPVToken_);
  // edm::ValueMap<float> trackQualityFromPV_ = iEvent.get(trackQualityFromPVToken_);
  for (const reco::PFCandidate& pfc : pfCandidates) {
    if (pfc.pt() < 0.8) continue;
    if (abs(pfc.eta()) > 3.0) continue;
    if (pfc.charge() == 0) continue;
    const reco::TrackRef& trkRef = pfc.trackRef();
    if (!trkRef) continue;
    if (trackSigmaT0FromBS_[trkRef] > 0.0) {
      TString suffix = "";
      if (abs(pfc.eta()) < 1.5) suffix += "BTL";
      else suffix += "ETL";
      histos1_["pfCandidate_bs_tval" + suffix]->Fill(trackT0FromBS_[trkRef] * 1000.0);
      histos1_["pfCandidate_bs_terr" + suffix]->Fill(trackSigmaT0FromBS_[trkRef] * 1000.0);
      histos1_["pfCandidate_bs_tsig" + suffix]->Fill(trackT0FromBS_[trkRef] / trackSigmaT0FromBS_[trkRef]);
      // histos1_["pfCandidate_bs_tqual" + suffix]->Fill(trackQualityFromBS_[trkRef]);
    }
    if (trackSigmaT0FromPV_[trkRef] > 0.0) {
      TString suffix = "";
      if (abs(pfc.eta()) < 1.5) suffix += "BTL";
      else suffix += "ETL";
      histos1_["pfCandidate_pv_tval" + suffix]->Fill(trackT0FromPV_[trkRef] * 1000.0);
      histos1_["pfCandidate_pv_terr" + suffix]->Fill(trackSigmaT0FromPV_[trkRef] * 1000.0);
      histos1_["pfCandidate_pv_tsig" + suffix]->Fill(trackT0FromPV_[trkRef] / trackSigmaT0FromPV_[trkRef]);
      // histos1_["pfCandidate_pv_tqual" + suffix]->Fill(trackQualityFromPV_[trkRef]);
    }
  }
}


void VertexNtuplizer::book1DHistograms(edm::Service<TFileService>& fs,
    std::map<TString, std::vector<float>>& vars,
    const TString& obj, std::map<TString, TH1F*>& histos) {

  for (const auto& iter : vars) {
    TString name = obj + "_" + iter.first;
    histos[name] = fs->make<TH1F>(name, name,
        iter.second[0], iter.second[1], iter.second[2]);
    histos[name]->Sumw2();
  }
}


void VertexNtuplizer::book2DHistograms(edm::Service<TFileService>& fs,
    std::map<TString, std::vector<float>>& vars,
    const TString& obj, std::map<TString, TH2F*>& histos) {

  for (const auto& iter : vars) {
    TString name = obj + "_" + iter.first;
    histos[name] = fs->make<TH2F>(name, name,
        iter.second[0], iter.second[1], iter.second[2],
        iter.second[3], iter.second[4], iter.second[5]);
    histos[name]->Sumw2();
  }
}


void VertexNtuplizer::beginJob() {}


void VertexNtuplizer::endJob() {

  unsigned int nEvents = histos1_["nGPs"]->GetEntries();

  std::cout << std::endl;
  std::cout << "number of events = " << nEvents << std::endl;

  std::cout << std::endl;
  std::cout << "mean number of GVs                              = " << histos1_["ngv"]->GetMean() << "+=" << histos1_["ngv"]->GetStdDev() << std::endl;
  std::cout << "mean number of GVs from B mother                = " << histos1_["ngvB"]->GetMean() << "+=" << histos1_["ngvB"]->GetStdDev() << std::endl;
  std::cout << "mean number of GVs from D mother                = " << histos1_["ngvD"]->GetMean() << "+=" << histos1_["ngvD"]->GetStdDev() << std::endl;
  // std::cout << "mean number of GVs with SIM match               = " << histos1_["ngvs"]->GetMean() << "+=" << histos1_["ngvs"]->GetStdDev() << std::endl;
  // std::cout << "mean number of GVs from B mother with SIM match = " << histos1_["ngvsB"]->GetMean() << "+=" << histos1_["ngvsB"]->GetStdDev() << std::endl;
  // std::cout << "mean number of GVs from D mother with SIM match = " << histos1_["ngvsD"]->GetMean() << "+=" << histos1_["ngvsD"]->GetStdDev() << std::endl;

  std::cout << std::endl;
  std::cout << "mean number of inclusive SVs        = " << histos1_["nsvInclusive"]->GetMean() << "+=" << histos1_["nsvInclusive"]->GetStdDev() << std::endl;
  std::cout << "mean number of inclusive SVs MTD PV = " << histos1_["nsvInclusiveMTDPV"]->GetMean() << "+=" << histos1_["nsvInclusiveMTDPV"]->GetStdDev() << std::endl;
  std::cout << "mean number of slimmed SVs          = " << histos1_["nsvSlimmed"]->GetMean() << "+=" << histos1_["nsvSlimmed"]->GetStdDev() << std::endl;
  std::cout << "mean number of slimmed SVs MTD PV   = " << histos1_["nsvSlimmedMTDPV"]->GetMean() << "+=" << histos1_["nsvSlimmedMTDPV"]->GetStdDev() << std::endl;

  std::cout << std::endl;
  std::cout << "inclusive SVs" << std::endl;
  std::cout << "Integrated efficiency for all GVs                              = " << histos1_["gv_svInclusive_pdgIdBin"]->GetEntries() / histos1_["gv_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs from B mother                = " << histos1_["gvB_svInclusive_pdgIdBin"]->GetEntries() / histos1_["gvB_pdgIdBin"]->GetEntries() << std::endl;
  std::cout << "Integrated efficiency for all GVs from D mother                = " << histos1_["gvD_svInclusive_pdgIdBin"]->GetEntries() / histos1_["gvD_pdgIdBin"]->GetEntries() << std::endl;
  // std::cout << "Integrated efficiency for all GVs SimTrack match               = " << histos1_["gvs_svInclusive_pdgIdBin"]->GetEntries() / histos1_["gvs_pdgIdBin"]->GetEntries() << std::endl;
  // std::cout << "Integrated efficiency for all GVs from B mother SimTrack match = " << histos1_["gvsB_svInclusive_pdgIdBin"]->GetEntries() / histos1_["gvsB_pdgIdBin"]->GetEntries() << std::endl;
  // std::cout << "Integrated efficiency for all GVs from D mother SimTrack match = " << histos1_["gvsD_svInclusive_pdgIdBin"]->GetEntries() / histos1_["gvsD_pdgIdBin"]->GetEntries() << std::endl;

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
