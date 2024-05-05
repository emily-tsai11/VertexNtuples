// -*- C++ -*-
//
// Package:    VertexNtuples/VertexNtuplizer
// Class:      VertexNtuplizer
//
/**\class VertexNtuplizer VertexNtuplizer.cc VertexNtuples/VertexNtuplizer/plugins/VertexNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emily Minyun Tsai
//         Created:  Sat, 04 May 2024 12:05:41 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "../interface/GenVertex.h"
#include "../interface/GenVertexCollectionBuilder.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class VertexNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit VertexNtuplizer(const edm::ParameterSet&);
    ~VertexNtuplizer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file

    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;

  #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    edm::ESGetToken<SetupData, SetupRecord> setupToken_;
  #endif

    GenVertexCollectionBuilder* gvc;

    TH1I* histo;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VertexNtuplizer::VertexNtuplizer(const edm::ParameterSet& iConfig) :
    tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"))),
    simTracksToken_(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simTracks"))) {

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  histo = fs->make<TH1I>("charge", "Charges", 2, -1, 1);

  gvc = new GenVertexCollectionBuilder(iConfig);
}

VertexNtuplizer::~VertexNtuplizer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void VertexNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  gvc->build(iEvent, genParticlesToken_, simTracksToken_);

  GenVertexCollection genVertices = gvc->getGenVertexCollection();
  GenVertexCollection genVerticesSimMatch = gvc->getGenVertexSimMatchCollection();
  GenVertexCollection genVerticesNoNu = gvc->getGenVertexNoNuCollection();
  GenVertexCollection genVerticesNoNuSimMatch = gvc->getGenVertexNoNuSimMatchCollection();

  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
    histo->Fill(track.charge());
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void VertexNtuplizer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void VertexNtuplizer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void VertexNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexNtuplizer);
