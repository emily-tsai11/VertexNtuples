#include "../interface/VertexMatcher.h"


VertexMatcher::VertexMatcher(const edm::ParameterSet& iConfig) {

  trkMatchDrCut_ = iConfig.getUntrackedParameter<double>("recoTrkMatchDrCut");
  trkMatchPtCut_ = iConfig.getUntrackedParameter<double>("recoTrkMatchPtCut");
  vtxMatchFrac_ = iConfig.getUntrackedParameter<double>("vtxMatchFrac");
  jetRadius_ = iConfig.getUntrackedParameter<double>("jetRadius");
}


// VertexMatcher::~VertexMatcher() {}


bool VertexMatcher::match(GenVertex& gv, SecondaryVertex& sv, MatchAlgo algo) {

  if (algo == TRACK) {
    return vtxTrackMatch(gv, sv);
  }
  else if (algo == MATRIX) {
    return vtxMatrixMatch(gv, sv);
  }
  std::cout << "WARNING in VertexMatcher::match(gv, sv): No matching algorithm specified. Returning false." << std::endl;
  return false;
}


bool VertexMatcher::match(SecondaryVertex& sv, RecoJet& rj) {

  return reco::deltaR(sv.eta(), sv.phi(), rj.eta(), rj.phi()) < jetRadius_;
}


bool VertexMatcher::vtxTrackMatch(GenVertex& gv, SecondaryVertex& sv) {

  float nmatch = 0;
  for (const reco::Candidate* dau : *(gv.daughters())) {
    for (unsigned int iTrk = 0; iTrk < sv.nTracks(); iTrk++) {
      bool match = true;
      if (reco::deltaR(sv.trkEta()->at(iTrk), sv.trkPhi()->at(iTrk), dau->eta(), dau->phi()) > trkMatchDrCut_) match = false;
      if (abs(dau->pt() - sv.trkPt()->at(iTrk)) / (dau->pt() + sv.trkPt()->at(iTrk)) > trkMatchPtCut_) match = false;
      if (match) {
        nmatch++;
        break;
      }
    }
  }
  return (nmatch/(float)gv.nDaughters()) >= vtxMatchFrac_;
}


bool VertexMatcher::vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv) {

  // Not yet implemented
  return false;
}


void VertexMatcher::fill(std::map<TString, TH1F*>& histos, TString gvPrefix,
    TString svPrefix, GenVertex& gv, SecondaryVertex& sv) {

  double xres = vertexntuples::xres(gv, sv);
  double yres = vertexntuples::yres(gv, sv);
  double zres = vertexntuples::zres(gv, sv);
  double xpull = vertexntuples::xpull(gv, sv);
  double ypull = vertexntuples::ypull(gv, sv);
  double zpull = vertexntuples::zpull(gv, sv);
  double dxy = vertexntuples::dxy(gv, sv);
  double d3d = vertexntuples::d3d(gv, sv);
  double dxyerr = vertexntuples::dxyErr(gv, sv);
  double d3derr = vertexntuples::d3dErr(gv, sv);
  double dxysig = dxy / dxyerr;
  double d3dsig = d3d / d3derr;

  histos[gvPrefix + "_xres"]->Fill(xres);
  histos[gvPrefix + "_yres"]->Fill(yres);
  histos[gvPrefix + "_zres"]->Fill(zres);
  histos[gvPrefix + "_xpull"]->Fill(xpull);
  histos[gvPrefix + "_ypull"]->Fill(ypull);
  histos[gvPrefix + "_zpull"]->Fill(zpull);
  histos[gvPrefix + "_matchdxy"]->Fill(dxy);
  histos[gvPrefix + "_matchd3d"]->Fill(d3d);
  histos[gvPrefix + "_matchdxyerr"]->Fill(dxyerr);
  histos[gvPrefix + "_matchd3derr"]->Fill(d3derr);
  histos[gvPrefix + "_matchdxysig"]->Fill(dxysig);
  histos[gvPrefix + "_matchd3dsig"]->Fill(d3dsig);

  histos[svPrefix + "_xres"]->Fill(xres);
  histos[svPrefix + "_yres"]->Fill(yres);
  histos[svPrefix + "_zres"]->Fill(zres);
  histos[svPrefix + "_xpull"]->Fill(xpull);
  histos[svPrefix + "_ypull"]->Fill(ypull);
  histos[svPrefix + "_zpull"]->Fill(zpull);
  histos[svPrefix + "_matchdxy"]->Fill(dxy);
  histos[svPrefix + "_matchd3d"]->Fill(d3d);
  histos[svPrefix + "_matchdxyerr"]->Fill(dxyerr);
  histos[svPrefix + "_matchd3derr"]->Fill(d3derr);
  histos[svPrefix + "_matchdxysig"]->Fill(dxysig);
  histos[svPrefix + "_matchd3dsig"]->Fill(d3dsig);
}
