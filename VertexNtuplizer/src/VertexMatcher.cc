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
  return (nmatch / (float) gv.nDaughters()) >= vtxMatchFrac_;
}


bool VertexMatcher::vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv) {

  // Not yet implemented
  return false;
}


void VertexMatcher::fill(std::map<TString, TH1F*>& histos1, std::map<TString, TH2F*>& histos2,
    TString gvPrefix, TString svPrefix, GenVertex& gv, SecondaryVertex& sv) {

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

  histos1[gvPrefix + "_xres"]->Fill(xres);
  histos1[gvPrefix + "_yres"]->Fill(yres);
  histos1[gvPrefix + "_zres"]->Fill(zres);
  histos1[gvPrefix + "_xpull"]->Fill(xpull);
  histos1[gvPrefix + "_ypull"]->Fill(ypull);
  histos1[gvPrefix + "_zpull"]->Fill(zpull);
  histos1[gvPrefix + "_matchdxy"]->Fill(dxy);
  histos1[gvPrefix + "_matchd3d"]->Fill(d3d);
  histos1[gvPrefix + "_matchdxyerr"]->Fill(dxyerr);
  histos1[gvPrefix + "_matchd3derr"]->Fill(d3derr);
  histos1[gvPrefix + "_matchdxysig"]->Fill(dxysig);
  histos1[gvPrefix + "_matchd3dsig"]->Fill(d3dsig);

  histos2[gvPrefix + "_matchdxy_dxy"]->Fill(dxy, gv.dxy());
  histos2[gvPrefix + "_matchdxy_d3d"]->Fill(dxy, gv.d3d());
  histos2[gvPrefix + "_matchd3d_dxy"]->Fill(d3d, gv.dxy());
  histos2[gvPrefix + "_matchd3d_d3d"]->Fill(d3d, gv.d3d());

  histos1[svPrefix + "_xres"]->Fill(xres);
  histos1[svPrefix + "_yres"]->Fill(yres);
  histos1[svPrefix + "_zres"]->Fill(zres);
  histos1[svPrefix + "_xpull"]->Fill(xpull);
  histos1[svPrefix + "_ypull"]->Fill(ypull);
  histos1[svPrefix + "_zpull"]->Fill(zpull);
  histos1[svPrefix + "_matchdxy"]->Fill(dxy);
  histos1[svPrefix + "_matchd3d"]->Fill(d3d);
  histos1[svPrefix + "_matchdxyerr"]->Fill(dxyerr);
  histos1[svPrefix + "_matchd3derr"]->Fill(d3derr);
  histos1[svPrefix + "_matchdxysig"]->Fill(dxysig);
  histos1[svPrefix + "_matchd3dsig"]->Fill(d3dsig);

  histos2[svPrefix + "_matchdxy_dxy"]->Fill(dxy, sv.dxy());
  histos2[svPrefix + "_matchdxy_d3d"]->Fill(dxy, sv.d3d());
  histos2[svPrefix + "_matchd3d_dxy"]->Fill(d3d, sv.dxy());
  histos2[svPrefix + "_matchd3d_d3d"]->Fill(d3d, sv.d3d());
}
