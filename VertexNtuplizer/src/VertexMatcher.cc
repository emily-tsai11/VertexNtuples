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
  // else if (algo == MATRIX) {
  //   return vtxMatrixMatch(gv, sv);
  // }
  std::cout << "WARNING in VertexMatcher::match(gv, sv): No matching algorithm specified. Returning false." << std::endl;
  return false;
}


bool VertexMatcher::match(SecondaryVertex& sv, RecoJet& rj) {

  return reco::deltaR(sv.eta(), sv.phi(), rj.eta(), rj.phi()) < jetRadius_;
}


bool VertexMatcher::vtxTrackMatch(GenVertex& gv, SecondaryVertex& sv) {

  float nmatch = 0;
  std::vector<bool> trkMatched(sv.nTracks(), false);
  for (const reco::Candidate* dau : *(gv.daughters())) {
    for (unsigned int iTrk = 0; iTrk < sv.nTracks(); iTrk++) {
      bool match = true;
      float matchdR = reco::deltaR(sv.trkEta()->at(iTrk), sv.trkPhi()->at(iTrk), dau->eta(), dau->phi());
      float matchPtResNorm = abs(dau->pt() - sv.trkPt()->at(iTrk)) / (dau->pt() + sv.trkPt()->at(iTrk));
      if (matchdR > trkMatchDrCut_) match = false;
      if (matchPtResNorm > trkMatchPtCut_) match = false;
      if (match && !trkMatched.at(iTrk)) {
        trkMatched.at(iTrk) = true;
        gv.addPtResNorm(matchPtResNorm);
        gv.addDeltaR(matchdR);
        sv.addPtResNorm(matchPtResNorm);
        sv.addDeltaR(matchdR);
        nmatch++;
        break;
      }
    }
  }
  return (nmatch / (float) gv.nDaughters()) >= vtxMatchFrac_;
}


// Not yet implemented
// bool VertexMatcher::vtxMatrixMatch(GenVertex& gv, SecondaryVertex& sv) {

//   return false;
// }


void VertexMatcher::fill(std::map<TString, TH1F*>& histos1, std::map<TString, TH2F*>& histos2,
    TString gvPrefix, TString svPrefix, const GenVertex& gv, const SecondaryVertex& sv) {

  float xres = vertexntuples::xres(gv, sv);
  float yres = vertexntuples::yres(gv, sv);
  float zres = vertexntuples::zres(gv, sv);
  float xpull = vertexntuples::xpull(gv, sv);
  float ypull = vertexntuples::ypull(gv, sv);
  float zpull = vertexntuples::zpull(gv, sv);
  float dxy = vertexntuples::dxy(gv, sv);
  float d3d = vertexntuples::d3d(gv, sv);
  float dxyerr = vertexntuples::dxyErr(gv, sv);
  float d3derr = vertexntuples::d3dErr(gv, sv);
  float dxysig = dxy / dxyerr;
  float d3dsig = d3d / d3derr;

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

  for (unsigned int iTrk = 0; iTrk < gv.dauMatchDeltaR()->size(); iTrk++) {
    histos1[gvPrefix + "_trk_deltaR"]->Fill(gv.dauMatchDeltaR()->at(iTrk));
    histos1[gvPrefix + "_trk_ptResNorm"]->Fill(gv.dauMatchPtResNorm()->at(iTrk));
    histos2[gvPrefix + "_trk_deltaR_ptResNorm"]->Fill(gv.dauMatchDeltaR()->at(iTrk), gv.dauMatchPtResNorm()->at(iTrk));
  }
  for (unsigned int iTrk = 0; iTrk < sv.trkMatchDeltaR()->size(); iTrk++) {
    histos1[svPrefix + "_trk_deltaR"]->Fill(sv.trkMatchDeltaR()->at(iTrk));
    histos1[svPrefix + "_trk_ptResNorm"]->Fill(sv.trkMatchPtResNorm()->at(iTrk));
    histos2[svPrefix + "_trk_deltaR_ptResNorm"]->Fill(sv.trkMatchDeltaR()->at(iTrk), sv.trkMatchPtResNorm()->at(iTrk));
  }
}
