#ifndef VertexNtuples_VertexNtuplizer_VertexCalculator_h
#define VertexNtuples_VertexNtuplizer_VertexCalculator_h


#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "TMath.h"


class GenVertex;
class SecondaryVertex;


namespace vertexntuples {

  // const float catchInfsAndBound(float input, float replaceValue,
  //     float lowerBound, float upperBound,
  //     float offset, bool useOffsets);

  const float dxy(const reco::Candidate* dau, const reco::Vertex& v);
  const float dxyErr(const reco::Candidate* dau, const reco::Vertex& v);
  const float dz(const reco::Candidate* dau, const reco::Vertex& v);
  const float dzErr(const reco::Candidate* dau, const reco::Vertex& v);
  const float d3d(const reco::Candidate* dau, const reco::Vertex& v);
  const float d3dErr(const reco::Candidate* dau, const reco::Vertex& v);

  const float d3d(const reco::TrackBaseRef& trkRef);
  const float d3dErr(const reco::TrackBaseRef& trkRef);

  const float dxy(const reco::CandidatePtr& trkPtr, const reco::Vertex& v);
  const float dz(const reco::CandidatePtr& trkPtr, const reco::Vertex& v);
  const float d3d(const reco::CandidatePtr& trkPtr, const reco::Vertex& v);
  const float d3dErr(const reco::CandidatePtr& trkPtr, const reco::Vertex& v);

  Measurement1D dxy(const reco::Vertex& v1, const reco::Vertex& v2);
  const float dz(const reco::Vertex& v1, const reco::Vertex& v2);
  const float dzErr(const reco::Vertex& v1, const reco::Vertex& v2);
  Measurement1D d3d(const reco::Vertex& v1, const reco::Vertex& v2);

  Measurement1D dxy(const reco::VertexCompositePtrCandidate& v1, const reco::Vertex& v2);
  const float dz(const reco::VertexCompositePtrCandidate& v1, const reco::Vertex& v2);
  Measurement1D d3d(const reco::VertexCompositePtrCandidate& v1, const reco::Vertex& v2);

  const float xres(const GenVertex& gv, const SecondaryVertex& sv);
  const float yres(const GenVertex& gv, const SecondaryVertex& sv);
  const float zres(const GenVertex& gv, const SecondaryVertex& sv);
  const float xpull(const GenVertex& gv, const SecondaryVertex& sv);
  const float ypull(const GenVertex& gv, const SecondaryVertex& sv);
  const float zpull(const GenVertex& gv, const SecondaryVertex& sv);

  const float dxy(const GenVertex& gv, const SecondaryVertex& sv);
  const float dxyErr(const GenVertex& gv, const SecondaryVertex& sv);
  const float dz(const GenVertex& gv, const SecondaryVertex& sv);
  const float dzErr(const GenVertex& gv, const SecondaryVertex& sv);
  const float d3d(const GenVertex& gv, const SecondaryVertex& sv);
  const float d3dErr(const GenVertex& gv, const SecondaryVertex& sv);
}


#endif
