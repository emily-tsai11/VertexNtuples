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

  // Do not use these -- use the above
  const float dxy(const float x1, const float x2, const float y1, const float y2);
  const float dxyErr(const float x1, const float x2, const float xerr1, const float xerr2,
      const float y1, const float y2, const float yerr1, const float yerr2);
  const float dz(const float z1, const float z2);
  const float dzErr(const float zerr1, const float zerr2);
  const float d3d(const float dxy, const float dz);
  const float d3dErr(const float dxy, const float dxyerr, const float dz, const float dzerr);
}


#endif
