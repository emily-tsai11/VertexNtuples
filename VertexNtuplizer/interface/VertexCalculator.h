#ifndef VertexNtuples_VertexNtuplizer_VertexCalculator_h
#define VertexNtuples_VertexNtuplizer_VertexCalculator_h


#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "TMath.h"


namespace vertexntuples {

  const double catchInfsAndBound(float input, float replaceValue,
      float lowerBound, float upperBound,
      float offset, bool useOffsets);

  const double dxy(const reco::Candidate* dau, const reco::Vertex& v);
  const double dz(const reco::Candidate* dau, const reco::Vertex& v);
  const double d3d(const reco::Candidate* dau, const reco::Vertex& v);
  const double dxyErr(const reco::Candidate* dau, const reco::Vertex& v);
  const double dzErr(const reco::Candidate* dau, const reco::Vertex& v);
  const double d3dErr(const reco::Candidate* dau, const reco::Vertex& v);

  const double dxy(const HepMC::FourVector& pos, const reco::Vertex& v);
  const double dz(const HepMC::FourVector& pos, const reco::Vertex& v);
  const double d3d(const HepMC::FourVector& pos, const reco::Vertex& v);
  const double dxyErr(const HepMC::FourVector& pos, const reco::Vertex& v);
  const double dzErr(const HepMC::FourVector& pos, const reco::Vertex& v);
  const double d3dErr(const HepMC::FourVector& pos, const reco::Vertex& v);

  const double d3d(const reco::TrackBaseRef& trkRef);
  const double d3dErr(const reco::TrackBaseRef& trkRef);

  Measurement1D dxy(const reco::Vertex& v1, const reco::Vertex& v2);
  Measurement1D d3d(const reco::Vertex& v1, const reco::Vertex& v2);
  const double dz(const reco::Vertex& v1, const reco::Vertex& v2);
  const double dzErr(const reco::Vertex& v1, const reco::Vertex& v2);

  // Prefer not to use these -- use the above
  const double dxy(const double x1, const double x2, const double y1, const double y2);
  const double dz(const double z1, const double z2);
  const double d3d(const double dxy, const double dz);
  const double dxyErr(const double x1, const double x2, const double xerr1, const double xerr2,
      const double y1, const double y2, const double yerr1, const double yerr2);
  const double dzErr(const double zerr1, const double zerr2);
  const double d3dErr(const double dxy, const double dxyerr, const double dz, const double dzerr);
}


#endif
