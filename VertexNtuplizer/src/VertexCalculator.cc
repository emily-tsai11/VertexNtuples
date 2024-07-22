#include "../interface/VertexCalculator.h"

#include "../interface/GenVertex.h"
#include "../interface/SecondaryVertex.h"


namespace vertexntuples {

  // const float catchInfsAndBound(float input, float replaceValue,
  //     float lowerBound, float upperBound,
  //     float offset, bool useOffset) {
  //   float returnValue = input;
  //   if (edm::isNotFinite(returnValue)) returnValue = replaceValue;
  //   if (returnValue < -1e32 || returnValue > 1e32) returnValue = replaceValue;
  //   if (returnValue + offset < lowerBound) returnValue = lowerBound;
  //   if (returnValue + offset > upperBound) returnValue = upperBound;
  //   if (useOffset) returnValue += offset;
  //   return returnValue;
  // }

  // ---------------------------------------------------------------------- //

  namespace {

    const float dxy(const float x1, const float x2, const float y1, const float y2) {
      const float dx = abs(x1 - x2);
      const float dy = abs(y1 - y2);
      return TMath::Sqrt(dx * dx + dy * dy);
    }

    const float dxyErr(const float x1, const float x2, const float xerr1, const float xerr2,
        const float y1, const float y2, const float yerr1, const float yerr2) {
      const float dx = abs(x1 - x2);
      const float dy = abs(y1 - y2);
      const float dxerr = TMath::Sqrt(xerr1 * xerr1 + xerr2 * xerr2);
      const float dyerr = TMath::Sqrt(yerr1 * yerr1 + yerr2 * yerr2);
      const float dx2err = 2 * dx * dxerr;
      const float dy2err = 2 * dy * dyerr;
      const float dxy2err = TMath::Sqrt(dx2err * dx2err + dy2err * dy2err);
      return 0.5 * dxy2err / dxy(x1, x2, y1, y2);
    }

    const float dz(const float z1, const float z2) {
      return abs(z1 - z2);
    }

    const float dzErr(const float zerr1, const float zerr2) {
      return TMath::Sqrt(zerr1 * zerr1 + zerr2 * zerr2);
    }

    const float d3d(const float dxy, const float dz) {
      return TMath::Sqrt(dxy * dxy + dz * dz);
    }

    const float d3dErr(const float dxy, const float dxyerr, const float dz, const float dzerr) {
      const float dxy2err = 2 * dxy * dxyerr;
      const float dz2err = 2 * dz * dzerr;
      const float d3d2err = TMath::Sqrt(dxy2err * dxy2err + dz2err * dz2err);
      return 0.5 * d3d2err / d3d(dxy, dz);
    }
  }

  // ---------------------------------------------------------------------- //

  const float dxy(const reco::Candidate* dau, const reco::Vertex& v) {
    return dxy(dau->vx(), v.x(), dau->vy(), v.y());
  }

  const float dxyErr(const reco::Candidate* dau, const reco::Vertex& v) {
    // 0.0 error for generated daughter
    return dxyErr(dau->vx(), v.x(), 0.0, v.xError(), dau->vy(), v.y(), 0.0, v.yError());
  }

  const float dz(const reco::Candidate* dau, const reco::Vertex& v) {
    return dz(dau->vz(), v.z());
  }

  const float dzErr(const reco::Candidate* dau, const reco::Vertex& v) {
    // 0.0 error for generated daughter
    return dzErr(0.0, v.zError());
  }

  const float d3d(const reco::Candidate* dau, const reco::Vertex& v) {
    const float dxyval = dxy(dau, v);
    const float dzval = dz(dau, v);
    return d3d(dxyval, dzval);
  }

  const float d3dErr(const reco::Candidate* dau, const reco::Vertex& v) {
    const float dxyval = dxy(dau, v);
    const float dxyerr = dxyErr(dau, v);
    const float dzval = dz(dau, v);
    const float dzerr = dzErr(dau, v);
    return d3dErr(dxyval, dxyerr, dzval, dzerr);
  }

  // ---------------------------------------------------------------------- //

  const float d3d(const reco::TrackBaseRef& trkRef) {
    return d3d(trkRef->dxy(), trkRef->dz());
  }

  const float d3dErr(const reco::TrackBaseRef& trkRef) {
    return d3dErr(trkRef->dxy(), trkRef->dxyError(), trkRef->dz(), trkRef->dzError());
  }

  // ---------------------------------------------------------------------- //

  const float dxy(const reco::CandidatePtr& trkPtr, const reco::Vertex& v) {
    return dxy(trkPtr->vx(), v.x(), trkPtr->vy(), v.y());
  }

  const float dz(const reco::CandidatePtr& trkPtr, const reco::Vertex& v) {
    return dz(trkPtr->vz(), v.z());
  }

  const float d3d(const reco::CandidatePtr& trkPtr, const reco::Vertex& v) {
    const float dxyval = dxy(trkPtr, v);
    const float dzval = dz(trkPtr, v);
    return d3d(dxyval, dzval);
  }

  const float d3dErr(const reco::CandidatePtr& trkPtr, const reco::Vertex& v) {
    const float dxyval = dxy(trkPtr, v);
    const float dzval = dz(trkPtr, v);
    return d3dErr(dxyval, trkPtr->dxyError(), dzval, trkPtr->dzError());
  }

  // ---------------------------------------------------------------------- //

  Measurement1D dxy(const reco::Vertex& v1, const reco::Vertex& v2) {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix cm1;
    v1.fill(cm1);
    reco::Vertex v1new(v1.position(), cm1);
    return dist.distance(v1new, v2);
  }

  const float dz(const reco::Vertex& v1, const reco::Vertex& v2) {
    return dz(v1.z(), v2.z());
  }

  const float dzErr(const reco::Vertex& v1, const reco::Vertex& v2) {
    return dzErr(v1.zError(), v2.zError());
  }

  Measurement1D d3d(const reco::Vertex& v1, const reco::Vertex& v2) {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix cm1;
    v1.fill(cm1);
    reco::Vertex v1new(v1.position(), cm1);
    return dist.distance(v1new, v2);
  }

  // ---------------------------------------------------------------------- //

  Measurement1D dxy(const reco::VertexCompositePtrCandidate& v1, const reco::Vertex& v2) {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix cm1;
    v1.fillVertexCovariance(cm1);
    reco::Vertex v1new(v1.vertex(), cm1);
    return dist.distance(v1new, v2);
  }

  const float dz(const reco::VertexCompositePtrCandidate& v1, const reco::Vertex& v2) {
    return dz(v1.vz(), v2.z());
  }

  Measurement1D d3d(const reco::VertexCompositePtrCandidate& v1, const reco::Vertex& v2) {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix cm1;
    v1.fillVertexCovariance(cm1);
    reco::Vertex v1new(v1.vertex(), cm1);
    return dist.distance(v1new, v2);
  }

  // ---------------------------------------------------------------------- //

  const float xres(const GenVertex& gv, const SecondaryVertex& sv) {
    return gv.x() - sv.x();
  }

  const float yres(const GenVertex& gv, const SecondaryVertex& sv) {
    return gv.y() - sv.y();
  }

  const float zres(const GenVertex& gv, const SecondaryVertex& sv) {
    return gv.z() - sv.z();
  }

  const float xpull(const GenVertex& gv, const SecondaryVertex& sv) {
    return (gv.x() - sv.x()) / sv.xErr();
  }

  const float ypull(const GenVertex& gv, const SecondaryVertex& sv) {
    return (gv.y() - sv.y()) / sv.yErr();
  }

  const float zpull(const GenVertex& gv, const SecondaryVertex& sv) {
    return (gv.z() - sv.z()) / sv.zErr();
  }

  // ---------------------------------------------------------------------- //

  const float dxy(const GenVertex& gv, const SecondaryVertex& sv) {
    return dxy(gv.x(), sv.x(), gv.y(), sv.y());
  }

  const float dxyErr(const GenVertex& gv, const SecondaryVertex& sv) {
    return dxyErr(gv.x(), sv.x(), gv.xErr(), sv.xErr(), gv.y(), sv.y(), gv.yErr(), sv.yErr());
  }

  const float dz(const GenVertex& gv, const SecondaryVertex& sv) {
    return dz(gv.z(), sv.z());
  }

  const float dzErr(const GenVertex& gv, const SecondaryVertex& sv) {
    return dzErr(gv.zErr(), sv.zErr());
  }

  const float d3d(const GenVertex& gv, const SecondaryVertex& sv) {
    const float dxyval = dxy(gv, sv);
    const float dzval = dz(gv, sv);
    return d3d(dxyval, dzval);
  }

  const float d3dErr(const GenVertex& gv, const SecondaryVertex& sv) {
    const float dxyval = dxy(gv, sv);
    const float dxyerr = dxyErr(gv, sv);
    const float dzval = dz(gv, sv);
    const float dzerr = dzErr(gv, sv);
    return d3dErr(dxyval, dxyerr, dzval, dzerr);
  }
}
