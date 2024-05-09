#include "../interface/VertexCalculator.h"


namespace vertexntuples {

  const double catchInfsAndBound(float input, float replaceValue,
      float lowerBound, float upperBound,
      float offset, bool useOffset) {

    float returnValue = input;
    if (edm::isNotFinite(returnValue)) returnValue = replaceValue;
    if (returnValue < -1e32 || returnValue > 1e32) returnValue = replaceValue;
    if (returnValue + offset < lowerBound) returnValue = lowerBound;
    if (returnValue + offset > upperBound) returnValue = upperBound;
    if (useOffset) returnValue += offset;
    return returnValue;
  }

  const double dxy(const reco::Candidate* dau, const reco::Vertex& v) {

    return dxy(dau->vx(), v.x(), dau->vy(), v.y());
  }

  const double dz(const reco::Candidate* dau, const reco::Vertex& v) {
 
    return dz(dau->vz(), v.z());
  }

  const double d3d(const reco::Candidate* dau, const reco::Vertex& v) {

    double dxyval = dxy(dau->vx(), v.x(), dau->vy(), v.y());
    double dzval = dz(dau->vz(), v.z());
    return d3d(dxyval, dzval);
  }

  const double dxyErr(const reco::Candidate* dau, const reco::Vertex& v) {

    // 0.0 error for generated daughter
    return dxyErr(dau->vx(), v.x(), 0.0, v.xError(), dau->vy(), v.y(), 0.0, v.yError());
  }

  const double dzErr(const reco::Candidate* dau, const reco::Vertex& v) {

    // 0.0 error for generated daughter
    return dzErr(0.0, v.zError());
  }

  const double d3dErr(const reco::Candidate* dau, const reco::Vertex& v) {

    // 0.0 error for generated daughter
    double dxyval = dxy(dau->vx(), v.x(), dau->vy(), v.y());
    double dxyerr = dxyErr(dau->vx(), v.x(), 0.0, v.xError(), dau->vy(), v.y(), 0.0, v.yError());
    double dzval = dz(dau->vz(), v.z());
    double dzerr = dzErr(0.0, v.zError());
    return d3dErr(dxyval, dxyerr, dzval, dzerr);
  }

  const double d3d(const reco::TrackBaseRef& trkRef) {

    return d3d(trkRef->dxy(), trkRef->dz());
  }

  const double d3dErr(const reco::TrackBaseRef& trkRef) {

    return d3dErr(trkRef->dxy(), trkRef->dxyError(), trkRef->dz(), trkRef->dzError());
  }

  Measurement1D dxy(const reco::Vertex& v1, const reco::Vertex& v2) {

    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix cm1;
    v1.fill(cm1);
    reco::Vertex v1new(v1.position(), cm1);
    return dist.distance(v1new, v2);
  }

  Measurement1D d3d(const reco::Vertex& v1, const reco::Vertex& v2) {

    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix cm1;
    v1.fill(cm1);
    reco::Vertex v1new(v1.position(), cm1);
    return dist.distance(v1new, v2);
  }

  const double dz(const reco::Vertex& v1, const reco::Vertex& v2) {

    return dz(v1.z(), v2.z());
  }

  const double dzErr(const reco::Vertex& v1, const reco::Vertex& v2) {

    return dzErr(v1.zError(), v2.zError());
  }

  const double dxy(const double x1, const double x2, const double y1, const double y2) {

    float dx = abs(x1 - x2);
    float dy = abs(y1 - y2);
    return TMath::Sqrt(dx * dx + dy * dy);
  }

  const double dz(const double z1, const double z2) {

    return abs(z1 - z2);
  }

  const double d3d(const double dxy, const double dz) {

    return TMath::Sqrt(dxy * dxy + dz * dz);
  }

  const double dxyErr(const double x1, const double x2, const double xerr1, const double xerr2,
      const double y1, const double y2, const double yerr1, const double yerr2) {

    float dx = abs(x1 - x2);
    float dy = abs(y1 - y2);
    float dxerr = TMath::Sqrt(xerr1 * xerr1 + xerr2 * xerr2);
    float dyerr = TMath::Sqrt(yerr1 * yerr1 + yerr2 * yerr2);
    float dx2err = 2 * dx * dxerr;
    float dy2err = 2 * dy * dyerr;
    float dxy2err = TMath::Sqrt(dx2err * dx2err + dy2err * dy2err);
    return 0.5 * dxy2err / dxy(x1, x2, y1, y2);
  }

  const double dzErr(const double zerr1, const double zerr2) {

    return TMath::Sqrt(zerr1 * zerr1 + zerr2 * zerr2);
  }

  const double d3dErr(const double dxy, const double dxyerr, const double dz, const double dzerr) {

    float dxy2err = 2 * dxy * dxyerr;
    float dz2err = 2 * dz * dzerr;
    float d3d2err = TMath::Sqrt(dxy2err * dxy2err + dz2err * dz2err);
    return 0.5 * d3d2err / d3d(dxy, dz);
  }
}
