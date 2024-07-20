#include "../interface/VertexCalculatorExtended.h"


namespace vertexntuples {

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
