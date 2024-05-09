#include "../interface/VertexCalculatorExtended.h"


namespace vertexntuples {

  const double xres(const GenVertex& gv, const SecondaryVertex& sv) {

    return gv.x() - sv.x();
  }

  const double yres(const GenVertex& gv, const SecondaryVertex& sv) {

    return gv.y() - sv.y();
  }

  const double zres(const GenVertex& gv, const SecondaryVertex& sv) {

    return gv.z() - sv.z();
  }

  const double xpull(const GenVertex& gv, const SecondaryVertex& sv) {

    return (gv.x() - sv.x()) / sv.xErr();
  }

  const double ypull(const GenVertex& gv, const SecondaryVertex& sv) {

    return (gv.y() - sv.y()) / sv.yErr();
  }

  const double zpull(const GenVertex& gv, const SecondaryVertex& sv) {

    return (gv.z() - sv.z()) / sv.zErr();
  }

  const double dxy(const GenVertex& gv, const SecondaryVertex& sv) {

    return dxy(gv.x(), sv.x(), gv.y(), sv.y());
  }

  const double dz(const GenVertex& gv, const SecondaryVertex& sv) {

    return dz(gv.z(), sv.z());
  }

  const double d3d(const GenVertex& gv, const SecondaryVertex& sv) {

    double dxyval = dxy(gv.x(), sv.x(), gv.y(), sv.y());
    double dzval = dz(gv.z(), sv.z());
    return d3d(dxyval, dzval);
  }

  const double dxyErr(const GenVertex& gv, const SecondaryVertex& sv) {

    return dxyErr(gv.x(), sv.x(), gv.xErr(), sv.xErr(),
        gv.y(), sv.y(), gv.yErr(), sv.yErr());
  }

  const double dzErr(const GenVertex& gv, const SecondaryVertex& sv) {

    return dzErr(gv.zErr(), sv.zErr());
  }

  const double d3dErr(const GenVertex& gv, const SecondaryVertex& sv) {

    double dxyval = dxy(gv.x(), sv.x(), gv.y(), sv.y());
    double dxyerr = dxyErr(gv.x(), sv.x(), gv.xErr(), sv.xErr(),
        gv.y(), sv.y(), gv.yErr(), sv.yErr());
    double dzval = dz(gv.z(), sv.z());
    double dzerr = dzErr(gv.zErr(), sv.zErr());
    return d3dErr(dxyval, dxyerr, dzval, dzerr);
  }
}
