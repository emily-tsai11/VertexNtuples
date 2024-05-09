#ifndef VertexNtuples_VertexNtuplizer_VertexCalculatorExtended_h
#define VertexNtuples_VertexNtuplizer_VertexCalculatorExtended_h


#include "VertexCalculator.h"
#include "GenVertex.h"
#include "SecondaryVertex.h"


// Made this extended file to solve circular dependency issues
namespace vertexntuples {

  const double xres(const GenVertex& gv, const SecondaryVertex& sv);
  const double yres(const GenVertex& gv, const SecondaryVertex& sv);
  const double zres(const GenVertex& gv, const SecondaryVertex& sv);
  const double xpull(const GenVertex& gv, const SecondaryVertex& sv);
  const double ypull(const GenVertex& gv, const SecondaryVertex& sv);
  const double zpull(const GenVertex& gv, const SecondaryVertex& sv);

  const double dxy(const GenVertex& gv, const SecondaryVertex& sv);
  const double dz(const GenVertex& gv, const SecondaryVertex& sv);
  const double d3d(const GenVertex& gv, const SecondaryVertex& sv);
  const double dxyErr(const GenVertex& gv, const SecondaryVertex& sv);
  const double dzErr(const GenVertex& gv, const SecondaryVertex& sv);
  const double d3dErr(const GenVertex& gv, const SecondaryVertex& sv);
}


#endif
