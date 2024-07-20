#ifndef VertexNtuples_VertexNtuplizer_VertexCalculatorExtended_h
#define VertexNtuples_VertexNtuplizer_VertexCalculatorExtended_h


#include "VertexCalculator.h"
#include "GenVertex.h"
#include "SecondaryVertex.h"


// Made this extended file to solve circular dependency issues
namespace vertexntuples {

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
