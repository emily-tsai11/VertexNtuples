#ifndef SECONDARY_VERTEX
#define SECONDARY_VERTEX


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/Vertex.h"


class SecondaryVertex {

  public:

    SecondaryVertex(const reco::Vertex& sv);
    // ~SecondaryVertex();

  private:

    //
};


#endif
