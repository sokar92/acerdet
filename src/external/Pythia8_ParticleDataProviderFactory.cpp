#include "Pythia8_ParticleDataProviderFactory.h"
using namespace AcerDet::external;

// pythia8 header files
#include <Pythia8/Pythia.h>
using namespace Pythia8;

#define PARTS_SCAN 4000000

ParticleDataProvider Pythia8_ParticleDataProviderFactory::create() {
	ParticleDataProvider newProvider;
	Pythia pythia;
	pythia.init( 2212, 2212, 14000.0); // must be properly initialized to access particle data
	
	for (int j=-PARTS_SCAN; j<=PARTS_SCAN; ++j) {
		const string& name = pythia.particleData.name(j);
		Int32_t cht  = pythia.particleData.chargeType(j);
		Real64_t ch = pythia.particleData.charge(j);
		
		if (name.size() > 0 && name[0] != ' ')
			newProvider.insertInfo(j, cht, ch, name); 
	}
	
	return newProvider;
}
