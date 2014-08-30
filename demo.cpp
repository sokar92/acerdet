#include <iostream>
#include <cstring>

// pythia8 header files
#include <Pythia8/Pythia.h>
#include <Pythia8/Pythia8ToHepMC.h>

using std::cout;
using std::endl;
using namespace HepMC;

#include "src/AcerDET.h"
using namespace AcerDet;
using namespace AcerDet::conf;
using namespace AcerDet::core;
using namespace AcerDet::io;

#include "src/include/external/HepMC_InputConverter.h"
#include "src/include/external/Pythia8_ParticleDataProviderFactory.h"
#include "src/include/external/Root_HistogramManager.h"

// ----------------------
// -- Dummy Database ----
// ----------------------
class DbDummy {
public:
	void store( const io::OutputRecord& rec ) const {
//		printf ("DB: storing ...\n");
//		printf ("\tcells: %d\n", (int)rec.Cells.size());
//		printf ("\tclusters: %d\n", (int)rec.Clusters.size());
//		printf ("\n");
	}
};

// -----------------------
// -- clean temp record --
// -----------------------
void resetRecord( io::OutputRecord& rec ) {
	rec.clear();
}

int main( int argc, char **argv ) {

	if( argc < 3 ) {
		cout << "Usage: " << argv[0] << " <config_file> <events_count> " << endl;
		return -1;
	}
	int events_limit = atoi( argv[2] );
	
	//
    // Initialize pythia
    //
    Pythia8::Pythia pythia;
    Pythia8::Event& event = pythia.event;
    Pythia8ToHepMC toHepMC;

    // Use this line to output all available options and their default values
    //pythia.settings.listAll(); exit(0);

    pythia.readFile(argv[1]);
    pythia.init();

    //
    // Initialize AcerDET
    //
	const std::string configFileName = "acerdet.dat";
	Configuration configuration = Configuration::fromFile( configFileName );
	IParticleDataProviderFactory *pdpFactory =
		new external::Pythia8_ParticleDataProviderFactory();
	IHistogramManager *histoManager =
		new external::Root_HistogramManager();
	
	AcerDET acerDet(
		configuration,
		pdpFactory,
		histoManager 
	);
	acerDet.printInfo();

	DbDummy db;
	OutputRecord oRec;
	
	for( int iEvent = 0; iEvent < events_limit; ++iEvent ) {
		printf("Event %d\n", iEvent);
		
		if( !pythia.next() ) {
			printf("continue ...\n");
			continue;
		}
		
		resetRecord(oRec);

		GenEvent *hepmc = new GenEvent();
		toHepMC.fill_next_event( event, hepmc );

		acerDet.analyseRecord(
			external::HepMC_InputConverter::convert( *hepmc ),
			oRec
		);

		db.store(oRec);
		delete hepmc;
	}

	acerDet.printResults();
	acerDet.storeHistograms( "myFile.root" );

	delete pdpFactory;
	delete histoManager;
	return 0;
}
