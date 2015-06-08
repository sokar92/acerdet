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

#include "external/HepMC_InputConverter.h"
#include "external/Pythia8_ParticleDataProviderFactory.h"
#include "external/Root_HistogramManager.h"
#include "external/Root_NTupleManager.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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
	external::Root_NTupleManager nTuple;
	
	AcerDET acerDet(
		configuration,
		pdpFactory,
		histoManager 
	);
	nTuple.init();
	acerDet.printInfo();

	DbDummy db;
	OutputRecord oRec;
	
	for( int iEvent = 0; iEvent < events_limit; ++iEvent ) {
//		printf("Event %d\n", iEvent);
		
		if( !pythia.next() ) {
			printf("continue ...\n");
			continue;
		}
		
		resetRecord(oRec);

		GenEvent *hepmc = new GenEvent();
		toHepMC.fill_next_event( event, hepmc );
//		printf (" ---- PRINTING --- \n");
//		hepmc->print();

		InputRecord iRec = external::HepMC_InputConverter::convert( *hepmc );
		acerDet.analyseRecord(iRec, oRec);

		nTuple.fill(iRec, oRec, 1.0, 1); // 1.0 by default - to change, 1 = processId
		
		delete hepmc;
	}

	acerDet.printResults();
	
	char inputname[1000] = "";
	sprintf( inputname, argv[1] );
	
	char outputname[1000] = "";
	strcat( outputname, inputname );
	strcat( outputname, ".root" );
	
	TFile scanMiniTreeFile( outputname, "recreate" );
	scanMiniTreeFile.cd();

	acerDet.storeHistograms();
	nTuple.write();

	delete pdpFactory;
	delete histoManager;
	return 0;
}
