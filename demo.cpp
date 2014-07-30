#include <iostream>

// pythia8 header files
#include <Pythia8/Pythia.h>
#include <Pythia8/Pythia8ToHepMC.h>

using std::cout;
using std::endl;
using namespace HepMC;

#include <AcerDET/AcerDET.h>
using namespace AcerDet;
using namespace AcerDet::conf;
using namespace AcerDet::core;
using namespace AcerDet::io;

// ----------------------
// -- Dummy Database ----
// ----------------------
class DbDummy {
public:
	void store( const io::OutputRecord& rec ) const {
		printf ("DB: storing ...\n");
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
	AcerDET acerDet( configuration );
	acerDet.printInfo();


	Histogram histo("testing", 2.0, 12.5, 10);
	histo.insert(2.0);
	histo.insert(2.1);
	histo.insert(5.0);
	histo.insert(4.8);
	histo.insert(9.0);
	histo.insert(12.4);
	histo.insert(12.5);
	histo.insert(-2.5);
	histo.insert(12.5001);
	histo.print( false );
	
	printf (" proper -> %d\n", histo.storedProperlyCount());
	printf (" all -> %d\n", histo.storedCount());
	
	DbDummy db;
	OutputRecord oRec;
	
	events_limit = 2;
	for( int iEvent = 0; iEvent < events_limit; ++iEvent ) {
		printf("Event %d\n", iEvent);
		
		if( !pythia.next() ) {
			printf("continue ...\n");
			continue;
		}
		
		resetRecord(oRec);

		GenEvent *hepmc = new GenEvent();
		toHepMC.fill_next_event( event, hepmc );

		acerDet.analyseRecord( InputReader::computeEvent( *hepmc ), oRec );

		db.store(oRec);
		delete hepmc;
	}

	acerDet.printResults();

	return 0;
}
