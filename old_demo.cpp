/**
 *  @brief Example of using pythia8185 with HepMC
 *
 *  Usage ./pythia8-example.exe <config_file> <events_count>
 *
 *  The selected process and its options are listed in config files.
 *  First events will automatically be printed by pythia on screen
 *
 *  @note See pythia8185.options.log for a full list of options
 */

#include <iostream>

// pythia8 header files
#include <Pythia8/Pythia.h>
#include <Pythia8/Pythia8ToHepMC.h>

#include "src/AcerDetParticle.h"
#include "src/AcerDetInputReader.cxx"

using std::cout;
using std::endl;
using namespace HepMC;

int main( int argc, char **argv ) {

	if( argc < 3 ) {
		cout << "Usage: " << argv[0] << " <config_file> <events_count> " << endl;
		return -1;
	}

	int events_limit = atoi( argv[2] );

	// -----------------------
	// -- Initialize pythia --
	// -----------------------
	Pythia8::Pythia pythia;
	Pythia8::Event& event = pythia.event;
	Pythia8ToHepMC toHepMC;

	// Use this line to output all available options and their default value
	//pythia.settings.listAll(); exit(0);

	pythia.readFile(argv[1]);
	pythia.init();

	AcerDetInputReader reader;

	int taupair_count = 0;
	int print_count   = 3;

	// ----------------------
	// -- Begin event loop --
	// ----------------------
	for( int iEvent = 0; iEvent < events_limit; ++iEvent ) {
	
		// -------------------------
		// -- Generate next event --
		// -------------------------
		if( !pythia.next() ) continue;

		GenEvent hepmc;
		toHepMC.fill_next_event( event, &hepmc );

		// -----------------------
		// -- Analyze the event --
		// -----------------------
		GenParticle *taup = NULL, *taum = NULL;
		for( GenEvent::particle_iterator i  = hepmc.particles_begin(); i != hepmc.particles_end(); ++i ) {
			if( abs((*i)->pdg_id()) != 15 || (*i)->status() > 2 ) continue;

			// Take only taus that have Z or Z' as their mother
			GenParticle *first_mother = *((*i)->production_vertex()->particles_in_const_begin());

			if( first_mother->pdg_id() == 23 || first_mother->pdg_id() == 32 ) {
				// Store their pointers
				if( (*i)->pdg_id() == -15) taup = *i;
				else                       taum = *i;
			}
		}

		if( taup && taum ) {

			++taupair_count;

			GenParticle *mother = *(taup->production_vertex()->particles_in_const_begin());

			AcerDetParticle mZ( mother->momentum().px(), mother->momentum().py(), mother->momentum().pz(), mother->momentum().e(), mother->pdg_id() );
			AcerDetParticle mp(   taup->momentum().px(),   taup->momentum().py(),   taup->momentum().pz(),   taup->momentum().e(),   taup->pdg_id() );
			AcerDetParticle mm(   taum->momentum().px(),   taum->momentum().py(),   taum->momentum().pz(),   taum->momentum().e(),   taum->pdg_id() );

			if( print_count > 0 ) {
				cout << endl << "BEFORE ANALYSIS: " << endl;
				mZ.print();
				mp.print();
				mm.print();
			}

			// 1) boost tau pair to rest frame of their mother

			mp.boostToRestFrame(mZ);
			mm.boostToRestFrame(mZ);

			// 2) Rotate tau pair to frame where taus are along Z axis
			//    (order of these lines is important)

			double phi = mp.getAnglePhi();

			mp.rotateXY( -phi );
			mm.rotateXY( -phi );

			double theta = mp.getAngleTheta();

			mp.rotateXZ( M_PI-theta );
			mm.rotateXZ( M_PI-theta );

			// Now taus are back-to-back on Z axis

			if( print_count > 0 ) {
				cout << "IN REST FRAME OF TAU PAIR: " << endl;
				mZ.print();
				mp.print();
				mm.print();
			}

			// Back to laboratory frame
			mm.rotateXZ( theta - M_PI );
			mp.rotateXZ( theta - M_PI );
			mm.rotateXY( phi   );
			mp.rotateXY( phi   );
			mm.boostFromRestFrame(mZ);
			mp.boostFromRestFrame(mZ);

			if( print_count > 0 ) {
				cout << "AFTER ANALYSIS: " << endl;
				mZ.print();
				mp.print();
				mm.print();
			}
			if( print_count > 0 ) --print_count;
		}
	}
	
	pythia.statistics();
	cout << "Tau pairs analyzed: " << taupair_count << endl;

	return 0;
}
