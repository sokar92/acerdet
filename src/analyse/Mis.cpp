#include "Mis.h"
#include <cstdio>

using namespace AcerDet::analyse;

Mis::Mis( const Configuration& config ) :
	PTMUMIN	( config.Muon.MinMomenta ),
	ETAMAX	( config.Muon.MaxEta ),
	ETCELL	( config.Tau.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),

	histo_reconstructed_pT			("Mis: reconstructed p_T", 0.0, 200.0, 50),
	histo_reconstructed_pT_cells	("Mis: reconstructed p_T + cells", 0.0, 200.0, 50),
	histo_pTmiss					("Mis: pTmiss", 0.0, 200.0, 50),
	histo_pTnu						("Mis: pT nu", 0.0, 200.0, 50)
{}

Mis::~Mis() {}

void Mis::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Mis   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t muon coverage\n");
	printf (" min. muon p_T %lf\n", PTMUMIN);
	printf (" max. muon eta %lf\n", ETAMAX);
	printf ("\t unused cells ...\n");
	printf (" smearing %s\n", KEYSME ? "on" : "off");
	printf (" cells threshold %lf\n", ETCELL);
	printf ("\t invisible particles ...\n");
	printf (" KF code for invis %d\n", KFINVS);
	printf ("\n");
}

void Mis::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	
    // new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	// znajdz poczatek danych
	Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].stateID != 21) {
			NSTOP = i-1;
			NSTART = i;
			break;
		}
	}
	
	//sum up reconstructed momenta
	Real64_t PXREC = 0.0, PYREC = 0.0;
	Real64_t PXSUM = 0.0, PYSUM = 0.0;
	Real64_t PXXCALO = 0.0, PYYCALO = 0.0;
	Real64_t SUMET = 0.0;
	
	// add jets
	for (int i=0; i<orecord.Jets.size(); ++i) {
		const JetData jet = orecord.Jets[i];
		PXREC	+= jet.pT * cos( jet.phi_rec );
		PYREC	+= jet.pT * sin( jet.phi_rec );
		PXXCALO	+= jet.pT * cos( jet.phi_rec );
		PYYCALO	+= jet.pT * sin( jet.phi_rec );
		SUMET	+= jet.pT;
	}
	
	// add non-used clusters
	for (int i=0; i<orecord.Clusters.size(); ++i) {
		const ClusterData& cluster = orecord.Clusters[i];
		// if (KCLU(I,5) != 0) { TODO remove unused from vector before
			PXREC	+= cluster.pT * cos( cluster.phi_rec );
			PYREC	+= cluster.pT * sin( cluster.phi_rec );
			PXXCALO	+= cluster.pT * cos( cluster.phi_rec );
			PYYCALO	+= cluster.pT * sin( cluster.phi_rec );
			SUMET	+= cluster.pT;
		// }
	}
	
	// add isolated muons
	for (int i=0; i<orecord.Muons.size(); ++i) {
		const PartData& muon = orecord.Muons[i];
		PXREC += muon.pT * cos( muon.phi ); //_rec );
		PYREC += muon.pT * sin( muon.phi ); //_rec );
	}
	
	// add non-isolated muons not added to clusters
	for (int i=0; i<orecord.NonisolatedMuons.size(); ++i) {
		const PartData& muon = orecord.NonisolatedMuons[i];
		/* TODO mark used before
		if (KMUOX(I,5) != 0) {
			PXREC += muon.pT * cos( muon.phi ); //_rec );
			PYREC += muon.pT * sin( muon.phi ); //_rec );
		} else {
			SUMET -= muon.pT;
		}*/
	}
	
	// add isolated electrons
	for (int i=0; i<orecord.Electrons.size(); ++i) {
		const PartData& ele = orecord.Electrons[i];
		PXREC	+= ele.pT * cos( ele.phi ); //_rec );
		PYREC	+= ele.pT * sin( ele.phi ); //_rec );
		PXXCALO	+= ele.pT * cos( ele.phi ); //_rec );
		PYYCALO	+= ele.pT * sin( ele.phi ); //_rec );
		SUMET	+= ele.pT;
	}
	
	// add isolated photons
	for (int i=0; i<orecord.Photons.size(); ++i) {
		const PartData& pho = orecord.Photons[i];
		PXREC	+= pho.pT * cos( pho.phi ); //_rec );
		PYREC	+= pho.pT * sin( pho.phi ); //_rec );
		PXXCALO	+= pho.pT * cos( pho.phi ); //_rec );
		PYYCALO	+= pho.pT * sin( pho.phi ); //_rec );
		SUMET	+= pho.pT;
    }
    
    // store pT in histo
    Real64_t ETREC = sqrt( pow(PXREC, 2) + pow(PYREC, 2) );
    histo_reconstructed_pT.insert( ETREC );
    
    // smear cells energy not used for reconstruction
    // remove cells below threshold
    // add momenta in cells not used for reconstruction
    for (int i=0; i<orecord.Cells.size(); ++i) {
		const CellData& cell = orecord.Cells[i];
		/* TODO
		if (cell.state != 0) {
			if (KEYSME) {
				PEI = cell.pT * cosh( cell.eta );
				SIGSME = RESHAD(PEI, cell.eta, CALOTH, cell.pT, 0.0);
				
				PSME = PEI * (1.0 + SIGSME) / cosh( cell.eta );
				if (PSME < ETCELL) 
					PSME = 0.0;
					
				PXSUM	+= PSME * cos( cell.phi );
				PYSUM	+= PSME * sin( cell.phi );
				PXXCALO	+= PSME * cos( cell.phi );
				PYYCALO	+= PSME * sin( cell.phi );
				SUMET	+= PSME;
			}
		}*/
	}
	
	PXSUM += PXREC;
	PYSUM += PYREC;
	
	Real64_t ETSUM = sqrt( pow(PXSUM, 2) + pow(PYSUM, 2) );
	histo_reconstructed_pT_cells.insert( ETSUM );
	
	Real64_t PXXMISS = -PXSUM;
	Real64_t PYYMISS = -PYSUM; // ? po co
	Real64_t PTMISS = sqrt( pow(PXXMISS, 2) + pow(PYYMISS, 2) );
	histo_pTmiss.insert( PTMISS );
	
	// sum up momenta  of neutrinos 
	Real64_t PXXNUES = 0.0;
	Real64_t PYYNUES = 0.0;
	for (int i=NSTART; i<parts.size(); ++i) {   
		if (parts[i].isNeutrino() || abs(parts[i].typeID) == KFINVS) {
			PXXNUES += parts[i].pX();
			PYYNUES += parts[i].pY();
		}
	}
	
	Real64_t PTNUES = sqrt( pow(PXXNUES, 2) + pow(PYYNUES, 2) );
	histo_pTnu.insert( PTNUES );
}

void Mis::printResults() const {
	printf ("***********************************\n");
	printf ("*                                 *\n");
	printf ("*     ***********************     *\n");
	printf ("*     ***   Output from   ***     *\n");
	printf ("*     ***  analyse::Miss  ***     *\n");
	printf ("*     ***********************     *\n");
	printf ("*                                 *\n");
	printf ("***********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo_reconstructed_pT			.print( true );
	histo_reconstructed_pT_cells	.print( true );
	histo_pTmiss					.print( true );
	histo_pTnu						.print( true );
}
