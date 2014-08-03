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
/*
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
	Real64_t SUMET = 0.0;
	
	// add jets
	for (int i=0; i<orecord.Jets.size(); ++i) {
		PXREC += PJET(I,5) * COS(PJET(I,4));
		PYREC += PJET(I,5) * SIN(PJET(I,4));
		PXXCALO += PJET(I,5) * COS(PJET(I,4));
		PYYCALO += PJET(I,5) * SIN(PJET(I,4));
		SUMET += PJET(I,5);
	}
	
	// add non-used clusters
	for (int i=0; i<orecord.Clusters.size(); ++i) {
		if (KCLU(I,5) != 0) {
			PXREC += PCLU(I,5) * COS(PCLU(I,4));
			PYREC += PCLU(I,5) * SIN(PCLU(I,4));
			PXXCALO += PCLU(I,5) * COS(PCLU(I,4));
			PYYCALO += PCLU(I,5) * SIN(PCLU(I,4));
			SUMET += PCLU(I,5);
		}
	}
	
	// add isolated muons
	for (int i=0; i<NMUO; ++i) {
		PXREC += PMUO(I,5) * COS(PMUO(I,4));
		PYREC += PMUO(I,5) * SIN(PMUO(I,4));
	}
	
	// add non-isolated muons not added to clusters
	for (int i=0; i<NMUOX; ++i) {
		if (KMUOX(I,5) != 0) {
			PXREC += PMUOX(I,5) * COS(PMUOX(I,4));
			PYREC += PMUOX(I,5) * SIN(PMUOX(I,4));
		} else {
			SUMET -= PMUOX(I,5);
		}
	}
	
	// add isolated electrons
	for (int i=0; i<NELE; ++i) {
		PXREC += PELE(I,5) * COS(PELE(I,4));
		PYREC += PELE(I,5) * SIN(PELE(I,4));
		PXXCALO += PELE(I,5) * COS(PELE(I,4));
		PYYCALO += PELE(I,5) * SIN(PELE(I,4));
		SUMET += PELE(I,5);
	}
	
	// add isolated photons
	for (int i=0; i<NPHO; ++i) {
		PXREC += PPHO(I,5) * COS(PPHO(I,4));
		PYREC += PPHO(I,5) * SIN(PPHO(I,4));
		PXXCALO += PPHO(I,5) * COS(PPHO(I,4));
		PYYCALO += PPHO(I,5) * SIN(PPHO(I,4));
		SUMET += PPHO(I,5);
    }
    
    Real64_t ETREC = sqrt( pow(PXREC, 2) + pow(PYREC, 2) );
    // CALL HF1(IDENT+11,ETREC,1.0)
    
    // smear cells energy not used for reconstruction
    // remove cells below threshold
    // add momenta in cells not used for reconstruction
    for (int i=0; i<orecord.Cells.size(); ++i) {
		const CellData& cell = orecord.Cells[i];
		
		if (cell.state != 0) {
			if (KEYSME) {
				PEI = PCELL(I,5) * COSH(PCELL(I,1));
				SIGSME = RESHAD(PEI,PCELL(I,1),CALOTH,PCELL(I,5), 0.0);
				
				PSME = PEI * (1.0 + SIGSME) / COSH(PCELL(I,1));
				if (PSME < ETCELL) 
					PSME = 0.0;
					
				PXSUM += PSME * COS(PCELL(I,2));
				PYSUM += PSME * SIN(PCELL(I,2));
				PXXCALO += PSME * COS(PCELL(I,2));
				PYYCALO += PSME * SIN(PCELL(I,2));
				SUMET += PSME;
			}
		}
	}
	
	PXSUM += PXREC;
	PYSUM += PYREC;
	
	ETSUM = sqrt( pow(PXSUM, 2) + pow(PYSUM, 2) );
	//CALL HF1(IDENT+12,ETSUM,1.0)
	
	PXXMISS = -PXSUM;
	PYYMISS = -PYSUM; // ? po co
	PTMISS = sqrt( pow(PXXMISS, 2) + pow(PYYMISS, 2) );
	// CALL HF1(IDENT+13,PTMISS,1.0)
	
	// sum up momenta  of neutrinos 
	PXXNUES = 0.0;
	PYYNUES = 0.0;
	for (int i=NSTART; i<parts.size(); ++i) {   
		if (abs(parts[i].typeID) == 12 || 
			abs(parts[i].typeID) == 14 ||
			abs(parts[i].typeID) == 16 ||
			abs(parts[i].typeID) == KFINVS)
		{
			PXXNUES += parts[i].pX();
			PYYNUES += parts[i].pY();
		}
	}
	
	PTNUES = sqrt( pow(PXXNUES, 2) + pow(PYYNUES, 2) );
	// CALL HF1(IDENT+21,PTNUES,1.0)
*/
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
