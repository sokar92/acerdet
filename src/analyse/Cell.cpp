#include "Cell.h"
#include "../core/Typedefs.h"
#include "../core/Functions.h"
#include <cstdio>

using namespace AcerDet::analyse;
using namespace AcerDet::io;

Cell::Cell( const Configuration& config ) :
	ETACEL	( config.Cell.RapidityCoverage ),
	PTMIN	( config.Cell.MinpT ),
	ETTHR	( config.Cell.MinEt ),
	CALOTH	( config.Cell.EtaTransition ),
	DBETA	( config.Cell.GranularityEta ),
	DBPHI	( config.Cell.GranularityPhi ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histo ("Cell: multiplicity", 0.0, 500.0, 50)
{}

Cell::~Cell() {}

void Cell::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::Cell   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ...\n");
	printf (" eta coverage  %lf\n", ETACEL);
	printf (" E_T_min cell thresh %lf\n", ETTHR);
	printf (" eta gran. transition %lf\n", CALOTH);
	printf (" gran in eta(central) %lf\n", DBETA);
	printf (" gran in phi(central) %lf\n", DBPHI);
	printf ("\n\t B field apply ....\n");
	printf (" B-field %s\n", KEYFLD ? "on" : "off");
	if (KEYFLD) {
		printf (" p_T min non looping %lf\n", PTMIN);
	}
	printf ("\n");
}

void Cell::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {

	// new event to compute
	IEVENT++;
	
	// corrections to much detector granularity
	Int32_t NBETA = round(2.0 * ETACEL / DBETA);
	Int32_t NBPHI = round(6.4 / DBPHI);
	Real64_t PTLRAT = 1.0 / pow(sinh(ETACEL), 2.0);
	
	// reference to particle container
	const vector<Particle>& parts = irecord.particles();
	
	// temporary cells container
	vector<CellData> tempCells;
	
	// Loop over all particles.
	// Find cell that was hit by given particle.
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (!part.isStable())
			continue;
			
		Real64_t DETPHI = 0.0;
		Real64_t ETA, PHI, PT, PZ;

		PT = part.pT();
		PZ = part.pZ();
		if (PT * PT <= PTLRAT * PZ * PZ)
			continue;

		if (part.type == PT_UNKNOWN || part.isNeutrino() || part.type == PT_MUON || part.type == KFINVS)
			continue;

		if (KEYFLD && part.getKuchge() != 0) {
			if (part.pT() < PTMIN)
				continue;
				
			Real64_t CHRG = part.getKuchge() / 3.0;
			DETPHI = CHRG * part.foldPhi();
		}
		
		PT = part.pT();
		ETA = part.getEta();
		PHI = saturatePi(part.getPhi() + DETPHI);
		
		Int32_t IETA, IPHI;
		if (abs(ETA) < CALOTH) {
			IETA = 1 + static_cast<Int32_t>( (ETA + ETACEL) / 2.0 / ETACEL * NBETA );
			IPHI = 1 + static_cast<Int32_t>( (PHI + PI) / 2.0 / PI * NBPHI );
		} else {
			IETA = 1 + 2 * static_cast<Int32_t>( (ETA + ETACEL) / 2.0 / ETACEL * NBETA / 2.0 );
			IPHI = 1 + 2 * static_cast<Int32_t>( (PHI + PI) / 2.0 / PI * NBPHI / 2.0 );
		}
		
		Int32_t cellID = NBPHI * IETA + IPHI;
		
		// Add to cell already hit
		Bool_t found = false;
		for (int j=0; j<tempCells.size(); ++j) {
			if (cellID == tempCells[j].cellID) {
				tempCells[j].hits++;		// new part hits this cell
				tempCells[j].pT += PT;		// summing pT of part hits 
				found = true;
				break;
			}
		}
		
		// Or book new cell
		if (!found) {
			CellData newCell;
			newCell.cellID = cellID;			// not used ID
			newCell.hits = 1;				// only single hit for now
			newCell.state = 2;				// type?
			newCell.pT = PT;				// pT from single hit
			
			if (abs(ETA) < CALOTH) {
				newCell.eta = 2.0 * ETACEL * (IETA - 1.0 + 0.5) / NBETA - ETACEL;
				newCell.phi = 2.0 * PI * (IPHI - 1.0 + 0.5) / NBPHI - PI;
			} else {
				newCell.eta = 2.0 * ETACEL * (IETA - 1.0 + 1.0) / NBETA - ETACEL;
				newCell.phi = 2.0 * PI * (IPHI - 1.0 + 1.0) / NBPHI - PI;
			}

			tempCells.push_back(newCell);
		}
	}

	// Remove cells below threshold and store cells-map in output record
	for (int j=0; j<tempCells.size(); ++j) {
		// enough pT to create new cell
		if (tempCells[j].pT > ETTHR) {
			orecord.Cells.push_back(tempCells[j]);
		}
	}

	// call histogram
	histo.insert( orecord.Cells.size() );
}

void Cell::printResults() const {
	printf ("***********************************\n");
	printf ("*                                 *\n");
	printf ("*     ***********************     *\n");
	printf ("*     ***   Output from   ***     *\n");
	printf ("*     ***  analyse::Cell  ***     *\n");
	printf ("*     ***********************     *\n");
	printf ("*                                 *\n");
	printf ("***********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo.print( true );
}
