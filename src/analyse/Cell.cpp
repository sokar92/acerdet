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

	// corrections to much detector granularity
	Int32_t NBETA = round(2.0 * ETACEL / DBETA);
	Int32_t NBPHI = round(6.4 / DBPHI);
	
	// new event to compute
	IEVENT++;
	Int32_t NDUM = 0;
	
	Real64_t PTLRAT = 1.0 / pow(sinh(ETACEL), 2.0);
	
	// Loop over all particles.
	// Find cell that was hit by given particle.
	const vector<Particle>& parts = irecord.particles();
	Int32_t N = parts.size();
	
	for (int i=0;i<N;++i) {
		const Particle& part = parts[i];
		
		if (!part.isStable())
			continue;

		if (part.pT() <= PTLRAT * part.pZ() * part.pZ())
			continue;
			
		Int32_t KC = kfcomp(part.typeID);
		if (KC == 0 || KC == 12 || KC == 13 || KC == 14 || KC == 16 || KC == KFINVS)
			continue;
			
		Real64_t DETPHI = 0.0;
		Real64_t ETA, PHI, THETA, PT;

		if (KEYFLD && kuchge(part.typeID) != 0) {
			if (part.pT() < PTMIN)
				continue;
				
			Real64_t CHRG = kuchge(part.typeID) / 3.0;
			DETPHI = CHRG * part.foldPhi();
		}
		
		PT = part.pT();
		ETA = part.getEta();
		PHI = saturatePi(part.getPhi() + DETPHI);
		
		Int32_t IETA, IPHI;
		if (abs(ETA) < CALOTH) {
			IETA = 1 + (Int32_t)( (ETA + ETACEL) / 2.0 / ETACEL * NBETA );
			IPHI = 1 + (Int32_t)( (PHI + PI) / 2.0 / PI * NBPHI );
		} else {
			IETA = 1 + 2 * (Int32_t)( (ETA + ETACEL) / 2.0 / ETACEL * NBETA / 2.0 );
			IPHI = 1 + 2 * (Int32_t)( (PHI + PI) / 2.0 / PI * NBPHI / 2.0 );
		}
		Int32_t IEPTH = NBPHI * IETA + IPHI;
		
		// Add to cell already hit
		Bool_t found = false;
		for (int j=0; j<orecord.vDum.size(); ++j) {
			if (IEPTH == orecord.vDum[j].K[2]) {
				orecord.vDum[j].K[3]++;				// new part hits this cell
				orecord.vDum[j].P[4] += PT;			// summing pT of part hits 
				found = true;
				break;
			}
		}
		
		// Or book new cell
		if (!found) {
			Dum newDum;
			newDum.K[2] = IEPTH;			// not used ID?
			newDum.K[3] = 1;				// only single hit for now
			newDum.K[4] = 2;				// type?
			newDum.P[4] = PT;				// pT from single hit
			
			if (abs(ETA) < CALOTH) {
				newDum.P[0] = 2.0 * ETACEL * (IETA - 1.0 + 0.5) / NBETA - ETACEL;
				newDum.P[1] = 2.0 * PI * (IPHI - 1.0 + 0.5) / NBPHI - PI;
			} else {
				newDum.P[0] = 2.0 * ETACEL * (IETA - 1.0 + 1.0) / NBETA - ETACEL;
				newDum.P[1] = 2.0 * PI * (IPHI - 1.0 + 1.0) / NBPHI - PI;
			}

			orecord.vDum.push_back(newDum);
		}
	}
	
	// Remove cells below threshold and store cells-map in output record
	for (int j=0; j<orecord.vDum.size(); ++j) {
		// enough pT to create new cell
		if (orecord.vDum[j].P[4] > ETTHR) {
			orecord.vCell.push_back(orecord.vDum[j]);
		}
	}
	
	// call histogram
	// CALL HF1(IDENT+1, Real64_t(NCELL), 1.0)
}
