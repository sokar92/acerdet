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
	Real64_t NBETA = round(2.0 * ETACEL / DBETA);
	Real64_t NBPHI = round(6.4 / DBPHI);
	
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

		if (part.pT() <= PTLRAT * part.pZ())
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
	
	/*

      INTEGER NBETA, NBPHI
      REAL DBETA, DBPHI, CALOTH
      REAL ETACEL, ETTHR
      REAL ETA,PT,PHI
      REAL PTMIN
      INTEGER IPHI,IETPH,IC,KC,IETA
      REAL THETA,PTLRAT,DETPHI,CHRG
      REAL FLDPHI
      INTEGER KUCHGE, KFCOMP
      REAL ANGLE
      INTEGER NDUM, KDUM
      REAL PDUM
      INTEGER NBINA,IDENT
      REAL TMAXA,TMINA

c..... corrections to much detector granularity
->      NBETA = NINT(2 * ETACEL / DBETA)
->      NBPHI = NINT(6.4 / DBPHI)

->      ELSEIF(MODE == 0) THEN

->      IEVENT = IEVENT + 1							// nowy event
->      NDUM = 0

C...Loop over all particles. Find cell that was hit by given particle. 
->	PTLRAT = 1. / SINH(ETACEL)**2 
->      DO I = from 1 to N 
->		IF(K(I,1) <= 0 || K(I,1) > 10) continue 					// SPECJALNA WLASNOCS CZASTEK - jak definiowac?  '->' stable_particles?
->		IF(P(I,1)**2 + P(I,2)**2 <= PTLRAT * P(I,3)**2) continue			// if (px^2 + py^2 <= ptlrat * pz^3) {
->		KC = KFCOMP(K(I,2)) 
->		IF(KC == 0 || KC == 12 || KC == 14 || KC == 16 || KC == 13 || KC == KFINVS) continue 
->		DETPHI = 0.0
->		IF(KEYFLD == true && KUCHGE(K(I,2)) != 0) 
->		THEN
->			PT = SQRT(P(I,1)**2+P(I,2)**2)						// pt = sqrt( px^2 + py^2 )
->			IF(PT < PTMIN) continue
->			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
->			PHI = ANGLE(P(I,1),P(I,2))
->			CHRG = KUCHGE(K(I,2))/3.
->			IF(ABS(P(I,3)/P(I,4)) < 1) 
->			THEN
->				THETA = ACOS(P(I,3)/SQRT(PT**2+P(I,3)**2))				// kat z iloczynu skalarnego (pz,p)
->			ELSE
->				IF(P(I,3) > 0.0) THETA = 0						// normalizacja z wysyceniem do przedzialu [0, pi]
->				IF(P(I,3) < 0.0) THETA = PI
->			ENDIF
->			DETPHI = CHRG * FLDPHI(PT,ETA,THETA)
->		ENDIF 
->		PT = SQRT(P(I,1)**2+P(I,2)**2)						// pt = sqrt( px^2 + py^2 )
->		ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
->		PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py
->		PHI = PHI + DETPHI
->		IF(PHI > PI) 
->		THEN								// normalizacja kata do [-pi,pi]
->			PHI = PHI - 2 * PI
->		ELSEIF(PHI < -PI) 
->		THEN
->			PHI = PHI + 2 * PI
->		ENDIF
->		IF(ABS(ETA) < CALOTH) THEN						// jest do tego funkcja
->			IETA = INT( (ETA+ETACEL) / 2.0 / ETACEL * NBETA) + 1 
->			IPHI = INT( (PHI+PI) / 2.0 / PI * NBPHI ) + 1 
->		ELSE
->			IETA = 2 * INT( (ETA+ETACEL) / 2.0 / ETACEL * NBETA / 2.0 ) + 1 
->			IPHI = 2 * INT( (PHI+PI) / 2.0 / PI * NBPHI / 2.0 ) + 1 
->		ENDIF
->		IETPH = NBPHI * IETA + IPHI

 
C...Add to cell already hit, or book new cell.
->		DO IC = from 1 to NDUM { 
->			IF(IETPH == KDUM(IC,3)) {
->				KDUM(IC,4) = KDUM(IC,4) + 1 
->				PDUM(IC,5) = PDUM(IC,5) + PT 
->				continue 
->     			}
->		} 

->		NDUM = NDUM + 1 
->		KDUM(NDUM,3) = IETPH 
->		KDUM(NDUM,4) = 1 
->		KDUM(NDUM,5) = 2 

->		IF(ABS(ETA) < CALOTH) 
->		THEN						// caloth z configa
->			PDUM(NDUM,1) = 2.0 * ETACEL*(FLOAT(IETA)-1.0+0.5)/NBETA-ETACEL
->			PDUM(NDUM,2) = 2.0 * PI*(FLOAT(IPHI)-1.0+.5)/NBPHI-PI
->		ELSE
->			PDUM(NDUM,1) = 2.0 * ETACEL*(FLOAT(IETA)-1.0+1.0)/NBETA-ETACEL
->			PDUM(NDUM,2) = 2.0 * PI*(FLOAT(IPHI)-1.0+1.0)/NBPHI-PI
->		ENDIF
->		PDUM(NDUM,5) = PT  
->	ENDDO

C...Remove cells below threshold  and store cells-map in common /CELLS/ 
->	NCELL = 0
->	DO IC = from 1 to NDUM { 
->		IF(PDUM(IC,5) > ETTHR) 
->		THEN						// spelniony warunek -> dodaj cellke do outrecordu 
->			NCELL = NCELL + 1 
->			KCELL(NCELL,3) = KDUM(IC,3) 
->			KCELL(NCELL,4) = KDUM(IC,4) 
->			KCELL(NCELL,5) = KDUM(IC,5) 
->			PCELL(NCELL,1) = PDUM(IC,1) 
->			PCELL(NCELL,2) = PDUM(IC,2) 
->			PCELL(NCELL,5) = PDUM(IC,5)
->		ENDIF 
->	}

c.....histogram NCELL
      CALL HF1(IDENT+1, REAL(NCELL), 1.0)						// i do histogramu
C>>>>

C
      ELSEIF(MODE == 1) THEN
C     =========================

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM AcerDET      '
      WRITE(NOUT,BXTXT) '              ACDCEL             '
      WRITE(NOUT,BXTXT) '*********************************'
      
      WRITE(NOUT,BXCLO)
      
      CALL HPRINT(IDENT+1)

      ENDIF
C     =====

      END
*/

}
