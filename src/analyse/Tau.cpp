#include "Tau.h"
#include <cstdio>

using namespace AcerDet::analyse;

Tau::Tau( const Configuration& config, IHistogramManager* histoMng ) :
	ETJET	( config.Jet.MinEnergy ),

	PTTAU	( config.Tau.MinpT ),
	ETATAU	( config.Tau.MaxEta ),
	RJTAU	( config.Tau.MinR ),
	PTFRAC	( config.Tau.MaxR ),

	KEYHID	( config.Flag.HistogramId ),
	KEYTAU	( config.Flag.TauJetsLabeling ),

	IEVENT	( 0 ),
	
	histoManager(histoMng),
	histoRegistered( false )
{}

Tau::~Tau() {}

void Tau::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Tau   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t... tau-jets labeling ...\n");
	printf ("labeling on/off %s\n", KEYTAU ? "on" : "off");
	if (KEYTAU) {
		printf ("\ttau-jets ...\n");
		printf ("min tau-had p_T %lf\n", PTTAU);
		printf ("max tau-had eta %lf\n", ETATAU);
		printf ("max R_tauj for tau-jets %lf\n", RJTAU);
		printf ("tau-had frac. of jet %lf\n", PTFRAC);
	}
	printf ("\n");
}

void Tau::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {


        int idhist = 900 + KEYHID;

	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "Tau: tau-jets multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "Tau: taus multiplicity", 10, 0.0, 10.0);
	}
	
	//histo_jets	("Tau: jets-multiplicity", 0.0, 10.0, 10),
	//histo_taus	("Tau: multiplicity", 0.0, 10.0, 10)



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
	
	// look for tau-jets
	Int32_t NTAU = 0, NJETTAU = 0;
	for (int i=NSTART; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_TAU) {
			// if there are still jets
			if (!orecord.Jets.isEmpty()) {
				Bool_t TAUJET = true;
				
				// choose only hadronic tau-decay
				IN = -1;
				if (!part.hasDaughter())
					continue;
				
				for (int j=part.daughters.first, j<=part.daughters.second; ++j) {
					if (abs(parts[j].typeID) == 11 || abs(parts[j].typeID) == 13)
						TAUJET = false;
					// TODO uzyj enumow do typow
					if (abs(parts[j].typeID) == 16)
						IN = j;
				}
				
				if (IN < 0)
					continue;
				
				Particle tmpPart = part;
				tmpPart.momentum -= parts[IN].momentum;
				
				//PT = sqrt(
				//	pow(part.pX() - parts[IN].pX(), 2) +
				//	pow(part.pY() - parts[IN].pY(), 2)
				//);
				PT = tmpPart.pT(); 
				if (PT < PTTAU) 
					TAUJET = false;
				
				
				//ETA = SIGN(LOG((SQRT(PT**2+(P(I,3)-P(IN,3))**2)+ABS(P(I,3)-P(IN,3)))/PT),P(I,3)-P(IN,3)) 
				ETA = tmpPArt.getEta();
				if (abs(ETA) > ETATAU) 
					TAUJET = false;
					
				//PHI = ANGLE(P(I,1)-P(IN,1),P(I,2)-P(IN,2))
				PHI = tmpPart.getPhi();
				if (TAUJET) 
					NTAU++;
				
				// mark tau-jet
				Real64_t DR = 100.0;
				for (int j=0; j<orecord.Jets.size(); ++j) {
					
					Real64_t DDR = sqrt(
						pow(ETA - PJET(II,3), 2) +
						pow(PHI - PJET(II,4), 2)
					);
					
					if (abs(PHI - PJET(II,4)) > PI)
						DDR = sqrt(
							pow(ETA - PJET(II,3), 2) +
							pow(abs(PHI - PJET(II,4)) - 2*PI, 2)
						);

					if (DDR < DR) {
						JETTAU = j;
						DDR = DR;
					}
				}
				
				if (DR > RJTAU) {
					TAUJET = false;
					JETTAU = 0;
				}
				
				if (TAUJET && abs(PJET(JETTAU,3)) < 2.5 && PT/PJET(JETTAU,5) > PTFRAC) {
					KJET(JETTAU,2) = K(I,2); // ?
					NJETTAU++;
				}
			}
		}
	}

	histoManager
		->insert(idhist+11,  NJETTAU);
	histoManager
		->insert(idhist+21,  NTAU);
	
*/
}

void Tau::printResults() const {
	printf ("**********************************\n");
	printf ("*                                *\n");
	printf ("*     **********************     *\n");
	printf ("*     ***  Output from   ***     *\n");
	printf ("*     ***  analyse::Tau  ***     *\n");
	printf ("*     **********************     *\n");
	printf ("*                                *\n");
	printf ("**********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	//histo_jets.print( true );
	//histo_taus.print( true );
}
