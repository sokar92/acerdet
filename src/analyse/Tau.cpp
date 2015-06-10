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
	
	histoManager( histoMng ),
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

void Tau::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weigth ) {
	Int32_t idhist = 900 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "Tau: tau-jets multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "Tau: hadronic taus multiplicity", 10, 0.0, 10.0);
	}

 	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	
	// look for tau-jets
	Int32_t NTAU = 0, NJETTAU = 0;
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (part.status == PS_DECAYED && part.type == PT_TAU) {
		  // if there are still jets
		  if (orecord.Jets.empty()) 
		    continue;
		  
		  Bool_t TAUJET_tagger = true;
		  
		  // choose only hadronic tau-decay
		  if (!part.hasDaughter())
		    continue;
		  Int32_t IN = -1;
		  for (int j=part.daughters.first; j<=part.daughters.second; ++j) {
		    // daughter is e or m
		    if (parts[j].type == PT_ELECTRON
			|| parts[j].type == PT_MUON)
		      TAUJET_tagger = false;
		    
		    // daughter is tau neutrino
		    if (parts[j].type == PT_NEUTRINO_TAU)
		      IN = j;
		  }
		  
		  // there is no tau neutrino in particle daughters set
		  if (IN < 0)
		    continue;

		  Particle tmpPart = part;
		  tmpPart.momentum -= parts[IN].momentum;  // tau - neutrino_tau
		  
		  if (tmpPart.pT() < PTTAU) 
		    TAUJET_tagger = false;
		  
		  if (abs(tmpPart.getEta()) > ETATAU) 
		    TAUJET_tagger = false;
		  
		  if (TAUJET_tagger) 
		    NTAU++;
		  
		  // mark tau-jet
		  Real64_t DR = 100.0;
		  Int32_t JETTAU = -1;
		  for (int j=0; j<orecord.Jets.size(); ++j) {
		    
		    Real64_t DDR = sqrt(
					pow(tmpPart.getEta() - orecord.Jets[j].eta_rec, 2) +
					pow(tmpPart.getPhi() - orecord.Jets[j].phi_rec, 2)
					);
		    
		    if (abs(tmpPart.getPhi() - orecord.Jets[j].phi_rec) > PI)
		      DDR = sqrt(
				 pow(tmpPart.getEta() - orecord.Jets[j].eta_rec, 2) +
				 pow(abs(tmpPart.getPhi() - orecord.Jets[j].phi_rec) - 2*PI, 2)
				 );
		    
		    if (DDR < DR) {
		      JETTAU = j;
		      DR = DDR;
		    }
		  }
		  
		  if (DR > RJTAU) {
		    TAUJET_tagger = false;
		    JETTAU = -1;
		  }
		  
		  if (TAUJET_tagger  && abs(orecord.Jets[JETTAU].eta_rec) < 2.5 ) {
		    orecord.Jets[JETTAU].type = TAU_JET;
		  }
		}
	}

	//count tau-jets
	for (int j=0; j<orecord.Jets.size(); ++j) {
          if( orecord.Jets[j].type == TAU_JET) NJETTAU++;
	}


	histoManager
		->insert(idhist+11, NJETTAU);
	histoManager
		->insert(idhist+21, NTAU);
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
}
