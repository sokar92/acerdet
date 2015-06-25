#include "CJet.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Functions.h"
using namespace AcerDet::core;

CJet::CJet( const Configuration& config, IHistogramManager *histoMng ) :
	ETJET	( config.Jet.MinEnergy ),
	RCONE	( config.Cluster.ConeR ),

	PTCMIN	( config.CJet.MinMomenta ),
	ETCMAX	( config.CJet.MaxEta ),
	RJC	( config.CJet.MaxRcj ),

	KEYHID	( config.Flag.HistogramId ),
	KEYBCL	( config.Flag.BCJetsLabeling ),

	IEVENT	( 0 ),

	histoManager(histoMng),
	histoRegistered( false )
{}

CJet::~CJet() {}

void CJet::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::CJet   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t... jets labeling ...\n");
	printf ("labeling on/off %s\n", KEYBCL ? "on" : "off");
	if (KEYBCL) {
		printf ("\tcjets ...\n");
		printf ("min c-quark p_T %lf\n", PTCMIN);
		printf ("max c-quark eta %lf\n", ETCMAX);
		printf ("max R_cj for c-jets %lf\n", RJC);
	}
	printf ("\n");
}

void CJet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	Int32_t idhist = 800 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "CJet: c-jets multiplicity"       , 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "CJet: c-quarks HP multiplicity"  , 10, 0.0, 10.0);
		histoManager 
			->registerHistogram(idhist+23, "CJet: delta r cjet-cquark HP"    , 50, 0.0,  1.0);
		histoManager
			->registerHistogram(idhist+24, "CJet: pTcjet/pTcquark HP"        , 50, 0.0,  2.0);
	}

	// do not use this algorithm
	if (!KEYBCL)
		return;

	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	
	// look for c-jets
	Int32_t NJETC = 0;
	for (int i=0; i<parts.size(); ++i) {
	  const Particle& part = parts[i];
	  
	  if (part.status == PS_CASCADE_QUARK && part.type == PT_CJET) {
	    // if there is a c-quark found before hadronization
	    // if there are still jets
	    if (orecord.Jets.empty()) continue; 
	    
	    Bool_t CJET_tagger = true;
	    
	    // and this c-quark is the last one in the FSR cascade
	    if (part.hasDaughter()) {
	      for (int j=part.daughters.first; j<=part.daughters.second; ++j) {
		if (parts[j].type == PT_CJET) 
		  CJET_tagger = false;
	      }
	    }
	    
	    if (!CJET_tagger) 
	      continue;
	    
	    if (part.pT() < PTCMIN) 
	      continue;
	    
	    if (abs(part.getEta()) > ETCMAX)
	      continue;
	    
	    // mark c-jet
	    Real64_t DR = 100.0;
	    Int32_t JETC = -1;
				
	    for (int j=0; j<orecord.Jets.size(); ++j) {
	      // if already marked as b-jet don't overwrite
	    
	      if (orecord.Jets[j].type != B_JET) {
		Real64_t DDR = sqrt(
				    pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
				    pow(part.getPhi() - orecord.Jets[j].phi_rec, 2)
				    );
		
		if (abs(part.getPhi() - orecord.Jets[j].phi_rec) > PI)
		  DDR = sqrt(
			     pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
			     pow(abs(part.getPhi() - orecord.Jets[j].phi_rec) - 2*PI, 2)
			     );
		
	      // select closest jet to c-jet tagger	      
		if (DDR < DR) { 
		  JETC = j;
		  DR = DDR;
		}
	      }
	    }
	    
	    // labell  c-jet
	    if (JETC >= 0 && DR <= RJC) {
	      orecord.Jets[JETC].type = C_JET;
	    }
	  }
	}

	//count c-jets
	for (int j=0; j<orecord.Jets.size(); ++j) {
          if( orecord.Jets[j].type == C_JET) NJETC++;
	}
	
	// store count in histogram
	histoManager
	  ->insert(idhist+11, NJETC, weight);
	
	// check partons
	Int32_t IQUAC = 0, ICJET = 0;
	for (int i=6; i<parts.size(); ++i) {
	  const Particle& part = parts[i];
	  
	  if (part.status == PS_HISTORY && part.type == PT_CJET) { 

	    if (abs(part.getEta()) < ETCMAX
		&& part.pT() > ETJET) {
	      
	      IQUAC++;
	      Real64_t DR = 18.0;
	      
	      for (int j=0; j<orecord.Jets.size(); ++j) {
		if (orecord.Jets[j].type == C_JET) {
		  
		  Real64_t DDR = sqrt(
				      pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
				      pow(part.getPhi() - orecord.Jets[j].phi_rec, 2)
				      );
		  
		  if (abs(part.getPhi() - orecord.Jets[j].phi_rec) > PI)
		    DDR = sqrt(
			       pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
			       pow(abs(part.getPhi() - orecord.Jets[j].phi_rec) - 2*PI, 2)
			       );
							
		  if (DDR < DR) { 
		    ICJET = j;
		    DR = DDR;
		  }
		}
	      }
	    }
	  }
	}

	// fill histogram
	histoManager
	  ->insert(idhist+21, IQUAC, weight );
	
	for (int i=0; i<orecord.Jets.size(); ++i) {
	  Real64_t PTREC = 0.0;
	  Real64_t DETRMIN = RCONE;
	  
	  if (orecord.Jets[i].type == C_JET) {
	    for (int j=0; j<parts.size(); ++j) {
	      const Particle& part = parts[j];
	      
	      if (isHardProcess(parts, j) )
		continue;
	      
	      if( part.type != PT_CJET) 
		continue;
	      
	      Real64_t DPHIA = abs(part.getPhi() - orecord.Jets[i].phi_rec);
	      if (DPHIA > PI) 
		DPHIA -= 2*PI;
	      
	      Real64_t DETR = sqrt(
				   pow(part.getEta() - orecord.Jets[i].eta_rec , 2) + 
				   pow(DPHIA, 2)
				   );
	      
	      if (DETR < DETRMIN) {
		PTREC = part.pT();
		DETRMIN = DETR;
	      }
	    }
	  }
	  
	  if (PTREC != 0) {
	    histoManager
	      ->insert(idhist + 23, DETRMIN, weight);
	    histoManager
	      ->insert(idhist + 24, orecord.Jets[i].pT / PTREC, weight);
	  }
	}
}

void CJet::printResults() const {
	if (KEYBCL) {
		printf ("***********************************\n");
		printf ("*                                 *\n");
		printf ("*     ***********************     *\n");
		printf ("*     ***   Output from   ***     *\n");
		printf ("*     ***  analyse::CJet  ***     *\n");
		printf ("*     ***********************     *\n");
		printf ("*                                 *\n");
		printf ("***********************************\n");
	
		printf (" Analysed records: %d\n", IEVENT);
	}
}
