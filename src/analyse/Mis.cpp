#include "Mis.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Smearing.h"
using namespace AcerDet::core;

Mis::Mis( const Configuration& config, IHistogramManager *histoMng ) :
	PTMUMIN	( config.Muon.MinMomenta ),
	ETAMAX	( config.Muon.MaxEta ),
	ETCELL	( config.Tau.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histoManager(histoMng),
	histoRegistered( false )
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

void Mis::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	Int32_t idhist = 600 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "Mis: reconstructed p_T ", 50, 0.0, 200.0);
		histoManager
			->registerHistogram(idhist+12, "Mis: reconstructed p_T + cells", 50, 0.0, 200.0);
		histoManager
			->registerHistogram(idhist+13, "Mis: reconstructed pTmiss", 50, 0.0, 200.0);
		histoManager
			->registerHistogram(idhist+21, "Mis: true p_T invisible", 50, 0.0, 200.0);
	}

    // new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	//sum up reconstructed momenta - firstly clear
	orecord.Miss.clear();

	// add jets
	for (int i=0; i<orecord.Jets.size(); ++i) {
		const JetData jet = orecord.Jets[i];
		orecord.Miss.PXREC += jet.pT * cos( jet.phi_rec );
		orecord.Miss.PYREC += jet.pT * sin( jet.phi_rec );
		orecord.Miss.PXXCALO += jet.pT * cos( jet.phi_rec );
		orecord.Miss.PYYCALO += jet.pT * sin( jet.phi_rec );
		orecord.Miss.SUMET += jet.pT;
	}
	
	// add non-used clusters
	for (int i=0; i<orecord.Clusters.size(); ++i) {
		const ClusterData& cluster = orecord.Clusters[i];
		if (!cluster.alreadyUsed) {
			orecord.Miss.PXREC += cluster.pT * cos( cluster.phi_rec );
			orecord.Miss.PYREC += cluster.pT * sin( cluster.phi_rec );
			orecord.Miss.PXXCALO += cluster.pT * cos( cluster.phi_rec );
			orecord.Miss.PYYCALO += cluster.pT * sin( cluster.phi_rec );
			orecord.Miss.SUMET += cluster.pT;
		}
	}
	
	// add isolated muons
	for (int i=0; i<orecord.Muons.size(); ++i) {
		const ObjectData& muon = orecord.Muons[i];
		orecord.Miss.PXREC += muon.pT * cos( muon.phi );
		orecord.Miss.PYREC += muon.pT * sin( muon.phi );
	}
	
	// add non-isolated muons not added to clusters
	for (int i=0; i<orecord.NonisolatedMuons.size(); ++i) {
		const ObjectData& muon = orecord.NonisolatedMuons[i];
		if (!muon.alreadyUsed) {
			orecord.Miss.PXREC += muon.pT * cos( muon.phi );
			orecord.Miss.PYREC += muon.pT * sin( muon.phi );
		} else {
			orecord.Miss.SUMET -= muon.pT;
		}
	}
	
	// add isolated electrons
	for (int i=0; i<orecord.Electrons.size(); ++i) {
		const ObjectData& ele = orecord.Electrons[i];
		orecord.Miss.PXREC += ele.pT * cos( ele.phi );
		orecord.Miss.PYREC += ele.pT * sin( ele.phi );
		orecord.Miss.PXXCALO += ele.pT * cos( ele.phi );
		orecord.Miss.PYYCALO += ele.pT * sin( ele.phi );
		orecord.Miss.SUMET += ele.pT;
	}
	
	// add isolated photons
	for (int i=0; i<orecord.Photons.size(); ++i) {
		const ObjectData& pho = orecord.Photons[i];
		orecord.Miss.PXREC += pho.pT * cos( pho.phi );
		orecord.Miss.PYREC += pho.pT * sin( pho.phi );
		orecord.Miss.PXXCALO += pho.pT * cos( pho.phi );
		orecord.Miss.PYYCALO += pho.pT * sin( pho.phi );
		orecord.Miss.SUMET += pho.pT;
    }
    
    // store pT in histo
    Real64_t ETREC = sqrt(
		pow(orecord.Miss.PXREC, 2) +
		pow(orecord.Miss.PYREC, 2)
	);
    
    histoManager
      ->insert(idhist+11, ETREC, weight );
    
    // smear cells energy not used for reconstruction
    // remove cells below threshold
    // add momenta in cells not used for reconstruction
    for (int i=0; i<orecord.Cells.size(); ++i) {
		const CellData& cell = orecord.Cells[i];
		if (cell.status != 0) {
			if (KEYSME) {
				Real64_t PEI = cell.pT * cosh( cell.eta );
				Real64_t SIGSME = Smearing::forHadron(PEI, cell.eta, CALOTH);
				
				Real64_t PSME = PEI * (1.0 + SIGSME) / cosh( cell.eta );
				if (PSME < ETCELL) 
					PSME = 0.0;
					
				orecord.Miss.PXSUM += PSME * cos( cell.phi );
				orecord.Miss.PYSUM += PSME * sin( cell.phi );
				orecord.Miss.PXXCALO += PSME * cos( cell.phi );
				orecord.Miss.PYYCALO += PSME * sin( cell.phi );
				orecord.Miss.SUMET += PSME;
			}
		}
	}
	
	orecord.Miss.PXSUM += orecord.Miss.PXREC;
	orecord.Miss.PYSUM += orecord.Miss.PYREC;

	Real64_t ETSUM = sqrt(
		pow(orecord.Miss.PXSUM, 2) +
		pow(orecord.Miss.PYSUM, 2)
	);
	
	histoManager
	  ->insert(idhist+12, ETSUM, weight );
	
	//PXXMISS, PYYMISS
	Real64_t PXXMISS = -orecord.Miss.PXSUM;
	Real64_t PYYMISS = -orecord.Miss.PYSUM;
	Real64_t PTMISS = sqrt(
		pow(PXXMISS, 2) +
		pow(PYYMISS, 2)
	);
	
	histoManager
	  ->insert(idhist+13, PTMISS, weight );
	
	// sum up momenta  of neutrinos 
	Real64_t PXXNUES = 0.0;
	Real64_t PYYNUES = 0.0;
	for (int i=0; i<parts.size(); ++i) {
		
		if (parts[i].status != PS_FINAL)
			continue;
		
		if (parts[i].isNeutrino()
		|| abs(parts[i].pdg_id) == KFINVS) {
			PXXNUES += parts[i].pX();
			PYYNUES += parts[i].pY();
		}
		
	}
	
	Real64_t PTNUES = sqrt(
		pow(PXXNUES, 2) +
		pow(PYYNUES, 2)
	);
	
	orecord.Miss.PXNUE = PXXNUES;
	orecord.Miss.PYNUE = PYYNUES;
	
	histoManager
	  ->insert(idhist+21, PTNUES, weight );
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
}
