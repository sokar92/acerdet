#include "Test_Histograms.h"
#include <cstdio>
#include <iostream>
using namespace AcerDet::test;

using namespace AcerDet::analyse;

#include "../core/Smearing.h"
#include "../core/Functions.h"
using namespace AcerDet::core;

Test_Histograms::Test_Histograms( const Configuration& config, IHistogramManager* histoMng ) :
	
	KEYHID	( config.Flag.HistogramId ),
	TEST	( config.Flag.Test ),

	IEVENT	( 0 ),
	
	histoManager(histoMng),
	histoRegistered( false )

{}

Test_Histograms::~Test_Histograms() {}

void Test_Histograms::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::Test   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("Test mode: %s\n", TEST ? "on" : "off");
	printf ("\n");
}

void Test_Histograms::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	if ( TEST ){
		Int32_t muon_idhist = 2200 + KEYHID;
		Int32_t ele_idhist = 2300 + KEYHID;
		Int32_t pho_idhist = 2400 + KEYHID;
		Int32_t jet_idhist = 2500 + KEYHID;
		Int32_t mis_idhist = 2600 + KEYHID;
		Int32_t bjet_idhist = 2700 + KEYHID;
		Int32_t cjet_idhist = 2800 + KEYHID;
		if (!histoRegistered) {
			histoRegistered = true;
				
			histoManager
				->registerHistogram(muon_idhist+53, "Muon: Higher Z boson mass ISOLATED. H->2Z->4mu", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(muon_idhist+54, "Muon: Lower Z boson mass ISOLATED. H->2Z->4mu", 110, 0.0, 110.0);	
			histoManager
				->registerHistogram(muon_idhist+56, "Muon: Higgs mass ISOLATED. H->2Z->4mu", 100, 100, 140);
			histoManager
				->registerHistogram(muon_idhist+63, "Muon: Higher Z boson mass HARD. H->2Z->4mu", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(muon_idhist+64, "Muon: Lower Z boson mass HARD. H->2Z->4mu", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(muon_idhist+66, "Muon: Higgs mass HARD. H->2Z->4mu", 100, 100.0, 140.0);

		
			histoManager
				->registerHistogram(ele_idhist+41, "Electron: Z boson mass HARD-ISOLATED. Z->2e", 100, -50.0, 50.0);
			histoManager
				->registerHistogram(ele_idhist+51, "Electron: Z boson mass ISOLATED. Z->2e", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(ele_idhist+53, "Electron: Higher Z boson mass ISOLATED. H->2Z->4e", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(ele_idhist+54, "Electron: Lower Z boson mass ISOLATED. H->2Z->4e", 110, 0.0, 110.0);	
			histoManager
				->registerHistogram(ele_idhist+56, "Electron: Higgs mass ISOLATED. H->2Z->4e", 100, 100, 140);
			histoManager
				->registerHistogram(ele_idhist+61, "Electron: Z boson mass HARD. Z->2e", 110, 0.0, 110.0);	
			histoManager
				->registerHistogram(ele_idhist+63, "Electron: Higher Z boson mass HARD. H->2Z->4e", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(ele_idhist+64, "Electron: Lower Z boson mass HARD. H->2Z->4e", 110, 0.0, 110.0);
			histoManager
				->registerHistogram(ele_idhist+66, "Electron: Higgs mass HARD. H->2Z->4e", 100, 100.0, 140.0);
			
			
			histoManager
				->registerHistogram(pho_idhist+51, "Photon: Higgs mass ISOLATED", 100, 100.0, 140.0);
			histoManager
				->registerHistogram(pho_idhist+61, "Photon: Higgs mass HP", 100, 100.0, 140.0);
			histoManager
				->registerHistogram(pho_idhist+71, "Photon: photon pt ISOLATED", 100, 0.0, 100.0);
			histoManager
				->registerHistogram(pho_idhist+81, "Photon: photon pt HP", 100, 0.0, 100.0);
			
			histoManager
				->registerHistogram(jet_idhist+51, "UJet: Higgs mass from ujets", 100, 80.0,  140.0);
			histoManager
				->registerHistogram(jet_idhist+61, "UJet: Higgs mass from uquarks HP", 100, 80.0, 140.0);
			
			histoManager
				->registerHistogram(mis_idhist+2, "Mis: Px_miss - Px_nu", 200, -100, 100);
			histoManager
				->registerHistogram(mis_idhist+7, "Mis: Py_miss - Py_nu", 200, -100, 100);
			histoManager
				->registerHistogram(mis_idhist+12, "Mis: Px_calo - Px_nu", 200, -100, 100);
			histoManager
				->registerHistogram(mis_idhist+17, "Mis: Py_calo - Py_nu", 200, -100, 100);
			
			histoManager
				->registerHistogram(bjet_idhist+51, "BJet: Higgs mass from bjets", 100, 80.0,  140.0);
			histoManager
				->registerHistogram(bjet_idhist+61, "BJet: Higgs mass from bquarks HP", 100, 80.0, 140.0);
			
			histoManager
				->registerHistogram(cjet_idhist+51, "CJet: Higgs mass from cjets", 100, 80.0,  140.0);
			histoManager
				->registerHistogram(cjet_idhist+61, "CJet: Higgs mass from cquarks HP", 100, 80.0, 140.0);
				
		
		}

		// create container for given event
		Test test;
	

		// new event to compute
		IEVENT++;

		// reference to particles container
		const vector<Particle>& parts = irecord.particles();
		
		// filling test container
		for (int j=0; j<parts.size(); ++j) {
		
			// pick only particles from hard process
			if ( isHardProcess(parts, j) 
						|| parts[j].status == PS_BEAM ||  parts[j].status == PS_HISTORY ) {
				test.PutHARD( parts[j] );
		  
			}
		}
		
		for (int j=0; j<orecord.Photons.size(); ++j) {
			test.ISOphoton.push_back( orecord.Photons[j] );
		}
		
		for (int j=0; j<orecord.Muons.size(); ++j) {
			test.ISOmuon.push_back( orecord.Muons[j] );
		}	
	
		for (int j=0; j<orecord.Electrons.size(); ++j) {
			test.ISOelectron.push_back( orecord.Electrons[j] );
		}
	
		for (int j=0; j<orecord.Jets.size(); ++j) {
			test.jet.push_back( orecord.Jets[j] );
		}
	
	
		Real64_t pxmiss = - orecord.Miss.PXREC;
		Real64_t pymiss = - orecord.Miss.PYREC;
		Real64_t pxnue = orecord.Miss.PXNUE;
		Real64_t pynue = orecord.Miss.PYNUE;
		Real64_t pxcalo = orecord.Miss.PXXCALO;
		Real64_t pycalo = orecord.Miss.PYYCALO;

		
		test.Dividejets();

		// routines for filling specific types of histograms
		test.Fill4part( test.HARDmuon, muon_idhist, weight, histoManager );
		test.Fill4object( test.ISOmuon, muon_idhist, weight, histoManager );
		test.Fill4part( test.HARDelectron, ele_idhist, weight, histoManager );
		test.Fill4object( test.ISOelectron, ele_idhist, weight, histoManager );
		test.Fill2e( ele_idhist, weight, histoManager );
		test.Fill2gamma( pho_idhist, weight, histoManager );

		test.Filljets( test.bquark, test.bjet, bjet_idhist, weight, histoManager );
		test.Filljets( test.cquark, test.cjet, cjet_idhist, weight, histoManager );
		test.Filljets( test.uquark, test.ujet, mis_idhist, weight, histoManager );
	
	
		histoManager
			->insert(mis_idhist+2, pxmiss - pxnue, weight);
		histoManager
			->insert(mis_idhist+7, pymiss - pynue, weight);
		
		histoManager
			->insert(mis_idhist+12, pxcalo - pxnue, weight);
		histoManager
			->insert(mis_idhist+17, pycalo - pynue, weight);
	}
}

void Test_Histograms::printResults() const {
	printf ("***********************************\n");
	printf ("*                                 *\n");
	printf ("*     ***********************     *\n");
	printf ("*     ***   Output from   ***     *\n");
	printf ("*     ***  analyse::Test  ***     *\n");
	printf ("*     ***********************     *\n");
	printf ("*                                 *\n");
	printf ("***********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
}
