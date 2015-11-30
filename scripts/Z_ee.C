#define Z_ee_cxx
#include "Z_ee.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Z_ee::Loop()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	TH1F *gen_rec_mass = new TH1F("shist12341" , "Electron: Z boson mass HARD-ISOLATED", 100, -50, 50 );
	TH1F *rec_mass = new TH1F("shist12351" , "Electron: Z boson mass ISOLATED", 100, 80, 140 ); 
	TH1F *gen_mass = new TH1F("shist12361" , "Electron: Z boson mass HARD", 100, 80, 128 );
   

	float  m_gen = 0, m_rec = 0, E_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0;	
	Long64_t gen_ele_event = 0;
	for (Long64_t jentry = 0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	GetEntry(jentry);
	m_rec = 0;
	m_gen = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
        pz_sum = 0;
	gen_ele_event = 0;


	for (Long64_t kentry = 0; kentry < part_pdgId->size(); kentry++){

	  if(( part_pdgId->at(kentry) == 11 || part_pdgId->at(kentry) == -11 ) && part_mother_pdgId->at(kentry) == 23 ){	
	  
	    E_sum += part_E->at(kentry);
	    px_sum += part_px->at(kentry);
	    py_sum += part_py->at(kentry);
	    pz_sum += part_pz->at(kentry);    
	    
	    gen_ele_event++;
	    
	  }
	}
	if (gen_ele_event == 2) {
	  m_gen = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	  gen_mass->Fill(m_gen);
	}

	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;

	if (ele_n == 2){

      	  E_sum = ele_E->at(0) + ele_E->at(1);
	  px_sum = ele_px->at(0) + ele_px->at(1);
	  py_sum = ele_py->at(0) + ele_py->at(1);
	  pz_sum = ele_pz->at(0) + ele_pz->at(1);	
	  m_rec = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	 
	  rec_mass->Fill(m_rec);
	  gen_rec_mass->Fill(m_rec-m_gen);

	}
   }

}
