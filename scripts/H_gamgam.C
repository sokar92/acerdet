#define H_gamgam_cxx
#include "H_gamgam.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>

void H_gamgam::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F *rec_mass = new TH1F("shist12451" , "Photon: higgs mass ISOLATED", 100, 110, 140 );
   TH1F *rec_mass_sel = new TH1F("shist12456" , "Photon: higgs mass ISOLATED SEL", 100, 110, 140 );
   TH1F *gen_mass = new TH1F("shist12461" , "Photon: higgs mass HP", 100, 122, 128 );
   TH1F *gen_mass_sel = new TH1F("shist12466" , "Photon: higgs mass HP SEL", 100, 122, 128 );
   TH1F *pho_pt_rec = new TH1F("shist12471" , "Photon: photon pt ISOLATED", 100, 0, 140 );
   TH1F *pho_pt_rec_sel = new TH1F("shist12476" , "Photon: photon pt ISOLATED SEL", 100, 0, 140 );
   TH1F *pho_pt_gen = new TH1F("shist12481" , "Photon: photon pt HP", 100, 0, 140 );
   TH1F *pho_pt_gen_sel = new TH1F("shist12486" , "Photon: photon pt HP SEL", 100, 0, 140 );
   vector<float> pt;
   vector<float> eta;
   float  m = 0, E_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0, pt_0 = 0, pt_a = 0, pt_b = 0, eta_a = 0, eta_b = 0, factor = 0, qualified = 0;	
   Long64_t nbytes = 0, gen_pho_event = 0, gen_pho = 0, rec_pho = 0, selected = 0;

   for (Long64_t jentry = 0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	GetEntry(jentry);
	m = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
        pz_sum = 0;
	gen_pho_event = 0;
	pt.clear();
	eta.clear();
	for (Long64_t kentry = 0; kentry < part_pdgId->size(); kentry++){
	  if( part_pdgId->at(kentry) == 22 && part_mother_pdgId->at(kentry) == 25 ){	
	    factor = 1.0;
	    E_sum += part_E->at(kentry);
	    px_sum += part_px->at(kentry);
	    py_sum += part_py->at(kentry);
	    pz_sum += part_pz->at(kentry);
	    pt_a = sqrt(pow(part_px->at(kentry),2)+pow(part_py->at(kentry),2));
	    if ( part_pz->at(kentry) < 0 ) factor = -1.0;
	    eta_a =  log( (sqrt(pow(pt_a,2)+pow(part_pz->at(kentry),2)) + (factor * part_pz->at(kentry))) / pt_a) * factor;
	    
	    pt.push_back(pt_a);
	    eta.push_back(eta_a);

	    gen_pho_event++;
	  
	  }
	}
	if (gen_pho_event == 2) {
	  m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	  pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
	  gen_mass->Fill(m);
	 
	  pho_pt_gen->Fill(pt.at(0));
	  pho_pt_gen->Fill(pt.at(1));
	       
	  if( ( (eta.at(0) < 2.37 && eta.at(0) > -2.37)  && (eta.at(1) < 2.37 && eta.at(1) > -2.37) ) && ( ( pt.at(0) > 0.35/m && pt.at(1) > 0.25/m ) || ( pt.at(0) > 0.25/m && pt.at(1) > 0.35/m ) ) ){
	    gen_mass_sel->Fill(m);
	   
	    pho_pt_gen_sel->Fill(pt.at(0));
	    pho_pt_gen_sel->Fill(pt.at(1));
	  }
	}
	
	
        
      gen_pho += gen_pho_event;
      rec_pho += pho_n;
      E_sum = 0;
      px_sum = 0;
      py_sum = 0;
      pz_sum = 0;

      if (pho_n == 2){
       	qualified++;

	  factor = 1;
      	  E_sum = pho_E->at(0) + pho_E->at(1);
	  px_sum = pho_px->at(0) + pho_px->at(1);
	  py_sum = pho_py->at(0) + pho_py->at(1);
	  pz_sum = pho_pz->at(0) + pho_pz->at(1);	
	  m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	  pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
	  rec_mass->Fill(m);


	  pt_a = sqrt(pow(pho_px->at(0),2)+pow(pho_py->at(0),2));
	  if ( pho_pz->at(0) < 0 ) factor = -1.0;
	  eta_a = log( (sqrt(pow(pt_a,2)+pow(pho_pz->at(0),2)) + (factor * pho_pz->at(0))) / pt_a) * factor;
	  factor = 1;
	  pt_b = sqrt(pow(pho_px->at(1),2)+pow(pho_py->at(1),2));
	  if ( pho_pz->at(1) < 0 ) factor = -1.0;
	  eta_b = log( (sqrt(pow(pt_b,2)+pow(pho_pz->at(1),2)) + (factor * pho_pz->at(1))) / pt_b) * factor;
      
	  pho_pt_rec->Fill(pt_a);
	  pho_pt_rec->Fill(pt_b);
	 
	  if( ( (eta_a < 2.37 && eta_a > -2.37)  && (eta_b < 2.37 && eta_b > -2.37) ) && ( ( pt_a > 0.35/m && pt_b > 0.25/m ) || ( pt_a > 0.25/m && pt_b > 0.35/m ) ) ){
	    rec_mass_sel->Fill(m);
	   
	    pho_pt_rec_sel->Fill(pt_a);
	    pho_pt_rec_sel->Fill(pt_b);
	    selected++;
	  }

      }
   }
   cout << " selection rate: " << (selected/qualified)*100.0 <<  "%" << endl;
}
