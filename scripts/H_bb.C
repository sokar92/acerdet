#define H_bb_cxx
#include "H_bb.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void H_bb::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   TH1F *pt_pt = new TH1F("shist10724" , "BJet: pTbjet/pTbquark HP", 100, 0, 2 );
   TH1F *rec_mass = new TH1F("shist12751" , "BJet: Higgs mass from bjets", 100, 0, 200 );
   TH1F *gen_mass = new TH1F("shist12761" , "BJet: Higgs mass from bquarks HP", 100, 122, 128 );


   vector<float> pt_1;
   vector<float> pt_2;
   vector<float> eta_1;
   vector<float> eta_2;
   vector<float> phi_1;
   vector<float> phi_2;
   vector<float> r;
   float  m = 0, E_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0, pt_0 = 0, pt_a = 0, pt_b = 0, eta_a = 0, eta_b = 0, phi_a = 0, r_a = 0, factor = 0, qualified = 0, b1 = 0, b2 = 0, b3 = 0, x1 = 0, x2 = 0;	
   Long64_t gen_pho_event = 0, gen_pho = 0, rec_pho = 0, selected = 0;

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
	pt_1.clear();
	pt_2.clear();
	eta_1.clear();
	eta_2.clear();
	phi_1.clear();
	phi_2.clear();
	r.clear();
	phi_a = 0;
	eta_1 = 0;
	pt_a = 0;

	for (Long64_t kentry = 0; kentry < part_pdgId->size(); kentry++){

	  if( (part_pdgId->at(kentry) == 5 || part_pdgId->at(kentry) == -5) && part_mother_pdgId->at(kentry) == 25 ){	
		  factor = 1.0;
		  E_sum += part_E->at(kentry);
		  px_sum += part_px->at(kentry);
		  py_sum += part_py->at(kentry);
		  pz_sum += part_pz->at(kentry);
		  pt_a = sqrt(pow(part_px->at(kentry),2)+pow(part_py->at(kentry),2));
		  if ( part_pz->at(kentry) < 0 ) factor = -1.0;
		  eta_a =  log( (sqrt(pow(pt_a,2)+pow(part_pz->at(kentry),2)) + (factor * part_pz->at(kentry))) / pt_a) * factor;
		  phi_a = TMath::ASinH( (part_py->at(kentry)) / pt_a );
	        
	  
		  pt_1.push_back(pt_a);
		  eta_1.push_back(eta_a);
		  phi_1.push_back(phi_a);
		 
		  gen_pho_event++;
		
		}
	}
	if (gen_pho_event == 2) {
	        m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
		pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
		gen_mass->Fill(m);
	        
	}
	
	
        
      gen_pho += gen_pho_event;
      rec_pho += pho_n;
      E_sum = 0;
      px_sum = 0;
      py_sum = 0;
      pz_sum = 0;
      qualified = 0;
      
      for( Long64_t zentry = 0; zentry < jet_n; zentry++ ){

	if (jet_pdgId->at(zentry) == 5 || jet_pdgId->at(zentry)== -5){
	  qualified++;
	  factor = 1;
      	  E_sum += jet_E->at(zentry);
	  px_sum += jet_px->at(zentry);
	  py_sum += jet_py->at(zentry);
	  pz_sum += jet_pz->at(zentry);	
	  pt_a = sqrt(pow(jet_px->at(zentry),2)+pow(jet_py->at(zentry),2));
	  if ( jet_pz->at(zentry) < 0 ) factor = -1.0;
	  eta_a =  log( (sqrt(pow(pt_a,2)+pow(jet_pz->at(zentry),2)) + (factor * jet_pz->at(zentry))) / pt_a) * factor;
	  phi_a = TMath::ASinH( (jet_py->at(zentry)) / pt_a );

	  
	  pt_2.push_back(pt_a);
	  eta_2.push_back(eta_a);
	  phi_2.push_back(phi_a);
	}
      }	  
      if ( qualified == 2 ){
	m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	rec_mass->Fill(m);
      }

      b1 = 0;
      b2 = 0;
      b3 = 0;

      if( eta_1.size() == 2 && eta_2.size() == 2 ){

	for ( Long64_t count1 = 0; count1<2; count1++ ){
	  for ( Long64_t count2 = 0; count2<2; count2++ ){
	    r_a = sqrt(pow(phi_1.at(count1)-phi_2.at(count2),2)+pow(eta_1.at(count1)-eta_2.at(count2),2));
	    r.push_back(r_a);
        
	  }
	}
	if ( r.at(0) > r.at(1) ) b1 = 1;
	if ( r.at(2) > r.at(3) ) b2 = 1;
	if ( r.at(b1) > r.at(2+b2) ) b3 = 1;
	x2 = b3;
	x1 = (b3 * b2) + ((1-b3)*b1);
	if ( r.at( 2*x2 + x1 ) < 0.4 ){
	  pt_pt->Fill ( pt_2.at(x1)/pt_1.at(x2) );
	}
	if ( r.at( 2*(1-x2) + (1-x1) ) < 0.4 ){
	  pt_pt->Fill ( pt_2.at(1-x1)/pt_1.at(1-x2) );
	}
      }

 
   }

}
