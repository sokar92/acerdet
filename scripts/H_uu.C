#define H_uu_cxx
#include "H_uu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void H_uu::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F *pt_pt = new TH1F("shist12541" , "UJet: pTujet/pTuquark HP", 100, 0, 2 );
   TH1F *rec_mass = new TH1F("shist12551" , "UJet: Higgs mass from ujets", 100, 0, 200 );
   TH1F *gen_mass = new TH1F("shist12561" , "UJet: Higgs mass from uquarks HP
", 100, 122, 128 );
   

   vector<float> E_2;
   vector<float> px_2;
   vector<float> py_2;
   vector<float> pz_2;
   vector<float> pt_1;
   vector<float> pt_2;
   vector<float> eta_1;
   vector<float> eta_2;
   vector<float> phi_1;
   vector<float> phi_2;
   vector<float> r;
   float  m = 0, E_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0, pt_0 = 0, pt_a = 0, pt_b = 0, eta_a = 0, eta_b = 0, phi_a = 0, r_a = 0, factor = 0, qualified = 0, x1 = 0, x2 = 0;	
   Long64_t gen_pho_event = 0, gen_pho = 0, rec_pho = 0, part_numb = 0, selected = 0, winner1 = 0, winner2 = 0, run_up1 = 0, run_up2 = 0;
  
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
     part_numb = part_pdgId->size();
     E_2.clear();
     px_2.clear();
     py_2.clear();
     pz_2.clear();
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

	 
       if( (part_pdgId->at(kentry) == 1 || part_pdgId->at(kentry) == -1) && part_mother_pdgId->at(kentry) == 25 ){	
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
     
     if ( jet_n != 0 ){
       for( Long64_t zentry = 0; zentry < jet_n; zentry++ ){
	  
	 if (jet_pdgId->at(zentry) == 0 ){
	   factor = 1;
	   E_sum += jet_E->at(zentry);
	   px_sum += jet_px->at(zentry);
	   py_sum += jet_py->at(zentry);
	   pz_sum += jet_pz->at(zentry);	
	   pt_a = sqrt(pow(jet_px->at(zentry),2)+pow(jet_py->at(zentry),2));
	   if ( jet_pz->at(zentry) < 0 ) factor = -1.0;
	   eta_a =  log( (sqrt(pow(pt_a,2)+pow(jet_pz->at(zentry),2)) + (factor * jet_pz->at(zentry))) / pt_a) * factor;
	   phi_a = TMath::ASinH( (jet_py->at(zentry)) / pt_a );
	    
	  
	   E_2.push_back(jet_E->at(zentry));
	   px_2.push_back(jet_px->at(zentry));
	   py_2.push_back(jet_py->at(zentry));
	   pz_2.push_back(jet_pz->at(zentry));
	   pt_2.push_back(pt_a);
	   eta_2.push_back(eta_a);
	   phi_2.push_back(phi_a);
	 }
       }	  


       if( eta_1.size() == 2 ){

	 for ( Long64_t count1 = 0; count1<2; count1++ ){
	   for ( Long64_t count2 = 0; count2<eta_2.size(); count2++ ){
	     r_a = sqrt(pow(phi_1.at(count1)-phi_2.at(count2),2)+pow(eta_1.at(count1)-eta_2.at(count2),2));
	     r.push_back(r_a);
	  
	   }
	 }

	 winner1 = 0;
	 winner2 = 0;
	 run_up1 = 0;
	 run_up2 = 0;
	 qualified = 0;
	 E_sum = 0;
	 px_sum = 0;
	 py_sum = 0;
	 pz_sum = 0;
	
	 for ( Long64_t count = 1; count<eta_2.size(); count++ ){
	   if ( r.at(winner1) > r.at(count) ) run_up1 = winner1, winner1 = count;
	   if ( r.at(winner2 + eta_2.size()) > r.at(count + eta_2.size()) ) run_up2 = winner2, winner2 = count;
	 }
	
	 if ( winner1 == winner2 ) {
	   if ( r.at(winner1) <= r.at(eta_2.size() + winner2) ) winner2 = run_up2;
	   if ( r.at(winner1) > r.at(eta_2.size() + winner2) ) winner1 = run_up1;
	 }
        
	 if ( r.at( winner1 ) < 0.4 ){
	   qualified++;
	   pt_pt->Fill ( pt_2.at(winner1)/pt_1.at(0) );
	   E_sum += E_2.at(winner1);
	   px_sum += px_2.at(winner1);
	   py_sum += py_2.at(winner1);
	   pz_sum += pz_2.at(winner1);
	 }
	 if ( r.at( winner2 + eta_2.size() ) < 0.4 ){
	   qualified++;
	   pt_pt->Fill ( pt_2.at(winner2)/pt_1.at(1) );
	   E_sum += E_2.at(winner2);
	   px_sum += px_2.at(winner2);
	   py_sum += py_2.at(winner2);
	   pz_sum += pz_2.at(winner2);
	 }
	 if ( qualified == 2 ){
	 m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	 rec_mass->Fill(m);
	 }


       }

     }
   }

}
