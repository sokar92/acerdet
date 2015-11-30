#define H_4mu_cxx
#include "H_4mu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void H_4mu::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F *gen_mass = new TH1F("shist12266" , "Muon: Higgs mass HP", 100, 80, 128 );
   TH1F *gen_mass_sel = new TH1F("shist12267" , "Muon: Higgs mass HP SEL", 100, 80, 128 );
   TH1F *rec_mass = new TH1F("shist12256" , "Muon: Higgs mass ISO", 100, 80, 140 );
   TH1F *rec_mass_sel = new TH1F("shist12257" , "Muon: Higgs mass ISO SEL", 100, 80, 140 );
 

   TH1F *z_gen_mass_i_a = new TH1F("shist122631" , "Muon: Higher Z boson mass HP I", 100, 0, 140 );
   TH1F *z_gen_mass_i_b = new TH1F("shist122641" , "Muon: Lower Z boson mass HP I", 100, 0, 140 );
   TH2F *z_gen_mass_i_ab = new TH2F("shist122131" , "Muon: Z bosons masses HP I", 100, 0, 140, 100, 0, 140 );
   
   TH1F *z_gen_mass_ii_a = new TH1F("shist122632" , "Muon: Higher Z boson mass HP II", 100, 0, 140 );
   TH1F *z_gen_mass_ii_b = new TH1F("shist122642" , "Muon: Lower Z boson mass HP II", 100, 0, 140 );
   TH2F *z_gen_mass_ii_ab = new TH2F("shist122132" , "Muon: Z bosons masses HP II", 100, 0, 140, 100, 0, 140 );
   
   TH1F *z_gen_mass_iii_a = new TH1F("shist122633" , "Muon: Higher Z boson mass HP III", 100, 0, 140 );
   TH1F *z_gen_mass_iii_b = new TH1F("shist122643" , "Muon: Lower Z boson mass HP III", 100, 0, 140 );
   TH2F *z_gen_mass_iii_ab = new TH2F("shist122133" , "Muon: Z bosons masses HP III", 100, 0, 140, 100, 0, 140 );
 
   TH1F *z_rec_mass_i_a = new TH1F("shist122531" , "Muon: Higher Z boson mass ISO I", 100, 0, 140 );
   TH1F *z_rec_mass_i_b = new TH1F("shist122541" , "Muon: Lower Z boson mass ISO I", 100, 0, 140 );
   TH2F *z_rec_mass_i_ab = new TH2F("shist122141" , "Muon: Z bosons masses ISO I", 100, 0, 140, 100, 0, 140 );
  
   TH1F *z_rec_mass_ii_a = new TH1F("shist122532" , "Muon: Higher Z boson mass ISO II", 100, 0, 140 );
   TH1F *z_rec_mass_ii_b = new TH1F("shist122542" , "Muon: Lower Z boson mass ISO II", 100, 0, 140 );
   TH2F *z_rec_mass_ii_ab = new TH2F("shist122142" , "Muon: Z bosons masses ISO II", 100, 0, 140, 100, 0, 140 );

   TH1F *z_rec_mass_iii_a = new TH1F("shist122533" , "Muon: Higher Z boson mass ISO III", 100, 0, 140 );
   TH1F *z_rec_mass_iii_b = new TH1F("shist122543" , "Muon: Lower Z boson mass ISO III", 100, 0, 140 );
   TH2F *z_rec_mass_iii_ab = new TH2F("shist122143" , "Muon: Z bosons masses ISO III", 100, 0, 140, 100, 0, 140 );

   vector<float> E;
   vector<float> px;
   vector<float> py;
   vector<float> pz;
   vector<float> pt;
   vector<float> muo_pt;
   vector<float> eta;
   vector<float> id_plus;
   vector<float> id_minus;
   vector<float> muo_id_plus;
   vector<float> muo_id_minus;
	vector<float> mass;
   float  m = 0, m_1 = 0, m_2 = 0 , E_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0, pt_0 = 0, pt_a = 0, pt_b = 0, eta_a = 0, eta_b = 0, factor = 0, qualified = 0, x1 = 0, x2 = 0;
	
   float  pt_min = 6, eta_max = 2.7;

   Long64_t nbytes = 0, nb = 0, gen_muo_event = 0, gen_muo = 0, rec_muo = 0, part_numb = 0, selected = 0, good_1 = 0, b1 = 0, b2 = 0, b3 = 0,
good_2 = 0;

   for (Long64_t jentry = 0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	GetEntry(jentry);
	m = 0;
	m_1 = 0;
	m_2 = 0;
	E_sum = 0;
	good_1 = 0;
	good_2 = 0;
	px_sum = 0;
	py_sum = 0;
        pz_sum = 0;
	gen_muo_event = 0;
	part_numb  = part_pdgId->size();
	pt.clear();
	eta.clear();
	E.clear();
	px.clear();
	py.clear();
	pz.clear();
	id_plus.clear();
	id_minus.clear();
	muo_id_plus.clear();
	muo_id_minus.clear();
	muo_pt.clear();

	for (Long64_t kentry = 0; kentry < part_pdgId->size(); kentry++){
	 
	  if(( part_pdgId->at(kentry) == 13 || part_pdgId->at(kentry) == -13 ) && part_mother_pdgId->at(kentry) == 23 ){	
		  factor = 1.0;
		  E_sum += part_E->at(kentry);
		  px_sum += part_px->at(kentry);
		  py_sum += part_py->at(kentry);
		  pz_sum += part_pz->at(kentry);
		  pt_a = sqrt(pow(part_px->at(kentry),2)+pow(part_py->at(kentry),2));
		  if ( part_pz->at(kentry) < 0 ) factor = -1.0;
		  eta_a =  log( (sqrt(pow(pt_a,2)+pow(part_pz->at(kentry),2)) + (factor * part_pz->at(kentry))) / pt_a) * factor;
		  
		  E.push_back(part_E->at(kentry));
		  px.push_back(part_px->at(kentry));
		  py.push_back(part_py->at(kentry));
		  pz.push_back(part_pz->at(kentry));
		  pt.push_back(pt_a);
		  eta.push_back(eta_a);
		  if (part_pdgId->at(kentry) == 13) id_plus.push_back(gen_muo_event);
		  if (part_pdgId->at(kentry) == -13) id_minus.push_back(gen_muo_event);


		  gen_muo_event++;
		  if ( pt_a > pt_min && eta_a > -eta_max && eta_a < eta_max ) good_1++;
		 
		}
	}

        
	
	if (gen_muo_event == 4) {
	 
	        m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
		pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
		gen_mass->Fill(m);

		if( good_1 == 4 ){
		  gen_mass_sel->Fill(m);

		}
	}
	
      gen_muo += gen_muo_event;
      rec_muo += muo_n;
      E_sum = 0;
      px_sum = 0;
      py_sum = 0;
      pz_sum = 0;
      if (muo_n == 4){
       	qualified++;

	  factor = 1;
	  
	  

	  for ( Long64_t count = 0; count<muo_n; count++ ){

	    if (muo_pdgId->at(count) == 13) muo_id_plus.push_back(count);
	    if (muo_pdgId->at(count) == -13) muo_id_minus.push_back(count);
	    muo_pt.push_back(sqrt(pow(muo_px->at(count),2)+pow(muo_py->at(count),2)));

	    E_sum += muo_E->at(count);
	    px_sum += muo_px->at(count);
	    py_sum += muo_py->at(count);
	    pz_sum += muo_pz->at(count);	
	  }
	  
	  
	  m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	  pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
	  rec_mass->Fill(m);

	  for (Long64_t check = 0; check < muo_n; check++) {
	    factor = 1;
	    pt_a = sqrt(pow(muo_px->at(check),2)+pow(muo_py->at(check),2));
	    
	    if ( muo_pz->at(check) < 0 ) factor = -1.0;
	    eta_a = log( (sqrt(pow(pt_a,2)+pow(muo_pz->at(check),2)) + (factor * muo_pz->at(check))) / pt_a) * factor;
	  
	    if ( pt_a > pt_min && eta_a > -eta_max && eta_a < eta_max ) good_2++;
	  }

	  if( good_2 == 4 ){
	    rec_mass_sel->Fill(m);
	    selected++;
	  }

      }

      ///////////////////////////////////////////////////////////////////////////////////
      //////////                        Method I                   //////////////////////
      ///////////////////////////////////////////////////////////////////////////////////
      if( gen_muo_event == 4 ){
	x1 = 0;
	x2 = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	if ( pt.at(id_plus.at(1)) > pt.at(id_plus.at(0)) ) x1 = 1;
	if ( pt.at(id_minus.at(1)) > pt.at(id_minus.at(0)) ) x2 = 1;

	E_sum = E.at(id_plus.at(x1)) + E.at(id_minus.at(x2));
	px_sum = px.at(id_plus.at(x1)) + px.at(id_minus.at(x2));
	py_sum = py.at(id_plus.at(x1)) + py.at(id_minus.at(x2));
	pz_sum = pz.at(id_plus.at(x1)) + pz.at(id_minus.at(x2));	
	m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = E.at(id_plus.at(1-x1)) + E.at(id_minus.at(1-x2));
	px_sum = px.at(id_plus.at(1-x1)) + px.at(id_minus.at(1-x2));
	py_sum = py.at(id_plus.at(1-x1)) + py.at(id_minus.at(1-x2));
	pz_sum = pz.at(id_plus.at(1-x1)) + pz.at(id_minus.at(1-x2));
	m_2 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);

	z_gen_mass_i_ab->Fill(m_1,m_2);
	z_gen_mass_i_a->Fill(m_1);
	z_gen_mass_i_b->Fill(m_2);
	
      }  


      if( muo_n == 4 && muo_id_plus.size() == 2 ){
	
	x1 = 0;
	x2 = 0;
	m_1 = 0;
	m_2 = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	if (muo_pt.at(muo_id_plus.at(1)) > muo_pt.at(muo_id_plus.at(0)) ) x1 = 1;
	if (muo_pt.at(muo_id_minus.at(1)) > muo_pt.at(muo_id_minus.at(0)) ) x2 = 1;

	E_sum = muo_E->at(muo_id_plus.at(x1)) + muo_E->at(muo_id_minus.at(x2));
	px_sum = muo_px->at(muo_id_plus.at(x1)) + muo_px->at(muo_id_minus.at(x2));
	py_sum = muo_py->at(muo_id_plus.at(x1)) + muo_py->at(muo_id_minus.at(x2));
	pz_sum = muo_pz->at(muo_id_plus.at(x1)) + muo_pz->at(muo_id_minus.at(x2));	
	m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = muo_E->at(muo_id_plus.at(1-x1)) + muo_E->at(muo_id_minus.at(1-x2));
	px_sum = muo_px->at(muo_id_plus.at(1-x1)) + muo_px->at(muo_id_minus.at(1-x2));
	py_sum = muo_py->at(muo_id_plus.at(1-x1)) + muo_py->at(muo_id_minus.at(1-x2));
	pz_sum = muo_pz->at(muo_id_plus.at(1-x1)) + muo_pz->at(muo_id_minus.at(1-x2));
	m_2 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);

	z_rec_mass_i_ab->Fill(m_1,m_2);
	z_rec_mass_i_a->Fill(m_1);
	z_rec_mass_i_b->Fill(m_2);
      }  
      
      ///////////////////////////////////////////////////////////////////////////////////
      //////////                         Method II                 //////////////////////
      ///////////////////////////////////////////////////////////////////////////////////
      b1 = 0;
      b2 = 0;
      b3 = 0;
      x1 = 0;
      x2 = 0;
      pt.clear();
      pt_a = 0;
      if( gen_muo_event == 4 ){

	for ( Long64_t count1 = 0; count1<2; count1++ ){
	  for ( Long64_t count2 = 0; count2<2; count2++ ){
	    pt_a = sqrt(pow(px.at(id_plus.at(count1))+px.at(id_minus.at(count2)),2)+pow(py.at(id_plus.at(count1))+py.at(id_minus.at(count2)),2));
	    pt.push_back(pt_a);

	  }
	}
	if ( pt.at(0) < pt.at(1) ) b1 = 1;
	if ( pt.at(2) < pt.at(3) ) b2 = 1;
	if ( pt.at(b1) < pt.at(2+b2) ) b3 = 1;
        x2 = b3;
	x1 = (b3 * b2) + ((1-b3)*b1);
	

	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = E.at(id_plus.at(x1)) + E.at(id_minus.at(x2));
	px_sum = px.at(id_plus.at(x1)) + px.at(id_minus.at(x2));
	py_sum = py.at(id_plus.at(x1)) + py.at(id_minus.at(x2));
	pz_sum = pz.at(id_plus.at(x1)) + pz.at(id_minus.at(x2));	
	m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = E.at(id_plus.at(1-x1)) + E.at(id_minus.at(1-x2));
	px_sum = px.at(id_plus.at(1-x1)) + px.at(id_minus.at(1-x2));
	py_sum = py.at(id_plus.at(1-x1)) + py.at(id_minus.at(1-x2));
	pz_sum = pz.at(id_plus.at(1-x1)) + pz.at(id_minus.at(1-x2));
	m_2 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);

	z_gen_mass_ii_ab->Fill(m_1,m_2);
	z_gen_mass_ii_a->Fill(m_1);
	z_gen_mass_ii_b->Fill(m_2);
	
      }


      if( muo_n == 4 && muo_id_plus.size() == 2 ){
	
	x1 = 0;
	x2 = 0;
	m_1 = 0;
	m_2 = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	b1 = 0;
	b2 = 0;
	b3 = 0;

	muo_pt.clear();
	pt_a = 0;
        
	for ( Long64_t count1 = 0; count1<2; count1++ ){
	  for ( Long64_t count2 = 0; count2<2; count2++ ){
	    pt_a = sqrt(pow(muo_px->at(muo_id_plus.at(count1))+muo_px->at(muo_id_minus.at(count2)),2)+pow(muo_py->at(muo_id_plus.at(count1))+muo_py->at(muo_id_minus.at(count2)),2));
	   muo_pt.push_back(pt_a);

	  }
	}
	if ( muo_pt.at(0) < muo_pt.at(1) ) b1 = 1;
	if ( muo_pt.at(2) < muo_pt.at(3) ) b2 = 1;
	if ( muo_pt.at(b1) < muo_pt.at(2+b2) ) b3 = 1;
        x2 = b3;
	x1 = (b3 * b2) + ((1-b3)*b1);
	
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = muo_E->at(muo_id_plus.at(x1)) + muo_E->at(muo_id_minus.at(x2));
	px_sum = muo_px->at(muo_id_plus.at(x1)) + muo_px->at(muo_id_minus.at(x2));
	py_sum = muo_py->at(muo_id_plus.at(x1)) + muo_py->at(muo_id_minus.at(x2));
	pz_sum = muo_pz->at(muo_id_plus.at(x1)) + muo_pz->at(muo_id_minus.at(x2));	
	m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = muo_E->at(muo_id_plus.at(1-x1)) + muo_E->at(muo_id_minus.at(1-x2));
	px_sum = muo_px->at(muo_id_plus.at(1-x1)) + muo_px->at(muo_id_minus.at(1-x2));
	py_sum = muo_py->at(muo_id_plus.at(1-x1)) + muo_py->at(muo_id_minus.at(1-x2));
	pz_sum = muo_pz->at(muo_id_plus.at(1-x1)) + muo_pz->at(muo_id_minus.at(1-x2));
	m_2 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);

	z_rec_mass_ii_ab->Fill(m_1,m_2);
	z_rec_mass_ii_a->Fill(m_1);
	z_rec_mass_ii_b->Fill(m_2);
	
      }

	
      ///////////////////////////////////////////////////////////////////////////////////
      //////////                         Method III                //////////////////////
      ///////////////////////////////////////////////////////////////////////////////////
      
      b1 = 0;
      b2 = 0;
      b3 = 0;
      x1 = 0;
      x2 = 0;
      mass.clear();
      if( gen_muo_event == 4 ){

	for ( Long64_t count1 = 0; count1<2; count1++ ){
	  for ( Long64_t count2 = 0; count2<2; count2++ ){
	    E_sum = 0;
	    px_sum = 0;
	    py_sum = 0;
	    pz_sum = 0;
	    m_1 = 0;
	    E_sum = E.at(id_plus.at(count1)) + E.at(id_minus.at(count2));
	    px_sum = px.at(id_plus.at(count1)) + px.at(id_minus.at(count2));
	    py_sum = py.at(id_plus.at(count1)) + py.at(id_minus.at(count2));
	    pz_sum = pz.at(id_plus.at(count1)) + pz.at(id_minus.at(count2));	
	    m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	    mass.push_back(m_1);

	  }
	}
	if ( mass.at(0) < mass.at(1) ) b1 = 1;
	if ( mass.at(2) < mass.at(3) ) b2 = 1;
	if ( mass.at(b1) < mass.at(2+b2) ) b3 = 1;
	x2 = b3;
	x1 = (b3 * b2) + ((1-b3)*b1);
	
	m_1 = mass.at( ( 2 * x2 ) + x1 );
	m_2 = mass.at( ( 2 * ( 1 - x2 ) + ( 1 - x1 ) ) );
	
	if ( m_1 > 50 && m_1 < 106 && m_2 > 12 && m_2 < 115 ) z_gen_mass_iii_ab->Fill(m_1,m_2),  z_gen_mass_iii_a->Fill(m_1), z_gen_mass_iii_b->Fill(m_2);
        
      }

      
      if( muo_n == 4 && muo_id_plus.size() == 2 ){
	
	x1 = 0;
	x2 = 0;
	m_1 = 0;
	m_2 = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	b1 = 0;
	b2 = 0;
	b3 = 0;

	mass.clear();
	m_1 = 0;
        
	for ( Long64_t count1 = 0; count1<2; count1++ ){
	  for ( Long64_t count2 = 0; count2<2; count2++ ){
	    E_sum = 0;
	    px_sum = 0;
	    py_sum = 0;
	    pz_sum = 0;
	    m_1 = 0;
	    E_sum = muo_E->at(muo_id_plus.at(count1)) + muo_E->at(muo_id_minus.at(count2));
	    px_sum = muo_px->at(muo_id_plus.at(count1)) + muo_px->at(muo_id_minus.at(count2));
	    py_sum = muo_py->at(muo_id_plus.at(count1)) + muo_py->at(muo_id_minus.at(count2));
	    pz_sum = muo_pz->at(muo_id_plus.at(count1)) + muo_pz->at(muo_id_minus.at(count2));	
	    m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	    mass.push_back(m_1);

	  }
	}
	if ( mass.at(0) < mass.at(1) ) b1 = 1;
	if ( mass.at(2) < mass.at(3) ) b2 = 1;
	if ( mass.at(b1) < mass.at(2+b2) ) b3 = 1;
        x2 = b3;
	x1 = (b3 * b2) + ((1-b3)*b1);
	
	
        m_1 = mass.at( ( 2 * x2 ) + x1 );
	m_2 = mass.at( ( 2 * ( 1 - x2 ) + ( 1 - x1 ) ) );
	
	if ( m_1 > 50 && m_1 < 106 && m_2 > 12 && m_2 < 115 ) z_rec_mass_iii_ab->Fill(m_1,m_2), z_rec_mass_iii_a->Fill(m_1), z_rec_mass_iii_b->Fill(m_2);


      }
     
      
   }
   cout << " selection rate: " << (selected/qualified)*100.0 <<  "%" << endl;
}
