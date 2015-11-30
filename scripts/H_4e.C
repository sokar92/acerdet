#define H_4e_cxx
#include "H_4e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void H_4e::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F *gen_mass = new TH1F("shist12366" , "Electron: Higgs mass HP", 100, 80, 128 );
   TH1F *gen_mass_sel = new TH1F("shist12367" , "Electron: Higgs mass HP SEL", 100, 80, 128 );
   TH1F *rec_mass = new TH1F("shist12356" , "Electron: Higgs mass ISO", 100, 80, 140 );
   TH1F *rec_mass_sel = new TH1F("shist12357" , "Electron: Higgs mass ISO SEL", 100, 80, 140 );
 
   TH1F *z_gen_mass_i_a = new TH1F("shist123631" , "Electron: Higher Z boson mass HP I", 100, 0, 140 );
   TH1F *z_gen_mass_i_b = new TH1F("shist123641" , "Electron: Lower Z boson mass HP I", 100, 0, 140 );
   TH2F *z_gen_mass_i_ab = new TH2F("shist123131" , "Electron: Z bosons masses HP I", 100, 0, 140, 100, 0, 140 );
   
   TH1F *z_gen_mass_ii_a = new TH1F("shist123632" , "Electron: Higher Z boson mass HP II", 100, 0, 140 );
   TH1F *z_gen_mass_ii_b = new TH1F("shist123642" , "Electron: Lower Z boson mass HP II", 100, 0, 140 );
   TH2F *z_gen_mass_ii_ab = new TH2F("shist123132" , "Electron: Z bosons masses HP II", 100, 0, 140, 100, 0, 140 );
   
   TH1F *z_gen_mass_iii_a = new TH1F("shist123633" , "Electron: Higher Z boson mass HP III", 100, 0, 140 );
   TH1F *z_gen_mass_iii_b = new TH1F("shist123643" , "Electron: Lower Z boson mass HP III", 100, 0, 140 );
   TH2F *z_gen_mass_iii_ab = new TH2F("shist123133" , "Electron: Z bosons masses HP III", 100, 0, 140, 100, 0, 140 );
 
   TH1F *z_rec_mass_i_a = new TH1F("shist123531" , "Electron: Higher Z boson mass ISO I", 100, 0, 140 );
   TH1F *z_rec_mass_i_b = new TH1F("shist123541" , "Electron: Lower Z boson mass ISO I", 100, 0, 140 );
   TH2F *z_rec_mass_i_ab = new TH2F("shist123141" , "Electron: Z bosons masses ISO I", 100, 0, 140, 100, 0, 140 );
  
   TH1F *z_rec_mass_ii_a = new TH1F("shist123532" , "Electron: Higher Z boson mass ISO II", 100, 0, 140 );
   TH1F *z_rec_mass_ii_b = new TH1F("shist123542" , "Electron: Lower Z boson mass ISO II", 100, 0, 140 );
   TH2F *z_rec_mass_ii_ab = new TH2F("shist123142" , "Electron: Z bosons masses ISO II", 100, 0, 140, 100, 0, 140 );

   TH1F *z_rec_mass_iii_a = new TH1F("shist123533" , "Electron: Higher Z boson mass ISO III", 100, 0, 140 );
   TH1F *z_rec_mass_iii_b = new TH1F("shist123543" , "Electron: Lower Z boson mass ISO III", 100, 0, 140 );
   TH2F *z_rec_mass_iii_ab = new TH2F("shist123143" , "Electron: Z bosons masses ISO III", 100, 0, 140, 100, 0, 140 );

   vector<float> E;
   vector<float> px;
   vector<float> py;
   vector<float> pz;
   vector<float> pt;
   vector<float> ele_pt;
   vector<float> eta;
   vector<float> id_plus;
   vector<float> id_minus;
   vector<float> ele_id_plus;
   vector<float> ele_id_minus;
   vector<float> mass;

   float  m = 0, m_1 = 0, m_2 = 0 , E_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0, pt_0 = 0, pt_a = 0, pt_b = 0, eta_a = 0, eta_b = 0, factor = 0, qualified = 0, x1 = 0, x2 = 0;
	
   float  pt_min = 7, eta_max = 2.47;
   Long64_t nbytes = 0, nb = 0, gen_ele_event = 0, gen_ele = 0, rec_ele = 0, part_numb = 0, selected = 0, good_1 = 0, b1 = 0, b2 = 0, b3 = 0,
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
	gen_ele_event = 0;
	part_numb  = part_pdgId->size();
	pt.clear();
	eta.clear();
	E.clear();
	px.clear();
	py.clear();
	pz.clear();
	id_plus.clear();
	id_minus.clear();
	ele_id_plus.clear();
	ele_id_minus.clear();
	ele_pt.clear();

	for (Long64_t kentry = 0; kentry < part_pdgId->size(); kentry++){
	
	  if(( part_pdgId->at(kentry) == 11 || part_pdgId->at(kentry) == -11 ) && part_mother_pdgId->at(kentry) == 23 ){	
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
		  if (part_pdgId->at(kentry) == 11) id_plus.push_back(gen_ele_event);
		  if (part_pdgId->at(kentry) == -11) id_minus.push_back(gen_ele_event);


		  gen_ele_event++;
		  if ( pt_a > pt_min && eta_a > -eta_max && eta_a < eta_max ) good_1++;
	
		}
	}

        
	
	if (gen_ele_event == 4) {
	 
	        m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
		pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
		gen_mass->Fill(m);
        
	
		if( good_1 == 4 ){
		  gen_mass_sel->Fill(m);
		}
	}
	
	

      gen_ele += gen_ele_event;
      rec_ele += ele_n;
      E_sum = 0;
      px_sum = 0;
      py_sum = 0;
      pz_sum = 0;
      if (ele_n == 4){
       	qualified++;

	  factor = 1;
	  
	  

	  for ( Long64_t count = 0; count<ele_n; count++ ){

	    if (ele_pdgId->at(count) == 11) ele_id_plus.push_back(count);
	    if (ele_pdgId->at(count) == -11) ele_id_minus.push_back(count);
	    ele_pt.push_back(sqrt(pow(ele_px->at(count),2)+pow(ele_py->at(count),2)));

	    E_sum += ele_E->at(count);
	    px_sum += ele_px->at(count);
	    py_sum += ele_py->at(count);
	    pz_sum += ele_pz->at(count);	
	  }
	  
	  
	  m = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	  pt_0 = sqrt(pow(px_sum,2)+pow(py_sum,2));
	  rec_mass->Fill(m);



	  for (Long64_t check = 0; check < ele_n; check++) {
	    factor = 1;
	    pt_a = sqrt(pow(ele_px->at(check),2)+pow(ele_py->at(check),2));
	    
	    if ( ele_pz->at(check) < 0 ) factor = -1.0;
	    eta_a = log( (sqrt(pow(pt_a,2)+pow(ele_pz->at(check),2)) + (factor * ele_pz->at(check))) / pt_a) * factor;
	  
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
      if( gen_ele_event == 4 ){
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
      


      if( ele_n == 4 && ele_id_plus.size() == 2 ){
	
	x1 = 0;
	x2 = 0;
	m_1 = 0;
	m_2 = 0;
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	if (ele_pt.at(ele_id_plus.at(1)) > ele_pt.at(ele_id_plus.at(0)) ) x1 = 1;
	if (ele_pt.at(ele_id_minus.at(1)) > ele_pt.at(ele_id_minus.at(0)) ) x2 = 1;

	E_sum = ele_E->at(ele_id_plus.at(x1)) + ele_E->at(ele_id_minus.at(x2));
	px_sum = ele_px->at(ele_id_plus.at(x1)) + ele_px->at(ele_id_minus.at(x2));
	py_sum = ele_py->at(ele_id_plus.at(x1)) + ele_py->at(ele_id_minus.at(x2));
	pz_sum = ele_pz->at(ele_id_plus.at(x1)) + ele_pz->at(ele_id_minus.at(x2));	
	m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = ele_E->at(ele_id_plus.at(1-x1)) + ele_E->at(ele_id_minus.at(1-x2));
	px_sum = ele_px->at(ele_id_plus.at(1-x1)) + ele_px->at(ele_id_minus.at(1-x2));
	py_sum = ele_py->at(ele_id_plus.at(1-x1)) + ele_py->at(ele_id_minus.at(1-x2));
	pz_sum = ele_pz->at(ele_id_plus.at(1-x1)) + ele_pz->at(ele_id_minus.at(1-x2));
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
      if( gen_ele_event == 4 ){

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


      if( ele_n == 4 && ele_id_plus.size() == 2 ){
	
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

	ele_pt.clear();
	pt_a = 0;
        
	for ( Long64_t count1 = 0; count1<2; count1++ ){
	  for ( Long64_t count2 = 0; count2<2; count2++ ){
	    pt_a = sqrt(pow(ele_px->at(ele_id_plus.at(count1))+ele_px->at(ele_id_minus.at(count2)),2)+pow(ele_py->at(ele_id_plus.at(count1))+ele_py->at(ele_id_minus.at(count2)),2));
	   ele_pt.push_back(pt_a);

	  }
	}
	if ( ele_pt.at(0) < ele_pt.at(1) ) b1 = 1;
	if ( ele_pt.at(2) < ele_pt.at(3) ) b2 = 1;
	if ( ele_pt.at(b1) < ele_pt.at(2+b2) ) b3 = 1;
        x2 = b3;
	x1 = (b3 * b2) + ((1-b3)*b1);
	
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = ele_E->at(ele_id_plus.at(x1)) + ele_E->at(ele_id_minus.at(x2));
	px_sum = ele_px->at(ele_id_plus.at(x1)) + ele_px->at(ele_id_minus.at(x2));
	py_sum = ele_py->at(ele_id_plus.at(x1)) + ele_py->at(ele_id_minus.at(x2));
	pz_sum = ele_pz->at(ele_id_plus.at(x1)) + ele_pz->at(ele_id_minus.at(x2));	
	m_1 = sqrt(E_sum*E_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum);
	E_sum = 0;
	px_sum = 0;
	py_sum = 0;
	pz_sum = 0;
	E_sum = ele_E->at(ele_id_plus.at(1-x1)) + ele_E->at(ele_id_minus.at(1-x2));
	px_sum = ele_px->at(ele_id_plus.at(1-x1)) + ele_px->at(ele_id_minus.at(1-x2));
	py_sum = ele_py->at(ele_id_plus.at(1-x1)) + ele_py->at(ele_id_minus.at(1-x2));
	pz_sum = ele_pz->at(ele_id_plus.at(1-x1)) + ele_pz->at(ele_id_minus.at(1-x2));
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
      if( gen_ele_event == 4 ){

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

      
      if( ele_n == 4 && ele_id_plus.size() == 2 ){
	
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
	    E_sum = ele_E->at(ele_id_plus.at(count1)) + ele_E->at(ele_id_minus.at(count2));
	    px_sum = ele_px->at(ele_id_plus.at(count1)) + ele_px->at(ele_id_minus.at(count2));
	    py_sum = ele_py->at(ele_id_plus.at(count1)) + ele_py->at(ele_id_minus.at(count2));
	    pz_sum = ele_pz->at(ele_id_plus.at(count1)) + ele_pz->at(ele_id_minus.at(count2));	
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
