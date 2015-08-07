#include "Test.h"
#include <cstdio>


using namespace AcerDet::test;

Test::Test() {}

Test::~Test() {}

void Test::Fill2e( Int32_t idhist, Real64_t weight, IHistogramManager* histoManager ) {
	
	Int32_t counter = 0;
	Real64_t mass_1 = 0, mass_2 = 0;
	
	if (HARDelectron.size() == 2){
		counter++;
		allparticlemomenta(HARDelectron);
		mass_1 = mass();
		
		PutHisto( mass(), mass2(), 61, idhist, weight, histoManager );
		
	}
	
	if (ISOelectron.size() == 2){
		counter++;
		allobjectmomenta(ISOelectron);
		mass_2 = mass();
		
		PutHisto( mass(), mass2(), 51, idhist, weight, histoManager );
		
	}

	if ( counter == 2 && mass_1 > 0 && mass_2 > 0 ) {
		histoManager
			->insert(idhist+41, mass_1 - mass_2, weight);
	}

}

void Test::Fill2gamma( Int32_t idhist, Real64_t weight, IHistogramManager* histoManager ) {
		
	if (HARDphoton.size() == 2){
		allparticlemomenta(HARDphoton);
		PutHisto( mass(), mass2(), 61, idhist, weight, histoManager );		
	}
	
	if (ISOphoton.size() == 2){
		allobjectmomenta(ISOphoton);	
		PutHisto( mass(), mass2(), 51, idhist, weight, histoManager );	
	}
	if ( ISOphoton.size() > 0 ){
		for (int j=0; j<ISOphoton.size(); ++j) {
			histoManager
				->insert(idhist+71, ISOphoton.at(j).pT, weight);
		}
	}
	if ( HARDphoton.size() > 0 ){
		for (int j=0; j<HARDphoton.size(); ++j) {
			histoManager
				->insert(idhist+81, HARDphoton.at(j).pT(), weight);
		}
	}
}

void Test::Fill4part( vector<Particle> part, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager ) {
	
	vector<Int32_t> k;
	Dividepart ( part );
	if (part.size() == 4 && particle1.size() == 2 && particle2.size() == 2 ){
		allparticlemomenta(part);
		
		PutHisto( mass(), mass2(), 66, idhist, weight, histoManager );
		groupparticlemasses();
		k = Findhighest();

		histoManager
			->insert(idhist+63+k.at(0), storage1.at(k.at(1)), weight);
		histoManager
			->insert(idhist+64-k.at(0), storage2.at(k.at(2)), weight);
	}
}

void Test::Fill4object( vector<ObjectData> newObject, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager ) {
	
	vector<Int32_t> k;
	Divideobject( newObject );
	if (newObject.size() == 4 && object1.size() == 2 && object2.size() == 2){

		allobjectmomenta(newObject);
		PutHisto( mass(), mass2(), 56, idhist, weight, histoManager );
			
		groupobjectmasses();
		k = Findhighest();
	
		histoManager
			->insert(idhist+53+k.at(0), storage1.at(k.at(1)), weight);
		histoManager
			->insert(idhist+54-k.at(0), storage2.at(k.at(2)), weight);
	}

}

void Test::Filljets( vector<Particle> quark, vector<JetData> jet, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager ) {

	vector<Int32_t> k;
	
	if ( quark.size() == 2 ){
	
		allparticlemomenta(quark);
		PutHisto( mass(), mass2(), 61, idhist, weight, histoManager );
		Dividepart ( quark );
		
		if ( jet.size() > 1 ){
			
			groupdeltaR( jet );
			k = Findlowest();
			
			if ( storage1.at(k.at(1)) < 0.4 && storage2.at(k.at(2)) < 0.4 ){
			
				selectjetmomenta(jet,k.at(1),k.at(2));
				PutHisto( mass(), mass2(), 51, idhist, weight, histoManager );	
			}
				
		}
	}
}	


void Test::AddMomenta( Real64_t pT, Real64_t eta, Real64_t phi ){
		Real64_t px = pT * cos(phi);
		Real64_t py = pT * sin(phi);
		Real64_t pz = pT * sinh(eta);
		Real64_t E  = pT * cosh(eta);
		Vector4f mom( px, py, pz, E );
		momentum += mom;
}

vector<Int32_t> Test::Findlowest() {
	vector<Int32_t> results;
	Int32_t winner1 = 0, winner2 = 0, run_up1 = 1, run_up2 = 1, group = 0; 
	for ( Int32_t count = 1; count<storage1.size(); count++ ){
	   if ( storage1.at(winner1) > storage1.at(count) ) run_up1 = winner1, winner1 = count;
	   if ( storage2.at(winner2) > storage2.at(count) ) run_up2 = winner2, winner2 = count;
	  
	 }
	
	 if ( winner1 == winner2 ) {
		 
		if ( winner1 == 0 && storage1.size() > 2){
			
			for ( Int32_t count = 2; count < storage1.size(); count++ ) {
				if ( storage1.at(run_up1) > storage1.at(count) ) run_up1 = count;
				if ( storage2.at(run_up2) > storage2.at(count) ) run_up2 = count;
			}
		}
		
		if ( storage1.at(winner1) <= storage2.at(winner2) ) winner2 = run_up2;
		if ( storage1.at(winner1) > storage2.at(winner2) ) winner1 = run_up1;
	 }
	 
	 if ( storage1.at(winner1) > storage2.at(winner2) ) group = 1;
	 results.push_back(group);
	 results.push_back(winner1);
	 results.push_back(winner2);
	 
	 return results;
 }
 
vector<Int32_t> Test::Findhighest() {
	vector<Int32_t> results;
	Int32_t winner1 = 0, winner2 = 0, run_up1 = 1, run_up2 = 1, group = 0; 
	
	for ( Int32_t count = 1; count<storage1.size(); count++ ){
	   if ( storage1.at(winner1) < storage1.at(count) ) run_up1 = winner1, winner1 = count;
	   if ( storage2.at(winner2) < storage2.at(count) ) run_up2 = winner2, winner2 = count;
	}
	
	if ( winner1 == winner2 ) {
		if ( winner1 == 0 && storage1.size() > 2){
			
			for ( Int32_t count = 2; count < storage1.size(); count++ ) {
				if ( storage1.at(run_up1) < storage1.at(count) ) run_up1 = count;
				if ( storage2.at(run_up2) < storage2.at(count) ) run_up2 = count;
			}
		}
		 
		if ( storage1.at(winner1) >= storage2.at(winner2) ) winner2 = run_up2;
		if ( storage1.at(winner1) < storage2.at(winner2) ) winner1 = run_up1;
	}
	 
	if ( storage1.at(winner1) < storage2.at(winner2) ) group = 1;
	results.push_back(group);
	results.push_back(winner1);
	results.push_back(winner2);
	 
	return results;
}
 	
void Test::groupobjectmasses(){
	
	clearstorage();
	Vector4f momzero;
	
	if ( ( object1.size() == 2 ) && ( object2.size() == 2 ) ){
		
		for ( Int32_t j = 0; j<object2.size(); j++ ){
			momentum = momzero;
			AddMomenta(object1.at(0).pT, object1.at(0).eta, object1.at(0).phi);
			AddMomenta(object2.at(j).pT, object2.at(j).eta, object2.at(j).phi);
			if ( mass2() > 0 ) storage1.push_back(mass());
			
			momentum = momzero;
			AddMomenta(object1.at(1).pT, object1.at(1).eta, object1.at(1).phi);
			AddMomenta(object2.at(j).pT, object2.at(j).eta, object2.at(j).phi);
			if ( mass2() > 0 ) storage2.push_back(mass());
		}
	}
	
}

void Test::groupparticlemasses(){
	
	clearstorage();
	Vector4f momzero;
	
	if ( ( particle1.size() == 2 ) && ( particle2.size() == 2 ) ){
		
		for ( Int32_t j = 0; j<particle2.size(); j++ ){
			momentum = momzero;
			AddMomenta(particle1.at(0).pT(), particle1.at(0).getEta(), particle1.at(0).getPhi());
			AddMomenta(particle2.at(j).pT(), particle2.at(j).getEta(), particle2.at(j).getPhi());
			if ( mass2() > 0 ) storage1.push_back(mass());
			
			momentum = momzero;
			AddMomenta(particle1.at(1).pT(), particle1.at(1).getEta(), particle1.at(1).getPhi());
			AddMomenta(particle2.at(j).pT(), particle2.at(j).getEta(), particle2.at(j).getPhi());
			if ( mass2() > 0 ) storage2.push_back(mass());
		}
	}
	
}

Real64_t Test::deltaR( Particle part, JetData jet ) {
	
	Real64_t DPhi = abs(part.getPhi() - jet.phi_rec); 
	if (DPhi > PI) DPhi -= 2*PI;
	return sqrt(pow(DPhi,2)+pow(part.getEta()-jet.eta_rec,2));
	
}

void Test::groupdeltaR(vector<JetData> jet){
	
	clearstorage();
	Vector4f momzero;
	
	if ( ( particle1.size() == 1 ) && ( particle2.size() == 1 ) ){
		
		for ( Int32_t j = 0; j<jet.size(); j++ ){
			storage1.push_back( deltaR( particle1.at(0), jet.at(j) ) );
			storage2.push_back( deltaR( particle2.at(0), jet.at(j) ) );
		}
	}
	
}

void Test::PutHisto( Real64_t value, Real64_t condition, Int32_t id, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager ) {
	
	if ( condition > 0 ){
		histoManager
			->insert(idhist+id, value, weight);
	}
}

void Test::selectjetmomenta( vector<JetData> jet, Int32_t it_1, Int32_t it_2 ){
		Vector4f momzero;
		momentum = momzero;
		AddMomenta(jet.at(it_1).pT, jet.at(it_1).eta_rec, jet.at(it_1).phi_rec);
		AddMomenta(jet.at(it_2).pT, jet.at(it_2).eta_rec, jet.at(it_2).phi_rec);
}

void Test::allobjectmomenta(vector<ObjectData> newObject){
		Vector4f momzero;
		momentum = momzero;
		if ( newObject.size() > 0 ){
			for ( Int32_t i = 0; i<newObject.size(); i++){
				AddMomenta(newObject.at(i).pT, newObject.at(i).eta, newObject.at(i).phi);
			}
		}
}

void Test::allparticlemomenta(vector<Particle> part){
		Vector4f momzero;
		momentum = momzero;
		if ( part.size() > 0 ){
			for ( Int32_t i = 0; i<part.size(); i++){
				AddMomenta(part.at(i).pT(), part.at(i).getEta(), part.at(i).getPhi());
			}
		}
}

void Test::PutHARD( Particle part ){
		
		if ( part.pdg_id == 1 || part.pdg_id == -1 ) uquark.push_back( part );
		else if ( part.type == PT_BJET ) bquark.push_back( part );
		else if ( part.type == PT_CJET ) cquark.push_back( part );
		else if ( part.type == PT_MUON ) HARDmuon.push_back( part );
		else if ( part.type == PT_PHOTON ) HARDphoton.push_back( part );
		else if ( part.type == PT_ELECTRON ) HARDelectron.push_back( part );
		
}

void Test::Dividepart( vector<Particle> part ){
	particle1.clear();
	particle2.clear();
	if ( part.size() > 0 ) {
		for ( Int32_t j = 0; j<part.size(); j++){
			if ( part.at(j).pdg_id > 0 ) particle1.push_back( part.at(j) );
			else if ( part.at(j).pdg_id < 0 ) particle2.push_back( part.at(j) );
		}
	}
}

void Test::Divideobject( vector<ObjectData> newObject ){
	object1.clear();
	object2.clear();
	if ( newObject.size() > 0 ) {
		for ( Int32_t j = 0; j<newObject.size(); j++){
			if ( newObject.at(j).pdg_id > 0 ) object1.push_back( newObject.at(j) );
			else if ( newObject.at(j).pdg_id < 0 ) object2.push_back( newObject.at(j) );
		}
	}
}				

void Test::Dividejets(){
	if ( jet.size() > 0 ) {
		for ( Int32_t j = 0; j<jet.size(); j++){
			if ( jet.at(j).type == C_JET) cjet.push_back( jet.at(j) );
			else if ( jet.at(j).type == B_JET) bjet.push_back( jet.at(j) );
			else if ( jet.at(j).type == UNKNOWN ) ujet.push_back( jet.at(j) );
		}
	}
}
