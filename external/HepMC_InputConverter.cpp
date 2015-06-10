#include "HepMC_InputConverter.h"
using namespace AcerDet::external;

#include <cstdio>

Vector4f vec4(const HepMC::FourVector& v) {
	return Vector4f(v.x(), v.y(), v.z(), v.e());
}

ParticleType HepMC_InputConverter::getParticleType(int code) {
	if (code == -4 || code == 4) return PT_CJET;
	if (code == -5 || code == 5) return PT_BJET;
	if (code == -11 || code == 11) return PT_ELECTRON;
	if (code == -12 || code == 12) return PT_NEUTRINO_ELE;
	if (code == -13 || code == 13) return PT_MUON;
	if (code == -14 || code == 14) return PT_NEUTRINO_MUO;
	if (code == -15 || code == 15) return PT_TAU;
	if (code == -16 || code == 16) return PT_NEUTRINO_TAU;
	if (code == -22 || code == 22) return PT_PHOTON;
	if (code == 23) return PT_BOSON_Z;
	if (code == 24) return PT_BOSON_W;
	if (code == 25) return PT_BOSON_H;
	
	return PT_UNKNOWN;
}

ParticleStatus HepMC_InputConverter::getParticleStatus(HepMC::GenParticle* gpart) {
	if (gpart->is_beam()) return PS_BEAM; 	     // status_code == 4
	if (gpart->is_undecayed()) return PS_FINAL;  // status_code == 1 -> final state
	if (gpart->has_decayed()) return PS_DECAYED; // status_code == 2 -> before hadronization
	if (gpart->status() == 3) return PS_HISTORY; // documentation line -> history
	if (abs(gpart->status()) >= 30) return PS_CASCADE_QUARK;
	if (20 <= abs(gpart->status()) && abs(gpart->status()) < 30) return PS_HP_QUARK;
	return PS_NULL;
}

Int32_t extractMother(HepMC::GenVertex* ptr) {
	if (ptr == NULL)
		return -1;
	
	const HepMC::GenParticle* p = *(ptr->particles_in_const_begin());
	if (p == NULL)
		return -1;
	return p->barcode();
}

Int32_t extractDaughter1(HepMC::GenVertex* ptr) {
	if (ptr == NULL)
		return -1;
	
	const HepMC::GenParticle* p = *(ptr->particles_out_const_begin());
	if (p == NULL)
		return -1;
	return p->barcode();
}

Int32_t extractDaughter2(HepMC::GenVertex* ptr) {
	if (ptr == NULL)
		return -1;
	
	if (*(ptr->particles_out_const_begin()) == NULL)
		return -1;
	
	Int32_t code = -1;
	for (vector<HepMC::GenParticle*>::const_iterator i = ptr->particles_out_const_begin(); i != ptr->particles_out_const_end(); i++)
		code = (*i)->barcode();
	return code;
}

InputRecord HepMC_InputConverter::convert( const GenEvent& event ) {
	vector<Particle> parts;
	
	for( GenEvent::particle_const_iterator iter = event.particles_begin(); iter != event.particles_end(); ++iter ) {
		GenParticle* gpart = *iter;
		
		core::Particle part;
		
		// particle tree hierarchy
		part.barcode = gpart->barcode();
		HepMC::GenVertex* prod = gpart->production_vertex();
		part.mother = extractMother(prod);
		HepMC::GenVertex* decay = gpart->end_vertex();
		part.daughters = make_pair(extractDaughter1(decay), extractDaughter2(decay));

		// state (named & id)
		part.status = getParticleStatus(gpart);
		part.statusID = gpart->status();
		
		// type (named & id)
		part.type = getParticleType(gpart->pdg_id());
		part.pdg_id = gpart->pdg_id();

		// check conversion from Pythia8 event record
		//    printf("barcode=%d, status=%d,  pdgid=%d, mother=%d, daughter1=%d, daughter2=%d\n", 
		//    part.barcode, part.statusID, part.pdg_id, part.mother, extractDaughter1(decay), extractDaughter2(decay)); 
		
		// momentum as Vector4
		part.momentum = vec4(gpart->momentum());
		
		// production vertex (optional) as Vector4
		if (prod != NULL) {
			part.production = vec4(prod->position());
		}
		
		parts.push_back(part);
		//cout << part; // to delete in release!
	}
	
	return InputRecord(parts);
}
