#include "InputReader.h"
using namespace AcerDet::io;

#include <cstdio>

Vector4f vec4(const HepMC::FourVector& v) {
	return Vector4f(v.x(), v.y(), v.z(), v.e());
}

ParticleType InputReader::getParticleType(int hepmc_code) {
	int code = abs(hepmc_code);
	PDGcode pdg = (PDGcode)code;

	switch(pdg) {
	case QUARK_C: return PT_CJET;
	case QUARK_B: return PT_BJET;
	case LEPT_ELECTRON: return PT_ELECTRON;
	case LEPT_MUON: return PT_MUON;
	case LEPT_TAU: return PT_TAU;
	default: return PT_UNKNOWN;
	}
}

ParticleState InputReader::getParticleStatus(int hepmc_status) {
	switch(hepmc_status) {
	case 1: return PS_FINAL;
	case 2: return PS_DECAYED;
	case 3: return PS_HISTORY;
	case 4: return PS_BEAM;
	}
	return PS_NULL;
	/*
	if (gpart->is_beam()) part.state = PS_BEAM; 			// status_code == 4
	else if (gpart->is_undecayed()) part.state = PS_FINAL;  // status_code == 1 -> final state
	else if (gpart->has_decayed()) part.state = PS_DECAYED; // status_code == 2 -> before hadronization
	else if (gpart->status() == 3) part.state = PS_HISTORY; // documentation line -> history
	*/
}

InputRecord InputReader::computeEvent( const GenEvent& event ) {
	vector<Particle> parts;
	
	for( GenEvent::particle_const_iterator iter = event.particles_begin(); iter != event.particles_end(); ++iter ) {
		GenParticle* gpart = *iter;
		
		core::Particle part;
		
		// particle tree hierarchy
		part.id = gpart->barcode();
		part.mother = -1;
		part.daughters = make_pair(-1,-1);
		
		// state (named & id)
		part.state = getParticleStatus(gpart->status());
		part.stateID = gpart->status();
		
		// type (named & id)
		part.type = getParticleType(gpart->pdg_id());
		part.typeID = gpart->pdg_id();
		
		// momentum as Vector4
		part.momentum = vec4(gpart->momentum());
		
		// production vertex (optional) as Vector4
		HepMC::GenVertex* prod = gpart->production_vertex();
		if (prod != NULL) {
			part.production = vec4(prod->position());
		}
		
		// angles
		//part.phi = gpart->polarization().phi();
		//part.theta = gpart->polarization().theta();
		
		parts.push_back(part);
		// cout << part; // to delete in release!
	}
	
	return InputRecord(parts);
}
