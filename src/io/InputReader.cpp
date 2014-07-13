#include "InputReader.h"
using namespace AcerDet::io;

#include <cstdio>

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
		
		part.id = gpart->barcode();
		part.mother = -1;
		part.daughters = make_pair(-1,-1);
		
		part.state = getParticleStatus(gpart->status());
		
		const HepMC::FourVector& momentum = gpart->momentum();
		part.px = momentum.px();
		part.py = momentum.py();
		part.pz = momentum.pz();
		part.e = momentum.e();
		
		HepMC::GenVertex* prod = gpart->production_vertex();
		if (prod != NULL) {
			const HepMC::FourVector& pv = prod->position();
			part.prod_x = pv.x();
			part.prod_y = pv.y();
			part.prod_z = pv.z();
			part.prod_time = pv.e();
		}
		
		part.phi = gpart->polarization().phi();
		part.theta = gpart->polarization().theta();
		
		parts.push_back(part);
		part.print();
	}
	
	InputRecord record(parts);
	return record;
}
