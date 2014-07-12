#include "InputReader.h"
using namespace AcerDet::io;

#include <cstdio>

InputRecord InputReader::computeEvent( const GenEvent& event ) {
	vector<Particle> parts;
	
	for( GenEvent::particle_const_iterator iter = event.particles_begin(); iter != event.particles_end(); ++iter ) {
		GenParticle* gpart = *iter;
		
		//printf ("Particle:\n");
		//printf ("\tid: %d\n", gpart->pdg_id());
		//printf ("\tstatus: %d\n", gpart->status());
		//printf ("\tpolarization phi: %f\n", gpart->polarization().phi());
		//printf ("\tpolarization theta: %f\n", gpart->polarization().theta());
		//printf ("\tgenerated mass: %f\n", gpart->generated_mass());
		//printf ("\tbarcode: %d\n", gpart->barcode());
		//printf ("\tis_undecayed: %s\n", gpart->is_undecayed() ? "Yes" : "No");
		//printf ("\thas_decayed: %s\n", gpart->has_decayed() ? "Yes" : "No");
		//printf ("\tis_beam: %s\n", gpart->is_beam() ? "Yes" : "No");
		//printf ("\tpx= %f\n", gpart->momentum().px());
		//printf ("\tpy= %f\n", gpart->momentum().py());
		//printf ("\tpz= %f\n", gpart->momentum().pz());
		//printf ("\te= %f\n", gpart->momentum().e());
		//gpart->print();
		//printf ("\n");
		
		core::Particle part;
		
		part.id = gpart->pdg_id();
		part.mother = 0;
		part.daughters = make_pair(0,0);
		
		part.state = PS_NULL;
		if (gpart->is_beam()) part.state = PS_BEAM; 			// status_code == 4
		else if (gpart->is_undecayed()) part.state = PS_FINAL;  // status_code == 1 -> final state
		else if (gpart->has_decayed()) part.state = PS_DECAYED; // status_code == 2 -> before hadronization
		else if (gpart->status() == 3) part.state = PS_HISTORY; // documentation line -> history
		
		const HepMC::FourVector& momentum = gpart->momentum();
		part.px = momentum.px();
		part.py = momentum.py();
		part.pz = momentum.pz();
		part.e = momentum.e();
		
		const HepMC::FourVector& prod = gpart->momentum();//gpart->production_vertex()->position();
		part.x = prod.x();
		part.y = prod.y();
		part.z = prod.z();
		
		part.phi = gpart->polarization().phi();
		part.theta = gpart->polarization().theta();
		
		parts.push_back(part);
		part.print();
	}
	
	InputRecord record(parts);
	return record;
}
