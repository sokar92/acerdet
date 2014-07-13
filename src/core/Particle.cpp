#include "Particle.h"
using namespace AcerDet::core;

#define ParticleNull (-1)

Particle::Particle() : 
	state(PS_NULL), 
	stateID(-1),
	type(PT_UNKNOWN),
	typeID(0), 
	px(0), py(0), pz(0), e(0), 
	phi(0), theta(0),
	prod_x(0), prod_y(0), prod_z(0), prod_time(0), 
	id(-1), mother(-1), 
	daughters(make_pair(-1,-1)) 
{}

Bool_t Particle::hasMother() const { 
	return mother != ParticleNull; 
}

Bool_t Particle::hasDaughter() const { 
	return daughters.first != ParticleNull; 
}

Int32_t Particle::daughtersCount() const {
	if (hasDaughter())
		return 1 + daughters.second - daughters.first;
	return 0;
}

Bool_t Particle::isStable() const {
	return 0 < stateID && stateID <= 10;
}

Bool_t Particle::isBeam() const {
	return state == PS_BEAM;
}

Bool_t Particle::isDecayed() const {
	return state == PS_DECAYED;
}

void Particle::print() const {
	printf ("Particle:\n");
	printf ("\tID: %d\n", id);
	printf ("\tType: %s\n", getTypeName().c_str());
	printf ("\tState: %s\n", getStateName().c_str());
	printf ("\tPolarization (phi, theta): (%f, %f)\n", phi, theta);
	printf ("\tMomentum (px,py,pz,e) = (%f, %f, %f, %f)\n", px, py, pz, e);
	printf ("\tProduction (x,y,z,t) = (%f, %f, %f, %f)\n", prod_x, prod_y, prod_z, prod_time);
	printf ("\n");
}

string Particle::getTypeName() const {
	switch(type) {
	case PT_JET: return string("Jet");
	case PT_BJET: return string("B-Jet");
	case PT_CJET: return string("C-Jet");
	case PT_CELL: return string("Cell");
	case PT_CLUSTER: return string("Cluster");
	case PT_MUON: return string("Muon");
	case PT_ELECTRON: return string("Electron");
	case PT_PHOTON: return string("Photon");
	case PT_TAU: return string("Tau");
	default: return string("Unknown");
	}
}

string Particle::getStateName() const {
	switch(state) {
	case PS_BEAM: return string("Beam");
	case PS_FINAL: return string("Final");
	case PS_DECAYED: return string("Decayed");
	case PS_HISTORY: return string("Historical");
	default: return string("Unknown");
	}
}
