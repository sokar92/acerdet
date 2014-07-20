#include "Particle.h"
#include "Functions.h"
using namespace AcerDet::core;

#define ParticleNull (-1)

Particle::Particle() : 
	state(PS_NULL), 
	stateID(-1),
	type(PT_UNKNOWN),
	typeID(0),
	momentum(),
	production(),
//	phi(0), theta(0), 
	id(-1), mother(-1), 
	daughters(make_pair(-1,-1)) 
{}

/*
 * Properties
 */

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

/*
 * Useful
 */
Real64_t Particle::pT() const {
	return sqrt(momentum.p.x * momentum.p.x + momentum.p.y * momentum.p.y);
}

/*
 * Angles
 */
 
Real64_t Particle::getPhi() const {
	return angle(momentum.p.x, momentum.p.y);
}

Real64_t Particle::getEta() const {
	return sign(log( (momentum.p.length() + abs(momentum.p.z)) / pT()), momentum.p.z);
}

Real64_t Particle::getTheta() const {
	if (abs(momentum.p.z / momentum.e) < 1.0) 
		return acos(momentum.p.z / momentum.p.length());
	return momentum.p.z > 0 ? 0 : PI;
}

/*
 * Private functions
 */

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
