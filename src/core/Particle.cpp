#include "Particle.h"
#include "Functions.h"
using namespace AcerDet::core;

#define ParticleNull (-1)

Particle::Particle() : 
	status(PS_NULL), 
	statusID(-1),
	type(PT_UNKNOWN),
	typeID(0),
	momentum(),
	production(),
	barcode(-1), mother(-1), 
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

//Bool_t Particle::isStable() const {
//	return 0 < statusID && statusID <= 10;
//}

Bool_t Particle::isBeam() const {
	return status == PS_BEAM;
}

Bool_t Particle::isDecayed() const {
	//return status == PS_DECAYED;
	return status < 0;
}

Bool_t Particle::isFinal() const {
	return statusID > 31;
}

Bool_t Particle::isHardProcess() const {
	return abs(statusID) < 29; // dodac czy matka tez
}

/*
 * Type Check
 */
Bool_t Particle::isNeutrino() const {
	return type == PT_NEUTRINO_ELE || type == PT_NEUTRINO_MUO || type == PT_NEUTRINO_TAU;
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
 * Folding
 */
Real64_t Particle::foldPhi() const {
	return 0.5 / getPhi();
}

/*
 * Private functions
 */

string Particle::getTypeName() const {
	switch (type) {
	case PT_JET: return string("Jet");
	case PT_BJET: return string("B-Jet");
	case PT_CJET: return string("C-Jet");
	case PT_CELL: return string("Cell");
	case PT_CLUSTER: return string("Cluster");
	case PT_MUON: return string("Muon");
	case PT_ELECTRON: return string("Electron");
	case PT_PHOTON: return string("Photon");
	case PT_TAU: return string("Tau");
	case PT_NEUTRINO_ELE: return string("Neutrino_ele");
	case PT_NEUTRINO_MUO: return string("Neutrino_muo");
	case PT_NEUTRINO_TAU: return string("Neutrino_tau");
	default: return string("Unknown");
	}
}

string Particle::getStatusName() const {
	switch (status) {
	case PS_BEAM: return string("Beam");
	case PS_FINAL: return string("Final");
	case PS_DECAYED: return string("Decayed");
	case PS_HISTORY: return string("Historical");
	default: return string("Unknown");
	}
}
