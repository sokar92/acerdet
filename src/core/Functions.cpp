#include "Functions.h"
using namespace AcerDet::core;

Bool_t AcerDet::core::isHardProcess(const vector<Particle>& parts, int i) {
  //    particle is hard if is final
  //    one can use it if the below check is done recursively, i.e
  //    requiring any accessor to be H, W, or Z boson
  //	if (parts[i].status != PS_FINAL)
  //		return false;

	// if mother exists should be H, W, or Z boson {23, 24, 25}, top-quarks or tau lepton
	if (parts[i].hasMother()) {
		Int32_t mother_index = parts[i].mother;
		ParticleType mother_type = parts[mother_index].type;
		return     mother_type == PT_BOSON_Z
			|| mother_type == PT_BOSON_W
			|| mother_type == PT_BOSON_H
			|| mother_type == PT_TopQuark
			|| mother_type == PT_TAU;
	}

	return true;
}
