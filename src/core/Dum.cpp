#include "Dum.h"
using namespace AcerDet::core;

Dum::Dum() {
	for (int i=0;i<5;++i) {
		K[i] = 0;
		P[i] = 0.0;
	}
}

Dum::Dum(const Dum& d) {
	for (int i=0;i<5;++i) {
		K[i] = d.K[i];
		P[i] = d.P[i];
	}
}

Dum& Dum::operator = (const Dum& d) {
	for (int i=0;i<5;++i) {
		K[i] = d.K[i];
		P[i] = d.P[i];
	}
	return *this;
}
