#include "Smearing.h"
#include <cmath>
using namespace AcerDet::core;

std::default_random_engine Smearing::generator;
std::normal_distribution<Real64_t> Smearing::distribution(0.0, 1.0);

Real64_t Smearing::forHadron(
	Real64_t ene,
	Real64_t eta,
	Real64_t caloth) 
{
	Real64_t aa, sigma;
	do {
		aa = distribution(generator);
		if (abs(eta) < caloth) sigma = aa * 0.5 / sqrt(ene);
		else sigma = aa / sqrt(ene);
	} while (1.0 + sigma <= 0.0);
	return sigma;
}

Real64_t Smearing::forElectron(Real64_t ene) {
	Real64_t aa, sigma;
	do {
		aa = distribution(generator);
		sigma = aa * 0.12 / sqrt(ene);
	} while (1.0 + sigma <= 0.0);
	return sigma;
}

Real64_t Smearing::forMuon(Real64_t pt) {
	Real64_t aa, sigma;
	do {
		aa = distribution(generator);
		sigma = 0.0005 * pt * aa;
	} while (1.0 + sigma <= 0.0);
	return sigma;
}

Real64_t Smearing::forPhoton(Real64_t ene) {
	Real64_t aa, sigma;
	do {
		aa = distribution(generator);
		sigma = aa * 0.1 / sqrt(ene); 
	} while(1.0 + sigma <= 0.0);
	return sigma;
}
