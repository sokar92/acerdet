#ifndef _AcerDet_Core_Particle_
#define _AcerDet_Core_Particle_

#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <ostream>
using namespace std;

#include "Typedefs.h"
#include "StlDefs.h"
#include "Consts.h"
#include "Vector4.h"

namespace AcerDet {
	namespace core {
		
		/**
		  A current state of particle
		*/
		enum ParticleStatus {
			PS_NULL = 0,      /*!< Undefined Particle state. Should not occur in valid program. */
			PS_BEAM,          /*!< TODO: detailed description  */
			PS_FINAL,         /*!< TODO: detailed description  */
			PS_DECAYED,       /*!< TODO: detailed description  */
			PS_HISTORY,       /*!< TODO: detailed description  */
			PS_CASCADE_QUARK, /*!< b or c quark outside of hard process. */
			PS_HP_QUARK       /*!< b or c quark from hard process. */
		};
		
		/**
		  Type of particle
		*/
		enum ParticleType {
			PT_UNKNOWN = 0,            /*!< Unknown Particle type. If this type occurs in some particle check your algorithm! */
			
			PT_JET,                    /*!< Parton representing Jet of grouped particles. */
			PT_BJET,         /* 5 */   /*!< Parton representing a Jet of b particles. */
			PT_CJET,         /* 4 */   /*!< Parton representing a Jet of c particles. */
			
			PT_CELL,                   /*!< Grouping parton representing a single Cell in detector. */
			PT_CLUSTER,                /*!< Grouping parton representing a Cluster aka. collection of consecutive Cells in detector. */
			
			PT_MUON,         /* 13 */  /*!< A muon particle. */
			PT_ELECTRON,     /* 11 */  /*!< An electron or positron (e- / e+). */
			PT_PHOTON,       /* 22 */  /*!< Nonmassive gamma particle. */
			PT_TAU,          /* 15 */  /*!< A tau particle. */
			
			PT_NEUTRINO_ELE, /* 12 */  /*!< An electron neutrino. */
			PT_NEUTRINO_MUO, /* 14 */  /*!< A muon neutrino. */
			PT_NEUTRINO_TAU, /* 16 */  /*!< A tau neutrino. */
			
			PT_BOSON_Z,      /* 23 */  /*!< Z0 boson. */
			PT_BOSON_W,      /* 24 */  /*!< W+ boson. */
			PT_BOSON_H       /* 25 */  /*!< Higgs boson. */
		};
		
		/**
		 * Common class representing partons.
		 * Stores it's type, status, four-vectors representing momentum and production vertex (if exists), mother and daughters ids.
		 */
		class Particle {
		public:
			ParticleStatus status; /*!< Enumerated particle status. */
			Int32_t statusID;      /*!< Original not converted particle status (comes from generator). */

			ParticleType type;     /*!< Enumerated particle type. (converted from particle PDGid). */
			Int32_t pdg_id;        /*!< Original not converted particle PDGid. (global convention). */
			
			Vector4f momentum;     /*!< A four-vector representing total particle momentum. */
			Vector4f production;   /*!< A four-vector representing particle production vertex coordinates (if such vertex exists). */

			Int32_t barcode;       /*!< Particle barcode. A unique id representing particle. */
			Int32_t mother;        /*!< If particle has mother, mother's barcode. Otherwise -1.*/
			pair<Int32_t,Int32_t> daughters; /*!< If particle has daughters, a range of daughter barcodes. Otherwise <-1,-1>. */

			/**
			 * Default particle constructor.
			 * Initialize particle data.
			 */
			Particle();

			/**
			 * Check if given particle has mother.
			 */
			Bool_t hasMother() const;

			/**
			 * Check if given particle has any daughters.
			 */
			Bool_t hasDaughter() const;

			/**
			 * Number of particle daughters.
			 */
			Int32_t daughtersCount() const;

			/**
			 * Check if given particle is beam.
			 */
			Bool_t isBeam() const;
			
			/**
			 * Check if given particle is a nautrino.
			 */
			Bool_t isNeutrino() const;

			/**
			* Print out to std::ostream basic information about given particle.
			* Prints out barcode, type, status, momentum and production vertex.
			*/			
			friend ostream& operator << (ostream& str, const Particle& p) {
				str << "Particle:" << endl;
				str << "\tBarcode: " << p.barcode << endl;
				str << "\tType: " << p.getTypeName() << "(" << p.pdg_id << ")" << endl;
				str << "\tStatus: " << p.getStatusName() << "(" << p.statusID << ")" << endl;
				str << "\tMomentum (px,py,pz,e) = " << p.momentum << endl;
				str << "\tProduction (x,y,z,t) = " << p.production << endl;
				str << endl;
				return str;
			}
			
			/**
			 * Length of momentum four-vector projected on XY plane.
			 * Returns sqrt(x*x + y*y), where [x,y,z,e] is a particle momentum.
			 */
			Real64_t pT() const;
			
			/**
			 * Particle momentum X-coordinate, where [x,y,z,e] is a particle momentum.
			 */
			inline Real64_t pX() const { return momentum.p.x; }
			
			/**
			 * Particle momentum Y-coordinate, where [x,y,z,e] is a particle momentum.
			 */
			inline Real64_t pY() const { return momentum.p.y; }
			
			/**
			 * Particle momentum Z-coordinate, where [x,y,z,e] is a particle momentum.
			 */
			inline Real64_t pZ() const { return momentum.p.z; }
			
			/**
			 * Particle momentum E-coordinate, where [x,y,z,e] is a particle momentum.
			 */
			inline Real64_t e() const { return momentum.e; }
			
			/**
			 * Particle momentum phi angle after projection on XY plane.
			 */
			Real64_t getPhi() const;
			
			/**
			 * Particle momentum eta angle (as a logarithm).
			 */
			Real64_t getEta() const;
			
			/**
			 * Particle momentum theta angle.
			 */
			Real64_t getTheta() const;
			
			/**
			 * Particle momentum phi angle folding.
			 * Returns 0.5 / getPhi()
			 */
			Real64_t foldPhi() const;

		private:
			string getTypeName() const;
			
			string getStatusName() const;
		};
	}
}

#endif
