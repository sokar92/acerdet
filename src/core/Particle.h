#ifndef _AcerDet_Core_Particle_
#define _AcerDet_Core_Particle_

#pragma once

#include <cmath>
#include <cstdio>
#include <string>
using namespace std;

#include "Typedefs.h"
#include "StlDefs.h"
#include "Consts.h"

namespace AcerDet {
	namespace core {
		
		/*
		  A state of particle
		*/
		enum ParticleState {
			PS_NULL = 0,
			PS_BEAM,
			PS_FINAL,
			PS_DECAYED,
			PS_HISTORY
		};
		
		/*
		  Type of particle - self descriptive
		*/
		enum ParticleType {
			PT_UNKNOWN = 0,
			PT_JET,
			PT_BJET,
			PT_CJET,
			PT_CELL,
			PT_CLUSTER,
			PT_MUON,
			PT_ELECTRON,
			PT_PHOTON,
			PT_TAU
		};
		
		/*
		 * Common class representing partons
		 */
		class Particle {
		public:
			ParticleState state;
			Int32_t stateID;
			
			ParticleType type;
			Int32_t typeID;
			
			Real64_t px, py, pz, e; //momentum
			Real64_t phi, theta; //angles
			Real64_t prod_x, prod_y, prod_z, prod_time;	//production
			
			Int32_t id,mother;
			pair<Int32_t,Int32_t> daughters;
			
			/*
			 * Default ctor
			 * Initialize particle data
			 */
			Particle();
			
			/*
			 * Check if given particle has mother
			 */
			Bool_t hasMother() const;
			
			/*
			 * Check if given particle has any daughters
			 */
			Bool_t hasDaughter() const;
			
			/*
			 * Number of particle daughters
			 */
			Int32_t daughtersCount() const;
			
			/*
			 * Check if particle is stable
			 */
			Bool_t isStable() const;
			
			/*
			 * Check if particle is beam
			 */
			Bool_t isBeam() const;
			
			/*
			 * Check if particle has decayed status
			 */
			Bool_t isDecayed() const;
			
			/*
			 * Print basic informations about particle
			 */
			void print() const;
			
		private:
			string getTypeName() const;
			
			string getStateName() const;
		};

		/*
		  Single particle class with basic methods such as boost, rotation,
		  and angle calculation.
		*/

/*
		public:
			T getPx() const { return _px; }
			T getPy() const { return _py; }
			T getPz() const { return _pz; }
			T getE () const { return _e; }
			
			T getPhi() const { return _phi; }
			T getTheta() const { return _theta; }
			
			Int32_t getID() const { return _id; }
			pair<Int32_t,Int32_t> getDaughters() const { return _idDaughter; }
			ParticleState getState() const { return _state; }
			
			Bool_t	isBeam() const { return _state == PS_BEAM; }

			void setPx(T px) { _px = px; }
			void setPy(T py) { _py = py; }
			void setPz(T pz) { _pz = pz; }
			void setE (T e ) { _e  = e;  }
			
			void setPhi(T phi) { _phi = phi; }
			void setTheta(T theta) { _theta = theta ; }
			
			void setID(Int32_t id) { _id = id; }
			void setMother(Int32_t mother) { _idMother = mother; }
			void setDaughters(Int32_t first, Int32_t last) { _idDaughter = make_pair(first,last); }
			void setState(ParticleState ps) { _state = ps; }

			// Invariant mass. If mass is negative then -sqrt(-mass) is returned
			T recalculated_mass()  const
			{
				T p2 = _px*_px + _py*_py + _pz*_pz;
				T e2 = _e*_e;
				T m2 = e2 - p2;

				if ( m2 < 0.0 ) return -sqrt( -m2 );
				return sqrt( m2 );
			}

			void print()
			{
				if (_id)	printf("%4d: %15.8e %15.8e %15.8e %15.8e  | %15.8e\n", getID(), getPx(), getPy(), getPz(), getE(), recalculated_mass());
				else 		printf(" SUM: %15.8e %15.8e %15.8e %15.8e  | %15.8e\n", getPx(), getPy(), getPz(), getE(), recalculated_mass());
			}
		public:
			double getAnglePhi();
			double getAngleTheta();
			void rotateXZ(double theta);
			void rotateXY(double theta);
			void boostAlongZ(double pz, double e);
			void boostToRestFrame(ParticleT<T> &p);
			void boostFromRestFrame(ParticleT<T> &p);
		};

		template<typename T>
		inline double ParticleT<T>::getAnglePhi()
		{
			// conventions as in ANGFI(X,Y) of tauola.f and PHOAN1 of photos.f
			// but now 'phi' in name define that it is rotation in px py

			double buf = 0.;

			if ( fabs(getPy()) < fabs(getPx()) ) {
				buf = atan( fabs(getPy()/getPx()) );
				if (getPx() < 0.) buf = M_PI - buf;
			} else buf = acos( getPx()/sqrt(getPx()*getPx() + getPy()*getPy()) );

			if ( getPy() < 0. ) buf = 2.*M_PI - buf;
			return buf;
		}

		template<typename T>
		inline double ParticleT<T>::getAngleTheta()
		{
			// conventions as in ANGXY(X,Y) of tauola.f or PHOAN2 of photos.f
			// but now 'theta' in name define that it is rotation in pz px
			// note that first argument of PHOAN2 was usually z

			double buf = 0.;

			if ( fabs(getPx()) < fabs(getPz()) ) {
				buf = atan( fabs(getPx()/getPz()) );
				if ( getPz() < 0. ) buf = M_PI - buf;
			}
			else buf = acos( getPz()/sqrt(getPz()*getPz() + getPx()*getPx()) );
			return buf;
		}

		template<typename T>
		inline void ParticleT<T>::rotateXZ(double theta)
		{
			// as PHORO2 of photos.f
			T pX = getPx();
			T pZ = getPz();
			setPx( cos(theta)*pX + sin(theta)*pZ);
			setPz(-sin(theta)*pX + cos(theta)*pZ);
		}

		template<typename T>
		inline void ParticleT<T>::rotateXY(double theta)
		{
			// as PHORO3 of photos.f
			T pX = getPx();
			T pY = getPy();
			setPx( cos(theta)*pX - sin(theta)*pY);
			setPy( sin(theta)*pX + cos(theta)*pY);
		}

		template<typename T>
		inline void ParticleT<T>::boostAlongZ(double p_pz, double p_e)
		{
			// as PHOBO3 of photos.f
			double m = sqrt(p_e*p_e-p_pz*p_pz);
			T buf_pz = getPz();
			T buf_e = getE();

			setPz((p_e *buf_pz + p_pz*buf_e)/m);
			setE ((p_pz*buf_pz + p_e *buf_e)/m);
		}

		template<typename T>
		inline void ParticleT<T>::boostToRestFrame(ParticleT<T> &p)
		{
			T p_len = sqrt(p.getPx()*p.getPx() + p.getPy()*p.getPy() + p.getPz()*p.getPz());
			double phi   = p.getAnglePhi();
			p.rotateXY(-phi );
			double theta = p.getAngleTheta();
			p.rotateXY( phi );

			//Now rotate coordinates to get boost in Z direction.
			rotateXY(-phi  );
			rotateXZ(-theta);
			boostAlongZ(-1*p_len, p.getE());
			rotateXZ( theta);
			rotateXY( phi  );
		}

		template<typename T>
		inline void ParticleT<T>::boostFromRestFrame(ParticleT<T> &p)
		{
			T p_len = sqrt(p.getPx()*p.getPx() + p.getPy()*p.getPy() + p.getPz()*p.getPz());
			double phi   = p.getAnglePhi();
			p.rotateXY(-phi );
			double theta = p.getAngleTheta();
			p.rotateXY( phi );

			//Now rotate coordinates to get boost in Z direction.
			rotateXY(-phi  );
			rotateXZ(-theta);
			boostAlongZ(p_len, p.getE());
			rotateXZ( theta);
			rotateXY( phi  );
		}
*/
	}
}

#endif
