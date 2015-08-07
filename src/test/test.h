#ifndef _AcerDet_test_Test_
#define _AcerDet_test_Test_

#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <ostream>
using namespace std;

#include "../core/IHistogramManager.h"
#include "../core/ObjectData.h"
#include "../core/Particle.h"
#include "../core/JetData.h"
using namespace AcerDet::core;

#include "../core/Typedefs.h"
#include "../core/StlDefs.h"
#include "../core/Consts.h"
#include "../core/Vector4.h"

namespace AcerDet {
	namespace test {
		
		
		//! Helper class for creating test histograms. 
		/**
		 * Helper class for creating test histograms. 
		 * Stores particle data needed to create specific types of hostograms.
		 */
		class Test {
		public:		    
			
			vector<JetData> jet;			/*!< Stores all jets from a given event  */
			vector<Particle> HARDphoton;	/*!< Stores hard photons from a given event  */
			vector<Particle> HARDelectron;	/*!< Stores hard electrons from a given event  */
			vector<Particle> HARDmuon;		/*!< Stores hard muons from a given event  */
			vector<Particle> bquark;		/*!< Stores b-type quarks from a given event  */
			vector<Particle> cquark;		/*!< Stores c-type from a given event  */
			vector<Particle> uquark;		/*!< Stores u-type from a given event  */
			vector<ObjectData> ISOphoton;	/*!< Stores isolated photons from a given event  */
			vector<ObjectData> ISOelectron;	/*!< Stores isolated electrons from a given event  */
			vector<ObjectData> ISOmuon;		/*!< Stores isolated muons from a given event  */
			vector<JetData> bjet;			/*!< Stores b-type jets from a given event  */
			vector<JetData> cjet;			/*!< Stores c-type jets from a given event  */
			vector<JetData> ujet;			/*!< Stores jets suspected to be u-type from a given event  */
			vector<Particle> particle1;		/*!< Used to differentiate hard particles into particles and antiparticles. Stores particles */
			vector<Particle> particle2;		/*!< Used to differentiate hard particles into particles and antiparticles. Stores antiparticles */

			/**
			 * Default test constructor.
			 * Initialize test data.
			 */
			Test();
			
			/**
			 * Destructor.
			 */
			~Test();
			
			/**
			 * Returns square of an invariant mass. Uses values stored in a Vector4f momentum. 
			 */
			inline Real64_t mass2() const { return ( pow(momentum.e , 2.0) - pow(momentum.p.x , 2.0) - pow( momentum.p.y , 2.0) - pow(momentum.p.z , 2.0)); }
			
			/**
			 * Returns an invariant mass.
			 */
			inline Real64_t mass() const { return sqrt( mass2() ); }
			
			/**
			 * Returns a delta r for pair quark-jet.
			 */
			Real64_t deltaR( Particle part, JetData jet );
			
			/**
			 * Routine for filling histograms specific for Z->2e process.
			 */
			void Fill2e(Int32_t idhist, Real64_t weight, IHistogramManager* histoManager);
			
			/**
			 * Routine for filling histograms specific for H->2gamma process.
			 */
			void Fill2gamma(Int32_t idhist, Real64_t weight, IHistogramManager* histoManager);
			
			/**
			 * Routine for filling histograms specific for H->2Z->4X process. Works on hard particles ( type Particle ).
			 */
			void Fill4part(vector <Particle> part, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager);
			
			/**
			 * Routine for filling histograms specific for H->2Z->4X process. Works on isolated particles ( type ObjectData ).
			 */
			void Fill4object(vector <ObjectData> newObject, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager);
			
			/**
			 * Routine for filling histograms specific for H->2q process.
			 */
			void Filljets( vector<Particle> quark, vector<JetData> jet, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager);
			
			/**
			 * Converts momentum to [E,px,py,pz] form and add it to Vector4f momentum.  
			 */
			void AddMomenta(Real64_t pT, Real64_t eta, Real64_t phi);
			
			/**
			 * Sums all momenta from given ObjectData type container.  
			 */
			void allobjectmomenta(vector<ObjectData> newObject);
			
			/**
			 * Sums all momenta from given Particle type container.  
			 */
			void allparticlemomenta(vector<Particle> part);
			
			/**
			 * Sums 2 chosen momenta ( of it_1 and it_2 position in the container ) from given JetData type container.  
			 */
			void selectjetmomenta( vector<JetData> jet, Int32_t it_1, Int32_t it_2 );
			
			/**
			 * Fills particle1 and particle2 containers. 
			 * Particles of pdg_id > 0 go to particle1 container,  particles of pdg_id < 0 to particle2.
			 */
			void Dividepart( vector<Particle> part );
			
			/**
			 * Fills object1 and object2 containers. 
			 * Particles of pdg_id > 0 go to object1 container,  particles of pdg_id < 0 to object2.
			 */
			void Divideobject( vector<ObjectData> newObject );
			
			/**
			 * Fills bjet, cjet and ujet containers. 
			 * Jets of type B_JET go to bjet container,  jets of type C_JET to cjet container and jets of type UNKNOWN to ujet container.
			 */
			void Dividejets();
			
			/**
			 * Used for H->2Z->4X processes. Fills storage1 and storage2 containers with masses of particle-antiparticle pairs. 
			 * Storage1 is associated with first particle in object1 container and storage2 with second particle in object1 container.
			 */
			void groupobjectmasses();   
			
			/**
			 * Used for H->2Z->4X processes. Fills storage1 and storage2 containers with masses of particle-antiparticle pairs. 
			 * Storage1 is associated with the first particle in particle1 container. 
			 * Storage2 is associated with the second particle in particle1 container. 
			 */
			void groupparticlemasses();
			
			/**
			 * Used for H->2q processes. Fills storage1 and storage2 containers with values of delta r of quark-jet pairs. 
			 * Storage1 is associated with the first quark in Particle type container. 
			 * Storage2 is associated with the second quark in Particle type container.
			 */
			void groupdeltaR( vector<JetData> jet );
			
			/**
			 *Clears content of storage1 and storage2 containers.
			 */
			void clearstorage() { storage1.clear(), storage2.clear(); }
			
			/**
			 * Returns a vector with 3 elements. The first element indicates 
			 * storage with the lowest overall element ( 0 for storage1 and 1 for storage2 ).
			 * The second and the third element indicate positions of the lowest element in each storage.
			 * If the lowest values for both storages have the same position then position of the lowest overall 
			 * element is returned where for the other storage the second lowest value is returned.
			 */
			vector<Int32_t> Findlowest(); 
			
			/**
			 * Returns a vector with 3 elements. The first element indicates 
			 * storage with the highest overall element ( 0 for storage1 and 1 for storage2 ).
			 * The second and the third element indicate positions of the highest element in each storage.
			 * If the highest values for both storages have the same position then position of the highest overall 
			 * element is returned where for the other storage the second highest value is returned.
			 */
			vector<Int32_t> Findhighest(); 
			
			/**
			 * Routine for filling Particle type containers.
			 */
			void PutHARD( Particle part );
			
			/**
			 * Fills given histogram with a given value under condition > 0 ( used for filling mass() with condition mass2() > 0 ).
			 */
			void PutHisto( Real64_t value, Real64_t condition, Int32_t id, Int32_t idhist, Real64_t weight, IHistogramManager* histoManager );
			
			
		private:
			Vector4f momentum;			/*!< Used to eveluate invariant mass. Stores 4-momentum. */
			vector<Real64_t> storage1;	/*!< Container */
			vector<Real64_t> storage2;	/*!< Container */
			vector<ObjectData> object1; /*!< Used to differentiate isolated particles into particles and antiparticles. Stores particles */
			vector<ObjectData> object2; /*!< Used to differentiate isolated particles into particles and antiparticles. Stores antiparticles */
			
			
		};
	}
}

#endif
