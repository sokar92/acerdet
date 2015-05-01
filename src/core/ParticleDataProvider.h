#ifndef _AcerDet_Core_ParticleDataProvider_
#define _AcerDet_Core_ParticleDataProvider_

#include <string>
#include <map>
using namespace std;

#include "Typedefs.h"

#pragma once

namespace AcerDet {
	namespace core {

		//! A mapping between particle PDGid and it's data.
		/**
		 * A mapping between particle PDGid and it's data.
		 * Record data contains information about particle charge and name.
		 */
		class ParticleDataProvider {
		private:

			/**
			 * Internal structure representing particle data record.
			 */
			struct Data {
				Int32_t chargeType; /*!< Particle charging type. Determines whether particle is charged or not. */
				Real64_t charge;    /*!< Particle charge (if charged). */
				string name;        /*!< A name of particle - mapping from symbol to names. */
			};
			
			map<Int32_t, Data> parts; /*!< Mapping between particle PDGid and it's data record. */

		public:
			/**
			 * Default constructor.
			 * Creates a new instance for ParticleDataProvider.
			 */
			ParticleDataProvider();
			
			/**
			 * Checks if record base contains info about particle described by PDGid.
			 * \param pdgId particle PDGid.
			 * \return true if record base contains info about requested PDGid.
			 */
			Bool_t containsInfo(Int32_t pdgId) const;
			
			/**
			 * Insert new data record to records base.
			 * \param pdgId particle PDGid.
			 * \param chgType particle charge type.
			 * \param chg particle charge (if charged).
			 * \param name particle name.
			 */
			void insertInfo(Int32_t pdgId, Int32_t chgType, Real64_t chg, const string& name);
			
			/**
			 * Maps particle PDGid to it's charge type.
			 * \param pdgId particle PDGid
			 * \return charge type of particle.
			 */
			Int32_t getChargeType(Int32_t pdgId) const;
			
			/**
			 * Maps particle PDGid to it's charge (if charged).
			 * \param pdgId particle PDGid
			 * \return charge of particle.
			 */
			Real64_t getCharge(Int32_t pdgId) const;
			
			/**
			 * Maps particle PDGid to it's name.
			 * \param pdgId particle PDGid
			 * \return name of particle.
			 */
			const string& getName(Int32_t pdgId) const;
		};
	}
}

#endif
