#ifndef _AcerDet_Core_ParticleDataProvider_
#define _AcerDet_Core_ParticleDataProvider_

#include <string>
#include <map>
using namespace std;

#include "Typedefs.h"

#pragma once

namespace AcerDet {
	namespace core {

		class ParticleDataProvider {
		private:
			struct Data {
				Int32_t chargeType; /*!< detailed description  */
				Real64_t charge; /*!< detailed description  */
				string name; /*!< detailed description  */
			};
			map<Int32_t, Data> parts; /*!< detailed description  */
		public:
			ParticleDataProvider();
			
			Bool_t containsInfo(Int32_t pdgId) const;
			void insertInfo(Int32_t pdgId, Int32_t chgType, Real64_t chg, const string& name);
			
			Int32_t getChargeType(Int32_t pdgId) const;
			Real64_t getCharge(Int32_t pdgId) const;
			const string& getName(Int32_t pdgId) const;
		};
	}
}

#endif
