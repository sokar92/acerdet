#ifndef _AcerDet_Core_IHistogramManager_
#define _AcerDet_Core_IHistogramManager_

#include <string>
using namespace std;

#include "Typedefs.h"

#pragma once

namespace AcerDet {
	namespace core {
		
		//! Abstract interface describing contracts for any external implementation of histogram manager supported by AcerDET.
		/**
		 * Abstract interface describing contracts for any external implementation of histogram manager
		 * supported by AcerDET.
		 */
		class IHistogramManager {
		public:
		
			/**
			 * Virtual destructor.
			 */
			virtual ~IHistogramManager() {}
			
			/**
			 * Initialize histogram manager.
			 */
			virtual void init() = 0;
			
			/**
			 * Create new instance of histogram and store it's reference in histogram manager.
			 * \param id new histogram id - should be unique in manager.
			 * \param title new histogram name.
			 * \param blocks new histogram resolution.
			 * \param minVal minimal value which can be stored in histogram.
			 * \param maxVal maximal value which can be stored in histogram.
			 */
			virtual void registerHistogram(
				Int32_t id,
				const string& title,
				Int32_t blocks,
				Real64_t minVal,
				Real64_t maxVal ) = 0;
				
			/**
			 * Add new weighted record to histogram.
			 * \param id histogram id to store record in.
			 * \param value value to store.
			 * \param weight weight of stored record.
			 */
			virtual void insert( Int32_t id, Real64_t value, Real64_t weight = 1.0) = 0;
			
			/**
			 * Dump all histogram data to file system.
			 */
			virtual void storeHistograms() = 0;
		};

	}
}

#endif
