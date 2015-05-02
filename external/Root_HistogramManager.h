#ifndef _AcerDet_External_Root_HistogramManager_
#define _AcerDet_External_Root_HistogramManager_

#include "../src/core/IHistogramManager.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace external {
		
		//! ROOT specialization for AcerDET core::IHistogramManager
		class Root_HistogramManager 
			: public IHistogramManager {
		public:
		
			/**
			 * Initialize histogram manager.
			 */
			void init();
			
			/**
			 * Create new instance of histogram and store it's reference in histogram manager.
			 * \param id new histogram id - should be unique in manager.
			 * \param title new histogram name.
			 * \param blocks new histogram resolution.
			 * \param minVal minimal value which can be stored in histogram.
			 * \param maxVal maximal value which can be stored in histogram.
			 */
			void registerHistogram(
				Int32_t id,
				const string& title,
				Int32_t blocks,
				Real64_t minVal,
				Real64_t maxVal );
			
			/**
			 * Add new weighted record to histogram.
			 * \param id histogram id to store record in.
			 * \param value value to store.
			 * \param weight weight of stored record.
			 */
			void insert( Int32_t id, Real64_t value, Real64_t weigth );
			
			/**
			 * Dump all histogram data to file system.
			 */
			void storeHistograms();
		};

	}
}

#endif
