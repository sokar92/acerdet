#ifndef _AcerDet_Analyse_Cell_
#define _AcerDet_Analyse_Cell_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../core/ParticleDataProvider.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace analyse {

		/**
		 * Algorithm class for finding cells
		 */
		class Cell {
		private:
			Real64_t ETACEL;
			Real64_t PTMIN;
			Real64_t ETTHR;
			Real64_t CALOTH;
			Real64_t DBETA;
			Real64_t DBPHI;
			
			Int32_t KEYHID;
			Bool_t  KEYFLD;
			Int32_t KFINVS;
			
			Int32_t IEVENT;
			
			IHistogramManager *histoManager;
			Bool_t histoRegistered;

			const ParticleDataProvider& partProvider;
			
		public:
			Cell( 
				const Configuration&,
				IHistogramManager*,
				const ParticleDataProvider& );

			~Cell();
			
			/**
			 * Print basic informations about class
			 * and configuration to stdout.
			 */
			void printInfo() const;
			
			/**
			 * Analyse input record.
			 */
			void analyseRecord( const InputRecord&, OutputRecord& ); 
			
			/**
			 * Prints statistics to stdout.
			 */
			void printResults() const;
		};

	}
}

#endif
