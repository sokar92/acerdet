#ifndef _AcerDet_Analyse_Calibration_
#define _AcerDet_Analyse_Calibration_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
#include "./IAnalysePhase.h"

namespace AcerDet {
	namespace analyse {

		//! Analyse calibration
		/*!
		 * Calibration
		 */
		class Calibration : public IAnalysePhase {
		private:
			Real64_t RCONE; /*!< see Configuration */
			
			Int32_t KEYHID; /*!< see Configuration */
			Bool_t KEYCAL;  /*!< see Configuration */
			
			Int32_t IEVENT; /*!< see Configuration */
			
			IHistogramManager *histoManager; /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;          /*!< Indicates whether histograms are already registered and might be used now. */
			
		public:
			
			/**
			 * Constructor.
			 * Creates new instance of Calibration.
			 * \param conf AcerDET initial Configuration
			 */
			Calibration( const Configuration&, IHistogramManager* );
			
			/**
			 * Destructor.
			 */
			~Calibration();
			
			/**
			 * Print information about Calibration algorithm class to standard output.
			 */
			void printInfo() const;
			
			/*!
			 * \param input record.
			 * \param output record.
			 * \param weight event weight.
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output,
				Real64_t weight );
			
			/**
			 * Print Calibration algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
