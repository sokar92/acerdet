#ifndef _AcerDet_Analyse_Calibration_
#define _AcerDet_Analyse_Calibration_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		//! Analyse calibration
		/*!
		 * Calibration
		 */
		class Calibration {
		private:
			Real64_t RCONE; /*!< see Configuration */
			
			Int32_t KEYHID; /*!< see Configuration */
			Bool_t KEYCAL;  /*!< see Configuration */
			
			Int32_t IEVENT; /*!< see Configuration */
			Int32_t IDENT;  /*!< see Configuration */
			
		public:
			
			/**
			 * Constructor.
			 * Creates new instance of Calibration.
			 * \param conf AcerDET initial Configuration
			 */
			Calibration( const Configuration& );
			
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
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );
		};

	}
}

#endif
