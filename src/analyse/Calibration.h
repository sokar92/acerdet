#ifndef _AcerDet_Analyse_Calibration_
#define _AcerDet_Analyse_Calibration_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		//! Title
		/*!
		 * Detailed description
		 */
		class Calibration {
		private:
			Real64_t RCONE; /*!< detailed description  */
			
			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYCAL; /*!< detailed description  */
			
			Int32_t IEVENT; /*!< detailed description  */
			Int32_t IDENT; /*!< detailed description  */
		public:
			Calibration( const Configuration& );
			~Calibration();
			
			//! Print information about Calibration algorithm class to standard output
			void printInfo() const;
			
			//! Analyse input record
			/*!
			 * \param input <input data desc>
			 * \param output <output data desc>
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );
		};

	}
}

#endif
