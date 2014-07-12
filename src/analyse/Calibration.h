#ifndef _AcerDet_Analyse_Calibration_
#define _AcerDet_Analyse_Calibration_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Calibration {
		private:
			double RCONE;
			
			int KEYHID;
			bool KEYCAL;
			
			int IEVENT;
			int IDENT;
		public:
			Calibration( const Configuration& );
			~Calibration();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
		};

	}
}

#endif
