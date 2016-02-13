#ifndef _AcerDet_Analyse_Phase_I_
#define _AcerDet_Analyse_Phase_I_

#pragma once

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {
		
		class IAnalysePhase {
		public:
			virtual ~IAnalysePhase() {}
			
			virtual void printInfo() const = 0;
			
			virtual void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output,
				Real64_t weight ) = 0;
			
			virtual void printResults() const = 0;
		};
		
	}
}

#endif
