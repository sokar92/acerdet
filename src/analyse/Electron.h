#ifndef _AcerDet_Analyse_Electron_
#define _AcerDet_Analyse_Electron_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Electron {
		private:
			double ETCLU;
			double RCONE;
			double PTLMIN;
			double ETAMAX;
			double RJE;
			double RISOLJ;
			double RDEP;
			double EDMAX;

			int KEYHID;
			bool KEYSME;

			int IEVENT;
			
			core::Histogram histo_isol;
			core::Histogram histo_hard;
			core::Histogram histo_sum;
		public:
			Electron( const Configuration& );
			~Electron();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
		};

	}
}

#endif
