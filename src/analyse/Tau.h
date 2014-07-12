#ifndef _AcerDet_Analyse_Tau_
#define _AcerDet_Analyse_Tau_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Tau {
		private:
			double ETJET;
			double PTTAU;
			double ETATAU;
			double RJTAU;
			double PTFRAC;

			int KEYHID;
			bool KEYTAU;

			int IEVENT;
			
			core::Histogram histo_jets;
			core::Histogram histo_taus;
		public:
			Tau( const Configuration& );
			~Tau();

			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
		};
		
	}
}

#endif
