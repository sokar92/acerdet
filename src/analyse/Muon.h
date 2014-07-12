#ifndef _AcerDet_Analyse_Muon_
#define _AcerDet_Analyse_Muon_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Muon {
		private:
			double ETCLU;
			double RCONE;
			double PTMUMIN;
			double ETAMAX;
			double RISOLJ;
			double RDEP;
			double EDMAX;
			double PTMUMINT;

			int KEYHID;
			bool KEYSME;

			int IEVENT;
			
			core::Histogram histo_nonisol; //10
			core::Histogram histo_isol; //11
			core::Histogram histo_hard; //21
			core::Histogram histo_sum; //31
		public:
			Muon( const Configuration& );
			~Muon();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
		};

	}
}

#endif
