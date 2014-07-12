#ifndef _AcerDet_Analyse_Mis_
#define _AcerDet_Analyse_Mis_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Mis {
		private:
			double PTMUMIN;
			double ETAMAX;
			double ETCELL;
			double CALOTH;
			
			int KEYHID;
			bool KEYSME;
			int KFINVS;
			
			int IEVENT;
			
			core::Histogram histo_reconstructed_pT; //11
			core::Histogram histo_reconstructed_pT_cells; //12
			core::Histogram histo_pTmiss; //13
			core::Histogram histo_pTnu;	//21
		
		public:
			Mis( const Configuration& );
			~Mis();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
		};

	}
}

#endif
