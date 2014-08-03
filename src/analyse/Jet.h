#ifndef _AcerDet_Analyse_Jet_
#define _AcerDet_Analyse_Jet_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Jet {
		private:
			Real64_t ETJET;
			Real64_t ETAJET;
			Real64_t RCONE;
			Real64_t PTMIN;
			Real64_t CALOTH;

			Int32_t KEYHID;
			Bool_t KEYSME;
			Bool_t KEYFLD;
			Int32_t KFINVS;

			Int32_t IEVENT;
			
			core::Histogram histo_bJets; //IDENT + 1
			core::Histogram histo_delta_phi; //IDENT + 11
			core::Histogram histo_delta_eta; //IDENT + 12
			core::Histogram histo_delta_barycenter; //IDENT + 13
			core::Histogram histo_delta_parton; //IDENT + 23
			core::Histogram histo_pT_bySum; //IDENT + 14
			core::Histogram histo_pT_byPart; //IDENT + 24
			
		public:
			Jet( const Configuration& );
			~Jet();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
