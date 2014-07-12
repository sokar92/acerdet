#ifndef _AcerDet_Library_Included_
#define _AcerDet_Library_Included_

// ---------------------------------------
// -- Include AcerDET 2.0 Configuration -- 
// ---------------------------------------
#include "conf/Configuration.h"

// ------------------------------------------
// --  Include AcerDET 2.0 core components --
// ------------------------------------------
#include "core/Histogram.h"

// -------------------------------------------
// --  Include AcerDET 2.0 input components --
// -------------------------------------------
#include "io/InputRecord.h"
#include "io/OutputRecord.h"
#include "io/InputReader.h"

// ---------------------------------------------
// --  Include AcerDET 2.0 analyse components --
// ---------------------------------------------
#include "analyse/BJet.h"
#include "analyse/Calibration.h"
#include "analyse/Cell.h"
#include "analyse/CJet.h"
#include "analyse/Cluster.h"
#include "analyse/Electron.h"
#include "analyse/Jet.h"
#include "analyse/Mis.h"
#include "analyse/Muon.h"
#include "analyse/Photon.h"
#include "analyse/Tau.h"

namespace AcerDet {

	class AcerDET {
	private:
		analyse::BJet			analyse_BJet;
		analyse::Calibration	analyse_Calibration;
		analyse::Cell			analyse_Cell;
		analyse::CJet			analyse_CJet;
		analyse::Cluster		analyse_Cluster;
		analyse::Electron		analyse_Electron;
		analyse::Jet			analyse_Jet;
		analyse::Mis			analyse_Mis;
		analyse::Muon			analyse_Muon;
		analyse::Photon			analyse_Photon;
		analyse::Tau			analyse_Tau;
	public:
		AcerDET( const conf::Configuration& );
		~AcerDET();
		
		void printInfo() const;
		
		void analyseRecord( const io::InputRecord&, io::OutputRecord& );
	};

}

#endif
