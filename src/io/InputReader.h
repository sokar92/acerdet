#ifndef _AcerDet_IO_InputReader_
#define _AcerDet_IO_InputReader_

// HepMC event record
#include <HepMC/GenEvent.h>
using namespace HepMC;

#include "InputRecord.h"

namespace AcerDet {
	namespace io {
		
		/*
		 */
		class InputReader {
		public:
			static InputRecord computeEvent(const GenEvent&);
		};

	}
}

#endif
