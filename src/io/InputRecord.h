#ifndef _AcerDet_IO_InputRecord_
#define _AcerDet_IO_InputRecord_

#include "../core/Particle.h"
using namespace AcerDet::core;

#include <vector>
using namespace std;

namespace AcerDet {
	namespace io {

		class InputRecord {
/*
		public:
			Int32_t IdProc;

			// acmcproc.inc
			// INTEGER IDPRO
			// COMMON /ACMCPROC/ IDPRO
			Int32_t IdPro;

			// acdjets.inc
			// INTEGER MAXJET
			// PARAMETER(MAXJET=20)
			// INTEGER NJETS,KFJETS
			// REAL PXJETS,PYJETS,PZJETS,EEJETS
			// COMMON /ACDJETS/ NJETS,KFJETS(MAXJET),
			// + PXJETS(MAXJET),PYJETS(MAXJET),PZJETS(MAXJET),
			// + EEJETS(MAXJET)

			// const Int32_t MAX_JET = 20;
			Int32_t NJets, KfJets;
			// Real32_t PxJets, PyJest, PzJets, EeJets;

			// acdlept.inc
			// INTEGER MAXLEPT
			// PARAMETER(MAXLEPT=12)
			// INTEGER NLEPT,KFLEPT
			// REAL PXLEPT,PYLEPT,PZLEPT,EELEPT
			// COMMON /ACDLEPTT/ NLEPT,KFLEPT(MAXLEPT),
			// +                 PXLEPT(MAXLEPT),PYLEPT(MAXLEPT),PZLEPT(MAXLEPT),
			// +                 EELEPT(MAXLEPT)
			
			// const Int32_t MAX_LEPT = 12;
			Int32_t NLept, KfLept;
			// Real32_t PxLept, PyLept, PzLept, EeLept;

			// acdpart.inc
			// INTEGER MAXPART
			// PARAMETER(MAXPART=40)
			// INTEGER NPART,KSPART,KFPART
			// REAL PXPART,PYPART,PZPART,EEPART
			// COMMON /ACDPART/ NPART,KSPART(MAXPART),KFPART(MAXPART),
			// +                 PXPART(MAXPART),PYPART(MAXPART),
			// +                 PZPART(MAXPART),EEPART(MAXPART)
			
			// const Int32_t MAX_PART = 40;
			Int32_t NPart, KsPart, KfPart;
			// Real32_t PxPart, PyPart, PzPart, EePart;
			
			// acdphot.inc
			// INTEGER MAXPHOT
			// PARAMETER(MAXPHOT=12)
			// INTEGER NPHOT,KFPHOT
			// REAL PXPHOT,PYPHOT,PZPHOT,EEPHOT
			// COMMON /ACDPHOT/ NPHOT,KFPHOT(MAXPHOT),
			// +                 PXPHOT(MAXPHOT),PYPHOT(MAXPHOT),PZPHOT(MAXPHOT),
			// +                 EEPHOT(MAXPHOT) 
			
			// const Int32_t MAX_PHOT = 12;
			Int32_t NPhot, KfPhot;
			// Real32_t PxPhot, PyPhot, PzPhot, EePhot;

			// acdmiss.inc
			// REAL PXMISS,PYMISS,PXNUES,PYNUES,PXCALO,PYCALO
			// COMMON /ACDMISS/PXMISS,PYMISS,PXNUES,PYNUES,PXCALO,PYCALO
			Real32_t PxMiss, PyMiss, PxNues, PyNues, PxCalo, PyCalo; 

			// acdnout.inc
			// INTEGER NINP,NOUT
			// COMMON /ACDNOUT/ NINP,NOUT
			Int32_t NInp, NOut;

			// acmcevent.inc
			// INTEGER N, K
			// REAL P, V
			// COMMON /ACMCEVENT/ N,K(10000,5),P(10000,5),V(10000,5)
			Int32_t N,K;
			Real32_t P,V;
*/
		public:
			vector<Particle> parts;

		public:
			/*
			 * Constructor
			 * creates an input record from given particle list
			 */
			InputRecord(const vector<Particle>&);
			
			/*
			 * Access the particle list
			 */
			const vector<Particle>& particles() const;

		};

	}
}

#endif
