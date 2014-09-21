#ifndef _AcerDet_Configuration_
#define _AcerDet_Configuration_

#pragma once

#include <string>
using namespace std;

namespace AcerDet {
	namespace conf {

		//! Global configuration
		/*!
		 * A set of parameters used by analyse algorithms.
		 */
		struct Configuration {
			
			//! Common flags
			/*!
			 * A set of commonly used parameters.
			 */
			struct _Flag {
				int HistogramId;	/*!< detailed description  */
				bool Smearing;		/*!< detailed description  */
				bool BField;		/*!< detailed description  */
				int SusyParticle;	/*!< detailed description  */

				bool BCJetsLabeling;	/*!< detailed description  */
				bool TauJetsLabeling;	/*!< detailed description  */
				bool JetCalibration;	/*!< detailed description  */

				_Flag();

				void read( const string& );
				string write( ) const;

			} Flag;

			//! Cell-specific flags
			/*!
			 * A set of parameters used by analyse::Cell algorithm.
			 */
			struct _Cell {
				double RapidityCoverage;	/*!< detailed description  */
				double MinpT;			/*!< detailed description  */
				double MinEt;			/*!< detailed description  */
				double EtaTransition;	/*!< detailed description  */
				double GranularityEta;	/*!< detailed description  */
				double GranularityPhi;	/*!< detailed description  */

				_Cell();

				void read( const string& );
				string write( ) const;

			} Cell;

			//! Cluster-specific flags
			/*!
			 * A set of parameters used by analyse::Cluster algorithm.
			 */
			struct _Cluster {
				double RapidityCoverage;	/*!< detailed description  */
				double ConeR;		/*!< detailed description  */
				double MinEt;		/*!< detailed description  */
				double MinEtInit;	/*!< detailed description  */

				_Cluster();

				void read( const string& );
				string write( ) const;

			} Cluster;

			//! Muon-specific flags
			/*!
			 * A set of parameters used by analyse::Muon algorithm.
			 */
			struct _Muon {
				double MinMomenta;	/*!< detailed description  */
				double MaxEta;		/*!< detailed description  */
				double MinIsolRlj;	/*!< detailed description  */
				double ConeR;		/*!< detailed description  */
				double MaxEnergy;	/*!< detailed description  */

				_Muon();

				void read( const string& );
				string write( ) const;

			} Muon;

			//! Photon-specific flags
			/*!
			 * A set of parameters used by analyse::Photon algorithm.
			 */
			struct _Photon {
				double MinMomenta;	/*!< detailed description  */
				double MaxEta;		/*!< detailed description  */
				double MinJetsRlj;	/*!< detailed description  */
				double MinIsolRlj;	/*!< detailed description  */
				double ConeR;		/*!< detailed description  */
				double MaxEnergy;	/*!< detailed description  */

				_Photon();

				void read( const string& );
				string write( ) const;

			} Photon;

			//! Electron-specific flags
			/*!
			 * A set of parameters used by analyse::Electron algorithm.
			 */
			struct _Electron {
				double MinMomenta;	/*!< detailed description  */
				double MaxEta;		/*!< detailed description  */
				double MinJetsRlj;	/*!< detailed description  */
				double MinIsolRlj;	/*!< detailed description  */
				double ConeR;		/*!< detailed description  */
				double MaxEnergy;	/*!< detailed description  */

				_Electron();

				void read( const string& );
				string write( ) const;

			} Electron;

			//! Jet-specific flags
			/*!
			 * A set of parameters used by analyse::Jet algorithm.
			 */
			struct _Jet {
				double RapidityCoverage;	/*!< detailed description  */
				double MinEnergy;	/*!< detailed description  */

				_Jet();

				void read( const string& );
				string write( ) const;

			} Jet;

			//! BJet-specific flags
			/*!
			 * A set of parameters used by analyse::BJet algorithm.
			 */
			struct _BJet {
				double MinMomenta;	/*!< detailed description  */
				double MaxEta;		/*!< detailed description  */
				double MaxRbj;		/*!< detailed description  */

				_BJet();

				void read( const string& );
				string write( ) const;

			} BJet;

			//! CJet-specific flags
			/*!
			 * A set of parameters used by analyse::CJet algorithm.
			 */
			struct _CJet {
				double MinMomenta;	/*!< detailed description  */
				double MaxEta;		/*!< detailed description  */
				double MaxRcj;		/*!< detailed description  */

				_CJet();

				void read( const string& );
				string write( ) const;

			} CJet;

			//! Tau-specific flags
			/*!
			 * A set of parameters used by analyse::Tau algorithm.
			 */
			struct _Tau {
				double MinpT;	/*!< detailed description  */
				double MaxEta;	/*!< detailed description  */
				double MinR;	/*!< detailed description  */
				double MaxR;	/*!< detailed description  */

				_Tau();

				void read( const string& );
				string write( ) const;

			} Tau;

			//! Mis-specific flags
			/*!
			 * A set of parameters used by analyse::Mis algorithm.
			 */
			struct _Misc {
				double MinEt;	/*!< detailed description  */

				_Misc();

				void read( const string& );
				string write( ) const;

			} Misc;

			//! Save configuration to specific file
			/*!
			 * \param configuration configuration to save
			 * \param fileName name of destination file
			 */
			static void save(
				const Configuration& configuration,
				const std::string& fileName );

			//! Read configuration from given file
			/*!
			 * \param fileName name of file to read from
			 */
			static Configuration fromFile(
				const std::string& fileName );

			//! Returns default configuration
			static Configuration getDefault();
		};

	}
}

#endif
