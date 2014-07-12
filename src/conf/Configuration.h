#ifndef _AcerDet_Configuration_
#define _AcerDet_Configuration_

#pragma once

#include <string>
using namespace std;

namespace AcerDet {
	namespace conf {

		struct Configuration {
			struct _Flag {
				int HistogramId;
				bool Smearing;
				bool BField;
				int SusyParticle;

				bool BCJetsLabeling;
				bool TauJetsLabeling;
				bool JetCalibration;

				_Flag();

				void read( const string& );
				string write( ) const;

			} Flag;

			struct _Cell {
				double RapidityCoverage;
				double MinpT;
				double MinEt;
				double EtaTransition;
				double GranularityEta;
				double GranularityPhi;

				_Cell();

				void read( const string& );
				string write( ) const;

			} Cell;

			struct _Cluster {
				double RapidityCoverage;
				double ConeR;
				double MinEt;
				double MinEtInit;

				_Cluster();

				void read( const string& );
				string write( ) const;

			} Cluster;

			struct _Muon {
				double MinMomenta;
				double MaxEta;
				double MinIsolRlj;
				double ConeR;
				double MaxEnergy;

				_Muon();

				void read( const string& );
				string write( ) const;

			} Muon;

			struct _Photon {
				double MinMomenta;
				double MaxEta;
				double MinJetsRlj;
				double MinIsolRlj;
				double ConeR;
				double MaxEnergy;

				_Photon();

				void read( const string& );
				string write( ) const;

			} Photon;

			struct _Electron {
				double MinMomenta;
				double MaxEta;
				double MinJetsRlj;
				double MinIsolRlj;
				double ConeR;
				double MaxEnergy;

				_Electron();

				void read( const string& );
				string write( ) const;

			} Electron;

			struct _Jet {
				double RapidityCoverage;
				double MinEnergy;

				_Jet();

				void read( const string& );
				string write( ) const;

			} Jet;

			struct _BJet {
				double MinMomenta;
				double MaxEta;
				double MaxRbj;

				_BJet();

				void read( const string& );
				string write( ) const;

			} BJet;

			struct _CJet {
				double MinMomenta;
				double MaxEta;
				double MaxRcj;

				_CJet();

				void read( const string& );
				string write( ) const;

			} CJet;

			struct _Tau {
				double MinpT;
				double MaxEta;
				double MinR;
				double MaxR;

				_Tau();

				void read( const string& );
				string write( ) const;

			} Tau;

			struct _Misc {
				double MinEt;

				_Misc();

				void read( const string& );
				string write( ) const;

			} Misc;

			/* save configuration in specific file */
			static void save( const Configuration&, const std::string& );

			/* reads configuration from given file */
			static Configuration fromFile( const std::string& );

			/* returns default configuration */
			static Configuration getDefault();
		};

	}
}

#endif
