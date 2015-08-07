#include "Calibration.h"
#include <cstdio>

using namespace AcerDet::analyse;

Calibration::Calibration( const Configuration& config, IHistogramManager *histoMng ) :
	RCONE	( config.Cluster.ConeR ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYCAL	( config.Flag.JetCalibration ),

	IEVENT	( 0 ),

	histoManager( histoMng ),
	histoRegistered( false )
{
}

Calibration::~Calibration() {
}

void Calibration::printInfo() const {
	// print out title
	printf ("********************************************\n");
	printf ("*                                          *\n");
	printf ("*     ********************************     *\n");
	printf ("*     ***   analyse::Calibration   ***     *\n");
	printf ("*     ********************************     *\n");
	printf ("*                                          *\n");
	printf ("********************************************\n");

	// print out basic params
	printf ("\n\t jets calibration ....\n");
	printf (" calibration %s\n", KEYCAL ? "on" : "off");
	printf ("\n");
}

void Calibration::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
  Int32_t idhist = 1000 + KEYHID;
  if (!histoRegistered) {
    histoRegistered = true;
    histoManager
      ->registerHistogram(idhist+1, "Calibration: calibration correction", 50, 0.0, 2.0);
    histoManager
      ->registerHistogram(idhist+11, "Calibration: pT jets before calibration", 50, 0.0, 100.0);
    histoManager
      ->registerHistogram(idhist+12, "Calibration: pT jets after calibration", 50, 0.0, 100.0);
  }
	
  // do not use this algorithm
  if (!KEYCAL)
    return;
  
  // new event to compute
  IEVENT++;

  for (int j=0; j<orecord.Jets.size(); ++j) {

    if( orecord.Jets[j].type ==  TAU_JET ) continue;

    // histogram store
    histoManager
      ->insert(idhist+11, orecord.Jets[j].pT, weight);
      
    Real64_t X = orecord.Jets[j].pT;
    Real64_t Corr = 1.0;
    
    if( X < 40.0 ){
      Real64_t A0 =   .437876;
      Real64_t A1 =   .636209;
      Real64_t A2 =  -.709003E-01;
      Real64_t A3 =  +.320892E-02;
      Real64_t A4 =  -.660921E-04;
      Real64_t A5 =   .512782E-06;
      Corr = A0 + A1*X + A2*X*X +A3*X*X*X+A4*X*X*X*X+A5*X*X*X*X*X;

    } else if ( X > 40.0 && X < 60.0 ){
      Real64_t      A0 = 286.581;
      Real64_t      A1 =-29.0227;
      Real64_t      A2 = 1.17481;
      Real64_t      A3 =-0.0236602;
      Real64_t      A4 = 0.237002E-03;
      Real64_t      A5 =-9.44527E-07;
      Corr = A0+A1*X+A2*X*X+A3*X*X*X+A4*X*X*X*X+A5*X*X*X*X*X;   
      
    } else if ( X < 200.0){
      Real64_t   A0 =   1.05187;
      Real64_t   A1 =   -.103818E-02;
      Real64_t   A2 =    .340883E-05;
      Corr = A0+A1*X+A2*X*X;
    }
    orecord.Jets[j].pT = Corr * orecord.Jets[j].pT;

    // histogram store
    histoManager
      ->insert(idhist+1, Corr, weight);

    histoManager
      ->insert(idhist+12, orecord.Jets[j].pT, weight);
     
  }

}

void Calibration::printResults() const {
	if (KEYCAL) {
		printf ("***********************************\n");
		printf ("*                                 *\n");
		printf ("*     ***********************     *\n");
		printf ("*     ***   Output from   ***     *\n");
		printf ("*     ***  analyse::Calib ***     *\n");
		printf ("*     ***********************     *\n");
		printf ("*                                 *\n");
		printf ("***********************************\n");
	
		printf (" Analysed records: %d\n", IEVENT);
	}
}
