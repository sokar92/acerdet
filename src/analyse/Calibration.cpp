#include "Calibration.h"
#include <cstdio>

using namespace AcerDet::analyse;

Calibration::Calibration( const Configuration& config ) :
	RCONE	( config.Cluster.ConeR ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYCAL	( config.Flag.JetCalibration ),

	IEVENT	( 0 ),
	IDENT	( 1000 + KEYHID )
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

void Calibration::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
/*
      IF(KEYCAL.NE.0) THEN
c histograms initialization
      NBINA = 50
      TMAXA = 300.0
      TMINA = 0.0
      CALL HBOOK1(IDENT+10,'ACDCAL: calibration entry$',
     #                                     NBINA,TMINA,TMAXA,0.0)       
      CALL HBOOK1(IDENT+11,'ACDCAL: calibration factor$',
     #                                     NBINA,TMINA,TMAXA,0.0)       

      NBINA = 50
      TMAXA =   2.0
      TMINA =   0.0
      DO I=0,8
      CALL HBOOK1(IDENT+20+I,'ACDCAL: pTjetCALIB/pTparton $',
     #                                     NBINA,TMINA,TMAXA,0.)       
      CALL HBOOK1(IDENT+30+I,'ACDCAL: pTjetNONCALIB/pTparton $',
     #                                     NBINA,TMINA,TMAXA,0.)       
      ENDDO
 
      ENDIF

      ELSEIF(MODE.EQ.0.AND.KEYCAL.NE.0) THEN
C     ========================================
      IEVENT = IEVENT+1

      DO I=1,NJET
        PTJETCRU(I)=PJET(I,5)
        X=PJET(I,5)
        IF(X.LT.10) THEN
           CORE = 0.0
        ELSEIF(X.LE.40.0.AND.X.GE.10.) THEN
           A0 =   1.2715
           A1 =   .12241
           A2 =  -.10480E-01
           A3 =  +.33310E-03
           A4 =  -.47454E-05
           A5 =   .25436E-07
           CORE = A0+A1*X+A2*X**2+A3*X**3+A4*X**4+A5*X**5
           CORE = CORE*1.006
        ELSEIF (X.GT.40.0.AND.X.LT.60.0) THEN
           A0 = 11.3579
           A1 =-0.993207
           A2 = 0.0388214
           A3 =-0.000754235
           A4 =7.22229E-06
           A5 =-2.71501E-08
           CORE = A0+A1*X+A2*X**2+A3*X**3+A4*X**4+A5*X**5                                        ELSEIF(X.LT.200.0)THEN
           A0 =   1.18
           A1 =   -.16672E-02
           A2 =    .44414E-05
           CORE = A0+A1*X+A2*X**2  
        ELSE
           XX =   200.0
           A0 =   1.18
           A1 =   -.16672E-02                                                            
           A2 =    .44414E-05   
           CORE = A0+A1*XX+A2*XX**2
        ENDIF
        CALL HF1(IDENT+11,X,CORE)
        CALL HF1(IDENT+10,X,1.0)
        PJET(I,5) = PJET(I,5) * CORE
      ENDDO

      DO IJET=1,NJET
      PTREC=0
      DETRMIN=RCONE
      DO I=7,N 
      IF(K(I,1).NE.21.OR.ABS(K(I,2)).GT.10) GOTO 1040
      PT =SQRT(P(I,1)**2+P(I,2)**2)
      ETA=SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
      PHI=ANGLE(P(I,1),P(I,2))
      DPHIA=ABS(PJET(IJET,4)-PHI) 
      IF(DPHIA.GT.PI) DPHIA=DPHIA-2*PI
      DETR=SQRT((ETA-PJET(IJET,3))**2+DPHIA**2)
      IF(DETR.GT.DETRMIN) GOTO 1040 
      PTREC = PT
      DETRMIN=DETR
 1040 CONTINUE
      ENDDO 
      IF(PTREC.NE.0) THEN
         CALL HF1(IDENT+20,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).LT.20)
     #        CALL HF1(IDENT+21,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.20.AND.PTJETCRU(IJET).LT.25)
     #        CALL HF1(IDENT+22,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.25.AND.PTJETCRU(IJET).LT.30)
     #        CALL HF1(IDENT+23,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.30.AND.PTJETCRU(IJET).LT.40)
     #        CALL HF1(IDENT+24,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.40.AND.PTJETCRU(IJET).LT.60)
     #        CALL HF1(IDENT+25,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.60.AND.PTJETCRU(IJET).LT.100)
     #        CALL HF1(IDENT+26,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.100.AND.PTJETCRU(IJET).LT.200)
     #        CALL HF1(IDENT+27,PJET(IJET,5)/PTREC,1.0)
         IF(PTJETCRU(IJET).GT.200)
     #        CALL HF1(IDENT+28,PJET(IJET,5)/PTREC,1.0)
         CALL HF1(IDENT+30,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).LT.30)
     #        CALL HF1(IDENT+31,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.20.AND.PTJETCRU(IJET).LT.25)
     #        CALL HF1(IDENT+32,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.25.AND.PTJETCRU(IJET).LT.30)
     #        CALL HF1(IDENT+33,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.30.AND.PTJETCRU(IJET).LT.40)
     #        CALL HF1(IDENT+34,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.40.AND.PTJETCRU(IJET).LT.60)
     #        CALL HF1(IDENT+35,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.60.AND.PTJETCRU(IJET).LT.100)
     #        CALL HF1(IDENT+36,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GE.100.AND.PTJETCRU(IJET).LT.200)
     #        CALL HF1(IDENT+37,PTJETCRU(IJET)/PTREC,1.0)
         IF(PTJETCRU(IJET).GT.200)
     #        CALL HF1(IDENT+38,PTJETCRU(IJET)/PTREC,1.0)
      ENDIF
      ENDDO
*/
}
