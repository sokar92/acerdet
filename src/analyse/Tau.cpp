#include "Tau.h"
#include <cstdio>

using namespace AcerDet::analyse;

Tau::Tau( const Configuration& config ) :
	ETJET	( config.Jet.MinEnergy ),

	PTTAU	( config.Tau.MinpT ),
	ETATAU	( config.Tau.MaxEta ),
	RJTAU	( config.Tau.MinR ),
	PTFRAC	( config.Tau.MaxR ),

	KEYHID	( config.Flag.HistogramId ),
	KEYTAU	( config.Flag.TauJetsLabeling ),

	IEVENT	( 0 ),
	
	histo_jets	("Tau: jets-multiplicity", 0.0, 10.0, 10),
	histo_taus	("Tau: multiplicity", 0.0, 10.0, 10)
{}

Tau::~Tau() {}

void Tau::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Tau   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t... tau-jets labeling ...\n");
	printf ("labeling on/off %s\n", KEYTAU ? "on" : "off");
	if (KEYTAU) {
		printf ("\ttau-jets ...\n");
		printf ("min tau-had p_T %lf\n", PTTAU);
		printf ("max tau-had eta %lf\n", ETATAU);
		printf ("max R_tauj for tau-jets %lf\n", RJTAU);
		printf ("tau-had frac. of jet %lf\n", PTFRAC);
	}
	printf ("\n");
}

void Tau::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Tau: analyse record\n");
	/*
      SUBROUTINE ACDTAU(MODE,LPAR,YPAR)

      INCLUDE 'acerdet_inc/acdnout.inc'
      INCLUDE 'acerdet_inc/acmcevent.inc'
      INCLUDE 'acerdet_inc/jetall.inc'
      INTEGER KEYHID, KEYTAU, IEVENT 
      INTEGER I, II, NSTART,NSTOP
      REAL  DR, DDR
      REAL ANGLE
      REAL ETA, PHI, PT, PTFRAC
      REAL RJTAU, PTTAU, ETATAU, ETJET
      INTEGER NTAU, NJETTAU, JETTAU, IN
      LOGICAL TAUJET
      INTEGER NBINA
      INTEGER IDENT
      REAl TMAXA,TMINA
      INCLUDE 'acerdet_inc/ftbxdx.inc'

      ELSEIF(MODE.EQ.0.AND.KEYTAU.NE.0) THEN
C     =======================================
      IEVENT = IEVENT+1

      NSTOP=0
      NSTART=1
      DO I=1, N
       IF(K(I,1).NE.21) THEN
           NSTOP = I-1
           NSTART= I
           GOTO 500
       ENDIF
      ENDDO
 500  CONTINUE

c....look for tau-jets
      NTAU    = 0
      NJETTAU = 0
      DO I=NSTART,N
      IF(ABS(K(I,2)).EQ.15) THEN
c....if there are still jets
        IF(NJET.GT.0) THEN
c-------
        TAUJET  = .TRUE.
c....   choose only hadronic tau-decay
        IN=0
        IF(K(I,4).EQ.0.OR.K(I,5).EQ.0) GOTO 9999
        DO II=K(I,4),K(I,5)
         IF(ABS(K(II,2)).EQ.11.OR.ABS(K(II,2)).EQ.13)TAUJET= .FALSE.
         IF(ABS(K(II,2)).EQ.16)IN=II
        ENDDO
        IF(IN.EQ.0) GOTO 9999
c....
        PT=SQRT((P(I,1)-P(IN,1))**2+(P(I,2)-P(IN,2))**2) 
        IF(PT.LT.PTTAU) TAUJET= .FALSE.
        ETA=SIGN(LOG((SQRT(PT**2+(P(I,3)-P(IN,3))**2)
     #     +ABS(P(I,3)-P(IN,3)))/PT),P(I,3)-P(IN,3)) 
        IF(ABS(ETA).GT.ETATAU) TAUJET=.FALSE.
        PHI=ANGLE(P(I,1)-P(IN,1),P(I,2)-P(IN,2))
        IF(TAUJET) NTAU=NTAU+1
c === mark tau-jet
        DR=100.0
        DO II=1,NJET
          DDR=SQRT((ETA-PJET(II,3))**2+(PHI-PJET(II,4))**2)
          IF(ABS(PHI-PJET(II,4)).GT.PI)
     #     DDR=SQRT(   (ETA-PJET(II,3))**2
     #            +(ABS(PHI-PJET(II,4))-2*PI)**2 )
          IF(DDR.LT.DR) JETTAU=II
          DR=MIN(DDR,DR)
        ENDDO
        IF(DR.GT.RJTAU) THEN
            TAUJET=.FALSE.
            JETTAU=0
        ENDIF
        IF(TAUJET.AND.ABS(PJET(JETTAU,3)).LT.2.5.AND.
     #     PT/PJET(JETTAU,5).GT.PTFRAC) THEN
            KJET(JETTAU,2)=K(I,2)
            NJETTAU=NJETTAU+1
        ENDIF
 
c-------
        ENDIF
c-------
       ENDIF
 9999  CONTINUE
       ENDDO
       CALL HF1(IDENT+11,REAL(NJETTAU),1.0)
       CALL HF1(IDENT+21,REAL(NTAU),1.0)

C......

C
      ELSEIF(MODE.EQ.1.AND.KEYTAU.NE.0) THEN
C     ========================================

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '   ACDTAU      : WINDOW A        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)


      CALL HPRINT(IDENT+11)
      CALL HPRINT(IDENT+21)

      ENDIF
C     =====

      END
*/
}

void Tau::printResults() const {
	printf ("**********************************\n");
	printf ("*                                *\n");
	printf ("*     **********************     *\n");
	printf ("*     ***  Output from   ***     *\n");
	printf ("*     ***  analyse::Tau  ***     *\n");
	printf ("*     **********************     *\n");
	printf ("*                                *\n");
	printf ("**********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo_jets.print( true );
	histo_taus.print( true );
}
