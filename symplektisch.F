c     =================================================================

      PROGRAM symplektisch

c     =================================================================
c     Dieses Programm loesst die Bewegungsgleichung 
c     \vec{a} = -G M \vec{r} / |r|^3
c     fuer die Bewegung um einen massereichen Zentralkoerper mit
c     einem symplektischen Verfahren 8.Ordnung. Wir rechnen der Einfachheit
c     halber mit dimensionslosen Groessen und setzen G=M=1
c     -----------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION :: r (1:3,4) ! Ortsvektor
      DOUBLE PRECISION :: v (1:3,4) ! Geschwindigkeitsvektor
      DOUBLE PRECISION :: r_alt (1:3,4) 
      DOUBLE PRECISION :: v_alt (1:3,4)
      DOUBLE PRECISION :: dr (1:3)
      DOUBLE PRECISION :: m (4)

      DOUBLE PRECISION, PARAMETER :: dt = 2d2 ! Zeitschritt
      
c     Konstanten fuer Yoshidas Verfahren 8.Ordnung
      DOUBLE PRECISION :: coeff (0:7) 

      DOUBLE PRECISION :: zeit

      INTEGER :: n_zeitschritt, n2  ! Laufindex
      INTEGER :: j,l,k1,k2


c     Koeffizienten fuer optimales Verfahren 6.Ordnung
c$$$      coeff(3) =-0.117767998417887d1
c$$$      coeff(2) = 0.235573213359357d0
c$$$      coeff(1) = 0.784513610477560d0
c$$$      coeff(0) = 1.0d0 - 2.0d0 * (coeff(1) + coeff(2) + coeff(3))
      
c     Koeffizienten fuer optimales Verfahren 8.Ordnung
      coeff(7) = 0.74167036435061295344822780d0
      coeff(6) =-0.40910082580003159399730010d0
      coeff(5) = 0.19075471029623837995387626d0  
      coeff(4) =-0.57386247111608226665638773d0
      coeff(3) = 0.29906418130365592384446354d0
      coeff(2) = 0.33462491824529818378495798d0
      coeff(1) = 0.31529309239676659663205666d0
      coeff(0) = 1.0d0 - 2.0d0 * (coeff(1) + coeff(2) +  coeff(3) +
     &     coeff(4) + coeff(5) + coeff (6) + coeff (7))


c     -----------------------------------------------------------------
c     Beginn des eigentlichen Progamms
c     -----------------------------------------------------------------
c     Anfangsbedingungen:
      r(:,:) = 0.0d0
      v(:,:) = 0.0d0

      r(1,1) = 0.0d0 ; r(2,1) = 0.0d0 ! Kepler 47A (1. Stern)
      v(1,1) = 0.0d0 ; v(2,1) = -0.347d4*5.47d0
      r(1,2) = 0.0836d0*(1.0d0+0.0234d0)*1.496d11 ; r(2,2) = 0.0d0 !Kepler 47B (2. Stern)
      v(1,2) = 0.0d0 ; v(2,2) = 1.0d4*5.47d0
      r(1,3) = 0.989d0*1.496d11 ; r(2,3) = 0.0d0 ! Kepler 47c (‰uﬂerer Planet)
      v(1,3) = 0.0d0 ; v(2,3) = 3.625d4
      r(1,4) = -0.2956d0*1.496d11+2.0d9 ; r(2,4) = 0.0d0 ! Kepler 47b (innerer Planet)
      v(1,4) = 0.0d0 ; v(2,4) = -6.344d4
      m(1) = 1.043*1.989d30; m(2) = 0.362*1.989d30
      m(3) = 0.5d-12; m(4) = 0.5d-12

      zeit = 0.0d0

      m(:) = m(:)*6.673848d-11

      OPEN (1, FILE = 'ausgabe.dat', FORM = 'FORMATTED')
      WRITE(1,'(4a16)') 'Zeit','x','y','z'
      CLOSE(1)
      
      PRINT *,SUM(m(:)*v(1,:)),SUM(m(:)*v(2,:)),SUM(m(:)*v(3,:))
      v(1,:)=v(1,:)-SUM(m(:)*v(1,:))/SUM(m(:))
      v(2,:)=v(2,:)-SUM(m(:)*v(2,:))/SUM(m(:))
      v(3,:)=v(3,:)-SUM(m(:)*v(3,:))/SUM(m(:))

      DO n_zeitschritt = 1, 40000
         DO n2 = 1, 50

         r_alt(:,:)=r(:,:)          ! FORTRAN erlaubt uns oft eine aehnlich kompakte
         v_alt(:,:)=v(:,:)          ! Notation wie die mathematische Vektorschreibweise

c     -----------------------------------------------------------------
c     Verfahren 8.Ordnung
         DO j = 3, 0, -1
c     Einzelschritte wie im Leapfrog-Verfahren, aber mit dt'= dt * c_j
            r(:,:) = r(:,:) + coeff(j) * 0.5d0 * dt * v(:,:)
            DO k1 = 1, 4
               DO k2 = 1, 4
                  IF (k2 .NE. k1) THEN
                     dr (:) = r(:,k2) - r(:,k1)
                     v(:,k1) = v(:,k1) + coeff(j) * dt * dr (:) * 
     &                    m (k2) /
     &                    SQRT(dr(1) ** 2 + dr(2) ** 2 + 
     &                    dr(3) ** 2) ** 3
                  END IF
               END DO
            END DO
            r(:,:) = r(:,:) + coeff(j) * 0.5d0 * dt * v(:,:)
         END DO
c     2.Teil: Das Ganze mit den Koeffizienten c_j in umgekehrter Folge
         DO j = 1, 3
            r(:,:) = r(:,:) + coeff(j) * 0.5d0 * dt * v(:,:)
            DO k1 = 1, 4
               DO k2 = 1, 4
                  IF (k2 .NE. k1) THEN
                     dr (:) = r(:,k2) - r(:,k1)
                     v(:,k1) = v(:,k1) + coeff(j) * dt * dr (:) * 
     &                    m (k2) /
     &                    SQRT(dr(1) ** 2 + dr(2) ** 2 + 
     &                    dr(3) ** 2) ** 3
                  END IF
               END DO
            END DO
            r(:,:) = r(:,:) + coeff(j) * 0.5d0 * dt * v(:,:)
         END DO
c     -----------------------------------------------------------------

         zeit = zeit + dt

         END DO

c     Ausgabe des Ortsvektors als Funktion der Zeit
         WRITE(1,'(4e16.7)') zeit, r(1,1), r(2,1), r(3,1)
         WRITE(2,'(4e16.7)') zeit, r(1,2), r(2,2), r(3,2)
         WRITE(3,'(4e16.7)') zeit, r(1,3), r(2,3), r(3,3)
         WRITE(4,'(4e16.7)') zeit, r(1,4), r(2,4), r(3,4)

      END DO

      CLOSE (1) !Ausgabedatei schlieﬂen
      CLOSE (2)
      CLOSE (3)
      CLOSE (4)

      END PROGRAM symplektisch

c     =================================================================
