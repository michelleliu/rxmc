      Subroutine Tailc(Ib,Mytype,Utail)
      Implicit None

      Include 'commons.inc'

C     Compute The Tail Correction Energy Of A Particle

      Integer Ib,Mytype,I
      Double Precision Utail,Rho(Ntype),R3,S3

      Utail = 0.0d0

      Do I=1,Ntype
         Rho(I) = Dble(Npboxtype(Ib,I))/Dble(Box(Ib)**3)
      Enddo

      R3 = Rcutsq**(1.0d0/2.0d0)
      R3 = Rcutsq*R3

      Do I=1,Ntype
         S3 = Sig(I,Mytype)**3
         R3 = S3*1.0d0/R3
         R3 = S3*R3*(R3*R3/3.0d0-R3)
         Utail = Utail + (8.0d0/3.0d0)*Pi*Rho(I)*Eps(I,Mytype)*R3
      Enddo

      Return
      End
