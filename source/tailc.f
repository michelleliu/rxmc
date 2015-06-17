      Subroutine Tailc(Ib,Mytype,Utail)
      Implicit None

      Include 'commons.inc'

C     Compute The Tail Correction Energy Of A Particle

C     THIS IS ONLY FOR PURE COMPONENT!! Otherwise Must Change Densities

      Integer Ib,Mytype
      Double Precision Utail,Rho,R3,S3

      Rho = Dble(Npboxtype(Ib,Mytype))/Dble(Box(Ib)**3.0d0)

      R3 = Rcutsq**(1.0d0/2.0d0)
      R3 = Rcutsq*R3

      S3 = Sig(Mytype,Mytype)**3.0d0
      R3 = S3*1.0d0/R3
      R3 = S3*R3*(R3*R3/3.0d0-1.0d0)

      Utail = (8.0d0/3.0d0)*4.0d0*Datan(1.0d0)*Rho*Eps(Mytype,Mytype)*R3

      Return
      End
