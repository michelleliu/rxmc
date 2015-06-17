      Subroutine Tailc(Ib,Mytype,Utail)
      Implicit None

      Include 'commons.inc'

C     Compute The Tail Correction Energy Of A Particle

      Integer Ib,Mytype,I,J
      Double Precision Utail,Vol,R3,S3

      Utail = 0.0d0

      Vol=Dble(Box(Ib)**3.0d0)

      R3 = Rcutsq**(1.0d0/2.0d0)
      R3 = Rcutsq*R3

      If(Mytype.Eq.-1) Then
         Do I=1,Ntype
            Do J=1,Nytype
               S3 = Sig(I,J)**3.0d0
               R3 = S3*1.0d0/R3
               R3 = R3*(R3*R3/3.0d0-1.0d0)
               Utail = Utail + (8.0d0/(3.0d0*Vol))*4.0d0*Datan(1.0d0)
     &              *Npboxtype(I)*Npboxtype(J)*Eps(I,J)*S3*R3
            Enddo
         Enddo
      Else If(Mytype.Gt.0.And.Mytype.Le.Ntype) Then
         Do I=1,Ntype
            S3 = Sig(I,Mytype)**3.0d0
            R3 = S3*1.0d0/R3
            R3 = R3*(R3*R3/3.0d0-1.0d0)
            If(Mytype.Ne.I) Then
               Utail = Utail + (16.0d0/(3.0d0*Vol))*4.0d0*Datan(1.0d0)
     &           *Npboxtype(I)*Eps(I,Mytype)*S3*R3
            Else
               Utail = Utail + (8.0d0/(3.0d0*Vol))*4.0d0*Datan(1.0d0)
     &           *(2.0d0*Npboxtype(I)-1.0d0)*Eps(I,Mytype)*S3*R3
            Endif

         Enddo
      Else
         Stop "Error Mytype Tailc !!!"
      Endif

      Return
      End
