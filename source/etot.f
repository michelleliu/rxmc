      Subroutine Etot(Ib,Vir,Upot)
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Total Of Box Ib

      Integer I,J,Ib
      Double Precision Vir,Upot,Dx,Dy,Dz,R2,S2,Bx,Hbx
      Double Precision Radius,Utail

      Upot = 0.0d0
      Vir  = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

      Do I=1,Npart-1

         If(Ibox(I).Eq.Ib) Then

            Do J=I+1,Npart

               If(Ibox(J).Eq.Ib) Then

C                 Minimum Image Convention
                  Dx = Rx(I)-Rx(J)
                  Dy = Ry(I)-Ry(J)
                  Dz = Rz(I)-Rz(J)

                  If (Dx.Gt.Hbx) Then
                     Dx = Dx - Bx
                  Elseif (Dx.Lt.-Hbx) Then
                     Dx = Dx + Bx
                  Endif

                  If (Dy.Gt.Hbx) Then
                     Dy = Dy - Bx
                  Elseif (Dy.Lt.-Hbx) Then
                     Dy = Dy + Bx
                  Endif

                  If (Dz.Gt.Hbx) Then
                     Dz = Dz - Bx
                  Elseif (Dz.Lt.-Hbx) Then
                     Dz = Dz + Bx
                  Endif

                  R2 = Dx**2 + Dy**2 + Dz**2

                  If(Potential.Eq.1) Then
                     If(R2.Lt.Rcutsq) Then
                        S2   = (Sig(Types(I),Types(J)))**2
                        R2   = S2*1.0d0/R2
                        R2   = R2*R2*R2
                        Upot = Upot + 4.0d0*Eps(Types(I),Types(J))*R2*(R2-1.0d0)
     &                        - Ecut(Types(I),Types(J))
                        Vir  = Vir  + 48.0d0*Eps(Types(I),Types(J))*R2*(R2-0.5d0)
                     Endif
                  Else If(Potential.Eq.2) Then
                     If(R2.Lt.Rcutsq) Then
                        S2   = (Sig(Types(I),Types(J)))**2
                        R2   = S2*1.0d0/R2
                        R2   = R2*R2*R2
                        Upot = Upot + 4.0d0*Eps(Types(I),Types(J))*R2*(R2-1.0d0)
                        Vir  = Vir  + 48.0d0*Eps(Types(I),Types(J))*R2*(R2-0.5d0)
                     Endif
                  Else If(Potential.Eq.0) Then
                     Radius  = 0.5d0*(Eps(Types(I),Types(I))+Eps(Types(J),Types(J)))
                     If(R2.Lt.Radius) Then
                        Upot = Upot + 1000000000000d0
                        Vir  = Vir + 1000000000000d0
                     Else
                        Upot = Upot
                        Vir  = Vir
                     Endif
                  Else
                     Stop "Error Potential!!!"
                  Endif
               Endif
            Enddo
         Endif
      Enddo

C     Tail correction of box Ib
      If(Potential.Eq.2) Then
        Call Tailc(Ib,-1,Utail)
        Upot = Upot + Utail
      Endif

      Return
      End
