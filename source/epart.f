      Subroutine Epart(Ib,Vir,Upot,Xi,Yi,Zi,Ipart,Mytype)
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Of Particle Ipart In Box Ib

      Integer I,Ib,Ipart,Mytype
      Double Precision Vir,Upot,Dx,Dy,Dz,R2,S2,Bx,Hbx,Xi,Yi,Zi
      Double Precision Radius,Utail
      Logical Overlap

      Upot = 0.0d0
      Vir  = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

      Overlap = .False.

      If(Potential.Eq.2.And.Ntype.Eq.1) Then
C        Add Tail Correction
         Call Tailc(Ib,1,Utail)
         Upot = Upot + Utail
      Endif
      Do I=1,Npart

         If(.Not.Overlap) Then
            If(Ibox(I).Eq.Ib.And.I.Ne.Ipart) Then

               Dx = Rx(I)-Xi
               Dy = Ry(I)-Yi
               Dz = Rz(I)-Zi

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
                     S2   = (Sig(Types(I),Mytype))**2
                     R2   = S2*1.0d0/R2
                     R2   = R2*R2*R2
                     Upot = Upot + 4.0d0*Eps(Types(I),Mytype)*R2*(R2-1.0d0)
     &                    - Ecut(Types(I),Mytype)
                     Vir  = Vir  + 48.0d0*Eps(Types(I),Mytype)*R2*(R2-0.5d0)
                  Endif
               Else If(Potential.Eq.2) Then
                  If(R2.Lt.Rcutsq) Then
                     S2   = (Sig(Types(I),Mytype))**2
                     R2   = S2*1.0d0/R2
                     R2   = R2*R2*R2
                     Upot = Upot + 4.0d0*Eps(Types(I),Mytype)*R2*(R2-1.0d0)
                     Vir  = Vir  + 48.0d0*Eps(Types(I),Mytype)*R2*(R2-0.5d0)
                  Endif
               Else If(Potential.Eq.0) Then
                  Radius  = 0.5d0*(Eps(Types(I),Types(I))+Eps(Mytype,Mytype))
                  If(R2.Lt.Radius) Then
                     Overlap = .True.
                     Upot = Upot + 1000000000000d0
                     Vir  = Vir + 1000000000000d0
                  Else
                     Upot = Upot
                     Vir  = Vir
                  Endif
               Else
                  Stop "Error Potential!"
               Endif
            Endif
         Endif
      Enddo

      Return
      End
