      Subroutine Epair(Ib,Vir,Upot,Xi,Yi,Zi,Xj,Yj,Zj,Mytype,Yourtype)
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Of Particle Ipart In Box Ib

      Integer Ib,Mytype,Yourtype
      Double Precision Vir,Upot,Dx,Dy,Dz,R2,S2,Bx,Hbx,Xi,Yi,Zi,Xj,Yj,Zj
      Double Precision Radius

      Upot = 0.0d0
      Vir  = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

      Dx = Xj-Xi
      Dy = Yj-Yi
      Dz = Zj-Zi

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
            S2   = (Sig(Yourtype,Mytype))**2 ! FIXME
            R2   = S2*1.0d0/R2
            R2   = R2*R2*R2
            Upot = Upot + 4.0d0*Eps(Yourtype,Mytype)*R2*(R2-1.0d0)
     &           - Ecut(Yourtype,Mytype)
            Vir  = Vir  + 48.0d0*Eps(Yourtype,Mytype)*R2*(R2-0.5d0)
         Endif
      Else If(Potential.Eq.2) Then
         Stop "Error Potential Tail Correction !!!"
      Else If(Potential.Eq.0) Then
         Radius  = 0.5d0*(Eps(Yourtype,Yourtype)+Eps(Mytype,Mytype))
         If(R2.Lt.Radius) Then
            Upot = Upot + 1000000000000d0
            Vir  = Vir + 1000000000000d0
         Else
            Upot = Upot
            Vir  = Vir
         Endif
      Else
         Stop "Error potential!!"
      Endif

      Return
      End
