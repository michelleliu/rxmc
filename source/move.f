      Subroutine Move(Av1,Av2,Delta)
      Implicit None

      Include 'commons.inc'

C     Displace A Randomly Selected Particle

      Logical Laccept
      Integer Ib,Ipart
      Double Precision Xi,Yi,Zi,Rxo,Ryo,Rzo,Randomnumber,Unew,Uold
     $     ,Virnew,Virold,Av1,Av2,Delta
     $     ,Enew,Vnew,DiffE,DiffV

      If(Npart.Eq.0) Return

      Ipart = 1 + Int(Dble(Npart)*Randomnumber())
      Ib    = Ibox(Ipart)

      Xi = Rx(Ipart) + (2.0d0*Randomnumber()-1.0d0)*Delta
      Yi = Ry(Ipart) + (2.0d0*Randomnumber()-1.0d0)*Delta
      Zi = Rz(Ipart) + (2.0d0*Randomnumber()-1.0d0)*Delta

C     Put Back In The Box

      If(Xi.Lt.0.0d0) Then
         Xi = Xi + Box(Ib)
      Elseif(Xi.Gt.Box(Ib)) Then
         Xi = Xi - Box(Ib)
      Endif

      If(Yi.Lt.0.0d0) Then
         Yi = Yi + Box(Ib)
      Elseif(Yi.Gt.Box(Ib)) Then
         Yi = Yi - Box(Ib)
      Endif

      If(Zi.Lt.0.0d0) Then
         Zi = Zi + Box(Ib)
      Elseif(Zi.Gt.Box(Ib)) Then
         Zi = Zi - Box(Ib)
      Endif

      Rxo = Rx(Ipart)
      Ryo = Ry(Ipart)
      Rzo = Rz(Ipart)

      Call Epart(Ib,Virold,Uold,Rx(Ipart),Ry(Ipart),Rz(Ipart),Ipart,Types(Ipart))
      Call Epart(Ib,Virnew,Unew,Xi,Yi,Zi
     &     ,Ipart,Types(Ipart))

      Rx(Ipart) = Xi
      Ry(Ipart) = Yi
      Rz(Ipart) = Zi

      Call Etot(Ib,Vnew,Enew)

      DiffE=Unew-Uold-Enew+Etotal(Ib)
      DiffV=Virnew-Virold-Vnew+Vtotal(Ib)
      If(DiffE.Gt.1e-9) Then
         Stop "ERROR Energy Difference !!! "
      Endif


      If((-Beta*(Enew-Etotal(Ib))>0)) Then
         Laccept=.True.
      Else
         Call Accept(Dexp(-Beta*(Enew-Etotal(Ib))),Laccept)
      Endif

      Av2 = Av2 + 1.0d0

C     Accept Or Reject

      If(Laccept) Then
         Av1 = Av1 + 1.0d0

         Etotal(Ib) = Enew
         Vtotal(Ib) = Vnew

      Else

C     Reject, Restore Old Configuration
         Rx(Ipart) = Rxo
         Ry(Ipart) = Ryo
         Rz(Ipart) = Rzo

      Endif

      Return
      End
