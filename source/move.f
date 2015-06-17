      Subroutine Move(Av1,Av2,Delta)
      Implicit None

      Include 'commons.inc'

C     Displace A Randomly Selected Particle

      Logical Laccept
      Integer Ib,Ipart
      Double Precision Rxtrial,Rytrial,Rztrial,Rxo,Ryo,Rzo,Randomnumber,Unew,Uold
     $     ,Virnew,Virold,Av1,Av2,Delta
     $     ,Vnew,Vold,Enew,Eold

      If(Npart.Eq.0) Return

      Ipart = 1 + Int(Dble(Npart)*Randomnumber())
      Ib    = Ibox(Ipart)

      Rxtrial = Rx(Ipart) + (2.0d0*Randomnumber()-1.0d0)*Delta
      Rytrial = Ry(Ipart) + (2.0d0*Randomnumber()-1.0d0)*Delta
      Rztrial = Rz(Ipart) + (2.0d0*Randomnumber()-1.0d0)*Delta

C     Put Back In The Box

      If(Rxtrial.Lt.0.0d0) Then
         Rxtrial = Rxtrial + Box(Ib)
      Elseif(Rxtrial.Gt.Box(Ib)) Then
         Rxtrial = Rxtrial - Box(Ib)
      Endif

      If(Rytrial.Lt.0.0d0) Then
         Rytrial = Rytrial + Box(Ib)
      Elseif(Rytrial.Gt.Box(Ib)) Then
         Rytrial = Rytrial - Box(Ib)
      Endif

      If(Rztrial.Lt.0.0d0) Then
         Rztrial = Rztrial + Box(Ib)
      Elseif(Rztrial.Gt.Box(Ib)) Then
         Rztrial = Rztrial - Box(Ib)
      Endif

      Rxo = Rx(Ipart)
      Ryo = Ry(Ipart)
      Rzo = Rz(Ipart)

      !Call Etot(Ib,Vold,Eold)

C     Put New Particle In Box
      Rx(Ipart) = Rxtrial
      Ry(Ipart) = Rytrial
      Rz(Ipart) = Rztrial

      Call Etot(Ib,Vnew,Enew)

!      Call Epart(Ib,Virold,Uold,Xi,Yi,Zi,Ipart,Types(Ipart))
!      Call Epart(Ib,Virnew,Unew,Rxtrial,Rytrial,Rztrial
!     &     ,Ipart,Types(Ipart))

      !If((-Beta*(Unew-Uold)>0)) Then
      If((-Beta*(Enew-Etotal(Ib))>0)) Then
         Laccept=.True.
      Else
         !Call Accept(Dexp(-Beta*(Unew-Uold)),Laccept)
         Call Accept(Dexp(-Beta*(Enew-Etotal(Ib))),Laccept)
      Endif

      Av2 = Av2 + 1.0d0

C     Accept Or Reject

      If(Laccept) Then
         Av1 = Av1 + 1.0d0

         Etotal(Ib) = Enew
         Vtotal(Ib) = Vnew
         !Etotal(Ib) = Etotal(Ib) + Unew   - Uold
         !Vtotal(Ib) = Vtotal(Ib) + Virnew - Virold

         !Rx(Ipart) = Rxtrial
         !Ry(Ipart) = Rytrial
         !Rz(Ipart) = Rztrial
      Else
         Rx(Ipart) = Rxo
         Ry(Ipart) = Ryo
         Rz(Ipart) = Rzo
      Endif

      Return
      End
