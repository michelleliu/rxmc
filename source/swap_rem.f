      Subroutine Swap_Rem(Av1,Av2,Press)
      Implicit None

C     Remove A Particle In The Grand Canonical Ensemble

      Include 'commons.inc'

      Logical Laccept
      Integer Ipart,Iremtype
      Double Precision Xi,Yi,Zi,Randomnumber,Uold,Virold,Av1,Av2
     $     ,Press
C     Av1 is number of accepted moves
C     Av2 is number of attempted moves

      If(Npart.Eq.0) Return

      Iremtype = 1 + Int(Dble(Ntype)*Randomnumber())
C     Test that this type exists
      If(Nparts(Iremtype).gt.0) Then
C        Uncomment for debugging:
C         write(6,*) 'Nparts(type)           : ',Nparts(Iremtype)
         Ipart = 1 + Int(Dble(Npart)*Randomnumber())
         Do While (Types(Ipart).ne.Iremtype)
            Ipart = 1 + Int(Dble(Npart)*Randomnumber())
         Enddo

         Xi = Rx(Ipart)
         Yi = Ry(Ipart)
         Zi = Rz(Ipart)

         Av2 = Av2 + 1.0d0

         Call Epart(1,Virold,Uold,Xi,Yi,Zi,Ipart,Types(Ipart))

         Call Accept(Dble(Npart)/((Box(1)**3)*Press*MolFrac(Types(Ipart))
     &        *Beta*Dexp(-Beta*Uold))
     &        ,Laccept)

C        Accept Or Reject

         If(Laccept) Then
            Av1 = Av1 + 1.0d0

            Etotal(1) = Etotal(1) - Uold
            Vtotal(1) = Vtotal(1) - Virold

            Npart    = Npart    - 1
            Npbox(1) = Npbox(1) - 1
            Npboxtype(1,Types(Ipart)) = Npboxtype(1,Types(Ipart)) - 1
            Nparts(Types(Ipart)) = Nparts(Types(Ipart)) - 1

C           Fill in hole left by removing particle
            Rx(Ipart)   = Rx(Npart+1)
            Ry(Ipart)   = Ry(Npart+1)
            Rz(Ipart)   = Rz(Npart+1)
            Ibox(Ipart) = Ibox(Npart+1)
            Types(Ipart)= Types(Npart+1)

         Endif
      Else
C        Increment number of attempted moves
         Av2 = Av2 + 1.0d0
      Endif

      Return
      End
