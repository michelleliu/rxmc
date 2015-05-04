      Subroutine Swap_Add(Av1,Av2,Press)
      Implicit None

C     Add Particle In The Grand-Canonical Ensemble

      Include 'commons.inc'

      Logical Laccept
      Integer Iaddtype
      Double Precision Xi,Yi,Zi,Av1,Av2,Unew,Virnew,Randomnumber,Press

C     Uncomment for debugging:
C      write(6,*) 'Npart                 : ',Npart
      If(Npart+1.Gt.Maxpart) Stop "Error Swapadd, Maxpart !!!"

      Av2 = Av2 + 1.0d0

      Xi = Randomnumber()*Box(1)
      Yi = Randomnumber()*Box(1)
      Zi = Randomnumber()*Box(1)

      Iaddtype = 1 + Int(Dble(Ntype)*Randomnumber())

      Call Epart(1,Virnew,Unew,Xi,Yi,Zi,0,Iaddtype)

      !!! FIXME this should be with partial pressure!!!
      Call Accept(Dexp(-Beta*Unew)*Beta*Press*MolFrac(Iaddtype)
     &     *(Box(1)**3)/Dble(1+Npart)
     &     ,Laccept)

C     Accept Or Reject

      If(Laccept) Then
         Av1 = Av1 + 1.0d0

         Etotal(1) = Etotal(1) + Unew
         Vtotal(1) = Vtotal(1) + Virnew

C        FIXME !!! section fixed but needs checking
         Npart = Npart + 1
         Nparts(Iaddtype) = Nparts(Iaddtype) + 1
         Types(Npart) = Iaddtype

         Rx(Npart)   = Xi
         Ry(Npart)   = Yi
         Rz(Npart)   = Zi
         Ibox(Npart) = 1

         Npbox(1) = Npbox(1) + 1
         Npboxtype(1,Iaddtype) = Npboxtype(1,Iaddtype) + 1
      Endif

      Return
      End
