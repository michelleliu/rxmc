      Subroutine Swap_Gibbs(Av1,Av2)
      Implicit None

C     Particle Swap In The Gibbs Ensemble

      Include 'commons.inc'

      Logical Laccept
      Integer Iadd,Idel,Ipart,Iswaptype
      Double Precision Xi,Yi,Zi,Randomnumber,Av1,Av2,Upotadd,Upotdel
     $     ,Viradd,Virdel

C     Select At Random Which Box To Add / Remove

      If(Randomnumber().Lt.0.5d0) Then
         Iadd = 1
         Idel = 2
      Else
         Iadd = 2
         Idel = 1
      Endif

      If(Npbox(Idel).Eq.0) Return

C     Select type of particle to swap
      Iswaptype = 1 + Int(Dble(Ntype)*Randomnumber())
      If(Npboxtype(Idel,Iswaptype).Eq.0) Return

C     Select Particle To Be Removed Random From Box Idel
C     Make sure particle is of type Iswaptype
 1    Continue
      Ipart = 1 + Int(Dble(Npart)*Randomnumber())
      If((Ibox(Ipart).Ne.Idel).OR.(Types(Ipart).Ne.Iswaptype)) Then
         Goto 1
      Endif

      Call Epart(Idel,Virdel,Upotdel,Rx(Ipart),Ry(Ipart),Rz(Ipart),Ipart,Types(Ipart))

      Xi = Randomnumber()*Box(Iadd)
      Yi = Randomnumber()*Box(Iadd)
      Zi = Randomnumber()*Box(Iadd)

      Call Epart(Iadd,Viradd,Upotadd,Xi,Yi,Zi,0,Iswaptype)

      Call Accept((Dble(Npboxtype(Idel,Iswaptype))*(Box(Iadd)**3)/
     &     ((Box(Idel)**3)*Dble(Npboxtype(Iadd,Iswaptype)+1)))*
     &     Dexp(-Beta*(Upotadd-Upotdel)),Laccept)

      Av2 = Av2 + 1.0d0

C     Accept Or Reject

      If(Laccept) Then
         Npboxtype(Iadd,Iswaptype) = Npboxtype(Iadd,Iswaptype) + 1
         Npboxtype(Idel,Iswaptype) = Npboxtype(Idel,Iswaptype) - 1
         Npbox(Iadd) = Npbox(Iadd) + 1
         Npbox(Idel) = Npbox(Idel) - 1

         Types(Ipart) = Iswaptype
         Ibox(Ipart) = Iadd

         Av1 = Av1 + 1.0d0
         Etotal(Iadd) = Etotal(Iadd) + Upotadd
         Etotal(Idel) = Etotal(Idel) - Upotdel

         Vtotal(Iadd) = Vtotal(Iadd) + Viradd
         Vtotal(Idel) = Vtotal(Idel) - Virdel

      Endif

      Return
      End
