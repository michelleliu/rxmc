      Subroutine Init
      Implicit None

      Include 'commons.inc'

C     Generate An Initial Configuration

      Logical Laccept
      Integer I,J,Ipart,Ib,itype,Itmp
      Double Precision Fac,Randomnumber,Xi,Yi,Zi,Rxtrial,Rytrial,Rztrial
     $     ,Unew,Uold,Virnew,Virold

      Fac = Box(1)**3/(Box(1)**3 + Box(2)**3)

      Npbox(1) = 0
      Npbox(2) = 0

      Do itype=1,Ntype
         Npboxtype(1,itype)=0
         Npboxtype(2,itype)=0
      Enddo

C     Generate Random Coordinates And
C     Place Particles In A Certain Box

      Itmp = 0
      Do itype=1,Ntype
         Do I=1,Nparts(itype)
            Itmp = Itmp+1

            If (Verbose) Then
               Write(6,*) itype,I,Itmp
            Endif

            If(Nbox.Eq.1) Then
               Ibox(Itmp) = 1
            Else
               If(Randomnumber().Lt.Fac) Then
                  Ibox(Itmp) = 1
               Else
                  Ibox(Itmp) = 2
               Endif
            Endif

            Types(Itmp) = itype

            Npbox(Ibox(Itmp)) = Npbox(Ibox(Itmp)) + 1
            Npboxtype(Ibox(Itmp),itype) = Npboxtype(Ibox(Itmp),itype) +
     &           1

            Rx(Itmp) = Box(Ibox(Itmp))*Randomnumber()
            Ry(Itmp) = Box(Ibox(Itmp))*Randomnumber()
            Rz(Itmp) = Box(Ibox(Itmp))*Randomnumber()
         Enddo
      Enddo

C     Monte Carlo Displacements To Remove Initial Overlaps

      Do J=1,Npart*50
         Ipart = 1 + Int(Dble(Npart)*Randomnumber())
         Ib    = Ibox(Ipart)

         Rxtrial = Rx(Ipart) + Randomnumber() - 0.5d0
         Rytrial = Ry(Ipart) + Randomnumber() - 0.5d0
         Rztrial = Rz(Ipart) + Randomnumber() - 0.5d0

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

         Xi = Rx(Ipart)
         Yi = Ry(Ipart)
         Zi = Rz(Ipart)

         Call Epart(Ib,Virold,Uold,Xi,Yi,Zi,Ipart,Types(Ipart))
         Call Epart(Ib,Virnew,Unew,Rxtrial,Rytrial,Rztrial
     &        ,Ipart,Types(Ipart))

         Write(6,*) Unew,Uold
         Call Accept(Dexp(-0.5d0*(Unew-Uold)),Laccept)

         If(Laccept) Then
            Rx(Ipart) = Rxtrial
            Ry(Ipart) = Rytrial
            Rz(Ipart) = Rztrial
         Endif
      Enddo

      Return
      End
