      Subroutine Sample(Ic,Iensemble)
      Implicit None

C     Sample Simple Averages

      Include 'commons.inc'

      Integer I,J,Iensemble,Ic,Nb
      Double Precision Av1(5,2),Av2,Avmu1(2),Avmu2(2),Randomnumber,Xi,Yi
     $     ,Zi,Up,Dummy

      Integer Maxbin
      Parameter (Maxbin = 2000)

      Double Precision Delta,Idelta,Avrho1(Maxbin,2),Avrho2

      Parameter (Delta  = 0.001d0)
      Parameter (Idelta = 1.0d0/Delta)

      Save Av1,Av2,Avmu1,Avmu2,Avrho1,Avrho2

      If(Ic.Eq.1) Then

C     Set Counters To Zero

         Av2    = 0.0d0
         Avrho2 = 0.0d0

         Do I=1,5
            Do J=1,2
               Av1(I,J) = 0.0d0
               Avmu1(J) = 0.0d0
               Avmu2(J) = 0.0d0
            Enddo
         Enddo

         Do I=1,Maxbin
            Do J=1,2
               Avrho1(I,J) = 0.0d0
            Enddo
         Enddo

      Elseif(Ic.Eq.2) Then

C     Sample Ensemble Averages

         Av2 = Av2 + 1.0d0

         Do J=1,2
            Av1(1,J) = Av1(1,J) + Etotal(J)
C            Write(6,*) 'Energy',Av1(1,J)
            If(Iensemble.Ne.3) Av1(2,J) = Av1(2,J) +
     &           Dble(Npbox(J))/(Beta*(Box(J)**3)) +
     &           Vtotal(J)/(3.0d0*(Box(J)**3))
            Av1(3,J) = Av1(3,J) + Dble(Npbox(J))/(Box(J)**3)
            Av1(4,J) = Av1(4,J) + Box(J)**3
            Av1(5,J) = Av1(5,J) + Dble(Npbox(J))
         Enddo

      Elseif(Ic.Eq.3) Then

C     Sample Chemical Potential

         If(Iensemble.Ne.3) Then

            If(Iensemble.Eq.4) Then
               Nb = 2
            Else
               Nb = 1
            Endif

            Do J=1,Nb
               Do I=1,Npart
                  Xi = Randomnumber()*Box(J)
                  Yi = Randomnumber()*Box(J)
                  Zi = Randomnumber()*Box(J)

C                 FIXME !!!!
                  Call Epart(J,Dummy,Up,Xi,Yi,Zi,0,1)

                  Avmu1(J) = Avmu1(J) +
     &                 ((Box(J)**3)/(Dble(Npbox(J)+1)))*Dexp(-Beta*Up)

                  Avmu2(J) = Avmu2(J) + 1.0d0
               Enddo
            Enddo
         Endif

C     Sample Density Distribution

         If(Iensemble.Eq.4.Or.Iensemble.Eq.2) Then
            Avrho2 = Avrho2 + 1.0d0

            Do J=1,Nbox
               I = 1 + Int(Idelta*Dble(Npbox(J))/(Box(J)**3))

               If(I.Le.Maxbin) Avrho1(I,J) = Avrho1(I,J) + 1.0d0
            Enddo
         Endif

      Else

C     Write Averages

         Write(6,*)
         Write(6,*)
         Write(6,*) 'Averages Box 1'
         Write(6,*)
         Write(6,*) '<E>                  : ',Av1(1,1)/Av2
         If(Iensemble.Ne.3) Write(6,*) '<P>                  : '
     &        ,Av1(2,1)/Av2
         Write(6,*) '<Rho>                : ',Av1(3,1)/Av2
         Write(6,*) '<V>                  : ',Av1(4,1)/Av2
         Write(6,*) '<N>                  : ',Av1(5,1)/Av2
         If(Iensemble.Ne.3) Then
            Write(6,*) '<Mu> (Widom)         : '
     &           ,-Dlog(Avmu1(1)/Avmu2(1))/Beta
            Write(6,*) '<Mu>_Excess (Widom)  : '
     &           ,-Dlog(Avmu1(1)/Avmu2(1))/Beta - Dlog(Av1(3,1)/Av2)/Beta
            Write(6,*) '<Mu>_IG              : ',Dlog(Av1(3,1)/Av2)/Beta
         Endif
         Write(6,*)
         Write(6,*)

         If(Iensemble.Eq.4) Then
            Write(6,*) 'Averages Box 2'
            Write(6,*)
            Write(6,*) '<E>                  : ',Av1(1,2)/Av2
            Write(6,*) '<P>                  : ',Av1(2,2)/Av2
            Write(6,*) '<Rho>                : ',Av1(3,2)/Av2
            Write(6,*) '<V>                  : ',Av1(4,2)/Av2
            Write(6,*) '<N>                  : ',Av1(5,2)/Av2
            Write(6,*) '<Mu> (Widom)         : '
     &           ,-Dlog(Avmu1(2)/Avmu2(2))/Beta
            Write(6,*) '<Mu>_Excess (Widom)  : '
     &           ,-Dlog(Avmu1(2)/Avmu2(2))/Beta - Dlog(Av1(3,2)/Av2)/Beta
            Write(6,*) '<Mu>_IG              : ',Dlog(Av1(3,2)/Av2)/Beta
            Write(6,*)
            Write(6,*)
         Endif

         If(Iensemble.Eq.4.Or.Iensemble.Eq.2) Then
            Open(21,File="Rho",Status="Unknown")
            Do I=1,Maxbin
               Write(21,'(3e20.10)') (Dble(I)-0.5d0)*Delta
     &              ,Avrho1(I,1)/Avrho2,Avrho1(I,2)/Avrho2
            Enddo
            Close(21)
         Endif
      Endif

      Return
      End
