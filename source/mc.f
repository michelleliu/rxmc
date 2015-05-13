      Program Mc
      Implicit None

C     Monte Carlo In The Npt, Nvt, Muvt, And Gibbs Ensemble

      Include 'commons.inc'

      Logical Linit
      Integer Ncycle,Ninit,Nmove,Iensemble,I,J,Icycle,Icycle2,Dd(2)
     $     ,itype,jtype,icounter,irxn,jrxn
      Double Precision Pdisp,Pswap,Pvol,Prxn,Temp,Deltax,Deltav,V,E,Rm
     $     ,Randomnumber,Dummy,Avv1,Avv2,Avs1,Avs2,Avd1,Avd2,Press
     $     ,SumFrac,Avr1,Avr2

C     Set Verbose to .True. for verbose debugging statements
      Verbose = .False.

      Avv1 = 0.0d0
      Avv2 = 0.0d0
      Avs1 = 0.0d0
      Avs2 = 0.0d0
      Avd1 = 0.0d0
      Avd2 = 0.0d0
      Avr1 = 0.0d0
      Avr2 = 0.0d0

      Open(21,File="Input",Status="Unknown")
      Read(21,*)
      Read(21,*) Ncycle,Ninit,Linit,Temp,Iensemble
      Read(21,*)
      Read(21,*) Pdisp,Pswap,Pvol,Prxn
      Read(21,*)
C     Read number of types and reactions
      Read(21,*) Nrxn,Ntype
      Read(21,*)
      Read(21,*) Deltax,Deltav,Press
      Read(21,*)
      Read(21,*) Box(1),Box(2)
C     Loop over types
      Read(21,*)
      Do itype=1,Ntype
         Read(21,*) Nparts(itype),MolFrac(itype)
      Enddo
      SumFrac = 0
      Do itype=1,Ntype
         SumFrac = SumFrac + MolFrac(itype)
      Enddo
      Do itype=1,Ntype
         MolFrac(itype) = MolFrac(itype)/SumFrac
      Enddo
      Read(21,*)
      Do itype=1,Ntype
         Do jtype=itype,Ntype
            Read(21,*) Eps(itype,jtype),Sig(itype,jtype)
            Eps(jtype,itype)=Eps(itype,jtype)
            Sig(jtype,itype)=Sig(itype,jtype)
         Enddo
      Enddo

      Write(6,*) "Particles"
      Do itype=1,Ntype
         Write(6,*) itype, Nparts(itype), MolFrac(itype)
      Enddo
      Write(6,*)
      Write(6,*) "LJ Parameters"
      Do itype=1,Ntype
         Do jtype=1,Ntype
            Write(6,*) itype,jtype,Eps(itype,jtype),Sig(itype,jtype)
         Enddo
      Enddo
      Write(6,*)
C     Read in reactions
      Read(21,*)
      Do irxn=1,Nrxn
         Read(21,*) Nreactants(irxn),Nproducts(irxn)

      Enddo
      Read(21,*)
      Do irxn=1,Nrxn
         Read(21,*)
         Do jrxn=1,Nreactants(irxn)
            Read(21,*) Reactants(irxn,jrxn),ReactantStoich(irxn,jrxn)
         Enddo
         Do jrxn=1,Nproducts(irxn)
            Read(21,*) Products(irxn,jrxn),ProductStoich(irxn,jrxn)
         Enddo
      Enddo
      Close(21)

      Npart = 0
      Do itype=1,Ntype
         Npart = Npart+Nparts(itype)
      Enddo

      If(Temp.Le.0.0d0)   Stop "Error Temperature !!!"
      If(Pdisp.Le.0.0d0)  Pdisp = 0.0d0
      If(Pvol.Le.0.0d0)   Pvol  = 0.0d0
      If(Pswap.Le.0.0d0)  Pswap = 0.0d0
      If(Prxn.Le.0.0d0)  Prxn = 0.0d0
      If(Deltax.Lt.0.0d0) Stop "Error Deltax !!"
      If(Deltav.Lt.0.0d0) Stop "Error Deltav !!"
      If(Ncycle.Lt.100)   Stop "Minimal 100 Cycles !!!"

      Beta       = 1.0d0/Temp
      Npbox(1)   = 0
      Npbox(2)   = 0

      Do itype=1,Ntype
         Npboxtype(1,itype)=0
         Npboxtype(2,itype)=0
      Enddo

      Rcutsq = (2.5d0)**2
      Do itype=1,Ntype
         Do jtype=1,Ntype
            Ecut(itype,jtype)   = (Rcutsq**3)/(Eps(itype,jtype)**6)
            Ecut(itype,jtype)   = 1.0d0/Ecut(itype,jtype)
            Ecut(itype,jtype)   = 4.0d0*Sig(itype,jtype)*
     &           Ecut(itype,jtype)*(Ecut(itype,jtype)-1.0d0)
         Enddo
      Enddo
      write(6,*) 'ecut                 : ',ecut(1,1)

      Write(6,*) 'Ncycle               : ',Ncycle
      Write(6,*) 'Ninit                : ',Ninit
      Write(6,*) 'Linit                : ',Linit
      Write(6,*) 'Temp                 : ',Temp
      Write(6,*) 'Iens                 : ',Iensemble
      Write(6,*) 'Deltax               : ',Deltax
      Write(6,*) 'Deltav               : ',Deltav

C     Iensemble = 1    Nvt      Ensemble
C     Iensemble = 2    Npt      Ensemble
C     Iensemble = 3    Muvt     Ensemble
C     Iensemble = 4    Gibbs    Ensemble
C     Iensemble = 5    Reaction Ensemble

      If(Iensemble.Eq.1) Then
         Pswap = 0.0d0
         Pvol  = 0.0d0
         Nbox  = 1

         Write(6,*)
         Write(6,*) 'Nvt Ensemble'
         Write(6,*)

      Elseif(Iensemble.Eq.2) Then
         Pswap = 0.0d0
         Nbox  = 1

         If(Press.Le.0.0d0) Press = 0.0d0

         Write(6,*)
         Write(6,*) 'Npt Ensemble'
         Write(6,*)
         Write(6,*) 'Pressure             : ',Press
         Write(6,*)

      Elseif(Iensemble.Eq.3) Then
         Pvol  = 0.0d0
         Nbox  = 1

         If(Press.Lt.1.0d-100) Stop "Fugacity Too Low In Muvt Ensemble
     &        !!! "

         Write(6,*)
         Write(6,*) 'Muvt Ensemble'
         Write(6,*)
         Write(6,*) 'Mu                   : ',Dlog(Beta*Press)/Beta
         Write(6,*) 'Fugacity             : ',Press
         Write(6,*)

      Elseif(Iensemble.Eq.4) Then
         Nbox  = 2

         Write(6,*)
         Write(6,*) 'Gibbs Ensemble'
         Write(6,*)

      Elseif(Iensemble.Eq.5) Then
         Nbox  = 1
         Pvol  = 0.0d0

         Write(6,*)
         Write(6,*) 'Reaction Ensemble'
         Write(6,*) 'Number of Reactions  : ',Nrxn
         Write(6,*)
         Do irxn=1,Nrxn
            Write(6,*) 'Reaction ',irxn
            Write(6,*) '      Reactants'
            Do jrxn=1,Nreactants(irxn)
               Write(6,*) Reactants(irxn,jrxn)
            Enddo
            Write(6,*) '      Products'
            Do jrxn=1,Nproducts(irxn)
               Write(6,*) Products(irxn,jrxn)
            Enddo
         Enddo
      Else
         Stop "Error Ensemble !!!"
      Endif

      Dummy = Pdisp + Pswap + Pvol + Prxn
      Pdisp = Pdisp/Dummy
      Pswap = Pswap/Dummy
      Pvol  = Pvol/Dummy
      Prxn  = Prxn/Dummy

      Write(6,*) 'Pdisp                : ',Pdisp
      Write(6,*) 'Pswap                : ',Pswap
      Write(6,*) 'Pvol                 : ',Pvol
      Write(6,*) 'Prxn                 : ',Prxn

C     Add up probabilities so they are spaced correctly
      Pswap = Pswap + Pvol
      Prxn  = Prxn + Pswap

      If(Linit) Then
         Write(6,*)
         Write(6,*) 'Generate Initial Coordinates'
         Write(6,*)

         If(Npart.Lt.0.Or.Npart.Gt.Maxpart) Stop "Error Npart !!!"
C        Do itype=1,Ntype
C           If((Nparts(itype)).Lt.0.Or.(Nparts(itype)).Gt.Maxpart) Then
C              Stop "Error Nparts !!!"
C           Endif
C        Enddo
         If(Min(Box(1),Box(2)).Lt.5.0d0) Stop "Boxes Too Small !!!"

         Call Init
      Else
         Npbox(1)=0
         Npbox(2)=0
         Do I=1,Ntype
            Nparts(I)=0
            Do J=1,Nbox
               Npboxtype(J,I)=0
            Enddo
         Enddo
         Write(6,*)
         Write(6,*) 'Read Coordinates From Disk'
         Write(6,*)

         Open(21,File="Coordold",Status="Unknown")
         Read(21,*) Box(1),Box(2)
         Read(21,*) Npart
         Read(21,*) Ntype
         Read(21,*)

         If(Npart.Lt.0.Or.Npart.Gt.Maxpart) Stop "Error Npart !!!"
         If(Min(Box(1),Box(2)).Lt.5.0d0) Stop "Boxes Too Small !!!"

         Do I=1,Npart
            Read(21,*)Rx(I),Ry(I),Rz(I)
     &           ,Types(I),Ibox(I)

            If(Ibox(I).Eq.2.And.Iensemble.Ne.4) Stop "Particle
     &           In Wrong Box !!!"

            Npbox(Ibox(I)) = Npbox(Ibox(I)) + 1
            Nparts(Types(I)) = Nparts(Types(I)) + 1
            Npboxtype(Ibox(I),Types(I)) =
     &           Npboxtype(Ibox(I),Types(I)) + 1
         Enddo
      Endif

      Write(6,*) 'Box 1                : ',Box(1)
      If(Nbox.Eq.2) Write(6,*) 'Box 2                : ',Box(2)
      Npart = 0
      Do itype=1,Ntype
         Npart = Npart+Nparts(itype)
      Enddo
      Write(6,*) 'Total Npart          : ',Npart
      Do itype=1,Ntype
         Write(6,*) 'Npart type       : ',itype,Nparts(itype)
      Enddo
      Write(6,*) 'Npbox 1              : ',Npbox(1)
      Do itype=1,Ntype
         Write(6,*) 'Npbox 1, type    : ',itype,Npboxtype(1,itype)
      Enddo
      If(Nbox.Eq.2) Write(6,*) 'Npbox 2              : ',Npbox(2)
      Do itype=1,Ntype
         Write(6,*) 'Npbox 2, type    : ',itype,Npboxtype(2,itype)
      Enddo
      Write(6,*) '<Rho1>               : ',Dble(Npbox(1))/(Box(1)**3)
      If(Nbox.Eq.2)
     &     Write(6,*) '<Rho2>               : ',
     &     Dble(Npbox(2))/(Box(2)**3)

      If(Deltax.Gt.0.25d0*Min(Box(1),Box(2))) Stop "Deltax Too Large
     &     !!!"

      Etotal(1) = 0.0d0
      Etotal(2) = 0.0d0
      Vtotal(1) = 0.0d0
      Vtotal(2) = 0.0d0

      Do I=1,Nbox
         Write(6,*) 'Calling Etot'
         Call Etot(I,V,E)

         Vtotal(I) = V
         Etotal(I) = E
      Enddo

      Write(6,*)
      Write(6,*) 'Initial Energy Box 1 : ',Etotal(1)
      If(Nbox.Eq.2) Write(6,*) 'Initial Energy Box 2 : ',Etotal(2)
      Write(6,*) 'Initial Virial Box 1 : ',Vtotal(1)
      If(Nbox.Eq.2) Write(6,*) 'Initial Virial Box 2 : ',Vtotal(2)
      Write(6,*)

C     Start Of The Simulation

      Write(6,*)
      Write(6,*)
      Write(6,*) 'The Simulation Is Running.....'
      Write(6,*)
      Write(6,*)

      If(Iensemble.Eq.1) Call Sample_Radial(1)
      Call Sample(1,Iensemble)

      Open(22,File="Traject.xyz",Status="Unknown")
      Open(23,File="Results",Status="Unknown")

C     Starting simulation loop

      Do Icycle=1,Ncycle

         Nmove = Max(20,Npart)

         Do Icycle2=1,Nmove

C     Select Trial Move At Random

            Rm = Randomnumber()

            If(Rm.Lt.Pvol.And.(Iensemble.Eq.2.Or.Iensemble.Eq.4)) Then

C     Volume Change

               If(Iensemble.Eq.2) Then
                  Call Volume_Npt(Avv1,Avv2,Deltav,Press)
               Else
                  Call Volume_Gibbs(Avv1,Avv2,Deltav)
               Endif

            Elseif(Rm.Lt.Pswap.And.(Iensemble.Eq.3.Or.Iensemble.Eq.4)) Then

C     Particle Swap

               If(Iensemble.Eq.4) Then
                  Call Swap_Gibbs(Avs1,Avs2)
               Else
                  If(Randomnumber().Lt.0.5d0) Then
                     Call Swap_Add(Avs1,Avs2,Press)
                  Else
                     Call Swap_Rem(Avs1,Avs2,Press)
                  Endif
               Endif

            Elseif(Rm.Lt.Prxn.And.(Iensemble.Eq.5)) Then

C     Reaction Move

C              FIXME !!!!!!
C              If(Randomnumber().Lt.0.5d0) Then
C                 Call React_Forward(Avr1,Avr2,Press,stuff)
C              Else
C                 Call React_Reverse(Avr1,Avr2,Press,stuff)
C              Endif

            Else

C     Particle Displacement

               Call Move(Avd1,Avd2,Deltax)

            Endif

            If(Icycle.Gt.Ninit) Call Sample(2,Iensemble)

         Enddo
C     FLOATING POINT ERROR FIXME
C     Write results
         If(Nbox.Eq.2.AND.Ntype.Eq.1) Then
            Write(23,'(7e20.10)')
     &           Dble(Icycle),
     &           Dble(Npbox(1))/(Box(1)**3),
     &           Dble(Npbox(2))/(Box(2)**3),
     &           Etotal(1)/Max(0.5d0,Dble(Npbox(1))),
     &           Etotal(2)/Max(0.5d0,Dble(Npbox(2))),
     &           Dble(Npbox(1))/(Beta*(Box(1)**3))
     &           + Vtotal(1)/(3.0d0*(Box(1)**3)),
     &           Dble(Npbox(2))/(Beta*(Box(2)**3))
     &           + Vtotal(2)/(3.0d0*(Box(2)**3))
         Else If(Nbox.Eq.2.AND.Ntype.Eq.2) Then
            Write(23,'(11e20.10)')
     &           Dble(Icycle),
     &           Dble(Npboxtype(1,1))/(Box(1)**3),
     &           Dble(Npboxtype(1,2))/(Box(1)**3),
     &           Dble(Npboxtype(2,1))/(Box(2)**3),
     &           Dble(Npboxtype(2,2))/(Box(2)**3),
     &           Box(1),
     &           Box(2),
     &           Etotal(1)/Max(0.5d0,Dble(Npbox(1))),
     &           Etotal(2)/Max(0.5d0,Dble(Npbox(2))),
     &           Dble(Npbox(1))/(Beta*(Box(1)**3))
     &           + Vtotal(1)/(3.0d0*(Box(1)**3)),
     &           Dble(Npbox(2))/(Beta*(Box(2)**3))
     &           + Vtotal(2)/(3.0d0*(Box(2)**3))
         Else If(Nbox.Eq.2.AND.Ntype.Eq.3) Then
            Write(23,'(11e20.10)')
     &           Dble(Icycle),
     &           Dble(Npboxtype(1,1))/(Box(1)**3),
     &           Dble(Npboxtype(1,2))/(Box(1)**3),
     &           Dble(Npboxtype(1,3))/(Box(1)**3),
     &           Dble(Npboxtype(2,1))/(Box(2)**3),
     &           Dble(Npboxtype(2,3))/(Box(2)**3),
     &           Etotal(1)/Max(0.5d0,Dble(Npbox(1))),
     &           Etotal(2)/Max(0.5d0,Dble(Npbox(2))),
     &           Dble(Npbox(1))/(Beta*(Box(1)**3))
     &           + Vtotal(1)/(3.0d0*(Box(1)**3)),
     &           Dble(Npbox(2))/(Beta*(Box(2)**3))
     &           + Vtotal(2)/(3.0d0*(Box(2)**3))
         Else
            Write(23,'(4e20.10)')
     &           Dble(Icycle),
     &           Dble(Npbox(1))/(Box(1)**3),
     &           Etotal(1)/Max(0.5d0,Dble(Npbox(1))),
     &           Dble(Npbox(1))/(Beta*(Box(1)**3))
     &           + Vtotal(1)/(3.0d0*(Box(1)**3))
         Endif

         If(Icycle.Gt.Ninit) Then
            If(Iensemble.Eq.1) Call Sample_Radial(2)
            Call Sample(3,Iensemble)
         Endif

         If(Mod(Icycle,Ncycle/100).Eq.0) Then
            Write(22,*) Npart
            Enddo
            Write(22,*)

C     Write trajectory
            Do J=1,Npart
               If(Ibox(J).Eq.2) Then
                  Dummy = 4.0d0*(Box(1) + 2.0d0)
               Else
                  Dummy = 0.0d0
               Endif
C              Write(22,'(A,3f15.5)') 'Ar  '
               Write(22,'(I2,3f15.5)') Types(J)
     &              ,4.0d0*Rx(J)+Dummy
     &              ,4.0d0*Ry(J),4.0d0*Rz(J)
            Enddo

         Endif
      Enddo

C     End of simulation loop
C     Print closing statements

      Close(22)
      Close(23)

      If(Iensemble.Eq.1) Call Sample_Radial(3)
      Call Sample(4,Iensemble)

      Write(6,*) 'Frac. Acc. Displ.    : ',Avd1/Max(0.5d0,Avd2)
      Write(6,*) 'Frac. Acc. Swap      : ',Avs1/Max(0.5d0,Avs2)
      Write(6,*) 'Frac. Acc. Volume    : ',Avv1/Max(0.5d0,Avv2)
      Write(6,*) 'Frac. Acc. Reaction  : ',Avr1/Max(0.5d0,Avr2)

      Open(21,File="Coordnew",Status="Unknown")
      Write(21,*) Box(1),Box(2)
      Write(21,*) Npart

      Do I=1,Npart
         Write(21,'(3e20.10,i5)') Rx(I),Ry(I),Rz(I),Ibox(I)
      Enddo
      Close(21)

C     Check Energy Calculation

      Do I=1,Nbox
         Call Etot(I,V,E)

         Write(6,*)
         Write(6,*) 'Box                  : ',I
         Write(6,*) 'Particles            : ',Npbox(I)
         Do itype=1,Ntype
            Write(6,*) 'Particles type         : ',itype,Npboxtype(I,itype)
         Enddo
         Write(6,*) 'Box                  : ',Box(I)
         Write(6,*) 'Volume               : ',Box(I)**3
         Write(6,*) 'Rho                  : ',Dble(Npbox(I))/(Box(I)**3)
         Write(6,*) 'Energy               : ',E
         Write(6,*) 'Energy Sim.          : ',Etotal(I)
         Write(6,*) 'Diff                 : ',Dabs(Etotal(I)-E)
         Write(6,*) 'Virial               : ',V
         Write(6,*) 'Virial Sim.          : ',Vtotal(I)
         Write(6,*) 'Diff                 : ',Dabs(Vtotal(I)-V)
         Write(6,*)
      Enddo

      Dd(1) = 0
      Dd(2) = 0

      Do I=1,Npart
         Dd(Ibox(I)) = Dd(Ibox(I)) + 1
      Enddo

      If(Dd(1).Ne.Npbox(1).Or.Dd(2).Ne.Npbox(2)) Then
         Write(6,*) 'Error !!!!'
         Write(6,*) Npart
         Write(6,*) Dd(1),Npbox(1)
         Write(6,*) Dd(2),Npbox(2)
         Stop
      Endif

      Stop
      End
