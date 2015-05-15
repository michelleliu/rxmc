      Subroutine React_Forward(Irxn)
      Implicit None

      Integer Listreactants(10),Listproducts(10)

C     Change Reactants Into Products In The Reaction Ensemble

      Include 'commons.inc'

      Logical Lreact,Laccept,Doubledipping
      Integer Ipart,Irxn,ireactant,jreactant,Ireacttype
     $     ,ireactantmolecule,reactantcounter
     $     ,iproduct,Iprodtype,Ib
      Double Precision Xadd(4),Yadd(4),Zadd(4),Xdel(4),Ydel(4),Zdel(4)
     $     ,Randomnumber,Upotdel(4),Upotadd(4),Virdel(4),Viradd(4)
     $     ,Av1,Av2,Press,DeltaUrxn,Kequil
C      Double Precision Xi,Yi,Zi,Av1,Av2,Unew,Virnew,Randomnumber,Press

C     Equilibrium constant
      Kequil = 1.0

      If(Randomnumber().Lt.0.5d0) Then
         Ib = 1
      Else
         Ib = 2
      Endif

      If(Npart.Eq.0) Return
      If(Nrxn.Eq.0) Return

C     Choose a reaction
      Irxn = 1 + Int(Dble(Nrxn*Randomnumber()))

C     Check that there are enough reactants in Box Ib
      Lreact=.True.
      Do ireactant=1,Nreactants(Irxn)
         If(Npboxtype(Ib,Reactants(Irxn,ireactant)).lt.ReactantStoich(Irxn,
     $        ireactant)) Then
            Lreact=.False.
         Endif
      Enddo

      If(Lreact) Then
         reactantcounter = 0
         Do ireactant=1,Nreactants(Irxn)
            Ireacttype = Reactants(Irxn,ireactant)
            Do ireactantmolecule=1,Reactants(Irxn,ireactant)
               reactantcounter = reactantcounter + 1
C              Choose Nreactants(Irxn) reactants from Box Ib
               Doubledipping=.False.
               Ipart = 1 + Int(Dble(Npart)*Randomnumber())
               If(reactantcounter.gt.1) Then
                  Do jreactant=1,reactantcounter-1
C                    Check that we are not using the same particle twice
                     If(Listreactants(jreactant).eq.Ipart) Then
                        Doubledipping=.True.
                     Endif
                  Enddo
               Endif
               Do While (Types(Ipart).ne.Ireacttype.Or.Doubledipping
     $              .Or.Ibox(Ipart).ne.Ib)
                  Doubledipping=.False.
                  Ipart = 1 + Int(Dble(Npart)*Randomnumber())
                  If(reactantcounter.gt.1) Then
                     Do jreactant=1,reactantcounter-1
                        If(Listreactants(jreactant).eq.Ipart) Then
                           Doubledipping=.True.
                        Endif
                     Enddo
                  Endif
               Enddo
            Enddo
            Listreactants(reactantcounter)=Ipart
         Enddo

         DeltaUrxn = Dble(0)
         Do ireactant=1,Nreactants(Irxn)
            Ipart = Listreactants(ireactant)
            Xdel(ireactant) = Rx(Ipart)
            Ydel(ireactant) = Ry(Ipart)
            Zdel(ireactant) = Rz(Ipart)
            Call Epart(Ib,Virdel(ireactant),Upotdel(ireactant)
     $           ,Xdel(ireactant),Ydel(ireactant),Zdel(ireactant)
     $           ,Ipart,Types(Ipart))
            DeltaUrxn = DeltaUrxn - Upotdel(ireactant)
         Enddo

         Do iproduct=1,Nproducts(Irxn)
            Ipart = Listproducts(iproduct)
            Xadd(iproduct) = Randomnumber()*Box(Ib)
            Yadd(iproduct) = Randomnumber()*Box(Ib)
            Zadd(iproduct) = Randomnumber()*Box(Ib)
            Call Epart(Ib,Virdel(iproduct),Upotadd(iproduct)
     $           ,Xadd(iproduct),Yadd(iproduct),Zadd(iproduct)
     $           ,Ipart,Types(Ipart))
            DeltaUrxn = DeltaUrxn + Upotadd(iproduct)
         Enddo
C     Acceptance?
      !! FIXME
      Call Accept((Dble(Npboxtype(Ib,Ireacttype))*(Box(Ib)**3)/
     &     ((Box(Ib)**3)*Dble(Npboxtype(Ib,Ireacttype)+1)))*
     &     Dexp(-Beta*(DeltaUrxn)),Laccept)

C     Increment number of attempted moves
         Av2 = Av2 + 1.0d0

      Else
C     Increment number of attempted moves
         Av2 = Av2 + 1.0d0
      Endif
      Return
      End
