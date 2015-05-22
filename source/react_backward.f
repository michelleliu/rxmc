      Subroutine React_Backward(Av1,Av2)
      Implicit None

      Integer Listproducts(10)

C     Change Products Into Reactants In The Reaction Ensemble

      Include 'commons.inc'

      Logical Lreact,Laccept,Doubledipping
      Integer I,Ipart,Irxn,iproduct,jproduct,Iprodtype,Ireacttype
     $     ,istoich,reactantcounter,productcounter
     $     ,ireactant,jreactant,Ib,HoleCounter
     $     ,numproducts,numreactants,OldNpart !! number of molecules
     $     ,Mytype,Yourtype,selectioncount,J,Typeadd(10)
      Double Precision Av1,Av2,DeltaUrxn,DeltaVirrxn,Kequil
     $     ,AccNumer,AccDenom,Xdel,Ydel,Zdel,Xadd(10),Yadd(10),Zadd(10)
     $     ,Randomnumber,Upot,Vir,Dummy

      Double Precision Xi,Yi,Zi,Xj,Yj,Zj
C      Double Precision Xi,Yi,Zi,Av1,Av2,Unew,Virnew,Randomnumber

C     Equilibrium constant
      Kequil = 1.0d0
      AccNumer = 1.0d0
      AccDenom = 1.0d0

C     Initialize Listproducts array to zero
      Do I=1,10
         Listproducts(I)=0
         Xadd(I)=0.0d0
         Yadd(I)=0.0d0
         Zadd(I)=0.0d0
         Typeadd(I)=-1
      Enddo

      If(Randomnumber().Lt.0.5d0) Then
         Ib = 1
      Else
         Ib = 2
      Endif

      If(Npart.Eq.0) Return
      If(Nrxn.Eq.0) Return

C     Choose a reaction
      Irxn = 1 + Int(Dble(Nrxn*Randomnumber()))

C     Check that there are enough pr0ducts in Box Ib
      Lreact=.True.
      Do iproduct=1,Nproducts(Irxn)
         If(Npboxtype(Ib,Products(Irxn,iproduct)).lt.ProductStoich(Irxn,
     $        iproduct)) Then
            Lreact=.False.
         Endif
      Enddo

      If(Lreact) Then
         productcounter = 0
C        productcounter indexes Listproducts
C        Listproducts contains particle IDs of pr0duct molecules

         Do iproduct=1,Nproducts(Irxn)
            Iprodtype = Products(Irxn,iproduct)
C           Iprodtype is temporary variable which holds type of current pr0duct

            Do istoich=1,ProductStoich(Irxn,iproduct)
               productcounter = productcounter + 1
               AccNumer = AccNumer *
     $              ( Npboxtype(Ib,Products(Irxn,iproduct)) -
     $              istoich + 1 )
               Doubledipping=.False.
               Ipart = 1 + Int(Dble(Npart)*Randomnumber())
C              Choose pr0duct from Box Ib
               If(productcounter.gt.1) Then
C              Don't do this loop the first time
                  Do jproduct=1,productcounter-1
C                 jproduct indexes Listproducts in the small loop
                     If(Listproducts(jproduct).eq.Ipart) Then
                        Doubledipping=.True.
C                    Check that we are not using the same particle twice
                     Endif
                  Enddo
               Endif

               selectioncount = 0

               Do While (Types(Ipart).ne.Iprodtype.Or.Doubledipping
     $              .Or.Ibox(Ipart).ne.Ib)
                  selectioncount = selectioncount + 1
C                  If(selectioncount.Gt.10000) Then
CC                    Write trajectory and quit
C                    Write(6,*) Npboxtype(Ib,Iprodtype)
C    $                    ," alleged particles of type ",Iprodtype
C    $                    ," in box ",Ib
C                    Do J=1,Npart
C                       If(Ibox(J).Eq.2) Then
C                          Dummy = 4.0d0*(Box(1) + 2.0d0)
C                       Else
C                          Dummy = 0.0d0
C                       Endif
C                       Write(6,'(I2,3f15.5)') Types(J)
C    &                       ,4.0d0*Rx(J)+Dummy
C    &                       ,4.0d0*Ry(J),4.0d0*Rz(J)
C                    Enddo
C                    Stop "Error not enough particles!"
                     !If(Types(Ipart).ne.Iprodtype) Then
                     !   Write(6,*) 'Wrong type rxmc'
                     !Endif
                     !If(Ibox(Ipart).ne.Ib) Then
                     !   Write(6,*) 'Wrong box rxmc'
                     !Endif
C                  Endif
                  Doubledipping=.False.
                  !Write(6,*) 'Choosing particle again...'
                  Ipart = 1 + Int(Dble(Npart)*Randomnumber())
                  If(productcounter.gt.1) Then
                     Do jproduct=1,productcounter-1
                        If(Listproducts(jproduct).eq.Ipart) Then
                           Doubledipping=.True.
                           If(selectioncount.Gt.100) Then
                              !Write(6,*) 'Double dipping'
                           Endif
                        Endif
                     Enddo
                  Endif
               Enddo
C              Now pr0duct should meet all requirements
               Listproducts(productcounter)=Ipart
            Enddo
         Enddo
         numproducts=productcounter

C        Choose new coordinates for re4ctant moledules
         !! FIXME choose random coordinates for now
         !! FIXME pick smarter coordinates later

         reactantcounter = 0
         Do ireactant=1,Nreactants(Irxn)
            Do istoich=1,ReactantStoich(Irxn,ireactant)
               reactantcounter = reactantcounter + 1
               AccDenom = AccDenom *
     $              ( Npboxtype(Ib,Reactants(Irxn,ireactant)) -
     $              istoich + 1 + ReactantStoich(Irxn,ireactant))
               Xadd(reactantcounter) = Randomnumber()*Box(Ib)
               Yadd(reactantcounter) = Randomnumber()*Box(Ib)
               Zadd(reactantcounter) = Randomnumber()*Box(Ib)
               Typeadd(reactantcounter)=Reactants(Irxn,ireactant)
            Enddo
         Enddo
         numreactants=reactantcounter

C        Calculate DeltaUrxn for acceptance criteria:
C           - Calculate Upotdel of pr0ducts (epart)
C           - Sum Upotdel of pr0ducts
C           - Calculate pairwise energies (epair) and subtract from Upotdel
C           - Calculate Upotadd of re4ctants (epart)
C           - Sum Upotadd of re4ctants
C           - Calculate pairwise energies (epair) and add to Upotadd
C           - DeltaUrxn = Upotadd - Upotdel

         DeltaUrxn = 0.0d0
         DeltaVirrxn = 0.0d0

C        Calculate Delta U contribution of pr0ducts
         Do iproduct=1,numproducts
            Ipart = Listproducts(iproduct)
            Xdel = Rx(Ipart)
            Ydel = Ry(Ipart)
            Zdel = Rz(Ipart)
            Call Epart(Ib,Vir,Upot,Xdel,Ydel,Zdel
     $           ,Ipart,Types(Listproducts(iproduct)))
            DeltaUrxn = DeltaUrxn - Upot
            DeltaVirrxn = DeltaVirrxn - Vir
            Do jproduct=1,numproducts
               If(iproduct.ne.jproduct) Then
                  Xi = Rx(Listproducts(iproduct))
                  Yi = Ry(Listproducts(iproduct))
                  Zi = Rz(Listproducts(iproduct))
                  Xj = Rx(Listproducts(jproduct))
                  Yj = Ry(Listproducts(jproduct))
                  Zj = Rz(Listproducts(jproduct))
                  Mytype = Types(Listproducts(iproduct))
                  Yourtype = Types(Listproducts(jproduct))
                  Call Epair(Ib,Vir,Upot,Xi,Yi,Zi,Xj,Yj,Zj
     $                 ,Mytype,Yourtype)
                  DeltaUrxn = DeltaUrxn + Upot
                  DeltaVirrxn = DeltaVirrxn + Vir
               Endif
            Enddo
         Enddo

C        Calculate Delta U contribution of re4ctants
         Do ireactant=1,numreactants
            Ipart = Npart + 1
            Xi = Xadd(ireactant)
            Yi = Yadd(ireactant)
            Zi = Zadd(ireactant)
            Call Epart(Ib,Vir,Upot,Xi,Yi,Zi,Ipart
     $           ,Typeadd(reactantcounter))
            DeltaUrxn = DeltaUrxn + Upot
            DeltaVirrxn = DeltaVirrxn + Vir
            Do jreactant=1,numreactants
               If(jreactant.ne.ireactant) Then
                  Xi = Xadd(ireactant)
                  Yi = Yadd(ireactant)
                  Zi = Zadd(ireactant)
                  Xj = Xadd(jreactant)
                  Yj = Yadd(jreactant)
                  Zj = Zadd(jreactant)
                  Mytype = Typeadd(ireactant)
                  Yourtype = Typeadd(jreactant)
                  Call Epair(Ib,Vir,Upot,Xi,Yi,Zi,Xj,Yj,Zj
     $                 ,Mytype,Yourtype)
                  DeltaUrxn = DeltaUrxn - Upot
                  DeltaVirrxn = DeltaVirrxn - Vir
               Endif
            Enddo
            Do iproduct=1,numproducts
               Xi = Xadd(ireactant)
               Yi = Yadd(ireactant)
               Zi = Zadd(ireactant)
               Xj = Rx(Listproducts(iproduct))
               Yj = Ry(Listproducts(iproduct))
               Zj = Rz(Listproducts(iproduct))
               Mytype = Typeadd(ireactant)
               Yourtype = Types(Listproducts(iproduct))
               Call Epair(Ib,Vir,Upot,Xi,Yi,Zi,Xj,Yj,Zj
     $              ,Mytype,Yourtype)
                  DeltaUrxn = DeltaUrxn - Upot
                  DeltaVirrxn = DeltaVirrxn - Vir
            Enddo
         Enddo

         !Write(6,*) '  Checkpoint 1'
C        Acceptance?
         !Write(6,*) AccNumer,AccDenom,Kequil,Beta,DeltaUrxn
         If(DeltaUrxn.lt.((-1/Beta)*Log(AccDenom/(AccNumer*Kequil)))) Then
            Laccept=.True.
         Else
            Call Accept((AccNumer/AccDenom)*Kequil*Dexp(-Beta*DeltaUrxn)
     &           ,Laccept)
         Endif
         !Write(6,*) '  Checkpoint 2'

C        Accept Or Reject

         If(Laccept) Then

            !Write(6,*) 'Backward reaction accepted'
            Av1 = Av1 + 1.0d0

            !Write(6,*) 'reacting molecules'
            !Do iproduct=1,numproducts
               !Write(6,*) Listproducts(iproduct)
            !Enddo

C           Bookkeeping numbers of things
            OldNpart = Npart

            Do J=1,Npart
               If(Types(J).eq.0) Then
                  Write(6,*) "Type error before bookkeeping!"
               Endif
            Enddo

            Do iproduct=1,Nproducts(Irxn)
               Iprodtype=Products(Irxn,iproduct)
               Npboxtype(Ib,Iprodtype) = Npboxtype(Ib,Iprodtype) -
     &              ProductStoich(Irxn,iproduct)
               Npbox(Ib) = Npbox(Ib) - ProductStoich(Irxn,iproduct)
               Npart = Npart - ProductStoich(Irxn,iproduct)
               Nparts(Iprodtype) = Nparts(Iprodtype) -
     &              ProductStoich(Irxn,iproduct)

            Enddo
            Do ireactant=1,Nreactants(Irxn)
               Ireacttype=Reactants(Irxn,ireactant)
               Npboxtype(Ib,Ireacttype) = Npboxtype(Ib,Ireacttype) +
     &              ReactantStoich(Irxn,ireactant)
               Npbox(Ib) = Npbox(Ib) + ReactantStoich(Irxn,ireactant)
               Npart = Npart + ReactantStoich(Irxn,ireactant)
               Nparts(Ireacttype) = Nparts(Ireacttype) +
     &              ProductStoich(Irxn,ireactant)
            Enddo
            Etotal(Ib) = Etotal(Ib) - DeltaUrxn
            Vtotal(Ib) = Vtotal(Ib) - DeltaVirrxn

C           All the bookkeeping
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ stuff below here is sketchy

C           Establish re4ctant indices in place of pr0ducts
            !Write(6,*) 'Reactant indices'
            If(numreactants.gt.numproducts) Then
               Do ireactant=1,numreactants
                  If(ireactant.le.numproducts) Then
                     Rx(Listproducts(ireactant)) = Xadd(ireactant)
                     Ry(Listproducts(ireactant)) = Yadd(ireactant)
                     Rz(Listproducts(ireactant)) = Zadd(ireactant)
                     Ibox(Listproducts(ireactant)) = Ib
                     Types(Listproducts(ireactant)) = Typeadd(ireactant)
                     !Write(6,*) Listproducts(ireactant),' filled in'
                  Else
                     Rx(OldNpart+ireactant-numproducts) = Xadd(ireactant)
                     Ry(OldNpart+ireactant-numproducts) = Yadd(ireactant)
                     Rz(OldNpart+ireactant-numproducts) = Zadd(ireactant)
                     Ibox(OldNpart+ireactant-numproducts) = Ib
                     Types(OldNpart+ireactant-numproducts) = Typeadd(ireactant)
                     !Write(6,*) OldNpart+ireactant-numproducts,' added'
                  Endif
               Enddo
            Else
               Do ireactant=1,numreactants
                  Rx(Listproducts(ireactant)) = Xadd(ireactant)
                  Ry(Listproducts(ireactant)) = Yadd(ireactant)
                  Rz(Listproducts(ireactant)) = Zadd(ireactant)
                  Ibox(Listproducts(ireactant)) = Ib
                  Types(Listproducts(ireactant)) = Typeadd(ireactant)
                  !Write(6,*) Listproducts(ireactant),' filled in'
                  Listproducts(ireactant)=-1
               Enddo
C           Fill in holes left by removing pr0duct particles
C           Possibly super ratchet
               HoleCounter = -1
            !If(numreactants.le.numproducts) Then
               !Do iproduct=(numreactants+1),numproducts
               Do iproduct=(numreactants+1),numproducts
                  HoleCounter = HoleCounter + 1
C                 Move to next particle if current index is a pr0duct
                  Do While ( ANY( Listproducts==(OldNpart-HoleCounter) ))
                     HoleCounter = HoleCounter + 1
                  Enddo
C                 Fill in holes in pr0duct list
                  If(Listproducts(iproduct).le.OldNpart) Then
                     Rx(Listproducts(iproduct)) = Rx(OldNpart - HoleCounter)
                     Ry(Listproducts(iproduct)) = Ry(OldNpart - HoleCounter)
                     Rz(Listproducts(iproduct)) = Rz(OldNpart - HoleCounter)
                     Ibox(Listproducts(iproduct)) = Ibox(OldNpart - HoleCounter)
                     Types(Listproducts(iproduct)) = Types(OldNpart - HoleCounter)
C                    Write(6,*) OldNpart-Holecounter,' in place of
C    &                    ',Listproducts(iproduct)
                  Endif
               Enddo
            Endif
C           Check types
            Do J=1,Npart
               If(Types(J).eq.0) Then
                  Write(6,*) "React backward type error!"
               Endif
            Enddo
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end of sketchy

C        Increment number of attempted moves
         Av2 = Av2 + 1.0d0

         Endif
      Else
         !Write(6,*) '  Not reacting'
C        Increment number of attempted moves
         Av2 = Av2 + 1.0d0
      Endif
      !Write(6,*) '  Done react backward, Av2 = ',Av2
      Return
      End
