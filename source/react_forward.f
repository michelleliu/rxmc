      Subroutine React_Forward(Av1,Av2)
      Implicit None

      Integer Listreactants(10)

C     Change Reactants Into Products In The Reaction Ensemble

      Include 'commons.inc'

      Logical Lreact,Laccept,Doubledipping
      Integer I,Ipart,Irxn,ireactant,jreactant,Ireacttype,Iprodtype
     $     ,istoich,productcounter,reactantcounter
     $     ,iproduct,jproduct,Ib,HoleCounter
     $     ,numreactants,numproducts,OldNpart
     $     ,Mytype,Yourtype,selectioncount,Typeadd(10)
      Double Precision Av1,Av2,DeltaUrxn,DeltaVirrxn,Kequil
     $     ,AccNumer,AccDenom,Xdel,Ydel,Zdel,Xadd(10),Yadd(10),Zadd(10)
     $     ,Randomnumber,Upot(5),Vir

      Double Precision Xi,Yi,Zi,Xj,Yj,Zj

C     Equilibrium constant
      Kequil = 1.0d0
      AccNumer = 1.0d0
      AccDenom = 1.0d0

C     Initialize Listreactants array to zero
      Do I=1,10
         Listreactants(I)=0
         Xadd(I)=0.0d0
         Yadd(I)=0.0d0
         Zadd(I)=0.0d0
         Typeadd(I)=-2
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
C        reactantcounter indexes Listreactants
C        Listreactants contains particle IDs of reactant molecules

         Do ireactant=1,Nreactants(Irxn)
            Ireacttype = Reactants(Irxn,ireactant)
C           Ireacttype is temporary variable which holds type of current reactant

            Do istoich=1,ReactantStoich(Irxn,ireactant)
               reactantcounter = reactantcounter + 1
               AccNumer = AccNumer *
     $              ( Npboxtype(Ib,Reactants(Irxn,ireactant)) -
     $              istoich + 1 )
               Doubledipping=.False.
               Ipart = 1 + Int(Dble(Npart)*Randomnumber())
C              Choose reactant from Box Ib
               If(reactantcounter.gt.1) Then
C              Don't do this loop the first time
                  Do jreactant=1,reactantcounter-1
C                 jreact indexes Listreactants in the small loop
                     If(Listreactants(jreactant).eq.Ipart) Then
                        Doubledipping=.True.
C                    Check that we are not using the same particle twice
                     Endif
                  Enddo
               Endif

               selectioncount = 0

               Do While (Types(Ipart).ne.Ireacttype.Or.Doubledipping
     $              .Or.Ibox(Ipart).ne.Ib)
                  selectioncount = selectioncount + 1
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
C              Now reactant should meet all requirements
               Listreactants(reactantcounter)=Ipart
            Enddo
         Enddo
         numreactants=reactantcounter

C        Choose new coordinates for product moledules
         !! FIXME choose random coordinates for now
         !! FIXME pick smarter coordinates later

         productcounter = 0
         Do iproduct=1,Nproducts(Irxn)
            Do istoich=1,ProductStoich(Irxn,iproduct)
               productcounter = productcounter + 1
               AccDenom = AccDenom *
     $              ( Npboxtype(Ib,Products(Irxn,iproduct)) -
     $              istoich + 1 + ProductStoich(Irxn,iproduct))
               Xadd(productcounter) = Randomnumber()*Box(Ib)
               Yadd(productcounter) = Randomnumber()*Box(Ib)
               Zadd(productcounter) = Randomnumber()*Box(Ib)
               Typeadd(productcounter)=Products(Irxn,iproduct)
            Enddo
         Enddo
         numproducts=productcounter

C        Calculate DeltaUrxn for acceptance criteria:
C           - Calculate Upotdel of reactants (epart)
C           - Sum Upotdel of reactants
C           - Calculate pairwise energies (epair) and subtract from Upotdel
C           - Calculate Upotadd of products (epart)
C           - Sum Upotadd of products
C           - Calculate pairwise energies (epair) and add to Upotadd
C           - DeltaUrxn = Upotadd - Upotdel

         DeltaUrxn = 0.0d0
         DeltaVirrxn = 0.0d0

C        Calculate Delta U contribution of reactants
         Do ireactant=1,numreactants
            Ipart = Listreactants(ireactant)
            Xdel = Rx(Ipart)
            Ydel = Ry(Ipart)
            Zdel = Rz(Ipart)
            Call Epart(Ib,Vir,Upot(1),Xdel,Ydel,Zdel
     $           ,Ipart,Types(Listreactants(ireactant)))
            DeltaUrxn = DeltaUrxn - Upot(1)
            DeltaVirrxn = DeltaVirrxn - Vir
            Do jreactant=1,numreactants
               If(ireactant.ne.jreactant) Then
                  Xi = Rx(Listreactants(ireactant))
                  Yi = Ry(Listreactants(ireactant))
                  Zi = Rz(Listreactants(ireactant))
                  Xj = Rx(Listreactants(jreactant))
                  Yj = Ry(Listreactants(jreactant))
                  Zj = Rz(Listreactants(jreactant))
                  Mytype = Types(Listreactants(ireactant))
                  Yourtype = Types(Listreactants(jreactant))
                  Call Epair(Ib,Vir,Upot(2),Xi,Yi,Zi,Xj,Yj,Zj
     $                 ,Mytype,Yourtype)
                  DeltaUrxn = DeltaUrxn + Upot(2)
                  DeltaVirrxn = DeltaVirrxn + Vir
               Endif
            Enddo
         Enddo

C        Calculate Delta U contribution of products
         Do iproduct=1,numproducts
            Ipart = Npart + 1
            Xi = Xadd(iproduct)
            Yi = Yadd(iproduct)
            Zi = Zadd(iproduct)
            Call Epart(Ib,Vir,Upot(3),Xi,Yi,Zi,Ipart
     $           ,Typeadd(productcounter))
            DeltaUrxn = DeltaUrxn + Upot(3)
            DeltaVirrxn = DeltaVirrxn + Vir
            Do jproduct=1,numproducts
               If(jproduct.ne.iproduct) Then
                  Xi = Xadd(iproduct)
                  Yi = Yadd(iproduct)
                  Zi = Zadd(iproduct)
                  Xj = Xadd(jproduct)
                  Yj = Yadd(jproduct)
                  Zj = Zadd(jproduct)
                  Mytype = Typeadd(iproduct)
                  Yourtype = Typeadd(jproduct)
                  Call Epair(Ib,Vir,Upot(4),Xi,Yi,Zi,Xj,Yj,Zj
     $                 ,Mytype,Yourtype)
                  DeltaUrxn = DeltaUrxn + Upot(4)
                  DeltaVirrxn = DeltaVirrxn + Vir
               Endif
            Enddo
            Do ireactant=1,numreactants
               Xi = Xadd(iproduct)
               Yi = Yadd(iproduct)
               Zi = Zadd(iproduct)
               Xj = Rx(Listreactants(ireactant))
               Yj = Ry(Listreactants(ireactant))
               Zj = Rz(Listreactants(ireactant))
               Mytype = Typeadd(iproduct)
               Yourtype = Types(Listreactants(ireactant))
               Call Epair(Ib,Vir,Upot(5),Xi,Yi,Zi,Xj,Yj,Zj
     $              ,Mytype,Yourtype)
                  DeltaUrxn = DeltaUrxn - Upot(5)
                  DeltaVirrxn = DeltaVirrxn - Vir
            Enddo
         Enddo

C        Acceptance?
         If(DeltaUrxn.lt.((-1/Beta)*Log(AccDenom/(AccNumer*Kequil)))) Then
            Laccept=.True.
         Else
            Call Accept((AccNumer/AccDenom)*Kequil*Dexp(-Beta*DeltaUrxn)
     &           ,Laccept)
         Endif

C        Accept Or Reject

         If(Laccept) Then

            Av1 = Av1 + 1.0d0

C           Bookkeeping numbers of things
            OldNpart = Npart

            Do ireactant=1,Nreactants(Irxn)
               Ireacttype=Reactants(Irxn,ireactant)
               Npboxtype(Ib,Ireacttype) = Npboxtype(Ib,Ireacttype) -
     &              ReactantStoich(Irxn,ireactant)
               Npbox(Ib) = Npbox(Ib) - ReactantStoich(Irxn,ireactant)
               Npart = Npart - ReactantStoich(Irxn,ireactant)
               Nparts(Ireacttype) = Nparts(Ireacttype) -
     &              ReactantStoich(Irxn,ireactant)

            Enddo
            Do iproduct=1,Nproducts(Irxn)
               Iprodtype=Products(Irxn,iproduct)
               Npboxtype(Ib,Iprodtype) = Npboxtype(Ib,Iprodtype) +
     &              ProductStoich(Irxn,iproduct)
               Npbox(Ib) = Npbox(Ib) + ProductStoich(Irxn,iproduct)
               Npart = Npart + ProductStoich(Irxn,iproduct)
               Nparts(Iprodtype) = Nparts(Iprodtype) +
     &              ReactantStoich(Irxn,ireactant)
            Enddo
            Etotal(Ib) = Etotal(Ib) + DeltaUrxn
            Vtotal(Ib) = Vtotal(Ib) + DeltaVirrxn

C           All the bookkeeping

C           Establish product indices in place of reactants
            !Write(6,*) 'Product indices'
            If(numproducts.ge.numreactants) Then
               Do iproduct=1,numproducts
                  If(iproduct.le.numreactants) Then
                     Rx(Listreactants(iproduct)) = Xadd(iproduct)
                     Ry(Listreactants(iproduct)) = Yadd(iproduct)
                     Rz(Listreactants(iproduct)) = Zadd(iproduct)
                     Ibox(Listreactants(iproduct)) = Ib
                     Types(Listreactants(iproduct)) = Typeadd(iproduct)
                     !Write(6,*) Listreactants(iproduct),' filled in'
                  Else
                     Rx(OldNpart+iproduct-numreactants) = Xadd(iproduct)
                     Ry(OldNpart+iproduct-numreactants) = Yadd(iproduct)
                     Rz(OldNpart+iproduct-numreactants) = Zadd(iproduct)
                     Ibox(OldNpart+iproduct-numreactants) = Ib
                     Types(OldNpart+iproduct-numreactants) = Typeadd(iproduct)
                     !Write(6,*) OldNpart+iproduct-numreactants,' added'
                  Endif
               Enddo
            Else
               Do iproduct=1,numproducts
                  Rx(Listreactants(iproduct)) = Xadd(iproduct)
                  Ry(Listreactants(iproduct)) = Yadd(iproduct)
                  Rz(Listreactants(iproduct)) = Zadd(iproduct)
                  Ibox(Listreactants(iproduct)) = Ib
                  Types(Listreactants(iproduct)) = Typeadd(iproduct)
                  !Write(6,*) Listreactants(iproduct),' filled in'
                  Listreactants(iproduct)=-1
               Enddo
C           Fill in holes left by removing reactant particles
               HoleCounter = -1
               Do ireactant=(numproducts+1),numreactants
                  HoleCounter = HoleCounter + 1
C                 Move to next particle if current index is a reactant
                  Do While ( ANY( Listreactants==(OldNpart-HoleCounter) ))
                     HoleCounter = HoleCounter + 1
                  Enddo
C                 Fill in holes in reactant list
                  If(Listreactants(ireactant).le.OldNpart) Then
                     Rx(Listreactants(ireactant)) = Rx(OldNpart - HoleCounter)
                     Ry(Listreactants(ireactant)) = Ry(OldNpart - HoleCounter)
                     Rz(Listreactants(ireactant)) = Rz(OldNpart - HoleCounter)
                     Ibox(Listreactants(ireactant)) = Ibox(OldNpart - HoleCounter)
                     Types(Listreactants(ireactant)) = Types(OldNpart - HoleCounter)
C                    Write(6,*) OldNpart-Holecounter,' in place of
C    &                    ',Listreactants(ireactant)
                  Endif
               Enddo
            Endif

C        Increment number of attempted moves
         Av2 = Av2 + 1.0d0

         Endif
      Else
C        Not reacting
C        Increment number of attempted moves
         Av2 = Av2 + 1.0d0
      Endif
      Return
      End
