C     Global Variables
C
C     Maxpart  = Maximum Number Of Particles
C     Npart    = Number Of Particles
C     Ntype    = Number Of Types
C     Npboxtype= Number Of Type in Box
C     MolFrac  = Fraction of particles in ideal gas
C     Ibox     = In Which Box Is A Particle ?
C     Nbox     = Number Of Boxes (2 For Gibbs Ensemble, 1 For The Other Ones)
C     Rx/Ry/Rz = Coordinates
C     Ecut     = Cut-Off Energy
C     Rcutsq   = Cut-Off Radius Squared
C     Beta     = 1/Temperature
C     Etotal   = Total Energy Of A Box
C     Vtotal   = Total Virial Of A Box
C
C     Nrxn           = Number Of Reactions (Reaction Ensemble MC only)
C     Nreactants     = Number Of Reactants
C     Nproducts      = Number Of Products
C     Reactants      = Names Of Reactants (one entry per reactant)
C     ReactantStoich = Stoichiometric Coefficient Of Reactants
C     ProductsStoich = Stoichiometric Coefficient Of Products
C
C     Maxtype
C     Maxrxn
C     Eps()
C     Sig()
C     Verbose
C     Nparts()
C     Types()
C
C     Nprint   = Trajectory Print Frequency

      Integer Maxpart,Maxtype,Maxrxn

      Logical Verbose

      Parameter (Maxpart = 1024)
      Parameter (Maxtype = 10)
      Parameter (Maxrxn = 5)

      Double Precision Eps(Maxtype,Maxtype),Sig(Maxtype,Maxtype)

      Double Precision Rx(Maxpart),Ry(Maxpart),Rz(Maxpart),Box(2)
     $     ,Ecut(Maxtype,Maxtype),Rcutsq,Beta,Etotal(2),Vtotal(2)
     $     ,MolFrac(Maxtype)

      Integer Npart,Nparts(Maxtype),Ibox(Maxpart),Types(Maxpart)
     $     ,Npbox(2),Npboxtype(2,Maxtype),Nbox,Nrxn,Ntype
     $     ,Nreactants(Maxrxn),Nproducts(Maxrxn)
     $     ,Reactants(Maxrxn,4),Products(Maxrxn,4)
     $     ,ReactantStoich(Maxrxn,4),ProductStoich(Maxrxn,4)

      Integer Nprint

      Common /Mycomms/ Rx,Ry,Rz,MolFrac
     $     ,Box,Ecut,Rcutsq,Beta,Etotal,Eps,Sig
     $     ,Vtotal,Npart,Ibox,Npbox,Npboxtype,Nbox
     $     ,Nprint
     $     ,Nrxn,Ntype,Nparts,Types,Verbose
     $     ,Nreactants,Nproducts
     $     ,Reactants,Products,ReactantStoich,ProductStoich
