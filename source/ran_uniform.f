      Function Randomnumber()
      Implicit None

C     Ran3 From Numerical Recipies

      Integer Mbig,Mseed,Mj,Mk,Ma(55),T1,I,Iff,Ii,Inext,Inextp,K

      Double Precision Randomnumber,Fac
      
      Parameter (Mbig=1000000000,Mseed=161803398,Fac=1.0D-9)

      Save Iff,Inext,Inextp,Ma
      Data Iff /0/

      If(Iff.Eq.0)Then
         Iff=1

         Call System_Clock(T1)
         Call System_Clock(T1)

         Mj=Abs(Mseed-Abs(T1))
         Mj=Mod(Mj,Mbig)

         Ma(55)=Mj
         Mk=1
         Do I=1,54
            Ii=Mod(21*I,55)
            Ma(Ii)=Mk
            Mk=Mj-Mk
            If(Mk.Lt.0)Mk=Mk+Mbig
            Mj=Ma(Ii)
         Enddo
         Do K=1,4
            Do I=1,55
               Ma(I)=Ma(I)-Ma(1+Mod(I+30,55))
               If(Ma(I).Lt.0)Ma(I)=Ma(I)+Mbig
            Enddo
         Enddo
         
         Inext=0
         Inextp=31
         
         Do T1=1,1000
            Inext=Inext+1
            If(Inext.Eq.56)Inext=1
            Inextp=Inextp+1
            If(Inextp.Eq.56)Inextp=1
            Mj=Ma(Inext)-Ma(Inextp)
            If(Mj.Lt.0)Mj=Mj+Mbig
            Ma(Inext)=Mj
         Enddo
      Endif

      Inext=Inext+1
      If(Inext.Eq.56)Inext=1
      Inextp=Inextp+1
      If(Inextp.Eq.56)Inextp=1
      Mj=Ma(Inext)-Ma(Inextp)
      If(Mj.Lt.0)Mj=Mj+Mbig
      Ma(Inext)=Mj
      Randomnumber=Dble(Mj)*Fac
      Return
      End
