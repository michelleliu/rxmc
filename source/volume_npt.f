      Subroutine Volume_Npt(Av1,Av2,Delta,Press)
      Implicit None

C     Volume Change In the Npt Ensemble

      Include 'commons.inc'

      Logical Laccept
      Integer I
      Double Precision Rxo(Maxpart),Ryo(Maxpart),Rzo(Maxpart),Boxold
     $     ,Volold,Volnew,Df,Unew,Virnew,Delta,Av1,Av2,Press,Randomnumber

      Boxold = Box(1)
      Volold = Boxold**3
      Volnew = Volold + (2.0d0*Randomnumber()-1.0d0)*Delta

      If(Volnew.Le.0.0d0) Return

C     Scale Boxlength

      Box(1) = Volnew**(1.0d0/3.0d0)
      Df     = Box(1)/Boxold

      If(Dsqrt(Rcutsq).Gt.0.5d0*Box(1))
     &     Stop "Error Volume; Rcut Is Too Smal !!!"

C     Transform Coordinates

      Do I=1,Npart
         Rxo(I) = Rx(I)
         Ryo(I) = Ry(I)
         Rzo(I) = Rz(I)

         Rx(I) = Rx(I)*Df
         Ry(I) = Ry(I)*Df
         Rz(I) = Rz(I)*Df
      Enddo

C     Compute New Energy

      Call Etot(1,Virnew,Unew)

      Call Accept(Dexp(-Beta*(Unew - Etotal(1) +
     &     Press*(Volnew-Volold) -
     &     Dble(Npart)*Dlog(Volnew/Volold)/Beta)),Laccept)

      Av2 = Av2 + 1.0d0

C     Accept Or Reject

      If(Laccept) Then
         Av1     = Av1 + 1.0d0

         Etotal(1) = Unew
         Vtotal(1) = Virnew
      Else

C     Reject, Restore Cordinates

         Box(1) = Boxold

         Do I=1,Npart
            Rx(I) = Rxo(I)
            Ry(I) = Ryo(I)
            Rz(I) = Rzo(I)
         Enddo
      Endif

      Return
      End
