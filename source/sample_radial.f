      Subroutine Sample_Radial(Ichoise)
      Implicit None

      Include 'commons.inc'

C     Samples The Radial Distribution Function

      Integer I,J,Maxx,Ichoise,A

      Parameter(Maxx = 500)

      Double Precision Ggt,Gg(Maxx),Delta,R2,Dx,Dy,Dz,Bx,Hbx

      Save Ggt,Gg,Delta

      Bx  = Box(1)
      Hbx = 0.5d0*Box(1)

      If(Ichoise.Eq.1) Then
         Do I=1,Maxx
            Gg(I) = 0.0d0
         Enddo

         Ggt   = 0.0d0
         Delta = Dble(Maxx-1)/Hbx
      
      Elseif(Ichoise.Eq.2) Then

         Ggt = Ggt + 1.0d0

C     Loop Over All Pairs

         Do I=1,Npart-1
            Do J=I+1,Npart

               Dx = Rx(I) - Rx(J)
               Dy = Ry(I) - Ry(J)
               Dz = Rz(I) - Rz(J)
 
               If (Dx.Gt.Hbx) Then
                  Dx = Dx - Bx
               Elseif (Dx.Lt. - Hbx) Then
                  Dx = Dx + Bx
               Endif
 
               If (Dy.Gt.Hbx) Then
                  Dy = Dy - Bx
               Elseif (Dy.Lt. - Hbx) Then
                  Dy = Dy + Bx
               Endif
 
               If (Dz.Gt.Hbx) Then
                  Dz = Dz - Bx
               Elseif (Dz.Lt. - Hbx) Then
                  Dz = Dz + Bx
               Endif

               R2 = Dsqrt(Dx*Dx + Dy*Dy + Dz*Dz)

C     Calculate In Which Bin This Interaction Is In
C
C     Delta = 1/Binsize
C     Maxx  = Number Of Bins

               If(R2.Lt.Hbx) Then
                  
                  A     = Idint(R2*Delta) + 1
                  Gg(A) = Gg(A)           + 2.0d0
               
               Endif
            Enddo
         Enddo

      Else

C     Write Results To Disk

         Ggt   = 1.0d0/(Ggt*Dble(Npart))
         Delta = 1.0d0/Delta

         Open(21,File="Radial",Status="Unknown")
         Do I=1,Maxx-1
            R2 = (16.0d0*Datan(1.0d0)/3.0d0)*
     &           (Dble(Npart-1)/(Bx**3))*(Delta**3)*
     &           ((Dble(I))**3 - (Dble(I-1)**3))
            
            Write(21,*) ((Dble(I)-0.5d0)*Delta),Gg(I)*Ggt/R2
         Enddo
         Close(21)
      Endif

      Return
      End
