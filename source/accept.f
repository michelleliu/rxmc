      Subroutine Accept(Factor,Laccept)
      Implicit None

C     Test Acceptance; Draw Random Number For This

      Double Precision Randomnumber,Factor
      Logical Laccept

      If(Factor.Gt.1.0d0) Then
         Laccept = .True.
      Elseif(Factor.Lt.1.0d-10) Then
         Laccept = .False.
      Else
         If(Randomnumber().Lt.Factor) Then
            Laccept = .True.
         Else
            Laccept = .False.
         Endif
      Endif

      Return
      End
