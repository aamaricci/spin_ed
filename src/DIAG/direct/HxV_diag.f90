  htmp = zero
  !
  !> H_Imp: Diagonal Elements
  do io=1,Ns
     htmp = htmp + Hdiag(io)*Sz(io)
  enddo
  !
  i = j
  hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
