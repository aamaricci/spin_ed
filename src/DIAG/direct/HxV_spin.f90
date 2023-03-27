  !H = sum_i H(i)*Sz(i) + sum_{i,j} Hij * (Sz(i)*Sz(j) + 0.5*Sp(i)Sm(j) + 0.5*Sm(i)*Sp(j))
  htmp = zero
  !
  !Sum_i H_i*Sz(i)
  do io=1,Ns
     htmp = htmp + Hdiag(io)*Sz(io)
  enddo
  !
  !Sum_ij H_ij Sz(i)*Sz(j)
  do io=1,Ns
     do jo=io,Ns
        if(Hij(io,jo)==zero)cycle
        htmp = htmp - Hij(io,jo)*Sz(io)*Sz(jo)
     enddo
  enddo
  i = j
  hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
  !
  !Sum_ij H_ij 0.5*Sp(i)Sm(j) + 0.5*Sm(i)*Sp(j)
  do io=1,Ns
     do jo=io,Ns
        if(Hij(io,jo)==zero)cycle
        htmp = Hij(io,jo)/2d0
        Jcondition = (Nup(jo)==1) .AND. (Nup(io)==0)
        if (Jcondition) then
           call Sminus(jo,m,k1)
           call Splus(io,k1,k2)
           i = binary_search(Hsector%H(1)%map,k2)              
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        !
        Jcondition = (Nup(io)==1) .AND. (Nup(jo)==0)
        if (Jcondition) then
           call Splus(jo,m,k1)
           call Sminus(io,k1,k2)
           i = binary_search(Hsector%H(1)%map,k2)
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo


