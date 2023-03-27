  !H = sum_i H(i)*Sz(i) + sum_{i,j} Hij * (Sz(i)*Sz(j) + 0.5*Sp(i)Sm(j) + 0.5*Sm(i)*Sp(j))
  do i=MpiIstart,MpiIend
     m   = Hsector%H(1)%map(i)
     Nup = Bdecomp(m,Ns)
     Sz  = Nup - 0.5d0
     !
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
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0d,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0d,htmp,i,i)
     end select
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
              j = binary_search(Hsector%H(1)%map,k2)              
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
           !
           Jcondition = (Nup(io)==1) .AND. (Nup(jo)==0)
           if (Jcondition) then
              call Splus(jo,m,k1)
              call Sminus(io,k1,k2)
              j = binary_search(Hsector%H(1)%map,k2)
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
        enddo
     enddo
     !
  enddo
