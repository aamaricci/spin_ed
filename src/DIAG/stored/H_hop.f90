  do i=MpiIstart,MpiIend
     m   = Hsector%H(1)%map(i)
     Nup = Bdecomp(m,Ns)
     Sz  = Nup - 0.5d0
     !
     htmp = zero
     !
     !> H_hop: Off-diagonal elements, i.e. non-local part. 
     !remark: io=jo cant have simultaneously n=0 and n=1 (Jcondition)
     !        so diagonal element (in H_local) are neglected
     !H = sum_{i,j} Hij * (Sz(i)*Sz(j) + 0.5*Sp(i)Sm(j) + 0.5*Sm(i)*Sp(j))
     do io=1,Ns
        do jo=1,Ns
           ! Sz(i)*Sz(j)
           Jcondition = (Hij(io,jo)/=zero)
           if (Jcondition) then
              htmp = htmp - 2d0*Hij(io,jo)*Sz(io)*Sz(jo)
           end if
           !
           ! 0.5*Sp(i)Sm(j)
           Jcondition = &
                (Hij(io,jo)/=zero) .AND. (Nup(jo)==1) .AND. (Nup(io)==0)
           if (Jcondition) then
              call Sminus(jo,m,k1)
              call Splus(io,k1,k2)
              j = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(Hij(io,jo))/2d0
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
           !
           ! 0.5*Sm(i)*Sp(j)
           Jcondition = &
                (Hij(io,jo)/=zero) .AND. (Nup(io)==1) .AND. (Nup(jo)==0)
           if (Jcondition) then
              call Splus(jo,m,k1)
              call Sminus(io,k1,k2)
              j = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(Hij(io,jo))/2d0
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
