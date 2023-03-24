  do i=MpiIstart,MpiIend
     mup = Hsector%H(1)%map(i)
     Nup = bdecomp(mup,Ns)
     Sz  = Nup - 0.5d0
     !
     htmp = zero
     !> H_Imp: Diagonal Elements
     do io=1,Ns
        htmp = htmp + Hdiag(io)*Sz(io)
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0d,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0d,htmp,i,i)
     end select
     !
  enddo

