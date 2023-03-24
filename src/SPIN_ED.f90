MODULE SPIN_ED
  USE ED_INPUT_VARS, only: &
       ed_read_input , &
       Nsites        , &
       Norb          , &
       temp          , &
       eps           , &
       wini          , &
       wfin          , &
       Lmats         , &
       Lreal         , &
       LOGfile   

  USE ED_GRAPH_MATRIX, only: &
       ed_Hij_init     => Hij_init     , &
       ed_Hij_add_link => Hij_add_link , &
       ed_Hij_read     => Hij_read     , &
       ed_Hij_get      => Hij_get      , &
       ed_Hij_local    => Hij_local    , &
       ed_Hij_info     => Hij_info     , &
       ed_Hij_write    => Hij_write


  USE ED_MAIN, only: &
       ed_init_solver , &
       ed_solve


END MODULE SPIN_ED

