c!
c!   *..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
c!   | .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
c!   | |                                                                                      | |
c!   | |                                                                                      | |
c!   | |                   HSA-DC -- Hybrid Simulated Annealing approach for                  | |
c!   | |                           global minimization of DC functions                        | | 
c!   | |                                                                                      | |
c!   | |                                 (version 1)                                          | |
c!   | |                                                                                      | |
c!   | |                                                                                      | |
c!   | |                       by Kaisa Joki and Adil M. Bagirov                              | |
c!   | |                                                                                      | |
c!   | |                           (Last modified June 2025)                                  | |
c!   | |                                                                                      | |       
c!   | |                                                                                      | |
c!   | |                                                                                      | |
c!   | |     The software is free for academic teaching and research purposes, but we         | |
c!   | |     ask you to refer the reference given below, if you use it.                       | |
c!   | |                                                                                      | |
c!   | |                                                                                      | |
c!   | .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
c!   *..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*  
c!   |                                                                                          |
c!   |                                                                                          |
c!   |   Utilizes the aggregate subgradient method (AGGSUB) from [2] as a local search method.  |
c!   |                                                                                          |
c!   |   Utilizes the augmented subgradient method (AUGSUB) from [3] during the last iteration  |
c!   |   round to further improve the solution.                                                 |
c!   |                                                                                          |
c!   |   Utilizes the quadratic programming problem algorithm from [4].                         |
c!   |                                                                                          |
c!   |                                                                                          |
c!   | .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
c!   *..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
c!   |                                                                                          |
c!   |                                                                                          |
c!   |   References:                                                                            |
c!   |                                                                                          |
c!   |   [1] Bagirov, A.M., Joki, K., Mäkelä, M.M, Taheri, S. "A hybrid simulated annealing     |
c!   |       approach for global minimization of DC functions". (manuscript)                    |
c!   |                                                                                          |
c!   |   [2] Bagirov, A.M., Taheri, S., Joki, K., Karmitsa, N., Mäkelä, M.M.: "Aggregate        |
c!   |       subgradient method for nonsmooth DC optimization". Optim. Lett. 15(1),             | 
c!   |       pp. 83–96 (2021).                                                                  | 
c!   |                                                                                          |
c!   |   [3] Bagirov, A.M., Hoseini Monjezi, N., Taheri, S.: "An augmented subgradient          |
c!   |       method for minimizing nonsmooth DC functions". Comput. Optim. Appl. 80(1),         |
c!   |       pp. 411–438 (2021).                                                                |
c!   |                                                                                          |
c!   |   [4] Wolfe, P.H.: "Finding the nearest point in a polytope". Math. Program. 11(2),      |
c!   |       pp. 128–149 (1976).                                                                |
c!   |                                                                                          |
c!   |                                                                                          |
c!   | .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
c!   *..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*  
c!
c!---------------------------------------------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision spoint5_in(20,1000)   ! Original starting points from the boxes [-5,5] (row is a starting point)
      double precision spoint5(1000,20)      ! Used starting points from the boxes [-5,5] (column is a starting point)
      double precision spoint10_in(20,1000)  ! Original starting points from the boxes [-10,10] (row is a starting point)
      double precision spoint10(1000,20)     ! Used starting points from the boxes [-10,10] (column is a starting point)
      double precision spoint100_in(20,1000) ! Original starting points from the boxes [-100,100] (row is a starting point)
      double precision spoint100(1000,20)    ! Used starting points from the boxes [-100,100] (column is a starting point)
      
      COMMON /cnselect/nselect                                         ! Label indicating the problem solved inside the selected group of test problems
      COMMON /cstart5/spoint5,/cstart10/spoint10,/cstart100/spoint100  ! Starting points
      COMMON /cnprob/nprob0                                            ! Specifies whether a specific function is used in the test problem selected from Group 2
      COMMON /cgr/ngroup                                               ! Selected group of test problems (1="Group 1", 2="Group 2" and 3="Group 3")

      COMMON /caugsub/maxaugsubit        ! If no improvement during 'maxaugsubit' consecutive iterations, then AUGSUB is terminated
      COMMON /caugsubtol/augsub_stop_tol ! The tolerance used in the improvement check in 'maxaugsubit' consecutive iterations during AUGSUB
      COMMON /caggsub/maxaggsubit        ! If no improvement during 'maxaggsubit' consecutive iterations, then AGGSUB is terminated
      COMMON /caggsubtol/aggsub_stop_tol ! The tolerance used in the improvement check in 'maxaggsubit' consecutive iterations during AGGSUB
      
      COMMON /cnmax/nxpmax,nxbmax        ! Only 'nxpmax' ('nxbmax') number of the newest points used from 'xp' ('xb') table in comparisons.
      COMMON /cnonp/nonp                 ! If nonp = 1, then table 'xp' contains points only from the current round of globsub(). In addition, 'nxpmax' is not used.
                                         ! If nonp = 0, then table 'xp' contains all points generated during the execution but only 'nxpmax' newest ones are used.
      COMMON /cnonb/nonb                 ! If nonb = 1, then no previous solutions are used from table 'xb'. In addition, 'nxbmax' is not used.
                                         ! If nonb = 0, then previous solutions are used from table 'xb' (NOTE: only 'nxbmax' newest ones are used).

      COMMON /ctlimit/tlimit             ! Time limit
      COMMON /cpen/pen                   ! Penalization parameter for violating box constraints
      
      COMMON /cnlimit/nlimit1,nlimit2    ! Labels to select the correct values for the number of eps-subgradients for f1 and f2 (1="light", 2="heavier", 3="heavy")
                                        
      COMMON /crttemp/Tdec,Tinc,TEMP_in  ! Decrease and increase parameter for the temperature and initial temperature
      COMMON /cnTEMPupdate/nTEMPupdate   ! Is temperature updated (value>0) or not (value=0)
      
      COMMON /cnway/nway1,nway2          ! Directions for f1 and f2
      
c!--------------------------------------------      
c! Select the correct group of test problems
c!--------------------------------------------      
c!      ngroup = 1   ! Group 1
c!      ngroup = 2   ! Group 2 
      ngroup = 3   ! Group 3

c!---------------------------------------------------------------------
c!  -- START: Some USER selected parameters are set here --
c!---------------------------------------------------------------------
 
c! --- TEMPERATURE --- 
      Tdec = 0.75D+00               ! Decrease parameter for temperature
      Tinc = 2.0D+00                ! Increase parameter for temperature
      TEMP_in = 1.0d+4              ! Initial temperature     
      nTEMPupdate = 1               ! If nTEMPupdate = 0 temperature is not updated and if nTEMPupdate >0 temperature is updated with a certain procedure at the end of each globsub(). 
 
c! --- BOX CONSTRAINTS ---
      pen = 1.0d+03                 ! Penalization parameter for the box constraints
      
c! --- LABELs for the NUMBER of eps-subgradients for f1 and f2 ---
      nlimit1 = 1                   ! 1="light"=min(100,2*m), 2="heavier"=min(200,2*m) and 3="heavy"=min(400,2*m)
      nlimit2 = 3                   ! 1="light"=min(100,2*m), 2="heavier"=min(200,2*m) and 3="heavy"=min(400,2*m)
                                    ! Note: m = "the dimension of the problem"
      
c! --- DIRECTIONS for f1 and f2 ---
      nway1 = 0   
      nway2 = 2   
      
c!    Options for directions:     
c!      0 = Unit directions in the order of components (always +1 and -1 for each component, the maximum number of directions with this option is 2*m)
c!      1 = Random unit directions in the order of components (always randomly generates whether +1 or -1 for each component, the maximum number of directions with this option is m)
c!      2 = Half random unit directions in the order of components and half completely random directions (the maximum number of directions with this option is 2*m)
c!      3 = Completely random directions (the maximum number of directions with this option is 4*m)
c!      4 = All unit directions in the order of components + some random ones (if nlimit > 2*m, the maximum number of directions with this option is 4*m)
c!      Note: the values of 'nlimit1' and 'nlimit2' can be at most equal to the maximum number of directions specified in the options above.
      
 
c! --- AUGSUB method ---
      maxaugsubit = 25              ! If no improvement during 'maxaugsubit' consecutive iterations, then AUGSUB is terminated
      augsub_stop_tol = 1.0D-05     ! The tolerance used in the improvement check in 'maxagsubit' consecutive iterations
    
c! --- AGGSUB method ---   
      maxaggsubit = 10              ! If no improvement during 'maxaggsubit' consecutive iterations and the objective value is worse than the current best one, then AGGSUB is terminated
      aggsub_stop_tol = 1.0d-02     ! The tolerance used in the improvement check in 'maxaggsubit' consecutive iterations

c! --- POINTS used in comparisons ---      
      nxpmax = 1000                 ! Only 'nxpmax' number of the newest points used from 'xp' table in comparisons
      nxbmax = 200                  ! Only 'nxbmax' number of the newest points used from 'xb' table in comparisons
      
      nonp = 1                      ! If nonp = 1, then table 'xp' contains points only from current round of globsub(). In addition, 'nxpmax' is not used
                                    ! If nonp = 0, then table 'xp' contains all points generated during execution but only 'nxpmax' newest ones are used
                                    
      nonb = 0                      ! If nonb = 1, then no previous solutions are used from table 'xb'. In addition, 'nxbmax' is not used
                                    ! If nonb = 0, then the previous solutions are used from table 'xb' (NOTE: only 'nxbmax' newest ones are used)
                                    
c! --- TIME LIMIT ---      
      tlimit = 7.2d+03              ! Time limit in seconds
      
c!---------------------------------------------------------------------
c!  -- END: Some USER selected parameters are set here --
c!---------------------------------------------------------------------

      SELECT CASE(ngroup)
         CASE(1) 
          open(42,file='HSA-DC_G1.txt')
        
         CASE(2) 
          open(42,file='HSA-DC_G2.txt')
          
         CASE(3) 
          open(42,file='HSA-DC_G3.txt')   
      END SELECT
      
      WRITE(42,*) 'HSA-DC method, Group', ngroup
      WRITE(42,*) '--------------------------------------------------'
      WRITE(42,*) 'Temperature decrese paramter =', Tdec
      WRITE(42,*) 'Temperature increase paramter =', Tinc
      WRITE(42,*) 'Initial temperature =',TEMP_in
      IF (nTEMPupdate==0) THEN 
         WRITE(42,*) 'Temperature not updated.'
      ELSE   
         WRITE(42,*) 'Temperature is updated.'
      END IF     
      WRITE(42,*) 'penalization parameter =',pen
      WRITE(42,*) 'tlimit=', tlimit
      WRITE(42,*) '--------------------------------------------------'  
      WRITE(42,*) 'Inside globsub()'      
      WRITE(42,*) 'Label for number of eps-subgrad. of f1 = ',nlimit1     
      WRITE(42,*) 'Label for number of eps-subgrad. of f2 = ',nlimit2     
      WRITE(42,*) 'Label for directions used for f1 = ',nway1     
      WRITE(42,*) 'Label for directions used for f2 = ',nway2     
      WRITE(42,*) '--------------------------------------------------'    
      WRITE(42,*) 'AUGSUB method'    
      WRITE(42,*) 'maxaugsubit =',maxaugsubit
      WRITE(42,*) 'augsub_stop_tol =',augsub_stop_tol
      WRITE(42,*) '--------------------------------------------------'   
      WRITE(42,*) 'AGGSUB method'    
      WRITE(42,*) 'maxaggsubit =',maxaggsubit
      WRITE(42,*) 'aggsub_stop_tol =',aggsub_stop_tol
      WRITE(42,*) '--------------------------------------------------'        
      WRITE(42,*) 'nxpmax=',nxpmax
      IF (nonp==1) THEN
         WRITE(42,*) 'Table xp contains points only from the current', 
     1 ' round of globsub(). In addition, nxpmax is not used.'
      ELSE 
         WRITE(42,*) 'Table xp contains all points generated during ', 
     1 ' the execution of the method, but only nxpmax newest',
     2 ' ones are used.'
      END IF      
      WRITE(42,*) '--------------------------------------------------'  
      WRITE(42,*) 'nxbmax=',nxbmax
      IF (nonb==1) THEN
         WRITE(42,*) 'No previous solutions are used from table xb.', 
     1 ' In addition, nxbmax is not used.'
      ELSE 
         WRITE(42,*) 'The previous solutions are used from table xb', 
     1 ' (NOTE: only nxbmax newest ones are used).'
      END IF          
      WRITE(42,*) '--------------------------------------------------'
      WRITE(42,*) 
      

      
c!-----------------------------------------           
c! Starting points       
c!-----------------------------------------      
      open(77,file='startpoints5.txt',status='old',form='formatted')
      open(78,file='startpoints10.txt',status='old',form='formatted')
      open(79,file='startpoints100.txt',status='old',form='formatted')
c!=====================================================================      
      do i=1,20
       read(77,*) (spoint5_in(i,k),k=1,1000)
       npoint=i
      end do
      WRITE(42,*) 'Input start5. Number of records: ',npoint
c!------------------------------------
      do i=1,20
       read(78,*) (spoint10_in(i,k),k=1,1000)
       npoint=i
      end do
      WRITE(42,*) 'Input start10. Number of records: ',npoint
c!------------------------------------
      do i=1,20
       read(79,*) (spoint100_in(i,k),k=1,1000)
       npoint=i
      end do
      WRITE(42,*) 'Input start100. Number of records: ',npoint
      WRITE(42,*)
c!====================================================================
      CLOSE(77)
      CLOSE(78)
      CLOSE(79)
c!====================================================================   
c! -- Transposes of the starting point matrices --  
c!    Note: After transposes a column is a starting point    
      do i = 1, 20
        do k = 1,1000
            spoint5(k,i) = spoint5_in(i,k)      
            spoint10(k,i) = spoint10_in(i,k)      
            spoint100(k,i) = spoint100_in(i,k)      
        end do
      end do
c!====================================================================
      SELECT CASE(ngroup)
c!------------------------------------
      CASE(1)  ! Group 1
      write(42,141) 
141   format('  Prob#','   n','              Init. val.'
     1 ,'      Final val.','    CPU','     f1-eval','   f2-eval'
     1 ,'   grad-f1','   grad-f2','   #esc','    #st.p','    #ntemp')
      write(42,*)  
c!------------------------------------
      CASE(2)  ! Group 2
      write(42,142) 
142   format('  Prob1#','  Prob2#','   n','              Init. val.'
     1 ,'      Final val.','    CPU','     f1-eval','   f2-eval'
     1 ,'   grad-f1','   grad-f2','   #esc','    #st.p','    #ntemp')
      write(42,*) 
c!------------------------------------
      CASE(3)  ! Group 3
      write(42,143) 
143   format('  Prob#','   n','              Init. val.'
     1 ,'      Final val.','    CPU','     f1-eval','   f2-eval'
     1 ,'   grad-f1','   grad-f2','   #esc','    #st.p','    #ntemp')
      write(42,*)
      END SELECT
c!====================================================================
  
c!====================================================================
      nprob=0        ! Counter telling how many problems solved so far
c!--------------------------------------------------------------------
c! Random number seed selection
c!--------------------------------------------------------------------
      call random_seed()    
c!--------------------------------------------------------------------   

c!--------------------------------------------------------------------
c!  Initial information about the test problem groups and their problems
c!--------------------------------------------------------------------
c!====================================================================
      SELECT CASE(ngroup)
c!------------------------------------
      CASE(1)  ! Group 1
c!------------------------------------
      do i=1,8
       nselect=i
       if(i.eq.1) then
           n=2
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.2) then
           n=4
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if    
       if(i.eq.3) then
           n=2
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.4) then
           n=3
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.5) then
           do j=1,6
            if(j.eq.1) n=2 
            if(j.eq.2) n=5 
            if(j.eq.3) n=10 
            if(j.eq.4) n=50
            if(j.eq.5) n=100
            if(j.eq.6) n=200
            nprob=nprob+1
            print *,nprob
            call esub(n)
           end do
       end if
       if(i.eq.6) then
           n=3
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.7) then
           do j=1,6
            if(j.eq.1) n=2 
            if(j.eq.2) n=5 
            if(j.eq.3) n=10 
            if(j.eq.4) n=50
            if(j.eq.5) n=100
            if(j.eq.6) n=200
            nprob=nprob+1
            print *,nprob
            call esub(n)
           end do
       end if
       if(i.eq.8) then
           do j=1,5
            if(j.eq.1) n=5 
            if(j.eq.2) n=10 
            if(j.eq.3) n=50
            if(j.eq.4) n=100
            if(j.eq.5) n=200
            nprob=nprob+1
            print *,nprob
            call esub(n)
           end do
       end if
      end do 
      
c!------------------------------------
      CASE(2)  ! Group 2
c!------------------------------------
      do i=1,20
       nselect=i
       if(i.eq.1) then
        do j=1,5
          if(j.eq.1) n=2
          if(j.eq.2) n=5
          if(j.eq.3) n=10
          if(j.eq.4) n=50
          if(j.eq.5) n=100
          nprob0=1
          nprob=nprob+1 
          print *,nprob
          call esub(n)
        end do
       end if
       if((i.eq.2).or.(i.eq.3)) then
        do j=1,6
          if(j.eq.1) n=2
          if(j.eq.2) n=5
          if(j.eq.3) n=10
          if(j.eq.4) n=50
          if(j.eq.5) n=100
          if(j.eq.6) n=200
          nprob0=0         
          nprob=nprob+1
          print *,nprob
          call esub(n)
        end do
       end if
       if(i.eq.4) then
        do j=1,5
          if(j.eq.1) n=2
          if(j.eq.2) n=5
          if(j.eq.3) n=10
          if(j.eq.4) n=50
          if(j.eq.5) n=100
          nprob0=1
          nprob=nprob+1
          print *,nprob
          call esub(n)
        end do
       end if       
       if(i.eq.5) then
        do j=1,6
          if(j.eq.1) n=2
          if(j.eq.2) n=5
          if(j.eq.3) n=10
          if(j.eq.4) n=50
          if(j.eq.5) n=100
          if(j.eq.6) n=200
          nprob0=0
          nprob=nprob+1
          print *,nprob
          call esub(n)
        end do
       end if
       if(i.eq.6) then
        do j=1,5
          if(j.eq.1) n=2
          if(j.eq.2) n=5
          if(j.eq.3) n=10
          if(j.eq.4) n=50
          if(j.eq.5) n=100
          nprob0=1         
          nprob=nprob+1
          print *,nprob
          call esub(n)
        end do
       end if
       if(i.ge.7) then
        do j=1,6
          if(j.eq.1) n=2
          if(j.eq.2) n=5
          if(j.eq.3) n=10
          if(j.eq.4) n=50
          if(j.eq.5) n=100
          if(j.eq.6) n=200
          nprob0=0
          nprob=nprob+1
          print *,nprob
          call esub(n)
        end do
       end if
      end do
c!------------------------------------
      CASE(3)  ! Group 3
c!------------------------------------
      do i=1,6
       if(i.eq.1) then
           nselect=1
           n=2
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.2) then
           nselect=2           
           do j=1,6
            if(j.eq.1) n=2 
            if(j.eq.2) n=5 
            if(j.eq.3) n=10 
            if(j.eq.4) n=50
            if(j.eq.5) n=100
            if(j.eq.6) n=200
            nprob=nprob+1
            print *,nprob
            call esub(n)
           end do
       end if       
       if(i.eq.3) then
           nselect=3
           n=2
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.4) then
           nselect=4          
           do j=1,6
            if(j.eq.1) n=2 
            if(j.eq.2) n=5 
            if(j.eq.3) n=10 
            if(j.eq.4) n=50
            if(j.eq.5) n=100
            if(j.eq.6) n=200
            nprob=nprob+1
            print *,nprob
            call esub(n)
           end do
       end if
       if(i.eq.5) then
           nselect=5
           n=2
           nprob=nprob+1
           print *,nprob
           call esub(n)
       end if
       if(i.eq.6) then
           nselect=6
           do j=1,6
            if(j.eq.1) n=2 
            if(j.eq.2) n=5 
            if(j.eq.3) n=10 
            if(j.eq.4) n=50
            if(j.eq.5) n=100
            if(j.eq.6) n=200
            nprob=nprob+1
            print *,nprob
            call esub(n)
           end do
       end if
      end do 

      END SELECT 
c!====================================================================    
      CLOSE(42)

      STOP
      END
c!====================================================================   


c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!                    START: HSA-DC method 
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

      subroutine esub(n)
      implicit double precision (a-h,o-z)
      INTEGER n                                   ! The dimension of the problem
      double precision xstart(n,20)               ! Table containing starting points (column is a starting point)
      double precision x(n), xl(n), xu(n)         ! Point x and its lower and upper bounds 
      double precision xbest(n), xbest1(n)        ! Help vectors 
      double precision xb(n,20000)                ! Table with all solutions (column is a solution)
      double precision xp(n,100000)               ! Table with all points (column is a point)
      double precision fb(1000)                   ! Objective value table 
      double precision fxb(20000)                 ! Objective value table of local solutions
      double precision spoint5(1000,20)           ! Starting points (column is a point)
      double precision spoint10(1000,20)          ! Starting points (column is a point)
      double precision spoint100(1000,20)         ! Starting points (column is a point)
         
      COMMON /cpen/pen,/cnp/np,/csm/dbox,sm
     1 ,/cnf1/nf1,/cnf2/nf2,/cngrad1/ngrad1,/cngrad2/ngrad2
     2 ,/cnselect/nselect,/cind2/ind2,/cnp2/np2,/cnp1/np1
     
      COMMON /cstart5/spoint5,/cstart10/spoint10,/cstart100/spoint100  ! Starting points
      
      COMMON /cnprob/nprob0              ! Specifies whether a specific function is used in the test problem from Group 2 (needed in some test problems)
      COMMON /cgr/ngroup                 ! The selected group of test problems (1="Group 1", 2="Group 2" and 3="Group 3")

      COMMON /clim/limit1,limit2         ! The maximum number of eps-subgradients for f1 and f2
      COMMON /cTEMP/TEMP                 ! The current temperature
      COMMON /cntemp/ntemp               ! The number of temperature updates
      COMMON /crttemp/Tdec,Tinc,TEMP_in  ! Decrease and increase parameter for the temperature and the initial temperature

      COMMON /ctlimit/tlimit             ! Time limit
      
      COMMON /cnlimit/nlimit1,nlimit2    ! The labels to select the correct values for the number of eps-subgradient for f1 and f2 (1="light", 2="heavier", 3="heavy")
      
      COMMON /cnewbest/newbest,newind    ! 'newbwst' is used to indicate whether globsub() already produces a better solution (0=no and 1=yes). 'newind' is the index of the solution 
 
      COMMON /cnmax/nxpmax,nxbmax   ! Only 'nxpmax' ('nxbmax') number of the newest points used from 'xp' ('xb') table in comparisons
      COMMON /cnonb/nonb            ! If nonb = 1, then no previous solutions are used from table 'xb'. In addition, 'nxbmax' is not used.
                                    ! If nonb = 0, then the previous solutions are used from table 'xb' (NOTE: only 'nxbmax' newest ones are used).
 
c!====================================================================
c!   Test problem parameters
c!====================================================================
      SELECT CASE(ngroup) 
c!------------------------------------
      CASE(1)  ! Group 1
c!------------------------------------
          SELECT CASE(nselect)  ! Test problem from Group 1 is selected
          
           CASE(1)                        ! P_1
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=2

           CASE(2)                        ! P_2
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=3

           CASE(3)                        ! P_3
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=7

           CASE(4)                        ! P_4
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=8

           CASE(5)                        ! P_5 
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=10

           CASE(6)                        ! P_6 
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=11

           CASE(7)                        ! P_7 
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=14

           CASE(8)                        ! P_8 
               do i=1,n
                 xl(i) = -5.0d+00
                 xu(i) = 5.0d+00
               end do
               nprob1=15
          END SELECT
c!------ Starting point ------  
          if(nselect.le.7) then
           do i=1,20
            do j=1,n
             xstart(j,i)=spoint100(j,i)
            end do
           end do
          end if 
          if(nselect.eq.8) then
           do i=1,20
            do j=1,n
             xstart(j,i)=spoint5(j,i)
            end do
           end do
          end if
c!------------------------------------
      CASE(2)  ! Group 2
c!------------------------------------
          SELECT CASE(nselect) ! Test problem from Group 2 is selected
           
           CASE(1)                        ! P_9
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=1
               nprob2=6
          
           CASE(2)                        ! P_10
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=1
               nprob2=8

           CASE(3)                        ! P_11
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=1
               nprob2=14

           CASE(4)                        ! P_12
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=2
               nprob2=7

           CASE(5)                        ! P_13
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=2
               nprob2=9
               
           CASE(6)                        ! P_14
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=3
               nprob2=7

           CASE(7)                        ! P_15
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=10
               nprob2=8

           CASE(8)                        ! P_16
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=13
               nprob2=3

           CASE(9)                        ! P_17
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=13
               nprob2=8

           CASE(10)                       ! P_18
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=13
               nprob2=10

           CASE(11)                       ! P_19
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=13
               nprob2=14

           CASE(12)                       ! P_20
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=13
               nprob2=17

           CASE(13)                       ! P_21
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=13
               nprob2=18

           CASE(14)                       ! P_22
               do i=1,n
                 xl(i) = -5.0d+00
                 xu(i) = 5.0d+00
               end do
               nprob1=16
               nprob2=1

           CASE(15)                       ! P_23
               do i=1,n
                 xl(i) = -5.0d+00
                 xu(i) = 5.0d+00
               end do
               nprob1=16
               nprob2=3

           CASE(16)                       ! P_24
               do i=1,n
                 xl(i) = -5.0d+00
                 xu(i) = 5.0d+00
               end do
               nprob1=16
               nprob2=14

           CASE(17)                       ! P_25
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=17
               nprob2=2

           CASE(18)                       ! P_26
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=17
               nprob2=4

           CASE(19)                       ! P_27
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=17
               nprob2=8

           CASE(20)                       ! P_28
               do i=1,n
                 xl(i) = -1.0d+02
                 xu(i) = 1.0d+02
               end do
               nprob1=18
               nprob2=8

          END SELECT
c!------ Starting point ------  
          if((nselect.le.13).or.(nselect.ge.17)) then
           do i=1,20
            do j=1,n
             xstart(j,i)=spoint100(j,i)
            end do
           end do
          end if 
          if((nselect.ge.14).and.(nselect.le.16)) then
           do i=1,20
            do j=1,n
             xstart(j,i)=spoint5(j,i)
            end do
           end do
          end if
c!------------------------------------
      CASE(3)  ! Group 3
c!------------------------------------
          SELECT CASE(nselect)
               
           CASE(1)                        ! P_29
               do i=1,n
                 xl(i) = -1.0d+01
                 xu(i) = 1.0d+01
               end do
               nprob1=1

           CASE(2)                        ! P_30
               do i=1,n
                 xl(i) = -1.0d+01
                 xu(i) = 1.0d+01
               end do
               nprob1=2

           CASE(3)                        ! P_31
               do i=1,n
                 xl(i) = -5.0d+00
                 xu(i) = 5.0d+00
               end do
               nprob1=3
               
           CASE(4)                        ! P_32
               do i=1,n
                 xl(i) = -dble(10)
                 xu(i) = dble(10)
               end do
               nprob1=4

           CASE(5)                        ! P_33
               do i=1,n
                 xl(i) = -1.0d+01
                 xu(i) = 1.0d+01
               end do
               nprob1=5
               
           CASE(6)                        ! P_34 
               do i=1,n
                 xl(i) = -1.0d+01
                 xu(i) = 1.0d+01
               end do
               nprob1=6

          END SELECT
c!------ Starting point ------      
          if((nselect.le.2).or.(nselect.ge.4)) then
           do i=1,20
            do j=1,n
             xstart(j,i)=spoint10(j,i)
            end do
           end do
          end if 
          if(nselect.eq.3) then
           do i=1,20
            do j=1,n
             xstart(j,i)=spoint5(j,i)
            end do
           end do
          end if
      END SELECT
c!====================================================================
c!    Some method parameters
c!====================================================================
      m = n                      ! Dimension of the problem  
 
! Maximum number of eps-subgradients for f1   
      SELECT CASE(nlimit1) 
        CASE(1)                  !Light version
         limit1 = min(100,2*m)   
        CASE(2)                  !Heavier version
         limit1 = min(200,2*m)               
        CASE(3)                  !Heavy version
         limit1 = min(400,2*m)           
      END SELECT     
      
! Maximum number of eps-subgradients for f2         
      SELECT CASE(nlimit2) 
        CASE(1)                  !Light version
         limit2 = min(100,2*m)   
        CASE(2)                  !Heavier version
         limit2 = min(200,2*m)               
        CASE(3)                  !Heavy version
         limit2 = min(400,2*m)           
      END SELECT           
        
      nbundle1 = min(500,5*m+3)     ! The bundle size for AGGSUB
      nbundle2 = min(70,m+3)        ! The bundle size for AUGSUB
      
      imax=1
c!====================================================================     
c!  Parameter initialization for some test problems (needed in Group 2)     
c!====================================================================
      if (ngroup==2) then  
         if(nprob0.eq.1) call param(m) 
      end if 
c!====================================================================
c!   Maximum L1-norm distance between points satisfying box constraints 
c!====================================================================
      dbox=0.0d+00
      do i=1,n
       dbox=dbox+abs(xu(i)-xl(i))
      end do
c!====================================================================
c!   Maximum distance between one component inside box constraints 
c!====================================================================     
      sm=0.0d+00
      do i=1,n
       smhelp = abs(xu(i)-xl(i))
       if (smhelp > sm) then
          sm = smhelp
       end if
      end do  
c!====================================================================
c!   Tolerances to compare points in 'xp' and 'xb' tables
c!====================================================================    
      dev1 = 1.0d-2
      IF (sm>100d+00) THEN 
         dev2 = 0.2d+00
      ELSE   
         dev2 = 0.1d+00
      END IF       
c!====================================================================
c!     Initialization of the method and a starting point
c!====================================================================
      do istart=1,20              ! New starting point is looked through
       print 49,istart 
 49    format(i20)
       nesc=0                     ! initialization of the number of escapes from a local solution
       np=1                       ! Initialization of the number of points in point list 'xp' (first one is added just below)
       np2=0                      ! initialization of the number of points in local solution list 'xb'
       npall = 0                  ! Initialization of the number of points during the method
       nf1=0                      ! Initialization of the function evaluation counter for f1
       nf2=0                      ! Initialization of the function evaluation counter for f2
       ngrad1=0                   ! Initialization of the subgradient evaluation counter for f1
       ngrad2=0                   ! Initialization of the subgradient evaluation counter for f2
       TEMP=TEMP_in               ! Initial temperature is set
       ntemp = 0                  ! The number of temperature updates
       call cpu_time(time1)       ! Start for CPU time 
       do i=1,n                   
         x(i)=xstart(i,istart)    ! Starting point is stored into the current point 'x_k'
         xbest(i)=x(i)            ! Starting point is stored to be the best known solution
         xp(i,np)=x(i)            ! Starting point is stored to be the first solution in table 'xp' (point is a column)
       end do
       call func1(x,xl,xu,m,f1)           
       call func2(x,m,f2)
       finit=f1-f2                ! Objective value at the starting point
       fbest=finit                ! Best objective value so far 
       print 44,fbest
c!====================================================================
c!   Local search phase 
c!====================================================================
       nbundle= nbundle1              ! The bundle size for AGGSUB
       call aggsub(x,xl,xu,xb,fxb,fbest,m,nbundle,fvalue)  ! AGGSUB method executed starting from 'x'

       np2=np2+1                      ! First local solution obtained for 'xb'
       do i=1,m
        xb(i,np2)=x(i)                ! First local solution added to the local solution list 'xb' (column is a solution)
       end do
       fxb(np2)= fvalue               ! Objective value of the local solution is added to the list 'fxb'

       print 44,fvalue
       print*
       
       if(fvalue.lt.fbest) then       ! Is the new value better than the current best known one?
        fbest=fvalue                  ! The best objective value is updated
        do i=1,n
          xbest(i)=x(i)               ! The best solution vector is updated
        end do
       end if
       istep=1                        ! The initialization of the iteration round counter
       fb(istep)=fbest                ! The first element added to the list of best objective values 'bf'
       
c!====================================================================
c!   Global search phase
c!====================================================================
 3     PRINT*, 'temperature', TEMP
       call globsub(x,xl,xu,xp,m)          ! New auxiliary points are calculated
       istep=istep+1                       ! Iteration round counter updated
       
       if(np1.eq.0) then
        fb(istep)=fbest
        go to 11                           ! No auxiliary points generated with globsub(x) -> 11: improve solution with AUGSUB and STOP
       end if      
       npall=npall+np1                     ! The number of all points updated
       
c!-------------------------------------------------------------
       fbest2=fbest                        ! The current best value of the objective
       fbest1=1.0d+30                      ! The best value of auxiliary points after AGGSUB method
       do i4=1,np1                         ! All auxiliary points are looked through
c! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
        IF (nf1 > 10 000 000) go to 11     ! Have we reached the limit of function value evaluations? -> 11: improve solution with AUGSUB and STOP
        call cpu_time(time2) 
        time3=time2-time1                  ! The used CPU time
        if(time3.gt.tlimit) go to 11       ! Time limit exceeded -> 11: improve solution with AUGSUB and STOP
c! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
        ind2=0                             ! If value is changed to 1 during AGGSUB then no new point is found when AGGSUB is executed
c!------ Starting point start -----        
          i1=np-np1+i4
          do i2=1,m
            x(i2)=xp(i2,i1)                ! Correct auxiliary point picked from the 'xp' as a starting point (point is a column)
          end do
c!------ Starting point end -----   
        nbundle= nbundle1
        call aggsub(x,xl,xu,xb,fxb,fbest,m,nbundle,fvalue)  ! AGGSUB executed for the auxiliary point to obtain a new solution  
c!  -- START: New solution found -- 
        if(ind2.eq.0) then                 ! If .TRUE. a new solution is found
c! ----------------------------------------------------   
c!        NEW LOCAL SOLUTION ?
c! ----------------------------------------------------
         IF (nonb==0) THEN                    ! Current solution is compared to the ones in the table 'xb'
           nxbmax1=min(nxbmax,np2)            ! Quarantees that the correct number of points from 'xb' is looked through
           ii = np2-nxbmax1+1                 ! The first point looked  
           DO j1 = ii, np2                    ! We check is the new local solution different from the previous ones
             diff1 = ABS(fvalue-fxb(j1))      ! Difference between objective function values
             IF (diff1<dev1) THEN             ! Does the new local solution have the same objective value as the previous one?
               DO j2 = 1, m 
                 diff2 = ABS(x(j2)-xb(j2,j1)) ! Difference between one of the components of local solutions
                 IF (diff2 > dev2) then       ! New local solution differs from the previous one
                    go to 13                  ! We move on to compare with the next point in 'xb'
                 END IF
                 IF (j2==m) then              ! New local solution is the same as the previous one
                    go to 14                  ! We move on to look a new auxiliary point
                 END IF
               END DO             
             END IF
  13       END DO
         END IF   
c!  ----------------------------------------------------  
         np2=np2+1                         ! One new local solution is added to the table 'xb'
       
         do k=1,m
          xb(k,np2)=x(k)                   ! New solution added to the table 'xb' (column is a solution)
         end do
         fxb(np2)=fvalue                   ! The objective function value of the solution is added to the list 'fxb'         
         if(fvalue.lt.fbest1) then         ! Has the new solution better value than the current best one of the auxiliary points?
          fbest1=fvalue                    ! The best objective value of the auxiliary points is updated
          do i=1,m
           xbest1(i) = x(i)                ! The best solution vector of the auxiliary points is updated
          end do
         end if
         if(fvalue.lt.fbest) then          ! Has the new solution better value than the current best one?
          fbest=fvalue                     ! The best objective value is updated
          do i=1,m
           xbest(i) = x(i)                 ! The best solution vector is updated
          end do
          nesc=nesc+1                                ! The counter for escaping local solutions is updated by one
          dif1=(fbest2-fbest)/(abs(fbest)+1.0d+00)   ! Relative error between the new best and previous best objective function value    
          print 44,fvalue  
          if(dif1.ge.1.0d-03) go to 12               ! Is relative error big enough? (In this case, there is a clear improvement in the objective function value)
         end if
        end if
c!  -- END: New solution found --   
  14   end do
  12   continue  
       fb(istep)=fbest                                          ! The current best objective function value added to the list 'fb' (i.e. the best result during this iteration round) 
       if(istep.gt.imax) then
        dif=(fb(istep-imax)-fb(istep))/(abs(fb(istep))+1.0d+00) ! Relative error between two consecutive objective function values
        if(dif.le.1.0d-03) go to 11                             ! Is the relative error small enough? -> 11: improve solution with AUGSUB and STOP
       end if
       do i=1,m
        x(i)=xbest1(i)                               ! Update of the solution vector 'x_k'
       end do
c! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
        IF (nf1 > 10 000 000) go to 11               ! Have we reached the limit of objective function value evaluations? -> 11: improve solution with AUGSUB and STOP
        call cpu_time(time2) 
        time3=time2-time1                            ! Used CPU time
        if(time3.gt.tlimit) go to 11                 ! Time limit exceeded -> 11: improve solution with AUGSUB and STOP
c! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
      
       go to 3                                       ! One iteration round finished -> Go back to globsub()      

c!====================================================================
c!   Improving the final result with AUGSUB method
c!====================================================================        
  11   continue  
  
       call cpu_time(time4)                         ! End for the CPU time 
       time5=time4-time1                            ! CPU time for the problem solved so far
       
       PRINT*, '------------------------------------------------'      
       PRINT*, 'Before Augsub' 
       PRINT*, 'nf1=',nf1, 'nf2=', nf2
       PRINT*, 'ngrad1=',ngrad1, 'ngrad2=', ngrad2
       PRINT*, 'CPU=',time5
       PRINT*, '------------------------------------------------'

       nbundle= nbundle2                            ! The bundle size for AUGSUB
       call augsub(xbest,xl,xu,m,nbundle,fbest)
       print 44,fbest
       call cpu_time(time4)                         ! End for the CPU time 
       time5=time4-time1                            ! CPU time for the problem solved (final)

       PRINT*, '------------------------------------------------'      
       PRINT*, 'After Augsub' 
       PRINT*, 'nf1=',nf1, 'nf2=', nf2
       PRINT*, 'ngrad1=',ngrad1, 'ngrad2=', ngrad2
       PRINT*, 'CPU=',time5
       PRINT*, '------------------------------------------------'
       Print* 
c!====================================================================
c!   Printing the result
c!====================================================================    
       SELECT CASE(ngroup)    
c!------------------------------------
       CASE(1)  ! Group 1
         WRITE(42,51) nprob1,n,finit,fbest,time5,nf1,nf2,ngrad1,ngrad2
     1    ,nesc,npall,ntemp
c!------------------------------------
       CASE(2)  ! Group 2
         WRITE(42,52) nprob1,nprob2,n,finit,fbest,time5,nf1,nf2,ngrad1
     1  ,ngrad2,nesc,npall,ntemp
c!------------------------------------
       CASE(3)  ! Group 3
         WRITE(42,51) nprob1,n,finit,fbest,time5,nf1,nf2,ngrad1,ngrad2
     1    ,nesc,npall,ntemp
       END SELECT
c!====================================================================
      END DO
      
  44  FORMAT(3f14.6)    
  51  FORMAT(i5,i6,f24.6,f15.6,f10.4,4i10,2i7,2i7)
  52  FORMAT(i7,i7,i6,f24.6,f15.6,f10.4,4i10,2i7,2i7)
  
      RETURN
      END
      
c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!                    END: HSA-DC method
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/       


c!====================================================================
c!            The value of the DC component f1
c!====================================================================
      subroutine func1(x,xl,xu,m,f1)
      implicit double precision (a-h,o-z)
      double precision x(m),xl(m),xu(m),fi(m)
      COMMON /cnf1/nf1,/cpen/pen,/cka1/ind
     1 ,/cnselect/nselect
      COMMON /cgr/ngroup            ! The selected group of test problems (1="Group 1", 2="Group 2" and 3="Group 3")
      nf1=nf1+1                     ! One new function value evaluation

c!====================================================================
      SELECT CASE(ngroup)
c!------------------------------------
      CASE(1)  ! Group 1
c!------------------------------------
      SELECT CASE(nselect)
        
        CASE(1)                       ! P_1
        
           R1=ABS(x(1) - 1.0d+00)
           R2=ABS(x(1)) - x(2)
           IF (R2 > 0.0d+00) THEN
               f1 = R1+2.0d+02*R2
             ELSE
               f1 = R1
           END IF

        CASE(2)                       ! P_2

           f1 = ABS(x(1)-1.0d+00)
           f1 = f1+ABS(x(3)-1.0d+00)                
           f1 = f1+4.95d+00*(ABS(x(2)+x(4)-2.0d+00))
           f1 = f1+10.1d+00*(ABS(x(2)-1.0d+00)+ABS(x(4)-1.0d+00))      
           R1= ABS(x(1))-x(2)
           IF (R1 > 0.0d+00) THEN
            f1 = f1+200.0d+00*R1
           END IF 
           R2=ABS(x(3)) - x(4)
           IF (R2 > 0.0d+00) THEN
            f1 = f1+180.0d+00*R2
           END IF  

        CASE(3)                       ! P_3
         
           R1=x(1)**2
           R2=x(2)**2
           R3=abs(x(2))

           f1=ABS(x(1)-1.0d+00)+2.0d+02*DMAX1(0.0d+00, ABS(x(1))-x(2))                 

           a1 = R1+R2+R3
           a2 = x(1) + R1 + R2 + R3-0.5d+00
           a3 = ABS(x(1)-x(2)) + R3 - 1.0d+00
           a4 = x(1) + R1 + R2
        
           f1 = f1+10.0d+00*DMAX1(a1,a2,a3,a4)  

        CASE(4)                       ! P_4

           f1 = 9.0d+00 - 8.0d+00*x(1) -6.0d+00*x(2) - 4.0d+00*x(3)  
           f1 = f1 + 2.0d+00*(ABS(x(1)) + ABS(x(2)) +ABS(x(3)))
           f1 = f1 + 2.0d+00*(2.0d+00*x(1)**2+x(2)**2+x(3)**2)

           a1 = x(1) + x(2) + 2.0d+00*x(3) - 3.0d+00   
           a2 = -x(1)
           a3 = -x(2)
           a4 = -x(3)

           f1 = f1 + 1.0d+01*DMAX1(0.0d+00, a1, a2, a3, a4)

        CASE(5)                       ! P_5

           f1 = 0.0d+00
           DO i = 1,m
             f1 = f1+x(i)**2                             
           END DO  

        CASE(6)                       ! P_6

         f1 = 4.0d+00*abs(x(1))+2.0d+00*abs(x(2))+2.2d+01*abs(x(3))
         f1 = f1-3.3d+01*x(1)+1.6d+01*x(2)-2.4d+01*x(3)
         f1=f1+1.0d+02*
     1    dmax1(0.0d+00,2.0d+00*abs(x(2))-3.0d+00*x(1)-7.0d+00)
         f1 = f1+1.0d+02*dmax1(0.0d+00,abs(x(3))-4.0d+00*x(1)-1.1d+01)

        CASE(7)                       ! P_7

           apu = 0.0d+00
           DO j =1,m
              apu = apu+x(j)/dble(j)
           END DO
        
           apu1 = ABS(apu) 
         
           DO i = 2,m
             apu = 0.0d+00
             DO j = 1,m
                apu = apu+x(j)/dble(i+j-1)
             END DO
             apu=ABS(apu)
             IF (apu > apu1) THEN 
                 apu1 = apu
             END IF
           END DO
           f1 = dble(m)*apu1

        CASE(8)                       ! P_8
               
          DO i = 1,m-1
           a1 = x(i)**4 + x(i+1)**2
           a2 = (2.0d+00-x(i))**2 + (2.0d+00-x(i+1))**2
           a3 = 2.0d+00*EXP(-x(i)+x(i+1))
           fi(i)=dmax1(a1,a2,a3)
          END DO
        
          apu2 = fi(1)
          ind = 1
          DO i =2,m-1
           IF (fi(i) > apu2) THEN 
              apu2 = fi(i)
              ind = i
           END IF
          END DO
          f1 = dble(m-1)*apu2

      END SELECT
c!------------------------------------
      CASE(2)  ! Group 2
c!------------------------------------
      SELECT CASE(nselect)
        
        CASE(1)                       ! P_9
        
           apu = x(1)**2
           DO i = 2,m
             apu1=x(i)**2
             IF (apu < apu1 ) THEN 
                apu = apu1
             END IF
           END DO
           f1 = dble(m+1)*apu

        CASE(2)                       ! P_10
        
           apu = x(1)**2
           DO i = 2,m
             apu1=x(i)**2
             IF (apu < apu1 ) THEN 
                apu = apu1
             END IF
           END DO
           f1 = dble(m+1)*apu

        CASE(3)                       ! P_11
        
           apu = x(1)**2
           DO i = 2,m
             apu1=x(i)**2
             IF (apu < apu1 ) THEN 
                apu = apu1
             END IF
           END DO
           f1 = dble(m+1)*apu

        CASE(4)                       ! P_12
        
           f1 = x(1)**2
           DO i = 2,m
             f1 =  dmax1(f1,x(i)**2)
           END DO 

        CASE(5)                       ! P_13
        
           f1 = x(1)**2
           DO i = 2,m
             f1 =  dmax1(f1,x(i)**2)
           END DO 

        CASE(6)                       ! P_14
        
           f1 = 0.0d+00
           DO i = 1,m
             f1 =  f1+x(i)**2                             
           END DO 

        CASE(7)                       ! P_15
        
          f1=0.0d+00
          DO i=1,m-1
            y = -x(i)-x(i+1)
            z = y+x(i)**2+x(i+1)**2-1.0d+00
            f1= f1+dmax1(y,z)
          END DO  

        CASE(8)                       ! P_16
        
           f1 = 0.0d+00
           f2 = 0.0d+00
           DO i = 1,m
             f1 = f1+ABS(x(i)) 
             f2 = f2+dmax1(2.0d+00*(x(i)**2-x(i)-1.0d+00),0.0d+00)
           END DO
           f1=f1+1.0d+01*f2

        CASE(9)                       ! P_17
        
           f1 = 0.0d+00
           f2 = 0.0d+00
           DO i = 1,m
             f1 = f1+ABS(x(i)) 
             f2 = f2+dmax1(2.0d+00*(x(i)**2-x(i)-1.0d+00),0.0d+00)
           END DO
           f1=f1+1.0d+01*f2
           
        CASE(10)                      ! P_18
        
           f1 = 0.0d+00
           f2 = 0.0d+00
           DO i = 1,m
             f1 = f1+ABS(x(i)) 
             f2 = f2+dmax1(2.0d+00*(x(i)**2-x(i)-1.0d+00),0.0d+00)
           END DO
           f1=f1+1.0d+01*f2

        CASE(11)                      ! P_19
        
           f1 = 0.0d+00
           f2 = 0.0d+00
           DO i = 1,m
             f1 = f1+ABS(x(i)) 
             f2 = f2+dmax1(2.0d+00*(x(i)**2-x(i)-1.0d+00),0.0d+00)
           END DO
           f1=f1+1.0d+01*f2

        CASE(12)                      ! P_20
        
           f1 = 0.0d+00
           f2 = 0.0d+00
           DO i = 1,m
             f1 = f1+ABS(x(i)) 
             f2 = f2+dmax1(2.0d+00*(x(i)**2-x(i)-1.0d+00),0.0d+00)
           END DO
           f1=f1+1.0d+01*f2

        CASE(13)                      ! P_21
        
           f1 = 0.0d+00
           f2 = 0.0d+00
           DO i = 1,m
             f1 = f1+ABS(x(i)) 
             f2 = f2+dmax1(2.0d+00*(x(i)**2-x(i)-1.0d+00),0.0d+00)
           END DO
           f1=f1+1.0d+01*f2

        CASE(14)                      ! P_22
        
           DO i = 1,m-1
               a1 = x(i)**4+x(i+1)**2
               a2 = (2.0d+00-x(i))**2+(2.0d+00-x(i+1))**2
               a3 = 2.0d+00*EXP(-x(i)+x(i+1))
               fi(i)=dmax1(a1,a2,a3)
           END DO
           
           apu1 = fi(1)
           DO i =2,m-1
            apu1=dmax1(apu1,fi(i))
           END DO
           f1 = dble(m-1)*apu1 

        CASE(15)                      ! P_23
        
           DO i = 1,m-1
               a1 = x(i)**4+x(i+1)**2
               a2 = (2.0d+00-x(i))**2+(2.0d+00-x(i+1))**2
               a3 = 2.0d+00*EXP(-x(i)+x(i+1))
               fi(i)=dmax1(a1,a2,a3)
           END DO
           
           apu1 = fi(1)
           DO i =2,m-1
            apu1=dmax1(apu1,fi(i))
           END DO
           f1 = dble(m-1)*apu1 

        CASE(16)                      ! P_24
        
           DO i = 1,m-1
               a1 = x(i)**4+x(i+1)**2
               a2 = (2.0d+00-x(i))**2+(2.0d+00-x(i+1))**2
               a3 = 2.0d+00*EXP(-x(i)+x(i+1))
               fi(i)=dmax1(a1,a2,a3)
           END DO

           apu1 = fi(1)
           DO i =2,m-1
            apu1=dmax1(apu1,fi(i))
           END DO
           f1 = dble(m-1)*apu1 
           
        CASE(17)                      ! P_25
        
          y=0.0d+00
          DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
          END DO  
          y=2.0d+00*y
          f1=dmax1(0.0d+00,y)           
           
        CASE(18)                      ! P_26
        
          y=0.0d+00
          DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
          END DO  
          y=2.0d+00*y
          f1=dmax1(0.0d+00,y)           

        CASE(19)                      ! P_27
        
          y=0.0d+00
          DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
          END DO  
          y=2.0d+00*y
          f1=dmax1(0.0d+00,y)           

        CASE(20)                      ! P_28
        
           f1 = 0.0d+00
           DO i = 1, m-1
              f1 = f1 + x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
           END DO        

      END SELECT
c!------------------------------------
      CASE(3)  ! Group 3      
c!------------------------------------
      SELECT CASE(nselect)
        
        CASE(1)                       ! P_29

         f1=2.5d-01*x(1)**4+1.0d-01*x(1)+5.0d-01*x(2)**2 

        CASE(2)                       ! P_30

          f1=2.5d+01*dble(m)
          do i=1,m
             f1=f1+x(i)**2
          end do

        CASE(3)                       ! P_31

          f1=1.0d+00/6.0d+00+ x(1)**6+4.0d+00*x(1)**2+4.0d+00*x(2)**4
     1       +dabs(x(1))

        CASE(4)                       ! P_32

          f1=0.0d+00
          DO i=2,m
            f1=f1+(x(i)-1.0d+00)**2+x(i-1)**2 +x(i)**2
          end DO

        CASE(5)                       ! P_33 

          f1=2.0d+00*(x(1)**2+x(2)**2)

        CASE(6)                       ! P_34 

          f1=0.0d+00
          DO i=1,m-1
           f1=f1+dmax1(x(i+1)-x(i)+1,x(i)**2)
          end DO  
          f1=2.0d+00*f1

      END SELECT
c!------------------------------------    
      END SELECT
c!==================================================================== 
      fp=0.0d+00
      do i=1,m
       fp=fp+dmax1(0.0d+00,xl(i)-x(i),x(i)-xu(i))
      end do
      f1=f1+pen*fp
c!====================================================================
      return
      end
c!====================================================================
    
   
c!====================================================================
c!            The value of the DC component f2
c!====================================================================
      subroutine func2(x,m,f2)
      implicit double precision (a-h,o-z)
      double precision x(m)
      COMMON /cnf2/nf2,/cka2/ka2,/cnselect/nselect
      COMMON /cgr/ngroup            ! The selected group of test problems (1="Group 1", 2="Group 2" and 3="Group 3")
      nf2=nf2+1                     ! One new function value evaluation

c!====================================================================
      SELECT CASE(ngroup)
c!------------------------------------
      CASE(1)  ! Group 1
c!------------------------------------
      SELECT CASE(nselect)
      
        CASE(1)                       ! P_1
        
             f2 = 1.0d+02*(ABS(x(1))-x(2))
             
        CASE(2)                       ! P_2
        
             f2 = 4.95d+00*ABS(x(2)-x(4))
             f2 = f2+90.0d+00*(ABS(x(3))-x(4))
             f2 = f2+100.0d+00*(ABS(x(1))-x(2)) 

        CASE(3)                       ! P_3

             f2 = 1.0d+01*(x(1)**2+x(2)**2+ABS(x(2)))
             f2 = f2 + 1.0d+02*(ABS(x(1))-x(2)) 

        CASE(4)                       ! P_4
 
             f2 = ABS(x(1)-x(2))+ABS(x(1)-x(3)) 

        CASE(5)                       ! P_5

             f2 = 0.0d+00
             DO i = 2,m
               f2 = f2+ABS(x(i)-x(i-1))
             END DO      

        CASE(6)                       ! P_6

          f2=2.0d+01*(-7.0d+00*x(1)+2.0d+00*abs(x(2))-1.8d+01)

        CASE(7)                       ! P_7

          f2 = 0.0d+00
          DO i=1,m
            apu =0.0d+00
            DO j=1,m
              apu =apu +x(j)/dble(i+j-1)
            end DO
            f2 = f2 +abs(apu)
          end DO 

        CASE(8)                       ! P_8

          f2 = 0.0d+00
          DO i = 1,m-1
           a1 = x(i)**4 + x(i+1)**2
           a2 = (2.0d+00-x(i))**2 + (2.0d+00-x(i+1))**2
           a3 = 2.0d+00*EXP(-x(i)+x(i+1))
           fi=dmax1(a1,a2,a3)
           f2 = f2+fi
          END DO

      END SELECT
c!------------------------------------
      CASE(2)  ! Group 2
c!------------------------------------
      SELECT CASE(nselect)
      
        CASE(1)                       ! P_9
        
          apu = ABS(summa(x,m,1)) 
          DO j = 2,20
            apu1 = ABS(summa(x,m,j))
            IF (apu <= apu1 ) THEN 
              apu = apu1
            END IF
          END DO
          f2 = 2.0d+01*apu

        CASE(2)                       ! P_10
        
          f2 = 0.0d+00
          DO i = 2,m
             f2 = f2+ABS(x(i)-x(i-1))
          END DO      

        CASE(3)                       ! P_11
        
          apu = 0.0d+00
          DO j =1,m
             apu = apu + x(j)/dble(j)
          END DO         
          apu1 = ABS(apu)          
          DO i = 2,m
             apu = 0.0d+00  
             DO j = 1,m
               apu = apu+x(j)/dble(i+j-1)
             END DO
             apu1=dmax1(apu1,ABS(apu))
          END DO         
          f2 = m*apu1 

        CASE(4)                       ! P_12
        
          f2 = 0.0d+00
          DO j = 1,20
             f2 = f2 + ABS(summa(x,m,j))
          END DO   

        CASE(5)                       ! P_13 
        
          f2 = 0.0d+00
          DO i=1,m
           apu =0.0d+00
           DO j=1,m
            apu =apu +x(j)/dble(i+j-1)
           end DO
           f2 = dmax1(f2,abs(apu))
          end DO 
     
        CASE(6)                       ! P_14 
        
          f2 = 0.0d+00
          DO j = 1,20
             f2 = f2 + ABS(summa(x,m,j))
          END DO   

        CASE(7)                       ! P_15 
        
          f2 = 0.0d+00
          DO i = 2,m
             f2 = f2+ABS(x(i)-x(i-1))
          END DO                          

        CASE(8)                       ! P_16 
        
           f2 = 0.0d+00
           DO i = 1,m
             f2 =  f2+x(i)**2                             
           END DO 

        CASE(9)                       ! P_17 
        
          f2 = 0.0d+00
          DO i = 2,m
             f2 = f2+ABS(x(i)-x(i-1))
          END DO                

        CASE(10)                      ! P_18
        
          f2=0.0d+00
          DO i=1,m-1
            y = -x(i)-x(i+1)
            z = y+x(i)**2+x(i+1)**2-1.0d+00
            f2= f2+dmax1(y,z)
          END DO  
          
        CASE(11)                      ! P_19
        
          apu = 0.0d+00
          DO j =1,m
             apu = apu + x(j)/dble(j)
          END DO         
          apu1 = ABS(apu)          
          DO i = 2,m
             apu = 0.0d+00  
             DO j = 1,m
               apu = apu+x(j)/dble(i+j-1)
             END DO
             apu1=dmax1(apu1,ABS(apu))
          END DO         
          f2 = m*apu1           

        CASE(12)                      ! P_20
        
          y=0.0d+00
          DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
          END DO  
          y=2.0d+00*y
          f2=dmax1(0.0d+00,y)

        CASE(13)                      ! P_21 
        
           f2 = 0.0d+00
           DO i = 1, m-1
              f2 = f2 + x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
           END DO  

        CASE(14)                      ! P_22 
        
           apu = x(1)**2
           DO i = 2,m
             apu1=x(i)**2
             IF (apu < apu1 ) THEN 
                apu = apu1
             END IF
           END DO
           f2 = dble(m+1)*apu

        CASE(15)                      ! P_23 
        
           f2 = 0.0d+00
           DO i = 1,m
             f2 =  f2+x(i)**2                             
           END DO 

        CASE(16)                      ! P_24 
        
          apu = 0.0d+00
          DO j =1,m
             apu = apu + x(j)/dble(j)
          END DO         
          apu1 = ABS(apu)          
          DO i = 2,m
             apu = 0.0d+00  
             DO j = 1,m
               apu = apu+x(j)/dble(i+j-1)
             END DO
             apu1=dmax1(apu1,ABS(apu))
          END DO         
          f2 = m*apu1           

        CASE(17)                      ! P_25 
        
          f2 = x(1)**2
          DO j = 2,m
            f2 = dmax1(f2,x(j)**2)
          END DO

        CASE(18)                      ! P_26
        
           apu = ABS(x(1))
           DO i = 2,m
            IF (apu < ABS(x(i))) THEN 
               apu = ABS(x(i))
            END IF
           END DO
           f2 = m*apu  

        CASE(19)                      ! P_27 
        
          f2 = 0.0d+00
          DO i = 2,m
             f2 = f2+ABS(x(i)-x(i-1))
          END DO                          
          
        CASE(20)                      ! P_28
        
          f2 = 0.0d+00
          DO i = 2,m
             f2 = f2+ABS(x(i)-x(i-1))
          END DO

      END SELECT
c!------------------------------------
      CASE(3)  ! Group 3
c!------------------------------------
      SELECT CASE(nselect)

        CASE(1)                       ! P_29

          f2=5.0d-01*x(1)**2
             
        CASE(2)                       ! P_30 
        
          f2=0.0d+00
          DO i=1,m
           f2=f2+dabs(x(i))
          end DO
          f2=1.0d+01*f2

        CASE(3)                       ! P_31 
        
          f2=2.1d+00*x(1)**4+4.0d+00*x(2)**2

        CASE(4)                       ! P_32 
        
          f2=0.0d+00
          DO i=2,m
            f2=f2+abs(x(i-1)+x(i))
          END DO
     
        CASE(5)                       ! P_33 
        
          f2=abs(x(1)+x(2))  

        CASE(6)                       ! P_34 
        
          f2=0.0d+00
          DO i=1,m-1
            f2=f2+x(i)**2+x(i+1)-x(i)+1
          END DO 

      END SELECT
c!------------------------------------    
      END SELECT
c!====================================================================
      return
      end
c!====================================================================


c!====================================================================
c!            The subgradient of the DC component f1
c!====================================================================
      subroutine gradient1(x,xl,xu,m,grad)
      implicit double precision (a-h,o-z)
      double precision x(m),grad(m),xl(m),xu(m)
     1 ,a(5),abs_sign(m),fi(m)
      integer ind1(m)  
      COMMON /cngrad1/ngrad1,/cpen/pen,/cka1/ind
     1 ,/cnselect/nselect 
      COMMON /cgr/ngroup            ! The selected group of test problems (1="Group 1", 2="Group 2" and 3="Group 3")
      ngrad1=ngrad1+1               ! One new subgradient evaluation

c!====================================================================
      SELECT CASE(ngroup)
c!------------------------------------
      CASE(1)  ! Group 1
c!------------------------------------
      SELECT CASE(nselect)
         
         CASE(1)                      ! P_1
         
           grad(1)=0.0d+00
           grad(2)=0.0d+00
           IF(x(1).lt.1.0d+00) THEN
               grad(1) = grad(1)-1.0d+00
             ELSE
               grad(1) = grad(1)+1.0d+00
           END IF
        
           IF((ABS(x(1))-x(2)).ge.0.0d+00) THEN
               IF(x(1).lt.0.0d+00) THEN 
                    grad(1) = grad(1)-2.0d+02
                  ELSE
                    grad(1) = grad(1)+2.0d+02  
               END IF                   
               grad(2) = grad(2)-2.0d+02
           END IF        

         CASE(2)                      ! P_2
         
             grad(1) = 0.0d+00
             grad(2) = 0.0d+00
             grad(3) = 0.0d+00
             grad(4) = 0.0d+00
        
             IF(x(1) <= 1.0d+00) THEN
                 grad(1) = -1.0d+00
               ELSE
                 grad(1) = 1.0d+00
             END IF
        
             IF((ABS(x(1))-x(2)) >= 0.0d+00) THEN
                  IF(x(1) <= 0.0d+00) THEN 
                      grad(1) = grad(1)-200.0d+00
                    ELSE
                      grad(1) = grad(1)+200.0d+00  
                  END IF                    
                  grad(2) = -200.0d+00
             END IF
        
             IF((ABS(x(3))-x(4)) >= 0.0d+00) THEN
                  IF(x(3) <= 0.0d+00) THEN 
                      grad(3) = -180.0d+00
                    ELSE
                      grad(3) = 180.0d+00  
                  END IF                    
                  grad(4) = -180.0d+00
             END IF
        
             IF(x(3) <= 1.0d+00) THEN
                    grad(3) = grad(3)-1.0d+00
                 ELSE
                    grad(3) = grad(3)+1.0d+00
             END IF         
        
             IF(x(2) <= 1.0d+00) THEN
                   grad(2) = grad(2)-10.1d+00  
                 ELSE
                   grad(2) = grad(2)+10.1d+00
             END IF         
        
             IF(x(4) <= 1.0d+00) THEN
                   grad(4) = grad(4)-10.1d+00  
                 ELSE
                   grad(4) = grad(4)+10.1d+00
             END IF 
        
             IF((x(2)+x(4)-2.0d+00) <= 0.0d+00 ) THEN
                   grad(2) = grad(2) - 4.95d+00
                   grad(4) = grad(4) - 4.95d+00
                 ELSE
                   grad(2) = grad(2) + 4.95d+00
                   grad(4) = grad(4) + 4.95d+00
             END IF 

         CASE(3)                      ! P_3

             grad(1) = 0.0d+00                          
             grad(2) = 0.0d+00

             IF((x(1)-1.0d+00) <= 0.0d+00 ) THEN 
                   grad(1) = -1.0d+00
                ELSE
                   grad(1) = 1.0d+00
             END IF 

             R1=x(1)**2
             R2=x(2)**2
             R3=abs(x(2))
        
             a(1) = R1 + R2 + R3
             a(2) = x(1) + R1 + R2 + R3 - 0.5d+00
             a(3) = ABS(x(1) - x(2)) + R2 - 1.0d+00
             a(4) = x(1) + R1 + R2                
        
             ind = 1
             DO i = 2,4
                 IF(a(ind) < a(i)) THEN 
                     ind = i
                 END IF 
             END DO
        
             IF (ind == 1) THEN
                   grad(1) = grad(1) + 20.0d+00*x(1)
                   grad(2) = grad(2) + 20.0d+00*x(2)
                   IF (x(2) <= 0.0d+00) THEN
                       grad(2) = grad(2) - 10.0d+00
                     ELSE
                       grad(2) = grad(2) + 10.0d+00
                   END IF
             END IF
        
             IF (ind == 2) THEN
                grad(1) = grad(1) + 10.0d+00 + 20.0d+00*x(1)
                grad(2) = grad(2) + 20.0d+00*x(2)
                IF (x(2) <= 0.0d+00) THEN
                    grad(2) = grad(2) - 10.0d+00
                  ELSE
                    grad(2) = grad(2) + 10.0d+00
                END IF
             END IF 
        
             IF (ind == 3) THEN
                IF ((x(1) - x(2)) <= 0.0d+00) THEN
                    grad(1) = grad(1)-10.0d+00 
                    grad(2) = grad(2)+10.0d+00
                  ELSE
                    grad(1) = grad(1)+10.0d+00
                    grad(2) = grad(2)-10.0d+00                   
                END IF
                IF (x(2) <= 0.0d+00) THEN
                      grad(2) = grad(2)-10.0d+00                   
                  ELSE
                      grad(2) = grad(2)+10.0d+00                 
                END IF
             END IF 
        
             IF (ind == 4) THEN
                 grad(1) = grad(1) + 10.0d+00 + 20.0d+00*x(1)
                 grad(2) = grad(2) + 20.0d+00*x(2)
             END IF
        
             IF ((ABS(x(1))-x(2)) >= 0.0d+00 ) THEN
                  grad(2) = grad(2) - 200.0d+00
                  IF (x(1) <= 0.0d+00) THEN
                         grad(1) = grad(1) - 200.0d+00
                     ELSE
                         grad(1) = grad(1) + 200.0d+00
                  END IF 
             END IF                      

         CASE(4)                      ! P_4

             grad(1) = -8.0d+00 + 8.0d+00*x(1)
             grad(2) = -6.0d+00 + 4.0d+00*x(2)
             grad(3) = -4.0d+00 + 4.0d+00*x(3)
        
             IF (x(1) <= 0.0d+00) THEN 
                  grad(1) = grad(1) - 2.0d+00
                ELSE
                  grad(1) = grad(1) + 2.0d+00
             END IF                 
        
             IF (x(2) <= 0.0d+00) THEN 
                  grad(2) = grad(2) - 2.0d+00
               ELSE
                  grad(2) = grad(2) + 2.0d+00
             END IF 
        
             IF (x(3) <= 0.0d+00 ) THEN 
                  grad(3) = grad(3) - 2.0d+00
                ELSE
                  grad(3) = grad(3) + 2.0d+00
             END IF                 
        
             a(1) = x(1) + x(2) + 2.0d+00*x(3) - 3.0d+00
             a(2) = -x(1)
             a(3) = -x(2)
             a(4) = -x(3)
             a(5) = 0.0d+00
        
             ind = 5
        
             DO i = 1,4
                IF (a(ind) < a(i)) THEN 
                   ind = i
                END IF 
             END DO
        
             IF (ind == 1) THEN
                grad(1) = grad(1) + 10.0d+00 
                grad(2) = grad(2) + 10.0d+00 
                grad(3) = grad(3) + 20.0d+00 
             END IF
        
             IF (ind == 2) THEN
                 grad(1) = grad(1) - 10.0d+00 
             END IF 
        
             IF (ind == 3) THEN
                 grad(2) = grad(2) - 10.0d+00 
             END IF 
        
             IF (ind == 4) THEN
                 grad(3) = grad(3) - 10.0d+00 
             END IF 

         CASE(5)                      ! P_5

             DO i=1,m
               grad(i) = 2.0d+00*x(i)  
             END DO

         CASE(6)                      ! P_6
             
          grad(1)=-33.0d+00
          grad(2)=16.0d+00
          grad(3)=-24.0d+00
        
          IF(x(1).ge.0.0d+00) THEN
              grad(1)=grad(1)+4.0d+00
            ELSE
              grad(1)=grad(1)-4.0d+00
          END IF
        
          IF(x(2).ge.0.0d+00) THEN
              grad(2)=grad(2)+2.0d+00
            ELSE
              grad(2)=grad(2)-2.0d+00
          END IF
        
          IF(x(3).ge.0.0d+00) THEN
              grad(3)=grad(3)+2.2d+01
            ELSE
              grad(3)=grad(3)-2.2d+01
          END IF
        
          apu =2.0d+00*abs(x(2))-3.0d+00*x(1)-7.0d+00
          
          IF(apu.ge.0.0d+00) THEN
               grad(1)=grad(1)-300.0d+00
               IF(x(2).ge.0.0d+00)THEN
                   grad(2)=grad(2)+200.0d+00
                 ELSE
                   grad(2)=grad(2)-200.0d+00
               END IF
          END IF
        
          apu=abs(x(3))-4.0d+00*x(1)-11.0d+00
          IF(apu.ge.0.0d+00) THEN
               grad(1)=grad(1)-4.0d+02
               IF(x(3).ge.0.0d+00)THEN
                  grad(3)=grad(3)+1.0d+02
                ELSE
                  grad(3)=grad(3)-1.0d+02
               END IF
          END IF        

         CASE(7)                      ! P_7

          apu = 0.0d+00
          DO j =1,m
             apu = apu + x(j)/dble(j)
          END DO
          IF (apu >= 0.0d+00) THEN 
                abs_sign(1) = 1.0d+00
            ELSE   
              abs_sign(1) = -1.0d+00
          END IF        
          apu1 = ABS(apu) 
          ind = 1
        
          DO i = 2,m
           apu = 0.0d+00
           DO j = 1,m
             apu = apu + x(j)/dble(i+j-1)
           END DO
           IF (apu >= 0.0d+00) THEN 
                abs_sign(i) = 1.0d+00
            ELSE   
                abs_sign(i) = -1.0d+00
           END IF             
           apu2=ABS(apu)
           IF (apu2 > apu1) THEN 
               apu1 = apu2
               ind = i
           END IF
          END DO
        
          DO j = 1,m
            grad(j) = abs_sign(ind)*dble(m)/dble(ind+j-1)
          END DO    

         CASE(8)                      ! P_8

            do i=1,m
             grad(i)=0.0d+00
            end do

           f1=-1.0d+30
           DO i = 1,m-1
            a1 = x(i)**4 + x(i+1)**2
            a2 = (2.0d+00-x(i))**2 + (2.0d+00-x(i+1))**2
            a3 = 2.0d+00*EXP(-x(i)+x(i+1))
            ff=dmax1(a1,a2,a3)
        
            IF(a1.eq.ff) ind1(i)=1
            IF(a2.eq.ff) ind1(i)=2
            IF(a3.eq.ff) ind1(i)=3
            if(f1.lt.ff) then
             f1=ff
             ind=i
            end if
           END DO
           ind2=ind1(ind)

           IF(ind2==1) THEN
              grad(ind) = 4.0d+00*x(ind)**3 
              grad(ind+1) = 2.0d+00*x(ind+1)
           END IF
           IF(ind2==2) THEN
              grad(ind) = -2.0d+00*(2.0d+00-x(ind))
              grad(ind+1) = -2.0d+00*(2.0d+00-x(ind+1))
           END IF
           IF(ind2==3) THEN
              grad(ind) = -2.0d+00*EXP(-x(ind)+x(ind+1))
              grad(ind+1) = 2.0d+00*EXP(-x(ind)+x(ind+1))
           END IF
        
           do i=1,m
            grad(i) = dble(m-1)*grad(i)
           end do

      END SELECT
c!------------------------------------
      CASE(2)  ! Group 2
c!------------------------------------
      SELECT CASE(nselect)
         
        CASE(1)                       ! P_9
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
              apu1=x(i)**2
              IF (apu < apu1 ) THEN 
                 apu = apu1
                 ind = i
              END IF
           END DO
           grad(ind) = 2.0d+00*dble(m+1)*x(ind)    

        CASE(2)                       ! P_10
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
              apu1=x(i)**2
              IF (apu < apu1 ) THEN 
                 apu = apu1
                 ind = i
              END IF
           END DO
           grad(ind) = 2.0d+00*dble(m+1)*x(ind)    

        CASE(3)                       ! P_11
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
              apu1=x(i)**2
              IF (apu < apu1 ) THEN 
                 apu = apu1
                 ind = i
              END IF
           END DO
           grad(ind) = 2.0d+00*dble(m+1)*x(ind)    

        CASE(4)                       ! P_12
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
             IF (apu < x(i)**2) THEN 
               apu = x(i)**2
               ind = i
             END IF
           END DO
           grad(ind) = 2.0d+00*x(ind)

        CASE(5)                       ! P_13
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
             IF (apu < x(i)**2) THEN 
               apu = x(i)**2
               ind = i
             END IF
           END DO
           grad(ind) = 2.0d+00*x(ind)           
           
         CASE(6)                      ! P_14
         
           DO i=1,m
             grad(i) = 2.0d+00*x(i)
           END DO 

         CASE(7)                      ! P_15
         
           grad(1) = 0.0d+00
           DO i=1,m-1
              grad(i+1) = 0.0d+00
              y = -x(i)-x(i+1)
              z = y+x(i)**2+x(i+1)**2-1.0d+00
              IF (y >= z) THEN
                  grad(i) = grad(i)-1.0d+00
                  grad(i+1) = -1.0d+00
                ELSE
                  grad(i) = grad(i)-1.0d+00+2.0d+00*x(i)
                  grad(i+1) = -1.0d+00+2.0d+00*x(i+1)
              ENDIF
           END DO   
           
         CASE(8)                      ! P_16 
         
           DO i = 1, m
              IF (x(i) < 0.0d+00) THEN 
                  grad(i) = -1.0d+00
                ELSE
                  grad(i) = 1.0d+00
              END IF 
              y=2.0d+00*(x(i)**2-x(i)-1.0d+00)
              IF (y >= 0.0d+00) THEN
                 grad(i) = grad(i)+2.0d+01*(2.0d+00*x(i)-1.0d+00)
              END IF
           END DO 

         CASE(9)                      ! P_17   
         
           DO i = 1, m
              IF (x(i) < 0.0d+00) THEN 
                  grad(i) = -1.0d+00
                ELSE
                  grad(i) = 1.0d+00
              END IF 
              y=2.0d+00*(x(i)**2-x(i)-1.0d+00)
              IF (y >= 0.0d+00) THEN
                 grad(i) = grad(i)+2.0d+01*(2.0d+00*x(i)-1.0d+00)
              END IF
           END DO 

         CASE(10)                     ! P_18  
         
           DO i = 1, m
              IF (x(i) < 0.0d+00) THEN 
                  grad(i) = -1.0d+00
                ELSE
                  grad(i) = 1.0d+00
              END IF 
              y=2.0d+00*(x(i)**2-x(i)-1.0d+00)
              IF (y >= 0.0d+00) THEN
                 grad(i) = grad(i)+2.0d+01*(2.0d+00*x(i)-1.0d+00)
              END IF
           END DO 

         CASE(11)                     ! P_19  
         
           DO i = 1, m
              IF (x(i) < 0.0d+00) THEN 
                  grad(i) = -1.0d+00
                ELSE
                  grad(i) = 1.0d+00
              END IF 
              y=2.0d+00*(x(i)**2-x(i)-1.0d+00)
              IF (y >= 0.0d+00) THEN
                 grad(i) = grad(i)+2.0d+01*(2.0d+00*x(i)-1.0d+00)
              END IF
           END DO 

         CASE(12)                     ! P_20 
         
           DO i = 1, m
              IF (x(i) < 0.0d+00) THEN 
                  grad(i) = -1.0d+00
                ELSE
                  grad(i) = 1.0d+00
              END IF 
              y=2.0d+00*(x(i)**2-x(i)-1.0d+00)
              IF (y >= 0.0d+00) THEN
                 grad(i) = grad(i)+2.0d+01*(2.0d+00*x(i)-1.0d+00)
              END IF
           END DO 

         CASE(13)                     ! P_21 
         
           DO i = 1, m
              IF (x(i) < 0.0d+00) THEN 
                  grad(i) = -1.0d+00
                ELSE
                  grad(i) = 1.0d+00
              END IF 
              y=2.0d+00*(x(i)**2-x(i)-1.0d+00)
              IF (y >= 0.0d+00) THEN
                 grad(i) = grad(i)+2.0d+01*(2.0d+00*x(i)-1.0d+00)
              END IF
           END DO 

        CASE(14)                      ! P_22
        
           do i=1,m
             grad(i) = 0.0d+00
           end do
           DO i = 1,m-1
            a1 = x(i)**4+x(i+1)**2
            a2 = (2.0d+00-x(i))**2+(2.0d+00-x(i+1))**2
            a3 = 2.0d+00*EXP(-x(i)+x(i+1))
            fi(i)=dmax1(a1,a2,a3)
            if(a1==fi(i)) ind1(i)=1
            if(a2==fi(i)) ind1(i)=2
            if(a3==fi(i)) ind1(i)=3
           END DO
         
           apu1 = fi(1)
           i1 = 1
           DO i=2,m-1
              IF(fi(i) > apu1) THEN 
                 apu1 = fi(i)
                 i1 = i
              END IF
           END DO
         
           if(ind1(i1)==1) then
              grad(i1) = 4.0d+00*x(i1)**3 
              grad(i1+1) = 2.0d+00*x(i1+1)
           end if
           if(ind1(i1)==2) then
              grad(i1) = -2.0d+00*(2.0d+00-x(i1))
              grad(i1+1) = -2.0d+00*(2.0d+00-x(i1+1))
           end if
           if(ind1(i1)==3) then
              grad(i1) = -2.0d+00*EXP(-x(i1)+x(i1+1)) 
              grad(i1+1) = 2.0d+00*EXP(-x(i1)+x(i1+1)) 
           end if
         
           do i=1,m
              grad(i) = dble(m-1)*grad(i)
           end do   

        CASE(15)                      ! P_23
        
           do i=1,m
             grad(i) = 0.0d+00
           end do
           DO i = 1,m-1
            a1 = x(i)**4+x(i+1)**2
            a2 = (2.0d+00-x(i))**2+(2.0d+00-x(i+1))**2
            a3 = 2.0d+00*EXP(-x(i)+x(i+1))
            fi(i)=dmax1(a1,a2,a3)
            if(a1==fi(i)) ind1(i)=1
            if(a2==fi(i)) ind1(i)=2
            if(a3==fi(i)) ind1(i)=3
           END DO
         
           apu1 = fi(1)
           i1 = 1
           DO i=2,m-1
              IF(fi(i) > apu1) THEN 
                 apu1 = fi(i)
                 i1 = i
              END IF
           END DO
         
           if(ind1(i1)==1) then
              grad(i1) = 4.0d+00*x(i1)**3 
              grad(i1+1) = 2.0d+00*x(i1+1)
           end if
           if(ind1(i1)==2) then
              grad(i1) = -2.0d+00*(2.0d+00-x(i1))
              grad(i1+1) = -2.0d+00*(2.0d+00-x(i1+1))
           end if
           if(ind1(i1)==3) then
              grad(i1) = -2.0d+00*EXP(-x(i1)+x(i1+1)) 
              grad(i1+1) = 2.0d+00*EXP(-x(i1)+x(i1+1)) 
           end if
         
           do i=1,m
              grad(i) = dble(m-1)*grad(i)
           end do   
           
        CASE(16)                      ! P_24
        
           do i=1,m
             grad(i) = 0.0d+00
           end do
           DO i = 1,m-1
            a1 = x(i)**4+x(i+1)**2
            a2 = (2.0d+00-x(i))**2+(2.0d+00-x(i+1))**2
            a3 = 2.0d+00*EXP(-x(i)+x(i+1))
            fi(i)=dmax1(a1,a2,a3)
            if(a1==fi(i)) ind1(i)=1
            if(a2==fi(i)) ind1(i)=2
            if(a3==fi(i)) ind1(i)=3
           END DO
         
           apu1 = fi(1)
           i1 = 1
           DO i=2,m-1
              IF(fi(i) > apu1) THEN 
                 apu1 = fi(i)
                 i1 = i
              END IF
           END DO
         
           if(ind1(i1)==1) then
              grad(i1) = 4.0d+00*x(i1)**3 
              grad(i1+1) = 2.0d+00*x(i1+1)
           end if
           if(ind1(i1)==2) then
              grad(i1) = -2.0d+00*(2.0d+00-x(i1))
              grad(i1+1) = -2.0d+00*(2.0d+00-x(i1+1))
           end if
           if(ind1(i1)==3) then
              grad(i1) = -2.0d+00*EXP(-x(i1)+x(i1+1)) 
              grad(i1+1) = 2.0d+00*EXP(-x(i1)+x(i1+1)) 
           end if
         
           do i=1,m
              grad(i) = dble(m-1)*grad(i)
           end do   
           
        CASE(17)                      ! P_25
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           y=0.0d+00
           DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
           END DO  
           y=2.0d+00*y
           if(y.ge.0.0d+00) then
               do i=1,m-1
                grad(i)=grad(i)+4.0d+00*x(i)
                grad(i+1)=grad(i+1)+4.0d+00*(x(i+1)-1.0d+00)+2.0d+00
               end do          
           end if           
           
        CASE(18)                      ! P_26
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           y=0.0d+00
           DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
           END DO  
           y=2.0d+00*y
           if(y.ge.0.0d+00) then
               do i=1,m-1
                grad(i)=grad(i)+4.0d+00*x(i)
                grad(i+1)=grad(i+1)+4.0d+00*(x(i+1)-1.0d+00)+2.0d+00
               end do          
           end if                      
           
        CASE(19)                      ! P_27
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           y=0.0d+00
           DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2+x(i+1)-1.0d+00
           END DO  
           y=2.0d+00*y
           if(y.ge.0.0d+00) then
               do i=1,m-1
                grad(i)=grad(i)+4.0d+00*x(i)
                grad(i+1)=grad(i+1)+4.0d+00*(x(i+1)-1.0d+00)+2.0d+00
               end do          
           end if           

        CASE(20)                      ! P_28
        
           do i=1,m
             grad(i)=0.0d+00
           end do
       
           DO i = 1,m-1
             grad(i)=grad(i)+2.0d+00*x(i) 
             grad(i+1)=grad(i+1)+2.0d+00*(x(i+1)-1.0d+00)+1.0d+00 
           END DO

      END SELECT
c!------------------------------------
      CASE(3)  ! Group 3
c!------------------------------------
      SELECT CASE(nselect)

         CASE(1)                      ! P_29
         
           grad(1)=x(1)**3+1.0d-01
           grad(2)=x(2) 

         CASE(2)                      ! P_30 
         
           DO i=1,m
            grad(i)=2.0d+00*x(i)
           END DO

         CASE(3)                      ! P_31 
         
           grad(1)=6.0d+00*x(1)**5+8.0d+00*x(1)
           grad(2)=1.6d+01*x(2)**3
           IF(x(1).ge.0.0d+00) THEN
              grad(1) = grad(1)+1.0d+00
             ELSE
              grad(1) = grad(1)-1.0d+00 
           END IF

         CASE(4)                      ! P_32 
         
           grad(1)=2.0d+00*x(1)
           DO i=2,m-1
            grad(i)=6.0d+00*x(i)-2.0d+00
           END DO
           grad(m)=4.0d+00*x(m)-2.0d+00

         CASE(5)                      ! P_33 
         
           grad(1)=4.0d+00*x(1)
           grad(2)=4.0d+00*x(2)

         CASE(6)                      ! P_34 
         
           do i=1,m
             grad(i)=0.0d+00
           end do
           DO i=1,m-1
             d1= x(i+1)-x(i)+1
             d2= x(i)**2
             IF(d1.ge.d2) THEN
                grad(i) = grad(i)-2.0d+00
                grad(i+1) = grad(i+1)+2.0d+00
              ELSE
                grad(i) = grad(i)+4.0d+00*x(i) 
             END IF 
           END DO

      END SELECT
c!------------------------------------    
      END SELECT    
c!====================================================================
      do i=1,m
       if(x(i).le.xl(i)) then
        grad(i)=grad(i)-pen
       end if
       if(x(i).ge.xu(i)) then
        grad(i)=grad(i)+pen
       end if
      end do
c!====================================================================
      return
      end
c!====================================================================


c!====================================================================
c!            The subgradient of the DC component f2
c!====================================================================
      subroutine gradient2(x,m,grad)
      implicit double precision (a-h,o-z)
      double precision x(m),grad(m),abs_sign(m),fi(m) 
      COMMON /cngrad2/ngrad2,/cka2/ka2,/cnselect/nselect 
      COMMON /cgr/ngroup            ! The selected group of test problems (1="Group 1", 2="Group 2" and 3="Group 3")
      ngrad2=ngrad2+1               ! One new subgradient evaluation

c!====================================================================
      SELECT CASE(ngroup)
c!------------------------------------
      CASE(1)  ! Group 1
c!------------------------------------
      SELECT CASE(nselect)
      
        CASE(1)                       !P_1
        
          IF(x(1).lt.0.0d+00) THEN
                grad(1) = -1.0d+02 
              ELSE
                grad(1) = 1.0d+02
          END IF
          grad(2) = -1.0d+02      


        CASE(2)                       ! P_2 
        
            IF(x(1) <= 0.0d+00) THEN
                grad(1) = -100.0d+00 
              ELSE
                grad(1) = 100.0d+00
            END IF
        
            grad(2) = -100.0d+00
        
            IF(x(3) <= 0.0d+00) THEN
                 grad(3) = -90.0d+00 
              ELSE
                 grad(3) = 90.0d+00
            END IF
        
            grad(4) = -90.0d+00             
        
            IF((x(2)-x(4)) <= 0.0d+00) THEN
                 grad(2) = grad(2) - 4.95d+00
                 grad(4) = grad(4) + 4.95d+00
              ELSE
                 grad(2) = grad(2) + 4.95d+00 
                 grad(4) = grad(4) - 4.95d+00                 
            END IF  

         CASE(3)                      ! P_3

            grad(1) = 0.0d+00
            grad(2) = 0.0d+00
        
            IF (x(2) <= 0.0d+00) THEN 
                 grad(2) = grad(2) + 20.0d+00*x(2) - 10.0d+00 
               ELSE
                 grad(2) = grad(2) + 20.0d+00*x(2) + 10.0d+00 
            END IF
        
            grad(2) = grad(2) - 100.0d+00
            grad(1) = grad(1) + 20.0d+00*x(1)
        
            IF (x(1) <= 0.0d+00 ) THEN 
                 grad(1) = grad(1) - 100.0d+00
               ELSE
                 grad(1) = grad(1) + 100.0d+00
            END IF

          CASE(4)                     ! P_4

            do i=1,3
               grad(i) = 0.0d+00
            end do               
            IF((x(1)-x(2)) <= 0.0d+00) THEN 
               grad(1) = grad(1) - 1.0d+00 
               grad(2) = grad(2) + 1.0d+00
             ELSE
               grad(1) = grad(1) + 1.0d+00
               grad(2) = grad(2) - 1.0d+00 
            END IF
        
            IF ((x(1)-x(3)) <= 0.0d+00) THEN 
               grad(1) = grad(1) - 1.0d+00 
               grad(3) = grad(3) + 1.0d+00
             ELSE
               grad(1) = grad(1) + 1.0d+00
               grad(3) = grad(3) - 1.0d+00     
            END IF   

         CASE(5)                      ! P_5
     
             DO i=1,m
                grad(i) = 0.0d+00 
             END DO
             DO i = 2,m
               IF ((x(i) - x(i-1)) <= 0.0d+00 ) THEN
                    grad(i-1) = grad(i-1) + 1.0d+00  
                    grad(i) = grad(i) - 1.0d+00  
                 ELSE
                    grad(i-1) = grad(i-1) - 1.0d+00   
                    grad(i) = grad(i) + 1.0d+00      
               END IF
             END DO   

         CASE(6)                      ! P_6

          grad(1)=-140.0d+00
          grad(2)=0.0d+00
          grad(3)=0.0d+00
          IF(x(2).ge.0.0d+00) THEN
               grad(2)=40.0d+00
             ELSE
               grad(2)=-40.0d+00
          END IF

         CASE(7)                      ! P_7

          DO i=1,m
            grad(i)=0.0d+00
          END DO  
        
          DO i=1,m
           apu =0.0d+00
           DO j=1,m
             apu =apu+x(j)/dble(i+j-1)
           END DO
           IF(apu.ge.0.0d+00) THEN
              DO j=1,m
               grad(j)=grad(j)+1.0d+00/dble(i+j-1)
              END DO
            ELSE
              DO j=1,m
                grad(j)=grad(j)-1.0d+00/dble(i+j-1)
              END DO
           END IF
          END DO      

         CASE(8)                      ! P_8

          do i=1,m
            grad(i) = 0.0d+00 
          end do         

          DO i = 1,m-1
           a1 = x(i)**4 + x(i+1)**2
           a2 = (2.0d+00-x(i))**2 + (2.0d+00 -x(i+1))**2
           a3 = 2.0d+00*EXP(-x(i)+x(i+1))
           ff=dmax1(a1,a2,a3)
           if(a1==ff) then
              grad(i) = grad(i)+4.0d+00*x(i)**3 
              grad(i+1) = grad(i+1)+2.0d+00*x(i+1)          
           end if
           if(a2==ff) then
              grad(i) = grad(i)-2.0d+00*(2.0d+00-x(i))
              grad(i+1) = grad(i+1)-2.0d+00*(2.0d+00-x(i+1))       
           end if
           if(a3==ff) then
              grad(i) = grad(i)-2.0d+00*EXP(-x(i)+x(i+1)) 
              grad(i+1) = grad(i+1)+2.0d+00*EXP(-x(i)+x(i+1))      
           end if
          END DO     

      END SELECT
c!------------------------------------
      CASE(2)  ! Group 2
c!------------------------------------
      SELECT CASE(nselect)
      
        CASE(1)                       ! P_9
        
          apu = ABS(summa(x,m,1)) 
          ind = 1

          DO j = 2,20
           apu1=ABS(summa(x,m,j))
           IF (apu <= apu1 ) THEN 
              apu = apu1
              ind = j
           END IF
          END DO

          DO i = 1,m
            grad(i) = (0.05d+00*ind)**dble(i-1)
          END DO             
    
          IF (summa(x,m,ind).le.0.0d+00 ) THEN
             do i=1,m 
                grad(i) = -20.0d+00*grad(i)
             end do   
           ELSE
             do i=1,m 
                grad(i) = 20.0d+00*grad(i)
             end do   
          END IF

        CASE(2)                       ! P_10                        

           do i=1,m
              grad(i) = 0.0d+00    
           end do

           DO i = 2,m
             IF ((x(i) - x(i-1)).le.0.0d+00) THEN
                 grad(i-1) = grad(i-1) + 1.0d+00 
                 grad(i) = grad(i) - 1.0d+00 
               ELSE
                 grad(i-1) = grad(i-1) - 1.0d+00  
                 grad(i) = grad(i) + 1.0d+00     
             END IF
           END DO 

        CASE(3)                       ! P_11
        
           do i=1,m
              grad(i) = 0.0d+00       
           end do
           apu = 0.0d+00
           DO j =1,m
              apu = apu + x(j)/dble(j)
           END DO
           IF (apu >= 0.0d+00) THEN 
                 abs_sign(1) = 1.0d+00
            ELSE   
                 abs_sign(1) = -1.0d+00
           END IF        
           apu1 = ABS(apu) 
           ind = 1
        
           DO i = 2,m
             apu = 0.0d+00
             DO j = 1,m
                apu = apu+x(j)/dble(i+j-1)
             END DO
             IF (apu >= 0.0d+00) THEN 
                 abs_sign(i) = 1.0d+00
              ELSE   
                 abs_sign(i) = -1.0d+00
             END IF
             apu2=ABS(apu)
             IF (apu2 > apu1) THEN 
                apu1 = apu2
                ind = i
             END IF
           END DO
         
           DO j = 1,m
              grad(j) = abs_sign(ind)*dble(m)/dble(ind+j-1)
           END DO 

        CASE(4)                       ! P_12
        
          do i=1,m
            grad(i) = 0.0d+00 
          end do 

          DO j = 1,20
            IF (summa(x,m,j).le.0.0d+00) THEN 
               DO i = 1,m
                  grad(i) = grad(i)-(0.05d+00*j)**(i-1)
               END DO                  
              ELSE 
               DO i = 1,m
                  grad(i) = grad(i)+(0.05d+00*j)**(i-1)
               END DO                   
            END IF
          END DO         

        CASE(5)                       ! P_13 
        
           apu =0.0d+00
           DO j=1,m
              apu=apu+x(j)/dble(j)
           END DO
           fi(1)=apu
           apu1=ABS(apu)  
           ind=1

           DO i=2,m
            apu =0.0d+00
            DO j=1,m
              apu=apu+x(j)/dble(i+j-1)
            END DO
            fi(i)=apu
            apu=ABS(apu)
            IF(apu1.lt.apu) THEN
             apu1=apu
             ind=i
            END IF
           END DO 

           IF(fi(ind).ge.0.0d+00) THEN
                DO j=1,m
                  grad(j)=1.0d+00/dble(ind+j-1)
                END DO
              ELSE
                DO j=1,m
                  grad(j)=-1.0d+00/dble(ind+j-1)
                END DO
            END IF

        CASE(6)                       ! P_14 
        
          do i=1,m
            grad(i) = 0.0d+00 
          end do 

          DO j = 1,20
            IF (summa(x,m,j).le.0.0d+00) THEN 
               DO i = 1,m
                  grad(i) = grad(i)-(0.05d+00*j)**(i-1)
               END DO                  
              ELSE 
               DO i = 1,m
                  grad(i) = grad(i)+(0.05d+00*j)**(i-1)
               END DO                   
            END IF
          END DO         

        CASE(7)                       ! P_15 
        
           do i=1,m
              grad(i) = 0.0d+00    
           end do

           DO i = 2,m
             IF ((x(i) - x(i-1)).le.0.0d+00) THEN
                 grad(i-1) = grad(i-1) + 1.0d+00 
                 grad(i) = grad(i) - 1.0d+00 
               ELSE
                 grad(i-1) = grad(i-1) - 1.0d+00  
                 grad(i) = grad(i) + 1.0d+00     
             END IF
           END DO 

         CASE(8)                      ! P_16 
         
           DO i=1,m
             grad(i) = 2.0d+00*x(i)
           END DO 

         CASE(9)                      ! P_17 
         
           do i=1,m
              grad(i) = 0.0d+00    
           end do

           DO i = 2,m
             IF ((x(i) - x(i-1)).le.0.0d+00) THEN
                 grad(i-1) = grad(i-1) + 1.0d+00 
                 grad(i) = grad(i) - 1.0d+00 
               ELSE
                 grad(i-1) = grad(i-1) - 1.0d+00  
                 grad(i) = grad(i) + 1.0d+00     
             END IF
           END DO 

        CASE(10)                      ! P_18
        
           grad(1) = 0.0d+00
           DO i=1,m-1
              grad(i+1) = 0.0d+00
              y = -x(i)-x(i+1)
              z = y+x(i)**2+x(i+1)**2-1.0d+00
              IF (y >= z) THEN
                  grad(i) = grad(i)-1.0d+00
                  grad(i+1) = -1.0d+00
                ELSE
                  grad(i) = grad(i)-1.0d+00+2.0d+00*x(i)
                  grad(i+1) = -1.0d+00+2.0d+00*x(i+1)
              ENDIF
           END DO    

        CASE(11)                      ! P_19
        
           do i=1,m
              grad(i) = 0.0d+00       
           end do
           apu = 0.0d+00
           DO j =1,m
              apu = apu + x(j)/dble(j)
           END DO
           IF (apu >= 0.0d+00) THEN 
                 abs_sign(1) = 1.0d+00
            ELSE   
                 abs_sign(1) = -1.0d+00
           END IF        
           apu1 = ABS(apu) 
           ind = 1
        
           DO i = 2,m
             apu = 0.0d+00
             DO j = 1,m
                apu = apu+x(j)/dble(i+j-1)
             END DO
             IF (apu >= 0.0d+00) THEN 
                 abs_sign(i) = 1.0d+00
              ELSE   
                 abs_sign(i) = -1.0d+00
             END IF
             apu2=ABS(apu)
             IF (apu2 > apu1) THEN 
                apu1 = apu2
                ind = i
             END IF
           END DO
         
           DO j = 1,m
              grad(j) = abs_sign(ind)*dble(m)/dble(ind+j-1)
           END DO 

        CASE(12)                      ! P_20
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           y=0.0d+00
           DO i=1,m-1
             y=y+x(i)**2+(x(i+1)-1.0d+00)**2-x(i+1)-1.0d+00
           END DO  
           y=2.0d+00*y
           if(y >= 0.0d+00) then
               do i=1,m-1
                grad(i)=grad(i)+4.0d+00*x(i)
                grad(i+1)=grad(i+1)+4.0d+00*(x(i+1)-1.0d+00)+2.0d+00
               end do          
           end if

        CASE(13)                      ! P_21
        
           do i=1,m
             grad(i)=0.0d+00
           end do
       
           DO i = 1,m-1
             grad(i)=grad(i)+2.0d+00*x(i) 
             grad(i+1)=grad(i+1)+2.0d+00*(x(i+1)-1.0d+00)+1.0d+00 
           END DO

        CASE(14)                      ! P_22
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
              apu1=x(i)**2
              IF (apu < apu1 ) THEN 
                 apu = apu1
                 ind = i
              END IF
           END DO
           grad(ind) = 2.0d+00*dble(m+1)*x(ind)   

         CASE(15)                     ! P_23 
         
           DO i=1,m
             grad(i) = 2.0d+00*x(i)
           END DO 

         CASE(16)                     ! P_24
         
           do i=1,m
              grad(i) = 0.0d+00       
           end do
           apu = 0.0d+00
           DO j =1,m
              apu = apu + x(j)/dble(j)
           END DO
           IF (apu >= 0.0d+00) THEN 
                 abs_sign(1) = 1.0d+00
            ELSE   
                 abs_sign(1) = -1.0d+00
           END IF        
           apu1 = ABS(apu) 
           ind = 1
        
           DO i = 2,m
             apu = 0.0d+00
             DO j = 1,m
                apu = apu+x(j)/dble(i+j-1)
             END DO
             IF (apu >= 0.0d+00) THEN 
                 abs_sign(i) = 1.0d+00
              ELSE   
                 abs_sign(i) = -1.0d+00
             END IF
             apu2=ABS(apu)
             IF (apu2 > apu1) THEN 
                apu1 = apu2
                ind = i
             END IF
           END DO
         
           DO j = 1,m
              grad(j) = abs_sign(ind)*dble(m)/dble(ind+j-1)
           END DO 

        CASE(17)                      ! P_25
        
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = x(1)**2
           ind = 1
           DO i = 2,m
             IF (apu < x(i)**2) THEN 
               apu = x(i)**2
               ind = i
             END IF
           END DO
           grad(ind) = 2.0d+00*x(ind)

         CASE(18)                     ! P_26
         
           do i=1,m
             grad(i)=0.0d+00
           end do
           apu = ABS(x(1))
           ind = 1
                
           DO i = 2,m
             IF (apu < ABS(x(i))) THEN 
              apu = ABS(x(i))
              ind = i
             END IF
           END DO

           IF (x(ind).le.0.0d+00) THEN 
               grad(ind) = - dble(m)
             ELSE 
               grad(ind) = dble(m)
           END IF  

        CASE(19)                      ! P_27 
        
           do i=1,m
              grad(i) = 0.0d+00    
           end do

           DO i = 2,m
             IF ((x(i) - x(i-1)).le.0.0d+00) THEN
                 grad(i-1) = grad(i-1) + 1.0d+00 
                 grad(i) = grad(i) - 1.0d+00 
               ELSE
                 grad(i-1) = grad(i-1) - 1.0d+00  
                 grad(i) = grad(i) + 1.0d+00     
             END IF
           END DO 

        CASE(20)                      ! P_28 
        
           do i=1,m
              grad(i) = 0.0d+00    
           end do

           DO i = 2,m
             IF ((x(i) - x(i-1)).le.0.0d+00) THEN
                 grad(i-1) = grad(i-1) + 1.0d+00 
                 grad(i) = grad(i) - 1.0d+00 
               ELSE
                 grad(i-1) = grad(i-1) - 1.0d+00  
                 grad(i) = grad(i) + 1.0d+00     
             END IF
           END DO 

      END SELECT
c!------------------------------------
      CASE(3)  ! Group 3
c!------------------------------------
      SELECT CASE(nselect)

        CASE(1)                       ! P_29
        
          grad(1)=x(1)
          grad(2)=0.0d+00

        CASE(2)                       ! P_30
        
          DO i=1,m
            IF(x(i).ge.0.0d+00) THEN
              grad(i) = 1.0d+01
             ELSE
              grad(i) = -1.0d+01 
            END IF
          END DO
       
        CASE(3)                       ! P_31
        
          grad(1)=8.4d+00*x(1)**3
          grad(2)=8.00d+00*x(2)

        CASE(4)                       ! P_32
        
          DO i=1,m
           grad(i)=0.0d+00
          END DO

          DO i=2,m
           d1=x(i-1)+x(i)
           IF(d1.ge.0.0d+00) THEN
             grad(i-1)=grad(i-1)+1.0d+00
             grad(i)=grad(i)+1.0d+00
            ELSE
             grad(i-1)=grad(i-1)-1.0d+00
             grad(i)=grad(i)-1.0d+00
           END IF
          END DO 
          
        CASE(5)                       ! P_33
    
          d1=x(1)+x(2)
          IF(d1.ge.0.0d+00) THEN
             grad(1)=1.0d+00
             grad(2)=1.0d+00
            ELSE
             grad(1)=-1.0d+00
             grad(2)=-1.0d+00
          END IF
          
        CASE(6)                       ! P_34

          grad(1)=2.0d+00*x(1)-1.0d+00
          DO i=2,m-1
            grad(i)=2.0d+00*x(i)
          END DO
          grad(m)=1.0d+00 

      END SELECT
c!------------------------------------    
      END SELECT    
c!====================================================================
      return
      end
c!====================================================================
      
      
c!====================================================================
c!            Calculation of some parameters (needed in some test problems from Group 2)
c!====================================================================
      subroutine param(m)
      implicit double precision (a-h,o-z)
      double precision tj1(100,20)
      common /ctj1/tj1
      do j=1,20
       TJ=5.0d-02*j
       do i=1,m
        tj1(i,j)=tj**(i-1)
       end do
      end do
c!====================================================================
      return
      end     
c!====================================================================
     
      
c!====================================================================   
c!             Information for some functions (needed in some test problems from Group 2)
c!====================================================================
      DOUBLE PRECISION FUNCTION SUMMA(x,n,j)
       implicit double precision (a-h,o-z)
       double precision x(n),tj1(100,20)
       common /ctj1/tj1
       xi=1.0d+00/dble(n)
       SUMMA=0.0d+00
       do i=1,n
        SUMMA=SUMMA+(x(i)-xi)*tj1(i,j)
       end do
      END FUNCTION
c!====================================================================
  
  
c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!           START: Aggregate subgradient method (AGGSUB)
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/     

      subroutine aggsub(x,xl,xu,xb,fxb,fbest,n,nbundle,f)
      implicit double precision (a-h,o-z)
      double precision x(n),xl(n),xu(n)
      double precision xb(n,20000),fxb(20000)
      
      call optimum(x,xl,xu,xb,fxb,fbest,n,nbundle)
      call func1(x,xl,xu,n,f1)
      call func2(x,n,f2)
      f=f1-f2
      
      return
      end
c!====================================================================


c!====================================================================
      subroutine optimum(x,xl,xu,xb,fxb,fbest,m,nbundle)
      implicit double precision (a-h,o-z)
      PARAMETER(maxit=1000)
      double precision fvalues(maxit)
      double precision x(m),x1(m),g(m),v(m)
     1 ,v1(m),grad2(m),xl(m),xu(m)
      double precision xb(m,20000),fxb(20000)
      common /csm/dbox,sm,/cnp2/np2,/cind2/ind2
      COMMON /cnmax/nxpmax,nxbmax           ! Only 'nxpmax' ('nxbmax') number of the newest points used from 'xp' ('xb') table in comparisons.
      COMMON /cnonb/nonb                    ! If nonb = 1, then no previous solutions are used from the table 'xb'. In addition, 'nxbmax' is not used.
                                            ! If nonb = 0, then the previous solutions are used from the table 'xb' (NOTE: only 'nxbmax' newest ones are used).      
      COMMON /caggsub/maxaggsubit           ! If no improvement during 'maxaggsubit' consecutive iterations, then AGGSUB is terminated.
      COMMON /caggsubtol/aggsub_stop_tol    ! The tolerance used in the improvement check in 'maxaggsubit' consecutive iterations during AGGSUB.                                         
c!====================================================================
     
      dist1=1.0d-07
      c1=-2.0d-01
      eps0=1.0d-07
      div=2.0d-01
      tau0=1.0d+00
      taumin=1.0d-06*tau0
      maxiter=maxit
      sdif=1.0d-08
      mturn=4
      niter=0
      
      dev1 = 1.0d-2
      IF (sm>100d+00) THEN 
         dev2 = 0.2d+00
      ELSE   
         dev2 = 0.1d+00
      END IF     
c!====================================================================
      call func1(x,xl,xu,m,fmax)
      call func2(x,m,fmin)
      f2=fmax-fmin                       ! Objective function value at the starting point
      nxchange = 0                       ! Solution x is still the same as the starting point
c!====================================================================
      IF (nonb==0) THEN                  ! Current solution is compared to the ones in the table 'xb'
       nxbmax1=min(nxbmax,np2)           ! Quarantees that the correct number of points from 'xb' is looked through
       istart = np2-nxbmax1+1            ! The index of the first point looked 
       do k=istart,np2                   ! Specific local solutions are looked through
        diff = ABS(f2-fxb(k))            ! Difference between objective function values
        IF (diff<dev1) THEN              ! Does the starting point have the same objective value as a local solution?
          DO i = 1, m   
            difx = ABS(x(i)-xb(i,k))     ! Difference between one of the components of the starting point and local solution
            IF (difx > dev2) then        ! Does the starting point differ from the local solution?
                go to 13                 ! Starting point differs from the local solution -> 13: Move on to compare with the next point in 'xb'
            END IF
            IF (i==m) then               ! Starting point is same as a local solution
               ind2=1                    ! This labels the case that the starting point 'x' is too close to a known local solution
               return
            END IF
          END DO              
        END IF     
  13   end do
      END IF
c!====================================================================
      tau=tau0/div
  1   tau=div*tau                          ! Tau update (i.e. decrease)

c! ----------------------------------------      
c!        Tau STOPPING CONDITION 
c! ----------------------------------------   
      IF(tau.lt.taumin) then               ! If .TRUE. we STOP since tau is too small
         return
      END IF 
      if(tau.le.1.0d-01) then              ! Executed only when tau is smaller than 0.1
c! -------------------------------------------------------------    
c!      COMPARISON with the previous local solutions
c! -------------------------------------------------------------
       IF (nonb==0 .and. nxchange==1) THEN ! The current solution has changed and it is compared to the ones in the table 'xb'
        nxchange = 0                       ! The current solution is now looked through so label for 'nxchange' is put back to 0  
        do k=istart,np2                    ! Specific local solutions are looked through
          diff = ABS(f2-fxb(k))            ! Difference between objective values
          IF (diff<dev1) THEN              ! Does the current point have the same objective value as a local solution?
            DO i = 1, m 
              difx = ABS(x(i)-xb(i,k))     ! Difference between one of the components of current point and local solution
              IF (difx > dev2) then        ! Does the current point differ from the local solution?
                go to 14                   ! Current point differs from the local solution -> 14: Move on to compare with the next point in 'xb'
              END IF
              IF (i==m) then               ! Current solution is the same as a local solution
               ind2=1                      ! This labels the case that current point 'x' is too close to a known local solution
               return
              END IF
            END DO            
          END IF    
  14    end do  
       END IF  
c!--------------------------------------------------------------------           
      end if  
      do i=1,m
       g(i)=1.0d+00/dsqrt(DBLE(m))
      end do
      nnew=0
c!====================================================================
   2  niter=niter+1    

c! ----------------------------------------      
c!        Extra STOPPING CONDITION 
c! ----------------------------------------      
        IF (niter .gt. maxaggsubit) THEN
          diff = fvalues(niter-maxaggsubit)-fvalues(niter-1)
          diff = diff /(abs(fvalues(niter-1))+1.0d+00)
          IF (diff < aggsub_stop_tol .and. fvalues(niter-1)>fbest) THEN 
             RETURN
          END IF 
        END IF    
c! --------------------------------------------      
c!       Iteration limit STOPPING CONDITION 
c! --------------------------------------------   
        IF(niter.gt.maxiter) then
          return
        END IF      
        nnew=nnew+1
        f1=f2
        fvalues(niter)=f1
c!--------------------------------------------------------------------           
        if (nnew.gt.mturn) then
         mturn2=niter-mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.sdif) GO TO 1
        end if
        if (nnew.GE.(2*mturn)) then
         mturn2=niter-2*mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.(1.0d+01*sdif)) GO TO 1
        end if
c!--------------------------------------------------------------------           
        do nsub=1,nbundle
            call subgrad(nsub,x,xl,xu,m,tau,g,v1,grad2)
            dotprod=0.0d+00
            do i=1,m
             dotprod=dotprod+v1(i)*v1(i)
            end do
            r=dsqrt(dotprod)
            IF(r.lt.eps0) GO TO 1
            toler=dmax1(eps0,dist1*dotprod)
            IF(nsub.eq.1) then
             do i=1,m
              v(i)=v1(i)
             end do
            END if
            IF(nsub.gt.1) then
             d1=0.0d+00
             do i=1,m
              d1=d1+(v(i)-v1(i))**2
             end do
             if(d1.gt.(1.0d+01*eps0)) then
                   d3=0.0d+00
                   do i=1,m
                    d3=d3+v1(i)*(v1(i)-v(i))
                   end do
                   dlambda=d3/d1
                   if(dlambda.lt.0.0d+00) dlambda=0.0d+00
                   if(dlambda.gt.1.0d+00) dlambda=1.0d+00
                   do i=1,m
                    v(i)=v1(i)+dlambda*(v(i)-v1(i))
                   end do
                else
                   do i=1,m
                    v(i)=v1(i)
                   end do
             end if
            END if
c!====================================================================
            r=0.d+00
            do i=1,m
             r=r+v(i)*v(i)
            end do
            r=dsqrt(r)
            if((nsub.gt.1).and.(r.lt.toler)) GO TO 1
c!====================================================================
            do i=1,m
             g(i)=-v(i)/r
             x1(i)=x(i)+tau*g(i)
            end do
c!====================================================================
            call func1(x1,xl,xu,m,fmax4)
            call func2(x1,m,fmin4)
            f4=fmax4-fmin4
            f3=(f4-f1)/tau
            decreas=c1*r
            if(f3.le.decreas) then !line search
                call armijo(x,xl,xu,m,g,f1,f5,f4,tau,sigma,r)
                f2=f5
                do i=1,m
                  x(i)=x(i)+sigma*g(i)
                end do
                tau=1.2d+00*tau
                nxchange = 1               ! Solution x is changed
                GO TO 2
             end if
        END do
c!====================================================================
      go to 1
      return
      end
c!====================================================================


c!====================================================================
c! Subroutine 'subgrad' calculates subgradients needed in AGGSUB
c!====================================================================
      subroutine subgrad(nsub,x,xl,xu,n,tau,g,v,grad2)
      implicit double precision (a-h,o-z)
      double precision x1(n),g(n),x(n),grad1(n)
     1 ,grad2(n),v(n),xl(n),xu(n)
     
      do k=1,n
        x1(k)=x(k)+tau*g(k)
      end do

      call gradient1(x1,xl,xu,n,grad1)
      if(nsub.eq.1) call gradient2(x,n,grad2)

      do i=1,n
       v(i)=grad1(i)-grad2(i)
      end do
      
      return
      end
c!====================================================================


c!====================================================================
c! Line search (Armijo-type) used in Step 5 of AGGSUB
c!====================================================================
      subroutine armijo(x,xl,xu,n,g,f1,f5,f4,tau,sigma,r)
      implicit double precision (a-h,o-z)
      double precision x(n),g(n),x1(n),xl(n),xu(n)
      
      sigma=tau
      f5=f4
      k=0
      sigma1=tau
  1   k=k+1
      IF(k.gt.20) RETURN
      sigma1=2.0d+00*sigma1
      do i=1,n
       x1(i)=x(i)+sigma1*g(i)
      end do
      call func1(x1,xl,xu,n,fmax)
      call func2(x1,n,fmin)
      f50=fmax-fmin
      f30=f50-f1+1.0d-02*sigma1*r
      if(f30.le.0.0d+00) then
       sigma=sigma1
       f5=f50
       GO TO 1
      end if
c!====================================================================    
      return
      end
c!====================================================================

c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!           END: Aggregate subgradient method (AGGSUB)
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/   


c!====================================================================
c! Global search phase
c!====================================================================
      subroutine globsub(x,xl,xu,xp,n)
      implicit double precision (a-h,o-z)
      double precision x(n),xl(n),xu(n),xp(n,100000)
      
      call optimum3(x,xl,xu,xp,n)
      
      return
      end
c!====================================================================


c!====================================================================
      subroutine optimum3(x,xl,xu,xp,m)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdg=1000)
      double precision w(m,4*m),w1(m,4*m)
      double precision w2(m,4*m),w3(m,4*m)
      double precision x(m),x1(m),xl(m),xu(m)
     1 ,v2(m),v(m),xp(m,100000),d(m)
      double precision w2old(m,20000)               ! All different eps-subgradient obtained for f2 during this globsub()
      double precision dirs(m,20000)                ! All different search direction generated durig the current globsub() execution
      double precision rdist(4*m)
      double precision prod(maxdg,maxdg),z(maxdg)   ! Used also in Wolfe
      INTEGER ij(maxdg)                             ! Used also in Wolfe
      common /cij/ij,/cvert/jvertex,/cz/z,/ckmin/kmin,/cnp/np
     1 ,/csm/dbox,sm,/cdeviation/deviation
     2 ,/cnselect/nselect,/cnp1/np1,/cnewbest/newbest,newind
     
      COMMON /cnmax/nxpmax,nxbmax        ! Only 'nxpmax' ('nxbmax') number of the newest points used from 'xp' ('xb') table in comparisons.
      COMMON /cnonp/nonp                 ! If nonp = 1, then table 'xp' contains points only from the current round of globsub(). In addition, 'nxpmax' is not used.
                                         ! If nonp = 0, then table 'xp' contains all points generated during execution but only 'nxpmax' newest ones are used.
      COMMON /cTEMP/TEMP                 ! The current value of the temperature
      COMMON /cntemp/ntemp               ! The number of temperature updates
      COMMON /crttemp/Tdec,Tinc,TEMP_in  ! Decrease and increase parameter for the temperature and initial temperature
      COMMON /cnTEMPupdate/nTEMPupdate   ! Is temperature updated (value>0) or not (value=0)

      COMMON /cnway/nway1,nway2          ! Directions for f1 and f2
      
c!====================================================================

      if(m.lt.50) then        ! Parameters when dimension < 50
        eps=1.0d-04*dbox
        nstep=15
        nstep0=5
      end if
      if(m.ge.50) then        ! Parameters when dimension >= 50
        eps=3.0d-05*dbox
        nstep=30
        nstep0=10
      end if 
      
      dist1=1.0d-10
      c1=-1.0d-04
      eps0=1.0d-10
      taumin=5.0d-01
      
      np1=0                   ! Initialization of the number of points generated during this globsub() 
      IF (nonp==1) THEN       ! Table 'xp' contains only points generated during this globsub()
         np=0                 ! Table 'xp' is reset
      END IF          
      
c!   Tolerances for the amount of accepted points (used in temperature update)
      IF (m .le. 10) THEN 
         tol_low = 0.5d+00
         tol_up = 0.75d+00
      ELSE IF (m == 50) THEN 
         tol_low = 0.2d+00
         tol_up = 0.4d+00
      ELSE IF (m .ge. 100) THEN 
         tol_low = 0.1d+00
         tol_up = 0.25d+00
      END IF    
   
      IF (m .lt. 100) then 
         deviation = sm*5d-3
      ELSE 
         deviation = sm*5d-3*2.0d+00  
      END IF  
      
      newbest = 0    ! No better point found during globsub(), at least not yet
      newind = 0     ! The index of the new best point
      ndir = 0       ! The number of all directions so far
      
c!====================================================================
      call func1(x,xl,xu,m,fmax)
      call func2(x,m,fmin)
      f1=fmax-fmin            ! Objective funtion value at the starting point 'x'
      faux = f1               
c!--------------------------------------------------------------------
      d1=1.0d+30
      do i=1,m
       d2=dmax1(x(i)-xl(i),xu(i)-x(i))
       d1=dmin1(d1,d2)
      end do
      if(nstep0.gt.0) then 
        step1=2.0/dble(nstep0)
      end if
      step=d1/dble(nstep)     
      nstep1=nstep0+nstep
c!====================================================================  
c!            DC component f1          
c!====================================================================
      tau=1.0d-03
      ind1 = 1
      call esubdif(ind1,nway1,tau,nsub1,x,xl,xu,m,w1)
      dd3=0.0d+00
      do i=1,nsub1
       do j=1,m
        dd3=dd3+ABS(w1(j,i))
       end do
      end do
      dd5=dd3/dble(nsub1)     
      dd5=dd5/dble(m)     
      nsub3=1                  ! The number of different eps-subgradients
      do k=1,m
       w3(k,1)=w1(k,1)
      end do
      do i=2,nsub1
        do j = 1,nsub3
            do k = 1,m
              dd4 = ABS(w1(k,i)-w3(k,j))
              if (dd4 > 1.0d-02*dd5) then 
                go to 9        ! Different subgradient than the one compared to in w3
              end if 
              if (k==m) then           
                go to 10       ! Same subgradient than already in w3
              end if
            end do            
 9      end do  
        nsub3=nsub3+1          ! One more different eps-subgradient
        do k=1,m
         w3(k,nsub3)=w1(k,i)
        end do         
 10   end do      

c!---- different eps-subgradients ----
      nsub1=nsub3
      do i=1,nsub1
       do k=1,m
        w1(k,i)=w3(k,i)        
       end do
      end do
      
      nall = 0        ! The initialization of the number of the points generated during this globsub()        
      
c!--------------------------------------------------------------------
c!    START: LOOP for different TAU values
c!--------------------------------------------------------------------
      do istep=1,nstep1
       if(istep.le.nstep0) then
        tau=taumin+(istep-1)*step1
        tau1=tau
       end if
       if(istep.gt.nstep0) then
        tau=tau1+(istep-nstep0)*step
       end if
c!====================================================================  
c!            DC component f2          
c!====================================================================
       ind2 = 2
       call esubdif(ind2,nway2,tau,nsub2,x,xl,xu,m,w2)
       dd3=0.0d+00
       do i=1,nsub2
        do j=1,m
         dd3=dd3+ABS(w2(j,i))
        end do
       end do
       dd5=dd3/dble(nsub2)
       dd5=dd5/dble(m)         
       nsub3=1
       do k=1,m
        w3(k,1)=w2(k,1)
       end do      
       do i=2,nsub2
        do j = 1,nsub3
            do k = 1,m
              dd4 = ABS(w2(k,i)-w3(k,j))
              if (dd4 > 1.0d-02*dd5) then 
                go to 11       ! Different subgradient than the one compared to in w3
              end if 
              if (k==m) then           
                go to 12       ! Same subgradient than already in w3
              end if
            end do            
 11     end do  
        nsub3=nsub3+1          ! One more different eps-subgradient
        do k=1,m
         w3(k,nsub3)=w2(k,i)
        end do         
 12    end do   

c!---- different eps-subgradients ----
       nsub2=nsub3
       do i=1,nsub2
        do k=1,m
         w2(k,i)=w3(k,i)
        end do
       end do
       
c!--------------------------------------------------------------
c!    calculating search directions
c!--------------------------------------------------------------
        nei1 = 0      ! The number of rejected points (simply bad directions)
        nei2 = 0      ! The number of rejected points (too close to already generated points)
        nei3 = 0      ! The number of rejected points (simulated annealing rejections)
        njoo1 = 0     ! The number of accepted points (decrease the objective function value)
        njoo2 = 0     ! The number of accepted points (simulated annealing acceptions)
        ndirold = 0   ! The number of old directions 
        ndirnew = 0   ! The number of new directions
        
       do nsub=1,nsub2
c!      Eps-subgradient of the second DC component     
        do i=1,m
         v2(i)=w2(i,nsub)
        end do       
c!      Has the previous eps-subgradient already yielded the search direction?     
        DO i = 1, ndir
          DO j = 1, m
             dif = ABS(w2old(j,i)-v2(j))
             IF (dif > 1.0d-04) THEN
                go to 124
             END IF
             IF (j == m) then 
                do k = 1, m
                  d(k) = dirs(k,i)
                end do 
                ndirold = ndirold+1   ! One more "previous" direction obtained        
                go to 125
             END IF
          END DO
 124    END DO

        ndirnew = ndirnew+1           ! One new direction   
        do k=1,nsub1
         dnorm=0.0d+00          
         do j=1,m
          v(j)=w1(j,k)-v2(j)
          dnorm=dnorm+v(j)**2
         end do         
         rdist(k)=dsqrt(dnorm)
         IF(rdist(k).lt.eps0) THEN 
            nei1 = nei1 + 1           ! One more bad direction
            GO TO 3                   ! This search direction is not a good one -> 3: Go to look the next search direction
         END IF  
        end do
        rmin=1.0d+30
        rmean=0.0d+00 
        do k=1,nsub1
         r=rdist(k)
         do j=1,m
          w(j,k)=w1(j,k)-v2(j)
         end do
         rmean=rmean+r
         IF(r.lt.rmin) then
            kmin=k
            rmin=r
         END if        
         do i=1,k-1
           prod(i,k)=0.0d+00
           do j=1,m
             prod(i,k)=prod(i,k)+w(j,i)*w(j,k)
           end do
           prod(k,i)=prod(i,k)
         end do
         prod(k,k)=r**2
        end do 
        toler=dmax1(eps0,dist1*rmean/dble(nsub1))
c!====================================================================
        call wolfe(nsub1,prod)
c!      "Search direction" (not normalized or with minus sign)  
        do i=1,m
           v(i)=0.0d+00
           do j=1,jvertex
            v(i)=v(i)+w(i,ij(j))*z(j)
           END do
        END do
c!====================================================================
        r=0.0d+00
        do i=1,m
          r=r+v(i)*v(i)
        end do
        r=dsqrt(r)
        if(r.lt.toler) then 
           nei1 = nei1 + 1        ! One more bad direction
           GO TO 3                ! This search direction is not a good one -> 3: Go to look the next search direction
        end if   
        ndir = ndir + 1           ! One new search direction
        do i1=1,m
          d(i1)=-v(i1)/r          ! New search direction
          dirs(i1,ndir)= d(i1)    ! New search direction stored also into the vector 'dirs'
          w2old(i1,ndir) = v2(i1) ! Corresponding eps-subgradient of the second DC component f2 yielding this direction
        end do
        
 125    continue 
        do i2=1,m
         x1(i2)=x(i2)+tau*d(i2)                ! The new point
         if(x1(i2).gt.xu(i2)) then             ! Outside the upper bound
           x1(i2)=xu(i2)
         end if   
         if(x1(i2).lt.xl(i2)) then             ! Outside the lower bound
           x1(i2)=xl(i2) 
         end if   
        end do
        
        nall = nall+1                     ! One more point generated during this globsub()

        istart = 1
        nxpmax1=max(nxpmax,np1)           ! Quarantees that all points gererated during the current run of globsub() are taken into account
        IF (np>nxpmax1) THEN
           istart = np-nxpmax1+1      
        END IF 

c! ---------------------------------------------------------------------------------- 
c!      COMPARISON with the points already generated during this run of globsub()   
c! ----------------------------------------------------------------------------------    
        do i=istart,np                   ! All 'nxpmax' newest points from the list 'xp' are looked through       
          DO j = 1, m   
            difx = ABS(x1(j)-xp(j,i))    ! Difference between one of the components of the new point and the previous point
            IF (difx > deviation) then   ! Does the new point differ from the previous point?
              go to 13                   ! New point differs from the previous one -> 13: Move on to compare with the next point in 'xp'
            END IF
            IF (j==m) then               ! Is the new point same as the previous one?
              nei2 = nei2 + 1            ! One more point that is the same as previously
              go to 3                    ! The new point is the same as the one obtained previously -> 3: Go to look the next search direction
            END IF
          END DO    
  13    end do
c! ----------------------------------------------------------------    
c!      COMPARISON with the current best objective function value
c! ----------------------------------------------------------------    
        call func1(x1,xl,xu,m,fmax4)  ! f1 at the new point      
        call func2(x1,m,fmin4)        ! f2 at the new point
        f4=fmax4-fmin4                ! Objective function value at the new point
c!      Is the new point giving better objective function value than the current iterate?        
        if(f4.lt.f1) then             
          np=np+1                     ! The number of all points increased by one
          np1=np1+1                   ! The number of points during globsub() increased by one
          DO i = 1,m
             xp(i,np) = x1(i)         ! The new point is added to the point list 'xp'
          END DO 
          IF (f4.lt.faux) THEN                    ! Is the function value at the current auxiliary point better than in the previous one? If no previous auxiliary points then faux=f1.
             newind = np                          ! The index of the new better auxiliary point
             faux = f4                            ! The best auxiliary point during globsub() decreasing the current value of the objective function
             dif=(f1-faux)/(abs(faux)+1.0d+00)    ! Relative error between the current best auxiliary point and the current iterate
             if(dif.le.1.0d-03) then              ! Is the relative error small enough? -> Enough improvement in the auxiliary point 
               newbest = 1                        ! A new better auxiliary point found 
             end if              
          END IF 
          njoo1 = njoo1 +1                        ! One more point accepted due to the better objective function value
          GO TO 3
        end if      

c!-----------------------------------
c!     SIMULATED ANNEALING step
c!-----------------------------------
        call random_number(beta)     ! Random number generated
        F42=F1-F4                    ! The difference between objective function values (this is always <0!)
        F14=dmin1(0.0d+00,F42)
        F14=F14/TEMP
        IF(F14.LE.-1.0D+01) THEN
           P = 0.0D+00
          ELSE
           P = DMIN1(1.0d+00,EXP(F14))         
        END IF
        IF (beta.LE.P) THEN          ! Is the point accepted?       
          np=np+1                    ! The number of all points increased by one
          np1=np1+1                  ! The number of points during globsub() increased by one
          DO i = 1,m
           xp(i,np) = x1(i)          ! The new point is added to the point list 'xp'
          END DO         
          njoo2 = njoo2 +1           ! One more point ACCEPTED for the fixed tau using the acceptance mechanism from the simulated annealing      
        ELSE  
          nei3 = nei3 + 1            ! One more point REJECTED for the fixed tau using the acceptance mechanism from the simulated annealing 
        END IF          
            
  3    end do 
       
      end do
c!-------------------------------------------------------------------
c!    END: LOOP for different TAU values
c!-------------------------------------------------------------------
      
c!---------------------------------
c!     UPDATE of the temperature
c!--------------------------------- 
      IF (nTEMPupdate>0) THEN              ! Update procedure used for the temperature
      
        procent = (1.0d+00*np1)/nall       ! The amount of new auxiliary points during this globsub()     
        IF (procent < tol_low) THEN 
            TEMP = Tinc*TEMP               ! Increase the temperature
        ELSE IF (procent > tol_up) THEN  
            TEMP = Tdec*TEMP               ! Decrease the temperature
        END IF      
        ntemp = ntemp + 1                  ! One more temperature update
      
      END IF 
      
      return
      end
c!====================================================================
      

c!====================================================================
c!  Subroutines 'Wolfe' and 'Equations' solves the quadratic
c!  programming problem, to find a descent direction 
c!  in Step 3, Algorithm 2.
c!====================================================================

c!====================================================================
      subroutine wolfe(ndg,prod)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdg=1000)
      common /w01/a,/cij/ij,/cvert/jvertex,/cz/z,/ckmin/kmin
      INTEGER ij(maxdg)
      double precision z(maxdg),z1(maxdg),a(maxdg,maxdg)
     1 ,prod(maxdg,maxdg)
     
      j9=0
      jmax=1000*ndg
      jvertex=1
      ij(1)=kmin
      z(1)=1.0d+00
c!========================================
c!  To calculate X
c!========================================
 1    r=0.0d+00
      do i=1,jvertex
       do j=1,jvertex
        r=r+z(i)*z(j)*prod(ij(j),ij(i))
       end do
      end do
      IF(ndg.eq.1) RETURN
c!========================================
c!  To calculate <X,P_J> and J
c!========================================
      t0=1.0d+30
      do i=1,ndg
        t1=0.0d+00
        do j=1,jvertex
          t1=t1+z(j)*prod(i,ij(j))
        end do
        if(t1.lt.t0) then
            t0=t1
            kmax=i
        end if
      end do
c!========================================
c!  First stopping criterion
c!========================================
      rm=prod(kmax,kmax)
      do j=1,jvertex
       rm=dmax1(rm,prod(ij(j),ij(j)))
      end do
      r2=r-1.0d-12*rm
      if(t0.gt.r2) RETURN
c!========================================
c!  Second stopping criterion
c!========================================
      do i=1,jvertex
       if(kmax.eq.ij(i)) RETURN
      end do
c!========================================
c! Step 1(e) from Wolfe's algorithm
c!========================================
      jvertex=jvertex+1
      ij(jvertex)=kmax
      z(jvertex)=0.0d+00
c!========================================
 2    do i=1,jvertex
       do j=1,jvertex
        a(j,i)=1.0d+00+prod(ij(j),ij(i))
       end do
      end do
      j9=j9+1
      if(j9.gt.jmax) RETURN
      call equations(jvertex,z1)
      do i=1,jvertex
       if(z1(i).le.1.0d-10) go to 3
      end do      
      do i=1,jvertex
       z(i)=z1(i)
      end do
      go to 1
  3   teta=1.0d+00
      do i=1,jvertex
       z5=z(i)-z1(i)
       if(z5.gt.1.0d-10) teta=dmin1(teta,z(i)/z5)
      end do
      do i=1,jvertex
       z(i)=(1.0d+00-teta)*z(i)+teta*z1(i)
       if(z(i).le.1.0d-10) then
           z(i)=0.0d+00
           kzero=i
       end if
      end do
      j2=0
      do i=1,jvertex
       IF(i.ne.kzero) then
           j2=j2+1
           ij(j2)=ij(i)
           z(j2)=z(i)
       END if
      end do
      jvertex=j2
      go to 2    
      
      return
      end
c!====================================================================


c!====================================================================
      subroutine equations(n,z1)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdg=1000)
      allocatable b(:,:)
      common /w01/a
      double precision a(maxdg,maxdg),z1(maxdg)
      allocate(b(maxdg,maxdg)) 
      
      do i=1,n
       do j=1,n
        b(j,i)=a(j,i)
       end do
       b(n+1,i)=1.0d+00
      end do
      do i=1,n
       r=b(i,i)
       if(r.le.1.0d-10) r=1.0d-10
       do j=i,n+1
        b(j,i)=b(j,i)/r
       end do
       do k=i+1,n+1
        do j=i+1,n
         b(k,j)=b(k,j)-b(k,i)*b(i,j)
        end do
       end do
      end do
      z1(n)=b(n+1,n)
      do i=1,n-1
        k=n-i
        z1(k)=b(n+1,k)
        do j=k+1,n
         z1(k)=z1(k)-b(j,k)*z1(j)
        END do
      end do
      z2=0.0d+00
      do i=1,n
       z2=z2+z1(i)
      end do
      do i=1,n
       z1(i)=z1(i)/z2
      end do
      deallocate(b)
      
      return
      end
c!=====================================================================


c!=====================================================================
c! Calculation of eps-subdifferentials
c!=====================================================================
      subroutine esubdif(ind,nway,tau,nsub,x,xl,xu,m,w)
      implicit double precision (a-h,o-z)
      double precision x(m), v1(m), d(m), x3(m)
     1 , w(m,4*m), xl(m), xu(m)
      COMMON /clim/limit1,limit2    ! the selected numbers of eps-subgradients for f1 and f2      
c!=====================================================================

      IF ((nway < 0) .or.(nway > 4)) THEN  
        nway = 0                           ! The selection of 'nway' was invalid. Default is set. 
      END IF
      
      IF ((nway==0) .or. (nway==2)) THEN  
         nlimit1 = MIN(limit1,2*m)
         nlimit2 = MIN(limit2,2*m)
      ELSE IF (nway==1) THEN 
         nlimit1 = MIN(limit1,m)
         nlimit2 = MIN(limit2,m)
      ELSE IF (nway>3) THEN 
         nlimit1 = MIN(limit1,4*m)
         nlimit2 = MIN(limit2,4*m)       
      END IF 
      
c!====================================================================
      m1=4*m                                           ! At most '4*dimension' eps-subgradients will be generated
      nsub=0                                           ! Initilization of the number of generated eps-subgradients
      do nb=1,m1                                       ! 'm1' subgradients generated at most
       if((nsub.ge.nlimit1).and.(ind.eq.1)) return     ! STOP if already enough eps-subgradients for f1
       if((nsub.ge.nlimit2).and.(ind.eq.2)) return     ! STOP if already enough eps-subgradients for f2
 
c!*******************************************************
c!    Initialization fo the direction used to calculate an eps-subgradient 
       do j=1,m
        d(j)=0.0d+00
       end do

c!------------------------------------------------------------------
c!       DIRECTIONS for f1 and f2 
c!------------------------------------------------------------------
       SELECT CASE (nway)
c!-----------------------------------------------------------------------------------------------------
c!       CASE 0: Unit directions in the order of components (always +1 and -1 for each component)
c!-----------------------------------------------------------------------------------------------------    
         CASE(0)
           n1=nb/2                          ! Floor of the number 'nb/2'
           n2=2*n1
           if(nb.eq.n2) then                ! For even 'nb'
             d(n1)=-1.0d+00                 ! Component value is -1
           end if
           if(nb.ne.n2) then                ! For odd 'nb'
             n3=(nb+1)/2                    
             d(n3)=1.0d+00                  ! Component value is +1
           end if
c!------------------------------------------------------------------------------------------------------------------------------------
c!       CASE 1: Random unit directions in the order of components (always randomly generates whether +1 or -1 for each component)
c!------------------------------------------------------------------------------------------------------------------------------------               
         CASE(1) 
          call random_number(beta)     ! Random number        
          if (beta .lt. 0.5d+00) then  
            d(nb)=1.0d+00              ! Component value is +1
          else  
            d(nb)=-1.0d+00             ! Component value is -1
          end if  
c!-------------------------------------------------------------------------------------------------------------------
c!       CASE 2: Half random unit directions in the order of components and half completely random directions
c!-------------------------------------------------------------------------------------------------------------------         
         CASE(2)
          n1=nb/2                          ! Floor of number 'nb/2'
          n2=2*n1
c!        Random direction      
          if(nb.eq.n2) then                ! For even 'nb'
 2          continue
            do j = 1,m
              call random_number(beta)     ! Random number        
              if (beta .lt. (1.0d+00/3.0d+00)) then 
                d(j) = -1.0d+00
              else if (beta .lt. (2.0d+00/3.0d+00)) then 
                d(j) = 0.0d+00
              else 
                d(j) = 1.0d+00
              end if    
            end do
       
            summa = 0.0d+00
            do j = 1, m
              summa = summa + ABS(d(j))
            end do 
       
            IF (summa < 1.0d+00) THEN 
              go to 2                   ! Zero vector -> New direction needs to be generated
            END IF  
c!-------------------------------------------------         
c! Note: the above direction is not normalized    
c!-------------------------------------------------                     
          end if
          
c!        Random unit direction            
          if(nb.ne.n2) then              ! For odd 'nb'
            n3=(nb+1)/2
            call random_number(beta)     ! Random number        
            if (beta .lt. 0.5d+00) then  
              d(n3)=1.0d+00              ! Component value is +1 
            else  
              d(n3)=-1.0d+00             ! Component value is -1
            end if          
         end if      
         
c!-------------------------------------------------
c!       CASE 3: Completely random directions
c!-------------------------------------------------          
         CASE(3)
 1        continue
          do j = 1,m
            call random_number(beta)     ! Random number        
            if (beta .lt. (0.25d+00)) then 
              d(j) = -1.0d+00
            else if (beta .lt. (0.75d+00)) then 
              d(j) = 0.0d+00
            else 
              d(j) = 1.0d+00
            end if    
          end do
       
          summa = 0.0d+00
          do j = 1, m
            summa = summa + ABS(d(j))
          end do 
       
          IF (summa < 1.0d+00) THEN 
            go to 1
          END IF 
c!-------------------------------------------------         
c! Note: the above direction is not normalized    
c!-------------------------------------------------                   

c!----------------------------------------------------------------------------------------------------        
c!       CASE 4: All unit directions in the order of components + some random ones (if nlimit > 2*m)       
c!----------------------------------------------------------------------------------------------------
         CASE(4)
c!        Unit direction 
          IF (nb <=2*m) THEN

            n1=nb/2                          ! Floor of number 'nb/2'
            n2=2*n1
            if(nb.eq.n2) then                ! For even 'nb'
              d(n1)=-1.0d+00                 ! Component value is -1
            end if
            if(nb.ne.n2) then                ! For odd 'nb'
              n3=(nb+1)/2
              d(n3)=1.0d+00                  ! Component value is +1
            end if

c!        Random direction                
          ELSE
 3          continue
            do j = 1,m
              call random_number(beta)     ! Random number        
              if (beta .lt. (0.25d+00)) then 
                d(j) = -1.0d+00
              else if (beta .lt. (0.75d+00)) then 
                d(j) = 0.0d+00
              else 
                d(j) = 1.0d+00
              end if    
            end do
       
            summa = 0.0d+00
            do j = 1, m
              summa = summa + ABS(d(j))
            end do 
       
            IF (summa < 1.0d+00) THEN  
              go to 3                     ! Zero vector -> New direction needs to be generated
            END IF 
          END IF    
c!-------------------------------------------------         
c! Note: the above direction is not normalized    
c!-------------------------------------------------    
   
       END SELECT
c!*******************************************************
              
       do j=1,m
         x3(j)=x(j)+tau*d(j)             ! Point were eps-subgradient is calculated
       end do
       
       if(ind.eq.1) then                 ! eps-subgradient of f1 computed
           call gradient1(x3,xl,xu,m,v1) ! subgradient of f1 computed at 'x3' (it is an eps-subgradient)
           nsub=nsub+1                   ! the number of eps-subgradients increased by one
           do j=1,m 
             w(j,nsub)=v1(j)             ! New eps-subgradient stored to 'w'
           end do
       end if

       if(ind.eq.2) then            ! eps-subgradient of f2 computed
           call gradient2(x3,m,v1)  ! subgradient of f2 computed at 'x3' (it is an eps-subgradient)
           nsub=nsub+1              ! the number of eps-subgradients increased by one
           do j=1,m
             w(j,nsub)=v1(j)        ! New eps-subgradient stored to 'w'
           end do          
       end if
      end do
c!-------------------------------------------------------------------
      return
      end
c!=====================================================================      


c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!           START: Augmented subgradient method (AUGSUB)
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

      subroutine augsub(x,xl,xu,n,nbundle,f)
      implicit double precision (a-h,o-z)
      double precision x(n),xl(n),xu(n)
      
      call optimum2(x,xl,xu,n,nbundle)
      call func1(x,xl,xu,n,f1)
      call func2(x,n,f2)
      f=f1-f2
      
      return
      end
c!=====================================================================      


c!=====================================================================      
      subroutine optimum2(x,xl,xu,m,nbundle)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdg=1000, maxit=10000)
      double precision x(m),x1(m),xl(m),xu(m)
      double precision g(m),v(m+1),grad2(m)
      double precision prod(maxdg,maxdg),z(maxdg)
      double precision w(m+1,maxdg)
      double precision fvalues(maxit)
      INTEGER ij(maxdg)   
      common /cij/ij,/cvert/jvertex,/cz/z,/ckmin/kmin
      
      COMMON /caugsub/maxaugsubit        ! If no improvement during 'maxaugsubit' consecutive iterations, then AUGSUB is terminated
      COMMON /caugsubtol/augsub_stop_tol ! The tolerance used in the improvement check in 'maxaugsubit' consecutive iterations during AUGSUB                                        
      
c!=====================================================================

      dist1=1.0d-07
      c1=-2.0d-01
      div=2.0d-01
      eps0=1.0d-07
      tau0=1.0d-01
      taumin=1.0d-06*tau0
      maxiter=maxit
      sdif=1.0d-05
      mturn=4
      niter=0
      
c!=====================================================================
      tau=tau0/div
      call func1(x,xl,xu,m,fmax)
      call func2(x,m,fmin)
      f2=fmax-fmin      
        
  1   tau=div*tau
  
c! ----------------------------------------      
c!        Tau STOPPING CONDITION 
c! ----------------------------------------   
      IF(tau.lt.taumin) THEN 
           RETURN
      END IF 
      do i=1,m
       g(i)=1.0d+00/dsqrt(DBLE(m))
      end do
      nnew=0
c!=====================================================================
   2  niter=niter+1  

c! ---------------------------------------------      
c!        Iteration limit STOPPING CONDITION 
c! ---------------------------------------------    
        IF(niter.gt.maxiter) THEN 
           RETURN
        END IF    
        
        nnew=nnew+1
        f1=f2
        fmax1=fmax
        fvalues(niter)=f1
c! ----------------------------------------      
c!        Extra STOPPING CONDITION 
c! ----------------------------------------         
        IF(niter .gt. maxaugsubit) THEN                ! Testing if nothing has changed during 'maxaugsubit' rounds
            fcomp1=fvalues(niter-maxaugsubit) 
            fcomp2=fvalues(niter) 
            IF (ABS(fcomp2-fcomp1)<augsub_stop_tol) THEN 
               RETURN
            END IF 
        END IF      
c!-------------------------------------------------------------------
        if (nnew.gt.mturn) then
         mturn2=niter-mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.sdif) GO TO 1
        end if
        if (nnew.GE.(2*mturn)) then
         mturn2=niter-2*mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.(1.0d+01*sdif)) GO TO 1
        end if
c!-------------------------------------------------------------------
        rold = 0.0d+00
        do nsub=1,nbundle
        
            call subgrad2(nsub,fmax1,x,xl,xu,m,tau,g,v,grad2)           
            
            dotprod=0.0d+00
            do i=1,m+1
             dotprod=dotprod+v(i)*v(i)
            end do
                       
            rnew=dsqrt(dotprod)                          ! New 'r' value
            difr = ABS(rnew-rold)                        ! Difference between the new and previous 'r' value.
            
            IF(rnew.lt.eps0 .or. difr<1.0d-8) GO TO 1    ! We exit if the new 'r' is too small or two consecutive 'r' values are the same      
            rold = rnew                                  ! New value will be stored as an old 'r' value
            
            IF(nsub.eq.1) then
                rmean=rold
                kmin=1
                rmin=rold
            END if
            IF(nsub.gt.1) then
                rmin=dmin1(rmin,rold)
                IF(rold.eq.rmin) kmin=nsub
                rmean=((nsub-1)*rmean+rold)/dble(nsub)
            END if
            toler=dmax1(eps0,dist1*rmean)
            do i=1,nsub-1
             prod(nsub,i)=0.0d+00
             do j=1,m+1
              prod(nsub,i)=prod(nsub,i)+w(j,i)*v(j)
             end do
             prod(i,nsub)=prod(nsub,i)
            end do
            prod(nsub,nsub)=dotprod 
c!=====================================================================
            do i=1,m+1
             w(i,nsub)=v(i)
            end do
            call wolfe(nsub,prod)
c!=====================================================================
            do i=1,m+1
             v(i)=0.0d+00
             do j=1,jvertex
              v(i)=v(i)+w(i,ij(j))*z(j)
             END do
            END do
c!=====================================================================
            r=0.d+00
            do i=1,m+1
             r=r+v(i)*v(i)
            end do
            r=dsqrt(r)
            if(r.lt.toler) GO TO 1
c!=====================================================================
            do i=1,m
              g(i)=-v(i)/r
              x1(i)=x(i)+tau*g(i)
            end do
c!=====================================================================
            call func1(x1,xl,xu,m,fmax4)
            call func2(x1,m,fmin4)
            f4=fmax4-fmin4           
            f3=(f4-f1)/tau
            decreas=c1*r
            
            if(f3.le.decreas) then
                  fmax=fmax4
                  call armijo2(x,xl,xu,m,g,f1,f5,fmax,f4,tau,sigma,r) 
                  f2=f5      
                  do i=1,m
                   x(i)=x(i)+sigma*g(i)
                  end do
                  tau=1.2d+00*tau
                  GO TO 2
            end if
        END do
c!=====================================================================
      go to 1
      return
      end
c!=====================================================================


c!=====================================================================
c! Subroutine subgrad2 calculates subgradients or discrete gradients needed in AUGSUB
c!=====================================================================
      subroutine subgrad2(nsub,fmax1,x,xl,xu,n,tau,g,v,grad2)
      implicit double precision (a-h,o-z)
      double precision x1(n),g(n),x(n),grad1(n)
     1 ,grad2(n),v(n+1),xl(n),xu(n)

      do k=1,n
        x1(k)=x(k)+tau*g(k)
      end do
      
      call gradient1(x1,xl,xu,n,grad1)
      call func1(x1,xl,xu,n,ftau)
      if(nsub.eq.1) call gradient2(x,n,grad2)

      alpha=ftau-fmax1
      do i=1,n
       v(i)=grad1(i)-grad2(i)
       alpha=alpha-tau*grad1(i)*g(i)
      end do
      v(n+1)=alpha
      
      return
      end
c!=====================================================================


c!===========================================================
c! Line search (Armijo-type), Step 5 of AUGSUB
c!===========================================================
      subroutine armijo2(x,xl,xu,n,g,f1,f5,f5max,f4,tau,sigma,r)
      implicit double precision (a-h,o-z)
      double precision x(n),g(n),x1(n),xl(n),xu(n)
      
      sigma=tau
      f5=f4
      k=0
      sigma1=tau
  1   k=k+1
      IF(k.gt.20) RETURN
      sigma1=2.0d+00*sigma1
      do i=1,n
       x1(i)=x(i)+sigma1*g(i)
      end do
      call func1(x1,xl,xu,n,fmax)
      call func2(x1,n,fmin)
      f50=fmax-fmin
      f30=f50-f1+1.0d-02*sigma1*r
      if(f30.le.0.0d+00) then
       sigma=sigma1
       f5=f50
       f5max=fmax
       GO TO 1
      end if
      
      return
      end
!=====================================================================
      
c!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!           END: Augmented subgradient method (AUGSUB)
c!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
