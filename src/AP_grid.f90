

program main
use LSS_ximu_tests
USE de_model_init
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz,num_MCMCchain,numomwstds,nlines
  
  
  real(rt) :: omegam, omstds(1000), wstds(1000), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    nowlnlike,nowweight,nowom,noww0,nowH0,nowwa,nowomk,APlnlikemin,&
    t0,t1,t2,dt, ommin, ommax, w0min, w0max, wamin, wamax
  integer :: iline, numom, numw0, numwa, iom, iw0, iwa
  real(rt), allocatable :: APlnlikes(:), smutabstds(:,:,:,:,:)
  character(charlen) :: outputMCMCfile, mcmcdir, nowchisqstr, fileindexstr, MCMCfilestr, suffixstr=''
  type(omwpar) :: nowpar
  integer,parameter :: model_wcdm=3, model_cpl=4, model_owcdm=5, model_ocpl=6, model_lcdm=7, model_olcdm=8
  integer :: nowmodel
  logical :: smutabstds_inited, debug=.false., print_allinfo=.false.


  ommin=0.1;  ommax=0.4; 
  w0min=-1.4; w0max=-0.3; 
  wamin=-2.5; wamax=0.9;
  numom=100; numw0=100; numwa=100  

!  nowmodel = model_wcdm;   wamin=0.0;wamax=0.0;numwa=1
  nowmodel = model_cpl; !ommin=0.31; ommax=0.31; numom=1
!  nowmodel = model_olcdm

 
  !ommin=0.1; ommax=0.2; numom=34; suffixstr = 'base1omws_om0.2600_w-1.0000_1'!_fixom0.31'
  !ommin=0.20303030303d0; ommax=0.3; numom=33; suffixstr = 'base1omws_om0.2600_w-1.0000_2'!_fixom0.31'
  !ommin=0.30303030303d0; ommax=0.4; numom=33; suffixstr = 'base1omws_om0.2600_w-1.0000_3'!_fixom0.31'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_ExcludeLastThreeBins_B'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_nbins35to40'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_smax40'
!  suffixstr = 'base1omws_om0.3100_w-1.0000_smax40'
!  suffixstr = 'base1omws_om0.1100_w-2.0000'
!  print_allinfo = .true.
  

  
!---------------------------------------------------------
  !--------------------------------
  ! Settings of the model

  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/'
  if(nowmodel .eq. model_wcdm) then! .or. nowmodel .eq. model_owcdm) then
  	  de_model_lab = de_wCDM_lab
          MCMCfilestr = 'base_w_AP'	  
  elseif(nowmodel .eq. model_cpl ) then
  	  de_model_lab = de_CPL_lab
  	  MCMCfilestr = 'base_w_wa_AP'
  else
          print *, 'Wrong model : ', nowmodel
          stop
  endif


!---------------------------------------------------------
  !--------------------------------
  ! Preparation for the compute of AP likelihood
  
  numomwstds = 1
!  omstds(1)  = 0.11_rt;  wstds(1)  = -2.0_rt
  omstds(1)  = 0.26_rt;  wstds(1)  = -1.00_rt
  omstds(2)  = 0.26_rt;  wstds(2)  = -0.60_rt
!  omstds(3)  = 0.26_rt;  wstds(3)  = -1.40_rt
!  omstds(4)  = 0.26_rt;  wstds(4)  = -0.60_rt
!  omstds(5)  = 0.31_rt;  wstds(5)  = -1.40_rt
!  omstds(6)  = 0.31_rt;  wstds(6)  = -0.60_rt

  print *, '(Begin) Load in necessary files.'
!  call system('sleep 0'); print *, 'Compute/output covmats...';call calc_covmats();call output_covmats()
  print *, '* Load in covmats:'
  call load_covmats()
  print *, '* Invert covmats:'
  call invert_covmats()
  print *, '* Compute systematic correction:'
  call calc_syscor()
  smutabstds_inited = .false.
  allocate(smutabstds(nbins_database,mubins_database,3,nz,numomwstds))
  ! End 
  !--------------------------------

!---------------------------------------------------------  
  !--------------------------------
  ! Scan of grid
  if(.true.) then

    fileindexstr = '_1.txt'
    if(trim(adjustl(suffixstr)).eq.'') &
      suffixstr = trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl(suffixstr))//trim(adjustl(fileindexstr))

    print *
    print *, '###################################################'
    print *, '** Compuate AP chisqs for: '
    print *, '   ommin, ommax = ', ommin, ommax
    print *, '   w0min, w0max = ', w0min, w0max
    print *, '   wamin, wamax = ', wamin, wamax
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    print *, '** outputfile name: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    
    nlines = numom*numw0*numwa
    allocate(APlnlikes(nlines))
    
    print *, '** Computing ', nlines, 'chisqs...'
    iline = 1
    call cpu_time(t0); t1=t0; dt = 60.0;

    do iom = 1, numom
    do iw0 = 1, numw0
    do iwa = 1, numwa
      
      ! Begin model dependent
      nowom=ommin + (ommax-ommin)/max(dble(numom-1),1.0d0)*(iom-1);
      noww0=w0min + (w0max-w0min)/max(dble(numw0-1),1.0d0)*(iw0-1); 
      nowwa=wamin + (wamax-wamin)/max(dble(numwa-1),1.0d0)*(iwa-1);
      nowH0=70.0;nowomk=0.0
      
      ! Values of parameters
      de_CP%Ob0hsq  =  0.02253
      de_CP%h	    =  nowH0 / 100.0
      de_CP%alpha   =  0.142125E+01
      de_CP%beta    =  0.325121E+01  
      de_CP%Odm0    =  nowom - de_CP%Ob0hsq/de_CP%h**2.0
      de_CP%Ok0	    =  nowomk
      de_CP%wcdm%w  =  noww0 !! model dependent
      de_CP%CPL%w0  =  noww0 !! model dependent
      de_CP%CPL%wa  =  nowwa !! model dependent
      call de_init()
      ! End model dependent
      
      ! DAs & Hzs
      do iz = 1, nz
        DAs(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
        Hs(iz) = 100.0 / de_inv_e(zeffs(iz)) 
        if(debug) then
          nowpar%omegam = 0.06_rt; nowpar%w=-1.5_rt
          DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
          Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
        endif
        !print *, 'Check DA: ', DAz_wcdm(nowpar,zeffs(iz)), DAs(iz)
        !print *, 'Check H:  ', Hz_wcdm(nowpar,zeffs(iz)), Hs(iz)
        !Hs(iz) = Hz_wcdm(nowpar, zeffs(iz))
      enddo
!      stop
      
      ! AP likelihood
      if(.true.) then
        call smu_ximu_CalcDAHChisqs(& 
          DAs, Hs, & ! Values of DA, H in six cosmologies
          omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
          smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
          chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
          chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
          weightedstds = .false., avg_counts = .false. &
          ) 
        APlnlikes(iline) = sum(chisqs_syscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
!        APlnlikes(iline) = sum(chisqs_nosyscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
      else
        APlnlikes(iline) = 0.0d0
      endif

      if(debug) then
        print *, 'chisqs_nosyscor(1,1):', real(chisqs_nosyscor(1,1,:)), real(sum(chisqs_nosyscor(1,1,:)))
        print *, 'chisqs_syscor(1,1):  ', real(chisqs_syscor(1,1,:)), real(sum(chisqs_syscor(1,1,:)))
        print *, 'chisqs_nosyscor_all:', real(chisqs_nosyscor_all)
        print *, 'chisqs_syscor_all:  ', real(chisqs_syscor_all)
        stop
      endif
      
      call cpu_time(t2)
      if (t2-t1.gt.dt.or.print_allinfo) then
        write(*,'(f10.1,A,i5,A,f4.1,A)') (t2-t0)/dt, ' minutes passed.   #-parameters = ', &
           iline, ' (',100*float(iline)/float(nlines),'%)'
        write(*,'(A,1x,6(f9.4))') '             Current set of wei / chi2 / APchi2 / par:  ', &
          APlnlikes(iline), nowom, nowH0/100.0, noww0, nowwa, nowomk
        t1=t2
      endif
      iline = iline +1
    enddo
    enddo
    enddo
    
    APlnlikemin = minval(APlnlikes(1:nlines))
    
    ! write the values to new file
    print *, '  Write   AP chisqs to file: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    open(unit=3294,file=outputMCMCfile,action='write')
    iline = 1
    do iom = 1, numom
    do iw0 = 1, numw0
    do iwa = 1, numwa
      ! Begin model dependent
      nowom=ommin + (ommax-ommin)/max(dble(numom-1),1.0d0)*(iom-1);
      noww0=w0min + (w0max-w0min)/max(dble(numw0-1),1.0d0)*(iw0-1); 
      nowwa=wamin + (wamax-wamin)/max(dble(numwa-1),1.0d0)*(iwa-1);
      nowH0=70.0;nowomk=0.0
      
      nowlnlike = APlnlikes(iline)
      nowweight = exp(APlnlikemin - nowlnlike)
      write(3294,'(7(e14.7,1x))') nowweight, nowlnlike, &
        nowom, nowH0/100.0, noww0, nowwa, nowomk !model dependent
      iline = iline+1
    enddo 
    enddo
    enddo     
    close(3294)
    deallocate(APlnlikes)
  endif     
  
 
!  nowpar%omegam = 0.06_rt; nowpar%w = -1.5_rt

!  if(.true.) &
!   call smu_ximu_CalcDAHChisqs(& 
!    DAs, Hs, & ! Values of DA, H in six cosmologies
!    omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
!    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
!    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
!    chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
!    weightedstds = .true. &
!    ) 
!   do i = 1, N1
!   do j = 1, N2
!     if ( mubins(i) .eq. 25 .and. j.eq.1) then
!       print *, '* mubin / mucut = ', mubins(i),mucuts(j)
!       print *, '  chisqs (no cor) = ', real(chisqs_nosyscor(i,j,1:nz-1)), '; ', real(sum(chisqs_nosyscor(i,j,1:nz-1)))
!       print *, '  chisqs (corred) = ', real(chisqs_syscor(i,j,1:nz-1)), '; ', real(sum(chisqs_syscor(i,j,1:nz-1)))
!     endif
!   enddo
!   enddo
end program
