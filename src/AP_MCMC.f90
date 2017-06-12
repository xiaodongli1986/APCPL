

program main
use LSS_ximu_tests
USE de_model_init
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz,num_MCMCchain,numomwstds,nlines
  
  
  integer :: iomcol=1,iw0col=1,iwacol=1,iomkcol=1,iH0col=1,maxcol, ifile,iline
  
  real(rt) :: omegam, omstds(1000), wstds(1000), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    tmpX(1000),nowlnlike,nowAPlnlike,nowweight,nowom,noww0,nowH0,nowwa,nowomk,APlnlikemin,&
    t0,t1,t2,dt
  real(rt), allocatable :: APlnlikes(:), smutabstds(:,:,:,:,:)
  character(charlen) :: inputMCMCfile, outputMCMCfile, mcmcdir, nowchisqstr, fileindexstr, MCMCfilestr, suffixstr=''
  type(omwpar) :: nowpar
  integer,parameter :: model_wcdm=3, model_cpl=4, model_owcdm=5, model_ocpl=6, model_lcdm=7, model_olcdm=8
  integer :: nowmodel
  logical :: smutabstds_inited, debug=.false., avg_counts = .false., print_allinfo=.false.


!  nowmodel = model_wcdm
  nowmodel = model_cpl
!  nowmodel = model_olcdm
  suffixstr = 'base1omws_om0.2600_w-1.0000_nbins26to27'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_ExcludeLastThreeBins_B'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_nbins35to40'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_smax40'
!  suffixstr = 'base1omws_om0.3100_w-1.0000_smax40'
!  suffixstr = 'base1omws_om0.1100_w-2.0000'
!  print_allinfo = .true.
  
!---------------------------------------------------------
  !--------------------------------
  ! Settings of the model
  
  if(nowmodel .eq. model_wcdm) then! .or. nowmodel .eq. model_owcdm) then
	  de_model_lab  = de_wcdm_lab
	  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
          MCMCfilestr = 'base_w_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'	  
          iw0col=5; iH0col=37; iomcol=39; 
  elseif(nowmodel .eq. model_cpl ) then
  	  de_model_lab = de_CPL_lab
  	  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
  	  MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'
  	  iw0col=5; iwacol=6; iH0col=38; iomcol=40; 
  elseif(nowmodel .eq. model_olcdm ) then
  	  de_model_lab = de_LCDM_lab
  	  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
  	  MCMCfilestr = 'base_omegak_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'
  	  iw0col=1; iwacol=1; iH0col=37; iomcol=39; iomkcol=5;
  else
          print *, 'Wrong model : ', nowmodel
          stop
  endif

!  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
  
  ! Special setting for CMB+BAO chain
  !if(.false.) then
  ! mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO/'
  ! MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO'; iw0col=5; iwacol=6; iH0col=36; iomcol=38; 
  !endif
!   MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA_post_lensing'

  maxcol=max(iw0col,iwacol,iH0col,iomcol,iomkcol)+2
  num_MCMCchain = 4

  ! End of settings
  !--------------------------------

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
  ! Loop of all MCMC files
  do ifile = 1, num_MCMCchain

    ! File names
    write(fileindexstr,*) ifile
    fileindexstr = '_'//trim(adjustl(fileindexstr))//'.txt'
    inputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//trim(adjustl(fileindexstr))
    if(avg_counts) fileindexstr = '_avg_counts'//trim(adjustl(fileindexstr)) 
    if(trim(adjustl(suffixstr)).eq.'') &
      suffixstr = trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl(suffixstr))//trim(adjustl(fileindexstr))

    print *
    print *, '###################################################'
    print *, '** Compuate AP chisqs from file: '
    print *, '   ', trim(adjustl(inputMCMCfile))
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    print *, '** outputfile name: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    
    ! Open file and compute likelihoods...
    call de_count_line_number (trim(adjustl(inputMCMCfile)), nlines); allocate(APlnlikes(nlines))
    print *, '** Computing ', nlines, 'chisqs...'
    open(unit=3293,file=inputMCMCfile,action='read')
    iline = 1
    call cpu_time(t0); t1=t0; dt = 60.0;
    ! Loop of all files
    do while(.true.)
      read(3293,*,end=100) tmpX(1:maxcol)
      
      ! Begin model dependent
      nowweight=tmpX(1); nowlnlike=tmpX(2); nowom=tmpX(iomcol+2); noww0=tmpX(iw0col+2); nowH0=tmpX(iH0col+2)
      nowwa=tmpX(iwacol+2); nowomk=tmpX(iomkcol+2)
      
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_owcdm .or. &
      	 nowmodel.eq.model_lcdm .or. nowmodel.eq.model_olcdm) &
      	    nowwa = 0.0; 
      	    
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_cpl   .or. &
         nowmodel.eq.model_lcdm)  &
            nowomk = 0.0; 
            
      if(nowmodel.eq. model_lcdm .or. nowmodel.eq.model_olcdm) &
            noww0 = 0
            
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
          weightedstds = .false., avg_counts = avg_counts &
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
        write(*,'(A,e12.4,1x,f10.3,1x,6(f9.4))') '             Current set of wei / chi2 / APchi2 / par:  ', &
          nowweight, nowlnlike, APlnlikes(iline), nowom, nowH0/100.0, noww0, nowwa, nowomk
        t1=t2
      endif
      iline = iline +1
      cycle
100   exit 
    enddo
    close(3293)
    
    APlnlikemin = minval(APlnlikes(1:nlines))
    
    ! write the values to new file
    print *, '  Write   AP chisqs to file: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    open(unit=3293,file=inputMCMCfile,action='read')
    open(unit=3294,file=outputMCMCfile,action='write')
    iline = 1
    do while(.true.)
      read(3293,*,end=101) tmpX(1:maxcol)
      nowweight=tmpX(1); nowlnlike=tmpX(2); nowom=tmpX(iomcol+2); noww0=tmpX(iw0col+2); nowH0=tmpX(iH0col+2) ! model dependent
      nowwa=tmpX(iwacol+2); nowomk=tmpX(iomkcol+2)
      
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_owcdm .or. &
      	 nowmodel.eq.model_lcdm .or. nowmodel.eq.model_olcdm) &
      	    nowwa = 0.0; 
      	    
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_cpl   .or. &
         nowmodel.eq.model_lcdm)  &
            nowomk = 0.0; 
            
      if(nowmodel.eq. model_lcdm .or. nowmodel.eq.model_olcdm) &
            noww0 = 0
            
      nowAPlnlike = APlnlikes(iline)
      nowweight = nowweight * exp(APlnlikemin - nowAPlnlike)
      write(3294,'(7(e14.7,1x))') nowweight, nowlnlike+nowAPlnlike, &
        nowom, nowH0/100.0, noww0, nowwa, nowomk !model dependent
      iline = iline+1
      cycle
101   exit
    enddo      
    close(3293); close(3294)
    deallocate(APlnlikes)
  enddo     
  
 
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
