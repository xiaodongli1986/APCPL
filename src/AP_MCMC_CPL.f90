

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
  character(charlen) :: inputMCMCfile, outputMCMCfile, mcmcdir, nowchisqstr, fileindexstr, MCMCfilestr
  type(omwpar) :: nowpar
  logical :: smutabstds_inited 



!---------------------------------------------------------
  !--------------------------------
  ! Settings of the model
  
  de_model_lab  = de_wcdm_lab

  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
!  MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA_post_lensing'
  MCMCfilestr = 'base_w_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'

  iw0col=5; iH0col=37; iomcol=39; maxcol=max(iw0col,iH0col,iomcol,iwacol,iomkcol)+2
  
  num_MCMCchain = 4

  noww0=-1.0; nowwa=0.0; nowomk=0.0    

  ! End of settings
  !--------------------------------

!---------------------------------------------------------
  !--------------------------------
  ! Preparation for the compute of AP likelihood
  
  numomwstds = 1
  omstds(1)  = 0.26_rt;  wstds(1)  = -1.40_rt
!  omstds(2)  = 0.31_rt;  wstds(2)  = -1.00_rt

  print *, '(Begin) Load in necessary files.'
  print *, '* Load in covmats:'
  call load_covmats('')
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
    outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))//'_no4OVER5fact'//&
      trim(adjustl(fileindexstr))
    print *
    print *, '###################################################'
    print *, '** Compuate AP chisqs from file: '
    print *, '   ', trim(adjustl(inputMCMCfile))
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    
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
      ! Values of parameters
      de_CP%Ob0hsq  =  0.02253
      de_CP%h	    =  0.711833E+00
      de_CP%alpha   =  0.142125E+01
      de_CP%beta    =  0.325121E+01  
      de_CP%Odm0    =  nowom - de_CP%Ob0hsq/de_CP%h**2.0
      de_CP%Ok0	    =  0.0
      de_CP%wcdm%w  =  noww0 !! model dependent
      call de_init()
      ! End model dependent
      
      ! DAs & Hzs
      do iz = 1, nz
        DAs(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
        Hs(iz) = 100.0 / de_inv_e(zeffs(iz)) 
        !nowpar%omegam = nowom; nowpar%w=noww0
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
          weightedstds = .false. &
          ) 
        APlnlikes(iline) = sum(chisqs_syscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
      else
        APlnlikes(iline) = 0.0d0
      endif
      
      call cpu_time(t2)
      if (t2-t1.gt.dt) then
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
      nowAPlnlike = APlnlikes(iline)
      nowweight = nowweight * exp(APlnlikemin - nowAPlnlike)
      write(3294,'(7(e13.5,1x))') nowweight, nowlnlike+nowAPlnlike, &
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
