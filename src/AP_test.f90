program main
use LSS_ximu_tests
implicit none
  
  integer :: i,j,k,i1,i2,iz, nom,nw, imins(10)
  integer, parameter :: numom=71, numw=45, numomwstds = 1
  real(rt) :: omlist(numom), wlist(numw), ommin, ommax, wmin, wmax, omstds(numomwstds), wstds(numomwstds), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    smutabstds(nbins_database,mubins_database,3,nz,numomwstds)
  character(charlen) :: outputdir
  type(omwpar) :: nowpar
  logical :: smutabstds_inited 
  
! Dense1 subscan
!  integer, parameter :: numom=25, numw=51
!  ommin= 0.06_rt; ommax= 0.66_rt;
!  wmin =-2.5_rt; wmax  = 0.0_rt
! Dense2 subscan
!  integer, parameter :: numom=71, numw=45
  ommin= 0.06_rt; ommax= 0.41_rt;
  wmin =-1.5_rt; wmax  = -0.4_rt

!  call check_load_files()
!  stop

  do i = 1, numom
    omlist(i) = ommin + (ommax-ommin)/dble(numom-1.0_rt)*dble(i-1)
  enddo
  do j = 1, numw
    wlist(j) = wmin + (wmax-wmin)/dble(numw-1.0_rt)*dble(j-1)
  enddo

  print *,  'List of omegam to be computed:'
  do i =1, numom
   write(*,'(f8.4)',advance='no'), real(omlist(i)) 
  enddo
  print *
  print *,  'List of w to be computed:'
  do i =1, numw
   write(*,'(f8.4)',advance='no'), real(wlist(i)) 
  enddo
  print *

!  stop
  if(.false.) then
    print *, 'Compute covmats...'
    call calc_covmats()
    print *, 'Output covmats...'
    call output_covmats('')
  endif

  print *, 'Load in covmats...'
  call load_covmats('')
  print *, 'Invert covmats...'
  call invert_covmats()
  print *, 'Compute systematic correction...'
  call calc_syscor()
!  print *, mubins
!  print *, mucuts
!  stop
!  do iz = 2, nz
!    print *, real(dintxi_syscor(1:mubins(1)-1,1,5,iz-1))
!  enddo
!  stop
!  call check_load_files()
  omstds(1)  = 0.26_rt;  wstds(1)  = -1.00_rt
!  omstds(2)  = 0.31_rt;  wstds(2)  = -1.00_rt
!  omstds(3)  = 0.26_rt;  wstds(3)  = -1.40_rt
!  omstds(4)  = 0.26_rt;  wstds(4)  = -0.60_rt
!  omstds(5)  = 0.31_rt;  wstds(5)  = -1.40_rt
!  omstds(6)  = 0.31_rt;  wstds(6)  = -0.60_rt
!  omstds(7)  = 0.31_rt;  wstds(7)  = -1.00_rt
  if(.true.) then
   call smu_ximu_CalcOmWChisqs(&
    omlist=omlist, numom=numom, wlist=wlist, numw=numw, & ! List of omegam, w
    !omlist=omlist(1:1), numom=1, wlist=wlist(1:1), numw=1, & ! List of omegam, w
    outputdir='/home/xiaodongli/LSS/2PCF_AP/chisqs', &
    baseoutputfile='Dense2ScanIOTest_170420_002', & ! Basic name of the outputfile
!    omstd=0.26_rt, wstd=-1.0_rt & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base)
!    omstd=0.31_rt, wstd=-1.0_rt & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base)
!    omstds=0.31_rt, wstd=-0.6_rt & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base)
    omstds = omstds(1:numomwstds), wstds=wstds(1:numomwstds), numomwstds=numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
    weightedstds = .false., avg_counts=.false. &
    )
  endif
  
  smutabstds_inited = .false.
  nowpar%omegam = 0.06_rt; nowpar%w = -1.5_rt
  do iz = 1, nz
    DAs(iz) = DAz_wcdm(nowpar, zeffs(iz))
    Hs(iz) = Hz_wcdm(nowpar, zeffs(iz))
  enddo
  if(.false.) &
   call smu_ximu_CalcDAHChisqs(& 
    DAs, Hs, & ! Values of DA, H in six cosmologies
    omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
    chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
    weightedstds = .true., avg_counts = .false. &
    ) 
   do i = 1, N1
   do j = 1, N2
     if ( mubins(i) .eq. 36 .and. j.eq.1) then
       print *, '* nbin / mubin = ', mubins(i),j
       print *, '  chisqs (no cor) = ', real(chisqs_nosyscor(i,j,1:nz-1))
       print *, '  chisqs (corred) = ', real(chisqs_syscor(i,j,1:nz-1))
     endif
   enddo
   enddo
end program
