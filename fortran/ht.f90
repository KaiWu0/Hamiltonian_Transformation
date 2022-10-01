! Hamiltonian Transformation for band structure calculation
!------------------------------------------------------------------------
PROGRAM ht
  USE kinds,              ONLY: dp
  USE constants,          ONLY: pi, eps4, eps6
  USE io_global,          ONLY: stdin, ionode, stdout
  USE mp_global,          ONLY: mp_startup
  USE io_files,           ONLY: prefix
  USE environment,        ONLY: environment_start, environment_end
  USE read_cards_module,  ONLY: read_cards
  USE noncollin_module,   ONLY: noncolin, npol
  USE ener,               ONLY: ef
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: outdir, eig_file
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: ios
  !
  COMPLEX(dp), ALLOCATABLE :: wfc(:, :, :)
  REAL(dp), ALLOCATABLE :: eigval(:,:)
  INTEGER :: nk123=0, nr123=0, time, num_threads=8, rank = -1
  ! 1:nbnd_discard bands and eigenvalues are discarded if window_min is set.
  INTEGER :: nbnd_discard
  ! nbnd_window=nbnd-nbnd_discard
  INTEGER :: nbnd_window, funtype=2
  REAL(dp) :: qr_eps=-1, window_min=0, power=3, gap=-1
  REAL(dp) :: search_start=-3, search_end=-20
  REAL(dp), PARAMETER :: ry2ev=13.6056980659
  LOGICAL :: delete_top_bands=.false., read_eig=.false.
  !
  NAMELIST /input_ht/ outdir, prefix, qr_eps, delete_top_bands,&
    num_threads, window_min, search_start, search_end,&
    power, gap, funtype, read_eig, eig_file, rank

  ! initialise environment
  CALL mp_startup()
  CALL environment_start('ht')

  ! read and check input variables
  IF (ionode) THEN
    !CALL input_from_file()
    !CALL get_environment_variable('ESPRESSO_TMPDIR', outdir)
    READ (stdin, input_ht, iostat=ios)
    IF (ios /= 0) CALL errore('ht', 'read input_ht error', abs(ios))
    IF (trim(outdir) == ' ') outdir = './'
    IF (qr_eps <= 0) qr_eps=4e-3
  END IF

  !! read previous calculation info
  CALL read_file()
  !! read nscf kpoint information, saved in input_parameters module
  CALL read_cards('PW',stdin)
  !
  window_min=window_min/ry2ev+ef
  search_start=search_start/ry2ev+ef
  search_end=search_end/ry2ev+ef

  ! main body of Hamiltonian Transformation
  CALL tick(time)
  CALL ht_init()
  print *, 'ht_init    : time = ', tock(time)
  FLUSH(stdout)
  CALL tick(time)
  CALL ht_main()
  print *, 'ht_main    : time = ', tock(time)

  !! clean up and exit
  CALL environment_end('ht')
  CALL stop_pp()
!------------------------------------------------------------------------
CONTAINS
  SUBROUTINE ht_init()
    !----------------------------------------------------------------------------
    !
    ! Prepare wave functions for Hamiltonian transformation.
    !
    USE wavefunctions,        ONLY: evc
    USE gvect,                ONLY: mill
    USE klist,                ONLY: ngk, igk_k, nkstot
    USE symm_base,            ONLY: nsym, sr
    USE io_files,             ONLY: iunwfc, nwordwfc
    USE control_flags,        ONLY: gamma_only, io_level
    USE buffers,              ONLY: get_buffer
    USE fft_scalar,           ONLY: cfft3d
    USE buffers,              ONLY: open_buffer, close_buffer
    USE start_k,              ONLY: nk1, nk2, nk3
    USE wvfct,                ONLY: et, nbnd, npwx
    !
    ! For debug
    USE klist,      ONLY: xk, nkstot
    USE symm_base,  ONLY: nsym, s, fft_fact, ft
    USE cell_base,  ONLY: at
    !
    IMPLICIT NONE
    !
    INTEGER :: ng(3), maxg(3), ming(3), ftau(3, nsym)
    INTEGER :: ik, ir, ibnd, ig, iq, npw, n1, n2, n3, ipol, jpol, isym
    INTEGER :: nblock=1, iblock
    LOGICAL :: exst
    COMPLEX(dp) :: d_spin(2,2,48)
    INTEGER, ALLOCATABLE :: nl(:,:), index_xk(:), index_sym(:)
    COMPLEX(dp), ALLOCATABLE :: phase(:)
    COMPLEX(dp), ALLOCATABLE :: temppsic(:), temppsic_nc(:,:), psic_nc(:,:)
    REAL(dp),ALLOCATABLE :: xk_cryst(:,:), xr_cryst(:,:)
    INTEGER, ALLOCATABLE :: rir(:, :)
    REAL(dp) :: block_value(2,nbnd), min_band, max_band
    !
    COMPLEX(dp), ALLOCATABLE :: tmp(:,:,:)
    !
    ! debug
    !
    !integer :: i
    !real(dp) :: gap=1,power=2,shift=0.2
    !real(dp),dimension(5,2) :: x,y,ex,dy
    !x = reshape((/ (dble(i)/5, i = -5, 4) /),shape(x))
    !gap=0.4
    !power=3
    !shift=0.6
    !y = fun(gap,power,shift,x)
    !dy = dfun(gap,power,shift,x)
    !y=newton_inv(gap,power,shift,y)
    !write(*,*) x
    !write(*,*) y
    !write(*,*) dy
    !
    ! Wavefunction can be saved to corser grid than in scf,
    ! symmetry has to be tackled.
    IF (gamma_only) THEN
      CALL errore('ht_init','Do not support gamma_only now.')
    END IF
    maxg = MAXVAL(mill(:, igk_k(1:ngk(1), 1)), 2)
    ming = MINVAL(mill(:, igk_k(1:ngk(1), 1)), 2)
    ng = maxg - ming + 3
#if defined(_OPENMP)
    CALL omp_set_num_threads(1)
#endif
    CALL real_space_fft_grid(ng, ftau)
    write(*,*) 'wavefunction grid:', ng
    nr123 = ng(1)*ng(2)*ng(3)
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nr123, npol), psic_nc(nr123, npol) )
    ELSE
       ALLOCATE( temppsic(nr123) )
    ENDIF
    !
    ALLOCATE(nl(npwx,nkstot))
    DO ik = 1, nkstot
      DO ig = 1, ngk(ik)
        n1 = mill(1, igk_k(ig,ik))
        n2 = mill(2, igk_k(ig,ik))
        n3 = mill(3, igk_k(ig,ik))
        IF (n1 < 0) n1 = n1 + ng(1)
        IF (n2 < 0) n2 = n2 + ng(2)
        IF (n3 < 0) n3 = n3 + ng(3)
        nl(ig, ik) = 1 + n1 + n2*ng(1) + n3*ng(1)*ng(2)
      END DO
    END DO
    !
    ! Analysis the band gap of system, the best window_min should in the gap.
    ! A block of bands is a set of entangled bands but isolated from other bands.
    ! block_value(1,i) is the min value of i_th block, block_value(2,i) is the max
    ! value of i_th block.
    DO ibnd=1,nbnd
      min_band=MINVAL(et(ibnd,:))
      max_band=MAXVAL(et(ibnd,:))
      IF (ibnd==1 .OR. min_band > block_value(2,nblock)+0.4/ry2ev) THEN
        IF (ibnd>1) nblock=nblock+1
        block_value(1,nblock)=min_band
        block_value(2,nblock)=max_band
      ELSE
        block_value(1,nblock)=MIN(block_value(1,nblock),min_band)
        block_value(2,nblock)=MAX(block_value(2,nblock),max_band)
      END IF
    END DO
    DO iblock=nblock-1,1,-1
      IF (block_value(2,iblock)<search_start) EXIT
    END DO
    IF (abs(window_min-ef)<1e-6 .AND. nbnd>50 .AND. iblock>0) THEN
      IF (block_value(1,iblock+1)<search_end) THEN
        ! Cannot find gap between search_start and search_end,
        window_min = (search_start+search_end)/2
      ELSE
        window_min = (block_value(2,iblock)+block_value(1,iblock+1))/2
      END IF
      write(*,*) window_min
      write(*,*) "Discard bands too far below fermi level."
    END IF
    IF (window_min>ef) write(*,*) "window_min larger than fermi level, do not&
      & discard bands."
    write(*,*) "Band blocks:"
    DO iblock=1,nblock
      write(*,'(I3, F10.2, F10.2)') iblock,(block_value(1,iblock)-ef)*ry2ev,&
      &(block_value(2,iblock)-ef)*ry2ev
    END DO
    print*,"==============================================================================="
    ! If system is too large, discard some lowest bands.
    IF (window_min<ef) THEN
      DO ibnd = 1,nbnd
        IF (MAXVAL(et(ibnd,:))>=window_min) EXIT
      END DO
      nbnd_discard = ibnd-1
      nbnd_window = nbnd - ibnd + 1
      IF (nbnd_window<10) CALL errore('ht_init','Too few bands left.&
        & Check window_min or set it to 0.',1)
      write(*,'("Discard ",I5," bands, ",I5, " bands left.")')&
        nbnd_discard, nbnd_window
    ELSE
      nbnd_discard = 0
      nbnd_window = nbnd
    END IF
    !
    !! Transform wavefunction to coarse real space grid.
    !! Similar to exxinit in exx.f90.
    !
    nk123 = nk1*nk2*nk3
    ALLOCATE(xk_cryst(3,nk123),index_xk(nk123),index_sym(nk123))
    ALLOCATE(wfc(nr123*npol, nbnd_window, nk123))
    CALL generate_k_grid(xk_cryst, index_xk, index_sym)
    ALLOCATE (rir(nr123, nsym), xr_cryst(3,nr123))
    CALL generate_r_grid(ng, rir, xr_cryst)
    !
    IF (noncolin) THEN
       DO isym = 1, nsym
          CALL find_u( sr(:,:,isym), d_spin(:,:,isym) )
       ENDDO
    ENDIF
    !
    CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
    ALLOCATE(phase(nr123),eigval(nbnd_window,nk123))
    ALLOCATE(tmp(nr123,nbnd_window,nk123))
    DO ik = 1, nkstot
      IF (nkstot > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      npw = ngk(ik)
      DO ibnd = 1, nbnd_window
        IF (noncolin) THEN
!$omp parallel do default(shared) private(ir) firstprivate(nr123)
          DO ir = 1, nr123
            temppsic_nc(ir, 1) = (0._DP, 0._DP)
            temppsic_nc(ir, 2) = (0._DP, 0._DP)
          END DO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik)
          DO ig = 1, npw
            temppsic_nc(nl(ig, ik), 1) = evc(ig, ibnd+nbnd_discard)
          END DO
          CALL cfft3d(temppsic_nc(:,1), ng(1), ng(2), ng(3), ng(1), ng(2), ng(3), 1, 1)
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,npwx)
          DO ig = 1, npw
            temppsic_nc(nl(ig, ik), 2) = evc(ig + npwx, ibnd+nbnd_discard)
          END DO
          CALL cfft3d(temppsic_nc(:,2), ng(1), ng(2), ng(3), ng(1), ng(2), ng(3), 1, 1)
        ELSE
!$omp parallel do default(shared) private(ir) firstprivate(nr123)
          DO ir = 1, nr123
            temppsic(ir) = (0._DP, 0._DP)
          END DO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik)
          DO ig = 1, npw
            temppsic(nl(ig, ik)) = evc(ig, ibnd+nbnd_discard)
          END DO
          CALL cfft3d(temppsic, ng(1), ng(2), ng(3), ng(1), ng(2), ng(3), 1, 1)
        END IF
        !
        DO iq = 1, nk123
          !
          IF (index_xk(iq) /= ik) CYCLE
          phase = EXP((0.0_dp,1.0_dp)*2*pi*(xk_cryst(1,iq)*xr_cryst(1,:)+&
            xk_cryst(2,iq)*xr_cryst(2,:)+xk_cryst(3,iq)*xr_cryst(3,:)))
          phase = phase/SQRT(dble(nr123))
          isym = ABS(index_sym(iq))
          eigval(ibnd,iq) = et(ibnd+nbnd_discard,ik)
          !
          IF (noncolin) THEN
!$omp parallel do default(shared) private(ir) firstprivate(npol,nr123)
            DO ir = 1, nr123
              !DIR$ UNROLL_AND_JAM (2)
              DO ipol = 1, npol
                psic_nc(ir, ipol) = (0._DP, 0._DP)
              END DO
            END DO
!$omp parallel do default(shared) private(ipol,jpol,ir) firstprivate(npol,isym,nr123) reduction(+:psic_nc)
            DO ir = 1, nr123
              !DIR$ UNROLL_AND_JAM (4)
              DO ipol = 1, npol
                DO jpol = 1, npol
                  psic_nc(ir, ipol) = psic_nc(ir, ipol) + CONJG(d_spin(jpol, ipol, isym))* &
                                      temppsic_nc(rir(ir, isym), jpol)
                END DO
              END DO
            END DO
            !
            IF (index_sym(iq) > 0) THEN
                ! sym. op. without time reversal: normal case
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,iq)
                DO ir = 1, nr123
                  wfc(ir, ibnd, iq) = psic_nc(ir, 1)
                  wfc(ir + nr123, ibnd, iq) = psic_nc(ir, 2)
                END DO
            ELSE
              ! sym. op. with time reversal: spin 1->2*, 2->-1*
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,iq)
                DO ir = 1, nr123
                  wfc(ir, ibnd, iq) = CONJG(psic_nc(ir, 2))
                  wfc(ir + nr123, ibnd, iq) = -CONJG(psic_nc(ir, 1))
                END DO
            END IF
!$omp parallel do default(shared) private(ir) firstprivate(nr123)
            DO ir = 1, nr123
              wfc(ir,ibnd, iq) = wfc(ir,ibnd, iq)*phase(ir)
              wfc(ir+nr123,ibnd, iq) = wfc(ir+nr123,ibnd, iq)*phase(ir)
            END DO
          ELSE ! noncolinear
!$omp parallel do default(shared) private(ir) firstprivate(isym,ibnd,iq)
            DO ir = 1, nr123
              IF (index_sym(iq) > 0) THEN
                wfc(ir, ibnd, iq) = temppsic(rir(ir, isym))
              ELSE
                wfc(ir, ibnd, iq) = CONJG(temppsic(rir(ir, isym)))
              END IF
            END DO
!$omp parallel do default(shared) private(ir) firstprivate(nr123)
            DO ir = 1, nr123
              wfc(ir,ibnd, iq) = wfc(ir,ibnd, iq)*phase(ir)
            END DO
            !
          END IF ! noncolinear
        END DO ! iq
      END DO ! ibnd
    END DO ! ik
    CALL close_buffer  ( iunwfc, 'KEEP' )
    !
    DEALLOCATE(xk_cryst,index_xk,index_sym,nl,rir,xr_cryst,phase)
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, psic_nc)
    ELSE
       DEALLOCATE( temppsic )
    ENDIF
  END SUBROUTINE ht_init
  !------------------------------------------------------------------------
  SUBROUTINE ht_main()
    USE start_k,          ONLY : nk1, nk2, nk3
    USE fft_scalar,       ONLY : cfft3d
    USE input_parameters, ONLY : xk, nkstot
    USE uspp,             ONLY : okvan
    !
    IMPLICIT NONE
    !
    REAL(DP), EXTERNAL :: DZNRM2
    INTEGER :: info, qrlwork, size(3), mindim, nbk, nkpath, tm
    INTEGER :: rk, i, j, k, ik, eiglwork, neig, lrwork, liwork, u, iq
    INTEGER :: Rpts(3,nk123), itop, ibtm
    INTEGER, ALLOCATABLE :: jpvt(:), isuppz(:), iwork(:)
    REAL(dp) :: shift, tmp_norm
    REAL(dp) :: trust_ev, xkpath(3,nkstot), epath(nbnd_window,nkstot)
    REAL(dp), ALLOCATABLE :: qrrwork(:), eigrwork(:), eigtmp(:)
    COMPLEX(dp), ALLOCATABLE :: tau(:), qrwork(:), R(:,:,:), tmp(:,:)
    COMPLEX(dp), ALLOCATABLE :: eigwork(:), eigvectors(:,:), Q(:,:)
    COMPLEX(dp), ALLOCATABLE, TARGET :: H(:,:,:), M(:,:,:)
    COMPLEX(dp), POINTER :: H1d(:), H2d(:,:), M2d(:,:)
    COMPLEX(dp) :: phase(nkstot,nk123)
    !
    real(dp),dimension(5,2) :: x,y,ex,dy
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_num_threads
#endif
    !
#if defined(_OPENMP)
    CALL omp_set_num_threads(num_threads)
#endif
    ! Rename xk and nkstot for nscf to avoid confusion.
    nkpath = nkstot
    xkpath = xk
    !
    IF (read_eig) THEN
      print *, 'read eigenvalues from file',trim(eig_file), ', unit is eV.'
      u = find_free_unit()
      OPEN(unit=u,file=trim(eig_file),status='old',action='read')
      DO j=1,nk123
        READ(u,*) (eigval(i,j),i=1,nbnd_window)
      END DO
      CLOSE(u)
      eigval=eigval/ry2ev
    END IF
    !
    shift = 1*MAXVAL(eigval(nbnd_window,:))+0*MINVAL(eigval(nbnd_window,:))
    IF (gap<0) gap=4*(MAXVAL(eigval(nbnd_window,:))-MINVAL(eigval(nbnd_window,:)))
    eigval = fun(gap,power,shift,eigval)
    !
    nbk = nbnd_window*nk123
    mindim = MIN(nr123*npol,nbk)
    size(1) = nr123*npol
    size(2) = nbk
    size(3) = 1
    wfc = reshape(wfc,size)
    !
    ! QR decomposition of wavefunction
    CALL tick(tm)
    ALLOCATE(jpvt(nbk))
    ALLOCATE(tau(mindim), qrwork(1))
    ALLOCATE(qrrwork(2*nbk))
    jpvt(:) = 0
    CALL ZGEQP3(nr123*npol,nbk,wfc,nr123*npol,jpvt,tau,qrwork,-1,qrrwork,info)
    IF(info/=0) CALL errore('ht_main','Error in computing the QR factorization',1)
    qrlwork = AINT(REAL(qrwork(1)))
    jpvt(:) = 0
    DEALLOCATE(qrwork)
    ALLOCATE(qrwork(qrlwork))
    CALL ZGEQP3(nr123*npol,nbk,wfc,nr123*npol,jpvt,tau,qrwork,qrlwork,qrrwork,info)
    IF(info/=0) CALL errore('ht_main','Error in computing the QR factorization',1)
    print *, 'QR         : time = ', tock(tm)
    FLUSH(stdout)
    !
    ! generate Q for debug
    !ALLOCATE(Q(size(1),size(2)))
    !Q=wfc(:,:,1)
    !CALL ZUNGQR(nr123*npol,mindim,mindim,Q,nr123*npol,tau,qrwork,qrlwork,info)
    !IF(info/=0) CALL errore('ht_main','Error in computing the QR factorization',1)
    !
    DO i = 2, mindim
      IF (ABS(wfc(i,i,1))<ABS(wfc(1,1,1))*qr_eps) EXIT
    END DO
    rk = i-1
    IF (rank>=0) THEN
      print *, "You have assigned rank to ", rank
      rk=rank
    ELSE
      print *, 'rank = ', rk
    END IF
    FLUSH(stdout)
    ALLOCATE(R(rk,nbk,1))
    IF (okvan) THEN
      print *,"-------------------------------------------------------------------------------"
      print *,"Warning: The pseudopotential is PAW or ultrasoft, so wavefunctions are"
      print *,"nonorthogonal. HT method can be generalized to nonorthogonal basis sets, but it"
      print *,"has not been implemented yet."
      print *,"Here we normalize wave functions and pretend they are orthogonal, this may"
      print *,"cause some small errors."
      print *,"-------------------------------------------------------------------------------"
    END IF
    DO i = 1, nbk
      IF (okvan) THEN
        tmp_norm = DZNRM2(MIN(i,rk),wfc(:,i,1),1)
      ELSE
        tmp_norm = 1
      END IF
      R(1:MIN(i,rk),jpvt(i),1) = wfc(1:MIN(i,rk),i,1)/tmp_norm
    END DO
    DEALLOCATE(jpvt,qrrwork,tau,qrwork,wfc)
    !
    size(1) = rk
    size(2) = nbnd_window
    size(3) = nk123
    R = reshape(R,size)
    ALLOCATE(H(nk123,rk,rk))
    ALLOCATE(tmp(rk,rk))
    H(:,:,:) = 0
    CALL tick(tm)
    DO i = 1, nk123
      DO j = 1, nbnd_window
        R(:,j,i) = R(:,j,i)*sqrt(-eigval(j,i))
      END DO
      tmp(:,:)=0
      CALL ZHERK('U','N',rk,nbnd_window,-1.0_dp,R(:,:,i),rk,0,tmp,rk)
      H(i,:,:)=tmp
    ! H is Hermitian upper triangular matrix, recover down triangular elements.
      DO j = 2, rk
        H(i,j,1:j-1) = CONJG(H(i,1:j-1,j))
      END DO
    END DO
    print *, 'build H    : time = ', tock(tm)
    FLUSH(stdout)
    DEALLOCATE(R)
    !
    ! Fourier interpolation. FFT
    CALL tick(tm)
    H1d(1:nk123*rk*rk) => H
    CALL cfft3d(H1d, nk1, nk2, nk3, nk1, nk2, nk3, rk*rk, -1)
    H2d(1:nk123,1:rk*rk) => H
    print *, 'First FFT  : time = ', tock(tm)
    FLUSH(stdout)
    ! Interpolation back.
    ALLOCATE(M(nkpath,rk,rk))
    DO ik = 0, nk123-1
      k = ik/(nk1*nk2)
      j = (ik-k*nk1*nk2)/nk1
      i = ik-k*nk1*nk2-j*nk1
      IF (i>=(nk1+1)/2) i=i-nk1
      IF (j>=(nk2+1)/2) j=j-nk2
      IF (k>=(nk3+1)/2) k=k-nk3
      Rpts(1,ik+1)=i
      Rpts(2,ik+1)=j
      Rpts(3,ik+1)=k
    END DO
    DO iq = 1, nkpath
      phase(iq,:) = EXP((0.0_dp,2.0_dp)*pi*&
        (xkpath(1,iq)*Rpts(1,:) + xkpath(2,iq)*Rpts(2,:) + xkpath(3,iq)*Rpts(3,:)))
    END DO
    M2d(1:nkpath,1:rk*rk) => M
    CALL tick(tm)
    CALL ZGEMM('N','N',nkpath,rk*rk,nk123,(0.5_dp,0.0_dp),phase,&
      nkpath,H2d,nk123,(0.0_dp,0.0_dp),M2d,nkpath)
    DO iq = 1, nkpath
      M(iq,:,:) = M(iq,:,:) + TRANSPOSE(CONJG(M(iq,:,:)))
    END DO
    DEALLOCATE(H)
    print *, 'Second FT  : time = ', tock(tm)
    FLUSH(stdout)
    !
    ! Calculate eigenvalues
    CALL tick(tm)
    lrwork = rk*24
    liwork = rk*10
    ALLOCATE(eigvectors(rk,rk),isuppz(2*rk),eigwork(1))
    ALLOCATE(eigrwork(lrwork),iwork(liwork),eigtmp(rk))
    CALL ZHEEVR('N','I','U',rk,tmp,rk,0.0_dp,0.0_dp,1,nbnd_window,1e-16,neig,eigtmp, &
      eigvectors,rk,isuppz,eigwork,-1,eigrwork,lrwork,iwork,liwork,info)
    IF(info/=0) CALL errore('ht_main','Error in eigvenvalue decomposition',1)
    eiglwork = AINT(REAL(eigwork(1)))
    DEALLOCATE(eigwork)
    ALLOCATE(eigwork(eiglwork))
    DO ik = 1, nkpath
      tmp=M(ik,:,:)
      CALL ZHEEVR('N','I','U',rk,tmp,rk,0.0_dp,0.0_dp,1,nbnd_window,1e-16,neig,eigtmp, &
        eigvectors,rk,isuppz,eigwork,eiglwork,eigrwork,lrwork,iwork,liwork,info)
      IF(info/=0) CALL errore('ht_main','Error in eigvenvalue decomposition',1)
      epath(:,ik) = eigtmp(1:nbnd_window)
    END DO
    DEALLOCATE(eigwork,eigrwork,M,eigvectors,iwork,isuppz,eigtmp,tmp)
    print *, 'Diagonalize: time = ', tock(tm)
    FLUSH(stdout)
    !
    epath = newton_inv(gap,power,shift,epath)
    !
    !delete top bands
    itop = 0
    IF (delete_top_bands) THEN
      trust_ev = MINVAL(eigval(nbnd_window,:))+1e-4+shift;
      DO itop = MIN(nbnd_window-1,10),0,-1
        IF (MAXVAL(epath(nbnd_window-itop,:))>trust_ev) EXIT
      END DO
      itop=itop+1
      print *, itop, 'top bands are deleted.'
      IF (nbnd_window<10) THEN
        print *, 'If top bands are distorted, increase nbnd in scf calculation.'
      END IF
    END IF
    IF (window_min<ef) THEN
      DO ibtm = nbnd_window-itop,1,-1
        IF (MINVAL(epath(ibtm,:))<window_min) EXIT
      END DO
      print *, ibtm+nbnd_discard, 'bottom bands are deleted.'
    END IF
    DEALLOCATE(eigval)
    ! Output eigenvalues.
    u = find_free_unit()
    OPEN(newunit=u,file='band.txt')
    DO i=1,nbnd_window-itop
      WRITE(u,"(*(F16.6))") epath(i,:)*ry2ev
    END DO
    CLOSE(u)
    !
    !----------------------------------------------------------------------------
  END SUBROUTINE ht_main
  !----------------------------------------------------------------------------
  SUBROUTINE real_space_fft_grid(ng, ftau)
    !
    !! Generate smallest wavefunction FFT grid compatible with
    !! rotation matrices. Similar to scale_sym_ops in ruotaijk.f90.
    !
    USE fft_param, ONLY: nfftx
    USE symm_base, ONLY: nsym, s, fft_fact, ft
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT)  :: ng(3)
    INTEGER, INTENT(OUT)  :: ftau(3, nsym)
    INTEGER :: i1, i2, i3, isym, min_n123, n1, n2, n3
    REAL(dp) :: ftn(3, nsym)
    !
    min_n123 = nfftx*nfftx*nfftx
    loop_1: DO i1 = ng(1), MIN(ng(1)*2 + 1, nfftx)
      ftn(1, :) = ft(1, 1:nsym)*i1
      IF (MAXVAL(ABS(ftn(1, :) - NINT(ftn(1, :)))) > eps4) CYCLE
      loop_2: DO i2 = ng(2), MIN(ng(2)*2 + 1, nfftx)
        ftn(2, :) = ft(2, 1:nsym)*i2
        IF (MAXVAL(ABS(ftn(2, :) - NINT(ftn(2, :)))) > eps4) CYCLE
        DO isym = 1, nsym
          IF (MOD(s(2, 1, isym)*i1, i2) /= 0 .OR. MOD(s(1, 2, isym)*i2, i1) /= 0) CYCLE loop_2
        END DO
        loop_3: DO i3 = ng(3), MIN(ng(3)*2 + 1, nfftx)
          ftn(3, :) = ft(3, 1:nsym)*i3
          IF (MAXVAL(ABS(ftn(3, :) - NINT(ftn(3, :)))) > eps4) CYCLE
          DO isym = 1, nsym
            IF (MOD(s(3, 1, isym)*i1, i3) /= 0 .OR. MOD(s(1, 3, isym)*i3, i1) /= 0) CYCLE loop_3
            IF (MOD(s(3, 2, isym)*i2, i3) /= 0 .OR. MOD(s(2, 3, isym)*i3, i2) /= 0) CYCLE loop_3
          END DO
          IF (i1*i2*i3 >= min_n123) CYCLE loop_2
          IF (MOD(i1, fft_fact(1)) == 0 .AND. MOD(i2, fft_fact(2)) == 0 &
              .AND. MOD(i3, fft_fact(3)) == 0) THEN
            n1 = i1
            n2 = i2
            n3 = i3
            ftau = NINT(ftn)
            min_n123 = i1*i2*i3
          END IF
        END DO loop_3
      END DO loop_2
    END DO loop_1
    IF (min_n123 == nfftx*nfftx*nfftx) THEN
      CALL errore('ht', 'Cannot find corse real space fft grid.', 1)
    END IF
    ng(1) = n1
    ng(2) = n2
    ng(3) = n3
  END SUBROUTINE real_space_fft_grid
  !----------------------------------------------------------------------------
  SUBROUTINE generate_k_grid(xk_cryst, index_xk, index_sym)
    !
    !! Find all k-points equivalent by symmetry to the points in the k-list.
    !! Similar to exx_grid_init in exx_base.f90, but here k-points must be
    !! exactly same before and after symmetry operations.
    !
    USE symm_base,        ONLY: nsym, s, no_t_rev
    USE cell_base,        ONLY: at
    USE klist,            ONLY: xk, nkstot
    USE start_k,          ONLY: nk1, nk2, nk3
    USE noncollin_module, ONLY: domag
    USE control_flags,    ONLY: noinv
    !
    IMPLICIT NONE
    !
    INTEGER :: isym, ik, iq, i, j, k
    INTEGER, INTENT(OUT), DIMENSION(nk123) :: index_xk, index_sym
    REAL(DP), INTENT(OUT) :: xk_cryst(3, nk123)
    LOGICAL :: k_found
    REAL(DP) :: sxk(3), xk_cryst_sym(3)
    !
    DO ik = 0, nk123-1
      k = ik/(nk1*nk2)
      j = (ik-k*nk1*nk2)/nk1
      i = ik-k*nk1*nk2-j*nk1
      IF (i>=(nk1+1)/2) i=i-nk1
      IF (j>=(nk2+1)/2) j=j-nk2
      IF (k>=(nk3+1)/2) k=k-nk3
      xk_cryst(1,ik+1)=dble(i)/nk1
      xk_cryst(2,ik+1)=dble(j)/nk2
      xk_cryst(3,ik+1)=dble(k)/nk3
    END DO
    qloop: DO iq = 1, nk123
      k_found = .FALSE.
      DO ik = 1, nkstot
        IF (k_found) EXIT
        xk_cryst_sym = MATMUL(transpose(at), xk(:, ik))
        DO isym = 1, nsym
          IF (k_found) EXIT
          sxk = MATMUL(s(:,:,isym),xk_cryst_sym)
          IF (SUM(ABS(sxk - xk_cryst(:,iq))) <= eps6) THEN
            k_found = .TRUE.
            index_xk(iq) = ik
            index_sym(iq) = isym
          ELSE IF (.NOT. noinv .AND. SUM(ABS(-sxk - xk_cryst(:,iq))) <= eps6) THEN
            k_found = .TRUE.
            index_xk(iq) = ik
            index_sym(iq) = -isym
          END IF
        END DO
      END DO
      IF (.NOT. k_found) THEN
        IF (noncolin .AND. domag .AND. .NOT. no_t_rev) THEN
          CALL errore('ht - generate_k_grid', "Cannot find required k-points, &
          &set no_t_rev = .true. in the scf calculation.", 1)
        ELSE
          CALL errore('ht - generate_k_grid', "Cannot find required k-points, &
          &reduce or close symmetry in the scf calculation.", 1)
        END IF
      END IF
    END DO qloop
  END SUBROUTINE generate_k_grid
  !----------------------------------------------------------------------------
  SUBROUTINE generate_r_grid(nr, rir, xr_cryst)
    !-----------------------------------------------------------------------
    !! Generate rotated r grid, refer to exx_set_symm in exx_base.f90.
    !
    USE symm_base, ONLY: nsym, s, ft
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nr(3)
    INTEGER, INTENT(out) :: rir(nr123,nsym)
    REAL(dp), INTENT(out) :: xr_cryst(3,nr123)
    !
    ! ... local variables
    !
    INTEGER :: isym, i, j, k, ri, rj, rk, ir, nr1, nr2, nr3
    INTEGER, allocatable :: ftau(:, :), s_scaled(:, :, :)
    !
    rir = 0
    nr1 = nr(1)
    nr2 = nr(2)
    nr3 = nr(3)
    ALLOCATE (ftau(3, nsym), s_scaled(3, 3, nsym))
    CALL scale_sym_ops(nsym, s, ft, nr1, nr2, nr3, s_scaled, ftau)
    DO isym = 1, nsym
      DO k = 1, nr3
        DO j = 1, nr2
          DO i = 1, nr1
            CALL rotate_grid_point(s_scaled(1, 1, isym), ftau(1, isym), &
                                   i, j, k, nr1, nr2, nr3, ri, rj, rk)
            ir = i + (j - 1)*nr1 + (k - 1)*nr1*nr2
            rir(ir, isym) = ri + (rj - 1)*nr1 + (rk - 1)*nr1*nr2
            IF (isym == 1) THEN
              xr_cryst(1,ir) = dble(i-1)/nr1
              xr_cryst(2,ir) = dble(j-1)/nr2
              xr_cryst(3,ir) = dble(k-1)/nr3
            END IF
          END DO
        END DO
      END DO
    END DO
    !
    DEALLOCATE (s_scaled, ftau)
    !
  END SUBROUTINE generate_r_grid
  !----------------------------------------------------------------------------
  FUNCTION fun(a,n,s,x) RESULT(y)
    !-----------------------------------------------------------------------
    !! Transform function which makes Hamiltonian decays faster in real space.
    !
    IMPLICIT NONE
    !
    REAL(dp), INTENT(IN) :: a, n, s, x(:,:)
    REAL(dp), DIMENSION(SIZE(x,1),SIZE(x,2)) :: y
    !
    y=x-s
    SELECT CASE (funtype)
    CASE (0)
      WHERE (y>=0)
        y = 0
      ELSEWHERE (y<=-a)
        y = y + a*(1.-1./n)
      ELSEWHERE
        y = -(-y)**n/n/(a**(n-1))
      END WHERE
    CASE (1)
      WHERE (y>=0)
        y = 0
      ELSEWHERE (y<=-a)
        y = y + a/2
      ELSEWHERE
        y = a*(exp(-n*n/4)-exp(-n*n*(a+2*y)*(a+2*y)/4/a/a))/2/n/sqrt(pi)&
          -a*erfc(n/2)/4+(a/4+y/2)*erfc(n*(a/2+y)/a)
      END WHERE
    CASE (2)
      WHERE (y>=0)
        y = 0
      ELSEWHERE (y<=-a)
        y = y + a/2
      ELSEWHERE
        y = a*(exp(-n*n/4)-exp(-n*n*(a+2*y)*(a+2*y)/4/a/a))/2/n/sqrt(pi)/erf(n/2)&
          +(a+2*y)*(erf(n/2)-erf(n*(0.5+y/a)))/4/erf(n/2)
      END WHERE
    CASE DEFAULT
      CALL errore('ht - fun',"Error: wrong funtype!",1)
    END SELECT
  END FUNCTION fun
  !----------------------------------------------------------------------------
  FUNCTION dfun(a,n,s,x) RESULT(y)
    !-----------------------------------------------------------------------
    !! Derivative of transform function.
    !
    IMPLICIT NONE
    !
    REAL(dp), INTENT(IN) :: a, n, s, x(:,:)
    REAL(dp), DIMENSION(SIZE(x,1),SIZE(x,2)) :: y
    !
    y=x-s
    SELECT CASE (funtype)
    CASE (0)
      WHERE (y>=0)
        y = 0
      ELSEWHERE (y<=-a)
        y = 1
      ELSEWHERE
        y = (-y/a)**(n-1)
      END WHERE
    CASE (1)
      WHERE (y>=0)
        y = 0
      ELSEWHERE (y<=-a)
        y = 1
      ELSEWHERE
        y = erfc(n*(a/2+y)/a)/2
      END WHERE
    CASE (2)
      WHERE (y>=0)
        y = 0
      ELSEWHERE (y<=-a)
        y = 1
      ELSEWHERE
        y = 0.5-erf(n*(0.5+y/a))/2/erf(n/2)
      END WHERE
    CASE DEFAULT
      CALL errore('ht - dfun',"Error: wrong funtype!",1)
    END SELECT
  END FUNCTION dfun
  !----------------------------------------------------------------------------
  FUNCTION newton_inv(a,n,s,y) RESULT(x)
    !-----------------------------------------------------------------------
    !! Newton's method to compute the inverse of transform function.
    !
    IMPLICIT NONE
    !
    REAL(dp), INTENT(IN) :: a, n, s, y(:,:)
    REAL(dp), DIMENSION(SIZE(y,1),SIZE(y,2)) :: x, dx, res, df
    INTEGER :: nitermax=50, niter
    REAL(dp) :: dxmax, resmax=1e-12
    !
    ! Here all y<=0. Otherwise we must set y>0 elements to 0.
    dxmax=a/2
    x=y+s
    res=1.0
    niter=0
    DO WHILE (MAXVAL(ABS(res))>resmax .AND. niter<nitermax)
      res=fun(a,n,s,x)-y
      df=dfun(a,n,s,x)
      dx=0
      WHERE(df/=0) dx=-res/df
      WHERE(dx> dxmax) dx= dxmax
      WHERE(dx<-dxmax) dx=-dxmax
      x = x + dx
      niter = niter + 1
    END DO
    IF (MAXVAL(ABS(res))>resmax) CALL errore('ht - newton_inv',&
          "Newton's method does not converge.",1)
  END FUNCTION newton_inv
  !----------------------------------------------------------------------------
  SUBROUTINE tick(t)
    integer, intent(OUT) :: t
    call system_clock(t)
  END SUBROUTINE tick
  !----------------------------------------------------------------------------
  REAL FUNCTION tock(t)
    ! returns time in seconds from now to time described by t 
    integer, intent(in) :: t
    integer :: now, clock_rate
    call system_clock(now,clock_rate)
    tock = real(now - t)/real(clock_rate)
  END FUNCTION tock
  !----------------------------------------------------------------------------
END PROGRAM ht
