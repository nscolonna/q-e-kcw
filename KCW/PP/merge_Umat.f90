!#define DEBUG

PROGRAM merge
 ! Utility program to merge different U matrices in a single U
 ! Usage: merge_Umat.x <what> <path/to/what1> <path/to/what2> ... <path/to/whatN>
 ! <what>:  U matrices or centres (seedname_u.mat or seedname_centres.xyz
 !
 IMPLICIT NONE
 integer :: nk1, nk2, nb1, nb2
 integer :: nb
 complex*16, allocatable :: U1(:,:,:)
 complex*16, allocatable :: U2(:,:,:)
 complex*16, allocatable :: U(:,:,:)
 real*16, allocatable :: k(:,:)
 real*16 :: re, im
 integer :: ik, ib, jb, ii, iib, jjb
 character(256) :: dum,arg, what
 character(256), allocatable :: filename(:)
 integer :: nargs, iarg, ni, nf, nb1_
 !
 nargs = command_argument_count()
 !write(*,*) nargs 
 IF (nargs == 0 ) STOP
 allocate (filename(nargs) )
 WRITE(*,*) "Merging ",nargs-1, " blocks" !, what
 !
 nb=0
 iarg = 1
 CALL get_command_argument( iarg, arg )
 what=trim(arg)
 WRITE(*,*) "what = ",  trim(seedname)
 IF ( trim(what) == "centres" ) THEN 
    WRITE(*,*) "Merging centres not implemented yet. GOING to STOP"
    STOP
 ENDIF
 !
 DO WHILE ( iarg < nargs )
   iarg=iarg +1
   CALL get_command_argument( iarg, arg )
   !write (*,*) trim(arg), seedname, 'u_mat'
   filename(iarg)=trim(arg) !//'/'//trim(seedname)//trim('_u.mat')
   open (27, file=filename(iarg))
   read(27,*) dum
   read(27,*) nk1, nb1, nb1
   close(27)
   WRITE(*,*) "block # ", iarg-1, "file = ", trim(filename(iarg)), " nbnd = ", nb1
   nb = nb + nb1
 ENDDO
 !
 write(*,*) "Final # of bands   = ", nb
 write(*,*) "Final # of kpoints = ", nk1
 ALLOCATE( U(nb,nb,nk1))
 ALLOCATE (k(3,nk1))
 U=cmplx(0.D0,0.D0)
 !
 nb1_=0
 do iarg = 2, nargs
   open (27, file=filename(iarg))
   WRITE(*,*) "Reading ", TRIM(filename(iarg))
   read(27,*) dum
   read(27,*) nk1, nb1, nb1
   !write(*, '(3I12)') nk1, nb1, nb1
   ni=1+nb1_
   nf=nb1+nb1_
   write(*,*) "istart = ", ni, " iend = ", nf !, " increment =", nb1_
   do ik = 1, nk1
      read (27, *)
      !write(*,*) 
      read  (27, '(f15.10,sp,f15.10,sp,f15.10)') k(:, ik)
      !write (* , '(f15.10,sp,f15.10,sp,f15.10)') k(:, ik)
      !write(*,'(f15.10,sp,f15.10,sp,f15.10)')  k(:,ik)
      read  (27, '(f15.10,sp,f15.10)') ((U(ib, jb, ik), ib=ni, nf), jb=ni, nf)
      !write (* , '(f15.10,sp,f15.10)') ((U(ib, jb, ik), ib=ni, nf), jb=ni, nf)
      !read (27, *) ((U(ib, jb, ik), ib=ni, nf), jb=ni, nf)
   enddo
   nb1_=nb1_+nb1
   close(27)
   !
 end do
 !
 write(*,*) "Final # of bands   = ", nb
 write(*,*) "Final # of kpoints = ", nk1
 !WRITE(*,'(2(f15.10,sp,f15.10))') U(1,1,1), U(nb,nb,nk1)
 WRITE(*,*) "Writing ", TRIM('wann_u.mat')
 open (40, file='wann_u.mat', form='formatted')
 write (40, *) "Written by merge.x" 
 write (40, '(3I12)') nk1, nb, nb
 do ik = 1, nk1
   write(40, *)
   write(40, '(f15.10,sp,f15.10,sp,f15.10)') k(:, ik)
   write(40, '(f15.10,sp,f15.10)') ((U(ib, jb, ik), ib=1, nb), jb=1, nb)
   !write(40, *) ((U(ib, jb, ik), ib=1, nb), jb=1, nb)
   !write(*, *) ((U(ib, jb, ik), ib=1, nb), jb=1, nb)
 enddo
 WRITE (*,'(/, " GAME OVER", /)') 
 !
 !
END program
