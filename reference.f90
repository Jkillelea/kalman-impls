program foo
    implicit none

    ! (cols, rows) format, different from C
    integer, parameter :: DIMS = 12
    integer, parameter :: JOBFILENO = 3
    integer, parameter :: DATAFILENO = 4
    character*256 :: msg, datafile
    integer :: idx, datalen, ierr
    real :: dt, tstamp, tstampprev, x(DIMS), F(DIMS, DIMS), H(DIMS, DIMS),     &
            P(DIMS, DIMS), Q(DIMS, DIMS), mu_e(DIMS), Sigma_e(DIMS, DIMS),     &
            mu_a(DIMS), Sigma_a(DIMS, DIMS), K(DIMS, DIMS), a(3), g(3),        &
            B(DIMS, DIMS), u(DIMS)

    B          = 0.0
    F          = 0.0
    H          = 0.0
    K          = 0.0
    P          = 0.0
    Q          = 0.0
    Sigma_a    = 0.0
    Sigma_e    = 0.0
    dt         = 0.0
    mu_a       = 0.0
    mu_e       = 0.0
    tstamp     = 0.0
    tstampprev = 0.0
    u          = 0.0
    x          = 0.0

    ! set identity
    do idx = 1,DIMS
        F(idx, idx) = 1.0
    enddo

    ! Measurement matrix is mostly identity
    do idx = 7,DIMS
        H(idx, idx) = 1.0
    enddo

    write (*, *) "F:"
    write (*, "(12(f6.2))") F
    write (*, *) "H:"
    write (*, "(12(f6.2))") H


    open(unit=JOBFILENO, file='jobctl.txt', form='formatted', status='old', iostat=ierr)
    if (ierr .ne. 0) then
        write (*, *) "failed to open input file ", 'jobctl.txt'
        stop
    endif

    read (JOBFILENO, '(a5)') msg
    read (JOBFILENO, '(i99)') datalen

    write (*, *) "msg: ", msg
    write (*, *) "value: ", datalen

    read (JOBFILENO, '(a5)') msg
    read (JOBFILENO, '(a99)') datafile

    write (*, *) "msg: ", msg
    write (*, *) "value: ", datafile

    open(unit=DATAFILENO, file=datafile, form='formatted', status='old', iostat=ierr)
    if (ierr .ne. 0) then
        write (*, *) "failed to open input file ", datafile
        stop
    endif

    do idx = 1,datalen+1
        ! data from accelerometer
        tstampprev = tstamp

        read (DATAFILENO, *) tstamp, mu_a(7:12)

        dt   = tstamp - tstampprev
        x    = matmul(F, x) + matmul(B, u)
        mu_e = matmul(H, x)
        P    = matmul(F, matmul(P, transpose(F))) + Q
        K    = matmul(P, transpose(H)) 

        write (*, *) 'dt:'
        write (*, "(f10.9)") dt
        write (*, *) 'x:'
        write (*, "(12f6.2)") x
        write (*, *) 'P:'
        write (*, "(12f6.2)") P
        write (*, *) 'mu_a:'
        write (*, "(12f6.2)") mu_a
        write (*, *) 'mu_e:'
        write (*, "(12f6.2)") mu_e

    enddo

    ! real    :: A(4, 3), B(3, 4)
    ! integer :: asize(2), nrows, ncols
    ! A = reshape((/2, 0, 0, 1,   &
    !               0, 0, 0, 0,   &
    !               0, 0, 0, 1/), &
    !           shape(A))

    ! B = reshape((/1, 0, 0,    &
    !               0, 1, 0,    &
    !               0, 0, 1,    &
    !               1, 0, 1 /), &
    !           shape(B))

    ! asize = shape(A)
    ! ncols = asize(1)
    ! nrows = asize(2)

    ! ! colum major
    ! write (*, *) "A (4 cols 3 rows)"
    ! write (*, "(4(f4.0))") A

    ! write (*, *) "B (3 cols 4 rows)"
    ! write (*, "(3(f4.0))") B

    ! write (*, *) "A*B = 4 rows, 4 cols"
    ! write (*, "(4(f4.0))") matmul(A, B)

    ! write (*, *) "B*A = 3 rows, 3 cols"
    ! write (*, "(3(f4.0))") matmul(B, A)

end program

