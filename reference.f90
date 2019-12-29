program foo
    implicit none

    ! (cols, rows) format, different from C
    integer, parameter :: DIMS = 12;
    real :: x(DIMS), F(DIMS, DIMS), H(DIMS, DIMS), P(DIMS, DIMS),       & 
            Q(DIMS, DIMS), mu_e(DIMS), Sigma_e(DIMS, DIMS), mu_a(DIMS), & 
            Sigma_a(DIMS, DIMS), K(DIMS, DIMS)

    write (*, *) "x"
    write (*, "(12(f4.0))") x
    write (*, *) "F"
    write (*, "(12(f4.0))") F

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

