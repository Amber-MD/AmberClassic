c----------------------------------------------------------------------
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
c
c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     March 31, 1993
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  DGETRF computes an LU factorization of a general M-by-N matrix A
c  using partial pivoting with row interchanges.
c
c  The factorization has the form
c     A = P * L * U
c  where P is a permutation matrix, L is lower triangular with unit
c  diagonal elements (lower trapezoidal if m > n), and U is upper
c  triangular (upper trapezoidal if m < n).
c
c  This is the right-looking Level 3 BLAS version of the algorithm.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          On entry, the M-by-N matrix to be factored.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
c                has been completed, but the factor U is exactly
c                singular, and division by zero will occur if it is used
c                to solve a system of equations.
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
c     ..
c     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
c
c     Determine the block size for this environment.
c
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
c
c        Use unblocked code.
c
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
c
c        Use blocked code.
c
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
c
c           Factor diagonal and subdiagonal blocks and test for exact
c           singularity.
c
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
c
c           Adjust INFO and the pivot indices.
c
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
c
c           Apply interchanges to columns 1:J-1.
c
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
c
            IF( J+JB.LE.N ) THEN
c
c              Apply interchanges to columns J+JB:N.
c
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
c
c              Compute block row of U.
c
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
c
c                 Update trailing submatrix.
c
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
c
c     End of DGETRF
c
      END
c***********************************************************************
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c
c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     March 31, 1993
c
c     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  DGETRS solves a system of linear equations
c     A * X = B  or  A_prime * X = B
c  with a general N-by-N matrix A using the LU factorization computed
c  by DGETRF.
c
c  Arguments
c  =========
c
c  TRANS   (input) CHARACTER*1
c          Specifies the form of the system of equations:
c          = "N":        A * X = B  (No transpose)
c          = "T":  A_prime * X = B  (Transpose)
c          = "C":  A_prime * X = B  (Conjugate transpose = Transpose)
c
c  N       (input) INTEGER
c          The order of the matrix A.  N >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
c          The factors L and U from the factorization A = P*L*U
c          as computed by DGETRF.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,N).
c
c  IPIV    (input) INTEGER array, dimension (N)
c          The pivot indices from DGETRF; for 1<=i<=N, row i of the
c          matrix was interchanged with row IPIV(i).
c
c  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
c          On entry, the right hand side matrix B.
c          On exit, the solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            NOTRAN
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
c
      IF( NOTRAN ) THEN
c
c        Solve A * X = B.
c
c        Apply row interchanges to the right hand sides.
c
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
c
c        Solve L*X = B, overwriting B with X.
c
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
c
c        Solve U*X = B, overwriting B with X.
c
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
c
c        Solve A_prime * X = B.
c
c        Solve U_prime * X = B, overwriting B with X.
c
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
c
c        Solve L_prime * X = B, overwriting B with X.
c
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
c
c        Apply row interchanges to the solution vectors.
c
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
c
      RETURN
c
c     End of DGETRS
c
      END
c***********************************************************************
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
c
c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     June 30, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  DGETF2 computes an LU factorization of a general m-by-n matrix A
c  using partial pivoting with row interchanges.
c
c  The factorization has the form
c     A = P * L * U
c  where P is a permutation matrix, L is lower triangular with unit
c  diagonal elements (lower trapezoidal if m > n), and U is upper
c  triangular (upper trapezoidal if m < n).
c
c  This is the right-looking Level 2 BLAS version of the algorithm.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          On entry, the m by n matrix to be factored.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -k, the k-th argument had an illegal value
c          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            J, JP
c     ..
c     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
c
c     Quick return if possible
c
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
c
      DO 10 J = 1, MIN( M, N )
c
c        Find pivot and test for singularity.
c
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
c
c           Apply the interchange to columns 1:N.
c
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
c
c           Compute elements J+1:M of J-th column.
c
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
c
         ELSE IF( INFO.EQ.0 ) THEN
c
            INFO = J
         END IF
c
         IF( J.LT.MIN( M, N ) ) THEN
c
c           Update trailing submatrix.
c
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
c
c     End of DGETF2
c
      END
c***********************************************************************
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
c
c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  DLASWP performs a series of row interchanges on the matrix A.
c  One row interchange is initiated for each of rows K1 through K2 of A.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.
c
c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          On entry, the matrix of column dimension N to which the row
c          interchanges will be applied.
c          On exit, the permuted matrix.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.
c
c  K1      (input) INTEGER
c          The first element of IPIV for which a row interchange will
c          be done.
c
c  K2      (input) INTEGER
c          The last element of IPIV for which a row interchange will
c          be done.
c
c  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
c          The vector of pivot indices.  Only the elements in positions
c          K1 through K2 of IPIV are accessed.
c          IPIV(K) = L implies rows K and L are to be interchanged.
c
c  INCX    (input) INTEGER
c          The increment between successive values of IPIV.  If IPIV
c          is negative, the pivots are applied in reverse order.
c
c =====================================================================
c
c     .. Local Scalars ..
      INTEGER            I, IP, IX
c     ..
c     .. External Subroutines ..
      EXTERNAL           DSWAP
c     ..
c     .. Executable Statements ..
c
c     Interchange row I with row IPIV(I) for each of rows K1 through K2.
c
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
c
      RETURN
c
c     End of DLASWP
c
      END
c***********************************************************************
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
c
c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
c
c     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
c     ..
c
c  Purpose
c  =======
c
c  ILAENV is called from the LAPACK routines to choose problem-dependent
c  parameters for the local environment.  See ISPEC for a description of
c  the parameters.
c
c  This version provides a set of parameters which should give good,
c  but not optimal, performance on many of the currently available
c  computers.  Users are encouraged to modify this subroutine to set
c  the tuning parameters for their particular machine using the option
c  and problem size information in the arguments.
c
c  This routine will not function correctly if it is converted to all
c  lower case.  Converting it to all upper case is allowed.
c
c  Arguments
c  =========
c
c  ISPEC   (input) INTEGER
c          Specifies the parameter to be returned as the value of
c          ILAENV.
c          = 1: the optimal blocksize; if this value is 1, an unblocked
c               algorithm will give the best performance.
c          = 2: the minimum block size for which the block routine
c               should be used; if the usable block size is less than
c               this value, an unblocked routine should be used.
c          = 3: the crossover point (in a block routine, for N less
c               than this value, an unblocked routine should be used)
c          = 4: the number of shifts, used in the nonsymmetric
c               eigenvalue routines
c          = 5: the minimum column dimension for blocking to be used;
c               rectangular blocks must have dimension at least k by m,
c               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
c          = 6: the crossover point for the SVD (when reducing an m by n
c               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
c               this value, a QR factorization is used first to reduce
c               the matrix to a triangular form.)
c          = 7: the number of processors
c          = 8: the crossover point for the multishift QR and QZ methods
c               for nonsymmetric eigenvalue problems.
c
c  NAME    (input) CHARACTER*(*)
c          The name of the calling subroutine, in either upper case or
c          lower case.
c
c  OPTS    (input) CHARACTER*(*)
c          The character options to the subroutine NAME, concatenated
c          into a single character string.  For example, UPLO = "U",
c          TRANS = "T", and DIAG = "N" for a triangular routine would
c          be specified as OPTS = "UTN".
c
c  N1      (input) INTEGER
c  N2      (input) INTEGER
c  N3      (input) INTEGER
c  N4      (input) INTEGER
c          Problem dimensions for the subroutine NAME; these may not all
c          be required.
c
c (ILAENV) (output) INTEGER
c          >= 0: the value of the parameter specified by ISPEC
c          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
c
c  Further Details
c  ===============
c
c  The following conventions have been used when calling ILAENV from the
c  LAPACK routines:
c  1)  OPTS is a concatenation of all of the character options to
c      subroutine NAME, in the same order that they appear in the
c      argument list for NAME, even if they are not used in determining
c      the value of the parameter specified by ISPEC.
c  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
c      that they appear in the argument list for NAME.  N1 is used
c      first, N2 second, and so on, and unused problem dimensions are
c      passed a value of -1.
c  3)  The parameter value returned by ILAENV is checked for validity in
c      the calling subroutine.  For example, ILAENV is used to retrieve
c      the optimal blocksize for STRTRI as follows:
c
c      NB = ILAENV( 1, "STRTRI", UPLO // DIAG, N, -1, -1, -1 )
c      IF( NB.LE.1 ) NB = MAX( 1, N )
c
c  =====================================================================
c
c     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
c     ..
c     .. Executable Statements ..
c
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
c
c     Invalid value for ISPEC
c
      ILAENV = -1
      RETURN
c
  100 CONTINUE
c
c     Convert NAME to upper case if the first character is lower case.
c
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
c
c        ASCII character set
c
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
c
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
c
c        EBCDIC character set
c
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
c
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
c
c        Prime machines:  ASCII+128
c
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
c
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
c
      GO TO ( 110, 200, 300 ) ISPEC
c
  110 CONTINUE
c
c     ISPEC = 1:  block size
c
c     In these examples, separate code is provided for setting NB for
c     real and complex.  We assume that NB will take the same value in
c     single or double precision.
c
      NB = 1
c
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
c
  200 CONTINUE
c
c     ISPEC = 2:  minimum block size
c
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
c
  300 CONTINUE
c
c     ISPEC = 3:  crossover point
c
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
c
  400 CONTINUE
c
c     ISPEC = 4:  number of shifts (used by xHSEQR)
c
      ILAENV = 6
      RETURN
c
  500 CONTINUE
c
c     ISPEC = 5:  minimum column dimension (not used)
c
      ILAENV = 2
      RETURN
c
  600 CONTINUE 
c
c     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
c
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
c
  700 CONTINUE
c
c     ISPEC = 7:  number of processors (not used)
c
      ILAENV = 1
      RETURN
c
  800 CONTINUE
c
c     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
c
      ILAENV = 50
      RETURN
c
c     End of ILAENV
c
      END
c***********************************************************************
      LOGICAL          FUNCTION LSAME( CA, CB )
c
c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
c
c     .. Scalar Arguments ..
      CHARACTER          CA, CB
c     ..
c
c  Purpose
c  =======
c
c  LSAME returns .TRUE. if CA is the same letter as CB regardless of
c  case.
c
c  Arguments
c  =========
c
c  CA      (input) CHARACTER*1
c  CB      (input) CHARACTER*1
c          CA and CB specify the single characters to be compared.
c
c =====================================================================
c
c     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
c     ..
c     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
c     ..
c     .. Executable Statements ..
c
c     Test if the characters are equal
c
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
c
c     Now test for equivalence if both characters are alphabetic.
c
      ZCODE = ICHAR( 'Z' )
c
c     Use "Z" rather than "A" so that ASCII can be detected on Prime
c     machines, on which ICHAR returns a value with bit 8 set.
c
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
c
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
c
c        ASCII is assumed - ZCODE is the ASCII code of either lower or
c        upper case "Z".
c
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
c
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
c
c        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
c        upper case "Z".
c
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
c
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
c
c        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
c        plus 128 of either lower or upper case "Z".
c
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
c
c     RETURN
c
c     End of LSAME
c
      END
c***********************************************************************
      SUBROUTINE XERBLA( SRNAME, INFO )
c
c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
c
c     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
c     ..
c
c  Purpose
c  =======
c
c  XERBLA  is an error handler for the LAPACK routines.
c  It is called by an LAPACK routine if an input parameter has an
c  invalid value.  A message is printed and execution stops.
c
c  Installers may consider modifying the STOP statement in order to
c  call system-specific exception-handling facilities. 
c
c  Arguments
c  =========
c
c  SRNAME  (input) CHARACTER*6
c          The name of the routine which called XERBLA.
c
c  INFO    (input) INTEGER
c          The position of the invalid parameter in the parameter list
c          of the calling routine.
c
c =====================================================================
c
c     .. Executable Statements ..
c
      WRITE( *, FMT = 9999 )SRNAME, INFO
c
      stop
c
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
c
c     End of XERBLA
c
      END
c***********************************************************************
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  DGEMM  performs one of the matrix-matrix operations
c
c     C := alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X_prime,
c
c  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = "N" or "n",  op( A ) = A.
c
c              TRANSA = "T" or "t",  op( A ) = A_prime.
c
c              TRANSA = "C" or "c",  op( A ) = A_prime.
c
c           Unchanged on exit.
c
c  TRANSB - CHARACTER*1.
c           On entry, TRANSB specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSB = "N" or "n",  op( B ) = B.
c
c              TRANSB = "T" or "t",  op( B ) = B_prime.
c
c              TRANSB = "C" or "c",  op( B ) = B_prime.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANSA = "N" or "n",  and is  m  otherwise.
c           Before entry with  TRANSA = "N" or "n",  the leading  m by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. When  TRANSA = "N" or "n" then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
c           n  when  TRANSB = "N" or "n",  and is  k  otherwise.
c           Before entry with  TRANSB = "N" or "n",  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in the calling (sub) program. When  TRANSB = "N" or "n" then
c           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION.
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
c     and  columns of  A  and the  number of  rows  of  B  respectively.
c
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And if  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( NOTB )THEN
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B + beta*C.
c
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
c
c           Form  C := alpha*A_prime*B + beta*C
c
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B_prime + beta*C
c
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
c
c           Form  C := alpha*A_prime*B_prime + beta*C
c
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of DGEMM .
c
      END
c***********************************************************************
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
c     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
c     ..
c
c  Purpose
c  =======
c
c  DGER   performs the rank 1 operation
c
c     A := alpha*x*y_prime + A,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and A is an m by n matrix.
c
c  Parameters
c  ==========
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of the matrix A.
c           M must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix A.
c           N must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  X      - DOUBLE PRECISION array of dimension at least
c           ( 1 + ( m - 1 )*abs( INCX ) ).
c           Before entry, the incremented array X must contain the m
c           element vector x.
c           Unchanged on exit.
c
c  INCX   - INTEGER.
c           On entry, INCX specifies the increment for the elements of
c           X. INCX must not be zero.
c           Unchanged on exit.
c
c  Y      - DOUBLE PRECISION array of dimension at least
c           ( 1 + ( n - 1 )*abs( INCY ) ).
c           Before entry, the incremented array Y must contain the n
c           element vector y.
c           Unchanged on exit.
c
c  INCY   - INTEGER.
c           On entry, INCY specifies the increment for the elements of
c           Y. INCY must not be zero.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
c           Before entry, the leading m by n part of the array A must
c           contain the matrix of coefficients. On exit, A is
c           overwritten by the updated matrix.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. LDA must be at least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 2 Blas routine.
c
c  -- Written on 22-October-1986.
c     Jack Dongarra, Argonne National Lab.
c     Jeremy Du Croz, Nag Central Office.
c     Sven Hammarling, Nag Central Office.
c     Richard Hanson, Sandia National Labs.
c
c
c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
c     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
c
c     Start the operations. In this version the elements of A are
c     accessed sequentially with one pass through A.
c
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
c
      RETURN
c
c     End of DGER  .
c
      END
c***********************************************************************
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  DTRSM  solves one of the matrix equations
c
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c
c  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A_prime.
c
c  The matrix X is overwritten on B.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry, SIDE specifies whether op( A ) appears on the left
c           or right of X as follows:
c
c              SIDE = "L" or "l"   op( A )*X = alpha*B.
c
c              SIDE = "R" or "r"   X*op( A ) = alpha*B.
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = "U" or "u"   A is an upper triangular matrix.
c
c              UPLO = "L" or "l"   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = "N" or "n"   op( A ) = A.
c
c              TRANSA = "T" or "t"   op( A ) = A_prime.
c
c              TRANSA = "C" or "c"   op( A ) = A_prime.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = "U" or "u"   A is assumed to be unit triangular.
c
c              DIAG = "N" or "n"   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = "L" or "l"  and is  n  when  SIDE = "R" or "r".
c           Before entry  with  UPLO = "U" or "u",  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = "L" or "l",  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = "U" or "u",  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = "L" or "l"  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = "R" or "r"
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain  the  right-hand  side  matrix  B,  and  on exit  is
c           overwritten by the solution matrix  X.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*inv( A )*B.
c
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*inv( A_prime )*B.
c
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*inv( A ).
c
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*inv( A_prime ).
c
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of DTRSM .
c
      END
c***********************************************************************
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
c***********************************************************************
