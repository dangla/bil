C COPYRIGHT (c) 1995 Timothy A. Davis and
C            Council for the Central Laboratory of the Research Councils
C Original date 30 November 1995
C## Initialization of XOUT and IOUT added to subroutine MA38DD.  20/1/98
C## In cases where MA38ED is not called this avoids a possibly wrong
C## decision on storage requirements after label 9000 in MA38DD.
C## Initialization of XRMAX moved to top of subroutine MA38DD.   14/1/99
C## Small modifications to commenst on iout/xout in MA38DD and MA38ED.
C#######################################################################
C## MA38: double precision version
C## user-callable routines are listed first
C#######################################################################

C 12th July 2004 Version 1.0.0. Version numbering added.
C 21st February 2005 Version 1.1.0. FD05 dependence changed to FD15.
C      Calls to SGEMV with N=0 avoided to prevent subscript errors.

        SUBROUTINE MA38ID (KEEP, CNTL, ICNTL)
        INTEGER            KEEP(20)
        DOUBLE PRECISION   CNTL(10)
        INTEGER            ICNTL(20)

C=== MA38ID ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Initialize user-controllable parameters to default values, and
C  non-user controllable parameters.  This routine is normally
C  called once prior to any call to MA38AD.
C
C  This routine sets the default control parameters.  We recommend
C  changing these defaults under certain circumstances:
C
C  (1) If you know that your matrix has nearly symmetric nonzero
C       pattern, then we recommend setting Icntl (6) to 1 so that
C       diagonal pivoting is preferred.  This can have a significant
C       impact on the performance for matrices that are essentially
C       symmetric in pattern.
C
C   (2) If you know that your matrix is not reducible to block
C       triangular form, then we recommend setting Icntl (4) to 0
C       so that MA38 does not try to permute the matrix to block
C       triangular form (it will not do any useful work and will
C       leave the matrix in its irreducible form).  The work saved
C       is typically small, however.
C
C   The other control parameters typically have less effect on overall
C   performance.

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be modified on installation to reflect the computer
C  or environment in which this package is installed (printing control,
C  maximum integer, block size, and machine epsilon in particular).  If
C  you, the installer, modify this routine, please comment out the
C  original code, and add comments (with date) to describe the
C  installation.  Do not delete any original code.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C               --------------------------------------------------------
C  Icntl:       An integer array of size 20.  Need not be set by
C               caller on input.  On output, it contains default
C               integer control parameters.
C
C  Icntl (1):   Fortran output unit for error messages.
C               Default: 6
C
C  Icntl (2):   Fortran output unit for diagnostic messages.
C               Default: 6
C
C  Icntl (3):   printing-level.
C               0 or less: no output
C               1: error messages only
C               2: error messages and terse diagnostics
C               3: as 2, and print first few entries of all input and
C                       output arguments.  Invalid and duplicate entries
C                       are printed.
C               4: as 2, and print all entries of all input and
C                       output arguments.  Invalid and duplicate entries
C                       are printed.  The entire input matrix and its
C                       factors are printed.
C               5: as 4, and print out information on the data
C                       structures used to represent the LU factors,
C                       the assembly DAG, etc.
C               Default: 2
C
C  Icntl (4):   whether or not to attempt a permutation to block
C               triangular form.  If equal to one, then attempt the
C               permutation.  If you know the matrix is not reducible
C               to block triangular form, then setting Icntl (4) to
C               zero can save a small amount of computing time.
C               Default: 1 (attempt the permutation)
C
C  Icntl (5):   the number of columns to examine during the global
C               pivot search.  A value less than one is treated as one.
C               Default: 4
C
C  Icntl (6):   if not equal to zero, then pivots from the diagonal
C               of A (or the diagonal of the block-triangular form) are
C               preferred.  If the nonzero pattern of the matrix is
C               basically symmetric, we recommend that you change this
C               default value to 1 so that pivots on the diagonal
C               are preferred.
C               Default: 0 (do not prefer the diagonal)
C
C  Icntl (7):   block size for the BLAS, controlling the tradeoff
C               between the Level-2 and Level-3 BLAS.  Values less than
C               one are treated as one.
C               Default: 16, which is suitable for the CRAY YMP.
C
C  Icntl (8):   number of steps of iterative refinement to perform.
C               Values less than zero are treated as zero.  The matrix
C               must be preserved for iterative refinement to be done
C               (job=1 in MA38AD or MA38BD).
C               Default: 0  (no iterative refinement)
C
C  Icntl (9 ... 20):  set to zero.  Reserved for future releases.

C               --------------------------------------------------------
C  Cntl:        A double precision array of size 10.
C               Need not be set by caller on input.  On output, contains
C               default double precision control parameters.
C
C  Cntl (1):    pivoting tradeoff between sparsity-preservation
C               and numerical stability.  An entry A(k,k) is numerically
C               acceptable if:
C                  abs (A(k,k)) >= Cntl (1) * max (abs (A(*,k)))
C               Values less than zero are treated as zero (no numerical
C               constraints).  Values greater than one are treated as
C               one (partial pivoting with row interchanges).
C               Default: 0.1
C
C  Cntl (2):    amalgamation parameter.  If the first pivot in a
C               frontal matrix has row degree r and column degree c,
C               then a working array of size
C                  (Cntl (2) * c) - by - (Cntl (2) * r)
C               is allocated for the frontal matrix.  Subsequent pivots
C               within the same frontal matrix must fit within this
C               working array, or they are not selected for this frontal
C               matrix.  Values less than one are treated as one (no
C               fill-in due to amalgamation).  Some fill-in due to
C               amalgamation is necessary for efficient use of the BLAS
C               and to reduce the assembly operations required.
C               Default: 2.0
C
C  Cntl (3):    Normally not modified by the user.
C               Defines the smallest positive number,
C               epsilon = Cntl (3), such that fl (1.0 + epsilon)
C               is greater than 1.0 (fl (x) is the floating-point
C               representation of x).  If the floating-point mantissa
C               is binary, then Cntl (3) is 2 ** (-b+1), where b
C               is the number of bits in the mantissa (including the
C               implied bit, if applicable).
C
C               Cntl (3) is only
C               used in MA38JD to compute the sparse backward error
C               estimates, Rinfo (7) and Rinfo (8), when
C               Icntl (8) > 0 (the default is Icntl (8) = 0,
C               so by default, Cntl (3) is not used).
C
C  Cntl (4 ... 10):  set to zero.  Reserved for future releases.

C               --------------------------------------------------------
C  Keep:        An integer array of size 20.
C               Need not be set by the caller.  On output, contains
C               integer control parameters that are (normally) non-user
C               controllable (but can of course be modified by the
C               "expert" user or library installer).
C
C  Keep (1 ... 5):  unmodified (see MA38AD or MA38BD for a description).
C
C  Keep (6):    Formerly, this held the largest representable positive
C               integer.  MA38 no longer needs this, and it is set to
C               zero by MA38ID and ignored elsewhere.
C
C  Keep (7) and Keep (8): A column is treated as "dense" if
C               it has more than
C               max (0, Keep(7), Keep(8)*int(sqrt(float(n))))
C               original entries.  "Dense" columns are treated
C               differently that "sparse" rows and columns.  Dense
C               columns are transformed into a priori contribution
C               blocks of dimension cdeg-by-1, where cdeg is the number
C               of original entries in the column.  Modifying these two
C               parameters can change the pivot order.
C               Default:  Keep (7) = 64
C               Default:  Keep (8) = 1
C
C  Keep (9 ... 20):  set to zero.  Reserved for future releases.

C## End of user documentation for MA38ID ###############################

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine

C       functions called:       FD15AD
        DOUBLE PRECISION FD15AD
        EXTERNAL         FD15AD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I
        DOUBLE PRECISION ZERO, TENTH, TWO
        PARAMETER (TENTH = 0.1D0, TWO = 2.0D0, ZERO = 0.0D0)

C  i:       loop index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

C       ----------------------------------------------------------------
C       integer control parameters:
C       ----------------------------------------------------------------

        ICNTL (1) = 6
        ICNTL (2) = 6
        ICNTL (3) = 2
        ICNTL (4) = 1
        ICNTL (5) = 4
        ICNTL (6) = 0
        ICNTL (7) = 16
        ICNTL (8) = 0

C       Icntl (9 ... 20) is reserved for future releases:
        DO 10 I = 9, 20
           ICNTL (I) = 0
10      CONTINUE

C       ----------------------------------------------------------------
C       double precision control parameters:
C       ----------------------------------------------------------------

        CNTL (1) = TENTH
        CNTL (2) = TWO

C       Set machine epsilon by call to FD15.
        CNTL (3) = FD15AD('E')

C       Cntl (4 ... 10) is reserved for future releases:
        DO 30 I = 4, 10
           CNTL (I) = ZERO
30      CONTINUE

C       ----------------------------------------------------------------
C       integer control parameters in Keep:
C       ----------------------------------------------------------------

C       Keep(6) once held largest integer, now not used.
        KEEP (6) = 0
        KEEP (7) = 64
        KEEP (8) = 1

C       Keep (9 ... 20) is reserved for future releases:
        DO 20 I = 9, 20
           KEEP (I) = 0
20      CONTINUE

        RETURN
        END

        SUBROUTINE MA38AD (N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     *                     INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)
        INTEGER          N, NE, JOB
        LOGICAL          TRANSA
        INTEGER          LVALUE, LINDEX
        DOUBLE PRECISION VALUE(LVALUE)
        INTEGER          INDEX(LINDEX), KEEP(20)
        DOUBLE PRECISION CNTL(10)
        INTEGER          ICNTL(20), INFO(40)
        DOUBLE PRECISION RINFO(20)

C=== MA38AD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Given a sparse matrix A, find a sparsity-preserving and numerically-
C  acceptable pivot order and compute the LU factors, PAQ = LU.  The
C  matrix is optionally preordered into a block upper triangular form
C  (BTF).  Pivoting is performed within each diagonal block to maintain
C  sparsity and numerical stability.  The method used to factorize the
C  matrix is an unsymmetric-pattern variant of the multifrontal method.
C  Most of the floating-point work is done in the Level-3 BLAS (dense
C  matrix multiply).  In addition, approximate degrees are used in the
C  Markowitz-style pivot search to reduce the symbolic overhead.  For
C  best performance, be sure to use an optimized BLAS library.
C
C  This routine is normally preceded by a call to MA38ID, to
C  initialize the default control parameters.  MA38ID need only be
C  called once.  A call to MA38AD can be followed by any number of
C  calls to MA38CD, which solves a linear system using the LU factors
C  computed by this routine.  A call to MA38AD can also be followed by
C  any number of calls to MA38BD, which factorizes another matrix with
C  the same nonzero pattern as the matrix factorized by MA38AD (but with
C  different numerical values).
C
C  For more information, see T. A. Davis and I. S. Duff, "An
C  unsymmetric-pattern multifrontal method for sparse LU factorization",
C  SIAM J. Matrix Analysis and Applications (to appear), also
C  technical report TR-94-038, CISE Dept., Univ. of Florida,
C  P.O. Box 116120, Gainesville, FL 32611-6120, USA.  The method used
C  here is a modification of that method, described in T. A. Davis,
C  "A combined unifrontal/multifrontal method for unsymmetric sparse
C  matrices," TR-94-005.  (Technical reports are available via WWW at
C  http://www.cis.ufl.edu/).  The appoximate degree update algorithm
C  used here has been incorporated into an approximate minimum degree
C  ordering algorithm, desribed in P. Amestoy, T. A. Davis, and I. S.
C  Duff, "An approximate minimum degree ordering algorithm", SIAM J.
C  Matrix Analysis and Applications (to appear, also TR-94-039).  The
C  approximate minimum degree ordering algorithm is implemented as MC47
C  in the Harwell Subroutine Library (MC47 is not called by
C  MA38).

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  Requires the BLAS (Basic Linear Algebra Subprograms).  Ideally, you
C  should use vendor-optimized BLAS for your computer.  If you do not
C  have them, you may obtain the Fortran BLAS from NETLIB.  Send email
C  to netlib@ornl.gov with the two-line message:
C               send index from blas
C               send blas.shar from blas
C
C  To permamently disable any diagnostic and/or error printing, see
C  the "INSTALLATION NOTE:" comments in MA38YD and MA38ZD.
C
C  To change the default control parameters, see the
C  "INSTALLATION NOTE:" comments in MA38ID.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C           ------------------------------------------------------------
C  n:       An integer variable.
C           Must be set by caller on input (not modified).
C           Order of the matrix.

C           ------------------------------------------------------------
C  ne:      An integer variable.
C           Must be set by caller on input (not modified).
C           Number of entries in input matrix.  Restriction:  ne => 1.

C           ------------------------------------------------------------
C  job:     An integer variable.
C           Must be set by caller on input (not modified).
C           If job=1, then a column-oriented form of the input matrix
C           is preserved, otherwise, the input matrix is overwritten
C           with its LU factors.  If iterative refinement is to done
C           in MA38CD, (Icntl (8) > 0), then job must be set to 1.

C           ------------------------------------------------------------
C  transa:  A logical variable.
C           Must be set by caller on input (not modified).
C           If false then A is factorized: PAQ = LU.  Otherwise, A
C           transpose is factorized:  PA'Q = LU.

C           ------------------------------------------------------------
C  lvalue:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Value array.  Restriction:  lvalue >= 2*ne
C           is required to convert the input form of the matrix into
C           the internal represenation.  lvalue >= ne + axcopy is
C           required to start the factorization, where axcopy = ne if
C           job = 1, or axcopy = 0 otherwise.  During factorization,
C           additional memory is required to hold the frontal matrices.
C           The internal representation of the matrix is overwritten
C           with the LU factors, of size (Keep (2) - Keep (1) + 1
C           + axcopy), on output.

C           ------------------------------------------------------------
C  lindex:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Index array.  Restriction: lindex >= 3*ne+2*n+1,
C           is required to convert the input form of the matrix into
C           its internal representation.  lindex >= wlen + alen + acopy
C           is required to start the factorization, where
C           wlen <= 11*n + 3*dn + 8 is the size of the workspaces,
C           dn <= n is the number of columns with more than d
C           entries (d = max (64, sqrt (n)) is the default),
C           alen <= 2*ne + 11*n + 11*dn + dne is the size of the
C           internal representation of the matrix, dne <= ne is the
C           number of entries in such columns with more than d entries,
C           and acopy = ne+n+1 if job = 1, or acopy = 0 otherwize.
C           During factorization, the internal representation of size
C           alen is overwritten with the LU factors, of size
C           luilen = (Keep (5) - Keep (3) + 1 - acopy) on output.
C           Additional memory is also required to hold the unsymmetric
C           quotient graph, but this also overwrites the input matrix.
C           Usually about 7*n additional space is adequate for this
C           purpose.  Just prior to the end of factorization,
C           lindex >= wlen + luilen + acopy is required.

C           ------------------------------------------------------------
C  Value:   A double precision array of size lvalue.
C           Must be set by caller on input.  Modified on output.  On
C           input, Value (1..ne) holds the original matrix in triplet
C           form.  On output, Value holds the LU factors, and
C           (optionally) a column-oriented form of the original matrix
C           - otherwise the input matrix is overwritten with the LU
C           factors.

C           ------------------------------------------------------------
C  Index:   An integer array of size lindex.
C           Must be set by caller on input.  Modified on output.  On
C           input, Index (1..2*ne) holds the original matrix in triplet
C           form.  On output, Index holds the LU factors, and
C           (optionally) a column-oriented form of the original matrix
C           - otherwise the input matrix is overwritten with the LU
C           factors.
C
C           On input the kth triplet (for k = 1...ne) is stored as:
C                       A (row,col) = Value (k)
C                       row         = Index (k)
C                       col         = Index (k+ne)
C           If there is more than one entry for a particular position,
C           the values are accumulated, and the number of such duplicate
C           entries is returned in Info (2), and a warning flag is
C           set.  However, applications such as finite element methods
C           naturally generate duplicate entries which are then
C           assembled (added) together.  If this is the case, then
C           ignore the warning message.
C
C           On output, the LU factors and the column-oriented form
C           of A (if preserved) are stored in:
C               Value (Keep (1)...Keep (2))
C               Index (Keep (3)...Keep (5))
C           where Keep (2) = lvalue, and Keep (5) = lindex.

C           ------------------------------------------------------------
C  Keep:    An integer array of size 20.
C
C           Keep (1 ... 5):  Need not be set by caller on input.
C               Modified on output.
C               Keep (1): LU factors start here in Value
C               Keep (2) = lvalue: LU factors end here in Value
C               Keep (3): LU factors start here in Index
C               Keep (4): LU factors needed for MA38BD start here
C                             in Index
C               Keep (5) = lindex: LU factors end here in Index
C
C           Keep (6 ... 8):  Must be set by caller on input (not
C               modified).
C               integer control arguments not normally modified by the
C               user.  See MA38ID for details, which sets the defaults.
C               Keep (6) was the largest representable positive
C               integer but is set to zero and is no longer used.
C               Keep (7) and Keep (8) determine the
C               size of d, where columns with more than d original
C               entries are treated as a priori frontal matrices.
C
C           Keep (9 ... 20): Unused.  Reserved for future releases.

C           ------------------------------------------------------------
C  Cntl:    A double precision array of size 10.
C           Must be set by caller on input (not modified).
C           real control arguments, see MA38ID for a description,
C           which sets the defaults. MA38AD uses Cntl (1) and Cntl (2).

C           ------------------------------------------------------------
C  Icntl:   An integer array of size 20.
C           Must be set by caller on input (not modified).
C           Integer control arguments, see MA38ID for a description,
C           which sets the defaults.  MA38AD uses Icntl (1..7).

C           ------------------------------------------------------------
C  Info:    An integer array of size 40.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of MA38AD.
C
C           Info (1): zero if no error occurred, negative if
C               an error occurred and the factorization was not
C               completed, positive if a warning occurred (the
C               factorization was completed).
C
C               These errors cause the factorization to terminate:
C
C               Error   Description
C               -1      n < 1 or n > maximum value
C               -2      ne < 1 or ne > maximum value
C               -3      lindex too small
C               -4      lvalue too small
C               -5      both lindex and lvalue are too small
C
C               With these warnings the factorization was able to
C               complete:
C
C               Error   Description
C               1       invalid entries
C               2       duplicate entries
C               3       invalid and duplicate entries
C               4       singular matrix
C               5       invalid entries, singular matrix
C               6       duplicate entries, singular matrix
C               7       invalid and duplicate entries, singular matrix
C
C               Subsequent calls to MA38BD and MA38CD can only be made
C               if Info (1) is zero or positive.  If Info (1)
C               is negative, then some or all of the remaining
C               Info and Rinfo arrays may not be valid.
C
C           Info (2): duplicate entries in A.  A warning is set
C               if Info (2) > 0.  However, the duplicate entries
C               are summed and the factorization continues.  Duplicate
C               entries are sometimes intentional - for finite element
C               codes, for example.
C
C           Info (3): invalid entries in A, indices not in 1..n.
C               These entries are ignored and a warning is set
C               in Info (1).
C
C           Info (4): zero.  Used by MA38BD only.
C
C           Info (5): entries in A after adding duplicates and
C               removing invalid entries.
C
C           Info (6): entries in diagonal blocks of A.
C
C           Info (7): entries in off-diagonal blocks of A.  Zero
C               if Info (9) = 1.
C
C           Info (8): 1-by-1 diagonal blocks.
C
C           Info (9): blocks in block-triangular form.
C
C           Info (10): entries below diagonal in L.
C
C           Info (11): entries below diagonal in U.
C
C           Info (12): entries in L+U+offdiagonal part.
C
C           Info (13): frontal matrices.
C
C           Info (14): garbage collections performed on Index, when
C               memory is exhausted.  Garbage collections are performed
C               to remove external fragmentation.  If Info (14) is
C               excessively high, performance can be degraded.  Try
C               increasing lindex if that occurs.  Note that external
C               fragmentation in *both* Index and Value is removed when
C               either is exhausted.
C
C           Info (15): garbage collections performed on Value.
C
C           Info (16): diagonal pivots chosen.
C
C           Info (17): numerically acceptable pivots found in A.
C               If less than n, then A is singular (or nearly so).
C               The factorization still proceeds, and MA38CD can still
C               be called.  The zero-rank active submatrix of order
C               n - Info (17) is replaced with the identity matrix
C               (assuming BTF is not in use).  If BTF is in use, then
C               one or more of the diagonal blocks are singular.
C
C           Info (18): memory used in Index.
C
C           Info (19): minimum memory needed in Index
C               (or minimum recommended).  If lindex is set to
C               Info (19) on a subsequent call, then a moderate
C               number of garbage collections (Info (14)) will
C               occur.
C
C           Info (20): memory used in Value.
C
C           Info (21): minimum memory needed in Value
C               (or minimum recommended).  If lvalue is set to
C               Info (21) on a subsequent call, then a moderate
C               number of garbage collections (Info (15)) will
C               occur.
C
C           Info (22): memory needed in Index for the next call to
C               MA38BD.
C
C           Info (23): memory needed in Value for the next call to
C               MA38BD.
C
C           Info (24): zero.  Used by MA38CD only.
C
C           Info (25 ... 40): reserved for future releases

C           ------------------------------------------------------------
C  Rinfo:   A double precision array of size 20.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of MA38AD.
C
C           Rinfo (1): total flop count in the BLAS
C
C           Rinfo (2): total assembly flop count
C
C           Rinfo (3): total flops during pivot search
C
C           Rinfo (4): Level-1 BLAS flops
C
C           Rinfo (5): Level-2 BLAS flops
C
C           Rinfo (6): Level-3 BLAS flops
C
C           Rinfo (7): zero.  Used by MA38CD only.
C
C           Rinfo (8): zero.  Used by MA38CD only.
C
C           Rinfo (9 ... 20): reserved for future releases

C=======================================================================
C  TO BE PRESERVED BETWEEN CALLS TO MA38AD, MA38BD, MA38CD:
C=======================================================================
C
C  When calling MA38CD to solve a linear system using the factors
C  computed by MA38AD or MA38BD, the following must be preserved:
C
C       n
C       Value (Keep (1)...Keep (2))
C       Index (Keep (3)...Keep (5))
C       Keep (1 ... 20)
C
C  When calling MA38BD to factorize a subsequent matrix with a pattern
C  similar to that factorized by MA38AD, the following must be
C  preserved:
C
C       n
C       Index (Keep (4)...Keep (5))
C       Keep (4 ... 20)
C
C  Note that the user may move the LU factors to a different position
C  in Value and/or Index, as long as Keep (1 ... 5) are modified
C  correspondingly.

C## End of user documentation for MA38AD ###############################

C=======================================================================
C  CODING CONVENTIONS:
C=======================================================================
C
C  This package is written in ANSI Fortran 77.  To make the code more
C  understandable, the following coding conventions are followed for all
C  routines in this package:
C
C  1) Large code blocks are delimited with [...] comments.
C
C  2) GOTO usage:
C       a) Goto's used to return if an error condition is found are
C          written as "GO TO 9000" or "GO TO 9010".
C       b) Goto's used to exit loops prematurely are written as "GO TO",
C          and have a target label of 2000 or less.
C       c) Goto's used to jump to the next iteration of a do loop or
C          while loop (or to implement a while loop) are written as
C          "GOTO".
C       No other goto's are used in this package.
C
C  This package uses the following CRAY compiler directives to help
C  in the vectorization of loops.  Each of them operate on the
C  do-loop immediately following the directive.  Other compilers
C  normally treat these directives as ordinary comments.
C
C       CFPP$ NODEPCHK L        disables data dependency check, and
C                               asserts that no recursion exists.
C       CFPP$ NOLSTVAL L        disables the saving of last values of
C                               transformed scalars (indexes or promoted
C                               scalars, especially those in array
C                               subscripts).  Asserts that values do not
C                               need to be the same as in the scalar
C                               version (for later use of the scalars).
C       CDIR$ SHORTLOOP         asserts that the loop count is always
C                               64 or less.

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine
C       subroutines called:     MA38ND, MA38YD, MA38KD, MA38DD
C       functions called:       MAX, MIN
        INTRINSIC MAX, MIN
        EXTERNAL  MA38DD, MA38KD, MA38ND, MA38YD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, NZ, LUX1, LUI1, IUSE, XUSE, LUIR1, NZOFF, NBLKS
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)

C  Location of LU factors:
C  -----------------------
C  lux1:    real part of LU factors placed in Value (lux1 ... lvalue)
C  lui1:    integer part of LU factors placed in Index (lui1 ... lindex)
C  luir1:   Index (luir1 ... lindex) must be preserved for MA38BD
C
C  Memory usage:
C  -------------
C  iuse:    current memory usage in Index
C  xuse:    current memory usage in Value
C
C  Matrix to factorize:
C  --------------------
C  nblks:   number of diagonal blocks (1 if BTF not used)
C  nzoff:   entries in off-diagonal part (0 if BTF not used)
C  nz:      entries in matrix after removing invalid/duplicate entries
C
C  Other:
C  ------
C  nmax:    largest permissible value of n
C  i:       general loop index
C  presrv:  true if original matrix to be preserved

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  clear informational output, and Keep array (except Keep (6..8)):
C-----------------------------------------------------------------------

        DO 10 I = 1, 40
           INFO (I) = 0
10      CONTINUE
        DO 20 I = 1, 20
           RINFO (I) = ZERO
20      CONTINUE
        KEEP (1) = 0
        KEEP (2) = 0
        KEEP (3) = 0
        KEEP (4) = 0
        KEEP (5) = 0

C-----------------------------------------------------------------------
C  print input arguments if requested
C-----------------------------------------------------------------------

        CALL MA38YD (1, 1,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)

C-----------------------------------------------------------------------
C  initialize and check inputs
C-----------------------------------------------------------------------

        IUSE = 0
        XUSE = 0
        INFO (5) = NE
        INFO (6) = NE
        IF (N .LT. 1) THEN
C          n is too small
           CALL MA38ND (1, ICNTL, INFO, -1, -1)
           GO TO 9000
        ENDIF
        IF (NE .LT. 1) THEN
C          ne is too small
           CALL MA38ND (1, ICNTL, INFO, -2, -1)
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  get memory for conversion to column form
C-----------------------------------------------------------------------

        NZ = NE
        IUSE = 2*N+1 + MAX (2*NZ, N+1) + NZ
        XUSE = 2*NZ
        INFO (18) = IUSE
        INFO (20) = XUSE
        INFO (19) = IUSE
        INFO (21) = XUSE
        IF (LINDEX .LT. IUSE) THEN
C          set error flag if out of integer memory:
           CALL MA38ND (1, ICNTL, INFO, -3, IUSE)
        ENDIF
        IF (LVALUE .LT. XUSE) THEN
C          set error flag if out of real memory:
           CALL MA38ND (1, ICNTL, INFO, -4, XUSE)
        ENDIF
        IF (INFO (1) .LT. 0) THEN
C          error return, if not enough integer and/or real memory:
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  convert to column-oriented form and remove duplicates
C-----------------------------------------------------------------------

        CALL MA38KD (N, NZ, TRANSA, VALUE, LVALUE, INFO, ICNTL,
     $     INDEX, LINDEX-(2*N+1), INDEX(LINDEX-2*N), INDEX(LINDEX-N), 1)
        IF (INFO (1) .LT. 0) THEN
C          error return, if all entries invalid (nz is now 0):
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  current memory usage:
C-----------------------------------------------------------------------

C       Index (1..n+1): column pointers.  input matrix is now in
C       Index (1..nz+n+1) and Value (1..nz)
C       col pattern: Index (n+1+ Index (col) ... n+1+ Index (col+1))
C       col values:  Value (     Index (col) ...      Index (col+1))
C       at this point, nz <= ne (nz = ne if there are no invalid or
C       duplicate entries; nz < ne otherwise).

        IUSE = NZ + (N+1)
        XUSE = NZ

C-----------------------------------------------------------------------
C  factorize
C-----------------------------------------------------------------------

        PRESRV = JOB .EQ. 1
        IF (PRESRV) THEN

C          -------------------------------------------------------------
C          keep a copy of the original matrix in column-oriented form
C          -------------------------------------------------------------

C          copy column pointers (Cp (1..n+1) = Ap (1..n+1))
           IUSE = IUSE + (N+1)
CFPP$ NODEPCHK L
           DO 30 I = 1, N+1
              INDEX (NZ+N+1+I) = INDEX (I)
30         CONTINUE

           CALL MA38DD (N, NZ, INDEX (NZ+N+2),
     $          VALUE (NZ+1), LVALUE-NZ,
     $          INDEX (NZ+2*N+3), LINDEX-(NZ+2*N+2),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX (N+2), VALUE, N, NZ, KEEP, NE)
           IF (INFO (1) .LT. 0) THEN
C             error return, if MA38DD fails
              GO TO 9000
           ENDIF
C          adjust pointers to reflect Index/Value, not II/XX:
           LUX1 = LUX1 + NZ
           LUI1 = LUI1 + (NZ+2*N+2)

C          move preserved copy of A to permanent place
           LUX1 = LUX1 - NZ
           LUI1 = LUI1 - (NZ+N+1)
           DO 40 I = NZ+N+1, 1, -1
              INDEX (LUI1+I-1) = INDEX (I)
40         CONTINUE
           DO 50 I = NZ, 1, -1
              VALUE (LUX1+I-1) = VALUE (I)
50         CONTINUE

        ELSE

C          -------------------------------------------------------------
C          do not preserve the original matrix
C          -------------------------------------------------------------

           CALL MA38DD (N, NZ, INDEX,
     $          VALUE, LVALUE,
     $          INDEX (N+2), LINDEX-(N+1),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX, VALUE, 0, 1, KEEP, NE)
           IF (INFO (1) .LT. 0) THEN
C             error return, if MA38DD fails
              GO TO 9000
           ENDIF
C          adjust pointers to reflect Index/Value, not II/XX:
           LUI1 = LUI1 + (N+1)
        ENDIF

C-----------------------------------------------------------------------
C  wrap-up
C-----------------------------------------------------------------------

        IF (TRANSA) THEN
           INDEX (LINDEX-6) = 1
        ELSE
           INDEX (LINDEX-6) = 0
        ENDIF
        INDEX (LINDEX-5) = NZOFF
        INDEX (LINDEX-4) = NBLKS
        IF (PRESRV) THEN
           INDEX (LINDEX-3) = 1
        ELSE
           INDEX (LINDEX-3) = 0
        ENDIF
        INDEX (LINDEX-2) = NZ
        INDEX (LINDEX-1) = N
        INDEX (LINDEX) = NE

C       do not need preserved matrix (n+1+nz), or off-diagonal entries
C       (nzoff) for MA38BD:
        LUIR1 = LUI1
        IF (PRESRV) THEN
C          do not need preserved matrix for MA38BD
           LUIR1 = LUIR1 + N+1 + NZ
        ENDIF
        IF (NBLKS .GT. 1) THEN
C          do not need off-diagonal part for MA38BD
           LUIR1 = LUIR1 + NZOFF
        ENDIF

C       save location of LU factors
        KEEP (1) = LUX1
        KEEP (2) = LVALUE
        KEEP (3) = LUI1
        KEEP (4) = LUIR1
        KEEP (5) = LINDEX

C       update memory usage information
        IUSE = LINDEX - LUI1 + 1
        XUSE = LVALUE - LUX1 + 1
        INFO (22) = INFO (22) + (LINDEX - LUIR1 + 1)

C-----------------------------------------------------------------------
C  print the output arguments if requested, and return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (INFO (1) .LT. 0) THEN
           KEEP (1) = 0
           KEEP (2) = 0
           KEEP (3) = 0
           KEEP (4) = 0
           KEEP (5) = 0
        ENDIF

        INFO (18) = MIN (LINDEX, MAX (INFO (18), IUSE))
        INFO (20) = MIN (LVALUE, MAX (INFO (20), XUSE))

        CALL MA38YD (1, 2,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)
        RETURN
        END

        SUBROUTINE MA38BD (N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     *                     INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)
        INTEGER          N, NE, JOB
        LOGICAL          TRANSA
        INTEGER          LVALUE, LINDEX
        DOUBLE PRECISION VALUE(LVALUE)
        INTEGER          INDEX(LINDEX), KEEP(20)
        DOUBLE PRECISION CNTL(10)
        INTEGER          ICNTL(20), INFO(40)
        DOUBLE PRECISION RINFO(20)

C=== MA38BD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Given a sparse matrix A, and a sparsity-preserving and numerically-
C  acceptable pivot order and symbolic factorization, compute the LU
C  factors, PAQ = LU.  Uses the sparsity pattern and permutations from
C  a prior factorization by MA38AD or MA38BD.  The matrix A should have
C  the same nonzero pattern as the matrix factorized by MA38AD or
C  MA38BD.  The matrix can have different numerical values.  No
C  variations are made in the pivot order computed by MA38AD.  If a
C  zero pivot is encountered, an error flag is set and the
C  factorization terminates.
C
C  This routine can actually handle any matrix A such that (PAQ)_ij can
C  be nonzero only if (LU)_ij is be nonzero, where L and U are the LU
C  factors of the matrix factorized by MA38AD.  If BTF (block triangular
C  form) is used, entries above the diagonal blocks of (PAQ)_ij can have
C  an arbitrary sparsity pattern.  Entries for which (LU)_ij is not
C  present, or those below the diagonal blocks are invalid and ignored
C  (a warning flag is set and the factorization proceeds without the
C  invalid entries).  A listing of the invalid entries can be printed.
C
C  This routine must be preceded by a call to MA38AD or MA38BD.
C  A call to MA38BD can be followed by any number of calls to MA38CD,
C  which solves a linear system using the LU factors computed by this
C  routine or by MA38AD.  A call to MA38BD can also be followed by any
C  number of calls to MA38BD.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C           ------------------------------------------------------------
C  n:       An integer variable.
C           Must be set by caller on input (not modified).
C           Order of the matrix.  Must be identical to the value of n
C           in the last call to MA38AD.

C           ------------------------------------------------------------
C  ne:      An integer variable.
C           Must be set by caller on input (not modified).
C           Number of entries in input matrix.  Normally not modified
C           since the last call to MA38AD.
C           Restriction:  1 <= ne < (Keep (4)) / 2

C           ------------------------------------------------------------
C  job:     An integer variable.
C           Must be set by caller on input (not modified).
C           If job=1, then a column-oriented form of the input matrix
C           is preserved, otherwise, the input matrix is overwritten
C           with its LU factors.  If iterative refinement is to done
C           (Icntl (8) > 0), then job must be set to 1.  Can be
C           the same, or different, as the last call to MA38AD.

C           ------------------------------------------------------------
C  transa:  A logical variable.
C           Must be set by caller on input (not modified).
C           If false then A is factorized: PAQ = LU.  Otherwise, A
C           transpose is factorized:  PA'Q = LU.  Normally the same as
C           the last call to MA38AD.

C           ------------------------------------------------------------
C  lvalue:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Value array.  Restriction:  lvalue >= 2*ne,
C           although a larger will typically be required to complete
C           the factorization.  The exact value required is computed
C           by the last call to MA38AD or MA38BD (Info (23)).
C           This value assumes that the ne, job, and transa parameters
C           are the same as the last call.  Some garbage collection may
C           occur if lvalue is set to Info (23), but usually not
C           much.  We recommend lvalue => 1.2 * Info (23).  The
C           lvalue parameter is usually the same as in the last call to
C           MA38AD, however.

C           ------------------------------------------------------------
C  lindex:  An integer variable.
C           Must be set by caller on input (not modified).
C           Size of the Index array.  Restriction:
C           lindex >= 3*ne+2*n+1 + (Keep (5) - Keep (4) + 1),
C           although a larger will typically be required to complete
C           the factorization.  The exact value required is computed
C           by the last call to MA38AD or MA38BD (Info (22)).
C           This value assumes that the ne, job, and transa parameters
C           are the same as the last call.  No garbage collection ever
C           occurs in the Index array, since MA38BD does not create
C           external fragmentation in Index.  The lindex parameter is
C           usually the same as in the last call to MA38AD, however.
C           Note that lindex >= Keep (5) is also required, since
C           the pattern of the prior LU factors reside in
C           Index (Keep (4) ... Keep (5)).

C           ------------------------------------------------------------
C  Value:   A double precision array of size lvalue.
C           Must be set by caller on input (normally from the last call
C           to MA38AD or MA38BD).  Modified on output.  On input,
C           Value (1..ne) holds the original matrix in triplet form.
C           On output, Value holds the LU factors, and (optionally) a
C           column-oriented form of the original matrix - otherwise
C           the input matrix is overwritten with the LU factors.

C           ------------------------------------------------------------
C  Index:   An integer array of size lindex.
C           Must be set by caller on input (normally from the last call
C           to MA38AD or MA38BD).  Modified on output.  On input,
C           Index (1..2*ne) holds the original matrix in triplet form,
C           and Index (Keep (4) ... Keep (5)) holds the pattern
C           of the prior LU factors.  On output, Index holds the LU
C           factors, and (optionally) a column-oriented form of the
C           original matrix - otherwise the input matrix is overwritten
C           with the LU factors.
C
C           On input the kth triplet (for k = 1...ne) is stored as:
C                       A (row,col) = Value (k)
C                       row         = Index (k)
C                       col         = Index (k+ne)
C           If there is more than one entry for a particular position,
C           the values are accumulated, and the number of such duplicate
C           entries is returned in Info (2), and a warning flag is
C           set.  However, applications such as finite element methods
C           naturally generate duplicate entries which are then
C           assembled (added) together.  If this is the case, then
C           ignore the warning message.
C
C           On input, and the pattern of the prior LU factors is in
C               Index (Keep (4) ... Keep (5))
C
C           On output, the LU factors and the column-oriented form
C           of A (if preserved) are stored in:
C               Value (Keep (1)...Keep (2))
C               Index (Keep (3)...Keep (5))
C           where Keep (2) = lvalue, and Keep (5) = lindex.

C           ------------------------------------------------------------
C  Keep:    An integer array of size 20.
C
C           Keep (1 ... 3):  Need not be set by caller on input.
C               Modified on output.
C               Keep (1): new LU factors start here in Value
C               Keep (2) = lvalue: new LU factors end here in Value
C               Keep (3): new LU factors start here in Index
C
C           Keep (4 ... 5): Must be set by caller on input (normally
C               from the last call to MA38AD or MA38BD). Modified on
C               output.
C               Keep (4):  On input, the prior LU factors start here
C               in Index, not including the prior (optionally) preserved
C               input matrix, nor the off-diagonal pattern (if BTF was
C               used in the last call to MA38AD).  On output, the new
C               LU factors needed for MA38BD start here in Index.
C               Keep (5):  On input, the prior LU factors end here in
C               Index.  On output, Keep (5) is set to lindex, which
C               is where the new LU factors end in Index
C
C           Keep (6 ... 8):  Unused.  These are used by MA38AD only.
C               Future releases may make use of them, however.
C
C           Keep (9 ... 20): Unused.  Reserved for future releases.

C           ------------------------------------------------------------
C  Cntl:    A double precision array of size 10.
C           Must be set by caller on input (not modified).
C           double precision control arguments, see MA38ID for a
C           description, which sets the default values.  The current
C           version of MA38BD does not actually use Cntl.  It is
C           included to make the argument list of MA38BD the same as
C           MA38AD.  MA38BD may use Cntl in future releases.

C           ------------------------------------------------------------
C  Icntl:   An integer array of size 20.
C           Must be set by caller on input (not modified).
C           Integer control arguments, see MA38ID for a description,
C           which sets the default values.  MA38BD uses Icntl (1),
C           Icntl (2), Icntl (3), and Icntl (7).

C           ------------------------------------------------------------
C  Info:    An integer array of size 40.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of MA38BD.
C
C           Info (1): zero if no error occurred, negative if
C               an error occurred and the factorization was not
C               completed, positive if a warning occurred (the
C               factorization was completed).
C
C               These errors cause the factorization to terminate:
C
C               Error   Description
C               -1      n < 1 or n > maximum value
C               -2      ne < 1 or ne > maximum value
C               -3      lindex too small
C               -4      lvalue too small
C               -5      both lindex and lvalue are too small
C               -6      prior pivot ordering no longer acceptable
C               -7      LU factors are uncomputed, or are corrupted
C
C               With these warnings the factorization was able to
C               complete:
C
C               Error   Description
C               1       invalid entries
C               2       duplicate entries
C               3       invalid and duplicate entries
C               4       singular matrix
C               5       invalid entries, singular matrix
C               6       duplicate entries, singular matrix
C               7       invalid and duplicate entries, singular matrix
C
C               Subsequent calls to MA38BD and MA38CD can only be made
C               if Info (1) is zero or positive.  If Info (1)
C               is negative, then some or all of the remaining
C               Info and Rinfo arrays may not be valid.
C
C           Info (2): duplicate entries in A.  A warning is set
C               if Info (2) > 0.  However, the duplicate entries
C               are summed and the factorization continues.  Duplicate
C               entries are sometimes intentional - for finite element
C               codes, for example.
C
C           Info (3): invalid entries in A, indices not in 1..n.
C               These entries are ignored and a warning is set in
C               Info (1).
C
C           Info (4): invalid entries in A, not in prior LU
C               factors.  These entries are ignored and a warning is
C               set in Info (1).
C
C           Info (5): entries in A after adding duplicates and
C               removing invalid entries.
C
C           Info (6): entries in diagonal blocks of A.
C
C           Info (7): entries in off-diagonal blocks of A.  Zero
C               if Info (9) = 1.
C
C           Info (8): 1-by-1 diagonal blocks.
C
C           Info (9): blocks in block-triangular form.
C
C           Info (10): entries below diagonal in L.
C
C           Info (11): entries below diagonal in U.
C
C           Info (12): entries in L+U+offdiagonal part.
C
C           Info (13): frontal matrices.
C
C           Info (14): zero.  Used by MA38AD only.
C
C           Info (15): garbage collections performed on Value.
C
C           Info (16): diagonal pivots chosen.
C
C           Info (17): numerically acceptable pivots found in A.
C               If less than n, then A is singular (or nearly so).
C               The factorization still proceeds, and MA38CD can still
C               be called.  The zero-rank active submatrix of order
C               n - Info (17) is replaced with the identity matrix
C               (assuming BTF is not in use).  If BTF is in use, then
C               one or more of the diagonal blocks are singular.
C               MA38BD can be called if the value of Info (17)
C               returned by MA38AD was less than n, but the order
C               (n - Info (17)) active submatrix is still replaced
C               with the identity matrix.  Entries residing in this
C               submatrix are ignored, their number is included in
C               Info (4), and a warning is set in Info (1).
C
C           Info (18): memory used in Index.
C
C           Info (19): memory needed in Index (same as Info (18)).
C
C           Info (20): memory used in Value.
C
C           Info (21): minimum memory needed in Value
C               (or minimum recommended).  If lvalue is set to
C               Info (21) on a subsequent call, then a moderate
C               number of garbage collections (Info (15)) will
C               occur.
C
C           Info (22): memory needed in Index for the next call to
C               MA38BD.
C
C           Info (23): memory needed in Value for the next call to
C               MA38BD.
C
C           Info (24): zero.  Used by MA38CD only.
C
C           Info (25 ... 40): reserved for future releases

C           ------------------------------------------------------------
C  Rinfo:   A double precision array of size 20.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of MA38BD.
C
C           Rinfo (1): total flop count in the BLAS
C
C           Rinfo (2): total assembly flop count
C
C           Rinfo (3): zero.  Used by MA38AD only.
C
C           Rinfo (4): Level-1 BLAS flops
C
C           Rinfo (5): Level-2 BLAS flops
C
C           Rinfo (6): Level-3 BLAS flops
C
C           Rinfo (7): zero.  Used by MA38CD only.
C
C           Rinfo (8): zero.  Used by MA38CD only.
C
C           Rinfo (9 ... 20): reserved for future releases

C=======================================================================
C  TO BE PRESERVED BETWEEN CALLS TO MA38AD, MA38BD, MA38CD:
C=======================================================================
C
C  When calling MA38CD to solve a linear system using the factors
C  computed by MA38AD or MA38BD, the following must be preserved:
C
C       n
C       Value (Keep (1)...Keep (2))
C       Index (Keep (3)...Keep (5))
C       Keep (1 ... 20)
C
C  When calling MA38BD to factorize a subsequent matrix with a pattern
C  similar to that factorized by MA38AD, the following must be
C  preserved:
C
C       n
C       Index (Keep (4)...Keep (5))
C       Keep (4 ... 20)
C
C  Note that the user may move the LU factors to a different position
C  in Value and/or Index, as long as Keep (1 ... 5) are modified
C  correspondingly.

C## End of user documentation for MA38BD ###############################

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine
C       subroutines called:     MA38ND, MA38YD, MA38KD, MA38PD
C       functions called:       MAX, MIN
        INTRINSIC MAX, MIN
        EXTERNAL MA38KD, MA38ND, MA38PD, MA38YD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, NZ, LUX1, LUI1, IUSE, XUSE, N1, NZ1, NBLKS,
     $          LIND2, LUIR1, LUSIZ, LUI2, RPERMP, CPERMP,
     $          OFFPP, LUBLPP, BLKPP, ON, NZOFF, IP2, IO, PRL
        LOGICAL PRESRV, BADLU
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C
C  Matrix to factorize:
C  --------------------
C  nz:      number of entries, after removing invalid/duplicate entries
C  presrv:  true if original matrix to be preserved
C
C  Memory usage:
C  -------------
C  iuse:    current memory usage in Index
C  xuse:    current memory usage in Value
C  lind2:   allocatable part of Index is (1..lind2)
C
C  Location and status of LU factors:
C  ----------------------------------
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for this call to MA38BD
C  lusiz:   size of Index (luir1..lui2), needed from prior LU factors
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors in Value (lux1...lvalue)
C  ip2:     pointer into trailing part of LU factors in Index
C  badlu:   if true, then LU factors are corrupted or not computed
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  on:      size of Offp array
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  nblks:   number of diagonal blocks
C  nz1:     number of entries when prior matrix factorize
C  n1:      N argument in MA38AD or MA38BD when prior matrix factorized
C
C  Other:
C  ------
C  i:       loop index

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)

C-----------------------------------------------------------------------
C  clear informational output, and the unneeded part of the Keep array
C-----------------------------------------------------------------------

        DO 10 I = 1, 40
           INFO (I) = 0
10      CONTINUE
        DO 20 I = 1, 20
           RINFO (I) = ZERO
20      CONTINUE
        KEEP (1) = 0
        KEEP (2) = 0
        KEEP (3) = 0

C-----------------------------------------------------------------------
C  print input arguments if requested
C-----------------------------------------------------------------------

        CALL MA38YD (2, 1,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)

C-----------------------------------------------------------------------
C  check input arguments
C-----------------------------------------------------------------------

        IUSE = 0
        XUSE = 0
        INFO (5) = NE
        INFO (6) = NE
        IF (N .LT. 1) THEN
C          n is too small
           CALL MA38ND (2, ICNTL, INFO, -1, -1)
           GO TO 9000
        ENDIF
        IF (NE .LT. 1) THEN
C          ne is too small
           CALL MA38ND (2, ICNTL, INFO, -2, -1)
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  get pointers to integer part of prior LU factors
C-----------------------------------------------------------------------

        LUIR1 = KEEP (4)
        LUI2 = KEEP (5)
        LUSIZ = LUI2 - LUIR1 + 1

        BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR. LUI2.GT.LINDEX
        IF (BADLU) THEN
           CALL MA38ND (2, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF
        IF (2*NE .GT. LUIR1) THEN
           CALL MA38ND (2, ICNTL, INFO, -2, LUIR1/2)
C          error return, ne is too large:
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  shift the prior LU factors down to the end of Index.  If Keep and
C  lindex are unmodified from the prior call to MA38AD, then
C  Keep (5) = lindex, and this shift is not performed.
C-----------------------------------------------------------------------

        IF (LUI2 .LT. LINDEX) THEN
           DO 30 I = LINDEX, LINDEX - LUSIZ + 1, -1
              INDEX (I) = INDEX (I - LINDEX + LUI2)
30         CONTINUE
           LUIR1 = LINDEX - LUSIZ + 1
           KEEP (5) = LINDEX
           KEEP (4) = LUIR1
        ENDIF

C-----------------------------------------------------------------------
C  get seven scalars (transa, nzoff, nblks, presrv, nz, n, ne) from LU
C-----------------------------------------------------------------------

C       ne1 = Index (lindex), not required for MA38BD
        N1 = INDEX (LINDEX-1)
        NZ1 = INDEX (LINDEX-2)
C       presr1 = Index (lindex-3) .ne. 0, not required for MA38BD
        NBLKS = INDEX (LINDEX-4)
C       nzoff1 = Index (lindex-5), not required for MA38BD
C       trans1 = Index (lindex-6) .ne. 0, not required for MA38BD

C-----------------------------------------------------------------------
C  get pointers to permutation vectors
C-----------------------------------------------------------------------

        RPERMP = (LINDEX-6) - N
        CPERMP = RPERMP - N
        IP2 = CPERMP - 1

C-----------------------------------------------------------------------
C  get pointers to block-triangular information, if BTF was used
C-----------------------------------------------------------------------

        IF (NBLKS .GT. 1) THEN

C          -------------------------------------------------------------
C          get pointers to BTF arrays
C          -------------------------------------------------------------

           OFFPP = CPERMP - (N+1)
           BLKPP = OFFPP - (NBLKS+1)
           LUBLPP = BLKPP - (NBLKS)
           IP2 = LUBLPP - 1
           ON = N

        ELSE

C          -------------------------------------------------------------
C          matrix was factorized as a single block, pass dummy arg.
C          -------------------------------------------------------------

           OFFPP = 1
           BLKPP = 1
           LUBLPP = 1
           ON = 0

        ENDIF

        BADLU = N .NE. N1 .OR. NZ1 .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
        IF (BADLU) THEN
           CALL MA38ND (2, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  get memory for conversion to column form
C-----------------------------------------------------------------------

        NZ = NE
        IUSE = 2*N+1 + MAX (2*NZ,N+1) + NZ + LUSIZ
        XUSE = 2*NZ
        INFO (18) = IUSE
        INFO (20) = XUSE
        INFO (19) = IUSE
        INFO (21) = XUSE
        INFO (23) = XUSE
        LIND2 = LUIR1 - 1
        IF (LINDEX .LT. IUSE) THEN
C          set error flag if out of integer memory
           CALL MA38ND (2, ICNTL, INFO, -3, IUSE)
        ENDIF
        IF (LVALUE .LT. XUSE) THEN
C          set error flag if out of real memory
           CALL MA38ND (2, ICNTL, INFO, -4, XUSE)
        ENDIF
        IF (INFO (1) .LT. 0) THEN
C          error return, if not enough integer and/or real memory
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  convert to column-oriented form and remove duplicates
C-----------------------------------------------------------------------

        CALL MA38KD (N, NZ, TRANSA, VALUE, LVALUE, INFO, ICNTL,
     $     INDEX, LIND2-(2*N+1), INDEX (LIND2-2*N), INDEX (LIND2-N), 2)
        IF (INFO (1) .LT. 0) THEN
C          error return, if all entries are invalid (nz is now 0)
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  current memory usage:
C-----------------------------------------------------------------------

C       Index (1..n+1): column pointers.  input matrix is now in
C       Index (1..nz+n+1) and Value (1..nz)
C       col pattern: Index (n+1+ Index (col) ... n+1+ Index (col+1))
C       col values:  Value (     Index (col) ...      Index (col+1))
C       at this point, nz <= ne (nz = ne if there are no invalid or
C       duplicate entries; nz < ne otherwise).
C       Pattern of prior LU factors and BTF arrays are in
C       Index (Keep (4) ... Keep (5))

        IUSE = NZ + (N+1) + LUSIZ
        XUSE = NZ

C-----------------------------------------------------------------------
C  refactorize
C-----------------------------------------------------------------------

        PRESRV = JOB .EQ. 1
        IF (PRESRV) THEN

C          -------------------------------------------------------------
C          keep a copy of the original matrix in column-oriented form
C          -------------------------------------------------------------

C          copy column pointers (Cp (1..n+1) = Ap (1..n+1))
           IUSE = IUSE + (N+1)
CFPP$ NODEPCHK L
           DO 40 I = 1, N+1
              INDEX (NZ+N+1+I) = INDEX (I)
40         CONTINUE

           CALL MA38PD (N, NZ, INDEX (NZ+N+2),
     $          VALUE (NZ+1), LVALUE-NZ,
     $          INDEX (NZ+2*N+3), LIND2-(NZ+2*N+2),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX (N+2), VALUE, N, NZ,
     $          INDEX (LUIR1), IP2 - LUIR1 + 1,
     $          INDEX (LUBLPP), INDEX (BLKPP), INDEX (OFFPP), ON,
     $          INDEX (CPERMP), INDEX (RPERMP), NE)
           IF (INFO (1) .LT. 0) THEN
C             error return, if MA38PD fails
              GO TO 9000
           ENDIF
C          adjust pointers to reflect Index/Value, not II/XX:
           LUX1 = LUX1 + NZ
           LUI1 = LUI1 + (NZ+2*N+2)

C          move preserved copy of A to permanent place
           LUX1 = LUX1 - (NZ)
           LUI1 = LUI1 - (NZ+N+1)
           DO 50 I = NZ+N+1, 1, -1
              INDEX (LUI1+I-1) = INDEX (I)
50         CONTINUE
           DO 60 I = NZ, 1, -1
              VALUE (LUX1+I-1) = VALUE (I)
60         CONTINUE

        ELSE

C          -------------------------------------------------------------
C          do not preserve the original matrix
C          -------------------------------------------------------------

           CALL MA38PD (N, NZ, INDEX,
     $          VALUE, LVALUE,
     $          INDEX (N+2), LIND2-(N+1),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX, VALUE, 0, 1,
     $          INDEX (LUIR1), IP2 - LUIR1 + 1,
     $          INDEX (LUBLPP), INDEX (BLKPP), INDEX (OFFPP), ON,
     $          INDEX (CPERMP), INDEX (RPERMP), NE)
           IF (INFO (1) .LT. 0) THEN
C             error return, if MA38PD fails
              GO TO 9000
           ENDIF
C          adjust pointers to reflect Index/Value, not II/XX:
           LUI1 = LUI1 + (N+1)

        ENDIF

C-----------------------------------------------------------------------
C  wrap-up
C-----------------------------------------------------------------------

        IF (TRANSA) THEN
           INDEX (LINDEX-6) = 1
        ELSE
           INDEX (LINDEX-6) = 0
        ENDIF
        INDEX (LINDEX-5) = NZOFF
        INDEX (LINDEX-4) = NBLKS
        IF (PRESRV) THEN
           INDEX (LINDEX-3) = 1
        ELSE
           INDEX (LINDEX-3) = 0
        ENDIF
        INDEX (LINDEX-2) = NZ
        INDEX (LINDEX-1) = N
        INDEX (LINDEX) = NE

C       save location of LU factors
        KEEP (1) = LUX1
        KEEP (2) = LVALUE
        KEEP (3) = LUI1
        KEEP (4) = LUIR1
        KEEP (5) = LINDEX

C       update memory usage information
        IUSE = LINDEX - LUI1 + 1
        XUSE = LVALUE - LUX1 + 1

C-----------------------------------------------------------------------
C  print the output arguments if requested, and return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (INFO (1) .LT. 0) THEN
           KEEP (1) = 0
           KEEP (2) = 0
           KEEP (3) = 0
           KEEP (4) = 0
           KEEP (5) = 0
        ENDIF

        INFO (18) = MIN (LINDEX, MAX (INFO (18), IUSE))
        INFO (19) = INFO (18)
        INFO (22) = INFO (19)
        INFO (20) = MIN (LVALUE, MAX (INFO (20), XUSE))

        CALL MA38YD (2, 2,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)
        RETURN
        END

        SUBROUTINE MA38CD (N, JOB, TRANSC, LVALUE, LINDEX, VALUE, INDEX,
     *                     KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)
        INTEGER          N, JOB
        LOGICAL          TRANSC
        INTEGER          LVALUE, LINDEX
        DOUBLE PRECISION VALUE(LVALUE)
        INTEGER          INDEX(LINDEX), KEEP(20)
        DOUBLE PRECISION B(N), X(N), W(*), CNTL(10)
        INTEGER          ICNTL(20), INFO(40)
        DOUBLE PRECISION RINFO(20)

C=== MA38CD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Given LU factors computed by MA38AD or MA38BD, and the
C  right-hand-side, B, solve a linear system for the solution X.
C
C  This routine handles all permutations, so that B and X are in terms
C  of the original order of the matrix, A, and not in terms of the
C  permuted matrix.
C
C  If iterative refinement is done, then the residual is returned in W,
C  and the sparse backward error estimates are returned in
C  Rinfo (7) and Rinfo (8).  The computed solution X is the
C  exact solution of the equation (A + dA)x = (b + db), where
C    dA (i,j)  <= max (Rinfo (7), Rinfo (8)) * abs (A(i,j))
C  and
C    db (i) <= max (Rinfo (7) * abs (b (i)),
C                   Rinfo (8) * maxnorm (A) * maxnorm (x computed))
C  Note that dA has the same sparsity pattern as A.
C  The method used to compute the sparse backward error estimate is
C  described in M. Arioli, J. W. Demmel, and I. S. Duff, "Solving
C  sparse linear systems with sparse backward error," SIAM J. Matrix
C  Analysis and Applications, vol 10, 1989, pp. 165-190.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C           ------------------------------------------------------------
C  n:       An integer variable.
C           Must be set by caller on input (not modified).
C           Must be the same as passed to MA38AD or MA38BD.

C           ------------------------------------------------------------
C  job:     An integer variable.
C           Must be set by caller on input (not modified).
C           What system to solve (see the transc argument below).
C           Iterative refinement is only performed if job = 0,
C           Icntl (8) > 0, and only if the original matrix was
C           preserved (job = 1 in MA38AD or MA38BD).

C           ------------------------------------------------------------
C  transc:  A logical variable.
C           Must be set by caller on input (not modified).
C           solve with L and U factors or with L' and U', where
C           transa was passed to MA38AD or MA38BD.
C
C           If transa = false, then PAQ = LU was performed,
C           and the following systems are solved:
C
C                               transc = false          transc = true
C                               ----------------        ----------------
C                  job = 0      solve Ax = b            solve A'x = b
C                  job = 1      solve P'Lx = b          solve L'Px = b
C                  job = 2      solve UQ'x = b          solve QU'x = b
C
C           If transa = true, then A was transformed prior to LU
C           factorization, and P(A')Q = LU
C
C                               transc = false          transc = true
C                               ----------------        ----------------
C                  job = 0      solve A'x = b           solve Ax = b
C                  job = 1      solve P'Lx = b          solve L'Px = b
C                  job = 2      solve UQ'x = b          solve QU'x = b
C
C           Other values of job are treated as zero.  Iterative
C           refinement can be done only when solving Ax=b or A'x=b.
C
C           The comments below use Matlab notation, where
C           x = L \ b means x = (L^(-1)) * b, premultiplication by
C           the inverse of L.

C           ------------------------------------------------------------
C  lvalue:  An integer variable.
C           Must be set by caller on input (not modified).
C           The size of Value.

C           ------------------------------------------------------------
C  lindex:  An integer variable.
C           Must be set by caller on input (not modified).
C           The size of Index.

C           ------------------------------------------------------------
C  Value:   A double precision array of size lvalue.
C           Must be set by caller on input (normally from last call to
C           MA38AD or MA38BD) (not modified).
C           The LU factors, in Value (Keep (1) ... Keep (2)).
C           The entries in Value (1 ... Keep (1) - 1) and in
C           Value (Keep (2) + 1 ... lvalue) are not accessed.

C           ------------------------------------------------------------
C  Index:   An integer array of size lindex.
C           Must be set by caller on input (normally from last call to
C           MA38AD or MA38BD) (not modified).
C           The LU factors, in Index (Keep (3) ... Keep (5)).
C           The entries in Index (1 ... Keep (3) - 1) and in
C           Index (Keep (5) + 1 ... lindex) are not accessed.

C           ------------------------------------------------------------
C  Keep:    An integer array of size 20.
C
C           Keep (1..5): Must be set by caller on input (normally from
C               last call to MA38AD or MA38BD) (not modified).
C               Layout of the LU factors in Value and Index

C           ------------------------------------------------------------
C  B:       A double precision array of size n.
C           Must be set by caller on input (not modified).
C           The right hand side, b, of the system to solve.

C           ------------------------------------------------------------
C  W:       A double precision array of size 2*n or 4*n.
C           Need not be set by caller on input.  Modified on output.
C           Workspace of size W (1..2*n) if Icntl (8) = 0, which
C           is the default value.  If iterative refinement is
C           performed, and W must be of size W (1..4*n) and the
C           residual b-Ax (or b-A'x) is returned in W (1..n).

C           ------------------------------------------------------------
C  X:       A double precision array of size n.
C           Need not be set by caller on input.  Modified on output.
C           The solution, x, of the system that was solved.  Valid only
C           if Info (1) is greater than or equal to 0.

C           ------------------------------------------------------------
C  Cntl:    A double precision array of size 10.
C           Must be set by caller on input (not modified).
C           real control parameters, see MA38ID for a description,
C           which sets the defaults.

C           ------------------------------------------------------------
C  Icntl:   An integer array of size 20.
C           Must be set by caller on input (not modified).
C           Integer control parameters, see MA38ID for a description,
C           which sets the defaults.  In particular, Icntl (8) is
C           the maximum number of steps of iterative refinement to be
C           performed.

C           ------------------------------------------------------------
C  Info:    An integer array of size 40.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of MA38CD.
C
C           Info (1) is the error flag.  If Info (1) is -7, then
C           the LU factors are uncomputed, or have been corrupted since
C           the last call to MA38AD or MA38BD.  No system is solved,
C           and X (1..n) is not valid on output.  If Info (1) is 8,
C           then iterative refinement was requested but cannot be done.
C           To perform iterative refinement, the original matrix must be
C           preserved (job = 1 in MA38AD or MA38BD) and Ax=b or A'x=b
C           must be solved (job = 0 in MA38CD).  Info (24) is the
C           steps of iterative refinement actually taken.

C           ------------------------------------------------------------
C  Rinfo:   A double precision array of size 20.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of MA38CD.
C
C           If iterative refinement was performed then
C           Rinfo (7) is the sparse error estimate, omega1, and
C           Rinfo (8) is the sparse error estimate, omega2.

C=======================================================================
C  TO BE PRESERVED BETWEEN CALLS TO MA38AD, MA38BD, MA38CD:
C=======================================================================
C
C  The following must be unchanged since the call to MA38AD or MA38BD
C  that computed the LU factors:
C
C       n
C       Value (Keep (1) ... Keep (2))
C       Index (Keep (3) ... Keep (5))
C       Keep (1 ... 20)

C## End of user documentation for MA38CD ###############################

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine
C       subroutines called:     MA38ND, MA38YD, MA38JD
C       functions called:       MAX
        INTRINSIC MAX
        EXTERNAL   MA38JD, MA38ND, MA38YD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER NBLKS, OFFIP, OFFXP, N1, NZ, NE, OFFPP, BLKPP, LUBLPP,
     $          APP, AN, ANZ, ON, LUI1, LUI2, LUX1, LUX2, AIP, AXP,
     $          CPERMP, RPERMP, NZOFF, IRSTEP, YP, LY, LW, SP, IP1, IP2,
     $          XP1, LUIR1, IO, PRL
        LOGICAL PRESRV, BADLU
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C
C  Location and status of LU factors:
C  ----------------------------------
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for a call to MA38BD
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors start in Value (lux1...)
C  lux2:    real part of LU factors end in Value (...lux1)
C  ip1:     pointer into leading part of LU factors in Index
C  ip2:     pointer into trailing part of LU factors in Index
C  xp1:     pointer into leading part of LU factors in Value
C  badlu:   if true, then LU factors are corrupted or not computed
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  app:     Ap (1..n+1) array located in Index (app...app+n)
C  axp:     Ax (1..nz) array located in Value (axp...axp+nz-1)
C  aip:     Ai (1..nz) array located in Index (aip...aip+nz-1)
C  an:      n if A is preserved, 1 otherwise
C  anz:     nz if A is preserved, 1 otherwise
C  offip:   Offi (1..nzoff) array loc. in Index (offip...offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array loc. in Value (offxp...offxp+nzoff-1)
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  on:      size of Offp (1..n+1):  n if nblks > 1, 1 otherwise
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  ...      seven scalars in Index (lui2-6...lui2):
C  nzoff:   number of entries in off-diagonal part
C  nblks:   number of diagonal blocks
C  presrv:  true if original matrix was preserved when factorized
C  nz:      entries in A
C  n1:      N argument in MA38AD or MA38BD when matrix factorized
C  ne:      NE argument in MA38AD or MA38BD when matrix factorized
C
C  Arrays allocated from W work array:
C  -----------------------------------
C  lw:      size of W
C  yp:      Y (1..n) located in W (yp...yp+n-1)
C  sp:      S (1..n) located in W (sp...sp+n-1)
C  ly:      size of Y and S
C
C  Other:
C  ------
C  irstep:  maximum number of iterative refinement steps to take

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)

C-----------------------------------------------------------------------
C  clear informational output
C-----------------------------------------------------------------------

        INFO (1) = 0
        INFO (24) = 0
        RINFO (7) = ZERO
        RINFO (8) = ZERO

C-----------------------------------------------------------------------
C  print input arguments if requested
C-----------------------------------------------------------------------

        IRSTEP = MAX (0, ICNTL (8))
        IF (IRSTEP .EQ. 0) THEN
           LW = 2*N
        ELSE
           LW = 4*N
        ENDIF
        CALL MA38YD (3, 1,
     $          N, NE, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, N, W, LW)

C-----------------------------------------------------------------------
C  get pointers to LU factors
C-----------------------------------------------------------------------

        LUX1 = KEEP (1)
        LUX2 = KEEP (2)
        LUI1 = KEEP (3)
        LUIR1 = KEEP (4)
        LUI2 = KEEP (5)
        BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1
     $     .OR. LUI2 .GT. LINDEX
     $     .OR. LUX1 .LE. 0 .OR. LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE
     $     .OR. LUI1 .LE. 0 .OR. LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
        IF (BADLU) THEN
           CALL MA38ND (3, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  get seven scalars (transa, nzoff, nblks, presrv, nz, n, ne) from LU
C-----------------------------------------------------------------------

        NE = INDEX (LUI2)
        N1 = INDEX (LUI2-1)
        NZ = INDEX (LUI2-2)
        PRESRV = INDEX (LUI2-3) .NE. 0
        NBLKS = INDEX (LUI2-4)
        NZOFF = INDEX (LUI2-5)
C       transa = Index (lui2-6) .ne. 0, we don't actually need this here

C-----------------------------------------------------------------------
C  get pointers to permutation vectors
C-----------------------------------------------------------------------

        RPERMP = (LUI2-6) - N
        CPERMP = RPERMP - N
        IP2 = CPERMP - 1
        XP1 = LUX1
        IP1 = LUI1

C-----------------------------------------------------------------------
C  get pointers to preserved column-oriented copy of input matrix
C-----------------------------------------------------------------------

        IF (PRESRV) THEN

C          -------------------------------------------------------------
C          original matrix preserved in Index (lui1..lui1+nz+n) and
C          Value (lux1..lux1+nz-1)
C          -------------------------------------------------------------

           APP = IP1
           AIP = APP + N+1
           IP1 = AIP + NZ
           AXP = XP1
           XP1 = AXP + NZ
           AN = N
           ANZ = NZ

        ELSE

C          -------------------------------------------------------------
C          original matrix not preserved, pass dummy argument to MA38JD
C          -------------------------------------------------------------

           APP = 1
           AIP = 1
           AXP = 1
           AN = 1
           ANZ = 1

        ENDIF

C-----------------------------------------------------------------------
C  get pointers to block-triangular information, if BTF was used
C-----------------------------------------------------------------------

        IF (NBLKS .GT. 1) THEN

C          -------------------------------------------------------------
C          get pointers to off-diagonal nonzeros, and BTF arrays
C          -------------------------------------------------------------

           OFFIP = IP1
           IP1 = IP1 + NZOFF
           OFFXP = XP1
           XP1 = XP1 + NZOFF
           OFFPP = CPERMP - (N+1)
           BLKPP = OFFPP - (NBLKS+1)
           LUBLPP = BLKPP - (NBLKS)
           IP2 = LUBLPP - 1
           ON = N

        ELSE

C          -------------------------------------------------------------
C          matrix was factorized as a single block, pass dummy arg.
C          -------------------------------------------------------------

           OFFIP = 1
           OFFXP = 1
           OFFPP = 1
           BLKPP = 1
           LUBLPP = 1
           ON = 1

        ENDIF

        BADLU = N .NE. N1 .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $     NBLKS .LE. 0 .OR. NBLKS .GT. N .OR.
     $     XP1 .GT. LUX2 .OR. NZOFF .LT. 0 .OR. IP1 .NE. LUIR1
        IF (BADLU) THEN
           CALL MA38ND (3, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  get the number of steps of iterative refinement
C-----------------------------------------------------------------------

        IF (IRSTEP .GT. 0 .AND. .NOT. PRESRV) THEN
C          original matrix not preserved (MA38AD/MA38BD job .ne. 1)
           CALL MA38ND (3, ICNTL, INFO, 8, 0)
           IRSTEP = 0
        ENDIF
        IF (IRSTEP .GT. 0 .AND. (JOB .EQ. 1 .OR. JOB .EQ. 2)) THEN
C          iterative refinement for Ax=b and A'x=b only (job = 0)
           CALL MA38ND (3, ICNTL, INFO, 8, 1)
           IRSTEP = 0
        ENDIF
        IF (IRSTEP .EQ. 0) THEN
C          pass a dummy argument as Y, which is not accessed in MA38JD
           YP = 1
           LY = 1
           SP = 1
           LW = 2*N
        ELSE
C          pass W (yp ... yp+n-1) as Y (1..n) to MA38JD
           YP = 2*N+1
           LY = N
           SP = 3*N+1
           LW = 4*N
        ENDIF

C-----------------------------------------------------------------------
C  solve; optional iterative refinement and sparse backward error
C-----------------------------------------------------------------------

        CALL MA38JD (N, JOB, TRANSC, LUX2-XP1+1, VALUE (XP1),
     $     IP2-LUIR1+1, INDEX (LUIR1), B, X,
     $     W, W (N+1), LY, W (YP), W (SP),
     $     CNTL, INFO, RINFO, INDEX (CPERMP), INDEX (RPERMP),
     $     AN, ANZ, INDEX (APP), INDEX (AIP), VALUE (AXP),
     $     ON, MAX (1, NZOFF), INDEX (OFFPP), INDEX (OFFIP),
     $     VALUE (OFFXP), NBLKS, INDEX (LUBLPP), INDEX (BLKPP), IRSTEP)

C-----------------------------------------------------------------------
C  print output arguments if requested
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        CALL MA38YD (3, 2,
     $          N, NE, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, N, W, LW)
        RETURN
        END


C#######################################################################
C##  non-user-callable routines:
C#######################################################################


        SUBROUTINE MA38DD (N, NZ, CP, XX, XSIZE, II, ISIZE, XTAIL,
     *          ITAIL, IUSE, XUSE, NZOFF, NBLKS, ICNTL, CNTL, INFO,
     *          RINFO, PRESRV, AP, AI, AX, AN, ANZ, KEEP, NE)
        INTEGER N, NZ, ISIZE, II(ISIZE), ICNTL(20), INFO(40),
     *          CP(N+1), XSIZE, XTAIL, ITAIL, IUSE, XUSE, AN, ANZ,
     *          AP(AN+1), AI(ANZ), KEEP(20), NZOFF, NBLKS, NE
        LOGICAL PRESRV
        DOUBLE PRECISION XX(XSIZE), CNTL(10), RINFO(20), AX(ANZ)

C=== MA38DD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Factorize an unsymmetric sparse matrix in column-form, optionally
C  permuting the matrix to upper block triangular form and factorizing
C  each diagonal block.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       nz:             entries in matrix, after removing duplicates
C                       and invalid entries.
C       ne:             number of triplets, unchanged from MA38AD
C       Cp (1..n+1):    column pointers of input matrix
C       presrv:         if true then preserve original matrix
C       xsize:          size of XX
C       isize:          size of II
C       iuse:           memory usage in Index on input
C       xuse:           memory usage in Value on input
C       Icntl:          integer control parameters, see MA38ID
C       Cntl:           real control parameters, see MA38ID
C       Keep (6..8):    integer control parameters, see MA38ID
C
C       if presrv is true:
C           an:                 = n, order of preserved matrix
C           anz:                = anz, order of preserved matrix
C           Ap (1..an+1):       column pointers of preserved matrix
C           Ai (1..nz):         row indices of preserved matrix
C           Ax (1..nz):         values of preserved matrix
C           II:                 unused on input
C           XX:                 unused on input
C       else
C           an:                 1
C           anz:                1
C           Ap:                 unused
C           Ai:                 unused
C           Ax:                 unused
C           II (1..nz):         row indices of input matrix
C           XX (1..nz):         values of input matrix

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       XX (xtail ... xsize), xtail:
C
C                       LU factors are located in XX (xtail ... xsize),
C                       including values in off-diagonal part if matrix
C                       was permuted to block triangular form.
C
C       II (itail ... isize), itail:
C
C                       LU factors are located in II (itail ... isize),
C                       including pattern, row and column permutations,
C                       block triangular information, etc.  See umf2fa
C                       for more information.
C
C       Info:           integer informational output, see MA38AD
C       Rinfo:          real informational output, see MA38AD
C
C       iuse:           memory usage in Index on output
C       xuse:           memory usage in Value on output
C
C       nzoff:          entries in off-diagonal part (0 if BTF not used)
C       nblks:          number of diagonal blocks (1 if BTF not used)

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38AD
C       subroutines called:     MA38ND, MA38HD, MA38ED, MA38MD
C       functions called:       MAX
        INTRINSIC MAX
        EXTERNAL MA38ED, MA38HD, MA38MD, MA38ND

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER KN, NZDIA, BLKPP, LUBLPP, P, OFFIP, XHEAD, ROW,
     $          OFFXP, OFFPP, IHEAD, K1, K2, BLK, PRP, P2, CPERMP,
     $          RPERMP, NSGLTN, NPIV, MNZ, NSYM, K, COL, RMAX, CMAX,
     $          TOTNLU, XRMAX, XRUSE
        LOGICAL TRYBTF, IOUT, XOUT
        DOUBLE PRECISION
     $          ZERO, ONE, A
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C  Allocated array pointers:
C  -------------------------
C  blkpp:   Blkp (1..nblks+1) array located in II (blkpp..blkp+nblks)
C  lublpp:  LUblkp (1..nblks) array loc. in II (lublpp..lublpp+nblks-1)
C  offip:   Offi (1..nzoff) array located in II (offip..offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array located in XX (offxp..offxp+nzoff-1)
C  offpp:   Offp (1..n+1) array located in II (offpp..offpp+n)
C  cpermp:  Cperm (1..n) array located in II (cpermp..cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in II (rpermp..rpermp+n-1)
C  prp:     Pr (1..n) work array located in II (prp..prp+n-1)
C
C  BTF information:
C  ----------------
C  k1:      starting index of diagonal block being factorized
C  k2:      ending index of diagonal block being factorized
C  kn:      the order of the diagonal block being factorized
C  blk:     block number of diagonal block being factorized
C  trybtf:  true if BTF is to be attempted (= Icntl (4) .eq. 1)
C  nzdia:   number of entries in diagonal blocks (= nz if BTF not used)
C  nsgltn:  number of 1-by-1 diagonal blocks ("singletons")
C  npiv:    number of numerically valid singletons
C  a:       numerical value of a singleton
C  mnz:     nzoff
C
C  Memory usage:
C  -------------
C  xhead:   XX (1..xhead-1) is in use, XX (xhead..xtail-1) is free
C  ihead:   II (1..ihead-1) is in use, II (ihead..itail-1) is free
C  iout:    true if MA38ED ran out of integer memory, but did not
C           set error flag
C  xout:    true if MA38FD ran out of real memory, but did not
C           set error flag
C
C  Estimated memory for MA38BD:
C  ----------------------------
C  rmax:    largest contribution block is cmax-by-rmax
C  cmax:       "         "         "    "   "   "  "
C  totnlu:  total number of LU arrowheads in all diagonal blocks
C  xrmax:   estimated maximum real memory usage for MA38BD
C  xruse:   estimated current real memory usage for MA38BD
C
C  Other:
C  ------
C  k:       loop index (kth pivot)
C  row:     row index
C  col:     column index
C  p:       pointer
C  p2:      pointer
C  nsym:    number of symmetric pivots chosen

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  get input parameters and initialize
C-----------------------------------------------------------------------

        NBLKS = 1
        NZOFF = 0
        NZDIA = NZ
        NSGLTN = 0
        NPIV = 0
        RMAX = 1
        CMAX = 1
        TOTNLU = 0
        XOUT = .FALSE.
        IOUT = .FALSE.
        XRMAX = 2*NE
        IF (PRESRV) THEN
C          original matrix is not in Cp/II/XX, but in Ap/Ai/Ax:
           IHEAD = 1
           XHEAD = 1
        ELSE
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1

C-----------------------------------------------------------------------
C  allocate permutation arrays: Cperm (1..n) and Rperm (1..n), and
C  seven scalars:  transa, nzoff, nblks, presrv, nz, n, ne
C  (in that order) at tail of II (in LU factors)
C-----------------------------------------------------------------------

        ITAIL = ITAIL - (2*N+7)
        IUSE = IUSE + (2*N+7)
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = INFO (18)
        CPERMP = ITAIL
        RPERMP = CPERMP + N
        IF (IHEAD .GT. ITAIL) THEN
C          error return, if not enough integer memory:
           GO TO 9000
        ENDIF

C-----------------------------------------------------------------------
C  Find permutations to block upper triangular form, if requested.
C-----------------------------------------------------------------------

        TRYBTF = ICNTL (4) .EQ. 1
        IF (TRYBTF) THEN

C          -------------------------------------------------------------
C          get workspace at tail of II of size 6n+2
C          -------------------------------------------------------------

           ITAIL = ITAIL - (N+1)
           OFFPP = ITAIL
           ITAIL = ITAIL - (5*N+1)
           P = ITAIL
           IUSE = IUSE + (6*N+2)
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = INFO (18)

C          -------------------------------------------------------------
           IF (PRESRV) THEN
C          find permutation, but do not convert to BTF form
C          -------------------------------------------------------------

              IF (IHEAD .GT. ITAIL) THEN
C                error return, if not enough integer memory:
                 GO TO 9000
              ENDIF
              CALL MA38HD (AX, ANZ, AI, ANZ, N, NZ, NZDIA, NZOFF,
     $           NBLKS, CP, II (CPERMP), II (RPERMP), II(P), II(P+N),
     $           II (P+2*N), II (P+3*N), II (P+4*N), II (OFFPP),
     $           PRESRV)

C          -------------------------------------------------------------
           ELSE
C          find permutation, convert to BTF form, and discard original
C          -------------------------------------------------------------

C             use additional size nz temporary workspace in II and XX
              IHEAD = IHEAD + NZ
              XHEAD = XHEAD + NZ
              IUSE = IUSE + NZ
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (19) = INFO (18)
              INFO (21) = INFO (20)
              IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
C                error return, if not enough integer and/or real memory:
                 GO TO 9000
              ENDIF
              CALL MA38HD (XX, 2*NZ, II, 2*NZ, N, NZ, NZDIA, NZOFF,
     $              NBLKS, CP, II (CPERMP), II (RPERMP), II(P), II(P+N),
     $              II (P+2*N), II (P+3*N), II (P+4*N), II (OFFPP),
     $              PRESRV)
C             deallocate extra workspace in II and XX
              IHEAD = IHEAD - NZ
              XHEAD = XHEAD - NZ
              IUSE = IUSE - NZ
              XUSE = XUSE - NZ
           ENDIF

C          -------------------------------------------------------------
C          deallocate workspace, and allocate BTF arrays if required
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1) THEN
C             replace (6*n+2) workspace at tail of II with
C             Blkp (1..nblks+1) and LUblkp (1..nblks), Offp (1..n+1)
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              ITAIL = LUBLPP
              IUSE = IUSE - (6*N+2) + (2*NBLKS+N+2)
           ELSE
C             The matrix is irreducible.  There is only one block.
C             Remove everything at tail of II, except
C             for the 2*n permutation vectors and the 7 scalars.
C             (transa, nzoff, nblks, presrv, nz, n, ne).
              ITAIL = (ISIZE + 1) - (2*N+7)
              IUSE = IUSE - (6*N+2)
           ENDIF

        ENDIF

C-----------------------------------------------------------------------
C current memory usage:
C-----------------------------------------------------------------------

C       if .not. presrv then
C               input matrix is now in II (1..nz) and XX (1..nz)
C               off-diagonal part: in II/XX (1..nzoff)
C                       col pattern: II (Offp (col) ... Offp (col+1))
C                       col values:  XX (Offp (col) ... Offp (col+1))
C               diagonal blocks: in II/XX (nzoff+1..nz)
C                       col pattern: II (Cp (col) ... Cp (col+1))
C                       col values:  XX (Cp (col) ... Cp (col+1))
C               total: nz+n+1 integers, nz reals
C       else
C               input matrix is now in Ai (1..nz) and Ax (1..nz),
C               in original (non-BTF) order:
C                       col pattern: Ai (Ap (col) ... Ap (col+1))
C                       col values:  Ax (Ap (col) ... Ap (col+1))
C               Cp is a size n+1 integer workspace
C               total: nz+2*(n+1) integers, nz reals
C
C       if (nblks > 1) then
C               at tail of II (in order): 2*nblks+n+2
C                       LUblkp (1..nblks)
C                       Blkp (1..nblks+1)
C                       Offp (1..n+1)
C               total: (2*nblks+n+2) integers
C
C       remainder at tail of II:
C               Cperm (1..n)
C               Rperm (1..n)
C               seven scalars: transa, nzoff, nblks, presrv, nz, n, ne
C
C   Grand total current memory usage (including II,XX,Cp,Ai,Ap,Ax):
C
C       presrv  nblks>1         integers, iuse =
C       F       F               nz+  (n+1)+(2*n+7)
C       F       T               nz+  (n+1)+(2*n+7)+(2*nblks+n+2)
C       T       F               nz+2*(n+1)+(2*n+7)
C       T       T               nz+2*(n+1)+(2*n+7)+(2*nblks+n+2)
C
C   real usage is xuse = nz

C       ----------------------------------------------------------------
C       get memory usage for next call to MA38BD
C       ----------------------------------------------------------------

        XRUSE = NZ

C-----------------------------------------------------------------------
C factorization
C-----------------------------------------------------------------------

        IF (NBLKS .EQ. 1) THEN

C          -------------------------------------------------------------
C          factorize the matrix as a single block
C          -------------------------------------------------------------

           CALL MA38ED (CP, N, II (CPERMP), II (RPERMP), NZOFF,
     $          ITAIL, XTAIL, XX, XSIZE, XUSE, II, ITAIL-1, IUSE,
     $          ICNTL, CNTL, INFO, RINFO, NBLKS,
     $          AP, AI, AX, PRESRV, 1, AN, ANZ, II, KEEP,
     $          RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
           IF (IOUT .OR. XOUT) THEN
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF
           IF (INFO (1) .LT. 0) THEN
C             error return, if error in MA38FD:
              GO TO 9010
           ENDIF
C          original matrix has been deallocated
           IHEAD = 1
           XHEAD = 1

C          -------------------------------------------------------------
C          make the index of the block relative to start of LU factors
C          -------------------------------------------------------------

           II (ITAIL) = 1

        ELSE

C          -------------------------------------------------------------
C          factorize the block-upper-triangular form of the matrix
C          -------------------------------------------------------------

           PRP = OFFPP
           IF (PRESRV) THEN
C             count the off-diagonal entries during factorization
              NZOFF = 0
C             compute temp inverse permutation in II (prp..prp+n-1)
CFPP$ NODEPCHK L
              DO 10 K = 1, N
                 II (PRP + II (RPERMP+K-1) - 1) = K
10            CONTINUE
           ENDIF

           DO 30 BLK = NBLKS, 1, -1

C             ----------------------------------------------------------
C             factorize the kn-by-kn block, A (k1..k2, k1..k2)
C             ----------------------------------------------------------

C             get k1 and k2, the start and end of this block
              K1 = II (BLKPP+BLK-1)
              K2 = II (BLKPP+BLK) - 1
              KN = K2-K1+1
              IF (.NOT. PRESRV) THEN
                 P = CP (K1)
                 CP (K2+1) = IHEAD
              ENDIF

              IF (KN .GT. 1) THEN

C                -------------------------------------------------------
C                factor the block (the block is not a singleton)
C                -------------------------------------------------------

                 CALL MA38ED (CP (K1), KN,
     $              II (CPERMP+K1-1), II (RPERMP+K1-1), NZOFF,
     $              ITAIL, XTAIL, XX, XTAIL-1, XUSE, II, ITAIL-1,
     $              IUSE, ICNTL, CNTL, INFO, RINFO, NBLKS,
     $              AP, AI, AX, PRESRV, K1, AN, ANZ, II (PRP), KEEP,
     $              RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
                 IF (IOUT .OR. XOUT) THEN
C                   error return, if not enough int. and/or real memory:
                    GO TO 9000
                 ENDIF
                 IF (INFO (1) .LT. 0) THEN
C                   error return, if error in MA38FD:
                    GO TO 9010
                 ENDIF
                 IF (PRESRV) THEN
                    IHEAD = 1
                    XHEAD = 1
                 ELSE
                    IHEAD = P
                    XHEAD = P
                 ENDIF

C                -------------------------------------------------------
C                save the location of the LU factors in LUbkp (blk)
C                -------------------------------------------------------

                 II (LUBLPP+BLK-1) = ITAIL

              ELSE

C                -------------------------------------------------------
C                get the value of singleton at A (k1,k1), if it exists
C                -------------------------------------------------------

                 A = ZERO
                 IF (PRESRV) THEN
C                   find the diagonal entry in the unpermuted matrix
                    COL = II (CPERMP + K1 - 1)
                    DO 20 P2 = AP (COL), AP (COL + 1) - 1
                       ROW = II (PRP + AI (P2) - 1)
                       IF (ROW .LT. K1) THEN
C                         entry in off-diagonal blocks
                          NZOFF = NZOFF + 1
                       ELSE
                          A = AX (P2)
                       ENDIF
20                  CONTINUE
                    IHEAD = 1
                    XHEAD = 1
                 ELSE IF (P .NE. IHEAD) THEN
                    A = XX (P)
                    IHEAD = P
                    XHEAD = P
                    IUSE = IUSE - 1
                    XUSE = XUSE - 1
                    XRUSE = XRUSE - 1
                 ENDIF

C                -------------------------------------------------------
C                store the 1-by-1 LU factors of a singleton
C                -------------------------------------------------------

                 NSGLTN = NSGLTN + 1
                 IF (A .EQ. ZERO) THEN
C                   the diagonal entry is either not present, or present
C                   but numerically zero.  This is a singular matrix,
C                   replace with 1-by-1 identity matrix.
                    A = ONE
                 ELSE
C                   increment pivot count
                    NPIV = NPIV + 1
                 ENDIF
                 XTAIL = XTAIL - 1
C                note: if the matrix is not preserved and nonsingular
C                then we will not run out of memory at this point.
                 XUSE = XUSE + 1
                 XRUSE = XRUSE + 1
                 XRMAX = MAX (XRMAX, XRUSE)
                 INFO (20) = MAX (INFO (20), XUSE)
                 INFO (21) = MAX (INFO (21), XUSE)
C                error return, if not enough real memory:
                 IF (XHEAD .GT. XTAIL) THEN
                    GO TO 9000
                 ENDIF
                 II (LUBLPP+BLK-1) = -XTAIL
                 XX (XTAIL) = A

              ENDIF

30         CONTINUE

C          -------------------------------------------------------------
C          make the index of each block relative to start of LU factors
C          -------------------------------------------------------------

CFPP$ NODEPCHK L
           DO 40 P = LUBLPP, LUBLPP + NBLKS - 1
              IF (II (P) .GT. 0) THEN
                 II (II (P)) = II (II (P)) - XTAIL + 1
                 II (P) = II (P) - ITAIL + 1
              ELSE
C                this is a singleton
                 II (P) = (-II (P)) - XTAIL + 1
              ENDIF
40         CONTINUE

C          -------------------------------------------------------------
C          allocate temporary workspace for Pr (1..n) at head of II
C          -------------------------------------------------------------

           PRP = IHEAD
           IHEAD = IHEAD + N
           IUSE = IUSE + N

C          -------------------------------------------------------------
C          allocate a single entry in case the LU factors are empty
C          -------------------------------------------------------------

           IF (NBLKS .EQ. N) THEN
C             otherwise, arrays in MA38BD and MA38CD would have
C             zero size, which can cause an address fault later on
              ITAIL = ITAIL - 1
              IUSE = IUSE + 1
              P2 = ITAIL
           ENDIF

C          -------------------------------------------------------------
C          allocate permanent copy of off-diagonal blocks
C          -------------------------------------------------------------

           ITAIL = ITAIL - NZOFF
           OFFIP = ITAIL
           XTAIL = XTAIL - NZOFF
           OFFXP = XTAIL
           IUSE = IUSE + NZOFF
           XUSE = XUSE + NZOFF
           XRUSE = XRUSE + NZOFF
           XRMAX = MAX (XRMAX, XRUSE)
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), IUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XUSE)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF

C          -------------------------------------------------------------
C          re-order the off-diagonal blocks according to pivot perm
C          -------------------------------------------------------------

C          use Cp as temporary work array:
           MNZ = NZOFF
           IF (PRESRV) THEN
              CALL MA38MD (CP, N, II (RPERMP), II (CPERMP), NZOFF,
     $          II (OFFPP), II (OFFIP), XX (OFFXP), II (PRP),
     $          ICNTL, AP, AI, AX, AN, ANZ, PRESRV, NBLKS, II (BLKPP),
     $          MNZ, 1, P)
           ELSE
              CALL MA38MD (CP, N, II (RPERMP), II (CPERMP), NZOFF,
     $          II (OFFPP), II (OFFIP), XX (OFFXP), II (PRP),
     $          ICNTL, AP, II, XX, 0, MNZ, PRESRV, 0, II(BLKPP),
     $          MNZ, 1, P)
           ENDIF
           IF (NBLKS .EQ. N) THEN
C             zero the only entry in the integer part of the LU factors
              II (P2) = 0
           ENDIF

C          -------------------------------------------------------------
C          deallocate Pr (1..n), and II/XX (1..nzoff) if present
C          -------------------------------------------------------------

           IHEAD = 1
           XHEAD = 1
           IUSE = IUSE - N
           IF (.NOT. PRESRV) THEN
              IUSE = IUSE - NZOFF
              XUSE = XUSE - NZOFF
           ENDIF

        ENDIF

C-----------------------------------------------------------------------
C  normal and error return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (IOUT .OR. IHEAD .GT. ITAIL) THEN
C          set error flag if not enough integer memory
           CALL MA38ND (1, ICNTL, INFO, -3, INFO (19))
        ENDIF
        IF (XOUT .OR. XHEAD .GT. XTAIL) THEN
C          set error flag if not enough real memory
           CALL MA38ND (1, ICNTL, INFO, -4, INFO (21))
        ENDIF

C       error return label, for error from MA38FD:
9010    CONTINUE

        INFO (4) = 0
        NZDIA = NZ - NZOFF
        INFO (5) = NZ
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (8) = NSGLTN
        INFO (9) = NBLKS
        INFO (12) = INFO (10) + INFO (11) + N + INFO (7)

C       Count the number of symmetric pivots chosen.  Note that some of
C       these may have been numerically unacceptable.
        NSYM = 0
        IF (INFO (1) .GE. 0) THEN
           DO 50 K = 1, N
              IF (II (CPERMP+K-1) .EQ. II (RPERMP+K-1)) THEN
C                this kth pivot came from the diagonal of A
                 NSYM = NSYM + 1
              ENDIF
50         CONTINUE
        ENDIF
        INFO (16) = NSYM

        INFO (17) = INFO (17) + NPIV
        RINFO (1) = RINFO (4) + RINFO (5) + RINFO (6)

        IF (INFO (1) .GE. 0 .AND. INFO (17) .LT. N) THEN
C          set warning flag if matrix is singular
           CALL MA38ND (1, ICNTL, INFO, 4, INFO (17))
        ENDIF

C       ----------------------------------------------------------------
C       Determine an upper bound on the amount of integer memory needed
C       (LINDEX) for a subsequent call to MA38BD.  If block-upper-
C       triangular-form is not in use (Info (9) is 1), then
C       this bound is exact.  If NE is higher in the call to MA38BD
C       than in the call to MA38AD, then add 3 integers for each
C       additional entry (including the 2 integers required for the
C       row and column indices of the additional triplet itself).
C       This estimate assumes that JOB and TRANSA are the same in
C       MA38AD and MA38BD.
C       ----------------------------------------------------------------

C       (Keep (5) - Keep (4) + 1), is added to Info (22)
C       in MA38AD, to complete the computation of the estimate.

        IF (PRESRV) THEN
           INFO (22) = MAX (3*NE+2*N+1, NE+3*N+2,
     $                           2*NZ+4*N+10+RMAX+3*CMAX+4*TOTNLU)
        ELSE
           INFO (22) = MAX (3*NE+2*N+1, NE+3*N+2, 2*NZ+3*N+2,
     $                             NZ+3*N+ 9+RMAX+3*CMAX+4*TOTNLU)
        ENDIF

C       ----------------------------------------------------------------
C       Approximate the amount of real memory needed (LVALUE) for a
C       subsequent call to MA38BD.  The approximation is an upper bound
C       on the bare minimum amount needed.  Some garbage collection may
C       occur, but MA38BD is guaranteed to finish if given an LVALUE of
C       size Info (23) and if the pattern is the same.  If NE is
C       higher in the call to MA38BD than in the call to MA38AD, then
C       add 2 reals for each additional entry (including the 1 real
C       required for the value of the additional triplet itself).
C       This estimate assumes that JOB and TRANSA are the same in
C       MA38AD and MA38BD.
C       ----------------------------------------------------------------

        INFO (23) = XRMAX
        RETURN
        END

        SUBROUTINE MA38ED (CP, N, CPERM, RPERM, NZOFF,
     *          ITAIL, XTAIL, XX, XSIZE, XUSE, II, ISIZE, IUSE,
     *          ICNTL, CNTL, INFO, RINFO, NBLKS,
     *          AP, AI, AX, PRESRV, K1, AN, ANZ, PR, KEEP,
     *          RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
        INTEGER XSIZE, ISIZE, N, ICNTL(20), INFO(40), XUSE, IUSE,
     *          ITAIL, XTAIL, II(ISIZE), CP(N+1), CPERM(N), NZOFF,
     *          AN, ANZ, RPERM(N), AI(ANZ), AP(AN+1), K1, PR(AN),
     *          NBLKS, KEEP(20), RMAX, CMAX, TOTNLU, XRMAX, XRUSE
        LOGICAL PRESRV, IOUT, XOUT
        DOUBLE PRECISION XX(XSIZE), CNTL(10), RINFO(20), AX(ANZ)

C=== MA38ED ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  MA38ED factorizes the n-by-n column-form matrix at the head of II/XX
C  or in Ap/Ai/Ax, and places its LU factors at the tail of II/XX.  The
C  input matrix overwritten if it is located in II/XX on input.  If
C  block-triangular-form (BTF) is in use, this routine factorizes a
C  single diagonal block.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix (or order of diagonal block
C                       if BTF is in use).
C       Cp (1..n+1):    column pointers for input matrix
C       nblks:          number of diagonal blocks in BTF form
C       isize:          size of II
C       xsize:          size of XX
C       k1:             first index of this matrix (1 if BTF not used)
C       Icntl:          integer control parameters, see MA38ID
C       Cntl:           real control parameters, see MA38ID
C       Keep (6..8):    integer control parameters, see MA38ID
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       rmax:           maximum ludegr seen so far (see MA38FD for info)
C       cmax:           maximum ludegc seen so far (see MA38FD for info)
C       totnlu:         total number of LU arrowheads constructed so far
C
C       if nblks>1 then:
C          Cperm (1..n):        col permutation to BTF
C          Rperm (1..n):        row permutation to BTF
C       else
C          Cperm (1..n):        undefined on input
C          Rperm (1..n):        undefined on input
C
C
C       presrv:         if true then input matrix is preserved
C
C       if presrv is true then:
C           an:                 order of preserved matrix (all blocks)
C           anz:                entries in preserved matrix
C           Ap (1..an+1):       column pointers for preserved matrix
C           Ai (1..anz):        row indices of preserved matrix
C           Ax (1..anz):        values of preserved matrix
C                               The preserved matrix is not in BTF form;
C                               it is in the orginal order.
C           if nblks > 1:
C               Pr (1..n):      inverse row permutations to BTF form
C               nzoff           entries in off-diagonal blocks
C                               seen so far
C           else
C               Pr (1..n):      undefined on input
C
C           II (1..isize):      undefined on input
C           XX (1..xsize):      undefined on input
C           Cp (1..n+1):        undefined on input
C
C       else, if presrv is false:
C           an:                         1
C           anz:                        1
C           II (1..Cp (1) - 1):         unused
C           II (Cp (1) ... Cp (n+1)-1): row indices of matrix to factor,
C                                       will be overwritten on output
C           II (Cp (n+1) ... isize):    unused on input
C
C           XX (1..Cp (1) - 1):         unused
C           XX (Cp (1) ... Cp (n+1)-1): values of matrix to factorize,
C                                       will be overwritten on output
C           XX (Cp (n+1) ... xsize):    unused on input
C                       If BTF is in use, then II and XX contain a
C                       single diagonal block.

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       XX (xtail ... xsize), xtail,  II (itail ... isize), itail:
C
C                       The LU factors of a single diagonal block.
C                       See MA38FD for a description.
C
C       II (cp1 ... itail-1):   undefined on output
C       XX (cp1 ... xtail-1):   undefined on output,
C                       where cp1 is equal to the value of Cp (1)
C                       if presrv is false, or cp1 = 1 if presrv is
C                       true.
C
C       Info:           integer informational output, see MA38AD
C       Rinfo:          real informational output, see MA38AD
C       Cperm (1..n):   the final col permutations, including BTF
C       Rperm (1..n):   the final row permutations, including BTF
C
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       rmax:           maximum ludegr seen so far (see MA38FD for info)
C       cmax:           maximum ludegc seen so far (see MA38FD for info)
C       totnlu:         total number of LU arrowheads constructed so far
C
C       if nblks>1 and presrv:
C           nzoff       entries in off-diagonal blocks seen so far
C
C       iout:           true if ran out of integer memory in MA38ED,
C         and if the corresponding error status has not
C         yet been set (the error status is set in the caller,
C         MA38DD)
C       xout:           true if ran out of real memory in MA38ED,
C         and if the corresponding error status has not
C         yet been set (the error status is set in the caller,
C         MA38DD)

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38DD
C       subroutines called:     MA38FD
C       functions called:       MAX, SQRT
        INTRINSIC MAX, SQRT
        EXTERNAL MA38FD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER CP1, PC, PEND, PCOL, CDEG, COL, CSIZ, NZ, XP, IP, IS, P,
     $          DN, DSIZ, WRKSIZ, I, CLEN, D1, D2, N2, ROW, CSCAL
        PARAMETER (CSCAL = 9)
        DOUBLE PRECISION XN

C  Original and expanded column-form:
C  ----------------------------------
C  cp1:     = value of Cp (1) on input
C  pc:      pointer to integer part of expanded column-form matrix
C  pend:    column col ends here in the input column-form matrix
C  pcol:    column col starts here in the input column-form matrix
C  cdeg:    degree (number of entries) in a column
C  clen:    number of original entries in a column (= degree, here,
C           but decreases in MA38FD)
C  csiz:    size of the integer part of an expanded column (cdeg+cscal)
C  cscal:   = 9, the number of scalars in column data structure
C  nz:      number of entries in the diagonal block being factorized
C  xp:      pointer to real part of the expanded column
C  ip:      pointer to integer part of the expanded column
C
C  Memory usage:
C  -------------
C  wrksiz:  size of integer workspace needed by MA38FD
C
C  "Dense" columns: (converted to a prior, or artificial, frontal mat.):
C  ----------------
C  d1:      = Keep (7), dense column control
C  d2:      = Keep (8), dense column control
C  dn:      number of "dense" columns
C  dsiz:    a column is "dense" if it has more than dsiz entries
C  xn:      = sqrt (real (n))
C  n2:      = int (sqrt (real (n)))
C
C  Other:
C  ------
C  row:     row index
C  col:     a column index
C  i:       loop index
C  p:       pointer

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IOUT = .FALSE.
        XOUT = .FALSE.

C-----------------------------------------------------------------------
C  count "dense" columns (they are treated as a priori frontal matrices)
C-----------------------------------------------------------------------

C       a column is "dense" if it has more than dsiz entries
        D1 = KEEP (7)
        D2 = KEEP (8)
        XN = N
        XN = SQRT (XN)
        N2 = XN
        DSIZ = MAX (0, D1, D2 * N2)
        DN = 0
        IF (PRESRV) THEN
           IF (NBLKS .EQ. 1) THEN
              DO 10 COL = 1, N
                 IF (AP (COL+1) - AP (COL) .GT. DSIZ) THEN
C                   this is a "dense" column
                    DN = DN + 1
                 ENDIF
10            CONTINUE
           ELSE
              DO 40 COL = 1, N
C                if col might be dense, check more carefully:
                 CDEG = AP (CPERM (COL) + 1)- AP (CPERM (COL))
                 IF (CDEG .GT. DSIZ) THEN
                    CDEG = 0
                    DO 20 P = AP (CPERM (COL)), AP (CPERM (COL) + 1) -1
                       ROW = PR (AI (P))
                       IF (ROW .GE. K1) THEN
                          CDEG = CDEG + 1
                          IF (CDEG .GT. DSIZ) THEN
C                            this is a "dense" column, exit out of loop
                             DN = DN + 1
                             GO TO 30
                          ENDIF
                       ENDIF
20                  CONTINUE
C                   loop exit label:
30                  CONTINUE
                 ENDIF
40            CONTINUE
           ENDIF
        ELSE
           DO 50 COL = 1, N
              IF (CP (COL+1) - CP (COL) .GT. DSIZ) THEN
C                this is a "dense" column
                 DN = DN + 1
              ENDIF
50         CONTINUE
        ENDIF

C-----------------------------------------------------------------------
C  get size of workspaces to allocate from II
C-----------------------------------------------------------------------

C       workspaces: WiR (n), WiC (n), WpR (n), WpC (n),
C       Wm (n), Head (n), Rp (n+dn), Wc (n+dn), Wr (n+dn), Wj (n)
        IF (NBLKS .EQ. 1) THEN
C          Rperm (1..n) is used as WiR (1..n), and
C          Cperm (1..n) is used as WiC (1..n) in MA38FD
           WRKSIZ = 8*N + 3*DN
        ELSE
           WRKSIZ = 10*N + 3*DN
        ENDIF

C-----------------------------------------------------------------------
C  construct the expanded column-form of the matrix or the diag. block
C-----------------------------------------------------------------------

        IF (PRESRV) THEN

C          -------------------------------------------------------------
C          allocate space for wrksiz workspace and nz+cscal*n
C          integers and nz reals for the expanded column-form matrix.
C          -------------------------------------------------------------

           CP1 = 1
           XP = 1
           IP = 1 + WRKSIZ
           IF (NBLKS .EQ. 1) THEN

C             ----------------------------------------------------------
C             construct copy of entire matrix
C             ----------------------------------------------------------

              NZ = ANZ
              IS = NZ + WRKSIZ + CSCAL*N
              IUSE = IUSE + IS
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (19) = MAX (INFO (19), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
              IOUT = IS .GT. ISIZE
              XOUT = NZ .GT. XSIZE
              IF (IOUT .OR. XOUT) THEN
C                error return, if not enough integer and/or real memory:
                 GO TO 9000
              ENDIF

              PC = IP
              DO 70 COL = 1, N
                 CP (COL) = PC - WRKSIZ
                 CDEG = AP (COL+1) - AP (COL)
                 CLEN = CDEG
                 CSIZ = CDEG + CSCAL
                 II (PC) = CSIZ
                 II (PC+1) = CDEG
                 II (PC+5) = 0
                 II (PC+6) = CLEN
                 II (PC+7) = 0
                 II (PC+8) = 0
                 II (PC+2) = XP
                 XP = XP + CDEG
                 PC = PC + CSCAL
                 P = AP (COL)
                 DO 60 I = 0, CDEG - 1
                    II (PC + I) = AI (P + I)
60               CONTINUE
                 PC = PC + CDEG
70            CONTINUE
              DO 80 P = 1, NZ
                 XX (P) = AX (P)
80            CONTINUE

           ELSE

C             ----------------------------------------------------------
C             construct copy of a single block in BTF form
C             ----------------------------------------------------------

C             check for memory usage during construction of block
              DO 100 COL = 1, N
                 PC = IP
                 CP (COL) = PC - WRKSIZ
                 IP = IP + CSCAL
                 IOUT = IP .GT. ISIZE
                 IF (IOUT) THEN
C                   error return, if not enough integer memory:
                    GO TO 9000
                 ENDIF
                 II (PC+2) = XP
                 CDEG = IP
                 DO 90 P = AP (CPERM (COL)), AP (CPERM (COL)+1)-1
                    ROW = PR (AI (P))
                    IF (ROW .GE. K1) THEN
                       IOUT = IP .GT. ISIZE
                       XOUT = XP .GT. XSIZE
                       IF (IOUT .OR. XOUT) THEN
C                         error return, if not enough memory
                          GO TO 9000
                       ENDIF
                       II (IP) = ROW - K1 + 1
                       XX (XP) = AX (P)
                       IP = IP + 1
                       XP = XP + 1
                    ELSE
C                      entry in off-diagonal part
                       NZOFF = NZOFF + 1
                    ENDIF
90               CONTINUE
                 CDEG = IP - CDEG
                 CLEN = CDEG
                 CSIZ = CDEG + CSCAL
                 II (PC) = CSIZ
                 II (PC+1) = CDEG
                 II (PC+5) = 0
                 II (PC+6) = CLEN
                 II (PC+7) = 0
                 II (PC+8) = 0
100           CONTINUE

              NZ = XP - 1
              IS = NZ + WRKSIZ + CSCAL*N
              IUSE = IUSE + IS
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (19) = MAX (INFO (19), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)

           ENDIF

C          -------------------------------------------------------------
C          get memory usage for next call to MA38BD
C          -------------------------------------------------------------

           XRUSE = XRUSE + NZ
           XRMAX = MAX (XRMAX, XRUSE)

        ELSE

C          -------------------------------------------------------------
C          allocate space for wrksiz workspace and additional cscal*n
C          space for the expanded column-form of the matrix.
C          -------------------------------------------------------------

           CP1 = CP (1)
           NZ = CP (N+1) - CP1
           PC = CP1 + WRKSIZ + (NZ+CSCAL*N)
           IUSE = IUSE + WRKSIZ + CSCAL*N
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), IUSE)
           IOUT = PC .GT. ISIZE+1
           IF (IOUT) THEN
C             error return, if not enough integer memory:
              GO TO 9000
           ENDIF

C          -------------------------------------------------------------
C          expand the column form in place and make space for workspace
C          -------------------------------------------------------------

           XP = NZ + 1
           IP = NZ + CSCAL*N + 1
           PEND = CP (N+1)
           DO 120 COL = N, 1, -1
              PCOL = CP (COL)
              DO 110 P = PEND-1, PCOL, -1
                 PC = PC - 1
                 II (PC) = II (P)
110           CONTINUE
              PC = PC - CSCAL
              CDEG = PEND - PCOL
              CLEN = CDEG
              PEND = PCOL
              CSIZ = CDEG + CSCAL
              IP = IP - CSIZ
              CP (COL) = IP
              II (PC) = CSIZ
              II (PC+1) = CDEG
              II (PC+5) = 0
              II (PC+6) = CLEN
              II (PC+7) = 0
              II (PC+8) = 0
              XP = XP - CDEG
              II (PC+2) = XP
120        CONTINUE
        ENDIF

C-----------------------------------------------------------------------
C  factorize the expanded column-form, with allocated workspaces
C-----------------------------------------------------------------------

        XP = CP1
        IP = CP1 + WRKSIZ

        IF (NBLKS .EQ. 1) THEN

C          pass Rperm and Cperm as the WiR and WiC arrays in MA38FD:
           CALL MA38FD (CP, NZ, N, 1, CPERM, RPERM, ITAIL, XTAIL,
     $          XX (XP), XSIZE-XP+1, II (IP), ISIZE-IP+1, ICNTL, CNTL,
     $          INFO, RINFO, .FALSE., IUSE, XUSE,
     $          RPERM, CPERM, II (CP1), II (CP1+N),
     $          II (CP1+2*N), II (CP1+3*N), II (CP1+4*N), II (CP1+5*N),
     $          II (CP1+6*N+DN), II (CP1+7*N+2*DN),
     $          DN, DSIZ, RMAX, CMAX, TOTNLU, XRMAX, XRUSE)

        ELSE

C          pass Cperm, Rperm, WiC and WiR as separate arrays, and
C          change Cperm and Rperm from the BTF permutations to the
C          final permutations (including BTF and numerical pivoting).
           CALL MA38FD (CP, NZ, N, N, CPERM, RPERM, ITAIL, XTAIL,
     $          XX (XP), XSIZE-XP+1, II (IP), ISIZE-IP+1, ICNTL, CNTL,
     $          INFO, RINFO, .TRUE., IUSE, XUSE,
     $          II (CP1), II (CP1+N), II (CP1+2*N), II (CP1+3*N),
     $          II (CP1+4*N), II (CP1+5*N), II (CP1+6*N), II (CP1+7*N),
     $          II (CP1+8*N+DN), II (CP1+9*N+2*DN),
     $          DN, DSIZ, RMAX, CMAX, TOTNLU, XRMAX, XRUSE)

        ENDIF

        IF (INFO (1) .LT. 0) THEN
C          error return, if error occured in MA38FD:
           RETURN
        ENDIF

C-----------------------------------------------------------------------
C  adjust tail pointers, and save pointer to numerical part of LU
C-----------------------------------------------------------------------

C       head = cp1
        IUSE = IUSE - WRKSIZ
        ITAIL = ITAIL + IP - 1
        XTAIL = XTAIL + XP - 1
        II (ITAIL) = XTAIL
        RETURN

C=======================================================================
C  Error return
C=======================================================================

C       error return label:
9000    CONTINUE
        RETURN
        END

        SUBROUTINE MA38FD (CP, NZ, N, PN, CPERM, RPERM, ITAIL, XTAIL,
     *          XX, XSIZE, II, ISIZE, ICNTL, CNTL, INFO, RINFO, PGIVEN,
     *          IUSE, XUSE, WIR, WIC, WPR, WPC, WM, HEAD,
     *          WJ, RP, WC, WR, DN, DSIZ,
     *          RMAX, CMAX, TOTNLU, XRMAX, XRUSE)
        INTEGER XSIZE, ISIZE, ICNTL(20), INFO(40), PN,
     *          ITAIL, XTAIL, NZ, N, II(ISIZE), CP(N+1), DN, DSIZ,
     *          RPERM(PN), CPERM(PN), WIR(N), WIC(N), WPR(N),
     *          WPC(N), WM(N), HEAD(N), RP(N+DN), WC(N+DN),
     *          WR(N+DN), IUSE, XUSE, WJ(N),
     *          RMAX, CMAX, TOTNLU, XRMAX, XRUSE
        LOGICAL PGIVEN
        DOUBLE PRECISION XX(XSIZE), CNTL(10), RINFO(20)

C=== MA38FD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  MA38FD factorizes the n-by-n input matrix at the head of II/XX
C  (in expanded column-form) and places its LU factors at the tail of
C  II/XX.  The input matrix is overwritten.   No BTF information is
C  used in this routine, except that the BTF permutation arrays are
C  modified to include the final permutations.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       Cp (1..n+1):    column pointers of expanded column-form,
C                       undefined on output
C       n:              order of input matrix
C       nz:             entries in input matrix
C       isize:          size of II
C       xsize:          size of XX
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       Icntl:          integer control parameters, see MA38ID
C       Cntl:           real control parameters, see MA38ID
C       Keep (6)        integer control parameter, see MA38ID
C       dn:             number of dense columns
C       dsiz:           entries required for col to be treated as dense
C       rmax:           maximum ludegr seen so far (see below)
C       cmax:           maximum ludegc seen so far (see below)
C       totnlu:         total number of LU arrowheads constructed so far
C       xrmax:          maximum real memory usage for MA38BD
C       xruse:          current real memory usage for MA38BD
C
C       pgiven:         true if Cperm and Rperm are defined on input
C       if pgiven then:
C          Cperm (1..pn):       col permutation to BTF, n = pn
C          Rperm (1..pn):       row permutation to BTF
C       else
C          Cperm (1..pn):       unaccessed pn = 1
C          Rperm (1..pn):       unaccessed
C
C       II (1..nz+cscal*n):             expanded column-form, see below
C       II (nz+cscal*n+1..isize):       undefined on input
C       XX (1..nz):                     expanded column-form, see below
C       XX (nz+1..xsize):               undefined on input

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       WiR (1..n)
C       WiC (1..n)
C       WpR (1..n)
C       WpC (1..n)
C       Wm (1..n)
C       Head (n)
C       Rp (1..n+dn)
C       Wr (1..n+dn)
C       Wc (1..n+dn)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       II (1..itail-1):        undefined on output
C       II (itail..isize):      LU factors of this matrix, see below
C       XX (1..xtail-1):        undefined on output
C       XX (xtail..xsize):      LU factors of this matrix, see below
C
C       Info:           integer informational output, see MA38AD
C       Rinfo:          real informational output, see MA38AD
C       if pgiven:
C          Cperm (1..n): the final col permutations, including BTF
C          Rperm (1..n): the final row permutations, including BTF
C
C       WiC (1..n):     row permutations, not including BTF
C       WiR (1..n):     column permutations, not including BTF
C
C       iuse:           memory usage in Index
C       xuse:           memory usage in Value
C       rmax:           maximum ludegr seen so far (see below)
C       cmax:           maximum ludegc seen so far (see below)
C       totnlu:         total number of LU arrowheads constructed so far
C       xrmax:          maximum real memory usage for MA38BD
C       xruse:          current real memory usage for MA38BD

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38ED
C       subroutines called:     MA38ND, MA38GD, DGEMV, DGEMM
C       functions called:       IDAMAX, ABS, MAX, MIN
        INTEGER IDAMAX
        INTRINSIC ABS, MAX, MIN
        EXTERNAL MA38GD, MA38ND, DGEMM, DGEMV, IDAMAX

C=======================================================================
C  DESCRIPTION OF DATA STRUCTURES:
C=======================================================================

C-----------------------------------------------------------------------
C  Column/element/arrowhead pointers:
C-----------------------------------------------------------------------
C
C  The Cp (1..n) array contains information about non-pivotal columns
C
C       p = Cp (j)
C       if (p = 0) then j is a pivotal column
C       else i is a non-pivotal column
C
C  The Rp (1..n) array contains information about non-pivotal rows,
C  unassembled frontal matrices (elements), and the LU arrowheads
C
C       p = Rp (i)
C       if (i > n) then
C          i is an artificial frontal matrix (a dense column)
C          if (p = 0) then i is assembled, else unassembled
C       else if (p = 0) then i is pivotal but not element/arrowhead
C       else if (Wc (i) >= 0 and Wc (i) <= n) then
C          i is a non-pivotal row
C       else if (Wc (i) = -(n+dn+2)) then
C          i is a pivotal row, an assembled element, and an LU arrowhead
C       else i an unassembled element

C-----------------------------------------------------------------------
C  Matrix being factorized:
C-----------------------------------------------------------------------
C
C    Each column is stored in II and XX:
C    -----------------------------------
C
C       if j is a non-pivotal column, pc = Cp (j):
C
C       csiz = II (pc) size of the integer data structure for col j,
C                        including the cscal scalars
C       cdeg = II (pc+1) degree of column j
C       cxp  = II (pc+2) pointer into XX for numerical values
C       next = II (pc+3) pointer to next block of memory in XX
C       prev = II (pc+4) pointer to previous block of memory in XX
C       celn = II (pc+5) number of elements in column j element list
C       clen = II (pc+6) number of original entries in column j
C       cnxt = II (pc+7) next column with same degree as col j
C       cprv = II (pc+8) previous column with same degree as col j
C       cep = (pc+9) pointer to start of the element list
C       II (cep ... cep + 2*celn - 1)
C                       element list (e,f) for the column
C       II (cep + 2*celn ... pc + csiz - clen - 1)
C                       empty
C       II (pc + csiz - clen ... pc + csiz - 1)
C                       row indices of original nonzeros in the column
C       XX (xp ... xp + clen - 1)
C                       numerical values of original nonzeros in the col
C
C       if cdeg = II (pc+1) = -(n+2), then this is a singular column
C       if cdeg = -1, then this column is deallocated
C
C    Each row is stored in II only:
C    ------------------------------
C
C       if i is a non-pivotal row, pr = Rp (i)
C
C       rsiz = II (pr) size of the integer data structure for row i,
C                        including the rscal scalars
C       rdeg = II (pr+1) degree of row i
C       reln = Wr (i) number of elements in row i element list
C       rlen = Wc (i) number of original entries in row i
C       rep  = (pr+2) pointer to start of the element list
C       II (rep ... rep + 2*reln - 1)
C                       element list (e,f) for the row
C       II (rep + 2*reln ... pr + rsiz - rlen - 1)
C                       empty
C       II (pr + rsiz - rlen ... pr + rsiz - 1)
C                       column indices of original nonzeros in the row
C
C       if rdeg = -1, then this row is deallocated

C-----------------------------------------------------------------------
C  Frontal matrices
C-----------------------------------------------------------------------
C
C   Each unassembled frontal matrix (element) is stored as follows:
C       total size: fscal integers, (fdimr*fdimc) reals
C
C       if e is an unassembled element, ep = Rp (e), and e is also
C       the first pivot row in the frontal matrix.
C
C       fluip  = II (ep)        pointer to LU arrowhead in II
C       fdimc  = II (ep+1)      column dimension of contribution block
C       fxp    = II (ep+2)      pointer to contribution block in XX
C       next   = II (ep+3)      pointer to next block in XX
C       prev   = II (ep+4)      pointer to previous block in XX
C       fleftr = II (ep+5)      number of unassembled rows
C       fleftc = II (ep+6)      number of unassembled columns
C       fextr = Wr (e) - w0     external row degree of the frontal mtx
C       fextc = Wc (e) - w0     external col degree of the frontal mtx
C       XX (fxp ... )
C               a 2-dimensional array, C (1..fdimc, 1..fdimr).
C               note that fdimr is not kept (it is not needed,
C               except for the current frontal).  If this is not the
C               current frontal matrix, then luip points to the
C               corresponding LU arrowhead, and the contribution block
C               is stored in C (1..ludegc, 1..ludegr) in the
C               C (1..fdimc, ...) array.
C
C               If memory is limited, garbage collection will occur.
C               In this case, the C (1..fdimc, 1..fdimr) array is
C               compressed to be just large enough to hold the
C               unassembled contribution block,
C               C (1..ludegc, 1..ludegr).

C-----------------------------------------------------------------------
C  Artificial frontal matrices
C-----------------------------------------------------------------------
C
C   An artificial frontal matrix is an original column that is treated
C   as a c-by-1 frontal matrix, where c is the number of original
C   nonzeros in the column.  Dense columns (c > dsiz) are treated this
C   way.  An artificial frontal matrix is just the same as a frontal
C   matrix created by the elimination of one or more pivots, except
C   that there is no corresponding LU arrowhead.  The row and column
C   patterns are stored in:
C
C       ep = Rp (e), where e = n+1 .. n+dn, where there are dn
C                    artificial frontal matrices.
C
C       lucp = (ep+9)   pointer to row pattern (just one column index)
C       lurp = (ep+8) pointer to column pattern (fdimc row indices)

C-----------------------------------------------------------------------
C  Current frontal matrix
C-----------------------------------------------------------------------
C
C  ffxp points to current frontal matrix (contribution block and LU
C  factors).  For example, if fflefc = 4, fflefr = 6, k = 3, and
C  gro = 2.0, then "x" is a term in the contribution block, "l" in L1,
C  "u" in U1, "L" in L2, "U" in U2, and "." is unused.  XX (fxp) is "X".
C  The first 3 pivot values (diagonal entries in U1) are 1,2, and 3.
C  For this frontal matrix, ffdimr = 12 (the number of columns), and
C  ffdimc = 8 (the number of rows).  The frontal matrix is
C  ffdimc-by-ffdimr
C
C                             |----------- col 1 of L1 and L2, etc.
C                             V
C       X x x x x x . . . L L L
C       x x x x x x . . . L L L
C       x x x x x x . . . L L L
C       x x x x x x . . . L L L
C       . . . . . . . . . . . .
C       U U U U U U . . . 3 l l         <- row 3 of U1 and U2
C       U U U U U U . . . u 2 l         <- row 2 of U1 and U2
C       U U U U U U . . . u u 1         <- row 1 of U1 and U2

C-----------------------------------------------------------------------
C  LU factors
C-----------------------------------------------------------------------
C
C   The LU factors are placed at the tail of II and XX.  If this routine
C   is factorizing a single block, then this decription is for the
C   factors of the single block:
C
C       II (itail):             xtail = start of LU factors in XX
C       II (itail+1):           nlu = number of LU arrowheads
C       II (itail+2):           npiv = number of pivots
C       II (itail+3):           maximum number of rows in any
C                               contribution block (max ludegc)
C       II (itail+4):           maximum number of columns in any
C                               contribution block (max ludegr)
C       II (itail+5..itail+nlu+4): LUp (1..nlu) array, pointers to each
C                               LU arrowhead, in order of their
C                               factorization
C       II (itail+nlu+5...isize):integer info. for LU factors
C       XX (xtail..xsize):      real values in LU factors
C
C   Each LU arrowhead is stored as follows:
C   ---------------------------------------
C
C       total size: (7 + ludegc + ludegr + nsons) integers,
C                   (luk**2 + ludegc*luk + luk*ludegc) reals
C
C       If e is an LU arrowhead, then luip = Rp (e), and luip >= itail.
C       When MA38FD returns, then luip is given by luip =
C       II (itail+s+1), where s = 1..nlu is the position of the LU
C       arrowhead in the LU factors (s=1,2,.. refers to the first,
C       second,.. LU arrowhead)
C
C       luxp   = II (luip) pointer to numerical LU arrowhead
C       luk    = II (luip+1) number of pivots in LU arrowhead
C       ludegr = II (luip+2) degree of last row of U (excl. diag)
C       ludegc = II (luip+3) degree of last col of L (excl. diag)
C       nsons  = II (luip+4) number of children in assembly DAG
C       ludimr = II (luip+5)
C       ludimc = II (luip+5) max front size is ludimr-by-ludimc,
C                       or zero if this LU arrowhead factorized within
C                       the frontal matrix of a prior LU arrowhead.
C       lucp   = (luip + 7)
C                       pointer to pattern of column of L
C       lurp   = lucp + ludegc
C                       pointer to patter of row of U
C       lusonp = lurp + ludegr
C                       pointer to list of sons in the assembly DAG
C       II (lucp ... lucp + ludegc - 1)
C                       row indices of column of L
C       II (lurp ... lurp + ludegr - 1)
C                       column indices of row of U
C       II (lusonp ... lusonp + nsons - 1)
C                       list of sons
C       XX (luxp...luxp + luk**2 + ludegc*luk + luk*ludegr - 1)
C                       pivot block (luk-by-luk) and the L block
C                       (ludegc-by-luk) in a single (luk+ludegc)-by-luk
C                       array, followed by the U block in a
C                       luk-by-ludegr array.
C
C   Pivot column/row pattern (also columns/rows in contribution block):
C       If the column/row index is negated, the column/row has been
C       assembled out of the frontal matrix into a subsequent frontal
C       matrix.  After factorization, the negative flags are removed,
C       and the row/col indices are replaced with their corresponding
C       index in the permuted LU factors.
C
C   List of sons:
C       1 <= son <= n:           son an LUson
C       n+1 <= son <= 2n:        son-n is an Uson
C       2n+n <= son <= 3n:       son-2n is a Lson
C       during factorzation, a son is referred to by its first
C       pivot column.  After factorization, they are numbered according
C       to their order in the LU factors.

C-----------------------------------------------------------------------
C  Workspaces:
C-----------------------------------------------------------------------
C
C   WiR (e):  link list of sons of the current element
C       WiR (e) = -1 means that e is not in the list.
C       WiR (e) = next+n+2 means that "next" is the element after e.
C       The end of the list is marked with WiR (e) = -(n+2).
C       sonlst points to the first element in the list, or 0 if
C       the sonlst is empty.
C
C   WiR (row), WiC (col):  used for pivot row/col offsets:
C
C       If WiR (row) >= 0 then the row is in the current
C       column pattern.  Similarly for WiC (col).
C
C       If WiR (row) is set to "empty" (<= -1), then
C       the row is not in the current pivot column pattern.
C
C       Similarly, if WiC (col) is set to -2, then the column is
C       not in the current pivot row pattern.
C
C       If WiC (col) = -1 then col is pivotal
C
C       After factorization, WiR/C holds the pivot permutations.
C
C   WpR/C (1..n):  the first part is used for the current frontal
C           matrix pattern.  During factorization, the last part holds
C           a stack of the row and column permutations (WpR/C (n-k+1)
C           is the k-th pivot row/column).
C
C   Head (1..n):        degree lists for columns.  Head (d) is the
C                       first column in list d with degree d.
C                       The cnxt and cprv pointers are stored in the
C                       column data structure itself.
C                       mindeg is the least non-empty list
C
C   Wm (1..n):          various uses
C   Wj (1..degc) or Wj (1..fdegc):      offset in pattern of a son

C-----------------------------------------------------------------------
C  Memory allocation in II and XX:
C-----------------------------------------------------------------------
C
C   II (1..ihead):      rows and columns of active submatrix, and
C                       integer information for frontal matrices.
C   XX (1..xhead):      values of original entries in columns of
C                       matrix, values of contribution blocks, followed
C                       by the current frontal matrix.
C
C   mhead:              a pointer to the first block in the head of
C                       XX.  Each block (a column or frontal matrix)
C                       contains a next and prev pointer for this list.
C                       If the list is traversed starting at mhead,
C                       then the pointers to the reals (cxp or fxp)
C                       will appear in strictly increasing order.
C                       Note that the next, prev, and real pointers
C                       are in II.  next and prev point to the next
C                       and previous block in II, and the real pointer
C                       points to the real part in XX.
C
C   mtail:              the end of the memory list.

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER SWPCOL, SWPROW, FDIMC, K0, COLPOS, ROWPOS, ROW2, RDEG2,
     $          P, I, J, FFROW, PIVROW, PIVCOL, LUDEGR, LUDEGC, E1,
     $          FXP, LURP, LUCP, IP, NEXT, FFLEFR, PC, MNEXT, MPREV,
     $          FFLEFC, FEDEGR, FEDEGC, K, XUDP, XDP, XSP, XLP, S, COL2,
     $          BESTCO, COL, E, ROW, COST, SRCHED, PR, F1, RSCAN, REP,
     $          KLEFT1, FFSIZE, FFXP, W0, FFDIMR, FFDIMC, KLEFT, XLDP,
     $          EP, SCAN1, SCAN2, SCAN3, SCAN4, NZL, NZU, DEGC, CEP
        INTEGER MINDEG, NSRCH, NPIV, ESON, LUIP1, DNZ, IWORST, WXP,
     $          NB, LUPP, NLU, NSONS, INEED, XNEED, LDIMC, LXP, RLEN2,
     $          RSIZ, LSONS, SONLST, XHEAD, IHEAD, DELN, DLEN,
     $          SLIST, XP, LUIP, RDEG, CDEG1, PFREE, XFREE, CDEG2,
     $          F, CDEG, MTAIL, MHEAD, RSIZ2, CSIZ2, IP2, MAXDR, MAXDC,
     $          XS, IS, LUXP, FSP, FLP, FDP, JJ, USONS, NDN, P2,
     $          CSIZ, CELN, CLEN, RELN, RLEN, UXP, PC2, PR2
        INTEGER CNXT, CPRV, CXP, FLUIP, LUSONP, FLEFTR, FLEFTC,
     $          FMAXR, FMAXC, SLOTS, LIMIT, RSCAL, CSCAL, FSCAL, EXTRA,
     $          FMAX, MINMEM, DUMMY1, DUMMY2, DUMMY3, DUMMY4, W1
        LOGICAL SYMSRC, PFOUND, MOVELU, OKCOL, OKROW, BETTER
        DOUBLE PRECISION
     $          TOLER, MAXVAL, RELPT, GRO, ONE, ZERO, X
        PARAMETER (ONE = 1.0D0, ZERO = 0.0D0,
     $          RSCAL = 2, CSCAL = 9, FSCAL = 7,
     $          MINMEM = 24)
        INTEGER INTZER

C  Current element and working array, C:
C  -------------------------------------
C  ffxp:    current working array is in XX (ffxp ... ffxp+ffsize-1)
C  ffsize:  size of current working array in XX
C  ffdimr:  row degree (number of columns) of current working array
C  ffdimc:  column degree (number of rows) of current working array
C  fflefr:  row degree (number of columns) of current contribution block
C  fflefc:  column degree (number of rows) of current contribution block
C  fmaxr:   max row degree (maximum front size is fmaxr-by-fmaxc)
C  fmaxc:   max col degree (maximum front size is fmaxr-by-fmaxc)
C  fedegr:  extended row degree
C  fedegc:  extended column degree
C  ffrow:   current element being factorized (a pivot row index)
C  pivrow:  current pivot row index
C  pivcol:  current pivot column index
C  e1:      first pivot row in the frontal matrix
C  gro:     frontal matrix amalgamation growth factor
C  usons:   pointer to a link list of Usons, in Wc, assembled this SCAN3
C  lsons:   pointer to a link list of Lsons, in Wr, assembled this SCAN4
C  sonlst:  pointer to a link list of sons, in WiR, of current element
C  swpcol:  the non-pivotal column to be swapped with pivot column
C  swprow:  the non-pivotal row to be swapped with pivot row
C  colpos:  position in WpR of the pivot column
C  rowpos:  position in WpC of the pivot row
C  k:       current pivot is kth pivot of current element
C  k0:      contribution block, C, has been updated with pivots 1..k0
C
C  LU arrowhead (a factorized element):
C  ------------------------------------
C  movelu:  true if a new LU arrowhead is to be created
C  luip:    current element is in II (luip ...)
C  luip1:   first element from current frontal matrix in II (luip1...)
C  ludegc:  degree of pivot column (excluding pivots themselves)
C  ludegr:  degree of pivot row (excluding pivots themselves)
C  lucp:    pattern of col(s) of current element in II (lucp...)
C  lurp:    pattern of row(s) of current element in II (lurp...)
C  lusonp:  list of sons of current element is in II (lusonp...)
C  nsons:   number of sons of current element
C  ldimc:   column dimension (number of rows) of [L1\U1 L2] block
C  luxp:    numerical values of LU arrowhead stored in XX (luxp ...)
C  lxp:     L2 block is stored in XX (lxp ...) when computed
C  uxp:     U2 block is stored in XX (uxp ...) when computed
C  nzu:     nonzeros above diagonal in U in current LU arrowhead
C  nzl:     nonzeros below diagonal in L in current LU arrowhead
C
C  Son, or element other than current element:
C  -------------------------------------------
C  e:       an element
C  eson:    an element
C  s:       a renumbered element (1..nlu) for MA38CD and MA38BD
C  ep:      frontal matrix integer data struct. in II (ep...ep+fscal-1)
C  fscal:   = 7, size of frontal matrix data structure
C  fluip:   LU arrowhead of e is in II (fluip ...)
C  fxp:     contribution block of son is in XX (fxp ...)
C  fdimc:   leading dimension of contribution block of e
C  lucp:    pattern of col(s) of e in II (lucp...)
C  lurp:    pattern of row(s) of e in II (lurp...)
C  ludegr:  row degree of contribution block of e
C  ludegr:  column degree of contribution block of e
C  maxdr:   maximum ludegr for any LU arrowhead, for MA38BD
C  maxdc:   maximum ludegc for any LU arrowhead, for MA38BD
C  degc:    compressed column offset vector of son is in Wj/Wm (1..degc)
C  fleftr:  remaining row degree (number of columns) of a contrib. block
C  fleftc:  remaining column degree (number of rows) of a contrib. block
C  xudp:    pointer to a column of a prior contribution block
C  xldp:    pointer to a row of a prior contribution block
C
C  Memory allocation:
C  ------------------
C  mhead:   head pointer for link list of blocks in XX
C  mtail:   tail pointer for link list of blocks in XX
C  mprev:   previous block, II (p+4), of the block located at p
C  mnext:   next block, II (p+3), of the block located at p
C  pfree:   II (pfree+2) is the largest known free block in XX
C  xfree:   size of largest known free block in XX
C  xhead:   XX (1..xhead-1) is in use, XX (xhead ..xtail-1) is free
C  xtail:   XX (xtail..xsize) is in use, XX (xhead ..xtail-1) is free
C  xneed:   bare minimum memory currently needed in XX
C  ihead:   II (1..ihead-1) is in use, II (ihead ..itail-1) is free
C  itail:   II (itail..isize) is in use, II (ihead ..itail-1) is free
C  ineed:   bare minimum memory currently needed in II
C  iworst:  worst possible current integer memory required
C  xs:      size of a block of memory in XX
C  is:      size of a block of memory in II
C  wxp:     pointer to a temporary workspace in XX (wxp ... )
C  slots:   number of slots added to element lists during garbage coll.
C  minmem:  smallest isize allowed
C
C  Wr and Wc flag arrays:
C  ----------------------
C  w0:      marker value for Wr (1..n) and Wc (1..n) arrays
C  w1:      used to test for possible integer overflow in flags
C  fmax:    largest row/col degree of an element seen so far
C
C  A column:
C  ---------
C  pc:      pointer to a column, in II (pc...)
C  pc2:     pointer to a column, in II (pc2...)
C  csiz:    size of integer data structure of a column
C  csiz2:   size of integer data structure of a column
C  cscal:   = 9, number of scalars in data structure of a column
C  cdeg:    degree of a column
C  cdeg1:   degree of a column
C  cdeg2:   degree of a column
C  celn:    number of elements in the element list of a column
C  clen:    number of original entries that remain in a column
C  cnxt:    next column with same degree as this column
C  cprv:    previous column with same degree as this column
C  cep:     pointer to the element list of a column
C  cxp:     pointer to the numerical values in a column
C  limit:   maximum size for row/col data structure (excl. scalars)
C
C  Dense columns:
C  --------------
C  dnz:     number of original entries that reside in "dense" columns
C  dn:      number of "dense" columns
C  ndn:     n + dn
C  extra:   number of extra slots to add to reconstructed "dense" cols
C
C  A row:
C  ------
C  pr:      pointer to a row, in II (pr...)
C  pr2:     pointer to a row, in II (pr2...)
C  rsiz:    size of integer data structure of a row
C  rsiz2:   size of integer data structure of a row
C  rscal:   = 2, number of scalars in data structure of a row
C  rdeg:    degree of a row
C  rdeg2:   degree of a row
C  reln:    number of elements in the element list of a row
C  rlen:    number of original entries that remain in a row
C  rlen2:   number of original entries that remain in a row
C  rep:     pointer to the element list of a row
C
C  Pivot search:
C  -------------
C  cost:    approximate Markowitz-cost of the current candidate pivot
C  bestco:  best approximate Markowitz-cost seen so far
C  srched:  number of non-singular candidates searched so far
C  mindeg:  minimum degree of columns in active submatrix
C  nsrch:   maximum number of columns to search
C  slist:   pointer to a link list of searched columns, in II
C  symsrc:  true if attempting to preserve symmetry
C  pfound:  true if pivot found during local search
C  okcol:   true if candidate pivot column is acceptable, so far
C  okrow:   true if candidate pivot row is acceptable, so far
C  toler:   pivot tolerance; abs(pivot) must be >= toler
C  maxval:  maximum absolute value in a candidate pivot column
C  relpt:   relative pivot tolerance (Cntl (1))
C  npiv:    number of pivots factorized so far, incl. current element
C  kleft:   number of rows/columns remaining in active submatrix
C  kleft1:  kleft - 1
C  better:  if true, then candidate is better than the prior candidate
C
C  Assembly:
C  ---------
C  f1:      degree prior to assembly next item
C  f:       offset into an element
C  rscan:   skip row assembly if more than rscan original entries
C  scan1:   start SCAN1 at WpC (scan1 ... fflefc)
C  scan2:   start SCAN2 at WpR (scan2 ... fflefr)
C  scan3:   start SCAN3 at WpR (scan3 ... fflefr)
C  scan4:   start SCAN4 at WpC (scan4 ... fflefc)
C  deln:    number of (e,f) tuples to delete from an element list
C  dlen:    number of original entries to delete from a row/col
C
C  Allocated arrays:
C  -----------------
C  lupp:    LUp (1..nlu) array located in II (lupp...lupp+nlu-1)
C  nlu:     number of LU arrowheads
C
C  Other:
C  ------
C  xdp:     destination pointer, into XX
C  xsp:     source pointer, into XX
C  xlp:     pointer into XX of location of last row/col in C
C  xp:      pointer into XX
C  ip:      pointer into II
C  ip2:     pointer into II
C  p2:      pointer into II
C  fsp:     source pointer, into XX
C  fsp:     destination pointer, into XX
C  flp:     last row/column in current contribution is in XX (flp...)
C  col,col2: a column index
C  row,row2: a row index
C  nb:      block size for tradeoff between Level-2 and Level-3 BLAS
C  p, i, j, k, x:  various uses
C  jj:      loop index
C  next:    next pointer, for a link list
C  dummy1:  dummy loop index for main factorization loop
C  dummy2:  dummy loop index for global pivot search loop
C  dummy3:  dummy loop index for outer frontal matrix factorization loop
C  dummy4:  dummy loop index for inner frontal matrix factorization loop

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C       ----------------------------------------------------------------
C       get control parameters and initialize various scalars
C       ----------------------------------------------------------------

        NSRCH = MAX (1, ICNTL (5))
        SYMSRC = ICNTL (6) .NE. 0
        NB = MAX (1, ICNTL (7))
        RELPT = MAX (ZERO, MIN (CNTL (1), ONE))
        GRO = MAX (ONE, CNTL (2))
        NDN = N + DN
        W0 = NDN + 2
C       currently: w0 = n+dn+2 < 2n+2
        KLEFT = N
        NPIV = 0
        NLU = 0
        MINDEG = 1
        FMAX = 1
        IHEAD = NZ + CSCAL*N + 1
        XHEAD = NZ + 1
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1
C       Cp (1) must equal 1, the first block
        XFREE = -1
        PFREE = 0
C       make sure integer space is at least of size minmem (simplifies
C       link list management and memory management)
        INFO (19) = MAX (INFO (19), IUSE+MINMEM)
        IF (IHEAD.GT.ITAIL.OR.ISIZE.LT.MINMEM.OR.XHEAD.GT.XTAIL) THEN
C          error return, if not enough integer and/or real memory:
           GO TO 9000
        ENDIF
        BESTCO = 0
        LIMIT = N + 2*NDN
        LSONS = NDN + 1
        USONS = NDN + 1

C       ----------------------------------------------------------------
C       initialize workspaces
C       ----------------------------------------------------------------

        DO 10 I = 1, N
           WIR (I) = -1
           WIC (I) = -2
           HEAD (I) = 0
           WC (I) = 0
           WR (I) = 0
10      CONTINUE

C       ----------------------------------------------------------------
C       initialize the link list for keeping track of real memory usage
C       ----------------------------------------------------------------

        MHEAD = 0
        MTAIL = 0
        DO 20 COL = 1, N
           PC = CP (COL)
           CLEN = II (PC+6)
           IF (CLEN .GT. 0) THEN
C             place the column in the link list of blocks in XX
              IF (MHEAD .EQ. 0) MHEAD = PC
              II (PC+4) = MTAIL
              II (PC+3) = 0
              IF (MTAIL .NE. 0) II (MTAIL+3) = PC
              MTAIL = PC
           ELSE
              II (PC+2) = 0
              II (PC+4) = 0
              II (PC+3) = 0
           ENDIF
20      CONTINUE

C       ----------------------------------------------------------------
C       convert dense columns to a-priori contribution blocks and
C       get the count of nonzeros in each row
C       ----------------------------------------------------------------

        E = N
        DNZ = 0
        DO 50 COL = 1, N
           PC = CP (COL)
           CLEN = II (PC+6)
           CEP = (PC+9)
           IF (CLEN .GT. DSIZ) THEN
C             this is a dense column - add to element list length
              DNZ = DNZ + CLEN
              DO 30 IP = CEP, CEP + CLEN - 1
                 ROW = II (IP)
                 WR (ROW) = WR (ROW) + 1
30            CONTINUE
C             convert dense column (in place) into a frontal matrix
              E = E + 1
              EP = PC
              RP (E) = EP
              FDIMC = CLEN
              FLEFTC = CLEN
              FLEFTR = 1
              II (EP+1) = FDIMC
              II (EP+5) = FLEFTR
              II (EP+6) = FLEFTC
              WR (E) = W0-1
              WC (E) = W0-1
              LURP = (EP+8)
              II (LURP) = COL
              FMAX = MAX (FMAX, FLEFTC)
           ELSE
C             this is a sparse column - add to orig entry length
              DO 40 IP = CEP, CEP + CLEN - 1
                 ROW = II (IP)
                 WC (ROW) = WC (ROW) + 1
40            CONTINUE
           ENDIF
50      CONTINUE

C       ----------------------------------------------------------------
C       get memory for row-oriented form, and dense column element lists
C       ----------------------------------------------------------------

        PR = IHEAD
        CSIZ = CSCAL + 2
        IS = (NZ + RSCAL*N + DNZ) + (DN * CSIZ)
        IHEAD = IHEAD + IS
        IUSE = IUSE + IS
        INEED = IUSE
        XNEED = XUSE
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .GT. ITAIL) THEN
C          error return, if not enough integer memory:
           GO TO 9000
        ENDIF

C       ----------------------------------------------------------------
C       if memory is available, add up to dsiz+6 extra slots in the
C       reconstructed dense columns to allow for element list growth
C       ----------------------------------------------------------------

        IF (DN .GT. 0) THEN
           EXTRA = MIN ((ITAIL - IHEAD) / DN, DSIZ + 6)
           CSIZ = CSIZ + EXTRA
           IS = DN * EXTRA
           IHEAD = IHEAD + IS
           IUSE = IUSE + IS
           INFO (18) = MAX (INFO (18), IUSE)
        ENDIF

C       ----------------------------------------------------------------
C       construct row pointers
C       ----------------------------------------------------------------

        DO 60 ROW = 1, N
           RP (ROW) = PR
           REP  = (PR+2)
           RELN = WR (ROW)
           RLEN = WC (ROW)
           RSIZ = 2*RELN + RLEN + RSCAL
           II (PR) = RSIZ
           RDEG = RELN + RLEN
           II (PR+1) = RDEG
           WM (ROW) = REP
           PR = PR + RSIZ
60      CONTINUE

C       ----------------------------------------------------------------
C       construct row element lists for dense columns
C       ----------------------------------------------------------------

        PC = PR
        DO 80 E = N+1, N+DN
           EP = RP (E)
           LUCP = (EP+9)
           FDIMC = II (EP+1)
CFPP$ NODEPCHK L
           DO 70 F = 0, FDIMC - 1
              ROW = II (LUCP+F)
              II (WM (ROW)    ) = E
              II (WM (ROW) + 1) = F
              WM (ROW) = WM (ROW) + 2
70         CONTINUE
C          re-construct dense columns as just an element list,
C          containing a single element tuple (e,f), where f = 0
           LURP = (EP+8)
           COL = II (LURP)
           CP (COL) = PC
           II (PC) = CSIZ
           CDEG = FDIMC
           II (PC+1) = CDEG
           II (PC+2) = 0
           II (PC+4) = 0
           II (PC+3) = 0
           II (PC+5) = 1
           II (PC+6) = 0
           II (PC+7) = 0
           II (PC+8) = 0
C          store the (e,0) tuple:
           CEP = (PC+9)
           II (CEP  ) = E
           II (CEP+1) = 0
           PC = PC + CSIZ
80      CONTINUE

C       ----------------------------------------------------------------
C       construct the nonzero pattern of the row-oriented form
C       ----------------------------------------------------------------

        DO 100 COL = 1, N
           PC = CP (COL)
           CEP = (PC+9)
           CLEN = II (PC+6)
CFPP$ NODEPCHK L
           DO 90 P = CEP, CEP + CLEN - 1
              ROW = II (P)
              II (WM (ROW)) = COL
              WM (ROW) = WM (ROW) + 1
90         CONTINUE
100     CONTINUE

C       count the numerical assembly of the original matrix
        RINFO (2) = RINFO (2) + NZ

C       ----------------------------------------------------------------
C       initialize the degree lists
C       ----------------------------------------------------------------

C       do so in reverse order to try to improve pivot tie-breaking
        DO 110 COL = N, 1, -1
           PC = CP (COL)
           CDEG = II (PC+1)
           IF (CDEG .LE. 0) THEN
C             empty column - remove from pivot search
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ELSE
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) II (CP (CNXT)+8) = COL
              HEAD (CDEG) = COL
           ENDIF
110     CONTINUE

C=======================================================================
C=======================================================================
C  MAIN FACTORIZATION LOOP [
C=======================================================================
C=======================================================================

        DO 1540 DUMMY1 = 1, N
C       (this loop is not indented due to its length)

C       ----------------------------------------------------------------
C       factorization is done if N pivots have been found
C       ----------------------------------------------------------------

        IF (NPIV .GE. N) THEN
           GO TO 2000
        ENDIF

C=======================================================================
C  Global pivot search, and initialization of a new frontal matrix [
C=======================================================================

        IF (MTAIL .NE. 0 .AND. II (MTAIL+6) .EQ. 0) THEN
C          tail block is free, delete it
           XP = II (MTAIL+2)
           XUSE = XUSE - (XHEAD - XP)
           XHEAD = XP
           IF (MTAIL .EQ. PFREE) THEN
              PFREE = 0
              XFREE = -1
           ENDIF
           MTAIL = II (MTAIL+4)
           IF (MTAIL .NE. 0) THEN
              II (MTAIL+3) = 0
           ELSE
C             singular matrix.  No columns or contribution blocks left.
              MHEAD = 0
           ENDIF
        ENDIF

C=======================================================================
C  Global pivot search:  find pivot row and column
C=======================================================================

        NSONS = 0
        SONLST = 0
        SRCHED = 0
        PIVCOL = 0
        SLIST = 0

        DO 255 DUMMY2 = 1, N

C          -------------------------------------------------------------
C          get col from column upper-bound degree list
C          -------------------------------------------------------------

           COL = 0
           DO 140 CDEG = MINDEG, N
              COL = HEAD (CDEG)
              IF (COL .NE. 0) THEN
C                exit out of loop if column found:
                 GO TO 150
              ENDIF
140        CONTINUE
           IF (COL .EQ. 0) THEN
C             exit out of loop if column not found (singular matrix):
              GO TO 260
           ENDIF
C          loop exit label:
150        CONTINUE
           PC = CP (COL)
           CNXT = II (PC+7)
           IF (CNXT .NE. 0) THEN
              II (CP (CNXT)+8) = 0
           ENDIF
           HEAD (CDEG) = CNXT
           MINDEG = CDEG

C          -------------------------------------------------------------
C          construct candidate column in Wm and XX (wxp..wxp+cdeg-1)
C          -------------------------------------------------------------

           XS = CDEG
C          use Wm (1..cdeg) for pattern [
C          use XX (wxp..wxp+xs-1) as workspace for values [

           IF (XS .GT. XTAIL-XHEAD) THEN

              INFO (15) = INFO (15) + 1
              INTZER = 0
              CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                     II, ISIZE, IHEAD, IUSE,
     $                     CP, RP, DN, N, WIR, WIC, WR, WC,
     $                     INTZER, INTZER, INTZER, INTZER, .FALSE.,
     $                     PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C             at this point, iuse = ineed and xuse = xneed
              PC = CP (COL)
           ENDIF

           WXP = XHEAD
           XHEAD = XHEAD + XS
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (XHEAD .GT. XTAIL) THEN
C             error return, if not enough real memory:
              GO TO 9000
           ENDIF

C          -------------------------------------------------------------
C          assemble the elements in the element list
C          -------------------------------------------------------------

           CDEG = 0
           CEP = (PC+9)
           CELN = II (PC+5)
           DO 190 IP = CEP, CEP + 2*CELN - 2, 2
              E = II (IP)
              F = II (IP+1)
              EP = RP (E)
              FDIMC = II (EP+1)
              FXP = II (EP+2)
              IF (E .LE. N) THEN
                 FLUIP = II (EP)
                 LUDEGC = II (FLUIP+3)
                 LUCP = (FLUIP + 7)
              ELSE
                 LUDEGC = FDIMC
                 LUCP = (EP+9)
              ENDIF
              XP = FXP + F * FDIMC
C             split into 3 loops so that they all vectorize on a CRAY
              CDEG1 = CDEG
              DO 160 P = LUCP, LUCP + LUDEGC - 1
                 ROW = II (P)
                 IF (ROW .GT. 0) THEN
                    IF (WIR (ROW) .LE. 0) THEN
                       CDEG = CDEG + 1
                       WM (CDEG) = ROW
                    ENDIF
                 ENDIF
160           CONTINUE
              DO 170 I = CDEG1+1, CDEG
                 ROW = WM (I)
                 WIR (ROW) = I
                 XX (WXP+I-1) = ZERO
170           CONTINUE
CFPP$ NODEPCHK L
              DO 180 J = 0, LUDEGC - 1
                 ROW = II (LUCP+J)
                 IF (ROW .GT. 0) THEN
                    XX (WXP + WIR (ROW) - 1) =
     $              XX (WXP + WIR (ROW) - 1) + XX (XP+J)
                 ENDIF
180           CONTINUE
190        CONTINUE

C          -------------------------------------------------------------
C          assemble the original entries in the column
C          -------------------------------------------------------------

           CDEG1 = CDEG
           CLEN = II (PC+6)
           CSIZ = II (PC)
           IP = PC + CSIZ - CLEN
           CXP = II (PC+2)
CFPP$ NODEPCHK L
           DO 200 I = 0, CLEN - 1
              ROW = II (IP+I)
              WM (CDEG+1+I) = ROW
              XX (WXP+CDEG+I) = XX (CXP+I)
200        CONTINUE
           CDEG = CDEG + CLEN

C          -------------------------------------------------------------
C          update the degree of this column (exact, not upper bound)
C          -------------------------------------------------------------

           II (PC+1) = CDEG

C          Wm (1..cdeg) holds the pattern of col being searched.
C          XX (wxp..wxp+cdeg-1) holds the numerical values of col being
C          searched.  WiR (Wm (1..cdeg1)) is 1..cdeg1.

C          -------------------------------------------------------------
C          find the maximum absolute value in the column
C          -------------------------------------------------------------

           MAXVAL = ABS (XX (WXP - 1 + IDAMAX (CDEG, XX (WXP), 1)))
           RINFO (3) = RINFO (3) + CDEG
           TOLER = RELPT * MAXVAL
           RDEG = N+1

C          -------------------------------------------------------------
C          look for the best possible pivot row in this column
C          -------------------------------------------------------------

           IF (CDEG .NE. 0 .AND. MAXVAL .GT. ZERO) THEN
              IF (SYMSRC) THEN
C                prefer symmetric pivots, if numerically acceptable
                 ROW = COL
                 ROWPOS = WIR (ROW)
                 IF (ROWPOS .LE. 0) THEN
C                   diagonal may be in original entries
                    DO 210 I = CDEG1 + 1, CDEG1 + CLEN
                       IF (WM (I) .EQ. ROW) THEN
                          ROWPOS = I
C                         exit out of loop if symmetric pivot found:
                          GO TO 220
                       ENDIF
210                 CONTINUE
C                   loop exit label:
220                 CONTINUE
                 ENDIF
                 IF (ROWPOS .GT. 0) THEN
C                   diagonal entry exists in the column pattern
                    X = ABS (XX (WXP-1+ROWPOS))
                    IF (X .GE. TOLER .AND. X .GT. ZERO) THEN
C                      diagonal entry is numerically acceptable
                       PR = RP (ROW)
                       RDEG = II (PR+1)
                    ENDIF
                 ENDIF
              ENDIF
              IF (RDEG .EQ. N+1) THEN
C                Continue searching - no diagonal found or sought for.
C                Minimize row degree subject to abs(value) constraints.
                 ROW = N+1
                 DO 230 I = 1, CDEG
                    ROW2 = WM (I)
                    PR = RP (ROW2)
                    RDEG2 = II (PR+1)
C                   among those numerically acceptable rows of least
C                   (upper bound) degree, select the row with the
C                   lowest row index
                    BETTER = RDEG2 .LT. RDEG .OR.
     $                      (RDEG2 .EQ. RDEG .AND. ROW2 .LT. ROW)
                    IF (BETTER) THEN
                       X = ABS (XX (WXP-1+I))
                       IF (X.GE.TOLER .AND. X.GT.ZERO) THEN
                          ROW = ROW2
                          RDEG = RDEG2
                          ROWPOS = I
                       ENDIF
                    ENDIF
230              CONTINUE
              ENDIF
           ENDIF

C          -------------------------------------------------------------
C          deallocate workspace
C          -------------------------------------------------------------

           XHEAD = XHEAD - XS
           XUSE = XUSE - XS
           XNEED = XNEED - XS
C          done using XX (wxp...wxp+xs-1) ]

C          -------------------------------------------------------------
C          reset work vector
C          -------------------------------------------------------------

           DO 240 I = 1, CDEG1
              WIR (WM (I)) = -1
240        CONTINUE

C          -------------------------------------------------------------
C          check to see if a pivot column was found
C          -------------------------------------------------------------

           IF (RDEG .EQ. N+1) THEN

C             ----------------------------------------------------------
C             no pivot found, column is zero
C             ----------------------------------------------------------

C             remove this singular column from any further pivot search
              CDEG = -(N+2)
              II (PC+1) = CDEG

           ELSE

C             ----------------------------------------------------------
C             save a list of the columns searched (with nonzero degrees)
C             ----------------------------------------------------------

              SRCHED = SRCHED + 1
              II (PC+7) = SLIST
              SLIST = COL

C             ----------------------------------------------------------
C             check if this is the best pivot seen so far
C             ----------------------------------------------------------

C             compute the true Markowitz cost without scanning the row
C             Wm (1..cdeg) holds pivot column, including pivot row index
C             Wm (rowpos) contains the candidate pivot row index
              COST = (CDEG - 1) * (RDEG - 1)
              IF (PIVCOL .EQ. 0 .OR. COST .LT. BESTCO) THEN
                 FFLEFC = CDEG
                 DO 250 I = 1, FFLEFC-1
                    WPC (I) = WM (I)
250              CONTINUE
C                remove the pivot row index from pivot column pattern
                 WPC (ROWPOS) = WM (FFLEFC)
                 PIVCOL = COL
                 PIVROW = ROW
                 BESTCO = COST
              ENDIF
           ENDIF

C          done using Wm (1..cdeg) for pattern ]
C          WpC (1..fflefc-1) holds pivot column (excl. pivot row index)

C          -------------------------------------------------------------
C          exit global pivot search if nsrch pivots have been searched
C          -------------------------------------------------------------

           IF (SRCHED .GE. NSRCH) THEN
              GO TO 260
           ENDIF

255     CONTINUE
C       exit label for loop 255:
260     CONTINUE

C=======================================================================
C  Quit early if no pivot found (singular matrix detected)
C=======================================================================

        IF (PIVCOL .EQ. 0) THEN
C          complete the column permutation vector in
C          WpC (n-npiv+1 ... n) in reverse order
           K = N - NPIV + 1
           DO 270 COL = 1, N
              IF (CP (COL) .NE. 0) THEN
C                this is a non-pivotal column
                 K = K - 1
                 WPC (K) = COL
                 CP (COL) = 0
              ENDIF
270        CONTINUE
C          complete the row permutation vector in
C          WpR (n-npiv+1 ... n) in reverse order
           K = N - NPIV + 1
           DO 280 ROW = 1, NDN
              IF (ROW .GT. N) THEN
C                this is an artificial frontal matrix
                 E = ROW
                 RP (E) = 0
              ELSE IF (RP (ROW) .NE. 0) THEN
                 RLEN = WC (ROW)
                 IF (RLEN .GE. 0 .AND. RLEN .LE. N) THEN
C                   this is a non-pivotal row
                    K = K - 1
                    WPR (K) = ROW
                    RP (ROW) = 0
                 ELSE IF (RLEN .NE. -(NDN+2)) THEN
C                   this is an unassembled element: convert to LU
                    E = ROW
                    EP = RP (ROW)
                    WR (E) = -(NDN+2)
                    WC (E) = -(NDN+2)
                    FLUIP = II (EP)
                    RP (E) = FLUIP
                 ENDIF
              ENDIF
280        CONTINUE
C          factorization is done, exit the main factorization loop:
           GO TO 2000
        ENDIF

C=======================================================================
C  Place the non-pivotal columns searched back in degree lists
C=======================================================================

        DO 300 I = 1, SRCHED
           COL = SLIST
           PC = CP (COL)
           SLIST = II (PC+7)
           IF (COL .NE. PIVCOL) THEN
              CDEG = II (PC+1)
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = COL
              ENDIF
              HEAD (CDEG) = COL
              MINDEG = MIN (MINDEG, CDEG)
           ENDIF
300     CONTINUE

C=======================================================================
C  Construct pivot row pattern
C=======================================================================

C       At this point, WiR (1..n) = -1 and WiC (1..n) is -2 for
C       nonpivotal columns and -1 for pivotal columns.
C       WiC (WpR (1..fflefr+1)) is set to zero in the code below.  It
C       will be set to the proper offsets in do 775, once ffdimc is
C       known (offsets are dependent on ffdimc, which is dependent on
C       fflefr calculated below, and the memory allocation).

C       ----------------------------------------------------------------
C       assemble the elements in the element list
C       ----------------------------------------------------------------

        PR = RP (PIVROW)
        FFLEFR = 0
        REP = (PR+2)
        RELN = WR (PIVROW)
        DO 330 IP = REP, REP + 2*RELN - 2, 2
           E = II (IP)
           EP = RP (E)
           IF (E .LE. N) THEN
              FLUIP = II (EP)
              LUCP = (FLUIP + 7)
              LUDEGR = II (FLUIP+2)
              LUDEGC = II (FLUIP+3)
              LURP = LUCP + LUDEGC
C             split into two loops so that they both vectorize on a CRAY
              F1 = FFLEFR
              DO 310 P = LURP, LURP + LUDEGR - 1
                 COL = II (P)
                 IF (COL .GT. 0) THEN
                    IF (WIC (COL) .EQ. -2) THEN
                       FFLEFR = FFLEFR + 1
                       WPR (FFLEFR) = COL
                    ENDIF
                 ENDIF
310           CONTINUE
              DO 320 I = F1+1, FFLEFR
                 WIC (WPR (I)) = 0
320           CONTINUE
           ELSE
C             this is an artifical element (a dense column)
              LURP = (EP+8)
              COL = II (LURP)
              IF (WIC (COL) .EQ. -2) THEN
                 FFLEFR = FFLEFR + 1
                 WPR (FFLEFR) = COL
                 WIC (COL) = 0
              ENDIF
           ENDIF
330     CONTINUE

C       ----------------------------------------------------------------
C       assemble the original entries in the pivot row
C       ----------------------------------------------------------------

        RSIZ = II (PR)
        RLEN = WC (PIVROW)
        DO 340 P = PR + RSIZ - RLEN, PR + RSIZ - 1
           COL = II (P)
           IF (WIC (COL) .EQ. -2) THEN
              FFLEFR = FFLEFR + 1
              WPR (FFLEFR) = COL
           ENDIF
340     CONTINUE
C       the exact degree of the pivot row is fflefr

C=======================================================================
C  Initialize the new frontal matrix
C=======================================================================

C       ffrow is the name of current frontal matrix
        FFROW = PIVROW
        E1 = PIVROW
        K = 1
        K0 = 0
        FFDIMR = MIN (KLEFT, INT (GRO * FFLEFR))
        FFDIMC = MIN (KLEFT, INT (GRO * FFLEFC))
        FMAXR = FFLEFR
        FMAXC = FFLEFC
        FFSIZE = FFDIMC * FFDIMR
        RSCAN = MAX (DSIZ, FFDIMR)

C       ----------------------------------------------------------------
C       compute the offsets for rows in the pivot column
C       and the offsets for columns in the pivot row
C       ----------------------------------------------------------------

        DO 350 I = 1, FFLEFC - 1
           WIR (WPC (I)) = I - 1
350     CONTINUE
        DO 360 I = 1, FFLEFR
           WIC (WPR (I)) = (I - 1) * FFDIMC
360     CONTINUE

C       ----------------------------------------------------------------
C       remove the pivot column index from the pivot row pattern
C       ----------------------------------------------------------------

        COL = WPR (FFLEFR)
        COLPOS = (WIC (PIVCOL)/FFDIMC)+1
        WPR (COLPOS) = COL
        WIC (COL) = WIC (PIVCOL)
        WIC (PIVCOL) = (FFDIMR - 1) * FFDIMC
        WIR (PIVROW) = FFDIMC - 1

C       ----------------------------------------------------------------
C       remove the pivot row/col from the nonzero count
C       ----------------------------------------------------------------

        FFLEFR = FFLEFR - 1
        FFLEFC = FFLEFC - 1

C       ----------------------------------------------------------------
C       allocate the working array, doing garbage collection if needed
C       also allocate space for a work vector of size ffdimc
C       ----------------------------------------------------------------

        IF (FFSIZE + FFDIMC .GT. XTAIL-XHEAD) THEN
           INFO (15) = INFO (15) + 1
           INTZER = 0
           CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                  II, ISIZE, IHEAD, IUSE,
     $                  CP, RP, DN, N, WIR, WIC, WR, WC,
     $                  INTZER, INTZER, INTZER, INTZER, .FALSE.,
     $                  PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C          at this point, iuse = ineed and xuse = xneed
        ENDIF

        FFXP = XHEAD
        XHEAD = XHEAD + FFSIZE
        WXP = XHEAD
        XHEAD = XHEAD + FFDIMC
        XUSE = XUSE + FFSIZE + FFDIMC
        XNEED = XNEED + FFSIZE + FFDIMC
        INFO (20) = MAX (INFO (20), XUSE)
        INFO (21) = MAX (INFO (21), XNEED)
        IF (XHEAD .GT. XTAIL) THEN
C          error return, if not enough real memory:
           GO TO 9000
        ENDIF

C       ----------------------------------------------------------------
C       get memory usage for next call to MA38BD
C       ----------------------------------------------------------------

        XRUSE = XRUSE + FFSIZE
        XRMAX = MAX (XRMAX, XRUSE)

C       ----------------------------------------------------------------
C       zero the working array
C       ----------------------------------------------------------------

C       zero the contribution block:
        DO 380 J = 0, FFLEFR - 1
           DO 370 I = 0, FFLEFC - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
370        CONTINUE
380     CONTINUE

C       zero the pivot row:
        DO 390 J = 0, FFLEFR - 1
           XX (FFXP + J*FFDIMC + FFDIMC-1) = ZERO
390     CONTINUE

C       zero the pivot column:
        DO 400 I = 0, FFLEFC - 1
           XX (FFXP + (FFDIMR-1)*FFDIMC + I) = ZERO
400     CONTINUE

C       zero the pivot entry itself:
        XX (FFXP + (FFDIMR-1)*FFDIMC + FFDIMC-1) = ZERO

C       ----------------------------------------------------------------
C       current workspace usage:
C       ----------------------------------------------------------------

C       WpC (1..fflefc):        holds the pivot column pattern
C                               (excluding the pivot row index)
C       WpC (fflefc+1 .. n-npiv):       unused
C       WpC (n-npiv+1 .. n):            pivot columns in reverse order
C
C       WpR (1..fflefr):        holds the pivot row pattern
C                               (excluding the pivot column index)
C       WpR (fflefr+1 .. n-npiv):       unused
C       WpR (n-npiv+1 .. n):            pivot rows in reverse order
C
C       C (1..ffdimr, 1..ffdimc):  space for the new frontal matrix.
C
C       C (i,j) is located at XX (ffxp+((i)-1)+((j)-1)*ffdimc)
C
C       WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               Also, WiR (pivrow) is ffdimc-1, the offset in C of
C               the pivot row itself.
C               Otherwise, WiR (1..n) is -1
C
C       WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               Also, WiC (pivcol) is (ffdimr-1)*ffdimc,
C               the offset in C of the pivot column itself.
C               Otherwise, WiC (1..n) is -2 for nonpivotal columns,
C               and -1 for pivotal columns

C       ----------------------------------------------------------------
C       remove the columns affected by this element from degree lists
C       ----------------------------------------------------------------

        DO 410 J = 1, FFLEFR
           PC = CP (WPR (J))
           CDEG = II (PC+1)
           IF (CDEG .GT. 0) THEN
              CNXT = II (PC+7)
              CPRV = II (PC+8)
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = CPRV
              ENDIF
              IF (CPRV .NE. 0) THEN
                 II (CP (CPRV)+7) = CNXT
              ELSE
                 HEAD (CDEG) = CNXT
              ENDIF
           ENDIF
410     CONTINUE

C=======================================================================
C  Initialization of new frontal matrix is complete ]
C=======================================================================

C=======================================================================
C  Assemble and factorize the current frontal matrix [
C=======================================================================

C       for first pivot in frontal matrix, do all scans
        SCAN1 = 0
        SCAN2 = 0
        SCAN3 = 0
        SCAN4 = 0

        DO 1395 DUMMY3 = 1, N
C       (this loop is not indented due to its length)

C=======================================================================
C  Degree update and numerical assembly [
C=======================================================================

        KLEFT1 = KLEFT - 1

C       ----------------------------------------------------------------
C       SCAN1:  scan the element lists of each row in the pivot col
C               and compute the external column degree for each frontal
C       ----------------------------------------------------------------

        ROW = PIVROW
        DO 440 J = SCAN1, FFLEFC
           IF (J .NE. 0) THEN
C             Get a row;  otherwise, scan the pivot row if j is zero.
              ROW = WPC (J)
           ENDIF
           PR = RP (ROW)
           REP = (PR+2)
           RELN = WR (ROW)
CFPP$ NODEPCHK L
           DO 430 P = REP, REP + 2*RELN - 2, 2
              E = II (P)
              IF (WC (E) .LT. W0) THEN
C                this is the first time seen in either scan 1 or 2:
                 EP = RP (E)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 WR (E) = FLEFTR + W0
                 WC (E) = FLEFTC + W0
              ENDIF
              WC (E) = WC (E) - 1
430        CONTINUE
440     CONTINUE

C       ----------------------------------------------------------------
C       SCAN2:  scan the element lists of each col in the pivot row
C               and compute the external row degree for each frontal
C       ----------------------------------------------------------------

        COL = PIVCOL
        DO 460 J = SCAN2, FFLEFR
           IF (J .NE. 0) THEN
C             Get a col;  otherwise, scan the pivot col if j is zero.
              COL = WPR (J)
           ENDIF
           PC = CP (COL)
           CELN = II (PC+5)
           CEP = (PC+9)
CFPP$ NODEPCHK L
           DO 450 P = CEP, CEP + 2*CELN - 2, 2
              E = II (P)
              IF (WR (E) .LT. W0) THEN
C                this is the first time seen in either scan 1 or 2:
                 EP = RP (E)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 WR (E) = FLEFTR + W0
                 WC (E) = FLEFTC + W0
              ENDIF
              WR (E) = WR (E) - 1
450        CONTINUE
460     CONTINUE

C       ----------------------------------------------------------------
C       SCAN3:  scan the element lists of each column in pivot row
C               do degree update for the columns
C               assemble effective Usons and LU-sons
C       ----------------------------------------------------------------

C       flag Usons in Wc (e) as scanned (all now unflagged) [
C       uses Wc (e) for the link list.  Wc (e) <= 0
C       means that e is in the list, the external column
C       degree is zero, and -(Wc (e)) is the next element in
C       the Uson list.

        COL = PIVCOL
        DO 700 JJ = SCAN3, FFLEFR

C          -------------------------------------------------------------
C          assemble and update the degree of a column
C          -------------------------------------------------------------

           IF (JJ .NE. 0) THEN
C             Get a col;  otherwise, scan the pivot col if jj is zero
              COL = WPR (JJ)
           ENDIF

C          -------------------------------------------------------------
C          compute the degree, and partition the element list into
C          two parts.  The first part are not LUsons or Usons, and
C          are not assembled.  The second part is assembled.
C          -------------------------------------------------------------

           CDEG = 0
           DELN = 0
           PC = CP (COL)
           CEP = (PC+9)
           CELN = II (PC+5)
           IP2 = CEP + 2*CELN - 2
           XUDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
           DO 470 IP = CEP, IP2, 2
              E = II (IP)
              IF (WC (E) .GT. W0) THEN
C                this element cannot be assembled
                    CDEG = CDEG + (WC (E) - W0)
              ELSE
C                delete this tuple from the element list
                 DELN = DELN + 1
                 WM (DELN) = IP
              ENDIF
470       CONTINUE

          IF (DELN .NE. 0) THEN

C             ----------------------------------------------------------
C             move the deleted tuples to the end of the element list
C             ----------------------------------------------------------

              P2 = IP2
              DO 480 I = DELN, 1, -1
                 E = II (WM (I)  )
                 F = II (WM (I)+1)
                 II (WM (I)  ) = II (P2  )
                 II (WM (I)+1) = II (P2+1)
                 II (P2  ) = E
                 II (P2+1) = F
                 P2 = P2 - 2
480           CONTINUE

C             ----------------------------------------------------------
C             assemble from LUsons and Usons (the deleted tuples)
C             ----------------------------------------------------------

              DO 670 IP = P2 + 2, IP2, 2

C                -------------------------------------------------------
C                this is an LUson or Uson.  If fextc < 0 then this has
C                already been assembled.
C                -------------------------------------------------------

                 E = II (IP)
                 IF (WC (E) .LT. W0) THEN
C                   go to next iteration if already assembled
                    GOTO 670
                 ENDIF

C                -------------------------------------------------------
C                get scalar info, add son to list if not already there
C                -------------------------------------------------------

                 EP = RP (E)
                 FDIMC = II (EP+1)
                 FXP = II (EP+2)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 IF (E .LE. N) THEN
                    FLUIP = II (EP)
                    LUDEGR = II (FLUIP+2)
                    LUDEGC = II (FLUIP+3)
                    LUCP = (FLUIP + 7)
                    LURP = LUCP + LUDEGC
                    IF (WIR (E) .EQ. -1) THEN
                       WIR (E) = SONLST - N - 2
                       SONLST = E
                       NSONS = NSONS + 1
                    ENDIF
                 ELSE
C                   an artificial frontal matrix
                    LUDEGR = 1
                    LUDEGC = FDIMC
                    LUCP = (EP+9)
                    LURP = (EP+8)
                 ENDIF

C                -------------------------------------------------------
                 IF (WR (E) .EQ. W0) THEN
C                this is an LUson - assemble an entire frontal matrix
C                -------------------------------------------------------

C                   ----------------------------------------------------
                    IF (LUDEGC .EQ. FLEFTC) THEN
C                   no rows assembled out of this frontal yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets) [
                       DO 490 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          WM (I+1) = WIR (ROW2)
490                    CONTINUE

C                      -------------------------------------------------
                       IF (LUDEGR .EQ. FLEFTR) THEN
C                      no rows or cols assembled out of frontal yet
C                      -------------------------------------------------

                          DO 510 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 500 I = 0, LUDEGC-1
                                XX (XDP + WM (I+1)) =
     $                          XX (XDP + WM (I+1)) +
     $                          XX (FXP + J*FDIMC + I)
500                          CONTINUE
510                       CONTINUE

C                      -------------------------------------------------
                       ELSE
C                      only cols have been assembled out of frontal
C                      -------------------------------------------------

                          DO 530 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             IF (COL2 .GT. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 520 I = 0, LUDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
520                             CONTINUE
                             ENDIF
530                       CONTINUE
                       ENDIF
C                      done using Wm (1..ludegc for offsets) ]

C                   ----------------------------------------------------
                    ELSE
C                   some rows have been assembled out of this frontal
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets) [
                       DEGC = 0
                       DO 540 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                          ENDIF
540                    CONTINUE

C                      -------------------------------------------------
                       IF (LUDEGR .EQ. FLEFTR) THEN
C                      only rows assembled out of this frontal
C                      -------------------------------------------------

                          DO 560 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 550 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
550                          CONTINUE
560                       CONTINUE

C                      -------------------------------------------------
                       ELSE
C                      both rows and columns assembled out of frontal
C                      -------------------------------------------------

                          DO 580 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             IF (COL2 .GT. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 570 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
570                             CONTINUE
                             ENDIF
580                       CONTINUE
                       ENDIF
C                      done using Wm (1..ludegc for offsets) ]
                    ENDIF

C                   ----------------------------------------------------
C                   deallocate the LUson frontal matrix
C                   ----------------------------------------------------

                    WR (E) = -(NDN+2)
                    WC (E) = -(NDN+2)
                    IF (E .LE. N) THEN
                       RP (E) = FLUIP
                       II (EP) = FSCAL
                       INEED = INEED - FSCAL
                    ELSE
                       RP (E) = 0
                       II (EP) = FDIMC + CSCAL
                       INEED = INEED - (FDIMC + CSCAL)
                    ENDIF
                    II (EP+1) = -1
                    II (EP+6) = 0
                    MPREV = II (EP+4)
                    MNEXT = II (EP+3)
                    XNEED = XNEED - LUDEGR*LUDEGC
                    IF (MNEXT .NE. 0 .AND. II (MNEXT+6) .EQ. 0) THEN
C                      next block is free - delete it
                       MNEXT = II (MNEXT+3)
                       II (EP+3) = MNEXT
                       IF (MNEXT .NE. 0) THEN
                          II (MNEXT+4) = EP
                       ELSE
                          MTAIL = EP
                       ENDIF
                    ENDIF
                    IF (MPREV .NE. 0 .AND. II (MPREV+6) .EQ. 0) THEN
C                      previous block is free - delete it
                       II (EP+2) = II (MPREV+2)
                       MPREV = II (MPREV+4)
                       II (EP+4) = MPREV
                       IF (MPREV .NE. 0) THEN
                          II (MPREV+3) = EP
                       ELSE
                          MHEAD = EP
                       ENDIF
                    ENDIF
C                   get the size of the freed block
                    IF (MNEXT .NE. 0) THEN
                       XS = II (MNEXT+2) - II (EP+2)
                    ELSE
                       XS = FFXP - II (EP+2)
                    ENDIF
                    IF (XS .GT. XFREE) THEN
C                      keep track of the largest free block
                       XFREE = XS
                       PFREE = EP
                    ENDIF

C                   ----------------------------------------------------
C                   get memory usage for next call to MA38BD
C                   ----------------------------------------------------

                    XRUSE = XRUSE - LUDEGR*LUDEGC

C                -------------------------------------------------------
                 ELSE IF (WR (E) - W0 .LE. FLEFTR/2) THEN
C                this is a Uson - assemble all possible columns
C                -------------------------------------------------------

C                   ----------------------------------------------------
C                   add to Uson list - to be cleared just after scan 3
C                   ----------------------------------------------------

                    WC (E) = -USONS
                    USONS = E

C                   ----------------------------------------------------
                    IF (LUDEGC .EQ. FLEFTC) THEN
C                   no rows assembled out of this Uson frontal yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets)
                       DO 590 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          WM (I+1) = WIR (ROW2)
590                    CONTINUE

                       DO 610 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             IF (WIC (COL2) .GE. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 600 I = 0, LUDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
600                             CONTINUE
C                               flag this column as assembled from Uson
                                II (LURP+J) = -COL2
                             ENDIF
                          ENDIF
610                    CONTINUE

C                   ----------------------------------------------------
                    ELSE
C                   some rows already assembled out of this Uson frontal
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
C                      use Wm (1..ludegc for offsets)
                       DEGC = 0
                       DO 620 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                          ENDIF
620                    CONTINUE

                       DO 640 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             IF (WIC (COL2) .GE. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 630 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
630                             CONTINUE
C                               flag this column as assembled from Uson
                                II (LURP+J) = -COL2
                             ENDIF
                          ENDIF
640                    CONTINUE

                    ENDIF

                    FLEFTR = WR (E) - W0
                    II (EP+5) = FLEFTR

C                -------------------------------------------------------
                 ELSE
C                this is a Uson - assemble just one column
C                -------------------------------------------------------

C                   get the offset, f, from the (e,f) tuple
                    F = II (IP+1)

C                   ----------------------------------------------------
                    IF (LUDEGC .EQ. FLEFTC) THEN
C                   no rows assembled out of this Uson yet
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 650 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          XX (XUDP + WIR (ROW2)) =
     $                    XX (XUDP + WIR (ROW2)) +
     $                    XX (FXP + F*FDIMC + I)
650                    CONTINUE

C                   ----------------------------------------------------
                    ELSE
C                   some rows already assembled out of this Uson
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 660 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN
                             XX (XUDP + WIR (ROW2)) =
     $                       XX (XUDP + WIR (ROW2)) +
     $                       XX (FXP + F*FDIMC + I)
                          ENDIF
660                    CONTINUE
                    ENDIF

C                   ----------------------------------------------------
C                   decrement count of unassembled cols in frontal
C                   ----------------------------------------------------

                    II (EP+5) = FLEFTR - 1
C                   flag the column as assembled from the Uson
                    II (LURP+F) = -COL
                 ENDIF

670           CONTINUE

C             ----------------------------------------------------------
C             update the count of (e,f) tuples in the element list
C             ----------------------------------------------------------

              II (PC+5) = II (PC+5) - DELN
              INEED = INEED - 2*DELN
           ENDIF

C          -------------------------------------------------------------
C          assemble the original column and update count of entries
C          -------------------------------------------------------------

           CLEN = II (PC+6)
           IF (CLEN .GT. 0) THEN
              CSIZ = II (PC)
              IP = PC + CSIZ - CLEN
              DLEN = 0
CFPP$ NODEPCHK L
              DO 680 I = 0, CLEN - 1
                 ROW = II (IP+I)
                 IF (WIR (ROW) .GE. 0) THEN
C                   this entry can be assembled and deleted
                    DLEN = DLEN + 1
                    WM (DLEN) = I
                 ENDIF
680           CONTINUE
              IF (DLEN .NE. 0) THEN
                 CXP = II (PC+2)
                 DO 690 J = 1, DLEN
                    I = WM (J)
                    ROW = II (IP+I)
C                   assemble the entry
                    XX (XUDP + WIR (ROW)) =
     $              XX (XUDP + WIR (ROW)) + XX (CXP+I)
C                   and delete the entry
                    II (IP +I) = II (IP +J-1)
                    XX (CXP+I) = XX (CXP+J-1)
690              CONTINUE
                 CLEN = CLEN - DLEN
                 CXP = CXP + DLEN
                 INEED = INEED - DLEN
                 XNEED = XNEED - DLEN
                 II (PC+6) = CLEN
                 IF (CLEN .NE. 0) THEN
                    II (PC+2) = CXP
                 ELSE
C                   deallocate the real portion of the column:
                    MPREV = II (PC+4)
                    MNEXT = II (PC+3)
                    IF (MNEXT .NE. 0 .AND. II (MNEXT+6) .EQ. 0) THEN
C                      next block is free - delete it
                       MNEXT = II (MNEXT+3)
                       II (PC+3) = MNEXT
                       IF (MNEXT .NE. 0) THEN
                          II (MNEXT+4) = PC
                       ELSE
                          MTAIL = PC
                       ENDIF
                    ENDIF
                    IF (MPREV .NE. 0 .AND. II (MPREV+6) .EQ. 0) THEN
C                      previous block is free - delete it
                       II (PC+2) = II (MPREV+2)
                       MPREV = II (MPREV+4)
                       II (PC+4) = MPREV
                       IF (MPREV .NE. 0) THEN
                          II (MPREV+3) = PC
                       ELSE
                          MHEAD = PC
                       ENDIF
                    ENDIF
                    IF (PC .EQ. MHEAD) THEN
C                      adjust the start of the block if this is head
                       II (PC+2) = 1
                    ENDIF
C                   get the size of the freed block
                    IF (MNEXT .NE. 0) THEN
                       XS = II (MNEXT+2) - II (PC+2)
                    ELSE
                       XS = FFXP - II (PC+2)
                    ENDIF
                    IF (XS .GT. XFREE) THEN
C                      keep track of the largest free block
                       XFREE = XS
                       PFREE = PC
                    ENDIF
                 ENDIF
              ENDIF
              CDEG = CDEG + CLEN
           ENDIF

C          -------------------------------------------------------------
C          compute the upper bound degree - excluding current front
C          -------------------------------------------------------------

           CDEG2 = II (PC+1)
           CDEG = MIN (KLEFT1 - FFLEFC, CDEG2, CDEG)
           II (PC+1) = CDEG

700     CONTINUE

C       ----------------------------------------------------------------
C       SCAN-3 wrap-up:  remove flags from assembled Usons
C       ----------------------------------------------------------------

C       while (usons .ne. ndn+1) do
710     CONTINUE
        IF (USONS .NE. NDN+1) THEN
           NEXT = -WC (USONS)
           WC (USONS) = W0
           USONS = NEXT
C       end while:
        GOTO 710
        ENDIF
C       done un-flagging usons, all now unflagged in Wc (e) ]

C       ----------------------------------------------------------------
C       SCAN4:  scan element lists of each row in the pivot column
C               do degree update for the rows
C               assemble effective Lsons
C       ----------------------------------------------------------------

C       flag Lsons in Wr (e) (all are now unflagged) [
C       uses Wr (e) for the link list.  Wr (e) <= 0 means
C       that e is in the list, the external row degree is zero, and
C       -(Wr (e)) is the next element in the Lson list.

        ROW = PIVROW
        DO 840 JJ = SCAN4, FFLEFC

C          -------------------------------------------------------------
C          assemble and update the degree of a row
C          -------------------------------------------------------------

           IF (JJ .NE. 0) THEN
C             Get a row;  otherwise, scan the pivot row if jj is zero
              ROW = WPC (JJ)
           ENDIF

C          -------------------------------------------------------------
C          compute the degree, and partition the element list into
C          two parts.  The first part are not LUsons or Lsons, and
C          are not assembled.  The second part is assembled.
C          -------------------------------------------------------------

           RDEG = 0
           DELN = 0
           PR = RP (ROW)
           REP = (PR+2)
           RELN = WR (ROW)
           IP2 = REP + 2*RELN - 2
CFPP$ NODEPCHK L
           DO 720 IP = REP, IP2, 2
              E = II (IP)
              IF (WR (E) .GT. W0) THEN
                 RDEG = RDEG + (WR (E) - W0)
              ELSE
                 DELN = DELN + 1
                 WM (DELN) = IP
              ENDIF
720        CONTINUE

           IF (DELN .NE. 0) THEN

C             ----------------------------------------------------------
C             move the deleted tuples to the end of the element list
C             ----------------------------------------------------------

              P2 = IP2
              DO 730 I = DELN, 1, -1
                 E = II (WM (I)  )
                 F = II (WM (I)+1)
                 II (WM (I)  ) = II (P2  )
                 II (WM (I)+1) = II (P2+1)
                 II (P2  ) = E
                 II (P2+1) = F
                 P2 = P2 - 2
730           CONTINUE

C             ----------------------------------------------------------
C             assemble from Lsons (the deleted tuples)
C             ----------------------------------------------------------

              DO 810 IP = P2 + 2, IP2, 2

C                -------------------------------------------------------
C                this is an LUson or Lson.  If fextr < 0 then this has
C                already been assembled.  All LUsons have already been
C                assembled (in SCAN3, above).
C                -------------------------------------------------------

                 E = II (IP)
                 IF (WR (E) .LT. W0) THEN
C                   go to next iteration if already assembled
                    GOTO 810
                 ENDIF

C                -------------------------------------------------------
C                get scalar info, add to son list if not already there
C                -------------------------------------------------------

                 EP = RP (E)
                 FDIMC = II (EP+1)
                 FXP = II (EP+2)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 IF (E .LE. N) THEN
                    FLUIP = II (EP)
                    LUDEGR = II (FLUIP+2)
                    LUDEGC = II (FLUIP+3)
                    LUCP = (FLUIP + 7)
                    LURP = LUCP + LUDEGC
                    IF (WIR (E) .EQ. -1) THEN
                       WIR (E) = SONLST - N - 2
                       SONLST = E
                       NSONS = NSONS + 1
                    ENDIF
                 ELSE
C                   an artificial frontal matrix
                    LUDEGR = 1
                    LUDEGC = FDIMC
                    LUCP = (EP+9)
                    LURP = (EP+8)
                 ENDIF

C                -------------------------------------------------------
                 IF (WC (E) - W0 .LE. FLEFTC/2) THEN
C                this is an Lson - assemble all possible rows
C                -------------------------------------------------------

C                   ----------------------------------------------------
C                   add to Lson list - to be cleared just after scan 4
C                   ----------------------------------------------------

                    WR (E) = -LSONS
                    LSONS = E

C                   compute the compressed column offset vector
C                   use Wm (1..ludegc for offsets) [
                    DEGC = 0
                    DO 740 I = 0, LUDEGC-1
                       ROW2 = II (LUCP+I)
                       IF (ROW2 .GT. 0) THEN
                          IF (WIR (ROW2) .GE. 0) THEN
C                            this row will be assembled in loop below
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
C                            flag the row as assembled from the Lson
                             II (LUCP+I) = -ROW2
                          ENDIF
                       ENDIF
740                 CONTINUE

C                   ----------------------------------------------------
                    IF (LUDEGR .EQ. FLEFTR) THEN
C                   no columns assembled out this Lson yet
C                   ----------------------------------------------------

                       DO 760 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                          DO 750 I = 1, DEGC
                             XX (XDP + WM (I)) =
     $                       XX (XDP + WM (I)) +
     $                       XX (FXP + J*FDIMC + WJ (I))
750                       CONTINUE
760                    CONTINUE

C                   ----------------------------------------------------
                    ELSE
C                   some columns already assembled out of this Lson
C                   ----------------------------------------------------

                       DO 780 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 770 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
770                          CONTINUE
                          ENDIF
780                    CONTINUE
                    ENDIF

C                   done using Wm (1..ludegc for offsets) ]
                    FLEFTC = WC (E) - W0
                    II (EP+6) = FLEFTC

C                -------------------------------------------------------
                 ELSE
C                this is an Lson - assemble just one row
C                -------------------------------------------------------

                    XLDP = FFXP + WIR (ROW)
C                   get the offset, f, from the (e,f) tuple
                    F = II (IP+1)

C                   ----------------------------------------------------
                    IF (LUDEGR .EQ. FLEFTR) THEN
C                   no columns assembled out this Lson yet
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 790 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          XX (XLDP + WIC (COL2)) =
     $                    XX (XLDP + WIC (COL2)) +
     $                    XX (FXP + J*FDIMC + F)
790                    CONTINUE

C                   ----------------------------------------------------
                    ELSE
C                   some columns already assembled out of this Lson
C                   ----------------------------------------------------

CFPP$ NODEPCHK L
                       DO 800 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             XX (XLDP + WIC (COL2)) =
     $                       XX (XLDP + WIC (COL2)) +
     $                       XX (FXP + J*FDIMC + F)
                          ENDIF
800                    CONTINUE
                    ENDIF

                    II (EP+6) = FLEFTC - 1
C                   flag the row as assembled from the Lson
                    II (LUCP+F) = -ROW
                 ENDIF

810           CONTINUE

C             ----------------------------------------------------------
C             update the count of (e,f) tuples in the element list
C             ----------------------------------------------------------

              WR (ROW) = WR (ROW) - DELN
              INEED = INEED - 2*DELN
           ENDIF

C          -------------------------------------------------------------
C          assemble the original row and update count of entries
C          -------------------------------------------------------------

           RLEN = WC (ROW)
           IF (RLEN .GT. 0) THEN
C             do not scan a very long row:
              IF (RLEN .LE. RSCAN) THEN
                 RSIZ = II (PR)
                 IP = PR + RSIZ - RLEN
                 DLEN = 0
CFPP$ NODEPCHK L
                 DO 820 P = IP, IP + RLEN - 1
                    COL = II (P)
                    IF (WIC (COL) .NE. -2) THEN
C                      this entry can be assembled and deleted
C                      if WiC (col) = -1, it is an older pivot col,
C                      otherwise (>=0) it is in the current element
                       DLEN = DLEN + 1
                       WM (DLEN) = P
                    ENDIF
820              CONTINUE
                 IF (DLEN .NE. 0) THEN
                    DO 830 J = 1, DLEN
C                      delete the entry
                       II (WM (J)) = II (IP+J-1)
830                 CONTINUE
                    RLEN = RLEN - DLEN
                    INEED = INEED - DLEN
                    WC (ROW) = RLEN
                 ENDIF
              ENDIF
              RDEG = RDEG + RLEN
           ENDIF

C          -------------------------------------------------------------
C          compute the upper bound degree - excluding current front
C          -------------------------------------------------------------

           RDEG2 = II (PR+1)
           RDEG = MIN (KLEFT1 - FFLEFR, RDEG2, RDEG)
           II (PR+1) = RDEG

840     CONTINUE

C       ----------------------------------------------------------------
C       SCAN-4 wrap-up:  remove flags from assembled Lsons
C       ----------------------------------------------------------------

C       while (lsons .ne. ndn+1) do
850     CONTINUE
        IF (LSONS .NE. NDN+1) THEN
           NEXT = -WR (LSONS)
           WR (LSONS) = W0
           LSONS = NEXT
C       end while:
        GOTO 850
        ENDIF
C       done un-flagging Lsons, all now unflagged in Wr (e) ]

C=======================================================================
C  Degree update and numerical assemble is complete ]
C=======================================================================

C=======================================================================
C  Factorize frontal matrix until next pivot extends it [
C=======================================================================

        DO 1324 DUMMY4 = 1, N
C       (this loop is not indented due to its length)

C       ----------------------------------------------------------------
C       Wc (e) = fextc+w0, where fextc is the external column
C               degree for each element (ep = Rp (e)) appearing in
C               the element lists for each row in the pivot column.
C               if Wc (e) < w0, then fextc is defined as II (ep+6)
C
C       Wr (e) = fextr+w0, where fextr is the external row
C               degree for each element (ep = Rp (e)) appearing in
C               the element lists for each column in the pivot row
C               if Wr (e) < w0, then fextr is defined as II (ep+5)
C
C       WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               WiR (pivrow) is the offset of the latest pivot row
C
C       WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               WiC (pivcol) is the offset of the latest pivot column
C
C       WpR (1..fflefr) is the pivot row pattern (excl pivot cols)
C       WpC (1..fflefc) is the pivot col pattern (excl pivot rows)
C       ----------------------------------------------------------------

C=======================================================================
C  Divide pivot column by pivot
C=======================================================================

C       k-th pivot in frontal matrix located in C(ffdimc-k+1,ffdimr-k+1)
        XDP = FFXP + (FFDIMR - K) * FFDIMC
        X = XX (XDP + FFDIMC - K)

C       divide C(1:fflefc,ffdimr-k+1) by pivot value
        X = ONE / X
        DO 870 P = XDP, XDP + FFLEFC-1
           XX (P) = XX (P) * X
870     CONTINUE
C       count this as a call to the Level-1 BLAS:
        RINFO (4) = RINFO (4) + FFLEFC

C=======================================================================
C  A pivot step is complete
C=======================================================================

        KLEFT = KLEFT - 1
        NPIV = NPIV + 1
        INFO (17) = INFO (17) + 1

C       ----------------------------------------------------------------
C       the pivot column is fully assembled and scaled, and is now the
C       (npiv)-th column of L. The pivot row is the (npiv)-th row of U.
C       ----------------------------------------------------------------

        WPR (N-NPIV+1) = PIVROW
        WPC (N-NPIV+1) = PIVCOL
        WIR (PIVROW) = -1
        WIC (PIVCOL) = -1

C       ----------------------------------------------------------------
C       deallocate the pivot row and pivot column
C       ----------------------------------------------------------------

        RLEN = WC (PIVROW)
        INEED = INEED - CSCAL - RSCAL - RLEN
        PR = RP (PIVROW)
        PC = CP (PIVCOL)
        II (PR+1) = -1
        II (PC+1) = -1
        RP (PIVROW) = 0
        CP (PIVCOL) = 0

C=======================================================================
C  Local search for next pivot within current frontal matrix [
C=======================================================================

        FEDEGC = FFLEFC
        FEDEGR = FFLEFR
        PFOUND = .FALSE.
        OKCOL = FFLEFC .GT. 0
        OKROW = .FALSE.

C       ----------------------------------------------------------------
C       find column of minimum degree in current frontal row pattern
C       ----------------------------------------------------------------

C       among those columns of least (upper bound) degree, select the
C       column with the lowest column index
        IF (OKCOL) THEN
           COLPOS = 0
           PIVCOL = N+1
           CDEG = N+1
C          can this be vectorized?  This is the most intensive
C          non-vector loop.
           DO 880 J = 1, FFLEFR
              COL = WPR (J)
              PC = CP (COL)
              CDEG2 = II (PC+1)
              BETTER = CDEG2 .GE. 0 .AND.
     $                (CDEG2 .LT. CDEG .OR.
     $                (CDEG2 .EQ. CDEG .AND. COL .LT. PIVCOL))
              IF (BETTER) THEN
                 CDEG = CDEG2
                 COLPOS = J
                 PIVCOL = COL
              ENDIF
880        CONTINUE
           OKCOL = COLPOS .NE. 0
        ENDIF

C=======================================================================
C  Assemble candidate pivot column in temporary workspace
C=======================================================================

        IF (OKCOL) THEN
           PC = CP (PIVCOL)
           CLEN = II (PC+6)
           OKCOL = FEDEGC + CLEN .LE. FFDIMC
        ENDIF

        IF (OKCOL) THEN

C          -------------------------------------------------------------
C          copy candidate column from current frontal matrix into
C          work vector XX (wxp ... wxp+ffdimc-1) [
C          -------------------------------------------------------------

           P = FFXP + (COLPOS - 1) * FFDIMC - 1
CFPP$ NODEPCHK L
           DO 890 I = 1, FFLEFC
              XX (WXP-1+I) = XX (P+I)
890        CONTINUE

C          -------------------------------------------------------------
C          update candidate column with previous pivots in this front
C          -------------------------------------------------------------

           IF (K-K0 .GT. 0 .AND. FFLEFC .NE. 0) THEN
              CALL DGEMV ('N', FFLEFC, K-K0,
     $          -ONE, XX (FFXP + (FFDIMR - K) * FFDIMC)        ,FFDIMC,
     $                XX (FFXP + (COLPOS - 1) * FFDIMC + FFDIMC - K), 1,
     $           ONE, XX (WXP)                                      , 1)
              RINFO (3) = RINFO (3) + 2*FFLEFC*(K-K0)
           ENDIF

C          -------------------------------------------------------------
C          Compute extended pivot column in XX (wxp..wxp-1+fedegc).
C          Pattern of pivot column is placed in WpC (1..fedegc)
C          -------------------------------------------------------------

C          assemble the elements in the element list
           CEP = (PC+9)
           CELN = II (PC+5)
           DO 930 IP = CEP, CEP + 2*CELN - 2, 2
              E = II (IP)
              F = II (IP+1)
              EP = RP (E)
              FLEFTC = II (EP+6)
              FDIMC = II (EP+1)
              FXP = II (EP+2)
              IF (E .LE. N) THEN
                 FLUIP = II (EP)
                 LUCP = (FLUIP + 7)
                 LUDEGC = II (FLUIP+3)
              ELSE
                 LUCP = (EP+9)
                 LUDEGC = FDIMC
              ENDIF
              XP = FXP + F * FDIMC
C             split into 3 loops so that they all vectorize on a CRAY
              F1 = FEDEGC
              DO 900 P = LUCP, LUCP + LUDEGC - 1
                 ROW = II (P)
                 IF (ROW .GT. 0) THEN
                    IF (WIR (ROW) .LT. 0) THEN
                       F1 = F1 + 1
                       WPC (F1) = ROW
                    ENDIF
                 ENDIF
900           CONTINUE
              OKCOL = F1 + CLEN .LE. FFDIMC
              IF (.NOT. OKCOL) THEN
C                exit out of loop if column too long:
                 GO TO 940
              ENDIF
              DO 910 I = FEDEGC+1, F1
                 ROW = WPC (I)
                 WIR (ROW) = I - 1
                 XX (WXP-1+I) = ZERO
910           CONTINUE
              FEDEGC = F1
CFPP$ NODEPCHK L
              DO 920 J = 0, LUDEGC - 1
                 ROW = II (LUCP+J)
                 IF (ROW .GT. 0) THEN
                    XX (WXP + WIR (ROW)) =
     $              XX (WXP + WIR (ROW)) + XX (XP+J)
                 ENDIF
920           CONTINUE
930        CONTINUE
C          loop exit label:
940        CONTINUE
        ENDIF

C=======================================================================
C  Find candidate pivot row - unless candidate pivot column is too long
C=======================================================================

        IF (OKCOL) THEN

C          -------------------------------------------------------------
C          assemble the original entries in the column
C          -------------------------------------------------------------

           CSIZ = II (PC)
           IP = PC + CSIZ - CLEN
           CXP = II (PC+2)
CFPP$ NODEPCHK L
           DO 950 I = 0, CLEN - 1
              ROW = II (IP+I)
              WIR (ROW) = FEDEGC + I
              WPC (FEDEGC+1+I) = ROW
              XX  (WXP+FEDEGC+I) = XX (CXP+I)
950        CONTINUE
           FEDEGC = FEDEGC + CLEN

C          -------------------------------------------------------------
C          update degree of candidate column - excluding current front
C          -------------------------------------------------------------

           CDEG = FEDEGC - FFLEFC
           II (PC+1) = CDEG

C          -------------------------------------------------------------
C          find the maximum absolute value in the column
C          -------------------------------------------------------------

           MAXVAL = ABS (XX (WXP-1 + IDAMAX (FEDEGC, XX (WXP), 1)))
           RINFO (3) = RINFO (3) + FEDEGC
           TOLER = RELPT * MAXVAL
           RDEG = N+1

C          -------------------------------------------------------------
C          look for the best possible pivot row in this column
C          -------------------------------------------------------------

           IF (MAXVAL .GT. ZERO) THEN
              IF (SYMSRC) THEN
C                prefer symmetric pivots, if numerically acceptable
                 PIVROW = PIVCOL
                 ROWPOS = WIR (PIVROW) + 1
                 IF (ROWPOS .GT. 0 .AND. ROWPOS .LE. FFLEFC) THEN
C                   diagonal entry exists in the column pattern
C                   also within the current frontal matrix
                    X = ABS (XX (WXP-1+ROWPOS))
                    IF (X.GE.TOLER .AND. X.GT.ZERO) THEN
C                      diagonal entry is numerically acceptable
                       PR = RP (PIVROW)
                       RDEG = II (PR+1)
                    ENDIF
                 ENDIF
              ENDIF
              IF (RDEG .EQ. N+1) THEN
C                Continue searching - no diagonal found or sought for.
C                Minimize row degree subject to abs(value) constraints.
                 PIVROW = N+1
                 DO 960 I = 1, FFLEFC
                    ROW2 = WPC (I)
                    PR = RP (ROW2)
                    RDEG2 = II (PR+1)
C                   among those numerically acceptable rows of least
C                   (upper bound) degree, select the row with the
C                   lowest row index
                    BETTER = RDEG2 .LT. RDEG .OR.
     $                      (RDEG2 .EQ. RDEG .AND. ROW2 .LT. PIVROW)
                    IF (BETTER) THEN
                       X = ABS (XX (WXP-1+I))
                       IF (X.GE.TOLER .AND. X.GT.ZERO) THEN
                          PIVROW = ROW2
                          RDEG = RDEG2
                          ROWPOS = I
                       ENDIF
                    ENDIF
960              CONTINUE
              ENDIF
           ELSE
C             remove this column from any further pivot search
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ENDIF
           OKROW = RDEG .NE. N+1
        ENDIF

C       done using XX (wxp...wxp+ffdimc-1) ]

C=======================================================================
C  If found, construct candidate pivot row pattern
C=======================================================================

        IF (OKROW) THEN

C          -------------------------------------------------------------
C          assemble the elements in the element list
C          -------------------------------------------------------------

           PR = RP (PIVROW)
           REP = (PR+2)
           RELN = WR (PIVROW)
           DO 990 IP = REP, REP + 2*RELN - 2, 2
              E = II (IP)
              EP = RP (E)
              IF (E .LE. N) THEN
                 FLUIP = II (EP)
                 LUCP = (FLUIP + 7)
                 LUDEGR = II (FLUIP+2)
                 LUDEGC = II (FLUIP+3)
                 LURP = LUCP + LUDEGC
                 FLEFTR = II (EP+5)
                 OKROW = FLEFTR .LE. FFDIMR
                 IF (.NOT. OKROW) THEN
C                   exit out of loop if row too long:
                    GO TO 1000
                 ENDIF
C                split into two loops so that both vectorize on a CRAY
                 F1 = FEDEGR
                 DO 970 P = LURP, LURP + LUDEGR - 1
                    COL = II (P)
                    IF (COL .GT. 0) THEN
                       IF (WIC (COL) .EQ. -2) THEN
                          F1 = F1 + 1
                          WPR (F1) = COL
                       ENDIF
                    ENDIF
970              CONTINUE
                 OKROW = F1 .LE. FFDIMR
                 IF (.NOT. OKROW) THEN
C                   exit out of loop if row too long:
                    GO TO 1000
                 ENDIF
                 DO 980 I = FEDEGR+1, F1
                    WIC (WPR (I)) = (I - 1) * FFDIMC
980              CONTINUE
                 FEDEGR = F1
              ELSE
C                this is an artificial element (a dense column)
                 LURP = (EP+8)
                 COL = II (LURP)
                 IF (WIC (COL) .EQ. -2) THEN
                    WIC (COL) = FEDEGR * FFDIMC
                    FEDEGR = FEDEGR + 1
                    WPR (FEDEGR) = COL
                    OKROW = FEDEGR .LE. FFDIMR
                    IF (.NOT. OKROW) THEN
C                      exit out of loop if row too long:
                       GO TO 1000
                    ENDIF
                 ENDIF
              ENDIF
990        CONTINUE
C          loop exit label:
1000       CONTINUE
        ENDIF

        IF (OKROW) THEN

C          -------------------------------------------------------------
C          assemble the original entries in the row
C          -------------------------------------------------------------

           RLEN = WC (PIVROW)
           IF (RLEN .GT. 0) THEN
              F1 = FEDEGR
              RSIZ = II (PR)
              P2 = PR + RSIZ
C             split into two loops so that they both vectorize on a CRAY
              DO 1010 P = P2 - RLEN, P2 - 1
                 COL = II (P)
                 IF (WIC (COL) .EQ. -2) THEN
C                   this entry cannot be assembled, do not delete
                    F1 = F1 + 1
                    WPR (F1) = COL
                 ENDIF
1010          CONTINUE
              RLEN2 = F1 - FEDEGR
              IF (RLEN2 .LT. RLEN) THEN
C                delete one or more entries in the row
                 DO 1020 I = FEDEGR+1, F1
                    II (P2 - F1 + I - 1) = WPR (I)
1020             CONTINUE
                 INEED = INEED - (RLEN - RLEN2)
                 WC (PIVROW) = RLEN2
              ENDIF

C             ----------------------------------------------------------
C             update the candidate row degree - excluding current front
C             ----------------------------------------------------------

              RDEG = F1 - FFLEFR
              II (PR+1) = RDEG

C             ----------------------------------------------------------
C             pivot is found if candidate pivot row is not too long
C             ----------------------------------------------------------

              OKROW = F1 .LE. FFDIMR
              IF (OKROW) THEN
                 DO 1030 I = FEDEGR+1, F1
                    WIC (WPR (I)) = (I - 1) * FFDIMC
1030             CONTINUE
                 FEDEGR = F1
              ENDIF

           ELSE

C             ----------------------------------------------------------
C             update the candidate row degree - excluding current front
C             ----------------------------------------------------------

              RDEG = FEDEGR - FFLEFR
              II (PR+1) = RDEG
           ENDIF
        ENDIF

C       ----------------------------------------------------------------
C       if pivot not found: clear WiR and WiC
C       ----------------------------------------------------------------

        PFOUND = OKROW .AND. OKCOL
        IF (.NOT. PFOUND) THEN
           MOVELU = K .GT. 0
           DO 1040 I = FFLEFR+1, FEDEGR
              WIC (WPR (I)) = -2
1040       CONTINUE
           FEDEGR = FFLEFR
           DO 1050 I = FFLEFC+1, FEDEGC
              WIR (WPC (I)) = -1
1050       CONTINUE
           FEDEGC = FFLEFC
        ELSE
           MOVELU = FEDEGC .GT. FFDIMC - K .OR. FEDEGR .GT. FFDIMR - K
        ENDIF

C       ----------------------------------------------------------------
C       WpR (1..fflefr)                 unextended pivot row pattern
C       WpR (fflefr+1 .. fedegr)        extended pattern, if pfound
C       WpR (fedegr+1 .. n-npiv)        empty space
C       WpR (n-npiv+1 .. n)             pivot row order
C
C       WpC (1..fflefc)                 unextended pivot column pattern
C       WpC (fflefc+1 .. fedegc)        extended pattern, if pfound
C       WpC (fedegc+1 .. n-npiv)        empty space
C       WpC (n-npiv+1 .. n)             pivot column order
C       ----------------------------------------------------------------

C=======================================================================
C  Local pivot search complete ]
C=======================================================================

C=======================================================================
C  Update contribution block: rank-nb, or if LU arrowhead to be moved
C=======================================================================

        IF (K-K0 .GE. NB .OR. MOVELU) THEN
           CALL DGEMM ('N', 'N', FFLEFC, FFLEFR, K-K0,
     $          -ONE, XX (FFXP + (FFDIMR - K) * FFDIMC), FFDIMC,
     $                XX (FFXP +  FFDIMC - K)          , FFDIMC,
     $           ONE, XX (FFXP)                        , FFDIMC)
           RINFO (6) = RINFO (6) + 2*FFLEFC*FFLEFR*(K-K0)
           K0 = K
        ENDIF

C=======================================================================
C  Move the LU arrowhead if no pivot found, or pivot needs room
C=======================================================================

        IF (MOVELU) THEN

C          allocate permanent space for the LU arrowhead
           LUDEGR = FFLEFR
           LUDEGC = FFLEFC
           XS = K*LUDEGC + K*LUDEGR + K*K
           IS = 7 + LUDEGC + LUDEGR + NSONS
           IF (IS .GT. ITAIL-IHEAD .OR. XS .GT. XTAIL-XHEAD) THEN
              IF (IS .GT. ITAIL-IHEAD) THEN
C                garbage collection because we ran out of integer mem
                 INFO (14) = INFO (14) + 1
              ENDIF
              IF (XS .GT. XTAIL-XHEAD) THEN
C                garbage collection because we ran out of real mem
                 INFO (15) = INFO (15) + 1
              ENDIF
              CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                     II, ISIZE, IHEAD, IUSE,
     $                     CP, RP, DN, N, WIR, WIC, WR, WC,
     $                     FFXP, FFSIZE, WXP, FFDIMC, .FALSE.,
     $                     PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C             at this point, iuse = ineed and xuse = xneed
           ENDIF

           ITAIL = ITAIL - IS
           LUIP = ITAIL
           IUSE = IUSE + IS
           INEED = INEED + IS
           XTAIL = XTAIL - XS
           LUXP = XTAIL
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), INEED)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF

C          -------------------------------------------------------------
C          get memory usage for next call to MA38BD
C          -------------------------------------------------------------

           XRUSE = XRUSE + XS
           XRMAX = MAX (XRMAX, XRUSE)

C          -------------------------------------------------------------
C          save the new LU arrowhead
C          -------------------------------------------------------------

C          save the scalar data of the LU arrowhead
           II (LUIP) = LUXP
           II (LUIP+1) = K
           II (LUIP+2) = LUDEGR
           II (LUIP+3) = LUDEGC
           II (LUIP+4) = NSONS
           II (LUIP+5) = 0
           II (LUIP+6) = 0
           E = FFROW
           IF (E .EQ. E1) THEN
C             this is the first LU arrowhead from this global pivot
              LUIP1 = LUIP
           ENDIF
           WR (E) = -(NDN+2)
           WC (E) = -(NDN+2)

C          save column pattern
           LUCP = (LUIP + 7)
           DO 1060 I = 0, LUDEGC-1
              II (LUCP+I) = WPC (I+1)
1060       CONTINUE

C          save row pattern
           LURP = LUCP + LUDEGC
           DO 1070 I = 0, LUDEGR-1
              II (LURP+I) = WPR (I+1)
1070       CONTINUE

C          add list of sons after the end of the frontal matrix pattern
C          this list of sons is for the refactorization (MA38BD) only.
           LUSONP = LURP + LUDEGR
           IP = LUSONP
           E = SONLST
C          while (e > 0) do
1080       CONTINUE
           IF (E .GT. 0) THEN
              EP = RP (E)
              IF (WC (E) .EQ. -(NDN+2)) THEN
C                LUson
                 II (IP) = E
              ELSE IF (WC (E) .EQ. W0) THEN
C                Uson
                 II (IP) = E + N
              ELSE IF (WR (E) .EQ. W0) THEN
C                Lson
                 II (IP) = E + 2*N
              ENDIF
              NEXT = WIR (E) + N + 2
              WIR (E) = -1
              E = NEXT
              IP = IP + 1
C          end while:
           GOTO 1080
           ENDIF
           NSONS = 0
           SONLST = 0

C          move the L1,U1 matrix, compressing the dimension from
C          ffdimc to ldimc.  The LU arrowhead grows on top of stack.
           LDIMC = K + LUDEGC
           XP = FFXP + (FFDIMR-1)*FFDIMC + FFDIMC-1
           DO 1100 J = 0, K-1
CFPP$ NODEPCHK L
              DO 1090 I = 0, K-1
                 XX (LUXP + J*LDIMC + I) = XX (XP - J*FFDIMC - I)
1090          CONTINUE
1100       CONTINUE

C          move L2 matrix, compressing dimension from ffdimc to ludegc+k
           IF (LUDEGC .NE. 0) THEN
              LXP = LUXP + K
              XP = FFXP + (FFDIMR-1)*FFDIMC
              DO 1120 J = 0, K-1
CFPP$ NODEPCHK L
                 DO 1110 I = 0, LUDEGC-1
                    XX (LXP + J*LDIMC + I) = XX (XP - J*FFDIMC + I)
1110             CONTINUE
1120          CONTINUE
           ENDIF

C          move the U2 block.
           IF (LUDEGR .NE. 0) THEN
              UXP = LUXP + K * LDIMC
              XP = FFXP + FFDIMC-1
              DO 1140 J = 0, LUDEGR-1
CFPP$ NODEPCHK L
                 DO 1130 I = 0, K-1
                    XX (UXP + J*K + I) = XX (XP + J*FFDIMC - I)
1130             CONTINUE
1140          CONTINUE
           ENDIF

C          one more LU arrowhead has been created
           NLU = NLU + 1
           NZU = (K*(K-1)/2) + K*LUDEGC
           NZL = (K*(K-1)/2) + K*LUDEGR
           INFO (10) = INFO (10) + NZL
           INFO (11) = INFO (11) + NZU

C          no more rows of U or columns of L in current frontal array
           K = 0
           K0 = 0

           IF (PFOUND) THEN

C             ----------------------------------------------------------
C             Place the old frontal matrix as the only item in the son
C             list, since the next "implied" frontal matrix will have
C             this as its son.
C             ----------------------------------------------------------

              NSONS = 1
              E = FFROW
              WIR (E) = - N - 2
              SONLST = E

C             ----------------------------------------------------------
C             The contribution block of the old frontal matrix is still
C             stored in the current frontal matrix, and continues (in a
C             unifrontal sense) as a "new" frontal matrix (same array
C             but with a new name, and the LU arrowhead is removed and
C             placed in the LU factors).  Old name is "ffrow", new name
C             is "pivrow".
C             ----------------------------------------------------------

              RP (E) = LUIP
              FFROW = PIVROW
           ENDIF
        ENDIF

C=======================================================================
C  Stop the factorization of this frontal matrix if no pivot found
C=======================================================================

C       (this is the only way out of loop 1395)
        IF (.NOT. PFOUND) THEN
C          exit out of loop 1395 if pivot not found:
           GO TO 1400
        ENDIF

C=======================================================================
C  Update the pivot column, and move into position as (k+1)-st col of L
C=======================================================================

        XSP = (COLPOS - 1) * FFDIMC
        XDP = (FFDIMR - K - 1) * FFDIMC
        FSP = FFXP + XSP
        FDP = FFXP + XDP

        IF (K-K0 .GT. 0 .AND. FFLEFC .NE. 0) THEN
           CALL DGEMV ('N', FFLEFC, K-K0,
     $          -ONE, XX (FDP + FFDIMC    ), FFDIMC,
     $                XX (FSP + FFDIMC - K), 1,
     $           ONE, XX (FSP             ), 1)
           RINFO (5) = RINFO (5) + 2*FFLEFC*(K-K0)
        ENDIF

        IF (FFLEFR .LT. FFDIMR - K) THEN

           XLP = (FFLEFR - 1) * FFDIMC
           IF (FFLEFR .EQ. COLPOS) THEN

C             ----------------------------------------------------------
C             move C(:,colpos) => C(:,ffdimr-k)
C             ----------------------------------------------------------

C             copy only what needs to be copied
C             column of the contribution block:
CFPP$ NODEPCHK L
              DO 1160 I = 0, FFLEFC - 1
                 XX (FDP+I) = XX (FSP+I)
1160          CONTINUE
C             column of the U2 block
CFPP$ NODEPCHK L
              DO 1170 I = FFDIMC - K, FFDIMC - 1
                 XX (FDP+I) = XX (FSP+I)
1170          CONTINUE

           ELSE

C             ----------------------------------------------------------
C             move C(:,colpos) => C(:,ffdimr-k)
C             move C(:,fflefr) => C(:,colpos)
C             ----------------------------------------------------------

              FLP = FFXP + XLP

C             copy only what needs to be copied
C             columns of the contribution block:
CFPP$ NODEPCHK L
              DO 1190 I = 0, FFLEFC - 1
                 XX (FDP+I) = XX (FSP+I)
                 XX (FSP+I) = XX (FLP+I)
1190          CONTINUE
C             columns of the U2 block:
CFPP$ NODEPCHK L
              DO 1200 I = FFDIMC - K, FFDIMC - 1
                 XX (FDP+I) = XX (FSP+I)
                 XX (FSP+I) = XX (FLP+I)
1200          CONTINUE

              SWPCOL = WPR (FFLEFR)
              WPR (COLPOS) = SWPCOL
              WIC (SWPCOL) = XSP
           ENDIF

           IF (FEDEGR .NE. FFLEFR) THEN
C             move column fedegr to column fflefr (pattern only)
              SWPCOL = WPR (FEDEGR)
              WPR (FFLEFR) = SWPCOL
              WIC (SWPCOL) = XLP
           ENDIF

        ELSE IF (COLPOS .NE. FFDIMR - K) THEN

C          -------------------------------------------------------------
C          swap C(:,colpos) <=> C (:,ffdimr-k)
C          -------------------------------------------------------------

C          swap only what needs to be swapped
C          columns of the contribution block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1220 I = 0, FFLEFC - 1
              X = XX (FDP+I)
              XX (FDP+I) = XX (FSP+I)
              XX (FSP+I) = X
1220       CONTINUE
C          columns of the U2 block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1230 I = FFDIMC - K, FFDIMC - 1
              X = XX (FDP+I)
              XX (FDP+I) = XX (FSP+I)
              XX (FSP+I) = X
1230       CONTINUE
           SWPCOL = WPR (FFDIMR - K)
           WPR (COLPOS) = SWPCOL
           WIC (SWPCOL) = XSP
        ENDIF

        WIC (PIVCOL) = XDP
        FEDEGR = FEDEGR - 1
        SCAN2 = FFLEFR
        FFLEFR = FFLEFR - 1

C=======================================================================
C  Move pivot row into position as (k+1)-st row of U, and update
C=======================================================================

        XSP = ROWPOS - 1
        XDP = FFDIMC - K - 1
        FSP = FFXP + XSP
        FDP = FFXP + XDP

        IF (FFLEFC .LT. FFDIMC - K) THEN

           XLP = FFLEFC - 1
           IF (FFLEFC .EQ. ROWPOS) THEN

C             ----------------------------------------------------------
C             move C(rowpos,:) => C(ffdimc-k,:)
C             ----------------------------------------------------------

C             copy only what needs to be copied
C             row of the contribution block:
CFPP$ NODEPCHK L
              DO 1250 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
1250          CONTINUE
C             row of the L2 block:
CFPP$ NODEPCHK L
              DO 1260 J = (FFDIMR - K - 1) * FFDIMC,
     $                    (FFDIMR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
1260          CONTINUE

           ELSE

C             ----------------------------------------------------------
C             move C(rowpos,:) => C(ffdimc-k,:)
C             move C(fflefc,:) => C(rowpos,:)
C             ----------------------------------------------------------

              FLP = FFXP + XLP
C             copy only what needs to be copied
C             rows of the contribution block:
CFPP$ NODEPCHK L
              DO 1280 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = XX (FLP+J)
1280          CONTINUE
C             rows of the L2 block:
CFPP$ NODEPCHK L
              DO 1290 J = (FFDIMR - K - 1) * FFDIMC,
     $                    (FFDIMR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = XX (FLP+J)
1290          CONTINUE
              SWPROW = WPC (FFLEFC)
              WPC (ROWPOS) = SWPROW
              WIR (SWPROW) = XSP
           ENDIF

           IF (FEDEGC .NE. FFLEFC) THEN
C             move row fedegc to row fflefc (pattern only)
              SWPROW = WPC (FEDEGC)
              WPC (FFLEFC) = SWPROW
              WIR (SWPROW) = XLP
           ENDIF

        ELSE IF (ROWPOS .NE. FFDIMC - K) THEN

C          -------------------------------------------------------------
C          swap C(rowpos,:) <=> C (ffdimc-k,:)
C          -------------------------------------------------------------

C          swap only what needs to be swapped
C          rows of the contribution block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1310 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC
              X = XX (FDP+J)
              XX (FDP+J) = XX (FSP+J)
              XX (FSP+J) = X
1310       CONTINUE
C          rows of the L2 block:
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1320 J = (FFDIMR - K - 1) * FFDIMC,
     $                 (FFDIMR - 1) * FFDIMC, FFDIMC
              X = XX (FDP+J)
           XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = X

1320       CONTINUE
           SWPROW = WPC (FFDIMC - K)
           WPC (ROWPOS) = SWPROW
           WIR (SWPROW) = XSP
        ENDIF

        WIR (PIVROW) = XDP
        FEDEGC = FEDEGC - 1
        SCAN1 = FFLEFC
        FFLEFC = FFLEFC - 1

        IF (K-K0 .GT. 0 .AND. FFLEFR .GT. 0) THEN
           CALL DGEMV ('T', K-K0, FFLEFR,
     $       -ONE, XX (FDP + 1)                    , FFDIMC,
     $             XX (FDP + (FFDIMR - K) * FFDIMC), FFDIMC,
     $        ONE, XX (FDP)                        , FFDIMC)
           RINFO (5) = RINFO (5) + 2*(K-K0)*FFLEFR
        ENDIF

C=======================================================================
C  Prepare for degree update and next local pivot search
C=======================================================================

C       ----------------------------------------------------------------
C       if only column pattern has been extended:
C               scan1:  new rows only
C               scan2:  no columns scanned
C               scan3:  all columns
C               scan4:  new rows only
C
C       if only row pattern has been extended:
C               scan1:  no rows scanned
C               scan2:  new columns only
C               scan3:  new columns only
C               scan4:  all rows
C
C       if both row and column pattern have been extended:
C               scan1:  new rows only
C               scan2:  new columns only
C               scan3:  all columns
C               scan4:  all rows
C
C       if no patterns have been extended:
C               scan1-4: none
C       ----------------------------------------------------------------

        IF (FEDEGC .EQ. FFLEFC) THEN
C          column pattern has not been extended
           SCAN3 = FFLEFR + 1
        ELSE
C          column pattern has been extended.
           SCAN3 = 0
        ENDIF

        IF (FEDEGR .EQ. FFLEFR) THEN
C          row pattern has not been extended
           SCAN4 = FFLEFC + 1
        ELSE
C          row pattern has been extended
           SCAN4 = 0
        ENDIF

C=======================================================================
C  Finished with step k (except for assembly and scaling of pivot col)
C=======================================================================

        K = K + 1

C       ----------------------------------------------------------------
C       exit loop if frontal matrix has been extended
C       ----------------------------------------------------------------

        IF (FEDEGR .NE. FFLEFR .OR. FEDEGC .NE. FFLEFC) THEN
           GO TO 1325
        ENDIF

1324    CONTINUE
C       exit label for loop 1324:
1325    CONTINUE

C=======================================================================
C  Finished factorizing while frontal matrix is not extended ]
C=======================================================================

C=======================================================================
C  Extend the frontal matrix [
C=======================================================================

C       ----------------------------------------------------------------
C       Zero the newly extended frontal matrix
C       ----------------------------------------------------------------

C       fill-in due to amalgamation caused by this step is
C       k*(fedegr-fflefr+fedegc-fflefc)

        DO 1350 J = FFLEFR, FEDEGR - 1
C          zero the new columns in the contribution block:
           DO 1330 I = 0, FEDEGC - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1330       CONTINUE
C          zero the new columns in U block:
           DO 1340 I = FFDIMC - K, FFDIMC - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1340       CONTINUE
1350    CONTINUE

CFPP$ NODEPCHK L
        DO 1380 I = FFLEFC, FEDEGC - 1
C          zero the new rows in the contribution block:
CFPP$ NODEPCHK L
           DO 1360 J = 0, FFLEFR - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1360       CONTINUE
C          zero the new rows in L block:
CFPP$ NODEPCHK L
           DO 1370 J = FFDIMR - K, FFDIMR - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1370       CONTINUE
1380    CONTINUE

C       ----------------------------------------------------------------
C       remove the new columns from the degree lists
C       ----------------------------------------------------------------

        DO 1390 J = FFLEFR+1, FEDEGR
           PC = CP (WPR (J))
           CDEG = II (PC+1)
           IF (CDEG .GT. 0) THEN
              CNXT = II (PC+7)
              CPRV = II (PC+8)
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = CPRV
              ENDIF
              IF (CPRV .NE. 0) THEN
                 II (CP (CPRV)+7) = CNXT
              ELSE
                 HEAD (CDEG) = CNXT
              ENDIF
           ENDIF
1390    CONTINUE

C       ----------------------------------------------------------------
C       finalize extended row and column pattern of the frontal matrix
C       ----------------------------------------------------------------

        FFLEFC = FEDEGC
        FFLEFR = FEDEGR
        FMAXR = MAX (FMAXR, FFLEFR + K)
        FMAXC = MAX (FMAXC, FFLEFC + K)

C=======================================================================
C  Done extending the current frontal matrix ]
C=======================================================================

1395    CONTINUE
C       exit label for loop 1395:
1400    CONTINUE

C=======================================================================
C  Done assembling and factorizing the current frontal matrix ]
C=======================================================================

C=======================================================================
C  Wrap-up:  complete the current frontal matrix [
C=======================================================================

C       ----------------------------------------------------------------
C       store the maximum front size in the first LU arrowhead
C       ----------------------------------------------------------------

        II (LUIP1+5) = FMAXR
        II (LUIP1+6) = FMAXC

C       one more frontal matrix is finished
        INFO (13) = INFO (13) + 1

C       ----------------------------------------------------------------
C       add the current frontal matrix to the degrees of each column,
C       and place the modified columns back in the degree lists
C       ----------------------------------------------------------------

C       do so in reverse order to try to improve pivot tie-breaking
        DO 1410 J = FFLEFR, 1, -1
           COL = WPR (J)
           PC = CP (COL)
C          add the current frontal matrix to the degree
           CDEG = II (PC+1)
           CDEG = MIN (KLEFT, CDEG + FFLEFC)
           IF (CDEG .GT. 0) THEN
              II (PC+1) = CDEG
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = COL
              ENDIF
              HEAD (CDEG) = COL
              MINDEG = MIN (MINDEG, CDEG)
           ENDIF
1410    CONTINUE

C       ----------------------------------------------------------------
C       add the current frontal matrix to the degrees of each row
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 1420 I = 1, FFLEFC
           ROW = WPC (I)
           PR = RP (ROW)
           RDEG = II (PR+1)
           RDEG = MIN (KLEFT, RDEG + FFLEFR)
           II (PR+1) = RDEG
1420    CONTINUE

C       ----------------------------------------------------------------
C       Reset w0 so that Wr (1..n) < w0 and Wc (1..n) < w0.
C       Also ensure that w0 + n would not cause integer overflow
C       ----------------------------------------------------------------

        W1 = W0 + FMAX + 1
        IF (W1 .LE. W0) THEN
           W0 = NDN+2
           DO 1430 E = 1, N+DN
              IF (WR (E) .GT. NDN) THEN
C                this is a frontal matrix
                 WR (E) = W0-1
                 WC (E) = W0-1
              ENDIF
1430       CONTINUE
        ELSE
           W0 = W1
        ENDIF

C       ----------------------------------------------------------------
C       deallocate work vector
C       ----------------------------------------------------------------

        XUSE = XUSE - FFDIMC
        XNEED = XNEED - FFDIMC
        XHEAD = XHEAD - FFDIMC

C       ----------------------------------------------------------------
C       get the name of this new frontal matrix, and size of
C       contribution block
C       ----------------------------------------------------------------

        E = FFROW
        XS = FFLEFR * FFLEFC
        FMAX = MAX (FMAX, FFLEFR, FFLEFC)

C       ----------------------------------------------------------------
C       get memory usage for next call to MA38BD
C       ----------------------------------------------------------------

        XRUSE = XRUSE - FFSIZE + XS

C       ----------------------------------------------------------------
C       if contribution block empty, deallocate and continue next step
C       ----------------------------------------------------------------

        IF (FFLEFR .LE. 0 .OR. FFLEFC .LE. 0) THEN
           RP (E) = LUIP
           XUSE = XUSE - FFSIZE
           XNEED = XNEED - FFSIZE
           XHEAD = FFXP
           DO 1440 I = 1, FFLEFR
              WIC (WPR (I)) = -2
1440       CONTINUE
           DO 1450 I = 1, FFLEFC
              WIR (WPC (I)) = -1
1450       CONTINUE
C          next iteration of main factorization loop 1540:
           GOTO 1540
        ENDIF

C       ----------------------------------------------------------------
C       prepare the contribution block for later assembly
C       ----------------------------------------------------------------

        IF (FSCAL .GT. ITAIL-IHEAD) THEN
           INFO (14) = INFO (14) + 1
           INTZER = 0
           CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                  II, ISIZE, IHEAD, IUSE,
     $                  CP, RP, DN, N, WIR, WIC, WR, WC,
     $                  FFXP, FFSIZE, INTZER, INTZER, .FALSE.,
     $                  PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C          at this point, iuse = ineed and xuse = xneed
        ENDIF

        EP = IHEAD
        IHEAD = IHEAD + FSCAL
        IUSE = IUSE + FSCAL
        INEED = INEED + FSCAL
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .GT. ITAIL) THEN
C          error return, if not enough integer memory:
C          (highly unlikely to run out of memory at this point)
           GO TO 9000
        ENDIF

        RP (E) = EP
        II (EP) = LUIP
        II (EP+5) = FFLEFR
        II (EP+6) = FFLEFC
        WR (E) = W0-1
        WC (E) = W0-1

C       count the numerical assembly
        RINFO (2) = RINFO (2) + XS

        IF (XS .LE. XFREE) THEN

C          -------------------------------------------------------------
C          compress and store the contribution block in a freed block
C          -------------------------------------------------------------

C          place the new block in the list in front of the free block
           XDP = II (PFREE+2)
           II (PFREE+2) = II (PFREE+2) + XS
           XFREE = XFREE - XS
           MPREV = II (PFREE+4)
           IF (XFREE .EQ. 0) THEN
C             delete the free block if its size is zero
              MNEXT = II (PFREE+3)
              PFREE = 0
              XFREE = -1
           ELSE
              MNEXT = PFREE
           ENDIF
           IF (MNEXT .NE. 0) THEN
              II (MNEXT+4) = EP
           ELSE
              MTAIL = EP
           ENDIF
           IF (MPREV .NE. 0) THEN
              II (MPREV+3) = EP
           ELSE
              MHEAD = EP
           ENDIF
           DO 1470 J = 0, FFLEFR - 1
CFPP$ NODEPCHK L
              DO 1460 I = 0, FFLEFC - 1
                 XX (XDP + J*FFLEFC + I) = XX (FFXP + J*FFDIMC + I)
1460          CONTINUE
1470       CONTINUE
           XHEAD = FFXP
           XUSE = XUSE - FFSIZE
           XNEED = XNEED - FFSIZE + XS
           FFDIMC = FFLEFC
           II (EP+1) = FFDIMC
           II (EP+2) = XDP
           II (EP+3) = MNEXT
           II (EP+4) = MPREV

        ELSE

C          -------------------------------------------------------------
C          deallocate part of the unused portion of the frontal matrix
C          -------------------------------------------------------------

C          leave the contribution block C (1..fflefc, 1..fflefr) at the
C          head of XX, with column dimension of ffdimc and in space
C          of size (fflefr-1)*ffdimc for the first fflefr columns, and
C          fflefc for the last column.
           XNEED = XNEED - FFSIZE + XS
           XS = FFSIZE - (FFLEFC + (FFLEFR-1)*FFDIMC)
           XHEAD = XHEAD - XS
           XUSE = XUSE - XS
           II (EP+1) = FFDIMC
           II (EP+2) = FFXP
           II (EP+3) = 0
           II (EP+4) = MTAIL
           IF (MTAIL .EQ. 0) THEN
              MHEAD = EP
           ELSE
              II (MTAIL+3) = EP
           ENDIF
           MTAIL = EP
        ENDIF

C       ----------------------------------------------------------------
C       add tuples to the amount of integer space needed - and add
C       limit+cscal to maximum need to account for worst-case possible
C       reallocation of rows/columns.  Required integer memory usage
C       is guaranteed not to exceed iworst during the placement of (e,f)
C       tuples in the two loops below.
C       ----------------------------------------------------------------

        INEED = INEED + 2*(FFLEFR+FFLEFC)
        IWORST = INEED + LIMIT + CSCAL
        INFO (19) = MAX (INFO (19), IWORST)
        INFO (18) = MAX (INFO (18), IWORST)

C       ----------------------------------------------------------------
C       place (e,f) in the element list of each column
C       ----------------------------------------------------------------

        DO 1500 I = 1, FFLEFR
           COL = WPR (I)
           PC = CP (COL)
           CELN = II (PC+5)
           CSIZ = II (PC)
           CLEN = II (PC+6)
C          clear the column offset
           WIC (COL) = -2

C          -------------------------------------------------------------
C          make sure an empty slot exists - if not, create one
C          -------------------------------------------------------------

           IF (2*(CELN+1) + CLEN + CSCAL .GT. CSIZ) THEN

C             ----------------------------------------------------------
C             no room exists - reallocate elsewhere
C             ----------------------------------------------------------

C             at least this much space is needed:
              IS = 2 * (CELN + 1) + CLEN
C             add some slots for growth: at least 8 tuples,
C             or double the size - whichever is larger (but with a total
C             size not larger than limit+cscal)
              IS = MIN (IS + MAX (16, IS), LIMIT)
              CSIZ2 = IS + CSCAL

C             ----------------------------------------------------------
C             make sure enough room exists: garbage collection if needed
C             ----------------------------------------------------------

              IF (CSIZ2 .GT. ITAIL-IHEAD) THEN
C                garbage collection:
                 INFO (14) = INFO (14) + 1
                 INTZER = 0
                 CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                        II, ISIZE, IHEAD, IUSE,
     $                        CP, RP, DN, N, WIR, WIC, WR, WC,
     $                        INTZER, INTZER, INTZER, INTZER, .TRUE.,
     $                        PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C                at this point, iuse+csiz2 <= iworst and xuse = xneed
                 PC = CP (COL)
                 CSIZ = II (PC)
              ENDIF

C             ----------------------------------------------------------
C             get space for the new copy
C             ----------------------------------------------------------

              PC2 = IHEAD
              IHEAD = IHEAD + CSIZ2
              IUSE = IUSE + CSIZ2
              INFO (18) = MAX (INFO (18), IUSE)
              IF (IHEAD .GT. ITAIL) THEN
C                error return, if not enough integer memory:
                 GO TO 9000
              ENDIF

C             ----------------------------------------------------------
C             make the copy, leaving hole in middle for element list
C             ----------------------------------------------------------

C             copy the cscal scalars, and the element list
CFPP$ NODEPCHK L
              DO 1480 J = 0, CSCAL + 2*CELN - 1
                 II (PC2+J) = II (PC+J)
1480          CONTINUE

C             copy column indices of original entries (XX is unchanged)
CFPP$ NODEPCHK L
              DO 1490 J = 0, CLEN - 1
                 II (PC2+CSIZ2-CLEN+J) = II (PC+CSIZ-CLEN+J)
1490          CONTINUE

              IF (CLEN .GT. 0) THEN
C                place the new block in the memory-list
                 MNEXT = II (PC2+3)
                 MPREV = II (PC2+4)
                 IF (MNEXT .NE. 0) THEN
                    II (MNEXT+4) = PC2
                 ELSE
                    MTAIL = PC2
                 ENDIF
                 IF (MPREV .NE. 0) THEN
                    II (MPREV+3) = PC2
                 ELSE
                    MHEAD = PC2
                 ENDIF
              ENDIF

              CP (COL) = PC2
              II (PC2) = CSIZ2

C             ----------------------------------------------------------
C             deallocate the old copy of the column in II (not in XX)
C             ----------------------------------------------------------

              II (PC+1) = -1
              II (PC+6) = 0
              PC = PC2
           ENDIF

C          -------------------------------------------------------------
C          place the new (e,f) tuple in the element list of the column
C          -------------------------------------------------------------

           CEP = (PC+9)
           II (CEP + 2*CELN  ) = E
           II (CEP + 2*CELN+1) = I - 1
           II (PC+5) = CELN + 1
1500    CONTINUE

C       ----------------------------------------------------------------
C       place (e,f) in the element list of each row
C       ----------------------------------------------------------------

        DO 1530 I = 1, FFLEFC
           ROW = WPC (I)
           PR = RP (ROW)
           RSIZ = II (PR)
           RELN = WR (ROW)
           RLEN = WC (ROW)
C          clear the row offset
           WIR (ROW) = -1

C          -------------------------------------------------------------
C          make sure an empty slot exists - if not, create one
C          -------------------------------------------------------------

           IF (2*(RELN+1) + RLEN + RSCAL .GT. RSIZ) THEN

C             ----------------------------------------------------------
C             no room exists - reallocate elsewhere
C             ----------------------------------------------------------

C             at least this much space is needed:
              IS = 2 * (RELN + 1) + RLEN
C             add some extra slots for growth - for at least 8
C             tuples, or double the size (but with a total size not
C             larger than limit+rscal)
              IS = MIN (IS + MAX (16, IS), LIMIT)
              RSIZ2 = IS + RSCAL

C             ----------------------------------------------------------
C             make sure enough room exists: garbage collection if needed
C             ----------------------------------------------------------

              IF (RSIZ2 .GT. ITAIL-IHEAD) THEN
C                garbage collection:
                 INFO (14) = INFO (14) + 1
                 INTZER = 0
                 CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                        II, ISIZE, IHEAD, IUSE,
     $                        CP, RP, DN, N, WIR, WIC, WR, WC,
     $                        INTZER, INTZER, INTZER, INTZER, .TRUE.,
     $                        PFREE, XFREE, MHEAD, MTAIL, SLOTS)
C                at this point, iuse+rsiz2 <= iworst and xuse = xneed
                 PR = RP (ROW)
                 RSIZ = II (PR)
              ENDIF

C             ----------------------------------------------------------
C             get space for the new copy
C             ----------------------------------------------------------

              PR2 = IHEAD
              IHEAD = IHEAD + RSIZ2
              IUSE = IUSE + RSIZ2
              INFO (18) = MAX (INFO (18), IUSE)
              IF (IHEAD .GT. ITAIL) THEN
C                error return, if not enough integer memory:
                 GO TO 9000
              ENDIF

C             ----------------------------------------------------------
C             make the copy, leaving hole in middle for element list
C             ----------------------------------------------------------

C             copy the rscal scalars, and the element list
CFPP$ NODEPCHK L
              DO 1510 J = 0, RSCAL + 2*RELN - 1
                 II (PR2+J) = II (PR+J)
1510          CONTINUE

C             copy the original entries
CFPP$ NODEPCHK L
              DO 1520 J = 0, RLEN - 1
                 II (PR2+RSIZ2-RLEN+J) = II (PR+RSIZ-RLEN+J)
1520          CONTINUE

              RP (ROW) = PR2
              II (PR2) = RSIZ2

C             ----------------------------------------------------------
C             deallocate the old copy of the row
C             ----------------------------------------------------------

              II (PR+1) = -1
              PR = PR2
           ENDIF

C          -------------------------------------------------------------
C          place the new (e,f) tuple in the element list of the row
C          -------------------------------------------------------------

           REP = (PR+2)
           II (REP + 2*RELN  ) = E
           II (REP + 2*RELN+1) = I - 1
           WR (ROW) = RELN + 1
1530    CONTINUE

C=======================================================================
C  Wrap-up of factorized frontal matrix is complete ]
C=======================================================================

1540    CONTINUE
C       exit label for loop 1540:
2000    CONTINUE

C=======================================================================
C=======================================================================
C  END OF MAIN FACTORIZATION LOOP ]
C=======================================================================
C=======================================================================

C=======================================================================
C  Wrap-up:  store LU factors in their final form [
C=======================================================================

C       ----------------------------------------------------------------
C       deallocate all remaining columns, rows, and frontal matrices
C       ----------------------------------------------------------------

        IUSE = IUSE - (IHEAD - 1)
        XUSE = XUSE - (XHEAD - 1)
        INEED = IUSE
        XNEED = XUSE
        IHEAD = 1
        XHEAD = 1

        IF (NLU .EQ. 0) THEN
C          LU factors are completely empty (A = 0).
C          Add one integer and one real, to simplify rest of code.
C          Otherwise, some arrays in MA38BD or MA38CD would have
C          zero size, which can cause an address fault.
           ITAIL = ISIZE
           XTAIL = XSIZE
           IUSE = IUSE + 1
           XUSE = XUSE + 1
           INEED = IUSE
           XNEED = XUSE
           IP = ITAIL
           XP = XTAIL
        ENDIF

C       ----------------------------------------------------------------
C       compute permutation and inverse permutation vectors.
C       use WiR/C for the row/col permutation, and WpR/C for the
C       inverse row/col permutation.
C       ----------------------------------------------------------------

        DO 2010 K = 1, N
C          the kth pivot row and column:
           ROW = WPR (N-K+1)
           COL = WPC (N-K+1)
           WIR (K) = ROW
           WIC (K) = COL
2010    CONTINUE
C       replace WpR/C with the inversion permutations:
        DO 2020 K = 1, N
           ROW = WIR (K)
           COL = WIC (K)
           WPR (ROW) = K
           WPC (COL) = K
2020    CONTINUE

        IF (PGIVEN) THEN
C          the input matrix had been permuted from the original ordering
C          according to Rperm and Cperm.  Combine the initial
C          permutations (now in Rperm and Cperm) and the pivoting
C          permutations, and place them back into Rperm and Cperm.
           DO 2030 ROW = 1, N
              WM (WPR (ROW)) = RPERM (ROW)
2030       CONTINUE
           DO 2040 ROW = 1, N
              RPERM (ROW) = WM (ROW)
2040       CONTINUE
           DO 2050 COL = 1, N
              WM (WPC (COL)) = CPERM (COL)
2050       CONTINUE
           DO 2060 COL = 1, N
              CPERM (COL) = WM (COL)
2060       CONTINUE
C       else
C          the input matrix was not permuted on input.  Rperm and Cperm
C          in MA38ED have been passed to this routine as WiR and WiC,
C          which now contain the row and column permutations.  Rperm and
C          Cperm in this routine (MA38FD) are not defined.
        ENDIF

C       ----------------------------------------------------------------
C       allocate nlu+3 integers for xtail, nlu, npiv and LUp (1..nlu)
C       ----------------------------------------------------------------

        IS = NLU + 5
        LUIP1 = ITAIL
        ITAIL = ITAIL - IS
        IUSE = IUSE + IS
        INEED = IUSE
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .LE. ITAIL) THEN

C          -------------------------------------------------------------
C          sufficient memory exist to finish the factorization
C          -------------------------------------------------------------

           II (ITAIL+1) = NLU
           II (ITAIL+2) = NPIV
           LUPP = ITAIL+5
           IF (NLU .EQ. 0) THEN
C             zero the dummy entries, if LU factors are empty
              II (IP) = 0
              XX (XP) = ZERO
           ENDIF

C          -------------------------------------------------------------
C          convert the LU factors into the new pivot order
C          -------------------------------------------------------------

           S = 0
           MAXDR = 1
           MAXDC = 1
           DO 2100 K = 1, N
              E = WIR (K)
              LUIP = RP (E)
              IF (LUIP .GT. 0) THEN
C                this is an LU arrowhead - save a pointer in LUp:
                 S = S + 1
C                update pointers to LU arrowhead relative to start of LU
                 II (LUPP+S-1) = LUIP - LUIP1 + 1
                 LUXP = II (LUIP)
                 II (LUIP) = LUXP - XTAIL + 1
C                convert the row and column indices to their final order
C                pattern of a column of L:
                 P = (LUIP + 7)
                 LUDEGC = II (LUIP+3)
                 MAXDC = MAX (MAXDC, LUDEGC)
                 DO 2070 J = 1, LUDEGC
                    II (P) = WPR (ABS (II (P)))
                    P = P + 1
2070             CONTINUE
C                pattern of a row of U:
                 LUDEGR = II (LUIP+2)
                 MAXDR = MAX (MAXDR, LUDEGR)
                 DO 2080 J = 1, LUDEGR
                    II (P) = WPC (ABS (II (P)))
                    P = P + 1
2080             CONTINUE
C                convert the LUsons, Usons, and Lsons:
                 NSONS = II (LUIP+4)
                 DO 2090 J = 1, NSONS
                    ESON = II (P)
                    IF (ESON .LE. N) THEN
C                      an LUson
                       II (P) = WM (ESON)
                    ELSE IF (ESON .LE. 2*N) THEN
C                      a Uson
                       II (P) = WM (ESON-N) + N
                    ELSE
C                      an Lson
                       II (P) = WM (ESON-2*N) + 2*N
                    ENDIF
                    P = P + 1
2090             CONTINUE
C                renumber this LU arrowhead
                 WM (E) = S
              ENDIF
2100       CONTINUE

           CMAX = MAX (CMAX, MAXDC)
           RMAX = MAX (RMAX, MAXDR)
           TOTNLU = TOTNLU + NLU

           II (ITAIL+3) = MAXDC
           II (ITAIL+4) = MAXDR

C          -------------------------------------------------------------
C          get memory usage for next call to MA38BD
C          -------------------------------------------------------------

           XRUSE = XRUSE - NZ
           RETURN
        ENDIF

C=======================================================================
C  LU factors are now stored in their final form ]
C=======================================================================

C=======================================================================
C  Error conditions
C=======================================================================

C       error return label:
9000    CONTINUE
        IF (IHEAD .GT. ITAIL .OR. ISIZE .LT. MINMEM) THEN
C          error return if out of integer memory
           CALL MA38ND (1, ICNTL, INFO, -3, INFO (19))
        ENDIF
        IF (XHEAD .GT. XTAIL) THEN
C          error return if out of real memory
           CALL MA38ND (1, ICNTL, INFO, -4, INFO (21))
        ENDIF
        RETURN
        END

        SUBROUTINE MA38GD (XX, XSIZE, XHEAD, XUSE,
     *          II, ISIZE, IHEAD, IUSE,
     *          CP, RP, DN, N, WIR, WIC, WR, WC,
     *          FFXP, FFSIZE, WXP, FFDIMC, DOSLOT,
     *          PFREE, XFREE, MHEAD, MTAIL, SLOTS)
        INTEGER N, DN, ISIZE, II(ISIZE), IHEAD, RP(N+DN),
     *          CP(N+1), WIR(N), WIC(N), XSIZE, XUSE,
     *          IUSE, XHEAD, FFXP, FFSIZE, WXP,
     *          FFDIMC, WR(N+DN), WC(N+DN), PFREE, XFREE, MHEAD,
     *          MTAIL, SLOTS
        LOGICAL DOSLOT
        DOUBLE PRECISION XX (XSIZE)

C=== MA38GD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Garbage collection for MA38FD.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       II/XX:          integer/real workspace, containing matrix being
C                       factorized and partially-computed LU factors
C       isize:          size of II
C       xsize:          size of XX
C       xhead:          XX (1..xhead) is in use (matrix, frontal mtc's)
C       xtail:          XX (xtail..xsize) is in use (LU factors)
C       xuse:           memory usage in Value
C       ihead:          II (1..ihead) is in use (matrix, frontal mtc's)
C       itail:          II (itail..isize) is in use (LU factors)
C       iuse:           memory usage in Index
C       Cp (1..n+1):    pointers to columns
C       Rp (1..n+dn):   pointers to rows, frontal matrices, and LU
C                       arrowheads
C       dn:             number of dense columns
C       n:              order of matrix
C       Icntl:          integer control parameters, see MA38ID
C       Wr (1..n):      see MA38FD
C       Wc (1..n):      see MA38FD
C       ffxp:           pointer to current contribution block
C       ffsize:         size of current contribution block
C       mhead:          pointer to first block in memory list
C       mtail:          pointer to last block in memory list
C       doslot:         true if adding slots
C       if doslot:
C           WiR (1..n)  if WiR (row) >= 0 then add (or keep) an extra
C                       slot in the row's element list
C           WiC (1..n)  if WiR (col) >= 0 then add (or keep) an extra
C                       slot in the col's element list

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       II/XX:          external fragmentation is removed at head
C       xhead:          XX (1..xhead) is in use, reduced in size
C       xuse:           memory usage in Value, reduced
C       ihead:          II (1..ihead) is in use, reduced in size
C       iuse:           memory usage in Index, reduced
C       pfree:          pointer to free block in memory list, set to 0
C       xfree:          size of free block in XX, set to -1
C       mhead:          pointer to first block in memory list
C       mtail:          pointer to last block in memory list
C       ffxp            current working array has been shifted
C       wxp             current work vector has been shifted

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38FD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER WHAT, FSIZ, ROW, COL, P, IDP, XDP, I, E, EP, FDIMC,
     $          LUDEGR, LUDEGC, J, PC, CELN, CLEN, RELN, RLEN,
     $          CSIZ1, CSIZ2, RSIZ1, RSIZ2, FLUIP, CXP, FXP, RDEG,
     $          CDEG, CSCAL, RSCAL, FSCAL
        PARAMETER (CSCAL = 9, RSCAL = 2, FSCAL = 7)
        LOGICAL SLOT

C  Compression:
C  ------------
C  what:    what this block of memory is (a row, column, etc.)
C  idp:     int. destination pointer, current block moved to II (idp...)
C  xdp:     real destination pointer, current block moved to XX (xdp...)
C  slot:    true if adding, or keeping, a size-2 hole in an element list
C
C  Columns:
C  --------
C  cscal:   = 9, the number of scalars in column data structure
C  celn:    number of (e,f) tuples in element list of a column
C  clen:    number of unassembled original entries in a column
C  cdeg:    degree of a column (number of entries, including fill-in)
C  cxp:     a column is in XX (cxp...) prior to compression
C  pc:      column is in II (pc ...) prior to compression
C  csiz1:   size of a column in II, prior to compression
C  csiz2:   size of a column in II, after compression
C  col:     column index
C
C  Rows:
C  -----
C  rscal:   = 2, the number of scalars in row data structure
C  reln:    number of (e,f) tuples in element list of a row
C  rlen:    number of unassembled original entries in a row
C  rsiz1:   size of a row in II, prior to compression
C  rsiz2:   size of a row in II, after compression
C  rdeg:    degree of a row (number of entries, including fill-in)
C  row:     row index
C
C  Frontal matrices:
C  -----------------
C  fscal:   = 7, the number of scalars in element data structure
C  fluip:   element is in II (fluip...) prior to compression
C  fxp:     a frontal matrix is in XX (fxp...) prior to compression
C  e:       an element
C  fdimc:   column dimension (number of rows) of a frontal matrix
C  ludegr:  row degree (number of columns) of a contribution block
C  ludegc:  column degree (number of rows) of a contribution block
C  fsiz:    size of an artificial frontal matrix
C  ep:      an artificial frontal matrix is in II (ep ...) prior to comp
C
C  Other:
C  ------
C  p:       pointer
C  i:       general loop index
C  j:       general loop index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        SLOTS = 0

C-----------------------------------------------------------------------
C   prepare the non-pivotal rows/cols and unassembled elements
C-----------------------------------------------------------------------

C       place the size of each block of memory at the beginning,
C       and mark the 2nd entry in each block with what it is

C       ----------------------------------------------------------------
C       mark the columns
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 10 COL = 1, N
           PC = CP (COL)
           IF (PC .NE. 0) THEN
C             this is a non-pivotal, non-null column
              CDEG = II (PC+1)
              CP (COL) = CDEG
              II (PC+1) = COL+N
           ENDIF
10      CONTINUE

C       ----------------------------------------------------------------
C       mark the rows and frontal matrices
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 20 ROW = 1, N
           P = RP (ROW)
           RLEN = WC (ROW)
           IF (P .EQ. 0) THEN
C             a pivotal row
              CONTINUE
           ELSE IF (RLEN .GE. 0 .AND. RLEN .LE. N) THEN
C             this is a non-pivotal, non-null row
              RDEG = II (P+1)
              RP (ROW) = RDEG
              II (P+1) = ROW+2*N
           ELSE IF (WR (ROW) .EQ. -(N+DN+2)) THEN
C             a pivotal row, and an assembled element
              CONTINUE
           ELSE
C             this is an unassembled frontal matrix
C             the size is implicitly fscal
              FDIMC = II (P+1)
              RP (ROW) = FDIMC
              II (P+1) = ROW
           ENDIF
20      CONTINUE

C       ----------------------------------------------------------------
C       mark the artificial frontal matrices
C       ----------------------------------------------------------------

CFPP$ NODEPCHK L
        DO 30 E = N+1, N+DN
           EP = RP (E)
           IF (EP .NE. 0) THEN
C             this is an unassembled artificial frontal matrix
C             the size is II (ep+1) + cscal
              FDIMC = II (EP+1)
              RP (E) = FDIMC
              II (EP+1) = E+2*N
           ENDIF
30      CONTINUE

C-----------------------------------------------------------------------
C  scan the link list and compress the reals
C-----------------------------------------------------------------------

        XDP = 1
        P = MHEAD
C       while (p .ne. 0) do
40      CONTINUE
        IF (P .NE. 0) THEN

           WHAT = II (P+1)

C          -------------------------------------------------------------
           IF (WHAT .GT. 3*N) THEN
C          -------------------------------------------------------------

C             this is an unassembled artificial frontal matrix
              E = WHAT - 2*N
              FXP = II (P+2)
              II (P+2) = XDP
CFPP$ NODEPCHK L
              DO 50 J = 0, RP (E) - 1
                 XX (XDP+J) = XX (FXP+J)
50            CONTINUE
              XDP = XDP + RP (E)

C          -------------------------------------------------------------
           ELSE IF (WHAT .EQ. -1 .OR. II (P+6) .EQ. 0) THEN
C          -------------------------------------------------------------

C             this is a real hole - delete it from the link list
              IF (II (P+4) .NE. 0) THEN
                 II (II (P+4)+3) = II (P+3)
              ELSE
                 MHEAD = II (P+3)
              ENDIF
              IF (II (P+3) .NE. 0) THEN
                 II (II (P+3)+4) = II (P+4)
              ELSE
                 MTAIL = II (P+4)
              ENDIF

C          -------------------------------------------------------------
           ELSE IF (WHAT .LE. N) THEN
C          -------------------------------------------------------------

C             this is an unassembled frontal matrix
              E = WHAT
              FXP = II (P+2)
              II (P+2) = XDP
              FLUIP = II (P)
              LUDEGR = II (FLUIP+2)
              LUDEGC = II (FLUIP+3)
              FDIMC = RP (E)
              IF (FDIMC .EQ. LUDEGC) THEN
C                contribution block is already compressed
CFPP$ NODEPCHK L
                 DO 60 I = 0, (LUDEGR * LUDEGC) - 1
                    XX (XDP+I) = XX (FXP+I)
60               CONTINUE
              ELSE
C                contribution block is not compressed
C                compress XX (fxp..) to XX (xdp..xdp+(ludegr*ludegc)-1)
                 DO 80 J = 0, LUDEGR - 1
CFPP$ NODEPCHK L
                    DO 70 I = 0, LUDEGC - 1
                       XX (XDP + J*LUDEGC + I) = XX (FXP + J*FDIMC + I)
70                  CONTINUE
80               CONTINUE
                 RP (E) = LUDEGC
              ENDIF
              XDP = XDP + LUDEGR*LUDEGC

C          -------------------------------------------------------------
           ELSE IF (WHAT .LE. 2*N) THEN
C          -------------------------------------------------------------

C             this is a column
              CXP = II (P+2)
              II (P+2) = XDP
              CLEN = II (P+6)
CFPP$ NODEPCHK L
              DO 90 J = 0, CLEN - 1
                 XX (XDP+J) = XX (CXP+J)
90            CONTINUE
              XDP = XDP + CLEN

C          -------------------------------------------------------------
           ENDIF
C          -------------------------------------------------------------

C          -------------------------------------------------------------
C          get the next item in the link list
C          -------------------------------------------------------------

           P = II (P+3)

C       end while:
        GOTO 40
        ENDIF

        PFREE = 0
        XFREE = -1

C       ----------------------------------------------------------------
C       shift the current working array (if it exists)
C       ----------------------------------------------------------------

        IF (FFXP .NE. 0) THEN
CFPP$ NODEPCHK L
           DO 100 I = 0, FFSIZE - 1
              XX (XDP+I) = XX (FFXP+I)
100        CONTINUE
           FFXP = XDP
           XDP = XDP + FFSIZE
        ENDIF

C       ----------------------------------------------------------------
C       shift the current work vector (if it exists)
C       ----------------------------------------------------------------

        IF (WXP .NE. 0) THEN
           WXP = XDP
           XDP = XDP + FFDIMC
        ENDIF

C-----------------------------------------------------------------------
C  scan from the top of integer memory (1) to bottom (ihead) and
C  compress the integers
C-----------------------------------------------------------------------

        P = 1
        IDP = P
C       while (p .lt. ihead) do:
110     CONTINUE
        IF (P .LT. IHEAD) THEN

           WHAT = II (P+1)

C          -------------------------------------------------------------
           IF (WHAT .GT. 3*N) THEN
C          -------------------------------------------------------------

C             this is an unassembled artificial frontal matrix
              E = WHAT - 2*N
              FSIZ = RP (E) + CSCAL
              II (P+1) = RP (E)
              RP (E) = IDP
CFPP$ NODEPCHK L
              DO 120 I = 0, FSIZ - 1
                 II (IDP+I) = II (P+I)
120           CONTINUE
C             shift pointers in the link list
              IF (II (IDP+4) .NE. 0) THEN
                 II (II (IDP+4)+3) = IDP
              ELSE
                 MHEAD = IDP
              ENDIF
              IF (II (IDP+3) .NE. 0) THEN
                 II (II (IDP+3)+4) = IDP
              ELSE
                 MTAIL = IDP
              ENDIF
              P = P + FSIZ
              IDP = IDP + FSIZ

C          -------------------------------------------------------------
           ELSE IF (WHAT .EQ. -1) THEN
C          -------------------------------------------------------------

C             this is a integer hole
              P = P + II (P)

C          -------------------------------------------------------------
           ELSE IF (WHAT .GE. 1 .AND. WHAT .LE. N) THEN
C          -------------------------------------------------------------

C             this is an unassembled frontal matrix (fscal integers)
              E = WHAT
              FDIMC = RP (E)
              II (P+1) = FDIMC
              RP (E) = IDP
CFPP$ NODEPCHK L
              DO 130 I = 0, FSCAL - 1
                 II (IDP+I) = II (P+I)
130           CONTINUE
C             shift pointers in the link list
              IF (II (IDP+4) .NE. 0) THEN
                 II (II (IDP+4)+3) = IDP
              ELSE
                 MHEAD = IDP
              ENDIF
              IF (II (IDP+3) .NE. 0) THEN
                 II (II (IDP+3)+4) = IDP
              ELSE
                 MTAIL = IDP
              ENDIF
              P = P + FSCAL
              IDP = IDP + FSCAL

C          -------------------------------------------------------------
           ELSE IF (WHAT .LE. 2*N) THEN
C          -------------------------------------------------------------

C             this is a non-pivotal column
              CSIZ1 = II (P)
              COL = WHAT - N
              CELN = II (P+5)
              CLEN = II (P+6)
              CSIZ2 = 2*CELN + CLEN + CSCAL
              SLOT = DOSLOT .AND. WIC (COL) .GE. 0 .AND. P .GE. IDP+2
              IF (SLOT) THEN
C                keep (or make) one extra slot for element list growth
                 CSIZ2 = CSIZ2 + 2
                 SLOTS = SLOTS + 2
              ENDIF
              CDEG = CP (COL)
              II (P+1) = CDEG
              CP (COL) = IDP
              II (P) = CSIZ2
C             copy the cscal scalars and the celn (e,f) tuples
CFPP$ NODEPCHK L
              DO 140 I = 0, CSCAL + 2*CELN - 1
                 II (IDP+I) = II (P+I)
140           CONTINUE
              IF (CLEN .GT. 0) THEN
C                shift pointers in the link list
                 IF (II (IDP+4) .NE. 0) THEN
                    II (II (IDP+4)+3) = IDP
                 ELSE
                    MHEAD = IDP
                 ENDIF
                 IF (II (IDP+3) .NE. 0) THEN
                    II (II (IDP+3)+4) = IDP
                 ELSE
                    MTAIL = IDP
                 ENDIF
              ENDIF
              P = P + CSIZ1 - CLEN
              IDP = IDP + CSCAL + 2*CELN
              IF (SLOT) THEN
C                skip past the slot
                 IDP = IDP + 2
              ENDIF
C             copy the clen original row indices
CFPP$ NODEPCHK L
              DO 150 I = 0, CLEN - 1
                 II (IDP+I) = II (P+I)
150           CONTINUE
              P = P + CLEN
              IDP = IDP + CLEN

C          -------------------------------------------------------------
           ELSE
C          -------------------------------------------------------------

C             this is a non-pivotal row
              RSIZ1 = II (P)
              ROW = WHAT - 2*N
              RELN = WR (ROW)
              RLEN = WC (ROW)
              RSIZ2 = 2*RELN + RLEN + RSCAL
              SLOT = DOSLOT .AND. WIR (ROW) .GE. 0 .AND. P .GE. IDP+2
              IF (SLOT) THEN
C                keep (or make) one extra slot for element list growth
                 RSIZ2 = RSIZ2 + 2
                 SLOTS = SLOTS + 2
              ENDIF
              RDEG = RP (ROW)
              II (P+1) = RDEG
              RP (ROW) = IDP
              II (P) = RSIZ2
C             copy the rscal scalars, and the reln (e,f) tuples
CFPP$ NODEPCHK L
              DO 160 I = 0, RSCAL + 2*RELN - 1
                 II (IDP+I) = II (P+I)
160           CONTINUE
              P = P + RSIZ1 - RLEN
              IDP = IDP + RSCAL + 2*RELN
              IF (SLOT) THEN
C                skip past the slot
                 IDP = IDP + 2
              ENDIF
C             copy the rlen original column indices
CFPP$ NODEPCHK L
              DO 170 I = 0, RLEN - 1
                 II (IDP+I) = II (P+I)
170           CONTINUE
              P = P + RLEN
              IDP = IDP + RLEN

C          -------------------------------------------------------------
           ENDIF
C          -------------------------------------------------------------

C          -------------------------------------------------------------
C          move to the next block
C          -------------------------------------------------------------

C       end while:
        GOTO 110
        ENDIF

C-----------------------------------------------------------------------
C  deallocate the unused space
C-----------------------------------------------------------------------

        IUSE = IUSE - (IHEAD - IDP)
        IHEAD = IDP
        XUSE = XUSE - (XHEAD - XDP)
        XHEAD = XDP
        RETURN
        END

        SUBROUTINE MA38HD (XX, XSIZE, II, ISIZE, N, NZ, NZDIA, NZOFF,
     *          NBLKS, CP, CPERM, RPERM, PR, PC,
     *          W, ZPERM, BP, OFFP,
     *          PRESRV)
        INTEGER N, NZ, ISIZE, II(ISIZE), NZDIA, NZOFF, NBLKS, CP(N+1),
     *          CPERM(N), RPERM(N), PR(N), PC(N), W(N), ZPERM(N),
     *          BP(N+1), OFFP(N+1), XSIZE
        LOGICAL PRESRV
        DOUBLE PRECISION XX (XSIZE)

C=== MA38HD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Find permutations to block triangular form:
C       1) permute the matrix so that it has a zero-free diagonal.
C       2) find strongly-connected components of the corresponding
C          graph.  Each diagonal block corresponds to exactly one
C          strongly-connected component.
C       3) convert the matrix to block triangular form, unless it is
C          to be preserved in its original form.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       presrv:         true if original matrix is to be preserved
C       n:              order of matrix
C       nz:             entries in matrix
C       isize:          size of II
C       xsize:          size of XX
C       Cp (1..n+1):    column pointers
C       XX (1..nz):     values
C       II (1..nz):     row indices
C       Icntl:          integer control arguments
C
C          input matrix in column form is in:
C          XX (1..nz), II (1..nz), n, nz, Cp (1..n+1), where
C               II (Cp(col) ... Cp(col+1)-1): row indices
C               XX (Cp(col) ... Cp(col+1)-1): values
C          if presrv is false then xsize and isize must be >= 2*nz
C          otherwise, xsize and isize must be >= nz

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       Pr (1..n), Pc (1..n), W (1..n), Zperm (1..n)

C======================================================================
C  OUTPUT:
C=======================================================================
C
C       nblks: number of blocks
C       if (nblks > 1):
C
C           Cperm (1..n), Rperm (1..n): permutation to block form:
C               Rperm (newrow) = oldrow
C               Cperm (newcol) = oldcol
C
C           Bp (n-nblks+1...n+1) holds the start/end of blocks 1..nblks
C
C           if (presrv is false) then
C
C              input matrix is converted to block-upper-tri. form,
C              using II/XX (nz+1..2*nz) as workspace.
C              nzdia: nonzeros in diagonal blocks
C              nzoff: nonzeros in off-diagonal blocks
C              (nz = nzdia + nzoff)
C
C              off-diagonal column-oriented form in XX/II (1..nzoff)
C              col is located in
C              XX/II (Offp (col) ... Offp (col+1)-1)
C
C              diagonal blocks now in XX/II (nzoff+1 .. nzoff+nzdia)
C              col is located in
C              XX/II (Cp (col) ... Cp (col+1)-1)
C
C       else, nblks=1: and no other output is generated.

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38DD
C       subroutines called:     MC21BD, MC13ED
        EXTERNAL MC21BD, MC13ED

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER COL, NDIAG, I, PO, PB, BLK, P, ROW, K1, K

C  ndiag:   number of entries on the diagonal
C  po:      pointer into off-diagonal part
C  pb:      pointer into diagonal blocks
C  blk:     block number for current diagonal block
C  k1:      column col is in diagonal block A (k1.., k1...)
C  k:       kth row/col in BTF form is Rperm(k)/Cperm(k) in input matrix
C  p:       pointer
C  row:     row index
C  col:     column index
C  i:       general loop index


C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        NZDIA = NZ
        NZOFF = 0

C-----------------------------------------------------------------------
C compute the length of each column
C-----------------------------------------------------------------------

        DO 10 COL = 1, N
           W (COL) = CP (COL+1) - CP (COL)
10      CONTINUE

C-----------------------------------------------------------------------
C find a column permutation for a zero-free diagonal
C-----------------------------------------------------------------------

        CALL MC21BD(N, II, NZ, CP, W, ZPERM, NDIAG, OFFP, CPERM, PR, PC)
C          MC21BD calling interface:
C          input:       n, II (1..nz), nz, Cp (n), W (n):
C                       n-by-n matrix, col is of length W (col),
C                       and its pattern is located in
C                       II (Cp (col) ... Cp (col)+W(col)-1)
C          output:      Zperm (n), the permutation, such that
C                       colold = Zperm (col), and ndiag (number of
C                       structural nonzeros on the diagonal.
C                       matrix is structurally singular if ndiag < n
C          workspace:   Offp, Cperm, Pr, Pc

C-----------------------------------------------------------------------
C  permute the columns of the temporary matrix to get zero-free diagonal
C-----------------------------------------------------------------------

        DO 20 COL = 1, N
           OFFP (COL) = CP (ZPERM (COL))
           W (COL) = CP (ZPERM (COL)+1) - CP (ZPERM (COL))
20      CONTINUE

C-----------------------------------------------------------------------
C  find a symmetric permutation into upper block triangular form
C  (that is, find the strongly-connected components in the graph).
C-----------------------------------------------------------------------

        CALL MC13ED(N, II, NZ, OFFP, W, RPERM, BP, NBLKS, CPERM, PR, PC)
C          MC13ED calling interface:
C          input:       n, II (1..nz), nz, Offp (n), W (n)
C                       n-by-n matrix, col of length W(col),
C                       in II (Offp(col) ... Offp(col)+W(col)-1), where
C                       this permuted matrix has a zero-free diagonal
C                       (unless the matrix is structurally singular).
C          output:      Rperm (n), Bp (n+1), nblks
C                       old = Rperm (new) is the symmetric permutation,
C                       there are nblks diagonal blocks, Bp (i) is
C                       the position in new order of the ith block.
C          workspace:   Cperm, Pr, Pc

C-----------------------------------------------------------------------
C  if more than one block, get permutations and block pointers,
C  and convert to block-upper-triangular form (unless matrix preserved)
C-----------------------------------------------------------------------

        IF (NBLKS .NE. 1) THEN

C          -------------------------------------------------------------
C          find the composite column permutation vector (Cperm):
C          -------------------------------------------------------------

           DO 30 COL = 1, N
              CPERM (COL) = ZPERM (RPERM (COL))
30         CONTINUE

C          -------------------------------------------------------------
C          convert to block-upper-triangular form, if not preserved
C          -------------------------------------------------------------

           IF (.NOT. PRESRV) THEN

C             ----------------------------------------------------------
C             find the inverse permutation vectors, Pr and Pc
C             ----------------------------------------------------------

              DO 40 K = 1, N
                 PC (CPERM (K)) = K
                 PR (RPERM (K)) = K
40            CONTINUE

C             ----------------------------------------------------------
C             construct flag array to determine if entry in block or not
C             ----------------------------------------------------------

              BP (NBLKS+1) = N+1
              DO 60 BLK = 1, NBLKS
                 DO 50 I = BP (BLK), BP (BLK+1)-1
                    W (I) = BP (BLK)
50               CONTINUE
60            CONTINUE

C             ----------------------------------------------------------
C             construct block-diagonal form in XX/II (nz+1..nz+nzdia)
C             ----------------------------------------------------------

C             These blocks are in a permuted order (according to Rperm
C             and Cperm).  The row indices in each block range from 1
C             to the size of the block.

              PB = NZ + 1
              DO 80 COL = 1, N
                 ZPERM (COL) = PB
                 K1 = W (COL)
CFPP$ NODEPCHK L
                 DO 70 P = CP (CPERM (COL)), CP (CPERM (COL)+1)-1
                    ROW = PR (II (P))
                    IF (W (ROW) .EQ. K1) THEN
C                      entry is in the diagonal block:
                       II (PB) = ROW - K1 + 1
                       XX (PB) = XX (P)
                       PB = PB + 1
                    ENDIF
70               CONTINUE
80            CONTINUE
C             Zperm (n+1) == pb  ( but Zperm (n+1) does not exist )
              NZDIA = PB - (NZ + 1)
              NZOFF = NZ - NZDIA

C             diagonal blocks now in XX/II (nz+1..nz+nzdia)
C             col is located in XX/II (Zperm (col) ... Zperm (col+1)-1)

C             ----------------------------------------------------------
C             compress original matrix to off-diagonal part, in place
C             ----------------------------------------------------------

C             The rows/cols of off-diagonal form correspond to rows/cols
C             in the original, unpermuted matrix.  They are permuted to
C             the final pivot order and stored in a row-oriented form,
C             after the factorization is complete (by MA38MD).

              PO = 1
              DO 100 COL = 1, N
                 OFFP (COL) = PO
                 K1 = W (PC (COL))
CFPP$ NODEPCHK L
                 DO 90 P = CP (COL), CP (COL+1)-1
                    ROW = PR (II (P))
                    IF (W (ROW) .NE. K1) THEN
C                      offdiagonal entry
                       II (PO) = II (P)
                       XX (PO) = XX (P)
                       PO = PO + 1
                    ENDIF
90               CONTINUE
100           CONTINUE
              OFFP (N+1) = PO

C             off-diagonal form now in XX/II (1..nzoff)
C             col is located in XX/II(Offp(col)..Offp(col+1)-1)

C             ----------------------------------------------------------
C             move block-diagonal part into place
C             ----------------------------------------------------------

              PB = NZ + 1
CFPP$ NODEPCHK L
              DO 110 I = 0, NZDIA - 1
                 II (PO+I) = II (PB+I)
                 XX (PO+I) = XX (PB+I)
110           CONTINUE
              DO 120 COL = 1, N
                 CP (COL) = ZPERM (COL) - NZDIA
120           CONTINUE
C             Cp (n+1) == nz+1  ( this is unchanged )

C             diagonal blocks now in XX/II (nzoff+1 .. nzoff+nzdia)
C             col is located in XX/II (Cp (col) ... Cp (col+1)-1)

           ENDIF

C          -------------------------------------------------------------
C          shift Bp (1 .. nblks+1) down to Bp (1+n-nblks .. n+1), which
C          then becomes the Blkp (1 .. nblks+1) array.
C          -------------------------------------------------------------

           BP (NBLKS+1) = N+1
CFPP$ NODEPCHK L
           DO 130 BLK = NBLKS + 1, 1, -1
              BP (BLK + (N-NBLKS)) = BP (BLK)
130        CONTINUE
        ENDIF

        RETURN
        END

        SUBROUTINE MA38JD (N, JOB, TRANSC, LUXSIZ, LUX,
     *          LUISIZ, LUI, B, X, R, Z, LY, Y, S, CNTL, INFO,
     *          RINFO, CPERM, RPERM, AN, ANZ, AP, AI, AX, ON,
     *          NZOFF, OFFP, OFFI, OFFX, NBLKS, LUBLKP, BLKP, IRSTEP)
        INTEGER N, JOB, LUXSIZ, LUISIZ, LUI(LUISIZ), LY, IRSTEP,
     *          INFO(40), CPERM(N), RPERM(N), AN,
     *          ANZ, AP(AN+1), AI(ANZ), ON, NZOFF, OFFP(ON+1),
     *          OFFI(NZOFF), NBLKS, LUBLKP(NBLKS), BLKP(NBLKS+1)
        LOGICAL TRANSC
        DOUBLE PRECISION LUX(LUXSIZ), B(N), X(N), R(N), Z(N), Y(LY),
     *          S(LY), CNTL(10), RINFO(20), AX(ANZ),
     *          OFFX(NZOFF)

C=== MA38JD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Solve a system, given LU factors, permutation arrays, original
C  matrix (if preserved), and off-diagonal blocks (if BTF was used).

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       job:            0: solve Ax=b, 1: solve Lx=b, 2: solve Ux=b
C       transc:         if true, solve with transposed factors instead
C       LUxsiz:         size of LUx
C       LUx (1..LUxsiz) real values in LU factors for each block
C       LUisiz:         size of LUi
C       LUi (1..LUisiz) integers in LU factors for each block
C       B (1..n):       right-hand-side
C       ly:             size of Y and S, ly=n if Y and S are used
C       Cntl:           real control parameters, see MA38ID
C       Icntl:          integer control parameters, see MA38ID
C       Cperm (1..n):   Q, column permutation array
C       Rperm (1..n):   P, row permutation array
C       presrv:         if true, then original matrix was preserved
C       nblks:          number of diagonoal blocks (1 if no BTF)
C       irstep:         maximum number of steps of iterative refinement
C
C       if presrv then
C           an:                 order of preserved matrix, n
C           anz:                number of entries in preserved matrix
C           Ap (1..an+1):       column pointers of preserved matrix
C           Ai (1..anz):        row indices of preserved matrix
C           Ax (1..anz):        values of preserved matrix
C           an, anz, Ap, Ai, Ax:        not accessed
C
C       if nblks > 1 then
C           on:                 n
C           nzoff:              number of off-diagonoal entries
C           Offp (1..n+1)       row pointers for off-diagonal part
C           Offi (1..nzoff):    column indices for off-diagonal part
C           Offx (1..nzoff):    values of off-diagonal part
C           LUBlkp (1..nblks):  pointers to LU factors of each block
C           Blkp (1..nblks+1):  index range of each block
C       else
C           on, nzoff, Offp, Offi, Offx, LUBlkp, Blkp:  not accessed

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       R (1..n), Z (1..n)
C       Y (1..ly), S (1..ly):   unaccessed if no iterative refinement

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       solution
C       Info:           integer informational output, see MA38ID
C       Rinfo:          real informational output, see MA38ID
C
C       if irsteps > 0 and presrv is true then
C           W (1..n):           residual
C           Rinfo (7):  sparse error estimate, omega1
C           Rinfo (8):  sparse error estimate, omega2

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38CD
C       subroutines called:     MA38ND, MA38LD, MA38TD, MA38UD, MA38SD
C       functions called:       IDAMAX, ABS, MAX
        INTRINSIC ABS, MAX
        INTEGER IDAMAX
        EXTERNAL MA38LD, MA38ND, MA38SD, MA38TD, MA38UD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER NLU, I, BLK, K1, K2, KN, P, STEP, NPIV, J
        DOUBLE PRECISION
     $          ZERO, ONE, XNORM, TAU, NCTAU, OMEGA1, OMEGA2, D1,
     $          D2, OMEGA, OMLAST, OM1LST, OM2LST, TWO, EPS, MAXEPS,
     $          THOSND, A, AXX, R2, X2, Y2, Z2
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $          MAXEPS = TWO ** (-15), THOSND = 1000.0D0)

C  LU factors:
C  -----------
C  blk:     current diagonal block
C  k1,k2:   current diagonal block is A (k1..k2, k1..k2)
C  kn:      size of diagonal block (= k2-k1+1)
C  nlu:     number of elements in the LU factors of a single diag block
C  npiv:    number of pivots in the LU factors of a single diag block
C
C  Iterative refinement and sparse backward error:
C  -----------------------------------------------
C  step:    number of steps of iterative refinement taken
C  xnorm:   ||x|| maxnorm of solution vector, x
C  tau:     threshold for selecting which estimate to use (1 or 2)
C  nctau:   1000*n*eps
C  eps:     largest positive value such that fl (1.0 + eps) = 1.0
C  omega1:  current sparse backward error estimate 1
C  omega2:  current sparse backward error estimate 2
C  d1:      divisor for omega1
C  d2:      divisor for omega2
C  omega:   omega1 + omega2
C  omlast:  value of omega from previous step
C  om1lst:  value of omega1 from previous step
C  om2lst:  value of omega2 from previous step
C  maxeps:  maximum value that eps is allowed to be
C  a:       value of an entry in A, A_ij
C  axx:     A_ij * x_j
C
C  Other:
C  ------
C  i,j:     loop indices
C  p:       pointer
C  r2:      R (i)
C  x2:      X (i)
C  y2:      Y (i)
C  z2:      Z (i)

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  initializations for sparse backward error
C-----------------------------------------------------------------------

        OMEGA = ZERO
        OMEGA1 = ZERO
        OMEGA2 = ZERO
        EPS = CNTL (3)
        IF (EPS .LE. ZERO .OR. EPS .GT. MAXEPS) THEN
C          eps is too small or too big: set to a large default value
           EPS = MAXEPS
        ENDIF
        NCTAU = THOSND * N * EPS

C-----------------------------------------------------------------------
C  get information on LU factorization if BTF was not used
C-----------------------------------------------------------------------

        IF (NBLKS .EQ. 1) THEN
C          p is 1, and LUi (p) is 1
           NLU = LUI (2)
           NPIV = LUI (3)
        ENDIF

C-----------------------------------------------------------------------
        IF (JOB .EQ. 1) THEN
C-----------------------------------------------------------------------

C          -------------------------------------------------------------
           IF (.NOT. TRANSC) THEN
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve P'Lx=b:  x = L \ Pb
C             ----------------------------------------------------------

              DO 10 I = 1, N
                 X (I) = B (RPERM (I))
10            CONTINUE
              IF (NBLKS .EQ. 1) THEN
                 CALL MA38LD (NLU, N, LUI(6), LUI(6+NLU), LUX,X,Z)
              ELSE
                 DO 20 BLK = 1, NBLKS
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .GT. 1) THEN
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38LD (NLU, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                    ENDIF
20               CONTINUE
              ENDIF

C          -------------------------------------------------------------
           ELSE
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve L'Px=b:  x = P' (L' \ b)
C             ----------------------------------------------------------

              DO 30 I = 1, N
                 R (I) = B (I)
30            CONTINUE
              IF (NBLKS .EQ. 1) THEN
                 CALL MA38TD (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,R,Z)
              ELSE
                 DO 40 BLK = 1, NBLKS
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .GT. 1) THEN
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38TD (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF
40               CONTINUE
              ENDIF
              DO 50 I = 1, N
                 X (RPERM (I)) = R (I)
50            CONTINUE

C          -------------------------------------------------------------
           ENDIF
C          -------------------------------------------------------------

C-----------------------------------------------------------------------
        ELSE IF (JOB .EQ. 2) THEN
C-----------------------------------------------------------------------

C          -------------------------------------------------------------
           IF (TRANSC) THEN
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve QU'x=b:  x = U' \ Q'b
C             ----------------------------------------------------------

              DO 60 I = 1, N
                 X (I) = B (CPERM (I))
60            CONTINUE
              IF (NBLKS .EQ. 1) THEN
                 CALL MA38SD (NLU, N, LUI(6), LUI(6+NLU), LUX,X,Z)
              ELSE
                 DO 100 BLK = 1, NBLKS
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .EQ. 1) THEN
                       X (K1) = X (K1) / LUX (LUBLKP (BLK))
                       R (K1) = X (K1)
                    ELSE
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38SD (NLU, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                       DO 70 I = K1, K2
                          R (I) = X (I)
70                     CONTINUE
                       CALL MA38TD (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF
                    DO 90 I = K1, K2
                       R2 = R (I)
                       DO 80 P = OFFP (I), OFFP (I+1)-1
                          X (OFFI (P)) = X (OFFI (P)) - OFFX (P) * R2
80                     CONTINUE
90                  CONTINUE
100              CONTINUE
              ENDIF

C          -------------------------------------------------------------
           ELSE
C          -------------------------------------------------------------

C             ----------------------------------------------------------
C             Solve UQ'x=b:  x = Q (U \ b)
C             ----------------------------------------------------------

              IF (NBLKS .EQ. 1) THEN
                 DO 110 I = 1, N
                    R (I) = B (I)
110              CONTINUE
                 CALL MA38UD (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,R,Z)
              ELSE
                 DO 150 BLK = NBLKS, 1, -1
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    DO 130 I = K1, K2
                       X2 = ZERO
                       DO 120 P = OFFP (I), OFFP (I+1)-1
                          X2 = X2 + OFFX (P) * R (OFFI (P))
120                    CONTINUE
                       X (I) = X2
130                 CONTINUE
                    IF (KN .EQ. 1) THEN
                       R (K1) = (B (K1) - X (K1)) / LUX (LUBLKP (BLK))
                    ELSE
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38LD (NLU, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                       DO 140 I = K1, K2
                          R (I) = B (I) - X (I)
140                    CONTINUE
                       CALL MA38UD (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF
150              CONTINUE
              ENDIF
              DO 160 I = 1, N
                 X (CPERM (I)) = R (I)
160           CONTINUE

C          -------------------------------------------------------------
           ENDIF
C          -------------------------------------------------------------

C-----------------------------------------------------------------------
        ELSE
C-----------------------------------------------------------------------

           DO 450 STEP = 0, IRSTEP

C             ----------------------------------------------------------
C             If transa was true in MA38AD or MA38BD, then C = A'.
C             Otherwise C = A.  In both cases, the factorization is
C             PCQ = LU, and C is stored in column-form in Ai,Ax,Ap if
C             it is preserved.
C             ----------------------------------------------------------

C             ----------------------------------------------------------
              IF (.NOT. TRANSC) THEN
C             ----------------------------------------------------------

C                -------------------------------------------------------
C                Solve Cx=b (step 0):
C                   x = Q (U \ L \ Pb)
C                and then perform iterative refinement (step > 0):
C                   x = x + Q (U \ L \ P (b-Cx))
C                -------------------------------------------------------

                 IF (STEP .EQ. 0) THEN
                    DO 170 I = 1, N
                       R (I) = B (RPERM (I))
170                 CONTINUE
                 ELSE
                    DO 180 I = 1, N
                       Z (I) = B (I)
180                 CONTINUE
                    DO 200 I = 1, N
                       X2 = X (I)
                       DO 190 P = AP (I), AP (I+1) - 1
                          Z (AI (P)) = Z (AI (P)) - AX (P) * X2
190                    CONTINUE
200                 CONTINUE
                    DO 210 I = 1, N
                       R (I) = Z (RPERM (I))
210                 CONTINUE
                 ENDIF
                 IF (NBLKS .EQ. 1) THEN
                    CALL MA38LD (NLU, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                    CALL MA38UD (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                 ELSE
                    DO 240 BLK = NBLKS, 1, -1
                       K1 = BLKP (BLK)
                       K2 = BLKP (BLK+1) - 1
                       KN = K2-K1+1
                       DO 230 I = K1, K2
                          R2 = R (I)
                          DO 220 P = OFFP (I), OFFP (I+1)-1
                             R2 = R2 - OFFX (P) * R (OFFI (P))
220                       CONTINUE
                          R (I) = R2
230                    CONTINUE
                       IF (KN .EQ. 1) THEN
                          R (K1) = R (K1) / LUX (LUBLKP (BLK))
                       ELSE
                          P = LUBLKP (BLK)
                          NLU = LUI (P+1)
                          NPIV = LUI (P+2)
                          CALL MA38LD (NLU, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                          CALL MA38UD (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                       ENDIF
240                 CONTINUE
                 ENDIF
                 IF (STEP .EQ. 0) THEN
                    DO 250 I = 1, N
                       X (CPERM (I)) = R (I)
250                 CONTINUE
                 ELSE
                    DO 260 I = 1, N
                       X (CPERM (I)) = X (CPERM (I)) + R (I)
260                 CONTINUE
                 ENDIF

C             ----------------------------------------------------------
              ELSE
C             ----------------------------------------------------------

C                -------------------------------------------------------
C                Solve C'x=b (step 0):
C                   x = P' (L' \ U' \ Q'b)
C                and then perform iterative refinement (step > 0):
C                   x = x + P' (L' \ U' \ Q' (b-C'x))
C                -------------------------------------------------------

                 IF (STEP .EQ. 0) THEN
                    DO 270 I = 1, N
                       R (I) = B (CPERM (I))
270                 CONTINUE
                 ELSE
                    DO 280 I = 1, N
                       Z (I) = B (I)
280                 CONTINUE
                    DO 300 I = 1, N
                       Z2 = Z (I)
                       DO 290 P = AP (I), AP (I+1) - 1
                          Z2 = Z2 - AX (P) * X (AI (P))
290                    CONTINUE
                       Z (I) = Z2
300                 CONTINUE
                    DO 310 I = 1, N
                       R (I) = Z (CPERM (I))
310                 CONTINUE
                 ENDIF
                 IF (NBLKS .EQ. 1) THEN
                    CALL MA38SD (NLU, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                    CALL MA38TD (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                 ELSE
                    DO 340 BLK = 1, NBLKS
                       K1 = BLKP (BLK)
                       K2 = BLKP (BLK+1) - 1
                       KN = K2-K1+1
                       IF (KN .EQ. 1) THEN
                          R (K1) = R (K1) / LUX (LUBLKP (BLK))
                       ELSE
                          P = LUBLKP (BLK)
                          NLU = LUI (P+1)
                          NPIV = LUI (P+2)
                          CALL MA38SD (NLU, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                          CALL MA38TD (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                       ENDIF
                       DO 330 I = K1, K2
                          R2 = R (I)
                          DO 320 P = OFFP (I), OFFP (I+1)-1
                             R (OFFI (P)) = R (OFFI (P)) - OFFX (P) * R2
320                       CONTINUE
330                    CONTINUE
340                 CONTINUE
                 ENDIF
                 IF (STEP .EQ. 0) THEN
                    DO 350 I = 1, N
                       X (RPERM (I)) = R (I)
350                 CONTINUE
                 ELSE
                    DO 360 I = 1, N
                       X (RPERM (I)) = X (RPERM (I)) + R (I)
360                 CONTINUE
                 ENDIF

C             ----------------------------------------------------------
              ENDIF
C             ----------------------------------------------------------

C             ----------------------------------------------------------
C             sparse backward error estimate
C             ----------------------------------------------------------

              IF (IRSTEP .GT. 0) THEN

C                xnorm = ||x|| maxnorm
                 XNORM = ABS (X (IDAMAX (N, X, 1)))

C                r (i) = (b-Ax)_i, residual (or A')
C                z (i) = (|A||x|)_i
C                y (i) = ||A_i||, maxnorm of row i of A (or A')
                 DO 370 I = 1, N
                    R (I) = B (I)
                    Z (I) = ZERO
                    Y (I) = ZERO
370              CONTINUE

                 IF (.NOT. TRANSC) THEN

C                   ----------------------------------------------------
C                   sparse backward error for Cx=b, C stored by column
C                   ----------------------------------------------------

                    DO 390 J = 1, N
                       X2 = X (J)
CFPP$ NODEPCHK L
                       DO 380 P = AP (J), AP (J+1) - 1
                          I = AI (P)
                          A = AX (P)
                          AXX = A * X2
                          R (I) = R (I) -     (AXX)
                          Z (I) = Z (I) + ABS (AXX)
                          Y (I) = Y (I) + ABS (A)
380                    CONTINUE
390                 CONTINUE

                 ELSE

C                   ----------------------------------------------------
C                   sparse backward error for C'x=b, C' stored by row
C                   ----------------------------------------------------

                    DO 410 I = 1, N
                       R2 = R (I)
                       Z2 = Z (I)
                       Y2 = Y (I)
CFPP$ NODEPCHK L
                       DO 400 P = AP (I), AP (I+1) - 1
                          J = AI (P)
                          A = AX (P)
                          AXX = A * X (J)
                          R2 = R2 -     (AXX)
                          Z2 = Z2 + ABS (AXX)
                          Y2 = Y2 + ABS (A)
400                    CONTINUE
                       R (I) = R2
                       Z (I) = Z2
                       Y (I) = Y2
410                 CONTINUE

                 ENDIF

C                -------------------------------------------------------
C                save the last iteration in case we need to reinstate it
C                -------------------------------------------------------

                 OMLAST = OMEGA
                 OM1LST = OMEGA1
                 OM2LST = OMEGA2

C                -------------------------------------------------------
C                compute sparse backward errors: omega1 and omega2
C                -------------------------------------------------------

                 OMEGA1 = ZERO
                 OMEGA2 = ZERO
                 DO 420 I = 1, N
                    TAU = (Y (I) * XNORM + ABS (B (I))) * NCTAU
                    D1 = Z (I) + ABS (B (I))
                    IF (D1 .GT. TAU) THEN
                       OMEGA1 = MAX (OMEGA1, ABS (R (I)) / D1)
                    ELSE IF (TAU .GT. ZERO) THEN
                       D2 = Z (I) + Y (I) * XNORM
                       OMEGA2 = MAX (OMEGA2, ABS (R (I)) / D2)
                    ENDIF
420              CONTINUE
                 OMEGA = OMEGA1 + OMEGA2
                 RINFO (7) = OMEGA1
                 RINFO (8) = OMEGA2

C                -------------------------------------------------------
C                stop the iterations if the backward error is small
C                -------------------------------------------------------

                 INFO (24) = STEP
                 IF (ONE + OMEGA .LE. ONE) THEN
C                   further iterative refinement will no longer improve
C                   the solution
                    RETURN
                 ENDIF

C                -------------------------------------------------------
C                stop if insufficient decrease in omega
C                -------------------------------------------------------

                 IF (STEP .GT. 0 .AND. OMEGA .GT. OMLAST / TWO) THEN
                    IF (OMEGA .GT. OMLAST) THEN
C                      last iteration better than this one, reinstate it
                       DO 430 I = 1, N
                          X (I) = S (I)
                          RINFO (7) = OM1LST
                          RINFO (8) = OM2LST
430                    CONTINUE
                    ENDIF
                    INFO (24) = STEP - 1
                    RETURN
                 ENDIF

C                -------------------------------------------------------
C                save current solution in case we need to reinstate
C                -------------------------------------------------------

                 DO 440 I = 1, N
                    S (I) = X (I)
440              CONTINUE

              ENDIF

450        CONTINUE

C-----------------------------------------------------------------------
        ENDIF
C-----------------------------------------------------------------------

        RETURN
        END

        SUBROUTINE MA38KD (N, NZ, TRANSA, XX, XSIZE, INFO, ICNTL,
     *                     II, ISIZE, W, WP, WHO)
        INTEGER ISIZE, II(ISIZE), N, NZ, W(N), WP(N+1), INFO(40),
     *          ICNTL(20), XSIZE,WHO
        DOUBLE PRECISION XX (XSIZE)
        LOGICAL TRANSA

C=== MA38KD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Convert input matrix (II,XX,n,nz) into from triplet form to column-
C  oriented form.  Remove invalid entries and duplicate entries.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              size of the matrix
C       nz:             number of nonzeros in the input matrix
C       transa:         if true then transpose input matrix
C       XX (1..nz):     values of triplet form
C       xsize:          size of XX, must be >= 2*nz
C       II (1..2*nz):   row and col indices of triplet form
C       isize:          size of II, must be >= max (2*nz,n+1) + nz
C       Icntl:          integer control parameters
C       who:            who called MA38KD, 1: MA38AD, 2: MA38BD
C
C       II must be at least of size (nz + max (2*nz, n+1))
C       XX must be at least of size (nz + max (  nz, n+1))
C
C       input triplet matrix:
C          if (transa) is false:
C               II (p)          row index, for p = 1..nz
C               II (nz+p)       col index
C               XX (p)          value
C          if (transa) is true:
C               II (p)          col index, for p = 1..nz
C               II (nz+p)       row index
C               XX (p)          value

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       nz:             number of nonzeros in the output matrix,
C                       after removing invalid entries, and summing up
C                       duplicate entries
C       II (n+2..nz+n+1): row indices in column-form
C       XX (1..nz):     values in column-form.
C       Info (1):       error flag
C       Info (3):       invalid entries
C       Info (2):       duplicate entries
C       Info (5):       remaining valid entries
C       Info (6):       remaining valid entries
C       Info (7):       0
C       Wp (1..n+1)     column pointers for column form
C       II (1..n+1)     column pointers for column form

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  MA38AD, MA38BD
C       subroutines called:     MA38ND, MA38ZD
C       functions called:       MAX
        INTRINSIC MAX
        EXTERNAL MA38ND, MA38ZD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER ROW, COL, PDEST, P, NZ1, PCOL, IP, XP, IO, PRL, NINVLD,
     $          NDUPL, I
        LOGICAL PR3

C  row:     row index
C  col:     column index
C  pdest:   location of an entry in the column-form, for dupl. removal
C  p:       pointer
C  nz1:     number of entries after removing invalid or duplic. entries
C  pcol:    column col starts here after duplicates removed
C  ip:      column-form copy of matrix placed in II (ip...ip+nz-1)
C  xp:      column-form copy of matrix placed in XX (xp...xp+nz-1)
C  ninvld:  number of invalid entries
C  ndupl:   number of duplicate entries
C  i:       a row index if transa is true, a column index otherwise
C  io:      I/O unit for warning messages (for invalid or dupl. entries)
C  prl:     printing level
C  pr3:     true if printing invalid and duplicate entries

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  get arguments and check memory sizes
C-----------------------------------------------------------------------

        IO = ICNTL (2)
        PRL = ICNTL (3)
        PR3 = PRL .GE. 3 .AND. IO .GE. 0

C-----------------------------------------------------------------------
C  count nonzeros in columns and check for invalid entries
C-----------------------------------------------------------------------

        NINVLD = 0
        NDUPL = 0
        DO 10 COL = 1, N
           W (COL) = 0
10      CONTINUE
        NZ1 = NZ
        DO 20 P = NZ, 1, -1
           ROW = II (P)
           COL = II (NZ+P)
           IF (ROW.LT.1.OR.ROW.GT.N.OR.COL.LT.1.OR.COL.GT.N) THEN
C             this entry is invalid - delete it
              IF (PR3) THEN
C                print the offending entry on the diagnostic I/O unit
                 CALL MA38ZD (WHO, 99, ROW, COL, XX(P), IO)
              ENDIF
              II (P)    = II (NZ1)
              II (NZ+P) = II (NZ+NZ1)
              XX (P)    = XX (NZ1)
              NZ1 = NZ1 - 1
           ELSE
              IF (TRANSA) THEN
C                factorizing A transpose
                 W (ROW) = W (ROW) + 1
              ELSE
C                factorizing A
                 W (COL) = W (COL) + 1
              ENDIF
           ENDIF
20      CONTINUE
        NINVLD = NZ - NZ1
        IF (NINVLD .NE. 0) THEN
C          invalid entries found - set warning flag and continue
           CALL MA38ND (WHO, ICNTL, INFO, 1, NINVLD)
        ENDIF

C-----------------------------------------------------------------------
C  convert triplet form to column-form
C-----------------------------------------------------------------------

        WP (1) = 1
        DO 30 I = 1, N
           WP (I+1) = WP (I) + W (I)
30      CONTINUE
        DO 40 I = 1, N
           W (I) = WP (I)
40      CONTINUE

C       ----------------------------------------------------------------
C       construct column-form in II (2*nz+1..3*nz) and XX (nz+1..2*nz)
C       ----------------------------------------------------------------

        IP = MAX (2*NZ, N+1)
        XP = NZ
        IF (TRANSA) THEN
           DO 50 P = 1, NZ1
              ROW = II (P)
              COL = II (NZ+P)
              II (IP + W (ROW)) = COL
              XX (XP + W (ROW)) = XX (P)
              W (ROW) = W (ROW) + 1
50         CONTINUE
        ELSE
           DO 60 P = 1, NZ1
              ROW = II (P)
              COL = II (NZ+P)
              II (IP + W (COL)) = ROW
              XX (XP + W (COL)) = XX (P)
              W (COL) = W (COL) + 1
60         CONTINUE
        ENDIF

C       ----------------------------------------------------------------
C       shift the matrix back to II (n+2..nz+n+1) and XX (n+2..nz+n+1)
C       ----------------------------------------------------------------

        NZ = NZ1
CFPP$ NODEPCHK L
        DO 70 P = 1, NZ
           II (N+1+P) = II (IP+P)
           XX (P) = XX (XP+P)
70      CONTINUE

C-----------------------------------------------------------------------
C  remove duplicate entries by adding them up
C-----------------------------------------------------------------------

        DO 80 ROW = 1, N
           W (ROW) = 0
80      CONTINUE
        PDEST = 1
        DO 100 COL = 1, N
           PCOL = PDEST
           DO 90 P = WP (COL), WP (COL+1)-1
              ROW = II (N+1+P)
              IF (W (ROW) .GE. PCOL) THEN
C                this is a duplicate entry
                 XX (W (ROW)) = XX (W (ROW)) + XX (P)
                 IF (PR3) THEN
C                   print the duplicate entry on the diagnostic I/O
C                   unit.  The row and column indices printed reflect
C                   the input matrix.
                    IF (TRANSA) THEN
                       CALL MA38ZD (WHO, 98, COL, ROW, XX (P), IO)
                    ELSE
                       CALL MA38ZD (WHO, 98, ROW, COL, XX (P), IO)
                    ENDIF
                 ENDIF
              ELSE
C                this is a new entry, store and record where it is
                 W (ROW) = PDEST
                 IF (PDEST .NE. P) THEN
                    II (N+1+PDEST) = ROW
                    XX (PDEST) = XX (P)
                 ENDIF
                 PDEST = PDEST + 1
              ENDIF
90         CONTINUE
           WP (COL) = PCOL
100     CONTINUE
        WP (N+1) = PDEST
        NZ1 = PDEST - 1
        NDUPL = NZ - NZ1
        IF (NDUPL .NE. 0) THEN
C          duplicate entries found - set warning flag and continue
           CALL MA38ND (WHO, ICNTL, INFO, 2, NDUPL)
        ENDIF
        NZ = NZ1

C-----------------------------------------------------------------------
C  save column pointers in II (1..n+1)
C-----------------------------------------------------------------------

        DO 110 COL = 1, N+1
           II (COL) = WP (COL)
110     CONTINUE

        INFO (2) = NDUPL
        INFO (3) = NINVLD
        INFO (5) = NZ
        INFO (6) = NZ
        INFO (7) = 0
        IF (NZ .EQ. 0) THEN
C          set error flag if all entries are invalid
           CALL MA38ND (WHO, ICNTL, INFO, -2, -1)
        ENDIF
        RETURN
        END

        SUBROUTINE MA38LD (NLU, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)

C=== MA38LD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  solves Lx = b, where L is the lower triangular factor of a matrix
C  (if BTF not used) or a single diagonal block (if BTF is used).
C  B is overwritten with the solution X.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       nlu:            number of LU arrowheads in the LU factors
C       npiv:           number of pivots found (normally n)
C       n:              order of matrix
C       LUp (1..nlu):   pointer to LU arrowheads in LUi
C       LUi ( ... ):    integer values of LU arrowheads
C       LUx ( ... ):    real values of LU arroheads
C       X (1..n):       the right-hand-side

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       the solution to Lx=b

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38JD
C       subroutines called:     DTRSV, DGEMV
        EXTERNAL DTRSV, DGEMV

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, K, S, LUIP, LUXP, LUK, LUDEGC, LUCP, LXP, ROW
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)

C  s:       an element, or LU arrowhead
C  k:       kth pivot
C  i:       ith row in L2 array in element s
C  luip:    integer part of s is in LUi (luip...)
C  luxp:    real part of s is in LUx (luxp...)
C  luk:     number of pivots in s
C  ludegc:  column degree of non-pivotal part of s
C  lucp:    pattern of column of s in LUi (lucp...lucp+ludegc-1)
C  lxp:     the ludegc-by-luk L2 block of s is in LUx (lxp...)
C  row:     row index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        K = 0
        DO 40 S = 1, NLU

C          -------------------------------------------------------------
C          get the s-th LU arrowhead (s = 1..nlu, in pivotal order)
C          -------------------------------------------------------------

           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LXP    = LUXP + LUK

           IF (LUK .EQ. 1) THEN

C             ----------------------------------------------------------
C             only one pivot, stride-1 sparse saxpy
C             ----------------------------------------------------------

              K = K + 1
C             L (k,k) is one
CFPP$ NODEPCHK L
              DO 10 I = 1, LUDEGC
                 ROW = LUI (LUCP+I-1)
C                col: k, L (row,col): LUx (lxp+i-1)
                 X (ROW) = X (ROW) - LUX (LXP+I-1) * X (K)
10            CONTINUE

           ELSE

C             ----------------------------------------------------------
C             more than one pivot
C             ----------------------------------------------------------

              CALL DTRSV ('L', 'N', 'U', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)
              DO 20 I = 1, LUDEGC
                 ROW = LUI (LUCP+I-1)
                 W (I) = X (ROW)
20            CONTINUE
              CALL DGEMV ('N', LUDEGC, LUK, -ONE,
     $           LUX (LXP), LUDEGC + LUK, X (K+1), 1, ONE, W, 1)
              DO 30 I = 1, LUDEGC
                 ROW = LUI (LUCP+I-1)
                 X (ROW) = W (I)
30            CONTINUE
              K = K + LUK
           ENDIF
40      CONTINUE
        RETURN
        END

        SUBROUTINE MA38MD (W, N, RPERM, CPERM, NZOFF,
     *          OFFP, OFFI, OFFX, PR,
     *          ICNTL, MP, MI, MX, MN, MNZ, PRESRV, NBLKS, BLKP,
     *          ONZ, WHO, NBELOW)
        INTEGER N, NZOFF, W(N+1), RPERM(N), CPERM(N), ONZ,
     *          OFFP(N+1), OFFI(ONZ), PR(N), ICNTL(20), MN, MNZ,
     *          MP(MN+1), MI(MNZ), NBLKS, BLKP(NBLKS+1), WHO, NBELOW
        LOGICAL PRESRV
        DOUBLE PRECISION OFFX(ONZ), MX(MNZ)

C=== MA38MD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Permute the off-diagonal blocks according to final pivot permutation.
C  This routine is called only if the block-triangular-form (BTF) is
C  used.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       Rperm (1..n):   the final row permutations, including BTF
C                       If i is the k-th pivot row, then Rperm (k) = i
C       Cperm (1..n):   the final column permutations, including BTF
C                       If j is the k-th pivot col, then Cperm (k) = j
C       Icntl:          integer control parameters, see MA38ID
C       Info:           integer informational parameters
C       who:            who called (1: MA38AD, 2: MA38BD)
C
C       if presrv is true then
C           mn:                 order of preserved matrix
C           mnz:                number of entries in preserved matrix
C           Mp (1..mn+1):       column pointers of preserved matrix
C           Mi (1..mnz):        row indices of preserved matrix
C           Mx (1..mnz):        values of preserved matrix
C           Blkp (1..nblks+1):  the index range of the blocks
C           nblks:              the number of diagonal blocks
C       else
C           mn:                 0
C           mnz:                nzoff
C           Mp:                 unaccessed
C           Offp (1..n+1):      column pointers for off-diagonal entries
C                               in original order
C           Mi (1..mnz):        the row indices of off-diagonal entries,
C                               in original order
C           Mx (1..mnz):        the values of off-diagonal entries,
C                               in original order
C           nblks:              0
C           Blkp (1..nblks+1):  unaccessed
C           nzoff:              number of entries in off-diagonal blocks

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       Offp (1..n+1):          row pointers for off-diagonal part
C       Offi (1..nzoff):        column indices in off-diagonal part
C       Offx (1..nzoff):        values in off-diagonal part
C       nzoff:                  number of entries in off-diagonal blocks
C       Pr (1..n):              inverse row permutation
C       nbelow:                 entries that are below the diagonal
C                               blocks (can only occur if who = 2)

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38DD, MA38QD, MA38PD
C       subroutines called:     MA38ZD
        EXTERNAL MA38ZD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER ROW, COL, P, BLK, K, K1, K2, IO, PRL
        LOGICAL PR3

C  row:     row index
C  col:     column index
C  p:       pointer
C  blk:     current diagonal block
C  k:       kth pivot
C  k1,k2:   current diaogonal block is A (k1..k2, k1..k2)
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C  pr3:     true if printing entries below diagonal blocks (MA38BD)

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)
        PR3 = PRL .GE. 3 .AND. IO .GE. 0

C-----------------------------------------------------------------------
C  compute inverse row permutation
C-----------------------------------------------------------------------

C       if original row i is the kth pivot row, then
C               Rperm (k) = i
C               Pr (i) = k
C       if original col j is the kth pivot col, then
C               Cperm (k) = j
CFPP$ NODEPCHK L
        DO 10 K = 1, N
           PR (RPERM (K)) = K
10      CONTINUE

C-----------------------------------------------------------------------
C  construct row-oriented pointers for permuted row-form
C-----------------------------------------------------------------------

        W (1) = 1
        DO 20 ROW = 2, N
           W (ROW) = 0
20      CONTINUE
        NBELOW = 0
        IF (PRESRV) THEN
           DO 50 BLK = 1, NBLKS
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              DO 40 COL = K1, K2
CFPP$ NODEPCHK L
                 DO 30 P = MP (CPERM (COL)), MP (CPERM (COL)+1)-1
                    ROW = PR (MI (P))
                    IF (ROW .LT. K1) THEN
C                      offdiagonal entry
                       W (ROW) = W (ROW) + 1
                    ELSE IF (ROW .GT. K2 .AND. WHO .EQ. 2) THEN
C                      This entry is below the diagonal block - invalid.
C                      This can only occur if who = 2 (MA38BD).
                       IF (PR3) THEN
C                         print the original row and column indices:
                          CALL MA38ZD (2, 96, MI (P), COL, MX (P), IO)
                       ENDIF
                       NBELOW = NBELOW + 1
                    ENDIF
30               CONTINUE
40            CONTINUE
50         CONTINUE
        ELSE
           DO 70 COL = 1, N
CFPP$ NODEPCHK L
              DO 60 P = OFFP (COL), OFFP (COL+1) - 1
                 ROW = PR (MI (P))
                 W (ROW) = W (ROW) + 1
60            CONTINUE
70         CONTINUE
        ENDIF
        DO 80 ROW = 2, N
           W (ROW) = W (ROW) + W (ROW-1)
80      CONTINUE
        W (N+1) = W (N)
C       W (row) now points just past end of row in Offi/x

C-----------------------------------------------------------------------
C  construct the row-oriented form of the off-diagonal values,
C  in the final pivot order.  The column indices in each row
C  are placed in ascending order (the access of Offi/Offx later on
C  does not require this, but it makes access more efficient).
C-----------------------------------------------------------------------

        IF (PRESRV) THEN
           DO 110 BLK = NBLKS, 1, -1
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              DO 100 COL = K2, K1, - 1
CFPP$ NODEPCHK L
                 DO 90 P = MP (CPERM (COL)), MP (CPERM (COL)+1)-1
                    ROW = PR (MI (P))
                    IF (ROW .LT. K1) THEN
C                      offdiagonal entry
                       W (ROW) = W (ROW) - 1
                       OFFI (W (ROW)) = COL
                       OFFX (W (ROW)) = MX (P)
                    ENDIF
90               CONTINUE
100           CONTINUE
110        CONTINUE
        ELSE
           DO 130 COL = N, 1, -1
CFPP$ NODEPCHK L
              DO 120 P = OFFP (CPERM (COL)), OFFP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 W (ROW) = W (ROW) - 1
                 OFFI (W (ROW)) = COL
                 OFFX (W (ROW)) = MX (P)
120           CONTINUE
130        CONTINUE
        ENDIF

C-----------------------------------------------------------------------
C  save the new row pointers
C-----------------------------------------------------------------------

        DO 140 ROW = 1, N+1
           OFFP (ROW) = W (ROW)
140     CONTINUE

        NZOFF = OFFP (N+1) - 1

        RETURN
        END

        SUBROUTINE MA38ND (WHO, ICNTL, INFO, ERROR, S)
        INTEGER WHO, ICNTL(20), INFO(40), ERROR, S

C=== MA38ND ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Print error and warning messages, and set error flags.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who             which user-callable routine called:
C                       1: MA38AD, 2: MA38BD, 3: MA38CD
C       Icntl (1):      I/O unit for error and warning messages
C       Icntl (3):      printing level
C       Info (1):       the error/warning status
C       error:          the applicable error (<0) or warning (>0).
C                       See MA38ZD for a description.
C       s:              the relevant offending value

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       Info (1):       the error/warning status

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  MA38KD, MA38AD, MA38DD, MA38ED, MA38FD,
C                               MA38BD, MA38PD, MA38RD, MA38CD, MA38JD
C       subroutines called:     MA38ZD
C       functions called:       MOD
        INTRINSIC MOD
        EXTERNAL MA38ZD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        LOGICAL BOTH
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)
        INTEGER IOERR, PRL

C  ioerr:   I/O unit for error and warning messages
C  prl:     printing level
C  both:    if true, then combine errors -3 and -4 into error -5

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IOERR = ICNTL (1)
        PRL = ICNTL (3)
        IF (ERROR .LT. 0) THEN
C          this is an error message
           BOTH = (INFO (1) .EQ. -3 .AND. ERROR .EQ. -4) .OR.
     $            (INFO (1) .EQ. -4 .AND. ERROR .EQ. -3)
           IF (BOTH) THEN
C             combine error -3 (out of integer memory) and error -4
C             (out of real memory)
              INFO (1) = -5
           ELSE
              INFO (1) = ERROR
           ENDIF
           IF (PRL .GE. 1) THEN
              CALL MA38ZD (WHO, ERROR, S, 0, ZERO, IOERR)
           ENDIF
        ELSE IF (ERROR .GT. 0) THEN
C          this is a warning message
           IF (INFO (1) .GE. 0) THEN
C             do not override a prior error setting, sum up warnings
              IF (MOD (INFO (1) / ERROR, 2) .EQ. 0) THEN
                 INFO (1) = INFO (1) + ERROR
              ENDIF
           ENDIF
           IF (PRL .GE. 2) THEN
              CALL MA38ZD (WHO, ERROR, S, 0, ZERO, IOERR)
           ENDIF
        ENDIF
        RETURN
        END

        SUBROUTINE MA38OD (XX, XSIZE, XHEAD, XUSE,
     *          LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     *          FFXP, FFSIZE, PFREE, XFREE)
        INTEGER LUI(*), NLU, FRDIMC(NLU+2), FRXP(NLU+2),
     *          FRNEXT(NLU+2), FRPREV(NLU+2), LUP(NLU),
     *          XSIZE, XUSE, XHEAD, FFXP, FFSIZE,
     *          PFREE, XFREE
        DOUBLE PRECISION XX(XSIZE)

C=== MA38OD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Garbage collection for MA38RD.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       XX:             real workspace, containing matrix being
C                       factorized and partially-computed LU factors
C       xsize:          size of XX
C       xhead:          XX (1..xhead) is in use (matrix, frontal mtc's)
C       xtail:          XX (xtail..xsize) is in use (LU factors)
C       xuse:           memory usage in Value
C       Icntl:          integer control parameters, see MA38ID
C       ffxp:           pointer to current contribution block
C       ffsize:         size of current contribution block
C       nlu:            number of LU arrowheads
C
C       FRdimc (1..nlu+2)       leading dimension of frontal matrices
C       FRxp (1..nlu+2)         pointer to frontal matrices in XX
C       FRnext (1..nlu+2)       pointer to next block in XX
C       FRprev (1..nlu+2)       pointer to previous block in XX
C       LUp (1..nlu)            pointer to LU arrowhead patters in LUi
C       LUi (*)                 pattern of LU factors

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       XX:             external fragmentation is removed at head
C       xhead:          XX (1..xhead) is in use, reduced in size
C       xuse:           memory usage in Value, reduced
C       pfree:          pointer to free block in memory list, set to 0
C       xfree:          size of free block in XX, set to -1
C       FRdimc          arrays for frontal matrices are compressed
C       FRxp            frontal matrices have been shifted
C       ffxp            current working array has been shifted

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38RD
C       functions called:       ABS
        INTRINSIC ABS

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER XDP, I, E, FDIMC, LUDEGR, LUDEGC, J, FLUIP, FXP,
     $          MHEAD, MTAIL

C  xdp:     real destination pointer, current block moved to XX (xdp...)
C  e:       an element
C  fdimc:   column dimension (number of rows) of a frontal matrix
C  ludegr:  row degree (number of columns) of a contribution block
C  ludegc:  column degree (number of rows) of a contribution block
C  fluip:   element is in LUi (fluip...)
C  fxp:     element is in XX (fxp...) prior to compression
C  mhead:   nlu+1, head pointer for contribution block link list
C  mtail:   nlu+2, tail pointer for contribution block link list
C  i:       general loop index
C  j:       general loop index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  scan the link list and compress the reals
C-----------------------------------------------------------------------

        MHEAD = NLU+1
        MTAIL = NLU+2
        XDP = FRXP (MHEAD)
        E = FRNEXT (MHEAD)

C       while (e .ne. mtail) do
10      CONTINUE
        IF (E .NE. MTAIL) THEN

           FDIMC = FRDIMC (E)

C          -------------------------------------------------------------
           IF (FDIMC .EQ. 0) THEN
C          -------------------------------------------------------------

C             this is a real hole - delete it from the link list

              FRNEXT (FRPREV (E)) = FRNEXT (E)
              FRPREV (FRNEXT (E)) = FRPREV (E)

C          -------------------------------------------------------------
           ELSE
C          -------------------------------------------------------------

C             this is an unassembled frontal matrix
              FXP = FRXP (E)
              FRXP (E) = XDP
              FLUIP = LUP (E)
              LUDEGR = ABS (LUI (FLUIP+2))
              LUDEGC = ABS (LUI (FLUIP+3))
              IF (FDIMC .EQ. LUDEGC) THEN
C                contribution block is already compressed
CFPP$ NODEPCHK L
                 DO 20 I = 0, (LUDEGR * LUDEGC) - 1
                    XX (XDP+I) = XX (FXP+I)
20               CONTINUE
              ELSE
C                contribution block is not compressed
C                compress XX (fxp..) to XX (xdp..xdp+(ludegr*ludegc)-1)
                 DO 40 J = 0, LUDEGR - 1
CFPP$ NODEPCHK L
                    DO 30 I = 0, LUDEGC - 1
                       XX (XDP + J*LUDEGC + I) = XX (FXP + J*FDIMC + I)
30                  CONTINUE
40               CONTINUE
                 FRDIMC (E) = LUDEGC
              ENDIF
              XDP = XDP + LUDEGR*LUDEGC

           ENDIF

C          -------------------------------------------------------------
C          get the next item in the link list
C          -------------------------------------------------------------

           E = FRNEXT (E)

C       end while:
        GOTO 10
        ENDIF

        FRXP (MTAIL) = XDP
        PFREE = 0
        XFREE = -1

C       ----------------------------------------------------------------
C       shift the current working array (if it exists)
C       ----------------------------------------------------------------

        IF (FFXP .NE. 0) THEN
CFPP$ NODEPCHK L
           DO 50 I = 0, FFSIZE - 1
              XX (XDP+I) = XX (FFXP+I)
50         CONTINUE
           FFXP = XDP
           XDP = XDP + FFSIZE
        ENDIF

C-----------------------------------------------------------------------
C  deallocate the unused space
C-----------------------------------------------------------------------

        XUSE = XUSE - (XHEAD - XDP)
        XHEAD = XDP
        RETURN
        END

        SUBROUTINE MA38PD (N, NZ, CP, XX, XSIZE, II, ISIZE, XTAIL,
     *          ITAIL, IUSE, XUSE, NZOFF, NBLKS, ICNTL, INFO,
     *          RINFO, PRESRV, AP, AI, AX, AN, ANZ, LUI, LUISIZ,
     *          LUBLKP, BLKP, OFFP, ON, CPERM, RPERM, NE)
        INTEGER N, NZ, ISIZE, II(ISIZE), ICNTL(20), INFO(40),
     *          CP(N+1), XSIZE, XTAIL, ITAIL, IUSE, XUSE, AN, ANZ,
     *          AP(AN+1), AI(ANZ), LUISIZ, LUI(LUISIZ), NBLKS,
     *          LUBLKP(NBLKS), BLKP(NBLKS+1), ON, OFFP(ON+1),
     *          CPERM(N), RPERM(N), NZOFF, NE
        LOGICAL PRESRV
        DOUBLE PRECISION XX(XSIZE), RINFO(20), AX(ANZ)

C=== MA38PD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Refactorize an unsymmetric sparse matrix in column-form, optionally
C  permuting the matrix to upper block triangular form and factorizing
C  each diagonal block.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n:              order of matrix
C       nz:             entries in matrix
C       Cp (1..n+1):    column pointers of input matrix
C       presrv:         if true then preserve original matrix
C       xsize:          size of XX
C       isize:          size of II
C       iuse:           memory usage in Index on input
C       xuse:           memory usage in Value on input
C       Icntl:          integer control parameters, see MA38ID
C       Cntl:           real control parameters, see MA38ID
C
C       if presrv is true:
C           an:                 = n, order of preserved matrix
C           anz:                = anz, order of preserved matrix
C           Ap (1..an+1):       column pointers of preserved matrix
C           Ai (1..nz):         row indices of preserved matrix
C           Ax (1..nz):         values of preserved matrix
C           II:                 unused on input
C           XX:                 unused on input
C       else
C           an:                 0
C           anz:                1
C           Ap:                 unused
C           Ai:                 unused
C           Ax:                 unused
C           II (1..nz):         row indices of input matrix
C           XX (1..nz):         values of input matrix
C
C       Information from prior LU factorization:
C
C       LUisiz:                 size of LUi
C       LUi (1..LUisiz):        patterns of LU factors, excluding
C                               prior preserved matrix (if it existed)
C                               and prior off-diagonal part (if it
C                               existed)
C       Cperm (1..n):           column permutations
C       Rperm (1..n):           row permutations
C       nblks:                  number of diagonal blocks for BTF
C       if nblks > 1:
C           LUblkp (1..nblks):  pointers to each diagonal LU factors
C           Blkp (1..nblks+1):  index range of diagonal blocks

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       XX (xtail ... xsize), xtail:
C
C                       LU factors are located in XX (xtail ... xsize),
C                       including values in off-diagonal part if matrix
C                       was permuted to block triangular form.
C
C       II (itail ... isize), itail:
C
C                       the off-diagonal nonzeros, if nblks > 1
C
C       Offp (1..n+1):  row pointers for off-diagonal part, if nblks > 1
C       Info:           integer informational output, see MA38AD
C       Rinfo:          real informational output, see MA38AD

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38BD
C       subroutines called:     MA38ND, MA38RD, MA38QD, MA38MD
C       functions called:       MAX
        INTRINSIC MAX
        EXTERNAL  MA38MD, MA38ND, MA38QD, MA38RD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, NZDIA, P, IHEAD, NSGLTN, NSYM, WP, ARIP, ARXP, NPIV,
     $          WRKSIZ, NLU, PRP, MC, MR, DUMMY1, DUMMY2, NZ2, K, BLK,
     $          K1, K2, KN, NZBLK, COL, ROW, PRL, IO, LUIP, MNZ, ARNZ,
     $          XHEAD, OFFIP, OFFXP, NOUTSD, NBELOW, NZORIG, XRMAX
        DOUBLE PRECISION
     $          ZERO, ONE, A
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C
C  Allocated array pointers:
C  -------------------------
C  wp:      W (1..n+1), or W (1..kn+1), work array located in II (wp...)
C  prp:     Pr (1..n) work array located in II (prp..prp+n-1)
C  arip:    Ari (1..nz) array located in II (arip..arip+nz-1)
C  arxp:    Arx (1..nz) array located in XX (arxp..arxp+nz-1)
C  offip:   Offi (1..nzoff) array located in II (offip..offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array located in XX (offxp..offip+nzoff-1)
C
C  Arrowhead-form matrix:
C  ----------------------
C  nz2:     number of entries in arrowhead matrix
C  nzblk:   number of entries in arrowhead matrix of a single block
C  arnz:    arrowhead form of blocks located in II/XX (1..arnz)
C
C  BTF information:
C  ----------------
C  k1:      starting index of diagonal block being factorized
C  k2:      ending index of diagonal block being factorized
C  kn:      the order of the diagonal block being factorized
C  blk:     block number of diagonal block being factorized
C  nsgltn:  number of 1-by-1 diagonal blocks ("singletons")
C  a:       numerical value of a singleton
C  mnz:     nzoff
C  noutsd:  entries in diagonal blocks, but not in LU (invalid)
C  nbelow:  entries below diagonal blocks (invalid)
C  nzoff:   entries above diagonal blocks (valid)
C  nzdia:   entries in diagonal blocks (valid)
C  nzorig:  total number of original entries
C
C  Memory usage:
C  -------------
C  xhead:   XX (1..xhead-1) is in use, XX (xhead..xtail-1) is free
C  ihead:   II (1..ihead-1) is in use, II (ihead..itail-1) is free
C  wrksiz:  total size of work arrays need in II for call to MA38RD
C  xrmax:   memory needed in Value for next call to MA38BD
C
C  Symbolic information and pattern of prior LU factors:
C  -----------------------------------------------------
C  nlu:     number of elements in a diagonal block
C  luip:    integer part of LU factors located in LUi (luip...)
C  mr,mc:   largest frontal matrix for this diagonal block is mc-by-mr
C
C  Other:
C  ------
C  k:       loop index (kth pivot)
C  i:       loop index
C  row:     row index
C  col:     column index
C  p:       pointer
C  nsym:    number of symmetric pivots chosen
C  dummy1:  argument returned by MA38QD, but not needed
C  dummy2:  argument returned by MA38QD, but not needed

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)
        NZORIG = NZ

        IF (PRESRV) THEN
C          original matrix is not in Cp/II/XX, but in Ap/Ai/Ax:
           IHEAD = 1
           XHEAD = 1
        ELSE
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF

        NZOFF = 0
        NZDIA = NZ
        NSGLTN = 0
        NPIV = 0
        NOUTSD = 0
        NBELOW = 0
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1

C-----------------------------------------------------------------------
C current memory usage:
C-----------------------------------------------------------------------

C       if .not. presrv then
C               input matrix is now in II (1..nz) and XX (1..nz)
C                       col pattern: II (Cp (col) ... Cp (col+1))
C                       col values:  XX (Cp (col) ... Cp (col+1))
C               total: nz+n+1 integers, nz reals
C       else
C               input matrix is now in Ai (1..nz) and Ax (1..nz)
C                       col pattern: Ai (Ap (col) ... Ap (col+1))
C                       col values:  Ax (Ap (col) ... Ap (col+1))
C               Cp is a size n+1 integer workspace
C               total: nz+2*(n+1) integers, nz reals
C
C       if (nblks > 1) then
C                       LUblkp (1..nblks)
C                       Blkp (1..nblks+1)
C                       Offp (1..n+1)
C               total: (2*nblks+n+2) integers
C
C       Cperm (1..n) and Rperm (1..n)
C               total: 2*n integers
C
C   Grand total current memory usage (including II,XX,Cp,Ai,Ap,Ax
C       and LUi):
C
C       presrv  nblks>1 integers, iuse =
C       F       F       LUisiz + nz+  (n+1)+(2*n+7)
C       F       T       LUisiz + nz+  (n+1)+(2*n+7)+(2*nblks+n+2)
C       T       F       LUisiz + nz+2*(n+1)+(2*n+7)
C       T       T       LUisiz + nz+2*(n+1)+(2*n+7)+(2*nblks+n+2)
C
C   real usage is xuse = nz

C-----------------------------------------------------------------------
C  get memory usage estimate for next call to MA38BD
C-----------------------------------------------------------------------

        XRMAX = 2*NE

C-----------------------------------------------------------------------
C  convert matrix into arrowhead format (unless BTF and preserved)
C-----------------------------------------------------------------------

        IF (NBLKS .GT. 1 .AND. PRESRV) THEN

C          -------------------------------------------------------------
C          BTF is to be used, and original matrix is to be preserved.
C          It is converted and factorized on a block-by-block basis,
C          using the inverse row permutation (computed and stored in
C          Offp (1..n)
C          -------------------------------------------------------------

           DO 10 K = 1, N
              OFFP (RPERM (K)) = K
10         CONTINUE

        ELSE

C          -------------------------------------------------------------
C          convert the entire input matrix to arrowhead form
C          -------------------------------------------------------------

C          -------------------------------------------------------------
C          allocate workspace: W (n+1), Pr (n), Ari (nz), Arx (nz)
C          -------------------------------------------------------------

           ITAIL = ITAIL - (2*N+1)
           IUSE = IUSE + 2*N+1
           PRP = ITAIL
           WP = PRP + N
           IUSE = IUSE + NZ
           XUSE = XUSE + NZ
           ARXP = XHEAD
           ARIP = IHEAD
           IHEAD = IHEAD + NZ
           XHEAD = XHEAD + NZ
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XUSE)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
C             error return, if not enough integer and/or real memory:
              GO TO 9000
           ENDIF

C          -------------------------------------------------------------
C          convert
C          -------------------------------------------------------------

           IF (NBLKS .EQ. 1) THEN
              IF (PRESRV) THEN
                 CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $              II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $              ICNTL, AP, BLKP, AI, AX, OFFP, ON, NZ,
     $              0, N, NZ2, I)
              ELSE
                 CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $              II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $              ICNTL, CP, BLKP, II, XX, OFFP, ON, NZ,
     $              0, N, NZ2, I)
              ENDIF
           ELSE
C             note that presrv is false in this case
              CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $           II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $           ICNTL, CP, BLKP, II, XX, OFFP, ON, NZ,
     $           0, N, NZ2, NBELOW)
           ENDIF

C          -------------------------------------------------------------
C          copy the arrowhead pointers from W (1..n+1) to Cp (1..n+1)
C          -------------------------------------------------------------

           DO 20 I = 1, N+1
              CP (I) = II (WP+I-1)
20         CONTINUE

C          -------------------------------------------------------------
C          deallocate W and Pr.  If not presrv deallocate Ari and Arx
C          -------------------------------------------------------------

           IUSE = IUSE - (2*N+1)
           IF (.NOT. PRESRV) THEN
C             Ari and Arx have been deallocated.
              XUSE = XUSE - NZ
              IUSE = IUSE - NZ
           ENDIF
           ITAIL = ISIZE + 1
           XTAIL = XSIZE + 1
           NZ = NZ2
           IHEAD = NZ + 1
           XHEAD = NZ + 1

        ENDIF

        INFO (5) = NZ
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (4) = NBELOW

C-----------------------------------------------------------------------
C  refactorization
C-----------------------------------------------------------------------

C       ----------------------------------------------------------------
C       if nblks=1
C          arrowhead form is now stored in II (1..nz) and XX (1..nz)
C          in reverse pivotal order (arrowhead n, n-1, ..., 2, 1).
C          The arrowhead form will be overwritten.
C       else if not presrv
C          off-diagonal part is in II (1..nzoff), XX (1..nzoff),
C          (with row pointers Offp (1..n+1)) followed by each diagonal
C          block (block 1, 2, ... nblks) in II/XX (nzoff+1...nz).
C          Each diagonal block is in arrowhead form, and in
C          reverse pivotal order (arrowhead k2, k2-1, ..., k1-1, k1).
C          The arrowhead form will be overwritten.
C       else (nblks > 1 and presrv)
C          II and XX are still empty.  Original matrix is in Ap, Ai,
C          and Ax.  Inverse row permutation (Pr) is in Offp (1..n).
C          The arrowhead form is not yet computed.
C       ----------------------------------------------------------------

        IF (NBLKS .EQ. 1) THEN

C          -------------------------------------------------------------
C          refactorize the matrix as a single block
C          -------------------------------------------------------------

           NLU = LUI (2)
           MC = LUI (4)
           MR = LUI (5)
           WRKSIZ = 2*N + MR + 3*MC + 4*(NLU+2)
           ITAIL = ITAIL - WRKSIZ
           IUSE = IUSE + WRKSIZ
           P = ITAIL
           INFO (18) = MAX (INFO (18), IUSE)
           IF (IHEAD .GT. ITAIL) THEN
C             error return, if not enough integer memory:
              GO TO 9000
           ENDIF

           CALL MA38RD (CP, NZ, N, XTAIL,
     $          XX, XSIZE, XUSE, II, CPERM, RPERM,
     $          ICNTL, INFO, RINFO, MC, MR,
     $          II (P), II (P+N), II (P+2*N), II (P+2*N+MR),
     $          II (P+2*N+MR+MC), II (P+2*N+MR+2*MC),
     $          II (P+2*N+MR+3*MC), II (P+2*N+MR+3*MC+(NLU+2)),
     $          II (P+2*N+MR+3*MC+2*(NLU+2)),
     $          II (P+2*N+MR+3*MC+3*(NLU+2)),
     $          NLU, LUI (6), LUI (NLU+6), NOUTSD,
     $          XRMAX)

           IF (INFO (1) .LT. 0) THEN
C             error return, if not enough real memory or bad pivot found
              GO TO 9010
           ENDIF

C          -------------------------------------------------------------
C          deallocate workspace and original matrix (reals already done)
C          -------------------------------------------------------------

           IUSE = IUSE - WRKSIZ - NZ
           ITAIL = ITAIL + WRKSIZ
           LUI (1) = 1
           IHEAD = 1
           XHEAD = 1

        ELSE

C          -------------------------------------------------------------
C          refactorize the block-upper-triangular form of the matrix
C          -------------------------------------------------------------

           IF (PRESRV) THEN
C             count the entries in off-diagonal part
              NZOFF = 0
           ENDIF

           DO 70 BLK = NBLKS, 1, -1

C             ----------------------------------------------------------
C             factorize the kn-by-kn block, A (k1..k2, k1..k2)
C             ----------------------------------------------------------

C             get k1 and k2, the start and end of this block
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              KN = K2-K1+1
              A = ZERO

C             ----------------------------------------------------------
C             get pointers to, or place the block in, arrowhead form
C             ----------------------------------------------------------

              IF (PRESRV) THEN

                 IF (KN .GT. 1) THEN

C                   ----------------------------------------------------
C                   convert a single block to arrowhead format, using
C                   the inverse row permutation stored in Offp
C                   ----------------------------------------------------

C                   ----------------------------------------------------
C                   compute nzblk, allocate II/XX (1..nzblk), W(1..kn+1)
C                   ----------------------------------------------------

                    NZBLK = 0
                    DO 40 K = K1, K2
                       COL = CPERM (K)
CFPP$ NODEPCHK L
                       DO 30 P = AP (COL), AP (COL+1) - 1
                          ROW = OFFP (AI (P))
                          IF (ROW .LT. K1) THEN
C                            entry in off-diagonal part
                             NZOFF = NZOFF + 1
                          ELSE IF (ROW .LE. K2) THEN
C                            entry in diagonal block
                             NZBLK = NZBLK + 1
                          ENDIF
30                     CONTINUE
40                  CONTINUE

                    ITAIL = ITAIL - (KN+1)
                    WP = ITAIL
                    IHEAD = NZBLK + 1
                    XHEAD = NZBLK + 1
                    IUSE = IUSE + NZBLK + KN+1
                    XUSE = XUSE + NZBLK
                    XRMAX = MAX (XRMAX, XUSE)
                    INFO (18) = MAX (INFO (18), IUSE)
                    INFO (20) = MAX (INFO (20), XUSE)
                    INFO (21) = MAX (INFO (21), XUSE)
                    IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
C                      error return, if not enough integer
C                      and/or real memory:
                       GO TO 9000
                    ENDIF

C                   ----------------------------------------------------
C                   convert blk from column-form in Ai/Ax to arrowhead
C                   form in II/XX (1..nzblk)
C                   ----------------------------------------------------

                    CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, OFFP,
     $                 II (WP), NBLKS, XX, II, DUMMY1, DUMMY2,
     $                 ICNTL, AP, BLKP, AI, AX, OFFP, 0, NZBLK,
     $                 BLK, KN, NZ2, I)

C                   ----------------------------------------------------
C                   copy the arrowhead pointers from W (1..kn+1)
C                   to Cp (k1 ... k2+1)
C                   ----------------------------------------------------

                    DO 50 I = 0, KN
                       CP (K1+I) = II (WP+I)
50                  CONTINUE

C                   Cp (k1) is nzblk + 1 and Cp (k2+1) is 1

C                   ----------------------------------------------------
C                   deallocate W (1..kn+1)
C                   ----------------------------------------------------

                    IUSE = IUSE - (KN+1)
                    ITAIL = ITAIL + (KN+1)

                 ELSE

C                   ----------------------------------------------------
C                   get the value of singleton at A (k1,k1) if it exists
C                   ----------------------------------------------------

C                   find the diagonal entry in the unpermuted matrix,
C                   and count the entries in the diagonal and
C                   off-diagonal blocks.
                    COL = CPERM (K1)
                    DO 60 P = AP (COL), AP (COL + 1) - 1
C                      inverse row permutation is stored in Offp
                       ROW = OFFP (AI (P))
                       IF (ROW .LT. K1) THEN
                          NZOFF = NZOFF + 1
                       ELSE IF (ROW .EQ. K1) THEN
                          A = AX (P)
C                      else
C                         this is an invalid entry, below the diagonal
C                         block.  It will be detected (and optionally
C                         printed) in the call to MA38MD below.
                       ENDIF
60                  CONTINUE

                    IHEAD = 1
                    XHEAD = 1
                 ENDIF

              ELSE

C                -------------------------------------------------------
C                the block is located in II/XX (Cp (k2+1) ... Cp (k1)-1)
C                and has already been converted to arrowhead form
C                -------------------------------------------------------

                 IF (BLK .EQ. 1) THEN
C                   this is the last block to factorize
                    CP (K2+1) = NZOFF + 1
                 ELSE
                    CP (K2+1) = CP (BLKP (BLK-1))
                 ENDIF

                 IHEAD = CP (K1)
                 XHEAD = IHEAD

                 IF (KN .EQ. 1) THEN
C                   singleton block in II/XX (Cp (k1+1) ... Cp (k1)-1)
                    IF (CP (K1) .GT. CP (K1+1)) THEN
                       A = XX (CP (K1) - 1)
                       IHEAD = IHEAD - 1
                       XHEAD = XHEAD - 1
                       IUSE = IUSE - 1
                       XUSE = XUSE - 1
                    ENDIF
                 ENDIF

              ENDIF

C             ----------------------------------------------------------
C             factor the block
C             ----------------------------------------------------------

              IF (KN .GT. 1) THEN

C                -------------------------------------------------------
C                The A (k1..k2, k1..k2) block is not a singleton.
C                block is now in II/XX (Cp (k2+1) ... Cp (k1)-1), in
C                arrowhead form, and is to be overwritten with LU
C                -------------------------------------------------------

                 ARNZ = CP (K1) - 1

C                if (presrv) then
C                   II/XX (1..arnz) holds just the current block, blk
C                else
C                   II/XX (1..arnz) holds the off-diagonal part, and
C                   blocks 1..blk, in that order.
C                endif

                 LUIP = LUBLKP (BLK)
C                luxp = LUi (luip), not needed for refactorization
                 NLU = LUI (LUIP+1)
C                npiv = LUi (luip+2), not needed for refactorization
                 MC = LUI (LUIP+3)
                 MR = LUI (LUIP+4)
                 WRKSIZ = 2*KN + MR + 3*MC + 4*(NLU+2)
                 ITAIL = ITAIL - WRKSIZ
                 IUSE = IUSE + WRKSIZ
                 P = ITAIL
                 INFO (18) = MAX (INFO (18), IUSE)
                 IF (IHEAD .GT. ITAIL) THEN
C                   error return, if not enough integer memory:
                    GO TO 9000
                 ENDIF

                 CALL MA38RD (CP (K1), ARNZ, KN, XTAIL,
     $                XX, XTAIL-1, XUSE, II, CPERM (K1), RPERM (K1),
     $                ICNTL, INFO, RINFO, MC, MR,
     $                II (P), II (P+KN), II (P+2*KN), II (P+2*KN+MR),
     $                II (P+2*KN+MR+MC), II (P+2*KN+MR+2*MC),
     $                II (P+2*KN+MR+3*MC), II (P+2*KN+MR+3*MC+(NLU+2)),
     $                II (P+2*KN+MR+3*MC+2*(NLU+2)),
     $                II (P+2*KN+MR+3*MC+3*(NLU+2)),
     $                NLU, LUI (LUIP+5), LUI (LUIP+NLU+5), NOUTSD,
     $                XRMAX)

                 IF (INFO (1) .LT. 0) THEN
C                   error return, if not enough real memory or bad pivot
                    GO TO 9010
                 ENDIF

C                -------------------------------------------------------
C                deallocate workspace and original matrix (reals
C                already deallocated in MA38RD)
C                -------------------------------------------------------

                 IUSE = IUSE - WRKSIZ
                 ITAIL = ITAIL + WRKSIZ
                 LUI (LUIP) = XTAIL
                 IUSE = IUSE - (IHEAD - CP (K2+1))
                 IHEAD = CP (K2+1)
                 XHEAD = IHEAD

              ELSE

C                -------------------------------------------------------
C                factor the singleton A (k1,k1) block, in a
C                -------------------------------------------------------

                 NSGLTN = NSGLTN + 1
                 IF (A .EQ. ZERO) THEN
C                   this is a singular matrix, replace with 1-by-1
C                   identity matrix.
                    A = ONE
                 ELSE
C                   increment pivot count
                    NPIV = NPIV + 1
                 ENDIF
                 XTAIL = XTAIL - 1
                 XUSE = XUSE + 1
                 XRMAX = MAX (XRMAX, XUSE)
                 INFO (20) = MAX (INFO (20), XUSE)
                 INFO (21) = MAX (INFO (21), XUSE)
C                note: if the matrix is not preserved and nonsingular
C                then we will not run out of memory
                 IF (XHEAD .GT. XTAIL) THEN
C                   error return, if not enough real memory:
                    GO TO 9000
                 ENDIF

C                -------------------------------------------------------
C                store the 1-by-1 LU factors
C                -------------------------------------------------------

                 XX (XTAIL) = A
                 LUBLKP (BLK) = -XTAIL

              ENDIF
70         CONTINUE

C          -------------------------------------------------------------
C          make the index of each block relative to start of LU factors
C          -------------------------------------------------------------
CFPP$ NODEPCHK L
           DO 80 BLK = 1, NBLKS
              IF (LUBLKP (BLK) .GT. 0) THEN
                 LUI (LUBLKP (BLK)) = LUI (LUBLKP (BLK)) - XTAIL + 1
              ELSE
C                this is a singleton
                 LUBLKP (BLK) = (-LUBLKP (BLK)) - XTAIL + 1
              ENDIF
80         CONTINUE

C          -------------------------------------------------------------
C          store the off-diagonal blocks
C          -------------------------------------------------------------

           IF (PRESRV) THEN

C             ----------------------------------------------------------
C             allocate temporary workspace for Pr (1..n) at head of II
C             ----------------------------------------------------------

              PRP = IHEAD
              IHEAD = IHEAD + N
              IUSE = IUSE + N

C             ----------------------------------------------------------
C             allocate permanent copy of off-diagonal blocks
C             ----------------------------------------------------------

              ITAIL = ITAIL - NZOFF
              OFFIP = ITAIL
              XTAIL = XTAIL - NZOFF
              OFFXP = XTAIL
              IUSE = IUSE + NZOFF
              XUSE = XUSE + NZOFF
              XRMAX = MAX (XRMAX, XUSE)
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
              IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
C                error return, if not enough integer and/or real memory:
                 GO TO 9000
              ENDIF

C             ----------------------------------------------------------
C             re-order the off-diagonal blocks according to pivot perm
C             ----------------------------------------------------------

C             use Cp as temporary work array:
              MNZ = NZOFF
              IF (NZOFF .EQ. 0) THEN
C                Offi and Offx are not accessed in MA38MD.  Set offip
C                and offxp to 1 (since offip = itail = isize+1, which
C                can generate an address fault, otherwise).
                 OFFIP = 1
                 OFFXP = 1
              ENDIF
              CALL MA38MD (CP, N, RPERM, CPERM, NZOFF,
     $             OFFP, II (OFFIP), XX (OFFXP), II (PRP),
     $             ICNTL, AP, AI, AX, AN, ANZ, PRESRV, NBLKS, BLKP,
     $             MNZ, 2, NBELOW)

C             ----------------------------------------------------------
C             deallocate Pr (1..n)
C             ----------------------------------------------------------

              IHEAD = 1
              XHEAD = 1
              IUSE = IUSE - N

           ELSE

C             off-diagonal entries are in II/XX (1..nzoff); shift down
C             to II/XX ( ... itail/xtail).  No extra memory needed.
              DO 90 I = NZOFF, 1, -1
                 II (ITAIL+I-NZOFF-1) = II (I)
                 XX (XTAIL+I-NZOFF-1) = XX (I)
90            CONTINUE
              IHEAD = 1
              XHEAD = 1
              ITAIL = ITAIL - NZOFF
              XTAIL = XTAIL - NZOFF
           ENDIF

        ENDIF

C       ----------------------------------------------------------------
C       clear the flags (negated row/col indices, and negated ludegr/c)
C       ----------------------------------------------------------------

        DO 100 I = 1, LUISIZ
           LUI (I) = ABS (LUI (I))
100     CONTINUE

C-----------------------------------------------------------------------
C  normal and error return
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE
        IF (IHEAD .GT. ITAIL) THEN
C          set error flag if not enough integer memory
           CALL MA38ND (2, ICNTL, INFO, -3, INFO (18))
        ENDIF
        IF (XHEAD .GT. XTAIL) THEN
C          set error flag if not enough real memory
           CALL MA38ND (2, ICNTL, INFO, -4, INFO (21))
        ENDIF

C       error return label, for error return from MA38RD:
9010    CONTINUE

C       Cp can now be deallocated in MA38BD:
        IUSE = IUSE - (N+1)

        INFO (4) = NOUTSD + NBELOW
        NZDIA = NZORIG - NZOFF - NOUTSD - NBELOW
        INFO (5) = NZOFF + NZDIA
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (8) = NSGLTN
        INFO (9) = NBLKS
        INFO (12) = INFO (10) + INFO (11) + N + INFO (7)

C       Count the number of symmetric pivots chosen.  Note that some
C       of these may have been numerically unacceptable.
        NSYM = 0
        DO 110 K = 1, N
           IF (CPERM (K) .EQ. RPERM (K)) THEN
C             this kth pivot came from the diagonal of A
              NSYM = NSYM + 1
           ENDIF
110     CONTINUE
        INFO (16) = NSYM

        INFO (17) = INFO (17) + NPIV
        RINFO (1) = RINFO (4) + RINFO (5) + RINFO (6)

C       set warning flag if entries outside prior pattern are present
        IF (INFO (4) .GT. 0) THEN
           CALL MA38ND (2, ICNTL, INFO, 1, -INFO (4))
        ENDIF

C       set warning flag if matrix is singular
        IF (INFO (1) .GE. 0 .AND. INFO (17) .LT. N) THEN
           CALL MA38ND (2, ICNTL, INFO, 4, INFO (17))
        ENDIF

C       ----------------------------------------------------------------
C       return memory usage estimate for next call to MA38BD
C       ----------------------------------------------------------------

        INFO (23) = XRMAX

        RETURN
        END

        SUBROUTINE MA38QD (PRESRV, N, NZ, CPERM, RPERM, PR,
     *          W, NBLKS, ARX, ARI, NZOFF, NZDIA,
     *          ICNTL, MP, BLKP, MI, MX, OFFP, ON, NZBLK,
     *          CBLK, KN, NZ2, NBELOW)
        INTEGER N, NZ, CPERM(N), RPERM(N), PR(N), KN, W(KN+1),
     *          NBLKS, NZBLK, ARI(NZBLK), NZOFF, NZDIA, MP(N+1),
     *          MI(NZ), ON, ICNTL(20), BLKP(NBLKS+1), NZ2,
     *          OFFP(ON+1), CBLK, NBELOW
        LOGICAL PRESRV
        DOUBLE PRECISION ARX(NZBLK), MX(NZ)

C=== MA38QD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Convert a column-oriented matrix into an arrowhead format.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       n               size of entire matrix
C       Mi (1..nz):     row indices of column form of entire matrix
C       Mx (1..nz):     values of column form of entire matrix
C       Mp (1..n+1)     column pointers for entire matrix
C       Cperm (1..n):   column permutations
C       Rperm (1..n):   row permutations
C
C       if nblks > 1 and presrv
C           cblk:               the block to convert
C           kn:                 the size of the block to convert
C       else
C           cblk:               0
C           kn                  n, size of input matrix

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       if nblks = 1 and not presrv
C
C           nzoff               0
C           nzdia               nz - (entries below in diagonal blocks)
C           nz2                 nzdia
C
C           Mi (1..nz2)         arrowheads for the diagonal block
C           Mx (1..nz2)
C           Ari, Arx            used as workspace
C           W (1..n+1)          pointer to each arrowhead in Mi/Mx
C
C           Offp                not accessed
C
C       if nblks = 1 and presrv
C
C           nzoff               0
C           nzdia               nz - (entries below in diagonal blocks)
C           nz2                 nzdia
C
C           Mi, Mx              not modified
C           Ari (1..nz2)        arrowheads for the diagonal block
C           Arx (1..nz2)
C           W (1..n+1)          pointer to each arrowhead in Ari/Arx
C
C           Offp                not accessed
C
C       else if nblks > 1 and not presrv
C
C           nzoff               number of entries in off-diagonal part
C           nzdia               number of entries in diagonal blocks
C                               (nz = nzoff + nzdia + entries below
C                               diagonal blocks)
C           nz2                 nzoff + nzdia
C
C           Mi (nzoff+1..nz2)   arrowheads for each diagonal block
C           Mx (nzoff+1..nz2)
C           Ari, Arx            used as workspace
C           W (1..n+1)          pointer to each arrowhead in Mi/Mx
C
C           Offp (1..n+1)       row pointers for off-diagonal part
C           Mi (1..nzoff)       col indices for off-diagonal part
C           Mx (1..nzoff)       values for off-diagonal part
C
C       else (nblks > 1 and presrv)
C
C           nzoff               0
C           nzdia               nonzeros in the diagonal block, cblk
C           nz2                 nzdia
C
C           Mi, Mx              not modified
C           Ari (1..nz2)        arrowheads for the diagonal block, cblk
C           Arx (1..nz2)
C           W (1..kn+1)         pointer to each arrowhead in Ari/Arx
C
C           Offp                not accessed

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38PD
C       subroutines called:     MA38MD
C       functions called:       MIN
        INTRINSIC MIN
        EXTERNAL MA38MD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, P, ROW, COL, BLK, BASE, K1, K2, K, B1, B2, K0

C  i:       loop index, arrowhead index
C  p:       pointer into column-form input matrix
C  row:     row index
C  col:     column index
C  blk:     current diagonal block
C  base:    where to start the construction of the arrowhead form
C  k1,k2:   current diagonal block is A (k1..k2, k1..k2)
C  k:       loop index, kth pivot
C  b1,b2:   convert blocks b1...b2 from column-form to arrowhead form
C  k0:      convert A (k0+1..., k0+1...) to arrowhead form

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  if entire matrix is to be converted, then create the off-diagonal
C  part in row-oriented form in Ari (1..nzoff) and Arx (1..nzoff) and
C  compute inverse row permutation.  Otherwise, the inverse row
C  permutation has already been computed.
C-----------------------------------------------------------------------

        NZOFF = 0
        NBELOW = 0
        IF (NBLKS .EQ. 1) THEN
           DO 10 K = 1, N
              PR (RPERM (K)) = K
10         CONTINUE
        ELSE IF (NBLKS .GT. 1 .AND. .NOT. PRESRV) THEN
           CALL MA38MD (W, N, RPERM, CPERM, NZOFF,
     $        OFFP, ARI, ARX, PR,
     $        ICNTL, MP, MI, MX, N, NZ, .TRUE., NBLKS, BLKP,
     $        NZ, 2, NBELOW)
        ENDIF

C-----------------------------------------------------------------------
C  construct the arrowhead form for the diagonal block(s)
C-----------------------------------------------------------------------

        DO 20 I = 1, KN+1
           W (I) = 0
20      CONTINUE

        BASE = NZOFF + 1

        IF (CBLK .NE. 0) THEN
C          convert just cblk
           K0 = BLKP (CBLK) - 1
           B1 = CBLK
           B2 = CBLK
        ELSE
C          convert all the block(s)
           K0 = 0
           B1 = 1
           B2 = NBLKS
        ENDIF

        DO 80 BLK = B1, B2

C          -------------------------------------------------------------
C          get the starting and ending indices of this diagonal block
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1) THEN
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
           ELSE
              K1 = 1
              K2 = N
           ENDIF

C          -------------------------------------------------------------
C          count the number of entries in each arrowhead
C          -------------------------------------------------------------

           DO 40 COL = K1, K2
              DO 30 P = MP (CPERM (COL)), MP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 IF (ROW .GE. K1 .AND. ROW .LE. K2) THEN
C                   this is in a diagonal block, arrowhead i
                    I = MIN (ROW, COL) - K0
                    W (I) = W (I) + 1
                 ENDIF
30            CONTINUE
40         CONTINUE

C          -------------------------------------------------------------
C          set pointers to point just past end of each arrowhead
C          -------------------------------------------------------------

           W (K2-K0+1) = W (K2-K0) + BASE
           DO 50 I = K2-K0, K1-K0+1, -1
              W (I) = W (I+1) + W (I-1)
50         CONTINUE
           W (K1-K0) = W (K1-K0+1)
C          W (i+1-k0) points just past end of arrowhead i in Ari/Arx

C          -------------------------------------------------------------
C          construct arrowhead form, leaving pointers in final state
C          -------------------------------------------------------------

           DO 70 COL = K1, K2
              DO 60 P = MP (CPERM (COL)), MP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 IF (ROW .GE. K1 .AND. ROW .LE. K2) THEN
                    IF (ROW .GE. COL) THEN
C                      diagonal, or lower triangular part
                       I = COL - K0 + 1
                       W (I) = W (I) - 1
                       ARI (W (I)) = ROW - K1 + 1
                       ARX (W (I)) = MX (P)
                    ELSE
C                      upper triangular part, flag by negating col
                       I = ROW - K0 + 1
                       W (I) = W (I) - 1
                       ARI (W (I)) = -(COL - K1 + 1)
                       ARX (W (I)) = MX (P)
                    ENDIF
                 ENDIF
60            CONTINUE
70         CONTINUE

           BASE = W (K1-K0)
           W (K2-K0+1) = 0
80      CONTINUE

        W (KN+1) = NZOFF + 1
        NZDIA = BASE - NZOFF - 1
        NZ2 = NZOFF + NZDIA

C       ----------------------------------------------------------------
C       if cblk = 0, the entire matrix has been converted:
C
C          W (i) now points just past end of arrowhead i in Ari/Arx
C          arrowhead i is located in Ari/Arx (W (i+1) ... W (i)-1),
C          except for the k2-th arrowhead in each block.  Those are
C          located in Ari/Arx (base ... W (k2) - 1), where base is
C          W (Blkp (blk-1)) if blk>1 or W (n+1) = nzoff + 1 otherwise.
C
C       otherwise, just one block has been converted:
C
C          W (i) now points just past end of arrowhead i in Ari/Arx,
C          where i = 1 is the first arrowhead of this block (not the
C          first arrowhead of the entire matrix).  Arrowhead i is
C          located in Ari/Arx (W (i+1) ... W (i)-1).
C          This option is used only if nblks>1 and presrv is true.
C       ----------------------------------------------------------------

C-----------------------------------------------------------------------
C  if not preserved, overwrite column-form with arrowhead form
C-----------------------------------------------------------------------

        IF (.NOT. PRESRV) THEN
           DO 90 I = 1, NZ
              MI (I) = ARI (I)
              MX (I) = ARX (I)
90         CONTINUE
        ENDIF

        RETURN
        END

        SUBROUTINE MA38RD (CP, NZ, N, XTAIL, XX, XSIZE, XUSE, ARI,
     *          CPERM, RPERM, ICNTL, INFO, RINFO, MC, MR,
     *          WIR, WIC, WPR, WPC, WM, WJ, FRDIMC, FRXP, FRNEXT,
     *          FRPREV, NLU, LUP, LUI, NOUTSD, XRMAX)
        INTEGER XSIZE, ICNTL(20), INFO(40), N, CPERM(N), RPERM(N),
     *          XTAIL, NZ, ARI(NZ), CP(N+1), MR, MC, NOUTSD,
     *          WIR(N), WIC(N), WPR(MR), XRMAX, WPC(MC), WM(MC),
     *          NLU, FRDIMC(NLU+2), FRXP(NLU+2), XUSE, WJ(MC),
     *          FRNEXT(NLU+2), FRPREV(NLU+2), LUP(NLU), LUI(*)
        DOUBLE PRECISION XX(XSIZE), RINFO(20)

C=== MA38RD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  MA38RD refactorizes the n-by-n input matrix at the head of XX
C  (in arrowhead form) and places its LU factors at the tail of
C  XX.  The input matrix is overwritten.   No BTF information is
C  used in this routine.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       Cp (1..n+1):    column pointers of arrowhead form
C       n:              order of input matrix
C       nz:             entries in input matrix
C       xsize:          size of XX
C       Icntl:          integer control parameters, see MA38ID
C       Cntl:           real control parameters, see MA38ID
C
C       Ari (1..nz):            arrowhead format of A
C       XX (1..nz):             arrowhead format of A, see below
C       XX (nz+1..xsize):       undefined on input, used as workspace
C
C       nlu:            number of LU arrowheads
C       LUp (1..nlu):   pointers to LU arrowheads in LUi
C       LUi (1.. ):     LU arrowheads
C
C       xuse:           memory usage in Value
C
C       noutsd:         entries not in prior LU pattern
C
C       Cperm (1..n):   column permutation
C       Rperm (1..n):   row permutation

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       WiR (1..n)
C       WiC (1..n)
C
C       WpR (1.. max ludegr)
C       WpC (1.. max ludegc)
C       Wm  (1.. max ludegc)
C       Wj  (1.. max ludegc)
C
C       FRdimc (1..nlu+2)
C       FRxp   (1..nlu+2)
C       FRnext (1..nlu+2)
C       FRprev (1..nlu+2)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       LUi (1..):              LU arrowheads, modified luxp pointers
C       XX (1..xtail-1):        undefined on output
C       XX (xtail..xsize):      LU factors of this matrix, see below
C
C       Info:           integer informational output, see MA38AD
C       Rinfo:          real informational output, see MA38AD
C
C       xuse:           memory usage in Value
C
C       noutsd:         entries not in prior LU pattern, incremented

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38PD
C       subroutines called:     MA38ND, MA38ZD, MA38OD, DGEMV,
C                               DGEMM, DTRSV, DTRSM
C       functions called:       ABS, MAX
        INTRINSIC ABS, MAX
        EXTERNAL MA38ND, MA38OD, MA38ZD, DGEMV, DGEMM, DTRSV, DTRSM

C=======================================================================
C  DESCRIPTION OF DATA STRUCTURES:
C=======================================================================

C-----------------------------------------------------------------------
C  Matrix being factorized:
C-----------------------------------------------------------------------
C
C  The input matrix is held in an arrowhead format.  For the kth pivot,
C  the nonzeros in the pivot row (A (k, k...n)) and pivot column
C  (A (k...n, k)) are stored in the kth arrowhead.  The kth arrowhead
C  is located in:
C       Ari (Cp (k+1) ... Cp (k)-1):    pattern
C       XX  (Cp (k+1) ... Cp (k)-1):    values
C
C  Suppose p is in the range Cp (k+1) to Cp (k)-1.  If Ari (p) is
C  greater than zero, then the entry is in row Ari (p), column k,
C  with value XX (p).  If Ari (p) is less than zero, then the entry is
C  in row k, column -Ari (p), with value XX (p).  The arrowheads are
C  stored in reverse order (arrowhead n, n-1, ... 2, 1) in Ari and XX.
C  Note that Cp (n+1) = 1 unless BTF is in use and the original matrix
C  is not preserved.   In all cases, the real part of the arrowhead
C  format (XX (Cp (n+1) ... Cp (1)-1)) is overwritten with the LU
C  factors.  The integer part (Ari (Cp (n+1) ... Cp (1)-1)) is not
C  overwritten, since MA38RD does not require dynamic allocation of
C  integer memory.

C-----------------------------------------------------------------------
C  Frontal matrices
C-----------------------------------------------------------------------
C
C   Each unassembled frontal matrix (element) is stored as follows:
C       total size: fscal integers, (fdimr*fdimc) reals
C
C       if e is an unassembled element, and not the current frontal
C       matrix:
C
C       fluip = LUp (e) pointer to LU arrowhead in II
C       fdimc = FRdimc (e)      column dimension of contribution block
C       fxp   = FRxp (e)        pointer to contribution block in XX
C       next  = FRnext (e)      pointer to next block in XX
C       prev  = FRprev (e)      pointer to previous block in XX
C       fdegr = abs (LUi (fluip+2))
C       fdegc = abs (LUi (fluip+2))
C       XX (fxp ... )
C               a 2-dimensional array, C (1..fdimc, 1..fdimr), where
C               fdimr = fdegr if the contribution block is compressed,
C               or fdimr = LUi (fluip+5) if not.  Note, however, that
C               fdimr is not needed.  The contribution block is stored
C               in C (1..fdegc, 1..fdegr) in the C (1..fdimc,...) array.
C
C               If memory is limited, garbage collection will occur.
C               In this case, the C (1..fdimc, 1..fdimr) array is
C               compressed to be just large enough to hold the
C               unassembled contribution block,
C               C (1..fdegc, 1..fdegr).

C-----------------------------------------------------------------------
C  Current frontal matrix
C-----------------------------------------------------------------------
C
C  ffxp points to current frontal matrix (contribution block and LU
C  factors).  For example, if fflefc = 4, fflefr = 6, luk = 3,
C  ffdimc = 8, ffdimr = 12, then "x" is a term in the contribution
C  block, "l" in L1, "u" in U1, "L" in L2, "U" in U2, and "." is unused.
C  XX (fxp) is "X". The first 3 pivot values (diagonal entries in U1)
C  are labelled 1, 2, and 3.  The frontal matrix is ffdimc-by-ffdimr.
C
C                   |----------- col 1 of L1 and L2, etc.
C                   V
C       X x x x x x L L L . . .
C       x x x x x x L L L . . .
C       x x x x x x L L L . . .
C       x x x x x x L L L . . .
C       U U U U U U 3 l l . . .         <- row 3 of U1 and U2
C       U U U U U U u 2 l . . .         <- row 2 of U1 and U2
C       U U U U U U u u 1 . . .         <- row 1 of U1 and U2
C       . . . . . . . . . . . .

C-----------------------------------------------------------------------
C  LU factors
C-----------------------------------------------------------------------
C
C   The LU factors are placed at the tail of XX.  If this routine
C   is factorizing a single block, then this description is for the
C   factors of the single block:
C
C       LUi (1..):      integer info. for LU factors
C       XX (xtail..xsize):      real values in LU factors
C
C   Each LU arrowhead (or factorized element) is stored as follows:
C   ---------------------------------------------------------------
C
C       total size: (7 + ludegc + ludegr + lunson) integers,
C                   (luk**2 + ludegc*luk + luk*ludegc) reals
C
C       If e is an LU arrowhead, then luip = LUp (e).
C
C       luxp   = LUi (luip) pointer to numerical LU arrowhead
C       luk    = LUi (luip+1) number of pivots in LU arrowhead
C       ludegr = LUi (luip+2) degree of last row of U (excl. diag)
C       ludegc = LUi (luip+3) degree of last col of L (excl. diag)
C       lunson = LUi (luip+4) number of children in assembly DAG
C       ffdimr = LUi (luip+5)
C       ffdimc = LUi (luip+6)
C                       max front size for this LU arrowhead is
C                       ffdimr-by-ffdimc, or zero if this LU arrowhead
C                       factorized within the frontal matrix of a prior
C                       LU arrowhead.
C       lucp   = (luip + 7)
C                       pointer to pattern of column of L
C       lurp   = lucp + ludegc
C                       pointer to patter of row of U
C       lusonp = lurp + ludegr
C                       pointer to list of sons in the assembly DAG
C       LUi (lucp ... lucp + ludegc - 1)
C                       row indices of column of L
C       LUi (lurp ... lurp + ludegr - 1)
C                       column indices of row of U
C       LUi (lusonp ... lusonp + lunson - 1)
C                       list of sons
C       XX (luxp...luxp + luk**2 + ludegc*luk + luk*ludegr - 1)
C                       pivot block (luk-by-luk) and the L block
C                       (ludegc-by-luk) in a single (luk+ludegc)-by-luk
C                       array, followed by the U block in a
C                       luk-by-ludegr array.
C
C   Pivot column/row pattern (also columns/rows in contribution block):
C       If the column/row index is negated, the column/row has been
C       assembled out of the frontal matrix into a subsequent frontal
C       matrix.  After factorization, the negative flags are removed.
C
C   List of sons:
C       1 <= son <= n:           son an LUson
C       n+1 <= son <= 2n:        son-n is an Uson
C       2n+n <= son <= 3n:       son-2n is a Lson

C-----------------------------------------------------------------------
C  Workspaces:
C-----------------------------------------------------------------------
C
C  WpC (1..ludegr):     holds the pivot column pattern
C                       (excluding the pivot row indices)
C
C  WpR (1..ludegr):     holds the pivot row pattern
C                       (excluding the pivot column indices)
C
C  WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               Otherwise, WiR (1..n) is < 0
C
C  WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               Otherwise, WiC (1..n) is < 0
C
C  Wm (1..degc) or Wm (1..fdegc):       a gathered copy of WiR
C  Wj (1..degc) or Wj (1..fdegc):       offset in pattern of a son

C-----------------------------------------------------------------------
C  Memory allocation in XX:
C-----------------------------------------------------------------------
C
C   XX (1..xhead):      values of original entries in arrowheads of
C                       matrix, values of contribution blocks, followed
C                       by the current frontal matrix.
C
C   mtail = nlu+2
C   mhead = nlu+1:      FRnext (mhead) points to the first contribution
C                       block in the head of XX.  The FRnext and FRprev
C                       arrays form a doubly-linked list.  Traversing
C                       the list from mhead to mtail gives the
C                       contribution blocks in ascending ordering of
C                       address (FRxp).  A block is free if FRdimc <= 0.
C                       The largest known free block in XX is pfree,
C                       located in
C                       XX (FRxp (pfree) ... FRxp (pfree) + xfree -1),
C                       unless pfree = 0, in which case no largest free
C                       block is known.

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER SWPCOL, SWPROW, FDIMC, K0, COLPOS, ROWPOS, PIVOT, FFPP,
     $          P, I, J, LUDEGR, LUDEGC, KPOS, SP, FFRP, FFCP, TYPE,
     $          FXP, LURP, LUCP, NEXT, FFLEFR, PREV, XHEAD, FDEGR,
     $          FFLEFC, K, XCDP, XDP, XSP, S, FDEGC, FLURP, FLUCP,
     $          COL, E, ROW, MHEAD, MTAIL, UXP, LUK, IO, FLUIP, LUSONP,
     $          FFSIZE, FFXP, FFDIMR, FFDIMC, XRDP, NPIV, NB, LUNSON,
     $          XNEED, LDIMR, LDIMC, LXP, PRL, XP, LUIP, PFREE, XFREE,
     $          XS, LUXP, FSP, FLP, FDP, DEGC, NZU, NZL, XRUSE
        LOGICAL PR3, ALLCOL, ALLROW
        DOUBLE PRECISION
     $          ONE, ZERO, X
        PARAMETER (ONE = 1.0D0, ZERO = 0.0D0)

C  Printing control:
C  -----------------
C  prl:     invalid entries printed if prl >= 3
C  io:      I/O unit for warning messages (printing invalid entries)
C  pr3:     true if invalid entries are to be printed when found
C
C  Current working array:
C  ----------------------
C  ffxp:    current working array is in XX (ffxp ... ffxp+ffsize-1)
C  ffsize:  size of current working array in XX
C  ffdimr:  row degree (number of columns) of current working array
C  ffdimc:  column degree (number of rows) of current working array
C  fflefr:  row degree (number of columns) of current contribution block
C  fflefc:  column degree (number of rows) of current contribution block
C  ffrp:     U2 block is in XX (ffrp ...)
C  ffcp:     L2 block is in XX (ffcp ...)
C  ffpp:     location in XX of the current pivot value
C
C  Current element:
C  ----------------
C  s:       current element being factorized
C  luip:    current element is in LUi (luip ...)
C  luk:     number of pivots in current element
C  ludegc:  degree of pivot column (excluding pivots themselves)
C  ludegr:  degree of pivot row (excluding pivots themselves)
C  ldimr:   row degree (number of columns) of current element
C  ldimc:   column degree (number of row) of current element
C  lucp:    pattern of col(s) of current element in LUi (lucp...)
C  lurp:    pattern of row(s) of current element in LUi (lurp...)
C  lusonp:  list of sons of current element is in LUi (lusonp...)
C  lunson:  number of sons of current element
C  sp:      pointer into list of sons of current element
C  luxp:    numerical values of LU arrowhead stored in XX (luxp ...)
C  lxp:     L2 block is stored in XX (lxp ...) when computed
C  uxp:     U2 block is stored in XX (uxp ...) when computed
C  nzu:     nonzeros above diagonal in U in current LU arrowhead
C  nzl:     nonzeros below diagonal in L in current LU arrowhead
C  swpcol:  the non-pivotal column to be swapped with pivot column
C  swprow:  the non-pivotal row to be swapped with pivot row
C  colpos:  position in WpR of the pivot column
C  rowpos:  position in WpC of the pivot row
C  kpos:    position in C to place pivot row/column
C  k:       current pivot is kth pivot of current element, k=1..luk
C  k0:      contribution block, C, has been updated with pivots 1..k0
C  npiv:    number of pivots factorized so far, excl. current element
C  pivot:   current pivot entry is A (pivot, pivot)
C  xcdp:    current pivot column is in XX (xcdp ...)
C  xrdp:    current pivot row is in XX (xrdp ...)
C
C  Son, or element other than current element:
C  -------------------------------------------
C  e:       an element other than s (a son of s, for example)
C  fluip:   LU arrowhead of e is in LUi (fluip ...)
C  fxp:     contribution block of son is in XX (fxp ...)
C  fdimc:   leading dimension of contribution block of a son
C  fdegr:   row degree of contribution block of son (number of columns)
C  fdegc:   column degree of contribution block of son (number of rows)
C  allcol:  true if all columns are present in son
C  allrow:  true if all rows are present in son
C  flucp:   pattern of col(s) of son in LUi (flucp...)
C  flurp:   pattern of row(s) of son in LUi (flurp...)
C  type:    an LUson (type = 1), Uson (type = 2) or Lson (type = 3)
C  degc:    compressed column offset vector of son is in Wj/Wm (1..degc)
C
C  Memory allocation:
C  ------------------
C  mhead:   nlu+1, head pointer for contribution block link list
C  mtail:   nlu+2, tail pointer for contribution block link list
C  prev:    FRprev (e) of the element e
C  next:    FRnext (e) of the element e
C  pfree:   FRxp (pfree) is the largest known free block in XX
C  xfree:   size of largest known free block in XX
C  xneed:   bare minimum memory currently needed in XX
C  xhead:   XX (1..xhead-1) is in use, XX (xhead ..) is free
C  xruse:   estimated memory needed in XX for next call to MA38BD,
C           assuming a modest number of garbage collections
C  xs:      size of a block of memory in XX
C
C  Other:
C  ------
C  xdp:     destination pointer, into XX
C  xsp:     source pointer, into XX
C  xp:      a pointer into XX
C  fsp:     source pointer, into XX
C  fsp:     destination pointer, into XX
C  flp:     last row/column in current contribution is in XX (flp...)
C  col:     a column index
C  row:     a row index
C  nb:      block size for tradeoff between Level-2 and Level-3 BLAS
C  p, i, j, x:  various uses

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

C       ----------------------------------------------------------------
C       get control parameters and initialize various scalars
C       ----------------------------------------------------------------

        IO = ICNTL (2)
        PRL = ICNTL (3)
        NB = MAX (1, ICNTL (7))
        NPIV = 0
        XHEAD = CP (1)
        XTAIL = XSIZE + 1
        XNEED = XUSE
        XRUSE = XUSE
        XRMAX = MAX (XRMAX, XRUSE)
        MHEAD = NLU+1
        MTAIL = NLU+2
        XFREE = -1
        PFREE = 0
        PR3 = PRL .GE. 3 .AND. IO .GE. 0

C       ----------------------------------------------------------------
C       initialize workspaces
C       ----------------------------------------------------------------

        DO 10 I = 1, N
           WIR (I) = -1
           WIC (I) = -1
10      CONTINUE

        DO 20 E = 1, NLU+2
           FRDIMC (E) = 0
           FRXP (E) = 0
           FRNEXT (E) = 0
           FRPREV (E) = 0
20      CONTINUE
        FRNEXT (MHEAD) = MTAIL
        FRPREV (MTAIL) = MHEAD
        FRXP (MHEAD) = XHEAD
        FRXP (MTAIL) = XHEAD

C       count the numerical assembly of the original matrix
        RINFO (2) = RINFO (2) + NZ

C       current working array is empty:
        FFLEFR = 0
        FFLEFC = 0
        FFSIZE = 0
        FFXP = XHEAD

C=======================================================================
C  Factorization [
C=======================================================================

        DO 600 S = 1, NLU

C=======================================================================
C  Get the next element to factorize
C=======================================================================

           LUIP = LUP (S)
           LUK = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUDEGR = LUI (LUIP+2)
           LUNSON = LUI (LUIP+4)
           LUCP = (LUIP + 7)
           LURP = LUCP + LUDEGC
           LUSONP = LURP + LUDEGR
           LDIMC = LUK + LUDEGC
           LDIMR = LUK + LUDEGR

C=======================================================================
C  Start new frontal matrix or merge with prior contribution block [
C=======================================================================

C          =============================================================
           IF (LUI (LUIP+6) .NE. 0) THEN
C          start new contribution block
C          =============================================================

C             ----------------------------------------------------------
C             clear the prior offsets
C             ----------------------------------------------------------

              DO 30 I = 1, FFLEFR
                 WIC (WPR (I)) = -1
30            CONTINUE
              DO 40 I = 1, FFLEFC
                 WIR (WPC (I)) = -1
40            CONTINUE

C             ----------------------------------------------------------
C             save prior contribution block (s-1), if it exists
C             ----------------------------------------------------------

              XS = FFLEFR * FFLEFC
              IF (FFSIZE .NE. 0) THEN
C                one more frontal matrix is finished
                 XNEED = XNEED - (FFSIZE - XS)
                 XRUSE = XRUSE - (FFSIZE - XS)
                 INFO (13) = INFO (13) + 1
C             else
C                prior contribution block does not exist
              ENDIF

              IF (FFLEFR .LE. 0 .OR. FFLEFC .LE. 0) THEN

C                -------------------------------------------------------
C                if prior contribution block nonexistent or empty
C                -------------------------------------------------------

                 XUSE = XUSE - (XHEAD - FRXP (MTAIL))
                 XHEAD = FRXP (MTAIL)

              ELSE

C                -------------------------------------------------------
C                prepare the prior contribution block for later assembly
C                -------------------------------------------------------

                 E = S - 1

C                count the numerical assembly
                 RINFO (2) = RINFO (2) + XS

                 IF (XS .LE. XFREE) THEN

C                   ----------------------------------------------------
C                   compress and store in a freed block
C                   ----------------------------------------------------

C                   place the new block in the list
                    XFREE = XFREE - XS
                    IF (PFREE .EQ. MTAIL) THEN
C                      place the new block at start of tail block
                       PREV = FRPREV (MTAIL)
                       NEXT = MTAIL
                       XDP = FRXP (MTAIL)
                       FRXP (MTAIL) = XDP + XS
                    ELSE
C                      place the new block at end of block
                       PREV = PFREE
                       NEXT = FRNEXT (PFREE)
                       XDP = FRXP (NEXT) - XS
                       IF (XFREE .EQ. 0 .AND. PFREE .NE. MHEAD) THEN
C                         delete the free block if its size is zero
                          PREV = FRPREV (PREV)
                          PFREE = 0
                          XFREE = -1
                       ENDIF
                    ENDIF
                    DO 60 J = 0, FFLEFR - 1
CFPP$ NODEPCHK L
                       DO 50 I = 0, FFLEFC - 1
                          XX (XDP+J*FFLEFC+I) = XX (FFXP+J*FFDIMC+I)
50                     CONTINUE
60                  CONTINUE
                    XUSE = XUSE - (XHEAD - FRXP (MTAIL))
                    XHEAD = FRXP (MTAIL)
                    FRXP (E) = XDP
                    FRDIMC (E) = FFLEFC

                 ELSE

C                   ----------------------------------------------------
C                   deallocate part of unused portion of frontal matrix
C                   ----------------------------------------------------

C                   leave the contribution block C (1:fflefc, 1:fflefr)
C                   at head of XX, with column dimension of ffdimc and
C                   space of size (fflefr-1)*ffdimc for the first
C                   fflefr columns, and fflefc for the last column.
                    XS = FFSIZE - (FFLEFC + (FFLEFR-1)*FFDIMC)
                    XHEAD = XHEAD - XS
                    XUSE = XUSE - XS
                    PREV = FRPREV (MTAIL)
                    NEXT = MTAIL
                    FRXP (MTAIL) = XHEAD
                    FRXP (E) = FFXP
                    FRDIMC (E) = FFDIMC
                 ENDIF

                 FRNEXT (PREV) = E
                 FRPREV (NEXT) = E
                 FRNEXT (E) = NEXT
                 FRPREV (E) = PREV

              ENDIF

              IF (PFREE .EQ. MTAIL) THEN
                 PFREE = 0
                 XFREE = -1
              ENDIF

C             ----------------------------------------------------------
C             allocate a new ffdimr-by-ffdimc frontal matrix
C             ----------------------------------------------------------

              FFDIMC = LUI (LUIP+6)
              FFDIMR = LUI (LUIP+5)
              FFSIZE = FFDIMR * FFDIMC
              FFXP = 0

C             ----------------------------------------------------------
C             allocate and zero the space, garbage collection if needed
C             ----------------------------------------------------------

              IF (FFSIZE .GT. XTAIL-XHEAD) THEN
                 INFO (15) = INFO (15) + 1
                 CALL MA38OD (XX, XSIZE, XHEAD, XUSE,
     $              LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $              FFXP, FFSIZE, PFREE, XFREE)
              ENDIF

              FFXP = XHEAD
              XHEAD = XHEAD + FFSIZE
              XUSE = XUSE + FFSIZE
              XNEED = XNEED + FFSIZE
              XRUSE = XRUSE + FFSIZE
              XRMAX = MAX (XRMAX, XRUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XNEED)
              IF (XHEAD .GT. XTAIL) THEN
C                error return, if not enough real memory:
                 GO TO 9000
              ENDIF

C             ----------------------------------------------------------
C             zero the frontal matrix
C             ----------------------------------------------------------

              DO 70 P = FFXP, FFXP + FFSIZE - 1
                 XX (P) = ZERO
70            CONTINUE

C             ----------------------------------------------------------
C             place pivot rows and columns in correct position
C             ----------------------------------------------------------

              DO 80 K = 1, LUK
                 WIC (NPIV + K) = (LDIMR - K) * FFDIMC
                 WIR (NPIV + K) =  LDIMC - K
80            CONTINUE

C             ----------------------------------------------------------
C             get the pivot row pattern of the new LU arrowhead
C             ----------------------------------------------------------

              DO 90 I = 0, LUDEGR - 1
                 COL = LUI (LURP+I)
                 WIC (COL) = I * FFDIMC
                 WPR (I+1) = COL
90            CONTINUE

C             ----------------------------------------------------------
C             get the pivot column pattern of the new LU arrowhead
C             ----------------------------------------------------------

              DO 100 I = 0, LUDEGC - 1
                 ROW = LUI (LUCP+I)
                 WIR (ROW) = I
                 WPC (I+1) = ROW
100           CONTINUE

C          =============================================================
           ELSE
C          merge with prior contribution block
C          =============================================================

C             ----------------------------------------------------------
C             prior block is located at XX (ffxp ... ffxp + ffsize - 1).
C             It holds a working array C (1..ffdimc, 1..ffdimr), with a
C             prior contribution block in C (1..fflefc, 1..fflefr).
C             The last pivot column pattern is WpC (1..fflefc), and
C             the last pivot row pattern is WpR (1..fflefr).  The
C             offsets WiR and WiC are:
C             WiR (WpC (i)) = i-1, for i = 1..fflefc, and -1 otherwise.
C             WiC (WpR (i)) = (i-1)*ffdimc, for i = 1..fflefr, else -1.
C             The prior LU arrowhead is an implicit LUson of the current
C             element (and is implicitly assembled into the same
C             frontal matrix).
C             ----------------------------------------------------------

C             ----------------------------------------------------------
C             zero the newly extended frontal matrix
C             ----------------------------------------------------------

C             zero the new columns in the contribution and LU blocks
C             C (1..ldimc, fflefr+1..ldimr) = 0
              DO 120 J = FFLEFR, LDIMR - 1
                 DO 110 I = 0, LDIMC - 1
                    XX (FFXP + J*FFDIMC + I) = ZERO
110              CONTINUE
120           CONTINUE

C             C (fflefc+1..ldimc, 1..fflefr) = 0
C             zero the new rows in the contribution and U blocks
              DO 140 I = FFLEFC, LDIMC - 1
CFPP$ NODEPCHK L
                 DO 130 J = 0, FFLEFR - 1
                    XX (FFXP + J*FFDIMC + I) = ZERO
130              CONTINUE
140           CONTINUE

C             ----------------------------------------------------------
C             move pivot rows and columns into correct position
C             ----------------------------------------------------------

              DO 220 K = 1, LUK

C                -------------------------------------------------------
C                kth pivot of frontal matrix, (npiv+k)th pivot of LU
C                -------------------------------------------------------

                 PIVOT = NPIV + K

C                -------------------------------------------------------
C                move the kth pivot column into position
C                -------------------------------------------------------

                 XSP = WIC (PIVOT)
                 KPOS = LDIMR - K + 1
                 XDP = (KPOS - 1) * FFDIMC
                 WIC (PIVOT) = XDP

                 IF (XSP .GE. 0) THEN
C                   pivot column is already in current frontal matrix,
C                   shift into proper position
                    COLPOS = (XSP / FFDIMC) + 1
                    FSP = FFXP + XSP
                    FDP = FFXP + XDP

                    IF (FFLEFR .LT. KPOS) THEN

                       IF (FFLEFR .EQ. COLPOS) THEN

C                         ----------------------------------------------
C                         move C(:,colpos) => C (:,kpos)
C                         C (:,colpos) = 0
C                         ----------------------------------------------
CFPP$ NODEPCHK L
                          DO 150 I = 0, LDIMC - 1
                             XX (FDP+I) = XX (FSP+I)
                             XX (FSP+I) = ZERO
150                       CONTINUE

                       ELSE

C                         ----------------------------------------------
C                         move C(:,colpos) => C (:,kpos)
C                         move C(:,fflefr) => C (:,colpos)
C                         C (:,fflefr) = 0
C                         ----------------------------------------------

                          FLP = FFXP + (FFLEFR - 1) * FFDIMC
CFPP$ NODEPCHK L
                          DO 160 I = 0, LDIMC - 1
                             XX (FDP+I) = XX (FSP+I)
                             XX (FSP+I) = XX (FLP+I)
                             XX (FLP+I) = ZERO
160                       CONTINUE

                          SWPCOL = WPR (FFLEFR)
                          WPR (COLPOS) = SWPCOL
                          WIC (SWPCOL) = XSP
                       ENDIF

                    ELSE IF (COLPOS .NE. KPOS) THEN

C                      -------------------------------------------------
C                      swap C (:,colpos) <=> C (:,kpos)
C                      -------------------------------------------------
CFPP$ NODEPCHK L
                       DO 180 I = 0, LDIMC - 1
                          X = XX (FDP+I)
                          XX (FDP+I) = XX (FSP+I)
                          XX (FSP+I) = X
180                    CONTINUE

                       SWPCOL = WPR (KPOS)
                       WPR (COLPOS) = SWPCOL
                       WIC (SWPCOL) = XSP
                    ENDIF

                    FFLEFR = FFLEFR - 1
                 ENDIF

C                -------------------------------------------------------
C                move the kth pivot row into position
C                -------------------------------------------------------

                 XSP = WIR (PIVOT)
                 KPOS = LDIMC - K + 1
                 XDP = (KPOS - 1)
                 WIR (PIVOT) = XDP

                 IF (XSP .GE. 0) THEN
C                   pivot row is already in current frontal matrix,
C                   shift into proper position
                    ROWPOS = XSP + 1
                    FSP = FFXP + XSP
                    FDP = FFXP + XDP

                    IF (FFLEFC .LT. KPOS) THEN

                       IF (FFLEFC .EQ. ROWPOS) THEN

C                         ----------------------------------------------
C                         move C(rowpos,:) => C (kpos,:)
C                         C (rowpos,:) = 0
C                         ----------------------------------------------
CFPP$ NODEPCHK L
                          DO 190 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC
                             XX (FDP+J) = XX (FSP+J)
                             XX (FSP+J) = ZERO
190                       CONTINUE

                       ELSE

C                         ----------------------------------------------
C                         move C(rowpos,:) => C (kpos,:)
C                         move C(fflefc,:) => C (rowpos,:)
C                         C (fflefc,:) = 0
C                         ----------------------------------------------

                          FLP = FFXP + (FFLEFC - 1)
CFPP$ NODEPCHK L
                          DO 200 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC
                             XX (FDP+J) = XX (FSP+J)
                             XX (FSP+J) = XX (FLP+J)
                             XX (FLP+J) = ZERO
200                       CONTINUE

                          SWPROW = WPC (FFLEFC)
                          WPC (ROWPOS) = SWPROW
                          WIR (SWPROW) = XSP
                       ENDIF

                    ELSE IF (ROWPOS .NE. KPOS) THEN

C                      -------------------------------------------------
C                      swap C (rowpos,:) <=> C (kpos,:)
C                      -------------------------------------------------
CFPP$ NODEPCHK L
                       DO 210 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC
                          X = XX (FDP+J)
                          XX (FDP+J) = XX (FSP+J)
                          XX (FSP+J) = X
210                    CONTINUE

                       SWPROW = WPC (KPOS)
                       WPC (ROWPOS) = SWPROW
                       WIR (SWPROW) = XSP
                    ENDIF

                    FFLEFC = FFLEFC - 1
                 ENDIF

220           CONTINUE

C             ----------------------------------------------------------
C             merge with pivot row pattern of new LU arrowhead
C             ----------------------------------------------------------

              I = FFLEFR
              DO 230 P = LURP, LURP + LUDEGR - 1
                 COL = LUI (P)
                 IF (WIC (COL) .LT. 0) THEN
                    WIC (COL) = I * FFDIMC
                    I = I + 1
                    WPR (I) = COL
                 ENDIF
230           CONTINUE

C             ----------------------------------------------------------
C             merge with pivot column pattern of new LU arrowhead
C             ----------------------------------------------------------

              I = FFLEFC
              DO 240 P = LUCP, LUCP + LUDEGC - 1
                 ROW = LUI (P)
                 IF (WIR (ROW) .LT. 0) THEN
                    WIR (ROW) = I
                    I = I + 1
                    WPC (I) = ROW
                 ENDIF
240           CONTINUE

           ENDIF

C=======================================================================
C  Done initializing frontal matrix ]
C=======================================================================

C=======================================================================
C  Assemble original arrowheads into the frontal matrix, and deallocate
C=======================================================================

C          -------------------------------------------------------------
C          current workspace usage:
C          -------------------------------------------------------------

C          WpC (1..ludegr):     holds the pivot column pattern
C                               (excluding the pivot row indices)
C
C          WpR (1..ludegr):     holds the pivot row pattern
C                               (excluding the pivot column indices)
C
C          C (1..ffdimr, 1..ffdimc):  space for the frontal matrix,
C               in XX (ffxp ... ffxp + ffsize - 1)
C
C          C (i,j) is located at XX (ffxp+((i)-1)+((j)-1)*ffdimc)
C
C          C (1..ludegc, 1..ludegr):            contribution block
C          C (ludegc+1..ludegc+luk, 1..ludegr):             U2 block
C          C (1..ludegc, ludegr+1..ludegr+luk):             L2 block
C          C (ludegc+1..ludegc+luk, ludegr+1..ludegr+luk):  L1\U1 block
C
C          WiR (row) >= 0 for each row in pivot column pattern.
C               offset into pattern is given by:
C               WiR (row) == offset - 1
C               Also, WiR (npiv+1 ... npiv+luk) is
C               ludegc+luk-1 ... ludegc, the offsets of the pivot rows.
C
C               Otherwise, WiR (1..n) is < 0
C
C          WiC (col) >= 0 for each col in pivot row pattern.
C               WiC (col) == (offset - 1) * ffdimc
C               Also, WiC (npiv+1 ... npiv+luk) is
C               ludegr+luk-1 ... ludegr, the offsets of the pivot rows.
C
C               Otherwise, WiC (1..n) is < 0

           DO 260 K = 1, LUK
              I = NPIV + K
              XCDP = FFXP + WIC (I)
              XRDP = FFXP + WIR (I)
              DO 250 P = CP (I+1), CP (I) - 1
                 J = ARI (P)
                 IF (J .GT. 0) THEN
C                   a diagonal entry, or lower triangular entry
C                   row = j, col = i
                    XP = XCDP + WIR (J)
                    IF (XP .LT. XCDP) THEN
C                      invalid entry - not in prior LU pattern
                       NOUTSD = NOUTSD + 1
                       IF (PR3) THEN
C                         get original row and column index and print it
                          ROW = RPERM (J)
                          COL = CPERM (I)
                          CALL MA38ZD (2, 97, ROW, COL, XX (P), IO)
                       ENDIF
                    ELSE
                       XX (XP) = XX (XP) + XX (P)
                    ENDIF
                 ELSE
C                   an upper triangular entry
C                   row = i, col = -j
                    XP = XRDP + WIC (-J)
                    IF (XP .LT. XRDP) THEN
C                      invalid entry - not in prior LU pattern
                       NOUTSD = NOUTSD + 1
                       IF (PR3) THEN
C                         get original row and column index and print it
                          ROW = RPERM (I)
                          COL = CPERM (-J)
                          CALL MA38ZD (2, 97, ROW, COL, XX (P), IO)
                       ENDIF
                    ELSE
                       XX (XP) = XX (XP) + XX (P)
                    ENDIF
                 ENDIF
250           CONTINUE
260        CONTINUE

C          deallocate the original arrowheads
           P = CP (NPIV + LUK + 1)
           XS = CP (NPIV + 1) - P
           FRXP (MHEAD) = P
           XNEED = XNEED - XS
           IF (XS .GT. XFREE) THEN
              XFREE = XS
              PFREE = MHEAD
           ENDIF

C=======================================================================
C  Assemble LUsons, Usons, and Lsons into the frontal matrix [
C=======================================================================

           DO 480 SP = LUSONP, LUSONP + LUNSON - 1

C             ----------------------------------------------------------
C             get the son and determine its type (LUson, Uson, or Lson)
C             ----------------------------------------------------------

              E = LUI (SP)
              IF (E .LE. N) THEN
C                LUson
                 TYPE = 1
              ELSE IF (E .LE. 2*N) THEN
C                Uson
                 E = E - N
                 TYPE = 2
              ELSE
C                Lson
                 E = E - 2*N
                 TYPE = 3
              ENDIF

C             ----------------------------------------------------------
C             if fdimc=0 this is the implicit LUson (already assembled)
C             ----------------------------------------------------------

              FDIMC = FRDIMC (E)
              IF (FDIMC .NE. 0) THEN

C                -------------------------------------------------------
C                get scalar info of the son (it needs assembling)
C                -------------------------------------------------------

                 FXP = FRXP (E)
                 FLUIP = LUP (E)
                 FDEGR = LUI (FLUIP+2)
                 FDEGC = LUI (FLUIP+3)
                 ALLCOL = FDEGR .GT. 0
                 ALLROW = FDEGC .GT. 0
                 FDEGR = ABS (FDEGR)
                 FDEGC = ABS (FDEGC)
                 FLUCP = (FLUIP + 7)
                 FLURP = FLUCP + FDEGC

C                use Wm (1..fdegc) for offsets:

C                -------------------------------------------------------
                 IF (TYPE .EQ. 1) THEN
C                this is an LUson - assemble an entire frontal matrix
C                -------------------------------------------------------

C                   ----------------------------------------------------
                    IF (ALLROW) THEN
C                   no rows assembled out of this LUson yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DO 270 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          WM (I+1) = WIR (ROW)
270                    CONTINUE

C                      -------------------------------------------------
                       IF (ALLCOL) THEN
C                      no rows or cols assembled out of LUson yet
C                      -------------------------------------------------

                          DO 290 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 280 I = 0, FDEGC-1
                                XX (XDP + WM (I+1)) =
     $                          XX (XDP + WM (I+1)) +
     $                          XX (FXP + J*FDIMC + I)
280                          CONTINUE
290                       CONTINUE

C                      -------------------------------------------------
                       ELSE
C                      some columns already assembled out of LUson
C                      -------------------------------------------------

                          DO 310 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             IF (COL .GT. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 300 I = 0, FDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
300                             CONTINUE
                             ENDIF
310                       CONTINUE

                       ENDIF

C                   ----------------------------------------------------
                    ELSE
C                   some rows already assembled out of LUson
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DEGC = 0
                       DO 320 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          IF (ROW .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                          ENDIF
320                    CONTINUE

C                      -------------------------------------------------
                       IF (ALLCOL) THEN
C                      some rows already assembled out of LUson
C                      -------------------------------------------------

                          DO 340 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 330 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
330                          CONTINUE
340                       CONTINUE

C                      -------------------------------------------------
                       ELSE
C                      rows and columns already assembled out of LUson
C                      -------------------------------------------------

                          DO 360 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             IF (COL .GT. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 350 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
350                             CONTINUE
                             ENDIF
360                       CONTINUE

                       ENDIF
                    ENDIF

C                   ----------------------------------------------------
C                   deallocate the LUson frontal matrix
C                   ----------------------------------------------------

                    FRDIMC (E) = 0
                    PREV = FRPREV (E)
                    NEXT = FRNEXT (E)
                    XNEED = XNEED - FDEGR*FDEGC
                    XRUSE = XRUSE - FDEGR*FDEGC

                    IF (FRDIMC (PREV) .LE. 0) THEN
C                      previous block is free - delete this block
                       FRNEXT (PREV) = NEXT
                       FRPREV (NEXT) = PREV
                       E = PREV
                       PREV = FRPREV (E)
                    ENDIF

                    IF (FRDIMC (NEXT) .LE. 0) THEN
C                      next block is free - delete this block
                       FRXP (NEXT) = FRXP (E)
                       IF (E .LE. NLU) THEN
                          FRNEXT (PREV) = NEXT
                          FRPREV (NEXT) = PREV
                       ENDIF
                       E = NEXT
                       NEXT = FRNEXT (E)
                       IF (FRNEXT (MHEAD) .EQ. MTAIL) THEN
C                         no blocks left except mhead and mtail
                          FRXP (MTAIL) = FRXP (MHEAD)
                       ENDIF
                    ENDIF

C                   get the size of the freed block
                    IF (NEXT .EQ. 0) THEN
C                      this is the mtail block
                       XS = FFXP - FRXP (E)
                    ELSE
                       XS = FRXP (NEXT) - FRXP (E)
                    ENDIF
                    IF (XS .GT. XFREE) THEN
C                      keep track of the largest free block
                       XFREE = XS
                       PFREE = E
                    ENDIF

C                -------------------------------------------------------
                 ELSE IF (TYPE .EQ. 2) THEN
C                Uson:  assemble all possible columns
C                -------------------------------------------------------

C                   ----------------------------------------------------
                    IF (ALLROW) THEN
C                   no rows assembled out of this Uson yet
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DO 370 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          WM (I+1) = WIR (ROW)
370                    CONTINUE

                       DO 390 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN
                             IF (WIC (COL) .GE. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 380 I = 0, FDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
380                             CONTINUE
C                               flag this column as assembled
                                LUI (FLURP+J) = -COL
                             ENDIF
                          ENDIF
390                    CONTINUE

C                   ----------------------------------------------------
                    ELSE
C                   some rows already assembled out of this Uson
C                   ----------------------------------------------------

C                      compute the compressed column offset vector
                       DEGC = 0
                       DO 400 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          IF (ROW .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                          ENDIF
400                    CONTINUE

                       DO 420 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN
                             IF (WIC (COL) .GE. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 410 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
410                             CONTINUE
C                               flag this column as assembled
                                LUI (FLURP+J) = -COL
                             ENDIF
                          ENDIF
420                    CONTINUE

                    ENDIF

C                   flag this element as missing some columns
                    LUI (FLUIP+2) = -FDEGR

C                -------------------------------------------------------
                 ELSE
C                Lson:  assemble all possible rows
C                -------------------------------------------------------

C                   compute the compressed column offset vector
                    DEGC = 0
                    DO 430 I = 0, FDEGC-1
                       ROW = LUI (FLUCP+I)
                       IF (ROW .GT. 0) THEN
                          IF (WIR (ROW) .GE. 0) THEN
C                            this row will be assembled in loop below
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
C                            flag this row as assembled
                             LUI (FLUCP+I) = -ROW
                          ENDIF
                       ENDIF
430                 CONTINUE

C                   ----------------------------------------------------
                    IF (ALLCOL) THEN
C                   no columns assembled out of this Lson yet
C                   ----------------------------------------------------

                       DO 450 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                          DO 440 I = 1, DEGC
                             XX (XDP + WM (I)) =
     $                       XX (XDP + WM (I)) +
     $                       XX (FXP + J*FDIMC + WJ (I))
440                       CONTINUE
450                    CONTINUE

C                   ----------------------------------------------------
                    ELSE
C                   some columns already assembled out of this Lson
C                   ----------------------------------------------------

                       DO 470 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 460 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
460                          CONTINUE
                          ENDIF
470                    CONTINUE

                    ENDIF

C                   flag this element as missing some rows
                    LUI (FLUIP+3) = -FDEGC

                 ENDIF

              ENDIF

480        CONTINUE

C=======================================================================
C  Done assemblying sons into the frontal matrix ]
C=======================================================================

C=======================================================================
C  Factorize the frontal matrix [
C=======================================================================

           K0 = 0
           FFLEFR = LDIMR
           FFLEFC = LDIMC
           FFCP = FFXP + FFLEFR * FFDIMC
           FFRP = FFXP + FFLEFC
           FFPP = FFXP + FFLEFC + FFLEFR * FFDIMC

           DO 500 K = 1, LUK

C             ----------------------------------------------------------
C             compute kth column of U1, and update pivot column
C             ----------------------------------------------------------

              IF (K-K0-2 .GT. 0) THEN
C                u1 = L1 \ u1.  Note that L1 transpose is stored, and
C                that u1 is stored with rows in reverse order.
                 CALL DTRSV ('U', 'N', 'U', K-K0-1,
     $                         XX (FFPP         ), FFDIMC,
     $                         XX (FFPP - FFDIMC), 1)
                 RINFO (5) = RINFO (5) + (K-K0-2)*(K-K0-1)
              ENDIF
              IF (K-K0-1 .GT. 0) THEN
C                l1 = l1 - L2*u1
                 CALL DGEMV ('N', FFLEFC, K-K0-1,
     $                   -ONE, XX (FFCP         ), FFDIMC,
     $                         XX (FFPP - FFDIMC), 1,
     $                    ONE, XX (FFCP - FFDIMC), 1)
                 RINFO (5) = RINFO (5) + 2*FFLEFC*(K-K0-1)
              ENDIF

              FFCP = FFCP - FFDIMC
              FFRP = FFRP - 1
              FFPP = FFPP - FFDIMC - 1
              FFLEFR = FFLEFR - 1
              FFLEFC = FFLEFC - 1

C             ----------------------------------------------------------
C             divide pivot column by pivot
C             ----------------------------------------------------------

C             k-th pivot in frontal matrix located in XX (ffpp)
              X = XX (FFPP)
              IF (ABS (X) .EQ. ZERO) THEN
C                error return, if pivot order from MA38AD not acceptable
                 GO TO 9010
              ENDIF
              X = ONE / X
              DO 490 P = FFCP, FFCP + FFLEFC - 1
                 XX (P) = XX (P) * X
490           CONTINUE
C             count this as a call to the Level-1 BLAS:
              RINFO (4) = RINFO (4) + FFLEFC
              INFO (17) = INFO (17) + 1

C             ----------------------------------------------------------
C             compute U1 (k0+1..k, k..ldimc) and
C             update contribution block: rank-nb, or if last pivot
C             ----------------------------------------------------------

              IF (K-K0 .GE. NB .OR. K .EQ. LUK) THEN
                 CALL DTRSM ('L', 'U', 'N', 'U', K-K0, FFLEFR, ONE,
     $                      XX (FFPP), FFDIMC,
     $                      XX (FFRP), FFDIMC)
                 CALL DGEMM ('N', 'N', FFLEFC, FFLEFR, K-K0,
     $                -ONE, XX (FFCP ), FFDIMC,
     $                      XX (FFRP ), FFDIMC,
     $                 ONE, XX (FFXP), FFDIMC)
                 RINFO (6) = RINFO (6) + FFLEFR*(K-K0-1)*(K-K0)
     $                                         + 2*FFLEFC*FFLEFR*(K-K0)
                 K0 = K
              ENDIF

500        CONTINUE

C=======================================================================
C  Done factorizing the frontal matrix ]
C=======================================================================

C=======================================================================
C  Save the new LU arrowhead [
C=======================================================================

C          allocate permanent space for the LU arrowhead
           XS = LUK*LUDEGC + LUK*LUDEGR + LUK*LUK

           IF (XS .GT. XTAIL-XHEAD) THEN
              INFO (15) = INFO (15) + 1
              CALL MA38OD (XX, XSIZE, XHEAD, XUSE,
     $              LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $              FFXP, FFSIZE, PFREE, XFREE)
           ENDIF

           XTAIL = XTAIL - XS
           LUXP = XTAIL
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           XRUSE = XRUSE + XS
           XRMAX = MAX (XRMAX, XRUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (XHEAD .GT. XTAIL) THEN
C             error return, if not enough real memory:
              GO TO 9000
           ENDIF

C          save the scalar data of the LU arrowhead
           LUI (LUIP) = LUXP

C          save column pattern (it may have been rearranged)
           DO 510 I = 0, LUDEGC-1
              LUI (LUCP+I) = WPC (I+1)
510        CONTINUE

C          save row pattern (it may have been rearranged)
           DO 520 I = 0, LUDEGR-1
              LUI (LURP+I) = WPR (I+1)
520        CONTINUE

C          move the L1,U1 matrix, compressing the dimension from
C          ffdimc to ldimc.  The LU arrowhead grows on top of stack.
           XP = FFXP + (LDIMR-1)*FFDIMC + LDIMC-1
           DO 540 J = 0, LUK-1
CFPP$ NODEPCHK L
              DO 530 I = 0, LUK-1
                 XX (LUXP + J*LDIMC + I) = XX (XP - J*FFDIMC - I)
530           CONTINUE
540        CONTINUE

C          move L2 matrix, compressing dimension from ffdimc to ldimc
           IF (LUDEGC .NE. 0) THEN
              LXP = LUXP + LUK
              XP = FFXP + (LDIMR-1)*FFDIMC
              DO 560 J = 0, LUK-1
CFPP$ NODEPCHK L
                 DO 550 I = 0, LUDEGC-1
                    XX (LXP + J*LDIMC + I) = XX (XP - J*FFDIMC + I)
550              CONTINUE
560           CONTINUE
           ENDIF

C          move the U2 block.
           IF (LUDEGR .NE. 0) THEN
              UXP = LUXP + LUK * LDIMC
              XP = FFXP + LDIMC-1
              DO 580 J = 0, LUDEGR-1
CFPP$ NODEPCHK L
                 DO 570 I = 0, LUK-1
                    XX (UXP + J*LUK + I) = XX (XP + J*FFDIMC - I)
570              CONTINUE
580           CONTINUE
           ENDIF

C          one more LU arrowhead has been refactorized
           NZU = (LUK*(LUK-1)/2) + LUK*LUDEGC
           NZL = (LUK*(LUK-1)/2) + LUK*LUDEGR
           INFO (10) = INFO (10) + NZL
           INFO (11) = INFO (11) + NZU

C          -------------------------------------------------------------
C          clear the pivot row and column offsets
C          -------------------------------------------------------------

           DO 590 PIVOT = NPIV + 1, NPIV + LUK
              WIR (PIVOT) = -1
              WIC (PIVOT) = -1
590        CONTINUE
           NPIV = NPIV + LUK

C=======================================================================
C  Done saving the new LU arrowhead ]
C=======================================================================

600     CONTINUE

C=======================================================================
C  Factorization complete ]
C=======================================================================

C=======================================================================
C  Wrap-up:  store LU factors in their final form
C=======================================================================

C       ----------------------------------------------------------------
C       Flag remaining arrowheads as invalid entries, if prior matrix
C       was singular.  Print them if requested.
C       ----------------------------------------------------------------

        IF (NPIV .LT. N) THEN
           IF (PR3) THEN
              DO 620 I = NPIV+1, N
                 DO 610 P = CP (I+1), CP (I) - 1
                    J = ARI (P)
                    IF (J .GT. 0) THEN
C                      a diagonal entry, or lower triangular entry
C                      get original row and column index
                       ROW = RPERM (J)
                       COL = CPERM (I)
                    ELSE
C                      an upper triangular entry
C                      get original row and column index
                       ROW = RPERM (I)
                       COL = CPERM (-J)
                    ENDIF
                    CALL MA38ZD (2, 95, ROW, COL, XX(P), IO)
610              CONTINUE
620           CONTINUE
           ENDIF
           NOUTSD = NOUTSD + (CP (NPIV+1) - CP (N+1))
        ENDIF

C       ----------------------------------------------------------------
C       deallocate all remaining input arrowheads and frontal matrices
C       ----------------------------------------------------------------

        IF (FFSIZE .NE. 0) THEN
           INFO (13) = INFO (13) + 1
        ENDIF
        XUSE = XUSE - (XHEAD - CP (N+1))
        XNEED = XUSE
        XHEAD = CP (N+1)

        IF (NLU .EQ. 0) THEN
C          LU factors are completely empty (A = 0).
C          Add one real, to simplify rest of code.
C          Otherwise, some arrays in MA38BD or MA38CD would have
C          zero size, which can cause an address fault.
           XTAIL = XSIZE
           XUSE = XUSE + 1
           XRUSE = XUSE
           XNEED = XUSE
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
        ENDIF

        IF (XHEAD .LE. XTAIL) THEN

C          -------------------------------------------------------------
C          sufficient memory to complete the factorization
C          -------------------------------------------------------------

           IF (NLU .EQ. 0) THEN
C             zero the dummy entry, although it won't be accessed:
              XX (XTAIL) = ZERO
           ENDIF

C          -------------------------------------------------------------
C          update pointers in LU factors
C          -------------------------------------------------------------

           DO 630 S = 1, NLU
              LUIP = LUP (S)
              LUXP = LUI (LUIP)
              LUI (LUIP) = LUXP - XTAIL + 1
630        CONTINUE

C          -------------------------------------------------------------
C          get memory usage estimate for next call to MA38BD
C          -------------------------------------------------------------

           XRUSE = XUSE
           XRMAX = MAX (XRMAX, XRUSE)
           RETURN

        ENDIF

C=======================================================================
C  Error conditions
C=======================================================================

C       error return label:
9000    CONTINUE
C       out of real memory
        CALL MA38ND (2, ICNTL, INFO, -4, INFO (21))
        RETURN

C       error return label:
9010    CONTINUE
C       original pivot order computed by MA38AD is no longer acceptable
        CALL MA38ND (2, ICNTL, INFO, -6, 0)
        RETURN
        END

        SUBROUTINE MA38SD (NLU, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)

C=== MA38SD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  solves U'x = b, where U is the upper triangular factor of a matrix
C  (if BTF not used) or a single diagonal block (if BTF is used).
C  B is overwritten with the solution X.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       nlu:            number of LU arrowheads in the LU factors
C       npiv:           number of pivots found (normally n)
C       n:              order of matrix
C       LUp (1..nlu):   pointer to LU arrowheads in LUi
C       LUi ( ... ):    integer values of LU arrowheads
C       LUx ( ... ):    real values of LU arroheads
C       X (1..n):       the right-hand-side

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       the solution to U'x=b

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38JD
C       subroutines called:     DTRSV, DGEMV
        EXTERNAL  DTRSV, DGEMV

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I, K, S, LUIP, LUXP, LUK, LUDEGR, LUDEGC, LURP, UXP,
     $          LUCP, ROW
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)

C  s:       an element, or LU arrowhead
C  k:       kth pivot
C  i:       ith column in U2' array in element s
C  luip:    s is in LUi (luip...)
C  luxp:    real part of s is in LUx (luxp...)
C  luk:     number of pivots in s
C  ludegc:  column degree of non-pivotal part of s
C  ludegr:  row degree of non-pivotal part of s
C  lucp:    pattern of column of s in LUi (lucp...lucp+ludegc-1)
C  lurp:    pattern of row of s in LUi (lurp...lurp+ludegr-1)
C  uxp:     the luk-by-ludegr U2 block of s is in LUx (uxp...)
C  row:     row index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        K = 0
        DO 40 S = 1, NLU

C          -------------------------------------------------------------
C          get s-th LU arrowhead (s = 1..nlu, in pivotal order)
C          -------------------------------------------------------------

           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGR = LUI (LUIP+2)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LURP   = LUCP + LUDEGC

           IF (LUK .EQ. 1) THEN

C             ----------------------------------------------------------
C             only one pivot, stride-1 sparse saxpy
C             ----------------------------------------------------------

              K = K + 1
C             divide by pivot, U (k,k): LUx (luxp)
              X (K) = X (K) / LUX (LUXP)
              UXP = LUXP + LUDEGC + 1
CFPP$ NODEPCHK L
              DO 10 I = 1, LUDEGR
                 ROW = LUI (LURP+I-1)
C                col: k, U (row,col): LUx (uxp+i-1)
                 X (ROW) = X (ROW) - LUX (UXP+I-1) * X (K)
10            CONTINUE

           ELSE

C             ----------------------------------------------------------
C             more than one pivot
C             ----------------------------------------------------------

              UXP = LUXP + LUK * (LUDEGC + LUK)
              CALL DTRSV ('U', 'T', 'N', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)
              DO 20 I = 1, LUDEGR
                 ROW = LUI (LURP+I-1)
                 W (I) = X (ROW)
20            CONTINUE
              IF(LUDEGR.GT.0)CALL DGEMV ('T', LUK, LUDEGR, -ONE,
     $           LUX (UXP), LUK, X (K+1), 1, ONE, W, 1)
              DO 30 I = 1, LUDEGR
                 ROW = LUI (LURP+I-1)
                 X (ROW) = W (I)
30            CONTINUE
              K = K + LUK

           ENDIF

40      CONTINUE
        RETURN
        END

        SUBROUTINE MA38TD (NLU, NPIV, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, NPIV, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)

C=== MA38TD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  solves L'x = b, where L is the lower triangular factor of a matrix
C  (if BTF not used) or a single diagonal block (if BTF is used).
C  B is overwritten with the solution X.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       nlu:            number of LU arrowheads in the LU factors
C       npiv:           number of pivots found (normally n)
C       n:              order of matrix
C       LUp (1..nlu):   pointer to LU arrowheads in LUi
C       LUi ( ... ):    integer values of LU arrowheads
C       LUx ( ... ):    real values of LU arroheads
C       X (1..n):       the right-hand-side

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n):

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       the solution to L'x=b

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38JD
C       subroutines called:     DTRSV, DGEMV
        EXTERNAL DTRSV, DGEMV

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER J, K, S, LUIP, LUXP, LUK, LUDEGC, LUCP, LXP, COL
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)

C  s:       an element, or LU arrowhead
C  k:       kth pivot
C  j:       jth column in L2' array in element s
C  luip:    s is in LUi (luip...)
C  luxp:    real part of s is in LUx (luxp...)
C  luk:     number of pivots in s
C  ludegc:  column degree of non-pivotal part of s
C  lucp:    pattern of column of s in LUi (lucp...lucp+ludegc-1)
C  lxp:     the ludegc-by-luk L2 block of s is in LUx (lxp...)
C  col:     column index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        K = NPIV
        DO 30 S = NLU, 1, -1

C          -------------------------------------------------------------
C          get s-th LU arrowhead (s = nlu..1, in reverse pivotal order)
C          -------------------------------------------------------------

           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LXP    = LUXP + LUK

           IF (LUK .EQ. 1) THEN

C             ----------------------------------------------------------
C             only one pivot, stride-1 sparse dot product
C             ----------------------------------------------------------

CFPP$ NODEPCHK L
              DO 10 J = 1, LUDEGC
                 COL = LUI (LUCP+J-1)
C                row: k, U (row,col): LUx (lxp+j-1)
                 X (K) = X (K) - LUX (LXP+J-1) * X (COL)
10            CONTINUE
C             L (k,k) is one
              K = K - 1

           ELSE

C             ----------------------------------------------------------
C             more than one pivot
C             ----------------------------------------------------------

              K = K - LUK
              DO 20 J = 1, LUDEGC
                 COL = LUI (LUCP+J-1)
                 W (J) = X (COL)
20            CONTINUE
              CALL DGEMV ('T', LUDEGC, LUK, -ONE,
     $           LUX (LXP), LUDEGC + LUK, W, 1, ONE, X (K+1), 1)
              CALL DTRSV ('L', 'T', 'U', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)

           ENDIF

30      CONTINUE
        RETURN
        END

        SUBROUTINE MA38UD (NLU, NPIV, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, NPIV, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)

C=== MA38UD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  solves Ux = b, where U is the upper triangular factor of a matrix
C  (if BTF not used) or a single diagonal block (if BTF is used).
C  B is overwritten with the solution X.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       nlu:            number of LU arrowheads in the LU factors
C       npiv:           number of pivots found (normally n)
C       n:              order of matrix
C       LUp (1..nlu):   pointer to LU arrowheads in LUi
C       LUi ( ... ):    integer values of LU arrowheads
C       LUx ( ... ):    real values of LU arroheads
C       X (1..n):       the right-hand-side

C=======================================================================
C  WORKSPACE:
C=======================================================================
C
C       W (1..n)

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       X (1..n):       the solution to Ux=b

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   MA38JD
C       subroutines called:     DTRSV, DGEMV
        EXTERNAL DTRSV, DGEMV

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER J, K, S, LUIP, LUXP, LUK, LUDEGR, LUDEGC, LURP, UXP,
     $          LUCP, COL
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)

C  s:       an element, or LU arrowhead
C  k:       kth pivot
C  j:       jth column in U2 array in element s
C  luip:    s is in LUi (luip...)
C  luxp:    real part of s is in LUx (luxp...)
C  luk:     number of pivots in s
C  ludegc:  column degree of non-pivotal part of s
C  ludegr:  row degree of non-pivotal part of s
C  lucp:    pattern of column of s in LUi (lucp...lucp+ludegc-1)
C  lurp:    pattern of row of s in LUi (lurp...lurp+ludegr-1)
C  uxp:     the luk-by-ludegr U2 block of s is in LUx (uxp...)
C  col:     column index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

        K = NPIV
        DO 30 S = NLU, 1, -1

C          -------------------------------------------------------------
C          get s-th LU arrowhead (s = nlu..1, in reverse pivotal order)
C          -------------------------------------------------------------

           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGR = LUI (LUIP+2)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LURP   = LUCP + LUDEGC
           UXP    = LUXP + LUK * (LUDEGC + LUK)

           IF (LUK .EQ. 1) THEN

C             ----------------------------------------------------------
C             only one pivot, stride-1 sparse dot product
C             ----------------------------------------------------------

CFPP$ NODEPCHK L
              DO 10 J = 1, LUDEGR
                 COL = LUI (LURP+J-1)
C                row: k, U (row,col): LUx (uxp+j-1)
                 X (K) = X (K) - LUX (UXP+J-1) * X (COL)
10            CONTINUE
C             divide by pivot, U (k,k): LUx (luxp)
              X (K) = X (K) / LUX (LUXP)
              K = K - 1

           ELSE

C             ----------------------------------------------------------
C             more than one pivot
C             ----------------------------------------------------------

              K = K - LUK
              DO 20 J = 1, LUDEGR
                 COL = LUI (LURP+J-1)
                 W (J) = X (COL)
20            CONTINUE
              IF(LUDEGR.GT.0) CALL DGEMV ('N', LUK, LUDEGR, -ONE,
     $           LUX (UXP), LUK, W, 1, ONE, X (K+1), 1)
              CALL DTRSV ('U', 'N', 'N', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)

           ENDIF

30      CONTINUE
        RETURN
        END

        SUBROUTINE MA38YD (WHO, WHERE,
     *          N, NE, JOB, TRANS, LVALUE, LINDEX, VALUE,
     *          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     *          B, X, LX, W, LW)
        INTEGER WHO, WHERE, N, NE, JOB, LVALUE, LINDEX, INDEX(LINDEX),
     *          KEEP(20), ICNTL(20), INFO(40), LX, LW
        DOUBLE PRECISION VALUE(LVALUE), CNTL(10), RINFO(20), B(LX),
     *          X(LX), W(LW)
        LOGICAL TRANS

C=== MA38YD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  print input/output arguments for MA38AD, MA38BD, and MA38CD

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be deleted on installation (replaced with a dummy
C  routine that just returns without printing) in order to completely
C  disable the printing of all input/output parameters.  To completely
C  disable all I/O, you can also replace the MA38ZD routine with a
C  dummy subroutine.  If you make this modification, please do
C  not delete any original code - just comment it out instead.  Add a
C  comment and date to your modifications.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who:            what routine called MA38YD:
C                       1: MA38AD, 2: MA38BD, 3: MA38CD
C       where:          called from where:
C                       1: entry of routine, else exit of routine
C       Icntl (3):      if < 3 then print nothing, if 3 then print
C                       terse output, if 4 print info and matrix
C                       values.  If 5, print everything.
C       Icntl (2):      I/O unit on which to print.  No printing
C                       occurs if < 0.
C
C       Parameters to print, see MA38AD, MA38BD, or MA38CD for
C       descriptions:
C
C           n, ne, job, trans, lvalue, lindex, Value, Index, Keep,
C           Icntl, Info, Rinfo, B, X, lx, W, lw

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       on Icntl (2) I/O unit only

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  MA38AD, MA38BD, MA38CD
C       functions called:       MIN
        INTRINSIC MIN

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        LOGICAL TRANSA, TRANSC, PRLU, BADLU, SGLTON, PRESRV, SYMBOL
        INTEGER IO, PRL, PRN, K, LUI1, LUI2, LUX1, LUX2, ROW, COL,
     $          FACNE, FACN, NZ, FACJOB, NBLKS, NZOFF, FACTRA, CPERMP,
     $          RPERMP, APP, AXP, AIP, OFFIP, OFFXP, LUBLPP, OFFPP,
     $          BLKPP, P1, P2, P, BLK, K1, K2, KN, LUIIP, LUXXP, NPIV,
     $          NLU, E, LUK, LUPP, LUIP, LUXP, LUDEGR, LUDEGC, LUNSON,
     $          LUSONP, LUCP, LURP, I, J, NZCOL, NZROW, UXP, SON,
     $          PRMAX, LUDIMR, LUDIMC, MAXDR, MAXDC, LUIR1, IP1, IP2,
     $          XP1
        DOUBLE PRECISION
     $          ONE
        PARAMETER (PRMAX = 10, ONE = 1.0D0)

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C  prn:     number of entries printed so far
C  prmax:   maximum number of entries to print if prl = 3
C  prlu:    true if printing LU factors
C
C  Location and status of LU factors:
C  ----------------------------------
C  transc:  TRANSC argument in MA38CD
C  transa:  TRANSA argument in MA38AD or MA38BD when matrix factorized
C  badlu:   true if LU factors uncomputed or corrupted
C  presrv:  true if original matrix was preserved when factorized
C  symbol:  true if only symbolic part of LU factors needed on input
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for a call to MA38BD
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors start in Value (lux1...)
C  lux2:    real part of LU factors end in Value (...lux1)
C  ip1:     pointer into leading part of LU factors in Index
C  ip2:     pointer into trailing part of LU factors in Index
C  xp1:     pointer into leading part of LU factors in Value
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  app:     Ap (1..n+1) array located in Index (app...app+n)
C  axp:     Ax (1..nz) array located in Value (axp...axp+nz-1)
C  aip:     Ai (1..nz) array located in Index (aip...aip+nz-1)
C  offip:   Offi (1..nzoff) array loc. in Index (offip...offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array loc. in Value (offxp...offxp+nzoff-1)
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  ...      seven scalars in Index (lui2-6...lui2):
C  factra:  0/1 if TRANSA argument was false/true in MA38AD or MA38BD
C  nzoff:   number of entries in off-diagonal part
C  nblks:   number of diagonal blocks
C  facjob:  JOB argument in MA38AD or MA38BD when matrix factorized
C  nz:      entries in A
C  facn:    N argument in MA38AD or MA38BD when matrix factorized
C  facne:   NE argument in MA38AD or MA38BD when matrix factorized
C
C  A single diagonal block and its LU factors:
C  -------------------------------------------
C  blk:     current diagonal block
C  k1,k2:   current diagonal is A (k1..k2, k1..k2)
C  kn:      order of current diagonal block (= k2-k1+1)
C  sglton:  true if current diagonal block is 1-by-1 (a singleton)
C  luiip:   LU factors of a diagonal block start in Index (luiip...)
C  luxxp:   LU factors of a diagonal block start in Value (luxxp...)
C  npiv:    number of pivots in a diagonal block (0 <= npiv <= kn)
C  nlu:     number of elements in a diagonal block
C  lupp:    LUp (1..nlu) array located in Index (lupp...lupp+nlu-1)
C
C  An element in the LU factors of a single diagonal block:
C  --------------------------------------------------------
C  e:       element
C  luk:     number of pivots in element e
C  luip:    integer part of element is in Index (luip...)
C  luxp:    real part of element e is in Value (luxp...)
C  ludegr:  row degree (number of columns) of U2 block in element e
C  ludegc:  column degree (number of rows) of L2 block in element e
C  lunson:  number of sons of element e in the assembly DAG
C  lusonp:  list of sons of element e in Index(lusonp...lusonp+lunson-1)
C  lucp:    column pattern (row indices) of L2 block in Index (lucp..)
C  lurp:    row pattern (column indices) of U2 block in Index (lurp..)
C  nzcol:   entries in a column of L, including unit diagonal
C  nzrow:   entries in a row of U, including non-unit diagonal
C  uxp:     a row of the U2 block located in Value (uxp...)
C  son:     a son of the element e
C  ludimr:  row dimension (number of columns) in frontal matrix
C  ludimc:  column dimension (number of rows) in frontal matrix
C  maxdr:   largest ludimr for this block
C  maxdc:   largest ludimc for this block
C
C  Other:
C  ------
C  row:     row index
C  col:     column index
C  k:       kth pivot, and general loop index
C  i, j:    loop indices
C  p:       pointer
C  p1:      column of A starts Ai/Ax (p1...), or row Offi/x (p1...)
C  p2:      column of A ends in Ai/Ax (...p2), or row Offi/x (...p2)

C=======================================================================
C  EXECUTABLE STATEMENTS:
C       if (printing disabled on installation) return
C=======================================================================

C-----------------------------------------------------------------------
C  get printing control parameters
C-----------------------------------------------------------------------

        IO  = ICNTL(2)
        PRL = ICNTL(3)
        IF (PRL .LT. 3 .OR. IO .LT. 0) THEN
C          printing has not been requested
           RETURN
        ENDIF

C-----------------------------------------------------------------------
C  who is this, and where.  Determine if LU factors are to be printed
C-----------------------------------------------------------------------

        IF (WHO .EQ. 1) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'MA38AD input:'
              PRLU = .FALSE.
           ELSE
              WRITE (IO, 6) 'MA38AD output:'
              PRLU = .TRUE.
           ENDIF
        ELSE IF (WHO .EQ. 2) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'MA38BD input:'
              PRLU = .TRUE.
           ELSE
              WRITE (IO, 6) 'MA38BD output:'
              PRLU = .TRUE.
           ENDIF
        ELSE IF (WHO .EQ. 3) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'MA38CD input:'
              PRLU = .TRUE.
           ELSE
              WRITE (IO, 6) 'MA38CD output:'
              PRLU = .FALSE.
           ENDIF
        ENDIF

C-----------------------------------------------------------------------
C  print scalar input arguments: n, ne, job, trans, lvalue, lindex
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1) THEN
           WRITE (IO, 1)  'Scalar arguments:'
           WRITE (IO, 1)  '   N:         ', N, ' : order of matrix A'
           IF (WHO .EQ. 3) THEN
C             MA38CD:
C             was A or A^T factorized?
              LUI2 = KEEP (5)
              TRANSA = .FALSE.
              IF (LUI2-6 .GE. 1 .AND. LUI2-6 .LE. LINDEX) THEN
                 TRANSA = INDEX (LUI2-6) .NE. 0
              ENDIF
              TRANSC = TRANS
              IF (.NOT. TRANSC) THEN
                 IF (JOB .EQ. 1) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve P''Lx=b'
                 ELSE IF (JOB .EQ. 2) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve UQ''x=b'
                 ELSE IF (.NOT. TRANSA) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve Ax=b (PAQ=LU was factorized)'
                 ELSE
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve A''x=b (PA''Q=LU was factorized)'
                 ENDIF
              ELSE
                 IF (JOB .EQ. 1) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve L''Px=b'
                 ELSE IF (JOB .EQ. 2) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve QU''x=b'
                 ELSE IF (.NOT. TRANSA) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve A''x=b (PAQ=LU was factorized)'
                 ELSE
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve Ax=b (PA''Q=LU was factorized)'
                 ENDIF
              ENDIF
              IF (TRANSC) THEN
                 WRITE (IO, 1)
     $           '   TRANSC:          .true. : see JOB above '
              ELSE
                 WRITE (IO, 1)
     $           '   TRANSC:         .false. : see JOB above '
              ENDIF
           ELSE
C             MA38AD or MA38BD:
              WRITE (IO, 1) '   NE:        ', NE,
     $        ' : entries in matrix A'
              IF (JOB .EQ. 1) THEN
                 WRITE (IO, 1) '   JOB:       ', JOB,
     $           ' : matrix A preserved'
              ELSE
                 WRITE (IO, 1) '   JOB:       ', JOB,
     $           ' : matrix A not preserved'
              ENDIF
              TRANSA = TRANS
              IF (TRANSA) THEN
                 WRITE (IO, 1)
     $           '   TRANSA:          .true. : factorize A transpose'
              ELSE
                 WRITE (IO, 1)
     $           '   TRANSA:         .false. : factorize A'
              ENDIF
           ENDIF
           WRITE (IO, 1) '   LVALUE:    ',LVALUE,
     $     ' : size of VALUE array'
           WRITE (IO, 1) '   LINDEX:    ',LINDEX,
     $     ' : size of INDEX array'
        ENDIF

C-----------------------------------------------------------------------
C  print control parameters:  Icntl, Cntl, and Keep (6..8)
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1) THEN
           WRITE (IO, 1)
     $     'Control parameters, normally initialized by MA38ID:'
           WRITE (IO, 1) '   ICNTL (1): ', ICNTL (1),
     $     ' : I/O unit for error and warning messages'
           WRITE (IO, 1) '   ICNTL (2): ', IO,
     $     ' : I/O unit for diagnostics'
           WRITE (IO, 1) '   ICNTL (3): ', PRL,
     $     ' : printing control'
           IF (WHO .EQ. 1) THEN
              IF (ICNTL (4) .EQ. 1) THEN
                 WRITE (IO, 1) '   ICNTL (4): ', ICNTL (4),
     $           ' : use block triangular form (BTF)'
              ELSE
                 WRITE (IO, 1) '   ICNTL (4): ', ICNTL (4),
     $           ' : do not permute to block triangular form (BTF)'
              ENDIF
              WRITE (IO, 1) '   ICNTL (5): ', ICNTL (5),
     $        ' : columns examined during pivot search'
              IF (ICNTL (6) .NE. 0) THEN
                 WRITE (IO, 1) '   ICNTL (6): ', ICNTL (6),
     $           ' : preserve symmetry'
              ELSE
                 WRITE (IO, 1) '   ICNTL (6): ', ICNTL (6),
     $           ' : do not preserve symmetry'
              ENDIF
           ENDIF
           IF (WHO .NE. 3) THEN
              WRITE (IO, 1) '   ICNTL (7): ', ICNTL (7),
     $        ' : block size for dense matrix multiply'
           ELSE
              WRITE (IO, 1) '   ICNTL (8): ', ICNTL (8),
     $        ' : maximum number of iterative refinement steps'
           ENDIF
           IF (WHO .EQ. 1) THEN
              WRITE (IO, 3) '   CNTL (1): ',CNTL (1),
     $        ' : relative pivot tolerance'
              WRITE (IO, 3) '   CNTL (2): ',CNTL (2),
     $        ' : frontal matrix growth factor'
              WRITE (IO, 1) '   KEEP (7):  ',KEEP(7),
     $        ' : dense row/col control, d1'
              WRITE (IO, 1) '   KEEP (8):  ',KEEP(8),
     $        ' : dense row/col control, d2'
           ELSE IF (WHO .EQ. 3) THEN
              WRITE (IO, 3) '   CNTL (3): ',CNTL(3),
     $        ' : machine epsilon'
           ENDIF
        ENDIF

C-----------------------------------------------------------------------
C  print the informational output
C-----------------------------------------------------------------------

        IF (WHERE .NE. 1) THEN
           WRITE (IO, 1) 'Output information:'
           IF (INFO (1) .LT. 0) THEN
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : error occurred.'
           ELSE IF (INFO (1) .GT. 0) THEN
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : warning occurred'
           ELSE
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : no error or warning occurred'
           ENDIF
           IF (WHO .NE. 3) THEN
              WRITE (IO, 1) '   INFO (2):  ', INFO (2),
     $        ' : duplicate entries in A'
              WRITE (IO, 1) '   INFO (3):  ', INFO (3),
     $        ' : invalid entries in A (indices not in 1..N)'
              WRITE (IO, 1) '   INFO (4):  ', INFO (4),
     $        ' : invalid entries in A (not in prior pattern)'
              WRITE (IO, 1) '   INFO (5):  ', INFO (5),
     $        ' : entries in A after summing duplicates'
              WRITE (IO, 1)
     $  '                             and removing invalid entries'
              WRITE (IO, 1) '   INFO (6):  ', INFO (6),
     $        ' : entries in diagonal blocks of A'
              WRITE (IO, 1) '   INFO (7):  ', INFO (7),
     $        ' : entries in off-diagonal blocks of A'
              WRITE (IO, 1) '   INFO (8):  ', INFO (8),
     $        ' : 1-by-1 diagonal blocks in A'
              WRITE (IO, 1) '   INFO (9):  ', INFO (9),
     $        ' : diagonal blocks in A (>1 only if BTF used)'
              WRITE (IO, 1) '   INFO (10): ', INFO (10),
     $        ' : entries below diagonal in L'
              WRITE (IO, 1) '   INFO (11): ', INFO (11),
     $        ' : entries above diagonal in U'
              WRITE (IO, 1) '   INFO (12): ', INFO (12),
     $        ' : entries in L + U + offdiagonal blocks of A'
              WRITE (IO, 1) '   INFO (13): ', INFO (13),
     $        ' : frontal matrices'
              WRITE (IO, 1) '   INFO (14): ', INFO (14),
     $        ' : integer garbage collections'
              WRITE (IO, 1) '   INFO (15): ', INFO (15),
     $        ' : real garbage collections'
              WRITE (IO, 1) '   INFO (16): ', INFO (16),
     $        ' : diagonal pivots chosen'
              WRITE (IO, 1) '   INFO (17): ', INFO (17),
     $        ' : numerically valid pivots found in A'
              WRITE (IO, 1) '   INFO (18): ', INFO (18),
     $        ' : memory used in INDEX'
              WRITE (IO, 1) '   INFO (19): ', INFO (19),
     $        ' : minimum memory needed in INDEX'
              WRITE (IO, 1) '   INFO (20): ', INFO (20),
     $        ' : memory used in VALUE'
              WRITE (IO, 1) '   INFO (21): ', INFO (21),
     $        ' : minimum memory needed in VALUE'
              WRITE (IO, 1) '   INFO (22): ', INFO (22),
     $        ' : memory needed in INDEX for next call to MA38BD'
              WRITE (IO, 1) '   INFO (23): ', INFO (23),
     $        ' : memory needed in VALUE for next call to MA38BD'
           ELSE
              WRITE (IO, 1) '   INFO (24): ', INFO (24),
     $        ' : steps of iterative refinement taken'
           ENDIF
           IF (WHO .NE. 3) THEN
              WRITE (IO, 3) '   RINFO (1):', RINFO (1),
     $        ' : total BLAS flop count'
              WRITE (IO, 3) '   RINFO (2):', RINFO (2),
     $        ' : assembly flop count'
              WRITE (IO, 3) '   RINFO (3):', RINFO (3),
     $        ' : pivot search flop count'
              WRITE (IO, 3) '   RINFO (4):', RINFO (4),
     $        ' : Level-1 BLAS flop count'
              WRITE (IO, 3) '   RINFO (5):', RINFO (5),
     $        ' : Level-2 BLAS flop count'
              WRITE (IO, 3) '   RINFO (6):', RINFO (6),
     $        ' : Level-3 BLAS flop count'
           ELSE IF (LW .EQ. 4*N) THEN
              WRITE (IO, 3) '   RINFO (7):', RINFO (7),
     $        ' : sparse error estimate omega1'
              WRITE (IO, 3) '   RINFO (8):', RINFO (8),
     $        ' : sparse error estimate omega2'
           ENDIF
        ENDIF

C-----------------------------------------------------------------------
C  print input matrix A, in triplet form, for MA38AD and MA38BD
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1 .AND. WHO .NE. 3) THEN

           IF (TRANSA) THEN
              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1) 'The input matrix A transpose:'
                 WRITE (IO, 1) '   VALUE (1 ... ',NE,
     $           ' ): numerical values'
                 WRITE (IO, 1) '   INDEX (1 ... ',NE,
     $           ' ): column indices'
                 WRITE (IO, 1) '   INDEX (',NE+1,' ... ',2*NE,
     $           ' ): row indices'
              ENDIF
              WRITE (IO, 1)
     $        'Input matrix A transpose (entry: row, column, value):'
           ELSE
              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1) 'The input matrix A:'
                 WRITE (IO, 1) '   VALUE (1 ... ',NE,
     $           ' ): numerical values'
                 WRITE (IO, 1) '   INDEX (1 ... ',NE,
     $           ' ): row indices'
                 WRITE (IO, 1) '   INDEX (',NE+1,' ... ',2*NE,
     $           ' ): column indices'
              ENDIF
              WRITE (IO, 1)
     $        'Input matrix A (entry: row, column, value):'
           ENDIF

           PRN = MIN (PRMAX, NE)
           IF (PRL .GE. 4) THEN
              PRN = NE
           ENDIF
           DO 20 K = 1, PRN
              IF (TRANSA) THEN
                 ROW = INDEX (K+NE)
                 COL = INDEX (K)
              ELSE
                 ROW = INDEX (K)
                 COL = INDEX (K+NE)
              ENDIF
              WRITE (IO, 2) K, ROW, COL, VALUE (K)
20         CONTINUE
           IF (PRN .LT. NE) THEN
              WRITE (IO, 7)
           ENDIF
        ENDIF

C-----------------------------------------------------------------------
C  print the LU factors:  MA38AD output, MA38BD input/output,
C                         and MA38CD input
C-----------------------------------------------------------------------

        IF (PRLU .AND. INFO (1) .LT. 0) THEN
           WRITE (IO, 1)
     $     'LU factors not printed because of error flag, INFO (1) ='
     $     , INFO (1)
           PRLU = .FALSE.
        ENDIF

        IF (PRLU) THEN

C          -------------------------------------------------------------
C          description of what must be preserved between calls
C          -------------------------------------------------------------

           LUX1 = KEEP (1)
           LUX2 = KEEP (2)
           LUI1 = KEEP (3)
           LUIR1 = KEEP (4)
           LUI2 = KEEP (5)

           XP1 = LUX1
           IP1 = LUI1
           IP2 = LUI2

C          -------------------------------------------------------------
C          on input to MA38BD, only the symbol information is used
C          -------------------------------------------------------------

           SYMBOL = WHO .EQ. 2 .AND. WHERE .EQ. 1

           IF (PRL .GE. 5) THEN
              IF (SYMBOL) THEN
                 WRITE (IO, 1)
     $           'KEEP (4...5) gives the location of LU factors'
                 WRITE (IO, 1)
     $           '   which must be preserved for calls to MA38BD: '
              ELSE
                 WRITE (IO, 1)
     $           'KEEP (1...5) gives the location of LU factors'
                 WRITE (IO, 1)
     $           '   which must be preserved for calls to MA38CD: '
                 WRITE (IO, 1) '      VALUE ( KEEP (1): ', LUX1,
     $           ' ... KEEP (2): ', LUX2,' )'
                 WRITE (IO, 1) '      INDEX ( KEEP (3): ', LUI1,
     $           ' ... KEEP (5): ', LUI2,' )'
                 WRITE (IO, 1) '   and for calls to MA38BD: '
              ENDIF
              WRITE (IO, 1) '      INDEX ( KEEP (4): ',LUIR1,
     $        ' ... KEEP (5): ', LUI2,' )'
           ENDIF

           BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR.
     $        LUI2 .GT. LINDEX
           IF (.NOT. SYMBOL) THEN
              BADLU = BADLU .OR. LUX1 .LE. 0 .OR.
     $        LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE .OR. LUI1 .LE. 0 .OR.
     $        LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
           ENDIF

C          -------------------------------------------------------------
C          get the 7 scalars, and location of permutation vectors
C          -------------------------------------------------------------

           IF (BADLU) THEN
C             pointers are bad, so these values cannot be obtained
              FACNE  = 0
              FACN   = 0
              NZ     = 0
              FACJOB = 0
              NBLKS  = 0
              NZOFF  = 0
              FACTRA = 0
           ELSE
              FACNE  = INDEX (LUI2)
              FACN   = INDEX (LUI2-1)
              NZ     = INDEX (LUI2-2)
              FACJOB = INDEX (LUI2-3)
              NBLKS  = INDEX (LUI2-4)
              NZOFF  = INDEX (LUI2-5)
              FACTRA = INDEX (LUI2-6)
           ENDIF

           PRESRV = FACJOB .NE. 0
           TRANSA = FACTRA .NE. 0
           RPERMP = (LUI2-6) - (FACN)
           CPERMP = RPERMP - (FACN)
           IP2 = CPERMP - 1

           IF (PRL .GE. 5) THEN
              IF (SYMBOL) THEN
                 WRITE (IO, 1) 'Layout of LU factors in INDEX:'
              ELSE
                 WRITE (IO, 1)
     $           'Layout of LU factors in VALUE and INDEX:'
              ENDIF
           ENDIF

C          -------------------------------------------------------------
C          get location of preserved input matrix
C          -------------------------------------------------------------

           IF (PRESRV) THEN
C             preserved column-form of original matrix
              APP = IP1
              AIP = APP + (FACN+1)
              IP1 = AIP + (NZ)
              AXP = XP1
              XP1 = XP1 + (NZ)
              IF (PRL .GE. 5 .AND. .NOT. SYMBOL) THEN
                 WRITE (IO, 1)'   preserved copy of original matrix:'
                 WRITE (IO, 1)'      INDEX ( ',APP,' ... ', AIP-1,
     $           ' ): column pointers'
                 WRITE (IO, 1)'      INDEX ( ',AIP,' ... ', IP1-1,
     $           ' ): row indices'
                 WRITE (IO, 1)'      VALUE ( ',AXP,' ... ', XP1-1,
     $           ' ): numerical values'
              ENDIF
           ELSE
              IF (PRL .GE. 5 .AND. .NOT. SYMBOL) THEN
                 WRITE (IO, 1) '   original matrix not preserved.'
              ENDIF
           ENDIF

           BADLU = BADLU .OR.
     $          N .NE. FACN .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
           IF (.NOT. SYMBOL) THEN
              BADLU = BADLU .OR. XP1 .GT. LUX2 .OR. NZOFF .LT. 0
           ENDIF
           IF (BADLU) THEN
              NBLKS = 0
           ENDIF

           IF (NBLKS .LE. 1) THEN

C             ----------------------------------------------------------
C             single block (or block triangular form not used),
C             or LU factors are corrupted
C             ----------------------------------------------------------

              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1)
     $           '   collection of elements in LU factors:'
                 WRITE (IO, 1) '      INDEX ( ',LUIR1,' ... ', IP2,
     $           ' ): integer data'
                 IF (.NOT. SYMBOL) THEN
                    WRITE (IO, 1) '      VALUE ( ',XP1,' ... ', LUX2,
     $              ' ): numerical values'
                 ENDIF
              ENDIF

           ELSE

C             ----------------------------------------------------------
C             block triangular form with more than one block
C             ----------------------------------------------------------

              OFFIP = IP1
              IP1 = IP1 + (NZOFF)
              OFFXP = XP1
              XP1 = XP1 + (NZOFF)
              OFFPP = CPERMP - (N+1)
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              IP2 = LUBLPP - 1
              BADLU = BADLU .OR. LUIR1 .GT. IP2
              IF (.NOT. SYMBOL) THEN
                 BADLU = BADLU .OR. IP1 .GT. IP2 .OR.
     $           XP1 .GT. LUX2 .OR. LUIR1 .NE. IP1
              ENDIF

              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1)
     $           '   matrix permuted to upper block triangular form.'
                 IF (NZOFF .NE. 0 .AND. .NOT. SYMBOL) THEN
                    WRITE (IO, 1)'   entries not in diagonal blocks:'
                    WRITE (IO, 1)'      INDEX ( ',OFFIP,' ... ',
     $              LUIR1-1, ' ): row indices'
                    WRITE (IO, 1)'      VALUE ( ',OFFXP,' ... ',
     $              XP1-1, ' ): numerical values'
                 ENDIF
                 WRITE (IO, 1)
     $  '   collection of elements in LU factors of diagonal blocks:'
                 IF (LUIR1 .LE. LUBLPP-1) THEN
                    WRITE (IO, 1) '      INDEX ( ',LUIR1,' ... ',
     $              IP2, ' ): integer data'
                 ENDIF
                 IF (XP1 .LE. LUX2 .AND. .NOT. SYMBOL) THEN
                    WRITE (IO, 1) '      VALUE ( ',XP1,' ... ', LUX2,
     $              ' ): numerical values'
                 ENDIF
                 WRITE (IO, 1) '   other block triangular data:'
                 WRITE (IO, 1) '      INDEX ( ',LUBLPP,' ... ',
     $           BLKPP-1, ' ): pointers to block factors'
                 WRITE (IO, 1) '      INDEX ( ', BLKPP,' ... ',
     $           OFFPP-1, ' ): index range of blocks'
                 IF (.NOT. SYMBOL) THEN
                    WRITE (IO, 1) '      INDEX ( ', OFFPP,' ... ',
     $              LUI2-7,' ): off-diagonal row pointers'
                 ENDIF
              ENDIF

           ENDIF

C          -------------------------------------------------------------
C          print location of permutation vectors and 7 scalars at tail
C          -------------------------------------------------------------

           IF (PRL .GE. 5) THEN
              WRITE (IO, 1)
     $        '   permutation vectors (start at KEEP(4)-2*N-6):'
              WRITE (IO, 1) '      INDEX ( ',CPERMP,' ... ',RPERMP-1,
     $        ' ): column permutations'
              WRITE (IO, 1) '      INDEX ( ',RPERMP,' ... ',LUI2-7,
     $        ' ): row permutations'
              WRITE (IO, 1) '   other data in INDEX: '
              WRITE (IO, 1) '      INDEX ( ',LUI2-6,' ): ', FACTRA,
     $        ' : TRANSA MA38AD/MA38BD argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2-5,' ): ', NZOFF,
     $        ' : entries in off-diagonal part'
              WRITE (IO, 1) '      INDEX ( ',LUI2-4,' ): ', NBLKS,
     $        ' : number of diagonal blocks'
              WRITE (IO, 1) '      INDEX ( ',LUI2-3,' ): ', FACJOB,
     $        ' : JOB MA38AD/MA38BD argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2-2,' ): ', NZ,
     $        ' : entries in original matrix'
              WRITE (IO, 1) '      INDEX ( ',LUI2-1,' ): ', FACN,
     $        ' : N MA38AD/MA38BD argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2  ,' ): ', FACNE,
     $        ' : NE MA38AD/MA38BD argument'
           ENDIF

           IF (.NOT. SYMBOL) THEN
              BADLU = BADLU .OR. IP1 .NE. LUIR1
           ENDIF
           IP1 = LUIR1
           IF (BADLU) THEN
              WRITE (IO, 1) 'LU factors uncomputed or corrupted.'
              PRESRV = .FALSE.
              NBLKS = 0
           ENDIF

C          -------------------------------------------------------------
C          copy of original matrix in column-oriented form
C          -------------------------------------------------------------

           IF (PRESRV .AND. .NOT. SYMBOL) THEN
              WRITE (IO, 8)
              WRITE (IO, 1)
     $        'Preserved copy of original matrix (stored by column),'
              WRITE (IO, 1) 'one entry per line (row index, value):'
              DO 40 COL = 1, N
                 P1 = INDEX (APP-1 + COL)
                 P2 = INDEX (APP-1 + COL+1) - 1
                 WRITE (IO, 1) '   col: ', COL
                 IF (PRL .EQ. 3) THEN
                    P2 = MIN (PRMAX, P2)
                 ENDIF
                 DO 30 P = P1, P2
                    WRITE (IO, 5) INDEX (AIP-1 + P), VALUE (AXP-1 + P)
30               CONTINUE
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN
C                   exit out of loop if done printing:
                    WRITE (IO, 7)
                    GO TO 50
                 ENDIF
40            CONTINUE
C             loop exit label:
50            CONTINUE
           ENDIF

C          -------------------------------------------------------------
C          entries in off-diagonal blocks, in row-oriented form
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1 .AND. .NOT. SYMBOL) THEN
              WRITE (IO, 8)
              WRITE (IO, 1)
     $        'Entries not in diagonal blocks (stored by row):'
              WRITE (IO, 1) 'one entry per line (column index, value):'
              IF (NZOFF .EQ. 0) THEN
                 WRITE (IO, 1) '   (none)'
              ENDIF
              DO 70 ROW = 1, N
                 P1 = INDEX (OFFPP-1 + ROW)
                 P2 = INDEX (OFFPP-1 + ROW+1) - 1
                 IF (P2 .GE. P1) THEN
                    WRITE (IO, 1) '   row: ', ROW
                    IF (PRL .EQ. 3) THEN
                       P2 = MIN (PRMAX, P2)
                    ENDIF
                    DO 60 P = P1, P2
                       WRITE (IO, 5)
     $                 INDEX (OFFIP-1 + P), VALUE (OFFXP-1 + P)
60                  CONTINUE
                 ENDIF
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN
C                   exit out of loop if done printing:
                    WRITE (IO, 7)
                    GO TO 80
                 ENDIF
70            CONTINUE
C             loop exit label:
80            CONTINUE
           ENDIF

C          -------------------------------------------------------------
C          LU factors of each diagonal block
C          -------------------------------------------------------------

           WRITE (IO, 8)
           IF (NBLKS .GT. 0) THEN
              IF (SYMBOL) THEN
                 WRITE (IO, 1) 'Nonzero pattern of prior LU factors:'
              ELSE
                 WRITE (IO, 1) 'LU factors:'
              ENDIF
           ENDIF
           PRN = 0
           DO 200 BLK = 1, NBLKS

C             ----------------------------------------------------------
C             print the factors of a single diagonal block
C             ----------------------------------------------------------

              IF (NBLKS .GT. 1) THEN
                 K1 = INDEX (BLKPP-1 + BLK)
                 K2 = INDEX (BLKPP-1 + BLK+1) - 1
                 KN = K2-K1+1
                 SGLTON = KN .EQ. 1
                 IF (SGLTON) THEN
C                   this is a singleton
                    LUXXP = XP1-1 + INDEX (LUBLPP-1 + BLK)
                 ELSE
                    LUIIP = IP1-1 + INDEX (LUBLPP-1 + BLK)
                 ENDIF
                 IF (BLK .GT. 1) THEN
                    WRITE (IO, 9)
                 ENDIF
              ELSE
                 SGLTON = .FALSE.
                 K1 = 1
                 K2 = N
                 KN = N
                 LUIIP = IP1
              ENDIF

              IF (SGLTON) THEN

C                -------------------------------------------------------
C                this is a singleton
C                -------------------------------------------------------

                 IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
C                   exit out of loop if done printing:
                    WRITE (IO, 7)
                    GO TO 210
                 ENDIF
                 PRN = PRN + 1
                 IF (SYMBOL) THEN
                    WRITE (IO, 1) 'Block: ', BLK,
     $              ' (singleton) at index : ', K1
                 ELSE
                    WRITE (IO, 4) 'Block: ', BLK,
     $              ' (singleton) at index : ', K1,' value: ',
     $              VALUE (LUXXP)
                 ENDIF
                 IF (PRL .GE. 5) THEN
                    WRITE (IO, 1) 'located in VALUE ( ', LUXXP,' )'
                 ENDIF

              ELSE

C                -------------------------------------------------------
C                this block is larger than 1-by-1
C                -------------------------------------------------------

                 LUXXP = XP1-1 + INDEX (LUIIP)
                 NLU = INDEX (LUIIP+1)
                 NPIV = INDEX (LUIIP+2)
                 MAXDC = INDEX (LUIIP+3)
                 MAXDR = INDEX (LUIIP+4)
                 LUPP = LUIIP+5
                 IF (NBLKS .GT. 1) THEN
                    WRITE (IO, 1) 'Block: ',BLK,' first index: ',K1,
     $              ' last index: ',K2
                 ENDIF
                 IF (PRL .GE. 5) THEN
                    WRITE (IO, 1) 'elements: ', NLU, ' pivots: ', NPIV
                    WRITE (IO, 1) 'largest contribution block: ',
     $                         MAXDC, ' by ', MAXDR
                    WRITE (IO, 1)'located in INDEX ( ',LUIIP,' ... )'
                    IF (.NOT. SYMBOL) THEN
                       WRITE (IO, 1) 'and in VALUE ( ',LUXXP,' ... )'
                    ENDIF
                 ENDIF
                 LUIIP = LUPP + NLU

C                Note: the indices of the LU factors of the block range
C                from 1 to kn, even though the kn-by-kn block resides in
C                A (k1 ... k2, k1 ... k2).
                 K = 0

                 DO 190 E = 1, NLU

C                   ----------------------------------------------------
C                   print a single element
C                   ----------------------------------------------------

                    LUIP = LUIIP-1 + INDEX (LUPP-1 + E)
                    LUXP = LUXXP-1 + INDEX (LUIP)
                    LUK  = INDEX (LUIP+1)
                    LUDEGR = INDEX (LUIP+2)
                    LUDEGC = INDEX (LUIP+3)
                    LUNSON = INDEX (LUIP+4)
                    LUDIMR = INDEX (LUIP+5)
                    LUDIMC = INDEX (LUIP+6)
                    LUCP = LUIP + 7
                    LURP = LUCP + LUDEGC
                    LUSONP = LURP + LUDEGR
                    IF (PRL .GE. 5) THEN
                       WRITE (IO, 1) '   e: ', E, ' pivots: ', LUK
                       WRITE (IO, 1) '   children in dag: ', LUNSON,
     $                 ' frontal matrix: ', LUDIMR, ' by ', LUDIMC
                       ENDIF

C                   ----------------------------------------------------
C                   print the columns of L
C                   ----------------------------------------------------

                    P = LUXP
                    DO 140 J = 1, LUK
                       COL = K+J
                       NZCOL = LUK-J+1+LUDEGC
                       WRITE (IO, 1) '      L, col: ', COL
                       PRN = PRN + 1
                       ROW = COL
                       IF (SYMBOL) THEN
                          WRITE (IO, 5) ROW
                       ELSE
C                         L is unit diagonal:
                          WRITE (IO, 5) ROW, ONE
                       ENDIF
                       P = P + 1
C                      pivot block
                       DO 120 I = J+1, LUK
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          ROW = K+I
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) ROW
                          ELSE
                             WRITE (IO, 5) ROW, VALUE (P)
                          ENDIF
                          P = P + 1
120                    CONTINUE
C                      L block
                       DO 130 I = 1, LUDEGC
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          ROW = INDEX (LUCP-1+I)
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) ROW
                          ELSE
                             WRITE (IO, 5) ROW, VALUE (P)
                          ENDIF
                          P = P + 1
130                    CONTINUE
                       P = P + J
140                 CONTINUE

C                   ----------------------------------------------------
C                   print the rows of U
C                   ----------------------------------------------------

                    UXP = LUXP + LUK*(LUDEGC+LUK)
                    DO 170 I = 1, LUK
                       ROW = K+I
                       NZROW = LUK-I+1+LUDEGR
                       WRITE (IO, 1) '      U, row: ', ROW
                       P = LUXP + (I-1) + (I-1) * (LUDEGC+LUK)
C                      pivot block
                       DO 150 J = I, LUK
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          COL = K+J
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) COL
                          ELSE
                             WRITE (IO, 5) COL, VALUE (P)
                          ENDIF
                          P = P + (LUDEGC+LUK)
150                    CONTINUE
                       P = UXP
C                      U block
                       DO 160 J = 1, LUDEGR
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          COL = INDEX (LURP-1+J)
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) COL
                          ELSE
                             WRITE (IO, 5) COL, VALUE (P)
                          ENDIF
                          P = P + LUK
160                    CONTINUE
                       UXP = UXP + 1
170                 CONTINUE

C                   ----------------------------------------------------
C                   print the sons of the element in the assembly DAG
C                   ----------------------------------------------------

                    IF (PRL .GE. 5) THEN
                       DO 180 I = 1, LUNSON
                          PRN = PRN + 1
                          SON = INDEX (LUSONP-1+I)
                          IF (SON .LE. KN) THEN
C                            an LUson
                             WRITE (IO, 1) '      LUson: ', SON
                          ELSE IF (SON .LE. 2*KN) THEN
C                            a Uson
                             WRITE (IO, 1) '      Uson:  ', SON-KN
                          ELSE
C                            an Lson
                             WRITE (IO, 1) '      Lson:  ', SON-2*KN
                          ENDIF
180                    CONTINUE
                    ENDIF

C                   ----------------------------------------------------
C                   increment count of pivots within this block
C                   ----------------------------------------------------

                    K = K + LUK
190              CONTINUE
              ENDIF
200        CONTINUE
C          loop exit label:
210        CONTINUE

C          -------------------------------------------------------------
C          row and column permutations
C          -------------------------------------------------------------

           IF (.NOT. BADLU) THEN
              PRN = MIN (PRMAX, N)
              IF (PRL .GE. 4) THEN
C                print all of Cperm and Rperm
                 PRN = N
              ENDIF
              WRITE (IO, 8)
              WRITE (IO, 1) 'Column permutations'
              DO 220 I = 1, PRN
                 WRITE (IO, 5) INDEX (CPERMP+I-1)
220           CONTINUE
              IF (PRN .LT. N) THEN
                 WRITE (IO, 7)
              ENDIF
              WRITE (IO, 8)
              WRITE (IO, 1) 'Row permutations'
              DO 230 I = 1, PRN
                 WRITE (IO, 5) INDEX (RPERMP+I-1)
230           CONTINUE
              IF (PRN .LT. N) THEN
                 WRITE (IO, 7)
              ENDIF
           ENDIF

        ENDIF

C-----------------------------------------------------------------------
C  print B (on input) or W and X (on output) for MA38CD
C-----------------------------------------------------------------------

        IF (WHO .EQ. 3) THEN
           WRITE (IO, 8)
           PRN = MIN (PRMAX, N)
           IF (PRL .GE. 4) THEN
C             print all of B, or W and X
              PRN = N
           ENDIF
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 1) 'W (1 ... ',LW,
     $        ' ), work vector: not printed'
              WRITE (IO, 1) 'B (1 ... ',N,' ), right-hand side: '
              DO 240 I = 1, PRN
                 WRITE (IO, 5) I, B (I)
240           CONTINUE
              IF (PRN .LT. N) THEN
                 WRITE (IO, 7)
              ENDIF
           ELSE
              IF (INFO (1) .LT. 0) THEN
                 WRITE (IO, 1) 'W (1 ... ',LW,' ), work vector, and'
                 WRITE (IO, 1) 'X (1 ... ',N, ' ), solution,'
                 WRITE (IO, 1)
     $           '   not printed because of error flag, INFO (1) = ',
     $           INFO (1)
              ELSE
                 IF (LW .EQ. 4*N) THEN
C                   MA38CD did iterative refinement
                    WRITE (IO, 1) 'W (1 ... ',N,' ), residual: '
                    DO 250 I = 1, PRN
                       WRITE (IO, 5) I, W (I)
250                 CONTINUE
                    IF (PRN .LT. N) THEN
                       WRITE (IO, 7)
                    ENDIF
                    WRITE (IO, 1) 'W (',N+1,' ... ',LW,
     $              ' ), work vector: not printed'
                 ELSE
C                   no iterative refinement
                    WRITE (IO, 1) 'W (1 ... ',LW,
     $              ' ), work vector: not printed'
                 ENDIF
                 WRITE (IO, 1) 'X (1 ... ',N,' ), solution: '
                 DO 260 I = 1, PRN
                    WRITE (IO, 5) I, X (I)
260              CONTINUE
                 IF (PRN .LT. N) THEN
                    WRITE (IO, 7)
                 ENDIF
              ENDIF
           ENDIF
        ENDIF

C-----------------------------------------------------------------------
C  who is this, and where:
C-----------------------------------------------------------------------

        IF (WHO .EQ. 1) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'end of MA38AD input '
           ELSE
              WRITE (IO, 6) 'end of MA38AD output'
           ENDIF
        ELSE IF (WHO .EQ. 2) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'end of MA38BD input '
           ELSE
              WRITE (IO, 6) 'end of MA38BD output'
           ENDIF
        ELSE IF (WHO .EQ. 3) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'end of MA38CD input '
           ELSE
              WRITE (IO, 6) 'end of MA38CD output'
           ENDIF
        ENDIF

        RETURN

C=======================================================================
C  FORMAT STATMENTS
C=======================================================================

1       FORMAT (' ', A, :, I12, :, A, :, I12, :,
     $               A, :, I12, :, A, :, I12, :, A, :, I12)
2       FORMAT (' ', I12, ': ', I12, ' ', I12, ' ', D11.4)
3       FORMAT (' ', A, D13.4, A)
4       FORMAT (' ', A, I12, A, I12, A, D11.4)
5       FORMAT (' ', I12, :, ': ', D11.4)
6       FORMAT (' ', 59('='), A)
7       FORMAT ('    ...')
8       FORMAT (' ', 79 ('-'))
9       FORMAT (' ', 79 ('.'))
        END

        SUBROUTINE MA38ZD (WHO, ERROR, I, J, X, IO)
        INTEGER WHO, ERROR, I, J, IO
        DOUBLE PRECISION X

C=== MA38ZD ============================================================
C
C  Unsymmetric-pattern multifrontal package (MA38). Double-precision.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  October 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C=======================================================================
C  NOT USER-CALLABLE

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Print error and warning messages for for MA38AD, MA38BD, and MA38CD.

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be deleted on installation (replaced with a dummy
C  routine that just returns without printing) in order to completely
C  disable the printing of all error and warning messages.  The error
C  and warning return flag (Info (1)) will not be affected.  To
C  completely disable all I/O, you can also replace the MA38YD routine
C  with a dummy subroutine.  If you make this modification, please do
C  not delete any original code - just comment it out instead.  Add a
C  comment and date to your modifications.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who:            what user-callable routine called MA38ZD:
C                       1: MA38AD, 2: MA38BD, 3: MA38CD
C       i, j, x:        the relevant offending value(s)
C       io:             I/O unit on which to print.  No printing
C                       occurs if < 0.
C       error:          the applicable error (<0) or warning (>0)
C                       Errors (<0) cause the factorization/solve to
C                       be terminated.  If an error occurs, a prior
C                       warning status is overwritten with the error
C                       status.
C
C  The following error codes are returned in Info (1) by MA38ND.
C  These errors cause the factorization or solve to terminate:
C
C  Where**      Error   Description
C
C  FA RF  -     -1      N < 1 or N > maximum value
C  FA RF  -     -2      NE < 1 or NE > maximum value
C  FA RF  -     -3      LINDEX too small
C  FA RF  -     -4      LVALUE too small
C  FA RF  -     -5      both LINDEX and LVALUE are too small
C   - RF  -     -6      prior pivot ordering no longer acceptable
C   - RF SO     -7      LU factors are uncomputed, or are corrupted
C
C  The following warning codes are returned in Info (1) by MA38ND.
C  The factorization or solve was able to complete:
C
C  FA RF  -     1       invalid entries
C  FA RF  -     2       duplicate entries
C  FA RF  -     3       invalid and duplicate entries
C  FA RF  -     4       singular matrix
C  FA RF  -     5       invalid entries, singular matrix
C  FA RF  -     6       duplicate entries, singular matrix
C  FA RF  -     7       invalid and duplicate entries, singular matrix
C   -  - SO     8       iterative refinement cannot be done
C
C  The following are internal error codes (not returned in Info (1))
C  for printing specific invalid or duplicate entries.  These codes are
C  for MA38KD, MA38MD, and MA38RD.  Invalid entries are ignored, and
C  duplicate entries are added together (and the factorization
C  continues).  Warning levels (1..7) will be set later by MA38ND,
C  above.
C
C  FA RF  -     99      invalid entry, out of range 1..N
C  FA RF  -     98      duplicate entry
C   - RF  -     97      invalid entry:  within a diagonal block, but not
C                       in the pattern of the LU factors of that block.
C   - RF  -     96      invalid entry:  below the diagonal blocks.  Can
C                       only occur if the matrix has been ordered into
C                       block-upper-triangular form.
C   - RF  -     95      invalid entry:  matrix is singular.  The
C                       remaining rank 0 submatrix yet to be factorized
C                       is replaced with the identity matrix in the LU
C                       factors.  Any entry that remains is ignored.

C ** FA: MA38AD, RF: MA38BD, SO: MA38CD

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C  Error or warning message printed on I/O unit

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  MA38ND, MA38KD, MA38PD, MA38RD

C=======================================================================
C  EXECUTABLE STATEMENTS:
C       if (printing disabled on installation) return
C=======================================================================

        IF (IO .LT. 0) THEN
C          printing of error / warning messages has not been requested
           RETURN
        ENDIF

        IF (WHO .EQ. 1) THEN

C          -------------------------------------------------------------
C          MA38AD error messages
C          -------------------------------------------------------------

           IF (ERROR .EQ. -1) THEN
              WRITE (IO, 1) 'MA38AD: N less than one.'
           ELSE IF (ERROR .EQ. -2) THEN
              WRITE (IO, 1) 'MA38AD: NE less than one.'
           ELSE IF (ERROR .EQ. -3) THEN
              WRITE (IO, 1)
     $        'MA38AD: LINDEX too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. -4) THEN
              WRITE (IO, 1)
     $        'MA38AD: LVALUE too small.  Must be greater than ', I

C          -------------------------------------------------------------
C          MA38AD cumulative warning messages
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 1) THEN
              WRITE (IO, 1) 'MA38AD: ', I,
     $        ' invalid entries ignored (out of range 1..N).'
           ELSE IF (ERROR .EQ. 2) THEN
              WRITE (IO, 1) 'MA38AD: ', I,' duplicate entries summed.'
           ELSE IF (ERROR .EQ. 4) THEN
              WRITE (IO, 1)
     $        'MA38AD: matrix is singular.  Only ', I, ' pivots found.'

C          -------------------------------------------------------------
C          MA38AD non-cumulative warning messages (internal error codes)
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 99) THEN
              WRITE (IO, 2)
     $        'MA38AD: invalid entry (out of range 1..N):', I, J, X
           ELSE IF (ERROR .EQ. 98) THEN
              WRITE (IO, 2)
     $        'MA38AD: duplicate entry summed:', I, J, X
           ENDIF

        ELSE IF (WHO .EQ. 2) THEN

C          -------------------------------------------------------------
C          MA38BD error messages
C          -------------------------------------------------------------

           IF (ERROR .EQ. -1) THEN
              WRITE (IO, 1) 'MA38BD: N less than one.'
           ELSE IF (ERROR .EQ. -2) THEN
              IF (I .LT. 0) THEN
                 WRITE (IO, 1) 'MA38BD: NE less than one.'
              ELSE
                 WRITE (IO, 1)
     $           'MA38BD: NE too large.  Must be less than ', I
              ENDIF
           ELSE IF (ERROR .EQ. -3) THEN
              WRITE (IO, 1)
     $        'MA38BD: LINDEX too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. -4) THEN
              WRITE (IO, 1)
     $        'MA38BD: LVALUE too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. -6) THEN
              WRITE (IO, 1) 'MA38BD: pivot order from MA38AD failed.'
           ELSE IF (ERROR .EQ. -7) THEN
              WRITE (IO, 1)
     $        'MA38BD: LU factors uncomputed or corrupted.'

C          -------------------------------------------------------------
C          MA38BD cumulative warning messages
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 1) THEN
              IF (I .GT. 0) THEN
                 WRITE (IO, 1) 'MA38BD: ', I,
     $           ' invalid entries ignored (out of range 1..N).'
              ELSE
                 WRITE (IO, 1) 'MA38BD: ',-I,
     $           ' invalid entries ignored (not in prior pattern).'
              ENDIF
           ELSE IF (ERROR .EQ. 2) THEN
              WRITE (IO, 1) 'MA38BD: ', I,' duplicate entries summed.'
           ELSE IF (ERROR .EQ. 4) THEN
              WRITE (IO, 1) 'MA38BD: matrix is singular.  Only ', I,
     $        ' pivots found.'

C          -------------------------------------------------------------
C          MA38BD non-cumulative warning messages (internal error codes)
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 99) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (out of range 1..N):', I, J, X
           ELSE IF (ERROR .EQ. 98) THEN
              WRITE (IO, 2)
     $        'MA38BD: duplicate entry summed:', I, J, X
           ELSE IF (ERROR .EQ. 97) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (not in pattern of prior factors)',
     $        I, J, X
           ELSE IF (ERROR .EQ. 96) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (below diagonal blocks):', I, J, X
           ELSE IF (ERROR .EQ. 95) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (prior matrix singular):', I, J, X
           ENDIF

        ELSE IF (WHO .EQ. 3) THEN

C          -------------------------------------------------------------
C          MA38CD error messages
C          -------------------------------------------------------------

           IF (ERROR .EQ. -7) THEN
              WRITE (IO, 1)
     $        'MA38CD: LU factors uncomputed or corrupted.'

C          -------------------------------------------------------------
C          MA38CD non-cumulative warning messages
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 8) THEN
              IF (I .EQ. 0) THEN
                 WRITE (IO, 1)
     $  'MA38CD: no iterative refinement: original matrix not preserved'
              ELSE
                 WRITE (IO, 1)
     $  'MA38CD: no iterative refinement: only for Ax=b or A''x=b'
              ENDIF
           ENDIF

        ENDIF

        RETURN

C=======================================================================
C  FORMAT STATMENTS
C=======================================================================

1       FORMAT (' ', A, :, I12, :, A)
2       FORMAT (' ', A, /,
     $          '    row: ', I12, ' col: ', I12,' value: ', D11.4)
        END


