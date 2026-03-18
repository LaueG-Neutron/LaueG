C ***************** This file must be compiled with R8 as the default *************

      SUBROUTINE NLSCON(N,M,MFIT,FCN,JAC,X,XSCAL,FI,FSCAL,RTOL,IOPT,
     $IERR,LIWK,IWK,LRWK,RWK)
C*    Begin Prologue NLSCON
      INTEGER N,M,MFIT
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N),FI(MFIT),FSCAL(MFIT)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*  Title
C
C     Numerical solution of nonlinear (NL) least squares (S)
C     problems with nonlinear constraints (CON), especially
C     designed for numerically sensitive problems.
C
C*  Written by        U. Nowak, L. Weimann 
C*  Purpose           Solution of highly nonlinear, optionally 
C                     constrained least squares problems
C*  Method            Damped affine invariant Gauss-Newton method
C                     (see references below)
C*  Category          K1b2b. - Nonlinear least squares approxi-
C                     mation with nonlinear constraints
C*  Keywords          Nonlinear least squares problems, 
C                     Gauss-Newton methods
C*  Version           2.3.2
C*  Revision          December 1993
C*  Latest Change     July 2000
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Lutz Weimann
C                     ZIB, Division Scientific Computing, 
C                          Department Scientific Software
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: weimann@zib.de
C
C*    References:
C
C     /1/ P. Deuflhard:
C         Newton Techniques for Highly Nonlinear Problems -
C         Theory and Algorithms.
C         Academic press Inc. (To be published)
C
C     /2/ U. Nowak, L. Weimann:
C         A Family of Newton Codes for Systems of Highly Nonlinear
C         Equations - Algorithm, Implementation, Application.
C         ZIB, Technical Report TR 90-10 (December 1990)
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from acquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Damped Gauss-Newton-algorithm with rank strategy for highly 
C     nonlinear least squares approximation problems (optionally
C     constrained, over- and underdetermined problems) - 
C     due to Ref.(1).
C
C     (The iteration is done by subroutine NCINT currently. NLSCON
C      itself does some house keeping and builds up workspace.)
C
C     The problem solved by this program looks as follows:
C
C     Denoting below the n-dimensional real space with IR(n),
C     the number of parameters to be estimated with N,
C     the number of measurements (data to be fitted) with MFIT, and
C     the number of equality constraints with MCON,
C     M := MCON + MFIT ,
C     let   F : IR(N) --> IR(MFIT) ,  G : IR(N) --> IR(MCON) 
C     be nonlinear functions,
C     FI     in IR(MFIT) the vector of measurement data and 
C     FSCAL  in IR(MFIT) the vector of measurement weights.
C
C     For M >= N, find a parameter vector X in IR(N), which 
C     minimizes  Sum (j=1,...,MFIT) ((Fj(X)-FIj)/FSCALj)**2 and
C     satisfies   G(X) = 0.
C     For M < N, find the parameter vector X in IR(N) with the 
C     smallest possible euclidian norm, 
C     which satisfies  F(X) = 0  and  G(X) = 0  . 
C
C     Jacobian approximation by numerical differences or user
C     supplied subroutine JAC.
C
C     The numerical solution of the arising linear least squares
C     problem is done by means of the subroutines DECCON and SOLCON
C     (QR decomposition with subcondition estimation, rank decision
C     and computation of the rank-deficient pseudoinverse) .
C     For special purposes these routines may be substituted.
C
C     A statistical a posteriori analysis of the parameter estimate
C     is optionally available.
C
C     This is a driver routine for the core solver NCINT.
C
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    External subroutines (to be supplied by the user)
C     =================================================
C 
C     (Caution: Arguments declared as (input) must not
C               be altered by the user subroutines ! )
C
C     FCN(N,M,MCON,X,F,IFAIL) 
C                      Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       M              Int    Number of measurement-vector
C                             components plus number of equality
C                             constraints (input)
C       MCON           Int    Number of equality constraints (input)
C       X(N)           Dble   Vector of parameters (input)
C       F(M)           Dble   Vector of equality constraints and
C                             measurement fitting values -
C                             the first MCON-components belonging
C                             to the constraints (output).
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             On input:  Has always value 0 (zero).
C                             On output: Indicates failure of FCN eval-
C                                uation, if having a value <= 2.
C                             If <0: NLSCON will be terminated with 
C                                    error code = 82, and IFAIL stored
C                                    to IWK(23).
C                             If =1: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced to it's half.
C                             If =2: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced by a reduct. factor, which 
C                                    must be output through F(1) by FCN,
C                                    and it's value must be >0 and < 1.
C                             Note, that if IFAIL = 1 or 2, additional
C                             conditions concerning the damping factor,
C                             e.g. the minimum damping factor or the
C                             bounded damping strategy may also influ-
C                             ence the value of the reduced damping 
C                             factor.
C
C
C     JAC(N,M,MCON,X,DFDX,IFAIL)
C                        Ext    Jacobian matrix subroutine
C       N                  Int    See parameter N of FCN above (input)
C       M                  Int    See parameter M of FCN above (input)
C       MCON               Int    See parameter MCON of FCN above 
C                                 (input)
C       X(N)               Dble   See parameter X of FCN above (input)
C       DFDX(M,N)          Dble   DFDX(i,k): partial derivative of
C                                 I-th component of FCN with respect
C                                 to X(k) (output)
C       IFAIL              Int    JAC evaluation-failure indicator. 
C                                 (output)
C                                 Has always value 0 (zero) on input.
C                                 Indicates failure of JAC evaluation
C                                 and causes termination of NLSCON,
C                                 if set to a negative value on output
C
C
C*    Input parameters of NLSCON
C     ==========================
C
C     N              Int    Number of parameters to be estimated
C     M              Int    Sum of number of measurement data and
C                           equality constraints
C     MFIT           Int    Number of measurement data (to be fitted)
C   * X(N)           Dble   Initial estimate of parameters
C   * XSCAL(N)       Dble   User scaling (lower threshold) of the 
C                           iteration vector X(N)
C     FI(MFIT)       Dble   Data obtained by measurements
C     FSCAL(MFIT)    Dble   User weighting vector of measurements
C   * RTOL           Dble   Required relative precision of
C                           solution components -
C                           RTOL.GE.EPMACH*TEN*N
C   * IOPT(50)       Int    Array of run-time options. Set to zero
C                           to get default values (details see below)
C
C*    Output parameters of NLSCON
C     ===========================
C
C   * X(N)           Dble   Solution parameters ( or final values,
C                           respectively )
C   * XSCAL(N)       Dble   After return with IERR.GE.0, it contains
C                           the latest internal scaling vector used
C                           After return with IERR.EQ.-1 in onestep-
C                           mode it contains a possibly adapted 
C                           (as described below) user scaling vector:
C                           If (XSCAL(I).LT. SMALL) XSCAL(I) = SMALL ,
C                           If (XSCAL(I).GT. GREAT) XSCAL(I) = GREAT .
C                           For SMALL and GREAT, see section machine
C                           constants below  and regard note 1.
C   * RTOL           Dble   Finally achieved (relative) accuracy
C                           The estimated absolute error of component i
C                           of x_out is approximately given by
C                             abs_err(i) = RTOL * XSCAL_out(i) ,
C                           where (approximately)
C                             XSCAL_out(i) = 
C                                max(abs(X_out(i)),XSCAL_in(i)).
C                           Note that RTOL_out may be greater than
C                           RTOL_in, but NLSCON claims 'solution found'
C                           - see IOPT(36).
C     IERR           Int    Return value parameter
C                           =-1 sucessfull completion of one iteration
C                               step, subsequent iterations are needed 
C                               to get a solution. (stepwise mode only)
C                           = 0 successfull completion of iteration
C                           > 0 see list of error messages below
C
C     Note 1.
C        The machine dependent values SMALL, GREAT and EPMACH are
C        gained from calls of the machine constants function D1MACH.
C        As delivered, this function is adapted to use constants 
C        suitable for all machines with IEEE arithmetic. If you use
C        another type of machine, you have to change the DATA state-
C        ments for IEEE arithmetic in D1MACH into comments and to 
C        uncomment the set of DATA statements suitable for your machine.
C
C*    Workspace parameters of NLSCON
C     ==============================
C
C     LIWK           Int    Declared dimension of integer
C                           workspace.
C                           Required minimum (for standard linear least
C                           squares problem solver): N+52
C   * IWK(LIWK)      Int    Integer Workspace.
C     LRWK           Int    Declared dimension of real workspace.
C                           Required minimum (for standard linear least
C                           squares problem solver and Jacobian computed 
C                           by numerical approximation - if the Jacobian
C                           is computed by a user subroutine JAC,
C                           decrease the expression noted below by N):
C                           (2*M+N)*N+8*M+10*N+MAX(M,N)+61
C   * RWK(LRWK)      Dble   Real Workspace
C
C     Note 2a.  A test on sufficient workspace is made. If this
C               test fails, IERR is set to 10 and an error-message
C               is issued from which the minimum of required
C               workspace size can be obtained.
C
C     Note 2b.  The first 50 elements of IWK and RWK are partially 
C               used as input for internal algorithm parameters (for
C               details, see below). In order to set the default values
C               of these parameters, the fields must be set to zero.
C               Therefore, it's recommanded always to initialize the
C               first 50 elements of both workspaces to zero.
C
C*   Options IOPT:
C    =============
C
C     Pos. Name   Default  Meaning
C
C       1  QSUCC  0        =0 (.FALSE.) initial call:
C                             NLSCON is not yet initialized, i.e. this
C                             is the first call for this nonlinear
C                             least squares problem.
C                             At successfull return with MODE=1,
C                             QSUCC is set to 1.
C                          =1 (.TRUE.) successive call:
C                             NLSCON is initialized already and is now
C                             called to perform one or more following
C                             Gauss-Newton-iteration steps.
C                             ATTENTION:
C                                Don't destroy the contents of
C                                IOPT(i) for 1 <= i <= 50 ,
C                                IWK(j)  for 1 <= j < NIWKFR and
C                                RWK(k)  for 1 <= k < NRWKFR.
C                                (Nevertheless, some of the options, e.g.
C                                 FCMIN, SIGMA, MPR..., can be modified
C                                 before successive calls.)
C       2  MODE   0        =0 Standard mode initial call:
C                             Return when the required accuracy for the
C                             iteration vector is reached. User defined
C                             parameters are evaluated and checked.
C                             Standard mode successive call:
C                             If NLSCON was called previously with
C                             MODE=1, it performs all remaining 
C                             iteration steps.
C                          =1 Stepwise mode:
C                             Return after one Gauss-Newton
C                             iteration step.
C       3  JACGEN 0        Method of Jacobian generation
C                          =0 Standard method is JACGEN=2
C                          =1 User supplied subroutine JAC will be 
C                             called to generate Jacobian matrix
C                          =2 Jacobian approximation by numerical
C                             differentation (see subroutine NCJAC)
C                          =3 Jacobian approximation by numerical
C                             differentation with feedback control
C                             (see subroutine NCJCF)
C       4..8               Reserved
C       9  ISCAL  0        Determines how to scale the iterate-vector:
C                          =0 The user supplied scaling vector XSCAL is
C                             used as a (componentwise) lower threshold
C                             of the current scaling vector
C                          =1 The vector XSCAL is always used as the
C                             current scaling vector
C      10                  Reserved
C      11  MPRERR 0        Print error messages
C                          =0 No output
C                          =1 Error messages
C                          =2 Warnings additionally
C                          =3 Informal messages additionally
C      12  LUERR  6        Logical unit number for error messages
C      13  MPRMON 0        Print iteration Monitor
C                          =0 No output
C                          =1 Standard output
C                          =2 Summary iteration monitor additionally
C                          =3 Detailed iteration monitor additionally
C                          =4,5,6 Outputs with increasing level addi-
C                             tional increasing information for code
C                             testing purposes. Level 6 produces
C                             in general extremely large output!
C      14  LUMON  6        Logical unit number for iteration monitor
C      15  MPRSOL 0        Print solutions
C                          =0 No output
C                          =1 Initial values and solution values
C                          =2 Intermediate iterates additionally
C      16  LUSOL  6        Logical unit number for solutions
C      17..18              Reserved
C      19  MPRTIM 0        Output level for the time monitor
C                          = 0 : no time measurement and no output
C                          = 1 : time measurement will be done and
C                                summary output will be written -
C                                regard note 4a.
C      20  LUTIM  6        Logical output unit for time monitor
C      21  QSTAT  0        Statistical Analysis of the final 
C                          least squares estimate:
C                          = 0 : Analysis will not be done
C                          = 1 : Analysis will be done,
C                                and certain results are stored
C                                to the RWK array (for details, see
C                                RWK description below)
C      22  MPRSTA 0        Printing of statistical Analysis for
C                          the final least squares estimate:
C                          = 0 : Printing will not be done
C                          = 1 : Printing will be done (and this
C                                implies QSTAT to be set to 1)
C      23..30              Reserved
C      31  NONLIN 3        Problem type specification
C                          =1 Linear problem
C                             Warning: If specified, no check will be
C                             done, if the problem is really linear, and
C                             NLSCON terminates unconditionally after
C                             one Gauss-Newton-iteration step.
C                          =2 Mildly nonlinear problem
C                          =3 Highly nonlinear problem
C                          =4 Extremely nonlinear problem
C      32  QRANK1 0        =0 (.FALSE.) Rank-1 updates by Broyden-
C                             approximation are inhibited.
C                          =1 (.TRUE.) Rank-1 updates by Broyden-
C                             approximation are allowed.
C      33..34              Reserved
C      35  QNSCAL 0        Inhibit automatic row scaling: 
C                          =0 (.FALSE.) Automatic row scaling of
C                             the linear system is activ: 
C                             Rows i=1,...,N will be divided by
C                             max j=1,...,N (abs(a(i,j))) 
C                          =1 (.TRUE.) No row scaling of the linear
C                             system. Recommended only for well row-
C                             scaled nonlinear least squares problems.
C      36  ITERM  0        Determines the iteration termination cri-
C                          terium to be chosen:
C                          =0 Iteration is terminated, if one of the
C                             stopping criteria used for ITERM=1 and
C                             ITERM=2 is satisfied (see below).
C                          =1 Iteration is terminated, if, from a
C                             statistical point of view, a resonable
C                             precision is achieved, i.e.
C                             simplified GN-correction < RTOL, and an
C                             estimate of the accuracy is available.
C                             Recommended to be used for incompatible
C                             problems.
C                          =2 Iteration is terminated, if 
C                             GN-correction < RTOL . Using this option
C                             may force 'to much' precision from
C                             the statistical point of view.
C      37                  Reserved
C      38  IBDAMP          Bounded damping strategy switch:
C                          =0 means currently always IBDAMP = off 
C                             (but may depend on the settings of other
C                              options in future versions)
C                          =1 means always IBDAMP = on 
C                          =2 means always IBDAMP = off 
C      39..45              Reserved
C      46..50              User options (see note 4b)

C     Note 3:
C         If NLSCON terminates with IERR=2 (maximum iterations)
C         or  IERR=3 (small damping factor), you may try to continue
C         the iteration by increasing NITMAX or decreasing FCMIN
C         (see RWK) and setting QSUCC to 1.
C
C     Note 4a:
C        The integrated time monitor calls the machine dependent
C        subroutine SECOND to get the current time stamp in form
C        of a real number (Single precision). As delivered, this
C        subroutine always return 0.0 as time stamp value. Refer
C        to the compiler- or library manual of the FORTRAN compiler
C        which you currently use to find out how to get the current
C        time stamp on your machine.
C
C     Note 4b:
C         The user options may be interpreted by the user replacable
C         routines NCSOUT, NCFACT, NCSOLV - the distributed version
C         of NCSOUT currently uses IOPT(46) as follows:
C         0 = standard plotdata output (may be postprocessed by a user-
C             written graphical program)
C         1 = plotdata output is suitable as input to the graphical
C             package GRAZIL (based on GKS), which has been developed
C             at ZIB. 
C
C
C*   Optional INTEGER input/output in IWK:
C    =======================================
C
C     Pos. Name          Meaning
C
C      1   NITER  IN/OUT Number of Gauss-Newton-iterations
C      2                 reserved
C      3   NCORR  IN/OUT Number of corrector steps
C      4   NFCN   IN/OUT Number of FCN-evaluations
C      5   NJAC   IN/OUT Number of Jacobian generations or
C                        JAC-calls
C      6                 reserved
C      7                 reserved
C      8   NFCNJ  IN/OUT Number of FCN-evaluations for Jacobian
C                        approximation
C      9   NREJR1 IN/OUT Number of rejected Gauss-Newton iteration
C                        steps done with a rank-1 computed Jacobian
C     10..11             Reserved
C     12   IDCODE IN/OUT Output: The 8 decimal digits program identi-
C                        fication number ppppvvvv, consisting of the
C                        program code pppp and the version code vvvv.
C                        Input: If containing a negative number,
C                        it will only be overwritten by the identi-
C                        fication number, immediately followed by
C                        a return to the calling program.      
C     13..15             Reserved
C     16   NIWKFR OUT    First element of IWK which is free to be used
C                        as workspace between Gauss-Newton iterations.
C                        For standard linear solvers:  51
C     17   NRWKFR OUT    First element of RWK which is free to be used
C                        as workspace between Newton iteration steps.
C                        For standard linear solver and numerically 
C                        approximated Jacobian computed by the 
C                        expression: N*(M+4)+4*M+61 
C                        If the Jacobian is computed by a user routine
C                        JAC, subtract N in this expression.
C     18   LIWKA  OUT    Length of IWK currently required
C     19   LRWKA  OUT    Length of RWK currently required
C     20..22             Reserved
C     23   IFAIL  OUT    Set in case of failure of NCFACT (IERR=80),
C                        N2SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
C                        to the nonzero IFAIL value returned by the 
C                        routine indicating the failure .
C     24..30             Reserved
C     31   NITMAX IN     Maximum number of permitted iteration
C                        steps (default: 50)
C     32   IRANK  IN     Initial rank
C                        =0 : full rank N
C                        =k with 0 < k < N : deficient rank assumed
C                           for the Jacobian in the starting point
C     33   NEW    IN/OUT Count of consecutive rank-1 updates
C     34   IFCCNT INTERN Count of consecutive special steps before 
C                        convergence stop of iteration
C     35..50             Reserved
C
C*   Optional REAL input/output in RWK:
C    ====================================
C
C     Pos. Name          Meaning
C
C      1..16             Reserved
C     17   CONV   OUT    The maximum norm of the latest ordinary 
C                        respective simplified (scaled) Gauss-Newton
C                        correction.
C     18   SUMX   OUT    Natural level (not Normx of printouts)
C                        of the current iterate, i.e. Sum(DX(i)**2),
C                        where DX = scaled Gauss-Newton correction.
C     19   DLEVF  OUT    Standard level (not Normf of printouts)
C                        of the current iterate, i.e. Norm2(F(X)),
C                        where F =  nonlinear model function.
C     20   FCBND  IN     Bounded damping strategy restriction factor
C                        (Default is 10)
C     21   FCSTRT IN     Damping factor for first Gauss-Newton iteration
C                        - overrides option NONLIN, if set (see note 5)
C     22   FCMIN  IN     Minimal allowed damping factor (see note 5)
C     23   SIGMA  IN     Broyden-approximation decision parameter
C                        Required choice: SIGMA.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs Broyden steps.
C                        Rank1 updates are inhibited, if 
C                        SIGMA.GT.1/FCMIN is set. (see note 5)
C     24                 Reserved
C     25   COND   IN     Maximum permitted subcondition for rank-
C                        decision by linear solver
C                        (Default: 1/epmach, epmach: relative
C                         machine precision) 
C     26   AJDEL  IN     Jacobian approximation without feedback:
C                        Relative pertubation for components
C                        (Default: sqrt(epmach*10), epmach: relative
C                         machine precision) 
C     27   AJMIN  IN     Jacobian approximation without feedback:
C                        Threshold value (Default: 0.0d0)
C                          The absolute pertubation for component k is
C                          computed by 
C                          DELX := AJDEL*max(abs(Xk),AJMIN)
C     28  ETADIF  IN     Jacobian approximation with feedback:
C                        Target value for relative pertubation ETA of X
C                        (Default: 1.0d-6)
C     29  ETAINI  IN     Jacobian approximation with feedback:
C                        Initial value for denominator differences
C                        (Default: 1.0d-6)
C     30                 Reserved
C     31  PREC    OUT    An estimate for the achieved relative accuracy.
C                        This number is only available, if IERR=0 or 
C                        IERR=1 and an estimate for the incompatibility
C                        factor kappa (SKAP=RWK(32)) can be made. If no
C                        meaningful information is at hand, a -1.0 is
C                        stored.
C     32  SKAP    OUT    An estimate of the incompatibility factor of
C                        the given nonlinear least squares problem.
C                        This number is only available, if IERR=0 or 
C                        IERR=1 and a certain further condition holds.
C                        If no meaningful information is at hand, a -1.0
C                        is stored.
C                        A small value of SKAP indicates that the given
C                        measurement data can be well fitted by the
C                        model, whereas a value SKAP.GE.0.5 may give a 
C                        hint to an insufficient modelling of the
C                        experiment from which the measurements
C                        originate.
C     33..49             Reserved
C     50  SIGMA2  OUT    Holds the estimated variance of the residual
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C     51..50+N
C         XL(N)          Holds the left bounds of the final parameters
C                        confidence intervals; 
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C     50+N+1..50+2*N
C         XR(N)          Holds the right bounds of the final parameters
C                        confidence intervals; 
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C     50+2*N+1..50+2*N+N*N         
C        VCV(N,N) OUT    Holds the columnwise stored correlation-matrix
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C
C     Note 5:
C       The default values of the internal parameters may be obtained
C       from the monitor output with at least IOPT field MPRMON set to 2
C       and by initializing the corresponding RWK-fields to zero. 
C
C*   Error messages:
C    ===============
C
C      1    Termination at stationary point (rank deficient Jacobian
C           and termination criterion fulfilled)
C      2    Termination after NITMAX iterations ( as indicated by
C           input parameter NITMAX=IWK(31) )
C      3    Termination, since damping factor became to small and
C           Jacobian rank was already reduced to zero
C     10    Integer or real workspace too small
C     20    Bad or inconsistent input to one or more of the 
C           dimensional parameters N,M and MFIT
C     21    Nonpositive value for RTOL supplied
C     22    Negative scaling value via vector XSCAL supplied
C     30    One or more fields specified in IOPT are invalid
C           (for more information, see error-printout)
C     80    Error signalled by linear solver routine NCFACT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used with standard routine NCFACT)
C     81    Error signalled by linear solver routine NCFIT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used with standard routine NCFIT)
C     82    Error signalled by user routine FCN (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     83    Error signalled by user routine JAC (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     180,182,183
C           see error codes 80,82,83, but the failure of NCFACT, FCN
C           or JAC occured when preparing the call of the statistics
C           subroutine STACON for the final iterate of NLSCON.
C
C     Note 6 : in case of failure:
C        -    use non-standard options
C        -    use another initial guess
C        -    or reformulate model
C        -    or turn to general optimization routine
C
C*    Machine dependent constants used:
C     =================================
C
C     DOUBLE PRECISION EPMACH  in  NLSCON, NCPCHK, NCINT
C     DOUBLE PRECISION GREAT   in  NLSCON, NCPCHK
C     DOUBLE PRECISION SMALL   in  NLSCON, NCPCHK, NCINT, NCSCAL
C
C*    Subroutines called: NCPCHK, NCINT, D1MACH
C
C     ------------------------------------------------------------
C*    End Prologue
C
C*    Summary of changes:
C     ===================
C      
C     2.2.1   91, June  3   Time monitor included
C     2.2.2   91, June  3   Bounded damping strategy implemented
C     2.2.3   91, July 26   AJDEL, AJMIN as RWK-options for JACGEN.EQ.2,
C                           ETADIF, ETAINI as RWK-opts. for JACGEN.EQ.3
C                           FCN-count changed for anal. Jacobian
C             91, Sept.     DECCON with new fail exit, for the case that
C                           the square root of a negative number will
C                           appear during pseudo inverse computation.
C                           (Occured, although theoretical impossible!)
C     2.2.6  91, Sept.  17  Damping factor reduction by FCN-fail imple-
C                           mented
C                Sept.  24  Meaning of option ITERM modified.
C     2.3    91, Dec.   20  New Release for CodeLib
C            92, June    2  Level of accepted simplified correction
C                           stored to RWK(IRWKI+4)
C     2.3.1  92, Oct.   13  Corrected errors in subr. NCJAC, NCJCF:
C                           dimensions of FX(M), FU(M),
C                           description of FCN in these subroutines
C     2.3.2  93, Nov.   24  Optional call of statistical analysis
C                           routine STACON to get statistics about
C                           the quality of the parameter estimate
C                           which has been computed by NLSCON.
C            00, July   12  RTOL output-value bug fixed
C 
C     ------------------------------------------------------------
C
C     PARAMETER (IRWKI=xx, LRWKI=yy)  
C     IRWKI: Start position of internally used RWK part
C     LRWKI: Length of internally used RWK part
C     (current values see parameter statement below)
C
C     INTEGER L4,L5,L51,L6,L61,L62,L7,L71,L8,L9,L10,L11,L12,L121,
C             L13,L14,L20
C     Starting positions in RWK of formal array parameters of internal
C     routine N1INT (dynamically determined in driver routine NLEQ1,
C     dependent on N and options setting)
C
C     Further RWK positions (only internally used)
C
C     Position  Name     Meaning
C
C     IRWKI              Reserved
C     IRWKI+1   FCA      Previous damping factor
C     IRWKI+(2..3)       Reserved
C     IRWKI+4   SUMXS    natural level of accepted simplified correction
C     IRWKI+(5..LRWKI-1) Free
C
C     Internal arrays stored in RWK (see routine NCINT for descriptions)
C
C     Position  Array         Type   Remarks
C
C     L4        AA(M,N)       Perm
C     L5        DX(N)         Perm  
C     L51       DXQ(N)        Perm 
C     L6        XA(N)         Perm
C     L61       FMODEL(M)     Perm
C     L62       F(M)          Perm
C     L7        FA(M)         Perm
C     L71       FW(M)         Perm
C     L72       ETA(N)        Perm   Only used for JACGEN=IOPT(3)=3
C     L8        A(M,N)        Temp
C     L9        QA(N,N)       Temp
C     L10       XW(N)         Temp
C     L11       DXQA(N)       Temp
C     L111      QU(M)         Temp
C     L112      RQ(M)         Temp
C     L113      RQKEEP(M)     Temp
C     L114      DELXQ(N)      Temp
C     L12       T1(N)         Temp
C     L121      T2(MAX(M,N))  Temp
C     L13       T3(M)         Temp   
C     L14                     Temp   Start position of array workspace 
C                                    needed for linear solver  
C     
C
      EXTERNAL D1MACH,NCINT
      INTRINSIC DBLE,MAX0,MIN0
      INTEGER IRWKI, LRWKI
      PARAMETER (IRWKI=51, LRWKI=10)  
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,IRANK,IBASFW,NITMAX,LUERR,LUMON,LUSOL,MAXMN,MCON,MINMN,
     $MPRERR,MPRMON,MPRSOL,M1,M2,NRWKFR,NRFRIN,NRW,NIWKFR,NIFRIN,NIW,
     $NONLIN,IRWKII,JACGEN,LUTIM,MPRTIM
      INTEGER L4,L5,L51,L6,L61,L62,L7,L71,L8,L9,L10,L11,L111,L12,L121,
     $L13,L14,L20,L112,L113,L114,L30,L31,L32,L72
      DOUBLE PRECISION COND,D1MACH,FC,FCMIN,PERCI,PERCR,TOLMIN
      DOUBLE PRECISION EPMACH,GREAT,SMALL
      LOGICAL QINIMO,QRANK1,QFCSTR,QSUCC,QBDAMP
      CHARACTER CHGDAT*20, PRODCT*8
C     Which version ?
      LOGICAL QVCHK
      INTEGER IVER
      PARAMETER( IVER=22112321 )
C
C     Version: 2.3.2             Latest change:
C     -----------------------------------------
C
      DATA      CHGDAT      /'July 12, 2000       '/
      DATA      PRODCT      /'NLSCON  '/
C*    Begin
      EPMACH = D1MACH(3)
      GREAT  = DSQRT(D1MACH(2)/TEN)
      SMALL  = D1MACH(6)
      TOLMIN = EPMACH*TEN*DBLE(N)
      IERR = 0
      QVCHK = IWK(12).LT.0
      IWK(12) = IVER
      IF (QVCHK) RETURN
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .EQ. 0) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C        Print iteration monitor?
      MPRMON = IOPT(13)
      LUMON = IOPT(14)
      IF (LUMON .LE. 0 .OR. LUMON .GT. 99) THEN
        LUMON = 6
        IOPT(14)=LUMON
      ENDIF
C        Print intermediate solutions?
      MPRSOL = IOPT(15)
      LUSOL = IOPT(16)
      IF (LUSOL .EQ. 0) THEN
        LUSOL = 6
        IOPT(16)=LUSOL
      ENDIF
C        Print time summary statistics?
      MPRTIM = IOPT(19)
      LUTIM = IOPT(20)
      IF (LUTIM .EQ. 0) THEN
        LUTIM = 6
        IOPT(20)=LUTIM
      ENDIF
      QSUCC = IOPT(1).EQ.1
      QINIMO = MPRMON.GE.1.AND..NOT.QSUCC
C     Print NLSCON heading lines
      IF(QINIMO)THEN
10000   FORMAT('    N L S C O N  *****  V e r s i o n  ',
     $         '2 . 3 . 2 ***',//,1X,'Gauss-Newton-Method for the ',
     $         'solution of nonlinear least squares problems',//)
        WRITE(LUMON,10000)
      ENDIF
C     Check input parameters and options
      CALL NCPCHK(N,M,MFIT,X,XSCAL,FI,FSCAL,RTOL,IOPT,
     $            IERR,LIWK,IWK,LRWK,RWK)
C     Exit, if any parameter error was detected till here
      IF (IERR.NE.0) RETURN 
      MCON = M-MFIT
      M1=M
      M2=M
      MAXMN=MAX0(M,N)
      JACGEN=IOPT(3)
      IF (JACGEN.EQ.0) JACGEN=2
      IOPT(3)=JACGEN
C     WorkSpace: RWK
      L30 = IRWKI
      L31=L30+N
      L32=L31+N
      IRWKII = L32+N*N
      L4=IRWKII+LRWKI
      L5=L4+M2*N
      L51=L5+N
      L6=L51+N
      L61=L6+N
      L62=L61+M
      L7=L62+M
      L71=L7+M
      L72=L71+M
      IF (JACGEN.NE.3) THEN
        L8=L72
      ELSE
        L8=L72+N
      ENDIF
      NRWKFR = L8
      L9=L8+M1*N
      L10=L9+N*N
      L11=L10+N
      L111=L11+N
      L112=L111+M
      L113=L112+M
      L114=L113+M
      L12=L114+N
      L121=L12+N
      L13=L121+MAXMN
      L14=L13+M
      NRW=L14-1
C     End WorkSpace at NRW
C     WorkSpace: IWK
      L20=51
      NIWKFR = L20
      NIW=L20-1
C     End WorkSpace at NIW
      IWK(16) = NIW+1
      IWK(17) = NRW+1
      NIFRIN = NIW+1
      NRFRIN = NRW+1
C
      IF(NRW.GT.LRWK.OR.NIW.GT.LIWK)THEN
        IERR=10
      ELSE
        IF(QINIMO)THEN
          PERCR = DBLE(NRW)/DBLE(LRWK)*100.0D0
          PERCI = DBLE(NIW)/DBLE(LIWK)*100.0D0
C         Print statistics concerning workspace usage
10050     FORMAT(' Real    Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//,
     $    ' Integer Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//)
          WRITE(LUMON,10050)LRWK,NRW,PERCR,LIWK,NIW,PERCI
        ENDIF
        IF(QINIMO)THEN
10051     FORMAT(/,' Number of parameters to be estimated (N) : ',I4,/,
     $    ' Number of data to fitted, e.g. observations (MFIT) : ',
     $    I4,/,' Number of equality constraints (MCON) : ',I4,//,
     $    ' Prescribed relative precision (RTOL) : ',D10.2,/)
          WRITE(LUMON,10051)N,MFIT,MCON,RTOL
10052     FORMAT(' The Jacobian is supplied by ',A)
          IF (JACGEN.EQ.1) THEN
            WRITE(LUMON,10052) 'a user subroutine'
          ELSE IF (JACGEN.EQ.2) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (without feedback strategy)'
          ELSE IF (JACGEN.EQ.3) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (feedback strategy included)'
          ENDIF
10057     FORMAT(' Automatic row scaling of the Jacobian is ',A,/)
          IF (IOPT(35).EQ.1) THEN
            WRITE(LUMON,10057) 'inhibited'
          ELSE
            WRITE(LUMON,10057) 'allowed'
          ENDIF
        ENDIF
        QRANK1=IOPT(32).EQ.1
        NONLIN=IOPT(31)
        IF (IOPT(38).EQ.1) QBDAMP = .TRUE.
        IF (IOPT(38).EQ.0 .OR. IOPT(38).EQ.2) QBDAMP = .FALSE.
        IF (QBDAMP) THEN
          IF (RWK(20).LT.ONE) RWK(20) = TEN
        ENDIF
        IF (QRANK1.AND.M.GT.N) THEN
10063     FORMAT(/,' Warning: Broyden-steps set to inhibited for ',
     $           'the overdetermined system')
          IF (MPRERR.GE.2) WRITE(LUERR,10063)
          QRANK1 = .FALSE.
          IOPT(32) = 0
        ENDIF
10064   FORMAT(' Rank-1 updates are ',A)
        IF (QINIMO) THEN
          IF (QRANK1) THEN
            WRITE(LUMON,10064) 'allowed'
          ELSE
            WRITE(LUMON,10064) 'inhibited'
          ENDIF
10065     FORMAT(' Problem is specified as being ',A)
          IF (NONLIN.EQ.1) THEN
            WRITE(LUMON,10065) 'linear'
          ELSE IF (NONLIN.EQ.2) THEN
            WRITE(LUMON,10065) 'mildly nonlinear'
          ELSE
            WRITE(LUMON,10065) 'highly nonlinear'
          ENDIF
10066     FORMAT(' Bounded damping strategy is ',A,:,/, 
     $           ' Bounding factor is ',D10.3)
          IF (QBDAMP) THEN
            WRITE(LUMON,10066) 'active', RWK(20)
          ELSE
            WRITE(LUMON,10066) 'off'
          ENDIF
        ENDIF
C       Maximum permitted number of iteration steps
        NITMAX=IWK(31)
        IF (NITMAX.LE.0) NITMAX=50
        IWK(31)=NITMAX
10068   FORMAT(' Maximum permitted number of iteration steps : ',
     $         I6)
        IF (QINIMO) WRITE(LUMON,10068) NITMAX
C       Initial damping factor for highly nonlinear problems
        QFCSTR=RWK(21).GT.ZERO
        IF (.NOT.QFCSTR) RWK(21)=1.0D-2
C       Minimal permitted damping factor
        IF (RWK(22).LE.ZERO) RWK(22)=1.0D-2
        FCMIN=RWK(22)
C       Broyden-update decision parameter SIGMA
        IF (RWK(23).LT.ONE) RWK(23)=2.0D0
        IF (.NOT.QRANK1) RWK(23)=10.0D0/FCMIN
C         Starting value of damping factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.2.AND..NOT.QFCSTR)THEN
C         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
C         for highly nonlinear problems
          FC = RWK(21)
        ENDIF
        RWK(21)=FC
C       Initial rank
        IRANK = IWK(32)
        MINMN=MIN0(M,N)
        IF (IRANK.LE.0.OR.IRANK.GT.MINMN) IWK(32) = MINMN
C       Maximum permitted sub condition number of matrix A
        COND = RWK(25)
        IF (COND.LT.ONE) COND = ONE/EPMACH
        RWK(25) = COND
        IF (MPRMON.GE.2.AND..NOT.QSUCC) THEN
10069     FORMAT(//,' Internal parameters:',//,
     $      ' Starting value for damping factor FCSTART = ',D9.2,/,
     $      ' Minimum allowed damping factor FCMIN = ',D9.2,/,
     $      ' Rank-1 updates decision parameter SIGMA = ',D9.2,/,
     $      ' Initial Jacobian pseudo-rank IRANK =',I6,/,
     $      ' Maximum permitted subcondition COND = ',D9.2)
          WRITE(LUMON,10069) RWK(21),FCMIN,RWK(23),IWK(32),COND
        ENDIF
        IF (.NOT.QSUCC) THEN
           IBASFW = L71-1+MCON
           DO 1 I=1,MFIT
             IF (FSCAL(I).GE.SMALL.AND.FSCAL(I).LE.GREAT) THEN
               RWK(IBASFW+I) = ONE/FSCAL(I)
             ELSE
               RWK(IBASFW+I) = ONE
               IF (FSCAL(I).NE.ZERO.AND.MPRERR.GE.2) THEN
                 WRITE(LUERR,10070) I,FSCAL(I)
10070            FORMAT(/,' Warning: Bad scaling value FSCAL(',I5,
     $                   ') = ',D10.2,' replaced by 1.0')
               ENDIF
             ENDIF
1          CONTINUE
        ENDIF
C       Store lengths of currently required workspaces
        IWK(18) = NIFRIN-1
        IWK(19) = NRFRIN-1
C
C       Initialize and start time measurements monitor
C
        IF ( IOPT(1).EQ.0 .AND. MPRTIM.NE.0 ) THEN
          CALL MONINI (' NLSCON',LUTIM)
          CALL MONDEF (0,'NLSCON')
          CALL MONDEF (1,'FCN')
          CALL MONDEF (2,'Jacobi')
          CALL MONDEF (3,'Lin-Fact')
          CALL MONDEF (4,'Lin-Sol')
          CALL MONDEF (5,'Output')
          CALL MONSRT ()
        ENDIF
C
C
        IERR=-1
C       If IERR is unmodified on exit, successive steps are required
C       to complete the Gauss-Newton iteration
        CALL NCINT(N,M,MFIT,FCN,JAC,X,XSCAL,FI,RTOL,NITMAX,NONLIN,
     $  IWK(32),IOPT,IERR,LRWK,RWK,NRFRIN,LIWK,IWK,NIFRIN,
     $  M1,M2,MAXMN,RWK(L30),RWK(L31),RWK(L32),
     $  RWK(L4),RWK(L8),RWK(L9),RWK(L5),RWK(L51),RWK(L6),RWK(L61),
     $  RWK(L62),RWK(L7),RWK(L72),RWK(L71),RWK(L10),RWK(L11),RWK(L111),
     $  RWK(L112),RWK(L113),RWK(L114),RWK(L12),RWK(L121),RWK(L13),
     $  RWK(21),RWK(22),RWK(23),RWK(IRWKII+1),COND,
     $  RWK(17),RWK(18),RWK(IRWKII+4),RWK(19),
     $  TOLMIN,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,IWK(1),
     $  IWK(3),IWK(4),IWK(5),IWK(8),IWK(9),IWK(33),IWK(34),QBDAMP)
C
      IF (MPRTIM.NE.0.AND.IERR.NE.-1.AND.IERR.NE.10) CALL MONEND
C
C       Free workspaces, so far not used between steps
        IWK(16) = NIWKFR
        IWK(17) = NRWKFR
      ENDIF
C     Print statistics
      IF (MPRMON.GE.1.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
10080   FORMAT(/, '   ******  Statistics * ', A8, ' *******', /,
     $            '   ***  Gauss-Newton iter.: ', I7,'  ***', /,
     $            '   ***  Corrector steps   : ', I7,'  ***', /,
     $            '   ***  Rejected rk-1 st. : ', I7,'  ***', /,
     $            '   ***  Jacobian eval.    : ', I7,'  ***', /,
     $            '   ***  Function eval.    : ', I7,'  ***', /,
     $            '   ***  ...  for Jacobian : ', I7,'  ***', /,
     $            '   *************************************', /)
        WRITE (LUMON,10080) PRODCT,IWK(1),IWK(3),IWK(9),IWK(5),
     $  IWK(4),IWK(8)
      ENDIF
C     Print workspace requirements, if insufficient
      IF (IERR.EQ.10) THEN
10090   FORMAT(///,20('*'),'Workspace Error',20('*'))
        IF (MPRERR.GE.1) WRITE(LUERR,10090)
        IF(NRW.GT.LRWK)THEN
10091     FORMAT(/,' Real Workspace dimensioned as',1X,I9,
     $    1X,'must be enlarged at least up to ',
     $    I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10091)LRWK,NRFRIN-1
        ENDIF
        IF(NIW.GT.LIWK)THEN
10092     FORMAT(/,' Integer Workspace dimensioned as ',
     $    I9,' must be enlarged at least up ',
     $    'to ',I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10092)LIWK,NIFRIN-1
        ENDIF
      ENDIF
C     End of subroutine NLSCON
      RETURN
      END
C
      SUBROUTINE NCPCHK(N,M,MFIT,X,XSCAL,FI,FSCAL,RTOL,IOPT,
     $IERR,LIWK,IWK,LRWK,RWK)
C*    Begin Prologue NCPCHK
      INTEGER N,M,MFIT
      DOUBLE PRECISION X(N),XSCAL(N),FI(MFIT),FSCAL(MFIT)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N C P C H K : Checking of input parameters and options
C                   for NLSCON.
C
C*    Parameters:
C     ===========
C
C     See parameter description in driver routine.
C
C*    Subroutines called: D1MACH
C
C*    Machine dependent constants used:
C     =================================
C
C     EPMACH = relative machine precision
C     GREAT = squareroot of maxreal divided by 10
C     SMALL = squareroot of "smallest positive machine number
C             divided by relative machine precision"
      DOUBLE PRECISION EPMACH,GREAT,SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
C
      EXTERNAL D1MACH
      INTRINSIC DBLE
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      INTEGER NUMOPT
      PARAMETER (NUMOPT=50)
      INTEGER IOPTL(NUMOPT),IOPTU(NUMOPT)
      DOUBLE PRECISION D1MACH,TOLMIN,TOLMAX,DEFSCL
      INTEGER I,MPRERR,LUERR,NONLIN
C
      DATA IOPTL /0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,
     $            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,0,0,0,
     $            -9999999,-9999999,-9999999,-9999999,-9999999/
      DATA IOPTU /1,1,3,0,0,0,0,0,1,0,3,99,6,99,3,99,0,0,1,99,
     $            1,1,0,0,0,0,0,0,0,0,4,1,0,0,1,
     $            2,0,2,0,0,0,0,0,0,0,
     $            9999999,9999999,9999999,9999999,9999999/
C
      EPMACH = D1MACH(3)
      GREAT  = DSQRT(D1MACH(2)/TEN)
      SMALL  = D1MACH(6)
      IERR = 0
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .LE. 0 .OR. LUERR .GT. 99) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C
C     Checking dimensional parameters N,M and MFIT
      IF ( N.LE.0 .OR. M.LE.0 .OR. MFIT.LT.0 .OR. MFIT.GT.M ) THEN
        IF (MPRERR.GE.1)  WRITE(LUERR,10011) N,M,MFIT
10011   FORMAT(/,' Error: Bad or inconsistent input to dimensional ',
     $         'parameters supplied',/,
     $         8X,'choose N and M positive, and MFIT.LE.M, MFIT ',
     $         'nonnegative',/,
     $         8X,'your input is: N = ',I5,' M = ',I5,' MFIT = ',I5)
        IERR = 20
      ENDIF
C
C     Problem type specification by user
      NONLIN=IOPT(31)
      IF (NONLIN.LE.0) NONLIN=3
      IOPT(31)=NONLIN
C
C     Checking and conditional adaption of the user-prescribed RTOL
      IF (RTOL.LE.ZERO) THEN
        IF (MPRERR.GE.1) 
     $      WRITE(LUERR,'(A)') '0Error: Nonpositive RTOL supplied'
        IERR = 21
      ELSE
        TOLMIN = EPMACH*TEN*DBLE(N)
        IF(RTOL.LT.TOLMIN) THEN
          RTOL = TOLMIN
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'increased ','smallest',RTOL
        ENDIF
        TOLMAX = 1.0D-1
        IF(RTOL.GT.TOLMAX) THEN
          RTOL = TOLMAX
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'decreased ','largest',RTOL
        ENDIF
10012   FORMAT(/,' Warning: User prescribed RTOL ',A,'to ',
     $         'reasonable ',A,' value RTOL = ',D11.2)
      ENDIF
C     
C     Test user prescribed accuracy and scaling on proper values
      IF (N.LE.0) RETURN 
      IF (NONLIN.EQ.3) THEN
        DEFSCL = RTOL
      ELSE
        DEFSCL = ONE
      ENDIF
      DO 10 I=1,N
        IF (XSCAL(I).LT.ZERO) THEN
          IF (MPRERR.GE.1) THEN 
            WRITE(LUERR,10013) I
10013       FORMAT(/,' Error: Negative value in XSCAL(',I5,') supplied')
          ENDIF
          IERR = 22
        ENDIF
        IF (XSCAL(I).EQ.ZERO) XSCAL(I) = DEFSCL
        IF ( XSCAL(I).GT.ZERO .AND. XSCAL(I).LT.SMALL ) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10014) I,XSCAL(I),SMALL
10014       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too small, ',
     $             'increased to',D9.2)
          ENDIF
          XSCAL(I) = SMALL
        ENDIF
        IF (XSCAL(I).GT.GREAT) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10015) I,XSCAL(I),GREAT
10015       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too big, ',
     $             'decreased to',D9.2)
          ENDIF
          XSCAL(I) = GREAT
        ENDIF
10    CONTINUE
C     Checks options
      DO 20 I=1,30
        IF (IOPT(I).LT.IOPTL(I) .OR. IOPT(I).GT.IOPTU(I)) THEN
          IERR=30
          IF (MPRERR.GE.1) THEN
            WRITE(LUERR,20001) I,IOPT(I),IOPTL(I),IOPTU(I)
20001       FORMAT(' Invalid option specified: IOPT(',I2,')=',I12,';',
     $             /,3X,'range of permitted values is ',I8,' to ',I8)
          ENDIF
        ENDIF
20    CONTINUE
C     End of subroutine NCPCHK
      RETURN
      END
C
      SUBROUTINE NCINT(N,M,MFIT,FCN,JAC,X,XSCAL,FI,RTOL,NITMAX,NONLIN,
     $IRANK,IOPT,IERR,LRWK,RWK,NRWKFR,LIWK,IWK,NIWKFR,M1,M2,MAXMN,
     $XL,XR,VCV,
     $AA,A,QA,DX,DXQ,XA,FMODEL,F,FA,ETA,FW,XW,DXQA,QU,RQ,RQKEEP,DELXQ,
     $T1,T2,T3,FC,FCMIN,SIGMA,FCA,COND,CONV,SUMX,SUMXS,DLEVF,TOLMIN,
     $MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,NCORR,NFCN,NJAC,
     $NFCNJ,NREJR1,NEW,IFCCNT,QBDAMP)
C*    Begin Prologue NCINT
      INTEGER N,M,MFIT
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N),FI(MFIT)
      DOUBLE PRECISION RTOL
      INTEGER NITMAX,NONLIN,IRANK
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER NRWKFR,LIWK
      INTEGER IWK(LIWK)
      INTEGER NIWKFR,M1,M2,MAXMN
      DOUBLE PRECISION VCV(N,N),XL(N),XR(N),AA(M2,N),A(M1,N),QA(N,N)
      DOUBLE PRECISION DX(N),DXQ(N),XA(N),FMODEL(M),F(M),FA(M),ETA(N)
      DOUBLE PRECISION XW(N),FW(M),DXQA(N),QU(M),RQ(M),RQKEEP(M),
     $                 DELXQ(N),T1(N),T2(MAXMN),T3(M)
      DOUBLE PRECISION FC,FCMIN,SIGMA,FCA,COND,CONV,SUMX,SUMXS,DLEVF,
     $                 TOLMIN
      INTEGER MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,
     $NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW,IFCCNT
      LOGICAL QBDAMP
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N C I N T : Core routine for NLSCON .
C     Damped Gauss-Newton-algorithm with rank-strategy for highly
C     nonlinear least squares estimation problems especially
C     designed for numerically sensitive problems.
C
C*    Parameters:
C     ===========
C
C       N,M,MFIT,FCN,JAC,X,XSCAL,FI,RTOL   
C                         See parameter description in driver routine
C
C       NITMAX     Int    Maximum number of allowed iterations
C       NONLIN     Int    Problem type specification
C                         (see IOPT-field NONLIN)
C       IRANK      Int    Initially proposed (in) and final (out) rank
C                         of Jacobian
C       IOPT       Int    See parameter description in driver routine
C       IERR       Int    See parameter description in driver routine
C       LRWK       Int    Length of real workspace
C       RWK(LRWK)  Dble   Real workspace array
C       NRWKFR     Int    First free position of RWK on exit 
C       LIWK       Int    Length of integer workspace
C       IWK(LIWK)  Int    Integer workspace array
C       NIWKFR     Int    First free position of IWK on exit 
C       M1         Int    Leading dimension of Jacobian array A
C                         for full case Jacobian: M
C                         (other matrix types are not yet implemented)
C       M2         Int    Leading dimension of Jacobian array AA
C                         for full case Jacobian: M
C       MAXMN      Int    Max(M,N) - length of some temporary
C                         workspace
C       VCV(N,N)   Dble   Correlation matrix, if IOPT(21)=1 or
C                         IOPT(22)=1 (out)
C       XL(N)     Dble    Confidence interval left sides,
C                         if IOPT(21)=1 or IOPT(22)=1 (out)
C       XR(N)     Dble    Confidence interval right sides,
C                         if IOPT(21)=1 or IOPT(22)=1 (out)
C       AA(M2,N)   Dble   Holds the originally computed Jacobian
C       A(M1,N)    Dble   Holds the Jacobian matrix (decomposed form
C                         after call of linear decomposition
C                         routine)
C       QA(N,N)    Dble   Holds the pseudo inverse in case of rank-
C                         deficiency
C       DX(N)      Dble   Current Gauss-Newton correction
C       DXQ(N)     Dble   Simplified Gauss-Newton correction J(k)*X(k+1)
C       XA(N)      Dble   Previous Gauss-Newton iterate
C       F(M)       Dble   Function (FCN) value minus measured values
C                         vector FI of current iterate
C       FA(M)      Dble   Holds the previous values of vector F(M)
C       ETA(N)     Dble   Jacobian approximation: updated scaled
C                         denominators
C       FW(M)      Dble   Scaling factors for rows of the system
C       XW(N)      Dble   Scaling factors for iteration vector
C       DXQA(N)    Dble   Previous simplified Gauss-Newton 
C                         correction J(k-1)*X(k)
C       QU(M)      Dble   Savespace for right hand side belonging
C                         to upper triangular linear system
C       RQ(M)      Dble   Gets the residuum of the linear problems
C                         solution in each iteration step :
C                         JAC(X(i-1))*DXQ(i)+F(i) for the i-th iterate
C       RQKEEP(M)  Dble   Keeps a copy of RQ(M) for restoring it in
C                         case of Jacobian recomputation (M.LE.N ,
C                         rejected rank-1 step, and Jacobian of
C                         previous iterate was rank-deficient).
C       DELXQ(N)   Dble   Gets the vector of parameters best fitting
C                         the linear problems residuum of the previous
C                         iterate
C       T1(N)      Dble   Workspace for linear solvers and internal
C                         subroutines
C       T2(MAXMN)  Dble   Workspace array for internal subroutines
C       T3(M)      Dble   Workspace array for internal subroutines
C       FC         Dble   Current Gauss-Newton iteration damping factor.
C       FCMIN      Dble   Minimum permitted damping factor. If
C                         FC becomes smaller than this value, one
C                         of the following may occur:
C                         a.    Recomputation of the Jacobian
C                               matrix by means of difference
C                               approximation (instead of Rank1
C                               update), if Rank1 - update
C                               previously was used
C                         b.    Rank reduction of Jacobi
C                               matrix ,  if difference
C                               approximation was used previously
C                               and Rank(A).NE.0
C                         c.    Fail exit otherwise
C       SIGMA      Dble   Decision parameter for rank1-updates.
C       FCA        Dble   Previous Gauss-Newton iteration damping factor.
C       COND       Dble   Maximum permitted subcondition for rank-
C                         decision by linear solver.
C       CONV       Dble   Scaled maximum norm of the Gauss-Newton-
C                         correction. Passed to RWK-field on output.
C       SUMX       Dble   Square of the natural level (see equal-
C                         named IOPT-output field)
C       SUMXS       Dble   Square of the "simplified" natural level
C                          (see equal-named RWK-internal field)
C       DLEVF      Dble   Square of the standard level (see equal-
C                         named IOPT-output field)
C       MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL :
C                         See description of equal named IOPT-fields
C                         in the driver subroutine
C       NITER,NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW :
C                         See description of equal named IWK-fields
C                         in the driver subroutine
C       QBDAMP     Logic  Flag, that indicates, whether bounded damping
C                         strategy is active:
C                         .true.  = bounded damping strategy is active
C                         .false. = normal damping strategy is active
C
C
C*    Internal double variables
C     =========================
C
C       AJDEL    See RWK(26) (num. diff. without feedback)
C       AJMIN    See RWK(27) (num. diff. without feedback)
C       COND1    Gets the subcondition of the linear system
C                as estimated by the linear solver (NCFACT)
C       CONVA    Holds the previous value of CONV .
C       DELQ     Gets the projection defect in case of rank-
C                deficiency.
C       DMUE     Temporary value used during computation of damping 
C                factors predictor.
C       EPDIFF   sqrt(10*epmach) (num. diff. with feedback)
C       ETADIF   See description of RWK(28) (num. diff. with feedback)
C       ETAINI   Initial value for all ETA-components (num. diff. fb.)
C       ETAMAX   Maximum allowed pertubation (num. diff. with feedback)
C       ETAMIN   Minimum allowed pertubation (num. diff. with feedback)
C       FCDNM    Used to compute the denominator of the damping 
C                factor FC during computation of it's predictor,
C                corrector and aposteriori estimate (in the case of
C                performing a Rank1 update) .
C       FCK2     Aposteriori estimate of FC.
C       FCMIN2   FCMIN**2 . Used for FC-predictor computation.
C       FCMINH   DSQRT(FCMIN).
C                Used in rank decision logical expression.
C       FCNUMP   Gets the numerator of the predictor formula for FC.
C       FCNMP2   Temporary used for predictor numerator computation.
C       FCNUMK   Gets the numerator of the corrector computation 
C                of FC .
C       SENS1    Gets the sensitivity of the Jacobian as
C                estimated by the linear solver NCFACT.
C       SKAP     In case of termination at stationary point:
C                incompatibility factor
C       SUMXA    Natural level of the previous iterate.
C       TH       Temporary variable used during corrector- and 
C                aposteriori computations of FC.
C
C*    Internal integer variables
C     ==========================
C
C     IFAIL      Gets the return value from subroutines called from
C                NCINT (NCFACT, NCFIT, FCN, JAC)
C     ISCAL      Holds the scaling option from the IOPT-field ISCAL      
C     MODE       Matrix storage mode (see IOPT-field MODE) 
C     NRED       Count of successive corrector steps
C     NILUSE     Gets the amount of IWK used by the linear solver
C     NRLUSE     Gets the amount of RWK used by the linear solver
C     NIWLA      Index of first element of IWK provided to the
C                linear solver
C     NRWLA      Index of first element of RWK provided to the
C                linear solver
C     LIWL       Holds the maximum amount of integer workspace
C                available to the linear solver
C     LRWL       Holds the maximum amount of real workspace
C                available to the linear solver
C
C*    Internal logical variables
C     ==========================
C
C     QGENJ      Jacobian updating technique flag:
C                =.TRUE.  : Call of analytical subroutine JAC or
C                           numerical differentiation
C                =.FALSE. : rank1- (Broyden-) update
C     QINISC     Iterate initial-scaling flag:
C                =.TRUE.  : at first call of NCSCAL
C                =.FALSE. : at successive calls of NCSCAL
C     QSUCC      See description of IOPT-field QSUCC.
C     QJCRFR     Jacobian refresh flag:
C                set to .TRUE. if damping factor gets too small
C                and Jacobian was computed by rank1-update. 
C                Indicates, that the Jacobian needs to be recomputed
C                by subroutine JAC or numerical differentation.
C     QLINIT     Initialization state of linear solver workspace:
C                =.FALSE. : Not yet initialized
C                =.TRUE.  : Initialized - NCFACT has been called at
C                           least one time.
C     QREPET     Operation mode flag for linear solver:
C                =.FALSE. : Normal operation (full rank matrix)
C                =.TRUE.  : Special operation in rank deficient case:
C                           Compute rank-deficient pseudo-inverse,
C                           partial recomputation when solving the
C                           linear system again prescribing a lower
C                           rank as before.
C     QNEXT      Set to .TRUE. to indicate success of the current
C                Gauss-Newton-step, i.e. : sucessfull monotonicity-test.
C     
C     QREDU      Set to .TRUE. to indicate that rank-reduction (or
C                refreshment of the Jacobian) is needed - if the
C                computed damping factor gets too small.
C     QSCALE     Holds the value of .NOT.QNSCAL. See description
C                of IOPT-field QNSCAL.
C
C*    Subroutines called:
C     ===================
C
C       NCFACT, NCFIT,  NCJAC,  NCJCF, NCLVLS, NCSCRF, NCPRJN,
C       NCRNK1, NCSOUT, NCPRV1, NCPRV2, NCSCAL,
C       MONON,  MONOFF
C
C*    Functions called:
C     =================
C
C       D1MACH, WNORM
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION EPMACH,SMALL
C 
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL NCFACT, NCFIT,  NCJAC,  NCJCF, NCLVLS, NCSCRF, NCPRJN,
     $         NCRNK1, NCSOUT, NCPRV1, NCPRV2, NCSCAL,
     $         MONON,  MONOFF, D1MACH, WNORM
      INTRINSIC DSQRT,DMIN1,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=10.0D0)
      INTEGER I,IFAIL,ISCAL,ITERM,K,MFOUT,MINMN,MODE,MCON,MCONH,ML,MU,
     $NRED,NILUSE,NRLUSE,NIWLA,NRWLA,LIWL,LRWL,L1,L2,IRANKA,IRANKC,
     $JACGEN,MODEFI,MPRTIM,IDUMMY
      DOUBLE PRECISION AJDEL,AJMIN,COND1,CONDCO,CONVA,DLEVXA,DELQ,
     $DMUE,D1MACH,DXANRM,DXNRM,WNORM,EPDIFF,ETAMIN,ETAMAX,ETAINI,
     $ETADIF,FCDNM,FCH,FCBND,FCBH,FCK2,FCMIN2,FCMINH,FCNUMP,FCNMP2,
     $FCNUMK,FCREDU,PREC,SENS1,SKAP,SUMXA,SUMXK,S1,TH,RSMALL,DUMMY(1)
      LOGICAL QGENJ,QINISC,QSUCC,QJCRFR,QLINIT,QREPET,QNEXT,QREDU,
     $        QSCALE,QRANK1,QMIXIO
      EPMACH = D1MACH(3)
      SMALL  = D1MACH(6)
C*    Begin
C       ----------------------------------------------------------
C       1 Initialization
C       ----------------------------------------------------------
C       1.1 Control-flags and -integers
        QSUCC = IOPT(1).EQ.1
        QSCALE = .NOT. IOPT(35).EQ.1
        QRANK1=IOPT(32).EQ.1
        QMIXIO = LUMON.EQ.LUSOL .AND. MPRMON.NE.0 .AND. MPRSOL.NE.0
        ISCAL = IOPT(9)
        ITERM = IOPT(36)
        MODE = IOPT(2)
        JACGEN = IOPT(3)
        MPRTIM = IOPT(19)
C       ----------------------------------------------------------
C       1.2 Derivated dimensional parameters
        MCON = M-MFIT
C       ----------------------------------------------------------
C       1.3 Derivated internal parameters
        MINMN = MIN0(M,N)
        FCMIN2 = FCMIN*FCMIN
        FCMINH = DSQRT(FCMIN)
        RSMALL = DSQRT(TEN*RTOL)
C       ----------------------------------------------------------
C       1.4 Adaption of input parameters, if necessary
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
C       ----------------------------------------------------------
C       1.5 Initial preparations
        QJCRFR = .FALSE.
        QLINIT = .FALSE.
        QREPET = .FALSE.
        IFAIL = 0
        FCBND = ZERO
        IF (QBDAMP) FCBND = RWK(20)
C       ----------------------------------------------------------
C       1.5.1 Numerical differentation related initializations
        IF (JACGEN.EQ.2) THEN
          AJDEL = RWK(26)
          IF (AJDEL.LE.SMALL) AJDEL = DSQRT(EPMACH*TEN)
          AJMIN = RWK(27)
        ELSE IF (JACGEN.EQ.3) THEN
          ETADIF = RWK(28)
          IF (ETADIF .LE. SMALL) ETADIF = 1.0D-6
          ETAINI = RWK(29)
          IF (ETAINI .LE. SMALL) ETAINI = 1.0D-6
          EPDIFF = DSQRT(EPMACH*TEN)
          ETAMAX = DSQRT(EPDIFF)
          ETAMIN = EPDIFF*ETAMAX
        ENDIF
C       ----------------------------------------------------------
C       1.5.2 Miscellaneous preparations of first iteration step
        IF (.NOT.QSUCC) THEN
          NITER = 0
          NCORR = 0
          NREJR1 = 0
          NFCN = 0
          NJAC = 0
          NFCNJ = 0
          IFCCNT = 0
          IRANKA = IRANK
          QGENJ = .TRUE.
          QINISC = .TRUE.
          FCA = FC
          FCK2 = FC
          CONV = ZERO
          IF (JACGEN.EQ.3) THEN
            DO 1520 L1=1,N
              ETA(L1)=ETAINI
1520        CONTINUE
          ENDIF
          DO 1521 L1=1,N
            XA(L1)=X(L1)
1521      CONTINUE
C         ------------------------------------------------------
C         1.6 Print monitor header
          IF(MPRMON.GE.2 .AND. .NOT.QMIXIO)THEN
16003       FORMAT(///,2X,74('*'))
            WRITE(LUMON,16003)
16004       FORMAT(/,8X,'It',7X,'Normf ',14X,'Normx ',6X,
     $             'Damp.Fct.',3X,'New',6X,'Rank')
            WRITE(LUMON,16004)
          ENDIF
C         --------------------------------------------------------
C         1.7 Startup step
C         --------------------------------------------------------
C         1.7.1 Computation of the residual vector
          IF (MPRTIM.NE.0) CALL MONON(1)
          CALL FCN(N,M,MCON,X,FMODEL,IFAIL)
          IF (MPRTIM.NE.0) CALL MONOFF(1)
          NFCN = NFCN+1
C     Exit, if ...
          IF (IFAIL.NE.0) THEN
            IERR = 82
            GOTO 4299
          ENDIF
          DO 1711 L1=1,MFIT
            F(L1+MCON)=FMODEL(L1+MCON)-FI(L1)
1711      CONTINUE
          DO 1712 L1=1,MCON
            F(L1)=FMODEL(L1)
1712      CONTINUE
        ELSE
          QINISC = .FALSE.
        ENDIF
C
C       Main iteration loop
C       ===================
C
C       Repeat
2       CONTINUE
C         --------------------------------------------------------
C         2 Startup of iteration step
          IF (.NOT.QJCRFR) THEN
C           ------------------------------------------------------
C           2.0 Scaling of variables X(N)
            CALL NCSCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,LRWK,RWK)
            QINISC = .FALSE.
            IF(NITER.NE.0)THEN
C             Preliminary pseudo-rank
              IRANKA = IRANK
              DO 210 L1=1,N
                DXQA(L1)=DXQ(L1)
210           CONTINUE
C             ----------------------------------------------------
C             2.1 Computation of linear residuum RQ(M)
              IF (M.GT.N.OR.IRANKA.NE.M) THEN
                DO 211 L1=MCON+1,M
                  S1=ZERO
                  DO 2111 L2=1,N
                    S1=S1+AA(L1,L2)*DXQA(L2)
2111              CONTINUE
                  RQ(L1)=S1+F(L1)
211             CONTINUE
                DO 212 L1=1,MCON
                  RQ(L1)=ZERO
212             CONTINUE
                DO 213 L1=1,M
                  RQKEEP(L1)=RQ(L1)
213             CONTINUE
                MFOUT=MFIT
              ELSE
                MFOUT=0
              ENDIF
C             ----------------------------------------------------
C             2.2 Aposteriori estimate of damping factor
              FCNUMP = ZERO
              DO 2201 L1=1,N
                FCNUMP=FCNUMP+(DX(L1)/XW(L1))**2
2201          CONTINUE
              TH = FC-ONE
              FCDNM = ZERO
              DO 2202 L1=1,N
                FCDNM=FCDNM+((DXQA(L1)+TH*DX(L1))/XW(L1))**2
2202          CONTINUE
              FCH = FC*FC*HALF*DSQRT(FCNUMP/FCDNM)
C             ------------------------------------------------------
C             2.2.1 Computation of the numerator of damping
C                   factor predictor
              FCNMP2 = ZERO
              DO 221 L1=1,N
                FCNMP2=FCNMP2+(DXQA(L1)/XW(L1))**2
221           CONTINUE
              FCNUMP = FCNUMP*FCNMP2
C             ----------------------------------------------------
C             2.2.2 Decision criterion for Jacobian updating
C                   technique
              QGENJ = FC.LT.FCA.AND.NEW.GT.0.OR.FCH.LT.FC*SIGMA
     $                .OR.IRANKA.NE.N.OR.M.GT.N
     $                .OR. .NOT.QRANK1
              FCA = FC
              FCK2 = FCA
C             IF(FC.GT.FCMINH) IRANK = MINMN
              IRANK = MINMN
              IF(NONLIN.GT.1) FC = DMIN1(ONE,FCH)
C             ----------------------------------------------------
C             2.2.3 Denominator for kappa (SKAP) estimate
              SUMXK = ZERO
              DO 223 L1=1,N
                SUMXK = SUMXK+(DX(L1)/XW(L1))**2
223           CONTINUE
            ENDIF
          ENDIF
          QJCRFR =.FALSE.
C         --------------------------------------------------------
C         2.3 Jacobian matrix (stored to array AA(M2,N))
C         --------------------------------------------------------
C         2.3.1 Jacobian generation by routine JAC or
C               difference approximation (If QGENJ.EQ..TRUE.)
C               - or -
C               Rank-1 update of Jacobian (If QGENJ.EQ..FALSE.)
          IF(QGENJ)THEN
            NEW = 0
            IF (JACGEN.EQ.1) THEN
               IF (MPRTIM.NE.0) CALL MONON(2)
               CALL JAC(N,M,MCON,X,AA,IFAIL)
               IF (MPRTIM.NE.0) CALL MONOFF(2)
            ELSE
              IF (MPRTIM.NE.0) CALL MONON(2)
              IF (JACGEN.EQ.3) 
     $          CALL NCJCF(FCN,N,M,MCON,M,X,FMODEL,AA,XW,ETA,ETAMIN,
     $                     ETAMAX,ETADIF,CONV,NFCNJ,T3,IFAIL)
              IF (JACGEN.EQ.2) 
     $          CALL NCJAC(FCN, N, M, MCON, M, X, FMODEL, AA, XW,
     $                     AJDEL, AJMIN, NFCNJ, T3, IFAIL)
              IF (MPRTIM.NE.0) CALL MONOFF(2)
            ENDIF
            NJAC = NJAC + 1
C     Exit, If ...
            IF (JACGEN.EQ.1 .AND. IFAIL.LT.0) THEN
              IERR = 83
              GOTO 4299
            ENDIF
            IF (JACGEN.NE.1 .AND. IFAIL.NE.0) THEN
              IERR = 82
              GOTO 4299
            ENDIF
          ELSE
            NEW = NEW+1
            CALL NCRNK1(N,M,XW,DX,F,FA,T1,T3,AA,FCA)
          ENDIF
C         --------------------------------------------------------
C         2.3.2 Copy matrix to Jacobian work array A(M2,N)
          DO 232 L1=1,N
            DO 2321 L2=1,M2
               A(L2,L1)=AA(L2,L1)
2321        CONTINUE
232       CONTINUE
C         --------------------------------------------------------
C         2.4 Prepare solution of the linear system
C         --------------------------------------------------------
C         2.4.1 internal column scaling of matrix A
          DO 2410 K=1,N
            S1 =-XW(K)
            DO 2411 L1=1,M
              A(L1,K)=A(L1,K)*S1
2411        CONTINUE
2410      CONTINUE
C         ------------------------------------------------------
C         2.4.2 Row scaling of matrix A(M,N)
          IF (QSCALE) THEN
            CALL NCSCRF(N,M,MCON,A,FW)
          ELSE
            DO 242 L1=1,M
              FW(L1)=ONE
242         CONTINUE
          ENDIF
C         --------------------------------------------------------
C         2.4.3 Save and scale values of F(M)
          DO 243 L1=1,M
            FA(L1)=F(L1)
            T3(L1)=F(L1)*FW(L1)
243       CONTINUE
C         --------------------------------------------------------
C         2.4.4 Scaling of the linear residuum RQ(M)
          IF (M.GT.N.OR.IRANKA.NE.M) THEN
            DO 244 L1=1,M
              RQ(L1)=RQ(L1)*FW(L1)
244         CONTINUE
          ENDIF
          QNEXT = .FALSE.
C         --------------------------------------------------------
C         3 Central part of iteration step
C
C         Pseudo-rank reduction loop
C         ==========================
C         DO (Until)
3         CONTINUE
C           --------------------------------------------------------
C           3.1 Solution of linear (MFIT,N)-least squares problem
C               with MCON equality constraints
C           --------------------------------------------------------
C           3.1.1 Decomposition of (N,N)-matrix A
            IF (.NOT.QLINIT) THEN
              NIWLA = IWK(18)+1
              NRWLA = IWK(19)+1
              LIWL = LIWK-NIWLA+1
              LRWL = LRWK-NRWLA+1
            ENDIF
            IF (QREPET) THEN
              IWK(NIWLA) = 1
            ELSE
              IWK(NIWLA) = 0
            ENDIF
            COND1 = COND
            MCONH = MCON
            IF (MPRTIM.NE.0) CALL MONON(3)
            CALL NCFACT(N,M1,MCONH,N,ML,MU,A,QA,COND1,IRANK,IOPT,IFAIL,
     $                  LIWL,IWK(NIWLA),NILUSE,LRWL,RWK(NRWLA),NRLUSE)
            IF (MPRTIM.NE.0) CALL MONOFF(3)
            IF (.NOT.QLINIT) THEN
              NIWKFR = NIWKFR+NILUSE
              NRWKFR = NRWKFR+NRLUSE
C             Store lengths of currently required workspaces
              IWK(18) = NIWKFR-1
              IWK(19) = NRWKFR-1
              QLINIT = .TRUE.
            ENDIF
C     Exit Repeat If ...
            IF(IFAIL.NE.0) THEN
              IERR = 80
              GOTO 4299
            ENDIF
            SENS1 = RWK(NRWLA)
            CONDCO = RWK(NRWLA+N+1)
            IRANKC = IWK(NIWLA+1)
C           --------------------------------------------------------
C           3.1.2 Solution of linear (N,N)-system
            IF (MPRTIM.NE.0) CALL MONON(4)
            CALL NCFIT(N,M1,MCONH,N,ML,MU,A,QA,T3,T2,IRANK,IOPT,IFAIL,
     $                 LIWL,IWK(NIWLA),IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
            IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
            IF(IFAIL.NE.0)  THEN
              IERR = 81
              GOTO 4299
            ENDIF
            IF(.NOT.QREPET.AND.IRANK.NE.0)THEN
              DO 312 L1=1,M
                QU(L1)=T3(L1)
312           CONTINUE
            ENDIF
C           --------------------------------------------------------
C           3.2 Evaluation of scaled natural level function SUMX
C               scaled maximum error norm CONV
C               evaluation of (scaled) standard level function
C               DLEVF ( DLEVF only, if MPRMON.GE.2 )
C               and computation of ordinary Gauss-Newton corrections DX(
C               N)
            CALL NCLVLS(N,M,T2,XW,F,DX,CONV,SUMX,DLEVF,MPRMON)
            DO 32 L1=1,N
              XA(L1)=X(L1)
32          CONTINUE
            SUMXA = SUMX
            DLEVXA = DSQRT(SUMXA/DBLE(FLOAT(N)))
            CONVA = CONV
            DXANRM = WNORM(N,DX,XW)
C           --------------------------------------------------------
C           3.3 A - priori estimate of damping factor FC
            QREDU = .FALSE.
            IF(NITER.NE.0.AND.NONLIN.NE.1)THEN
CWei;              IF(NEW.EQ.0.AND.(IRANK.EQ.N.OR.IRANKA.EQ.N).OR.
CWei;      *           QREPET)THEN
              IF(NEW.EQ.0.OR.QREPET)THEN
C               ------------------------------------------------------
C               3.3.1 Comp. of the denominator of a-priori estimate
                IF (M.GT.N.OR.IRANKA.NE.M) THEN
                  CALL NCFIT(N,M1,MCONH,N,ML,MU,A,QA,RQ,DELXQ,IRANK,
     $                       IOPT,IFAIL,LIWL,IWK(NIWLA),IDUMMY,LRWL,
     $                       RWK(NRWLA),IDUMMY)
C     Exit Repeat If ...
                  IF(IFAIL.NE.0)  THEN
                    IERR = 81
                    GOTO 4299
                  ENDIF
                  FCDNM = ZERO
                  DO 3312 L1=1,N
                    FCDNM=FCDNM+((DX(L1)-DXQ(L1))/XW(L1)-DELXQ(L1))**2
3312              CONTINUE
                ELSE
                  FCDNM = ZERO
                  DO 3313 L1=1,N
                    FCDNM=FCDNM+((DX(L1)-DXQ(L1))/XW(L1))**2
3313              CONTINUE
                ENDIF
                IF(IRANK.NE.N)THEN
C                 ------------------------------------------------
C                 3.3.2 Rank-deficient case (if previous rank
C                           was full) computation of the projected
C                       denominator of a-priori estimate
                  DO 332 L1=1,N
                    T1(L1)=DXQA(L1)/XW(L1)
332               CONTINUE
C                 Norm of projection of reduced component T1(N)
                  CALL NCPRJN(N,IRANK,DELQ,T1,RWK(NRWLA+1),T2,QA,
     $                       IWK(NIWLA+2))
                  FCDNM = FCDNM-DELQ
                ENDIF
                FCDNM = FCDNM*SUMX
C               ------------------------------------------------------
C               3.3.3 New damping factor
                IF(FCDNM.GT.FCNUMP*FCMIN2)THEN
                  DMUE = FCA*DSQRT(FCNUMP/FCDNM)
                  FC = DMIN1(DMUE,ONE)
                  IF (DMUE.GT.TEN) IFCCNT = IFCCNT + 1
                ELSE
                  FC = ONE
C$Test-begin
                  DMUE = -1.0D0
C$Test-end
                  IF (FCA.GE.ONE) IFCCNT = IFCCNT+1
                ENDIF
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,33201) IFCCNT, FC, FCA, DMUE, FCNUMP,
     $                               FCDNM
33201             FORMAT (/, ' +++ apriori estimate +++', /,
     $                   ' IFCCNT = ', I10, 8X,'  FC     = ', D18.10, /,
     $                   ' FCA    = ', D18.10, '  DMUE   = ', D18.10, /,
     $                   ' FCNUMP = ', D18.10, '  FCDNM  = ', D18.10, /,
     $                       ' ++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
CWei            FC = DMAX1(FC,FCMIN)
                IF (QBDAMP) THEN
                  FCBH = FCA*FCBND
                  IF (FC.GT.FCBH) THEN
                    FC = FCBH
                    IF (MPRMON.GE.4)
     $                WRITE(LUMON,*)' *** incr. rest. act. (a prio) ***'
                  ENDIF
                  FCBH = FCA/FCBND
                  IF (FC.LT.FCBH) THEN
                    FC = FCBH
                    IF (MPRMON.GE.4)
     $                WRITE(LUMON,*)' *** decr. rest. act. (a prio) ***'
                  ENDIF
                ENDIF
              ENDIF
              QREDU = FC.LT.FCMIN
            ENDIF
            QREPET = .FALSE.
            IF(.NOT.QREDU)THEN
C             --------------------------------------------------------
C             3.4 Save natural level for later computations of
C                 corrector and print iterate
              FCNUMK = SUMX
              IF (MPRMON.GE.2) THEN
                IF (MPRTIM.NE.0) CALL MONON(5)
                CALL NCPRV1(DLEVF,DLEVXA,FCA,NITER,NEW,IRANK,MPRMON,
     $                      LUMON,QMIXIO)
                IF (MPRTIM.NE.0) CALL MONOFF(5)
              ENDIF
              NRED = 0
C             Damping-factor reduction loop
C             ================================
C             DO (Until)
34            CONTINUE
C               ------------------------------------------------------
C               3.5 Preliminary new iterate
                DO 35 L1=1,N
                  X(L1)=XA(L1)+DX(L1)*FC
35              CONTINUE
C               -----------------------------------------------------
C               3.5.2 Exit, if problem is specified as being linear
C     Exit Repeat If ...
                IF( NONLIN.EQ.1 )THEN
                  IERR = 0
                  GOTO 4299
                ENDIF
C               ------------------------------------------------------
C               3.6.1 Computation of the residual vector
                IF (MPRTIM.NE.0) CALL MONON(1)
                CALL FCN(N,M,MCON,X,FMODEL,IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(1)
                NFCN = NFCN+1
C     Exit, if ...
                IF (IFAIL.LT.0) THEN
                  IERR = 82
                  GOTO 4299
                ENDIF
                IF(IFAIL.EQ.1 .OR. IFAIL.EQ.2) THEN
                  IF (IFAIL.EQ.1) THEN
                    FCREDU = HALF
                  ELSE
                    FCREDU = F(1)
C     Exit, If ...
                    IF (FCREDU.LE.0 .OR. FCREDU.GE.1) THEN
                      IERR = 83
                      GOTO 4299
                    ENDIF
                  ENDIF
                  IF (MPRMON.GE.2) THEN
36101               FORMAT(8X,I2,' FCN could not be evaluated ',
     $                     13X,F5.3,4X,I2,6X,I4)
                    WRITE(LUMON,36101)NITER,FC,NEW,IRANK
                  ENDIF
                  FCH = FC
                  FC = FCREDU*FC
                  IF (FCH.GT.FCMIN) FC = DMAX1(FC,FCMIN)
                  IF (QBDAMP) THEN
                    FCBH = FCH/FCBND
                    IF (FC.LT.FCBH) THEN
                      FC = FCBH
                      IF (MPRMON.GE.4) WRITE(LUMON,*)
     $                   ' *** decr. rest. act. (FCN redu.) ***'
                    ENDIF
                  ENDIF
                  IF (FC.LT.FCMIN) THEN
                    IERR = 3
                    GOTO 4299
                  ENDIF  
C     Break DO (Until) ...
                  GOTO 3109
                ENDIF
                DO 3611 L1=1,MFIT
                  F(L1+MCON)=FMODEL(L1+MCON)-FI(L1)
3611            CONTINUE
                DO 3612 L1=1,MCON
                  F(L1)=FMODEL(L1)
3612            CONTINUE
                DO 3613 L1=1,M
                  T3(L1)=F(L1)*FW(L1)
3613            CONTINUE
C               ------------------------------------------------------
C               3.6.2 Solution of linear (MFIT,N)-least squares problem
C                     with MCON equality constraints
                IF (QREPET) THEN
                  IWK(NIWLA) = 1
                ELSE
                  IWK(NIWLA) = 0
                ENDIF
                IF (MPRTIM.NE.0) CALL MONON(4)
                CALL NCFIT(N,M1,MCONH,N,ML,MU,A,QA,T3,T2,IRANK,IOPT,
     $                      IFAIL,LIWL,IWK(NIWLA),IDUMMY,LRWL,
     $                      RWK(NRWLA),IDUMMY)
                IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
                IF(IFAIL.NE.0)  THEN
                  IERR = 81
                  GOTO 4299
                ENDIF
C               ------------------------------------------------------
C               3.6.3 Evaluation of scaled natural level function
C                     SUMX
C                     scaled maximum error norm CONV and evaluation
C                     of (scaled) standard level function DLEVF
                CALL NCLVLS(N,M,T2,XW,F,DXQ,CONV,SUMX,DLEVF,MPRMON)
                DXNRM = WNORM(N,DXQ,XW)
C               ------------------------------------------------------
C               3.6.4 Convergence test
C     Exit Repeat If ...
                IF (ITERM.EQ.0) THEN
                  IF( (DXNRM.LE.RTOL .AND. IFCCNT.GE.3) .OR. 
     $                DXANRM.LE.RTOL) THEN
                    IERR = 0
                    GOTO 4299
                  ENDIF
                ELSE IF (ITERM.EQ.1) THEN
                  IF( DXNRM.LE.RTOL .AND. IFCCNT.GE.3 ) THEN
                    IERR = 0
                    GOTO 4299
                  ENDIF
                ELSE IF (ITERM.EQ.2) THEN
                  IF( DXANRM.LE.RTOL ) THEN
                    IERR = 0
                    GOTO 4299
                  ENDIF
                ENDIF
C           
                FCA = FC
C               ------------------------------------------------------
C               3.7 Natural monotonicity test
                IF(SUMX.GT.SUMXA)THEN
C                 ----------------------------------------------------
C                 3.8 Output of iterate
                  IF(MPRMON.GE.3) THEN
                    IF (MPRTIM.NE.0) CALL MONON(5)
                    CALL NCPRV2(DLEVF,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                          NITER,MPRMON,LUMON,QMIXIO,'*')
                    IF (MPRTIM.NE.0) CALL MONOFF(5)
                  ENDIF
C                 ----------------------------------------------------
C                 3.9 Evaluation of reduced damping factor
                  TH = FCA-ONE
                  FCDNM = ZERO
                  DO 39 L1=1,N
                    FCDNM=FCDNM+((DXQ(L1)+TH*DX(L1))/XW(L1))**2
39                CONTINUE
                  FC = FCA*FCA*HALF*DSQRT(FCNUMK/FCDNM)
                  IF (QBDAMP) THEN
                    FCBH = FCA/FCBND
                    IF (FC.LT.FCBH) THEN
                      FC = FCBH
                      IF (MPRMON.GE.4) WRITE(LUMON,*) 
     $                    ' *** decr. rest. act. (a post) ***'
                    ENDIF
                  ENDIF
                  NCORR = NCORR+1
                  NRED = NRED+1
                  IFCCNT = 0
C                 ----------------------------------------------------
C                 3.10 Rank reduction, if damping factor to small
                  QREDU  = FC.LT.FCMIN.OR.NEW.GT.0.AND.NRED.GT.1
                ELSE
                  QNEXT = .TRUE.
                ENDIF
3109          CONTINUE
              IF(.NOT.(QNEXT.OR.QREDU)) GOTO  34
C             UNTIL ( expression - negated above)
C             End of damping-factor reduction loop
C           =======================================
            ENDIF
            IF(QREDU)THEN
C             ------------------------------------------------------
C             3.11 Restore former values for repeting step
C                  step
              NREJR1 = NREJR1+1
              DO 3111 L1=1,N
                X(L1)=XA(L1)
3111          CONTINUE
              DO 3112 L1=1,M
                F(L1)=FA(L1)
3112          CONTINUE
              DO 3113 L1=1,N
                DXQ(L1)=DXQA(L1)
3113          CONTINUE
              IF(MPRMON.GE.2)THEN
31130           FORMAT(8X,I2,' Not accepted damping ',
     $                 'factor',13X,F5.3,4X,I2,6X,I4)
                WRITE(LUMON,31130)NITER,FC,NEW,IRANK
              ENDIF
              IFCCNT = 0
              FCA = FCK2
              IF(NITER.EQ.0)THEN
                FC = FCMIN
              ENDIF
              IF(NEW.GT.0)THEN
                QGENJ = .TRUE.
                QJCRFR = .TRUE.
                QREDU = .FALSE.
                FC = FCH
                IRANK = MINMN
                DO 3114 L1=1,M
                  RQ(L1)=RQKEEP(L1)
3114            CONTINUE
              ELSE
C               ------------------------------------------------
C               3.12 Pseudo-rank reduction
                QREPET = .TRUE.
                DO 42 L1=1,M
                  T3(L1)=QU(L1)
42              CONTINUE
                IRANK = IRANK-1
                IF(IRANK.EQ.0)THEN
                  IERR = 3
                  GOTO 4299
                ENDIF
              ENDIF
            ENDIF
          IF(.NOT.(.NOT.QREDU)) GOTO  3
C         UNTIL ( expression - negated above)
C
C         End of pseudo-rank reduction loop
C         =================================
          IF (QNEXT) THEN
C           ------------------------------------------------------
C           4 Preparations to start the following iteration step
C           ------------------------------------------------------
C           4.1 Print values
            IF(MPRMON.GE.3) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL NCPRV2(DLEVF,DSQRT(SUMX/DBLE(FLOAT(N))),FC,NITER+1,
     $                    MPRMON,LUMON,QMIXIO,'*')
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
C           Print the natural level of the current iterate and return
C           it in one-step mode
            SUMXS = SUMX
            SUMX = SUMXA
            IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL NCSOUT(N,MFOUT,XA,RQKEEP(MCON+1),2,IOPT,RWK,LRWK,IWK,
     $                    LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL NCSOUT(N,MFOUT,XA,DUMMY,1,IOPT,RWK,LRWK,IWK,
     $                    LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            NITER = NITER+1
C     Exit Repeat If ...
            IF(NITER.GE.NITMAX)THEN
              IERR = 2
              GOTO 4299
            ENDIF
            FCA = FC
C           ------------------------------------------------------
C           4.2 Return, if in one-step mode
C Exit Subroutine If ...
            IF (MODE.EQ.1) THEN
              IWK(18)=NIWLA-1
              IWK(19)=NRWLA-1
              IOPT(1)=1
              RETURN
            ENDIF
          ENDIF
        GOTO 2
C       End Repeat
4299    CONTINUE
C
C       End of main iteration loop
C       ==========================
C       ----------------------------------------------------------
C       9 Exits
C       ----------------------------------------------------------
C       9.1 Solution exit
        IF(IERR.EQ.0)THEN
          IF (NONLIN.NE.1) THEN
            DO 910 L1=1,N
              X(L1)=X(L1)+DXQ(L1)
910         CONTINUE
            IF(IRANK.LT.MINMN)THEN
              IERR = 1
            ENDIF
C           Print final monitor output
            IF(MPRMON.GE.2) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL NCPRV2(DLEVF,DSQRT(SUMX/DBLE(FLOAT(N))),FC,NITER+1,
     $                    MPRMON,LUMON,QMIXIO,'*')
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            IF(MPRMON.GE.1.AND.IERR.EQ.0) THEN
91001         FORMAT(///' Solution of nonlinear least squares problem ',
     $        'obtained within ',I3,' iteration steps')
              WRITE(LUMON,91001) NITER+1
            ENDIF
          ELSE
            IF(MPRMON.GE.1) THEN
91002         FORMAT(///' Least squares solution of linear system ',
     $        'of equations obtained by NLSCON',//,' No estimate ',
     $        'available for the achieved relative accuracy')
                WRITE(LUMON,91002)
            ENDIF
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2 Fail exit messages
C       ----------------------------------------------------------
C       9.2.1 Termination at stationary point
        IF(IERR.EQ.1.AND.MPRERR.GE.1)THEN
92101     FORMAT(/,' Iteration terminates at stationary point',/)
          WRITE(LUERR,92101)
        ENDIF
        RWK(31) = -ONE
        RWK(32) = -ONE
        IF((IERR.EQ.0.OR.IERR.EQ.1).AND.NONLIN.NE.1.AND.MPRMON.GE.1)THEN
          IF (SUMXK .LT. TOLMIN) THEN
            SKAP = -ONE
          ELSE IF (SUMXA .LT. TOLMIN) THEN
            SKAP = DSQRT(TOLMIN)
          ELSE
            SKAP = DSQRT(SUMXA/SUMXK)
          ENDIF
          IF (SKAP.GE.0) THEN
92202       FORMAT(/,' Incompatibility factor kappa',D10.3,2X,/)
            WRITE(LUMON,92202)SKAP
          ELSE 
92203       FORMAT(/,' Incompatibility factor kappa not available')
            WRITE(LUMON,92203)
          ENDIF
          IF ( IERR.EQ.0 .AND. SKAP.LT.ONE ) THEN
92204       FORMAT(/,' Achieved relative accuracy ',D10.3,/)
            IF (SKAP.GE.0) THEN
              PREC=DMAX1((SKAP/(ONE-SKAP))*
     $             DSQRT(SUMXA/DBLE(FLOAT(N))),EPMACH)
            ELSE 
              PREC=EPMACH
            ENDIF
            WRITE(LUMON,92204) PREC
          ENDIF
          RWK(31) = PREC
          RWK(32) = SKAP
        ENDIF
        RTOL = RWK(31)
C       ----------------------------------------------------------
C       9.2.2 Termination after more than NITMAX iterations
        IF(IERR.EQ.2.AND.MPRERR.GE.1)THEN
92210     FORMAT(/,' Iteration terminates after NITMAX ',
     $    '=',I3,'  Iteration steps')
          WRITE(LUERR,92210)NITMAX
        ENDIF
C       ----------------------------------------------------------
C       9.2.3 Gauss-Newton method fails to converge
        IF(IERR.EQ.3.AND.MPRERR.GE.1)THEN
92301     FORMAT(/,' Gauss-Newton method fails to converge')
          WRITE(LUERR,92301)
        ENDIF
C       ----------------------------------------------------------
C       9.2.5 Error exit due to linear solver routine NCFACT
        IF(IERR.EQ.80.AND.MPRERR.GE.1)THEN
92501     FORMAT(/,' Error ',I5,' signalled by linear solver NCFACT')
          WRITE(LUERR,92501) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.6 Error exit due to linear solver routine NCFIT
        IF(IERR.EQ.81.AND.MPRERR.GE.1)THEN
92601     FORMAT(/,' Error ',I5,' signalled by linear solver NCFIT')
          WRITE(LUERR,92601) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.7 Error exit due to fail of user function FCN
        IF(IERR.EQ.82.AND.MPRERR.GE.1)THEN
92701     FORMAT(/,' Error ',I5,' signalled by user function FCN')
          WRITE(LUERR,92701) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.7 Error exit due to fail of user function JAC
        IF(IERR.EQ.83.AND.MPRERR.GE.1)THEN
92801     FORMAT(/,' Error ',I5,' signalled by user function JAC')
          WRITE(LUERR,92801) IFAIL
        ENDIF
        IF(IERR.GE.80.AND.IERR.LE.83) IWK(23) = IFAIL
        IF ((IERR.EQ.82.OR.IERR.EQ.83).AND.NITER.LE.1.AND.MPRERR.GE.1)
     $  THEN
          WRITE (LUERR,92810)
92810     FORMAT(' Try to find a better initial guess for the solution')
        ENDIF
C       ----------------------------------------------------------
C       9.3 Common exit
        IF(MPRMON.GE.1)THEN
93001     FORMAT(/,3X,'Subcondition ( 1,',I4,') of ',A,' part',
     $    1X,D10.3)
93003     FORMAT(/,3X,'Subcondition ( ',I4,',',I4,') of ',A,' part',
     $    1X,D10.3)
93002     FORMAT(/,3X,'Sensitivity ( lsq )',1X,D10.3,2X,/)
          IF (MCON.GT.0) THEN
            WRITE(LUMON,93001) IRANKC,'constrained',CONDCO
            WRITE(LUMON,93003) IRANKC+1,IRANK,'least squares',COND1
          ELSE
            WRITE(LUMON,93001)IRANK,'least squares',COND1
          ENDIF
          WRITE(LUMON,93002)SENS1
        ENDIF
        SUMXS = SUMX
        SUMX = SUMXA
        IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
          IF (MPRTIM.NE.0) CALL MONON(5)
          CALL NCSOUT(N,MFOUT,XA,RQKEEP(MCON+1),2,IOPT,RWK,LRWK,IWK,
     $                LIWK,MPRSOL,LUSOL)
          IF (MPRTIM.NE.0) CALL MONOFF(5)
        ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
          IF (MPRTIM.NE.0) CALL MONON(5)
          CALL NCSOUT(N,MFOUT,XA,DUMMY,1,IOPT,RWK,LRWK,IWK,
     $                LIWK,MPRSOL,LUSOL)
          IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
        NITER = NITER+1
        IF(MPRSOL.GE.1)THEN
C         Print Solution or final iteration vector
          IF(IERR.EQ.0)THEN
             MODEFI = 3
          ELSE
             MODEFI = 4
          ENDIF
          IF (MPRTIM.NE.0) CALL MONON(5)
          CALL NCSOUT(N,0,X,DUMMY,MODEFI,IOPT,RWK,LRWK,IWK,
     $                LIWK,MPRSOL,LUSOL)
          IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
C       Return the latest internal scaling to XSCAL
        DO 93 I=1,N
          XSCAL(I)=XW(I)
93      CONTINUE
C       ----------------------------------------------------------
C       9.4 Optional computation of the statistical analysis
C           for the least squares problem solution
        IF (IOPT(21).EQ.1 .OR. IOPT(22).EQ.1) THEN
          IOPT(21) = 1
          IF (MPRTIM.NE.0) CALL MONON(1)
          CALL FCN(N,M,MCON,X,FMODEL,IFAIL)
          IF (MPRTIM.NE.0) CALL MONOFF(1)
          NFCN = NFCN+1
C     Exit, if ...
          IF (IFAIL.NE.0) THEN
            IERR = 182
            GOTO 9900
          ENDIF
          DO 941 L1=1,MFIT
            F(L1+MCON)=FMODEL(L1+MCON)-FI(L1)
941       CONTINUE
          DO 942 L1=1,MCON
            F(L1)=FMODEL(L1)
942       CONTINUE
          IF (JACGEN.EQ.1) THEN
             IF (MPRTIM.NE.0) CALL MONON(2)
             CALL JAC(N,M,MCON,X,AA,IFAIL)
             IF (MPRTIM.NE.0) CALL MONOFF(2)
          ELSE
            IF (MPRTIM.NE.0) CALL MONON(2)
            IF (JACGEN.EQ.3) 
     $        CALL NCJCF(FCN,N,M,MCON,M,X,FMODEL,AA,XW,ETA,ETAMIN,
     $                   ETAMAX,ETADIF,CONV,NFCNJ,T3,IFAIL)
            IF (JACGEN.EQ.2) 
     $        CALL NCJAC(FCN, N, M, MCON, M, X, FMODEL, AA, XW,
     $                   AJDEL, AJMIN, NFCNJ, T3, IFAIL)
            IF (MPRTIM.NE.0) CALL MONOFF(2)
          ENDIF
          NJAC = NJAC + 1
C     Exit, If ...
          IF (JACGEN.EQ.1 .AND. IFAIL.LT.0) THEN
            IERR = 183
            GOTO 9900
          ENDIF
          IF (JACGEN.NE.1 .AND. IFAIL.NE.0) THEN
            IERR = 182
            GOTO 9900
          ENDIF
          DO 943 L1=1,N
            DO 9430 L2=1,M2
               A(L2,L1)=AA(L2,L1)
9430        CONTINUE
943       CONTINUE
          IWK(NIWLA) = 0
          COND1 = COND
          MCONH = MCON
          IRANK =MINMN
          IF (MPRTIM.NE.0) CALL MONON(3)
          CALL NCFACT(N,M1,MCONH,N,ML,MU,A,QA,COND1,IRANK,IOPT,IFAIL,
     $                LIWL,IWK(NIWLA),NILUSE,LRWL,RWK(NRWLA),NRLUSE)
          IF (IFAIL.NE.0) THEN
            IERR = 180
            GOTO 9900
          ENDIF
          IRANKC = IWK(NIWLA+1)
          IF (MPRTIM.NE.0) CALL MONOFF(3)
          CALL STACON (N,MCON,M,MFIT,N,MCON,FI,X,
     $                 FMODEL(MCON+1),A,IRANKC,IRANK,RWK(NRWLA+1),
     $                 IWK(NIWLA+2),QA,RWK(NRWLA+N+1),IERR,VCV,
     $                 RWK(NRWLA+2*N+1),T3,XL,XR,RWK(50),
     $                 MPRMON,LUMON)
        ENDIF
        RETURN
9900    CONTINUE
        IF (MPRERR.GT.0) THEN
          IF (IERR.EQ.180) WRITE(LUERR,99001) 'DECCON'
          IF (IERR.EQ.182) WRITE(LUERR,99001) 'FCN'
          IF (IERR.EQ.183) WRITE(LUERR,99001) 'JAC'
        ENDIF
99001   FORMAT(/,1X,'Computation of the statistical analysis skipped',/,
     $           1X,'since ',A,' failed for the final Gauss-Newton ',
     $              'iterate')
C       End of exits
C       End of subroutine NCINT
      RETURN
      END
C
      SUBROUTINE NCSCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,LRWK,RWK)
C*    Begin Prologue SCALE
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N),XA(N),XW(N)
      INTEGER ISCAL
      LOGICAL QINISC
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C    
C     S C A L E : To be used in connection with NLSCON .
C       Computation of the internal scaling vector XW used for the
C       Jacobian matrix, the iterate vector and it's related
C       vectors - especially for the solution of the linear system
C       and the computations of norms to avoid numerical overflow.
C
C*    Input parameters
C     ================
C
C     N         Int     Number of unknowns
C     X(N)      Dble    Current iterate
C     XA(N)     Dble    Previous iterate
C     XSCAL(N)  Dble    User scaling passed from parameter XSCAL
C                       of interface routine NLSCON
C     ISCAL     Int     Option ISCAL passed from IOPT-field
C                       (for details see description of IOPT-fields)
C     QINISC    Logical = .TRUE.  : Initial scaling
C                       = .FALSE. : Subsequent scaling
C     LRWK      Int     Length of real workspace
C     RWK(LRWK) Dble    Real workspace (see description above)
C
C*    Output parameters
C     =================
C
C     XW(N)     Dble   Scaling vector computed by this routine
C                      All components must be positive. The follow-
C                      ing relationship between the original vector
C                      X and the scaled vector XSCAL holds:
C                      XSCAL(I) = X(I)/XW(I) for I=1,...N
C
C*    Subroutines called: D1MACH
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL D1MACH
      INTRINSIC DABS,DMAX1
      INTEGER L1
      DOUBLE PRECISION D1MACH,HALF
      PARAMETER (HALF=0.5D0)
      SMALL  = D1MACH(6)
C*    Begin
      DO 1 L1=1,N
        IF (ISCAL.EQ.1) THEN
          XW(L1) = XSCAL(L1)
        ELSE
          XW(L1)=DMAX1(XSCAL(L1),(DABS(X(L1))+DABS(XA(L1)))*HALF,SMALL)
        ENDIF
1     CONTINUE
C     End of subroutine NCSCAL
      RETURN
      END
C
      SUBROUTINE NCSCRF(N,M,MCON,A,FW)
C*    Begin Prologue SCROWF
      INTEGER N,M,MCON
      DOUBLE PRECISION A(M,N),FW(M)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S C R O W F : Row Scaling of a (M,N)-matrix in full storage
C                   mode
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N           Int    Number of columns of the matrix
C       M           Int    Number of rows of the matrix -
C                          leading dimension of A
C       MCON        Int    Number of rows to be automatically rescaled
C                          (the rows 1 to MCON will be rescaled)
C     * A(M,N)      Dble   Matrix to be scaled
C     * FW(M)       Dble   Scaling vector 
C
C*    Output parameters
C     =================
C
C       A(M,N)      Dble   The scaled matrix
C       FW(M)       Dble   Row scaling factors - FW(i) contains
C                          the factor by which the i-th row of A
C                          has been multiplied -
C                          unaltered in components MCON+1,...,M
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER J,K
      DOUBLE PRECISION S1,S2
C*    Begin
      DO 1 K=1,MCON
        S1=ZERO
        DO 2 J=1,N
          S2=DABS(A(K,J))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 J=1,N
            A(K,J)=A(K,J)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
      DO 4 K=MCON+1,M
        S1 = FW(K)
        DO 5 J=1,N
          A(K,J)=A(K,J)*S1
5       CONTINUE
4     CONTINUE
C     End of subroutine N1SCRF
      RETURN
      END
C
      SUBROUTINE NCFACT(N,M,MCON,LDAINV,ML,MU,A,AINV,COND,IRANK,IOPT,
     $IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK)
C*    Begin Prologue FACT
      INTEGER N,M,MCON,LDAINV,ML,MU
      DOUBLE PRECISION A(M,N),AINV(LDAINV,N)
      DOUBLE PRECISION COND
      INTEGER IRANK
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LAIWK,LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     F A C T : Call linear algebra subprogram for factorization of
C               a (N,N)-matrix with rank decision and casual compu-
C               tation of the rank deficient pseudo-inverse matrix
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C     N             Int    Number of parameters to be estimated
C     M             Int    Number of observations + equal. constraints
C     MCON          Int    Number of equality constraints
C     LDAINV        Int    Leading dimension of the matrix array AINV
C     ML            Int    Lower bandwidth of the matrix (only for
C                          banded systems)
C     MU            Int    Upper bandwidth of the matrix (only for
C                          banded systems)
C   * A(M,N)        Dble   Matrix storage.
C   * COND          Dble   On Input, COND holds the maximum permitted
C                          subcondition for the prescribed rank
C                          On Output, it holds the estimated subcon-
C                          dition of A
C     IOPT(50)      Int    Option vector passed from NLSCON
C
C*    Output parameters
C     =================
C
C     AINV(LDAINV,N) Dble   If matrix A is rank deficient, this array
C                           holds the pseudo-inverse of A
C     IFAIL          Int    Error indicator returned by this routine:
C                           = 0 matrix decomposition successfull
C                           =10 supplied (integer) workspace too small
C
C*    Workspace parameters
C     ====================
C
C     LIWK           Int    Length of integer workspace passed to this
C                           routine (In)
C     IWK(LIWK)      Int    Integer Workspace supplied for this routine
C     LAIWK          Int    Length of integer Workspace used by this 
C                           routine (out)       
C     LRWK           Int    Length of real workspace passed to this
C                           routine (In)                  
C     RWK(LRWK)      Dble   Real Workspace supplied for this routine
C     LARWK          Int    Length of real Workspace used by this 
C                           routine (out)
C
C*    Subroutines called:  DECCON
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL DECCON
      INTRINSIC DABS
      INTEGER IREPET, MPRERR, LUERR
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C*    Begin
      MPRERR=IOPT(11)
      LUERR=IOPT(12)
      LAIWK = N+2
      LARWK = 2*N+1
      IF (LIWK.GE.LAIWK.AND.LRWK.GE.LARWK) THEN
        IREPET = -IWK(1)
        IF (IREPET.EQ.0)  IWK(2) = MCON
        CALL DECCON(A,M,N,MCON,M,N,IWK(2),IRANK,COND,RWK(2),IWK(3),
     $              IREPET,AINV,RWK(N+2),IFAIL)
        IF (IFAIL.EQ.-2 .AND. MPRERR.GT.0) WRITE(LUERR,10001)
10001   FORMAT(1X,
     $       'DECCON failed to compute rank-deficient QR-decomposition',
     $        /)
        IF(IRANK.NE.0)THEN
          RWK(1) = DABS(RWK(IRANK+2))
        ELSE
          COND = ONE
          RWK(1) = ZERO
        ENDIF
      ELSE
        IFAIL = 10
10002   FORMAT(/,' Insuffient workspace for linear solver,',
     $         ' at least needed more needed : ',/,
     $         ' ',A,' workspace : ',I4)
        IF (LIWK.LT.LAIWK.AND.MPRERR.GT.0) 
     $    WRITE(LUERR,10002) 'Integer',LAIWK-LIWK
        IF (LRWK.LT.LARWK.AND.MPRERR.GT.0) 
     $    WRITE(LUERR,10002) 'Double',LARWK-LRWK
      ENDIF
      RETURN
      END
C
      SUBROUTINE NCFIT(N,M,MCON,LDAINV,ML,MU,A,AINV,B,Z,IRANK,IOPT,
     $IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK)
C*    Begin Prologue FIT
      INTEGER N,M,MCON,LDAINV,ML,MU
      DOUBLE PRECISION A(M,N),AINV(LDAINV,N)
      DOUBLE PRECISION B(M),Z(N)
      INTEGER IRANK
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK,LAIWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     F I T : Call linear algebra subprogram for (least squares) 
C             solution of the linear system A*Z = B
C
C*    Parameters
C     ==========
C
C     N,M,MCON,LDAINV,ML,MU,A,AINV,IRANK,IOPT,IFAIL,LIWK,IWK,LAIWK,
C     LRWK,RWK,LARWK :
C                        See description for subroutine NCFACT.          
C     B          Dble    In:  Right hand side of the linear system
C                        Out: Rhs. transformed to the upper trian-
C                             gular part of the linear system
C     Z          Dble    Out: Solution of the linear system
C
C     Subroutines called: SOLCON
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL SOLCON
      INTEGER IREPET
C*    Begin
      IREPET = -IWK(1)
      CALL SOLCON(A,M,N,MCON,M,N,Z,B,IWK(2),IRANK,RWK(2),IWK(3),
     $            IREPET,AINV,RWK(N+2))
      IFAIL = 0
      RETURN
      END
C
      SUBROUTINE NCLVLS(N,M,DX1,XW,F,DXQ,CONV,SUMX,DLEVF,MPRMON)
C*    Begin Prologue LEVELS
      INTEGER N,M,MPRMON
      DOUBLE PRECISION DX1(N),XW(N),F(M),DXQ(N)
      DOUBLE PRECISION CONV,SUMX,DLEVF
C     ------------------------------------------------------------
C
C*    Summary :
C
C     L E V E L S : To be used in connection with NLSCON .
C     provides descaled solution, error norm and level functions
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N              Int    Number of parameters to be estimated
C       M              Int    Number of measurements + eq. constraints
C       DX1(N)         Dble   array containing the scaled Gauss-Newton
C                             correction
C       XW(N)          Dble   Array containing the scaling values
C       F(M)           Dble   Array containing the residuum
C
C*    Output parameters
C     =================
C
C       DXQ(N)         Dble   Array containing the descaled Gauss-Newton
C                             correction
C       CONV           Dble   Scaled maximum norm of the Gauss-Newton
C                             correction
C       SUMX           Dble   Scaled natural level function value
C       DLEVF          Dble   Standard level function value (only
C                             if needed for print)
C       MPRMON         Int    Print information parameter (see
C                             driver routine NLSCON )
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1
      DOUBLE PRECISION S1
C*    Begin
C     ------------------------------------------------------------
C     1 Descaling of solution DX1 ( stored to DXQ )
      DO 1 L1=1,N
        DXQ(L1)=DX1(L1)*XW(L1)
1     CONTINUE
C     ------------------------------------------------------------
C     2 Evaluation of scaled natural level function SUMX and
C       scaled maximum error norm CONV
      CONV = ZERO
      DO 20 L1=1,N
        S1 = DABS(DX1(L1))
        IF(S1.GT.CONV) CONV=S1
20    CONTINUE
      SUMX = ZERO
      DO 21 L1=1,N
        SUMX = SUMX+DX1(L1)**2
21    CONTINUE
C     ------------------------------------------------------------
C     3 Evaluation of (scaled) standard level function DLEVF (only
C       if needed for print)
      IF(MPRMON.GE.2)THEN
        DLEVF = ZERO
        DO 3 L1=1,M
          DLEVF = DLEVF+F(L1)**2
3       CONTINUE
        DLEVF = DSQRT(DLEVF/DBLE(FLOAT(M)))
      ENDIF
C     End of subroutine NCLVLS
      RETURN
      END
C
      SUBROUTINE NCJAC (FCN, N, M, MCON, LDA, X, FX, A, YSCAL, AJDEL,
     $                  AJMIN, NFCN, FU, IFAIL)
C* Begin Prologue NCJAC
      EXTERNAL FCN
      INTEGER N, M, MCON, LDA
      DOUBLE PRECISION X(N), FX(M), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(M)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Evaluation of a dense Jacobian matrix using finite difference
C  approx. adapted for use in nonlinear least squares solver NLSCON
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C     FCN(N,M,MCON,X,F,IFAIL) 
C                      Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       M              Int    Number of measurement-vector
C                             components plus number of equality
C                             constraints (input)
C       MCON           Int    Number of equality constraints (input)
C       X(N)           Dble   Vector of parameters (input)
C       F(M)           Dble   Vector of equality constraints and
C                             measurement fitting values -
C                             the first MCON-components belonging
C                             to the constraints (output).
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             Whenever a nonzero value is returned
C                             by FCN routine NCJAC is terminated
C                             immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of columns of the Jacobian
C  M          Int     Number of rows of the Jacobian
C  MCON       Int     Number of equality constraints (MCON.LE.M)
C                     (only passed to FCN)
C  LDA        Int     Leading dimension of A (LDA .GE. M)
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(M)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
C  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(M)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called Subroutines:
C  -------------------
C
C     EXTERNAL FCN
      INTRINSIC DABS, DMAX1, DSIGN
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K
      DOUBLE PRECISION U, W
C
C* Begin
C
      IFAIL = 0
      DO 1 K = 1,N
         W = X(K)
         U = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
         X(K) = W + U
C
         CALL FCN (N, M, MCON, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99
C
         X(K) = W
         DO 11 I = 1,M
            A(I,K) = (FU(I) - FX(I)) / U  
 11      CONTINUE
 1    CONTINUE
C
99    CONTINUE
      RETURN
C
C
C* End of NCJAC
C
      END
      SUBROUTINE NCJCF (FCN, N, M, MCON, LDA, X, FX, A, YSCAL, ETA,
     $     ETAMIN, ETAMAX, ETADIF, CONV, NFCN, FU, IFAIL)
C* Begin Prologue NCJCF
      EXTERNAL FCN
      INTEGER N, M, MCON, LDA
      DOUBLE PRECISION X(N), FX(M), A(LDA,N), YSCAL(N), ETA(N),
     $     ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(M)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Approx. of dense Jacobian matrix for nonlinear least squares solver
C  NLSCON with feed-back control of discretization and rounding errors
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   January 1991
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C     FCN(N,M,MCON,X,F,IFAIL) 
C                      Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       M              Int    Number of measurement-vector
C                             components plus number of equality
C                             constraints (input)
C       MCON           Int    Number of equality constraints (input)
C       X(N)           Dble   Vector of parameters (input)
C       F(M)           Dble   Vector of equality constraints and
C                             measurement fitting values -
C                             the first MCON-components belonging
C                             to the constraints (output).
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             Whenever a nonzero value is returned
C                             by FCN routine NCJCF is terminated
C                             immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of columns of the Jacobian
C  M          Int     Number of rows of the Jacobian
C  MCON       Int     Number of equality constraints (MCON.LE.M)
C                     (only passed to FCN)
C  LDA        Int     Leading dimension of A (LDA .GE. M)
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(M)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  ETA(N)     Dble *  Array containing the scaled denominator
C                     differences
C  ETAMIN     Dble    Minimum allowed scaled denominator
C  ETAMAX     Dble    Maximum allowed scaled denominator
C  ETADIF     Dble    DSQRT (1.1*EPMACH)
C                     EPMACH = machine precision
C  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  ETA(N)     Dble *  Scaled denominator differences adjusted
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(M)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
C     EXTERNAL FCN
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT
C
C* Constants
C  ---------
C
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0,
     $           ZERO   = 0.0D0)
C
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K, IS
      DOUBLE PRECISION FHI, HG, U, SUMD, W
      LOGICAL QFINE
C
C* Begin
C
      DO 1 K = 1,N
         IS = 0
C        DO (Until)
 11         CONTINUE
            W = X(K)
            U = DSIGN (ETA(K)*YSCAL(K), X(K))
            X(K) = W + U
            CALL FCN (N, M, MCON, X, FU, IFAIL)
            NFCN = NFCN + 1
C           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
            X(K) = W
            SUMD = ZERO
            DO 111 I = 1,M
               HG = DMAX1 (DABS (FX(I)), DABS (FU(I)))
               FHI = FU(I) - FX(I)
               IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
               A(I,K) = FHI / U
 111        CONTINUE
            SUMD = DSQRT (SUMD / DBLE(M))
            QFINE = .TRUE.
            IF (SUMD .NE. ZERO .AND. IS .EQ. 0)THEN
               ETA(K) = DMIN1 (ETAMAX,
     $              DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
               IS = 1
               QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
            ENDIF
            IF (.NOT.(QFINE)) GOTO  11
C        UNTIL ( expression - negated above)
 1    CONTINUE
C
C     Exit from DO-loop
 99   CONTINUE
C
      RETURN
C
C* End of subroutine NCJCF
C
      END
      SUBROUTINE NCRNK1(N,M,XW,DX,F,FA,DXJ,DXF,A,FCA)
C*    Begin Prologue RANK1
      INTEGER N,M
      DOUBLE PRECISION FCA
      DOUBLE PRECISION A(M,N),DX(N),XW(N),F(M),FA(M),DXJ(N),
     $DXF(M)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     R A N K 1 : To be used in connection with NLSCON .
C     provides Rank-1 updates of Jacobian matrix A
C
C*    Input parameters
C     ================
C
C       N          Int    Number of columns of the Jacobian
C       M          Int    Number of rows of the Jacobian
C       XW(N)      Dble   Array containing the scaling factors
C       DX(N)      Dble   Last (unrelaxed) Gauss-Newton correction
C       F(M)       Dble   FCN(x(k)),  with x(k)= current iterate
C       FA(M)      Dble   FCN(x(k-1)),  with x(k-1)= previous
C                         iterate
C       FCA        Dble   Previous damping factor
C
C*    Output parameters:
C     ==================
C
C       A(M,N)     Dble   Array to contain the approximated
C                         Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C
C*    Workspace parameters:
C     =====================
C
C       DXJ(N)     Dble   For evaluation of broyden update -
C                         stores (deltaxk)t /( FCA
C                         * norm(deltaxk)**2 )
C       DXF(M)     Dble   For evaluation of Broyden update -
C                         stores(F-(1-FCA)*FA)
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION DNM,S1
      INTEGER K,L1
C*    Begin
        DO 10 L1=1,N
          DXJ(L1)=DX(L1)/XW(L1)
10      CONTINUE
        DNM = ZERO
        DO 20 L1=1,N
          DNM = DNM+DXJ(L1)**2
20      CONTINUE
        DO 30 L1=1,N
          DXJ(L1)=DXJ(L1)/XW(L1)
30      CONTINUE
        DNM = DNM*FCA
        IF(DNM.NE.ZERO)THEN
          S1 = FCA-ONE
          DO 41 L1=1,M
            DXF(L1)=F(L1)+FA(L1)*S1
41        CONTINUE
          DO 42 K=1,N
            S1 = DXJ(K)/DNM
            DO 421 L1=1,M
              A(L1,K) = A(L1,K)+DXF(L1)*S1
421         CONTINUE
42        CONTINUE
        ENDIF
C       End of subroutine NCRNK1
        RETURN
      END
C
      SUBROUTINE NCPRJN(N,IRANK,DEL,U,D,V,QE,PIVOT)
C*    Begin Prologue PRJCTN
      INTEGER IRANK,N
      INTEGER PIVOT(N)
      DOUBLE PRECISION DEL
      DOUBLE PRECISION QE(N,N)
      DOUBLE PRECISION U(N),D(N),V(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     P R J C T N :
C     To be used in connection with either DECOMP/SOLVE or 
C     DECCON/SOLCON .
C     Provides the projection to the appropriate subspace in case
C     of rank - reduction
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N              Int    Number of parameters to be estimated
C       IRANK                 Pseudo rank of decomposed Jacobian
C                             matrix
C       U(N)           Dble   Scaled Gauss-Newton correction
C       D(N)           Dble   Diagonal elements of upper
C                             triangular matrix
C       QE(N,N)        Dble   Part of pseudoinverse Jacobian
C                             matrix ( see QA of DECCON )
C       PIVOT(N)       Dble   Pivot vector resulting from matrix
C                             decomposition (DECCON)
C       V(N)           Dble   Real work array
C
C*    Output parameters
C     =================
C
C       DEL            Dble   Defekt
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,I,IRK1
      DOUBLE PRECISION S,SH
C*    Begin
      DO 1 I=1,N
        V(I)=U(PIVOT(I))
1     CONTINUE
      IRK1 = IRANK+1
      DEL = ZERO
      DO 2 I=IRK1,N
        SH = 0.0
        DO 21 L1=1,I-1
          SH = SH+QE(L1,I)*V(L1)
21      CONTINUE
        S =(V(I)-SH)/D(I)
        DEL = S*S+DEL
        V(I)=S
2     CONTINUE
C     End of subroutine NCPRJN
      RETURN
      END
C
      SUBROUTINE NCPRV1(DLEVF,DLEVX,FC,NITER,NEW,IRANK,MPRMON,LUMON,
     $                  QMIXIO)
C*    Begin Prologue NCPRV1
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,NEW,IRANK,MPRMON,LUMON
      LOGICAL QMIXIO
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N C P R V 1 : Printing of intermediate values (Type 1 routine)
C
C     Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of NCINT
C     FC,NITER,NEW,IRANK,MPRMON,LUMON
C                  See parameter descr. of subroutine NCINT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level
      IF(QMIXIO)THEN
1       FORMAT(2X,74('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',14X,'Normx ',18X,'New',6X,'Rank')
        IF (MPRMON.GE.3) WRITE(LUMON,2)
3       FORMAT(8X,'It',7X,'Normf ',14X,'Normx ',6X,'Damp.Fct.',3X,'New',
     $         6X,'Rank')
        IF (MPRMON.EQ.2) WRITE(LUMON,3)
      ENDIF
4     FORMAT(6X,I4,5X,D14.7,6X,D10.3,15X,I2,6X,I4)
      IF (MPRMON.GE.3.OR.NITER.EQ.0) 
     $  WRITE(LUMON,4) NITER,DLEVF,DLEVX,NEW,IRANK
5     FORMAT(6X,I4,5X,D14.7,6X,D10.3,6X,F5.3,4X,I2,6X,I4)
      IF (MPRMON.EQ.2.AND.NITER.NE.0) 
     $  WRITE(LUMON,5) NITER,DLEVF,DLEVX,FC,NEW,IRANK
      IF(QMIXIO)THEN
6       FORMAT(2X,74('*'))
        WRITE(LUMON,6)
      ENDIF
C     End of subroutine NCPRV1
      RETURN
      END
C
      SUBROUTINE NCPRV2(DLEVF,DLEVX,FC,NITER,MPRMON,LUMON,QMIXIO,
     $                  CMARK)
C*    Begin Prologue NCPRV2
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,MPRMON,LUMON
      LOGICAL QMIXIO
      CHARACTER*1 CMARK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N C P R V 2 : Printing of intermediate values (Type 2 routine)
C
C*    Parameters
C     ==========
C
C     DLEVF,DLEVX    See descr. of internal double variables of NCINT
C     FC,NITER,MPRMON,LUMON
C                  See parameter descr. of subroutine NCINT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C     CMARK Char*1    Marker character to be printed before DLEVX
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level, and damping
C     factor
      IF(QMIXIO)THEN
1       FORMAT(2X,74('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',14X,'Normx ',6X,'Damp.Fct.')
        WRITE(LUMON,2)
      ENDIF
3     FORMAT(6X,I4,5X,D14.7,4X,A1,1X,D10.3,6X,F5.3,19X,F5.3)
      WRITE(LUMON,3)NITER,DLEVF,CMARK,DLEVX,FC
      IF(QMIXIO)THEN
4       FORMAT(2X,74('*'))
        WRITE(LUMON,4)
      ENDIF
C     End of subroutine NCPRV2
      RETURN
      END
C
      SUBROUTINE NCSOUT(N,MFIT,X,RQ,MODE,IOPT,RWK,NRW,IWK,NIW,MPRINT,
     $                  LUOUT)
C*    Begin Prologue SOLOUT
      INTEGER N,MFIT
      DOUBLE PRECISION X(N),RQ(MFIT)
      INTEGER NRW
      INTEGER MODE
      INTEGER IOPT(50)
      DOUBLE PRECISION RWK(NRW)
      INTEGER NIW
      INTEGER IWK(NIW)
      INTEGER MPRINT,LUOUT
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L O U T : Printing of iterate (user customizable routine)
C
C*    Input parameters
C     ================
C
C     N           Int Number of equations/unknowns
C     X(N)     Dble   Iterate vector
C     RQ(MFIT) Dble   Linear residuum (without zero components
C                     corresponding to the equality constraints)
C     MODE            =1 This routine is called before the first
C                        Gauss-Newton iteration step
C                     =2 This routine is called with an intermedi-
C                        ate iterate X(N)
C                     =3 This is the last call with the solution
C                        vector X(N)
C                     =4 This is the last call with the final, but
C                        not solution vector X(N)
C     IOPT(50)    Int The option array as passed to the driver
C                     routine (elements 46 to 50 may be used
C                     for user options)
C     MPRINT      Int Solution print level 
C                     (see description of IOPT-field MPRINT)
C     LUOUT       Int the solution print unit 
C                     (see description of see IOPT-field LUSOL)
C
C
C*    Workspace parameters
C     ====================
C
C     NRW, RWK, NIW, IWK    see description in driver routine
C
C*    Use of IOPT by this routine
C     ===========================
C
C     Field 46:       =0 Standard output
C                     =1 GRAZIL suitable output
C
C     ------------------------------------------------------------
C*    End Prologue
      LOGICAL QGRAZ,QNORM
C*    Begin
      QNORM = IOPT(46).EQ.0
      QGRAZ = IOPT(46).EQ.1
      IF(QNORM) THEN
1        FORMAT('  ',A,' data:',/)
         IF (MODE.EQ.1) THEN
101        FORMAT('  Start data:',/,'  N =',I5,//,
     $            '  Format: iteration-number, (x(i),i=1,...N), ',
     $            'Normf , Normx ',/)
           WRITE(LUOUT,101) N
           WRITE(LUOUT,1) 'Initial'
         ELSE IF (MODE.EQ.3) THEN
           WRITE(LUOUT,1) 'Solution'
         ELSE IF (MODE.EQ.4) THEN
           WRITE(LUOUT,1) 'Final'
         ENDIF
2        FORMAT(' ',I5)
C        WRITE          NITER
         WRITE(LUOUT,2) IWK(1)
3        FORMAT((12X,3(D18.10,1X)))
         WRITE(LUOUT,3)(X(L1),L1=1,N)
C        WRITE          DLEVF,  DLEVX
         WRITE(LUOUT,3) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
         IF (MODE.EQ.2.AND.MPRINT.GE.3) THEN
301        FORMAT(/,'    Residuum for current iteration :')
           IF (MFIT.NE.0) THEN
             WRITE(LUOUT,301)
             WRITE(LUOUT,3) (RQ(L1),L1=1,MFIT)
           ENDIF
         ENDIF
         IF(MODE.EQ.1.AND.MPRINT.GE.2) THEN
           WRITE(LUOUT,1) 'Intermediate'
         ELSE IF(MODE.GE.3) THEN
           WRITE(LUOUT,1) 'End'
         ENDIF
      ENDIF
      IF(QGRAZ) THEN
        IF(MODE.EQ.1) THEN
10        FORMAT('&name com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,10)(I,I=1,N+2)
15        FORMAT('&def  com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,15)(I,I=1,N+2)
16        FORMAT(6X,': X=1, Y=',I3)
          WRITE(LUOUT,16) N+2
        ENDIF
20      FORMAT('&data ',I5)
C        WRITE          NITER
        WRITE(LUOUT,20) IWK(1) 
21      FORMAT((6X,4(D18.10)))
        WRITE(LUOUT,21)(X(L1),L1=1,N)
C        WRITE          DLEVF,  DLEVX
        WRITE(LUOUT,21) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
        IF(MODE.GE.3) THEN
30        FORMAT('&wktype 3111',/,'&atext x ''iter''')
          WRITE(LUOUT,30)
35        FORMAT('&vars = com',I3.3,/,'&atext y ''x',I3,'''',
     $           /,'&run')
          WRITE(LUOUT,35) (I,I,I=1,N)
36        FORMAT('&vars = com',I3.3,/,'&atext y ''',A,'''',
     $           /,'&run')
          WRITE(LUOUT,36) N+1,'Normf ',N+2,'Normx '
C39     FORMAT('&stop')
C       WRITE(LUOUT,39)
        ENDIF
      ENDIF
C     End of subroutine NCSOUT
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION WNORM(N,Z,XW)
      INTEGER N
      DOUBLE PRECISION Z(N), XW(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     W N O R M : Return the norm to be used in exit (termination)
C                 criteria
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     Z(N)     Dble  The vector, of which the norm is to be computed
C     XW(N)    Dble  The scaling values of Z(N)
C
C*    Output
C     ======
C
C     WNORM(N,Z,XW)  Dble  The mean square root norm of Z(N) subject
C                          to the scaling values in XW(N):
C                          = Sqrt( Sum(1,...N)((Z(I)/XW(I))**2) / N )
C
C     ------------------------------------------------------------
C*    End Prologue
      INTEGER I
      DOUBLE PRECISION S
C*    Begin
      S = 0.0D0
      DO 10 I=1,N
        S = S + ( Z(I)/XW(I) ) ** 2
10    CONTINUE
      WNORM = DSQRT( S / DBLE(FLOAT(N)) )
C     End of function WNORM
      RETURN
      END
C
      SUBROUTINE STACON (NDCL,MCODCL,MDCL,MFIT,N,MCON,Y,X,
     $                   YMODEL,A,IRANKC,IRANK,D,IPIV,AH,V,
     $                   IERR,VCV,RINV,RES,XL,XR,SIGMA2,
     $                   MPRMON,LUMON)
      INTEGER NDCL,MCODCL,MDCL,MFIT,N,MCON
      DOUBLE PRECISION Y(MFIT),X(N),YMODEL(MFIT),A(MDCL,NDCL)
      INTEGER IRANKC,IRANK
      DOUBLE PRECISION D(N)
      INTEGER IPIV(N)
      DOUBLE PRECISION AH(NDCL,N),V(N)
      INTEGER IERR
      DOUBLE PRECISION VCV(NDCL,N),RINV(NDCL,N),RES(MFIT),XL(N),XR(N)
      DOUBLE PRECISION SIGMA2
      INTEGER MPRMON,LUMON
C----------------------------------------------------------------------
C
C  Statistical analysis of constrained linear least squares estimates.
C  Computation of covariance matrix, correlation coefficients  and
C  confidence intervals for final parameter estimates.
C
C  To be used in connection with DECCON and SOLCON.
C
C----------------------------------------------------------------------
C
C  Date of latest change:  November 24, '93
C  By: U. Nowak, L. Weimann
C
C***********************************************************************
C
C   Input parameters
C   ----------------
C
C      NDCL          Int   Declared number of columns of a (see below)
C      MCODCL        Int
C      MDCL          Int   Declared number of rows of a (see below)
C      MFIT          Int   Number of measurements (observations)
C      N             Int   Number of parameters
C      MCON          Int   Number of equality constraints
C      Y(MFIT)       Dble  The vector of measurements
C      X(N)          Dble  Final (estimated) parameters
C      YMODEL(MFIT)  Dble  The model function vector evaluated for
C                          the final parameter estimates.
C      A(MDCL,NDCL)  Dble  Matrix of the model (unscaled Jacobian)
C      IRANKC        Int   Rank of the matrix constrains part
C      IRANK         Int   Rank of the matrix least squares part
C      D(N)          Dble  Diagonal elements of decomposed matrix A
C      IPIV(N)       Int   Column interchanges performed by 'DECCON'
C      AH(NDCL,NDCL) Dble  The rank deficient pseudo inverse computed
C                          by DECCON (if IRANK.LT.N)
C      V(N)          Dble  Work array
C      MPRMON        Int   The printing level:
C                          =0 : no printing will be done
C                          =1 : Information will be printed
C      LUMON         Int   The print unit number
C
C  Output parameters
C  -----------------
C
C      IERR          Int    Error indicator
C                           =0 : no error occured
C      RES(MFIT)     Dble   The residuum YMODEL-Y
C      RINV(NDCL,N)  Dble   The matrix of the correlation coefficients
C                           (lower triangle only) - only, if MPRMON.GT.0
C      VCV(N,N)      Dble   Correlation matrix
C      SIGMA2        Dble   Estimated variance of residual
C      XL(N)         Dble   The left bounds of the confidence intervals
C                           associated to the final parameter estimate
C                           vector X(N)
C      XR(N)         Dble   The right bounds of the confidence intervals
C                           associated to the final parameter estimate
C                           vector X(N)
C
C***********************************************************************
C
      DOUBLE PRECISION ZERO,ONE,HUNDRE
      PARAMETER ( ZERO=0.0D0, ONE=1.0D0, HUNDRE=100.0D0 )
      INTEGER I,IH,J,JP1,KRED,L,LH,M,NH1,NH2,NM1
      DOUBLE PRECISION H,FNMA,SUM
      DOUBLE PRECISION FISH15(47)
C
C  FISH15: Array containing upper 5% values of fisher(1,l)-distribution
C  (L=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
C  26,27,28,29,30,32,34,36,38,40,42,44,46,48,50,60,70,80,90,100,200,300)
C
      DATA FISH15/161.0D0,18.51D0,10.13D0,7.71D0,6.61D0,5.99D0,5.59D0,
     $ 5.32D0,5.12D0,4.96D0,4.84D0,4.75D0,4.67D0,4.6D0,4.54D0,4.49D0,
     $ 4.45D0,4.41D0,4.38D0,4.35D0,4.32D0,4.3D0,4.28D0,4.26D0,4.24D0,
     $ 4.23D0,4.21D0,4.2D0,4.18D0,4.17D0,4.15D0,4.13D0,4.11D0,4.1D0,
     $ 4.08D0,4.07D0,4.06D0,4.05D0,4.04D0,4.03D0,4.0D0,3.98D0,3.96D0,
     $ 3.95D0,3.94D0,3.89D0,3.88D0/
C
C  Estimate variance and standard deviation of residual
C 
      IF (MPRMON.GT.0) WRITE(LUMON,1200)
C
      SUM=ZERO
      DO 10 I=1,MFIT  
         RES(I)=YMODEL(I)-Y(I)
         SUM=SUM+RES(I)**2
10    CONTINUE
      SUM=SUM/DBLE(MFIT-(IRANK-MCON))
C
      SIGMA2=SUM
      IF (MPRMON.GT.0) THEN
        WRITE(LUMON,1000)
        WRITE(LUMON,1005) SIGMA2,DSQRT(SIGMA2)
      ENDIF
C
C Computation of covariance matrix
C---------------------------------
C
      IF(N.NE.1) THEN
C
C  Apply pseudo-inverse to cholesky decomposition of
C  variance covariance matrix of data errors. We assume, that
C  this matrix is sigma*I
C
        M=MFIT+MCON
        DO 3300 J=1,N
           DO 3010 I=1,N
             XL(I)=ZERO
3010       CONTINUE
           XL(J)=ONE
           IF(J.LE.MCON) XL(J)=ZERO
           KRED=-1
           CALL SOLCON (A,MDCL,NDCL,MCON,M,N,XR,XL,IRANKC,IRANK,D,
     1                IPIV,KRED,AH,V)
           DO 3250 I=1,N
              RINV(I,J)=XR(I)
3250       CONTINUE
3300    CONTINUE
C
C  Product RINV*RINV^T
        DO 3500 I=1,N
           DO 3500 J=1,N
              SUM=ZERO
              DO 3400 L=1,N
                 SUM=SUM+RINV(I,L)*RINV(J,L)
3400          CONTINUE
           VCV(I,J)=SIGMA2*SUM
3500    CONTINUE
      ELSE
        VCV(1,1)=SIGMA2/(D(1)*D(1))
      ENDIF
C
C  Pretty print out of the covariance matrix
C  -----------------------------------------
      IF (MPRMON.GT.0) THEN
      WRITE(LUMON,1010)
      NH2=0
250   CONTINUE
        NH1=NH2+1
        NH2=NH1+6
        IF(NH2.GT.N) NH2=N
        IF (MPRMON.GT.0) WRITE(LUMON,1011) (L,L=NH1,NH2)
        DO 260 I=NH1,N
          IH=I
          IF(IH.GT.NH2) IH=NH2
          IF (MPRMON.GT.0) WRITE(LUMON,1012) I,(VCV(I,J),J=NH1,IH)
260     CONTINUE
      IF(NH2.NE.N) GOTO 250
      ENDIF
C
C  Computation and pretty printout of the correlation coefficients
C  ---------------------------------------------------------------
      DO 300 I=1,N
        V(I)=DSQRT(VCV(I,I))
        IF (V(I).NE.ZERO) THEN
          RINV(I,I)=ONE
        ELSE
          RINV(I,I)=ZERO
        ENDIF
300   CONTINUE
C
      IF (MPRMON.GT.0) THEN
        IF(N.NE.1) THEN
          NM1=N-1
          DO 320 J=1,NM1
            JP1=J+1
            IF (V(J).EQ.ZERO) THEN
              DO 305 I=JP1,N
                RINV(I,J)=ZERO
305           CONTINUE
            ELSE
              DO 310 I=JP1,N
                IF (V(I).NE.ZERO) THEN
                  RINV(I,J)=VCV(I,J)/(V(I)*V(J))
                ELSE
                  RINV(I,J)=ZERO
                ENDIF
310           CONTINUE
            ENDIF
320       CONTINUE
C
          WRITE(LUMON,1020)
          NH2=0
350       CONTINUE
            NH1=NH2+1
            NH2=NH1+9
            IF(NH2.GT.N) NH2=N
            WRITE(LUMON,1021) (L,L=NH1,NH2)
            DO 360 I=NH1,N
              IH=I
              IF(IH.GT.NH2) IH=NH2
              WRITE(LUMON,1022) I,(RINV(I,J),J=NH1,IH)
360         CONTINUE
          IF(NH2.NE.N) GOTO 350
        ENDIF
      ENDIF
C
C  Standard error in parameters
C  ----------------------------
      IF (MPRMON.GT.0) THEN
        WRITE(LUMON,1030)
        DO 400 I=1,N
          IF (X(I).NE.ZERO) THEN
            XR(I)=DABS(HUNDRE*V(I)/X(I))
          ELSE
            XR(I)=DABS(HUNDRE*V(I))
          ENDIF
400     CONTINUE
        WRITE(LUMON,1060) (I,X(I),V(I),XR(I),I=1,N)
      ENDIF
C
      L=MFIT-(N-MCON)
      IF(L.GT.30) GOTO 410
      LH=L
      GOTO 450
410   IF(L.GT.50) GOTO 420
      LH=15+L/2
      GOTO 450
420   IF(L.GT.100) GOTO 430
      LH=35+L/10
      GOTO 450
430   IF(L.GE.300) GOTO 440
      LH=45+L/200
      GOTO 450
440   LH=47
450   CONTINUE
C
C  Associated confidence intervals
C  -------------------------------
C
      IF (MPRMON.GT.0) WRITE(LUMON,1050)
      FNMA=FISH15(LH)
      H=DSQRT(FNMA)
      IF (MPRMON.GT.0) WRITE(LUMON,1085) FNMA
      DO 460 I=1,N
        XR(I)=H*V(I)
        V(I)=HUNDRE*XR(I)/X(I)
        XL(I)=X(I)-XR(I)
        XR(I)=X(I)+XR(I)
        IF (MPRMON.GT.0) WRITE(LUMON,1062) I,XL(I),XR(I)
460   CONTINUE
C
      RETURN
C  --------------------------------------------------------------------
1200  FORMAT(///,10X,'under the assumptions of the classical linear',
     $          ' model:',/,11X,52(' '),///)
1000  FORMAT(/,3X,'Best unbiased estimate of variance and',
     $          ' standard deviation of residuals:',/,3X,71('-'))
1005  FORMAT(/,3X,'sigma2 =',D10.3,5X,'sigma =',D10.3)
1010  FORMAT(///,3X,'Covariance matrix of parameters',/,3X,17('-'))
1011  FORMAT(/,7(I10))
1012  FORMAT(1X,I3,7D10.2)
1020  FORMAT(///,3X,'Correlation coefficients',/,3X,24('-'))
1021  FORMAT(/,3X,10(I7))
1022  FORMAT(1X,I3,10F7.2)
1030  FORMAT(///,3X,'Standard deviation of parameters',/,3X,32('-'),/,
     $       5X,'No.',2X,'Estimate',11X,'sigma(X)')
1035  FORMAT(3X,I4,D14.3,D14.3)
1050  FORMAT(///,3X,'Independent confidence intervals',/,3X,32('-'))
1060  FORMAT(3X,I4,2X,D10.3,3X,'+/-',2X,D10.3,5X,'=',F8.2,' %')
1062  FORMAT(3X,I4,'  ( ',D10.3,' , ',D10.3,' )')
1075  FORMAT(3X,'rho =',D10.2,3X,'(',F6.1,' %)',//)
1085  FORMAT(3X,'   (on 95%-probability level using ',
     $       'F-distribution  F(alfa,1,m-n)=',F6.2,')',/)
C
      END
C*    End package
C
C
C*    Group  Linear Solver subroutines (Code DECCON/SOLCON)
C
      SUBROUTINE DECCON(A,NROW,NCOL,MCON,M,N,IRANKC,IRANK,COND,D,PIVOT,
     *KRED,AH,V,IERR)
C*    Begin Prologue DECCON
      INTEGER IRANKC,IRANK,MCON
      INTEGER M,N,NROW,NCOL,KRED
      INTEGER PIVOT(NCOL)
      DOUBLE PRECISION COND
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION D(NCOL),V(NCOL)
      INTEGER IERR
C     ------------------------------------------------------------
C
C*  Title
C
C*    Deccon - Constrained Least Squares QR-Decomposition
C
C*  Written by        P. Deuflhard, U. Nowak, L. Weimann 
C*  Purpose           Solution of least squares problems, optionally
C                     with equality constraints.
C*  Method            Constrained Least Squares QR-Decomposition
C                     (see references below)
C*  Category          D9b1. -  Singular, overdetermined or
C                              underdetermined systems of linear 
C                              equations, generalized inverses. 
C                              Constrained Least Squares solution
C*  Keywords          Linear Least Square Problems, constrained, 
C                     QR-decomposition, pseudo inverse.
C*  Version           1.2
C*  Revision          December 1993
C*  Latest Change     December 1993
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin
C                     Heilbronner Str. 10
C                     D-10711 Berlin-Wilmersdorf
C                     phone 0049+30+89604-0, 
C                     telefax 0049+30+89604-125
C*  Contact           Lutz Weimann 
C                     ZIB, Numerical Software Development 
C                     phone: 0049+30+89604-185 ;
C                     e-mail:
C                     RFC822 notation: weimann@sc.zib-berlin.de
C                     X.400: C=de;A=dbp;P=zib-berlin;OU=sc;S=Weimann
C
C*    References:
C     ===========
C
C       /1/ P.Deuflhard, V.Apostolescu:
C           An underrelaxed Gauss-Newton method for equality
C           constrained nonlinear least squares problems.
C           Lecture Notes Control Inform. Sci. vol. 7, p.
C           22-32 (1978)
C       /2/ P.Deuflhard, W.Sautter:
C           On rank-deficient pseudoinverses.
C           J. Lin. Alg. Appl. vol. 29, p. 91-111 (1980)
C    
C*    Related Programs:     SOLCON
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from acquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Constrained QR-decomposition of (M,N)-system  with
C     computation of pseudoinverse in case of rank-defeciency .
C     First MCON rows belong to equality constraints.
C
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    Input parameters
C     ================
C
C       A(NROW,NCOL) Dble   Array holding the (M,N)-Matrix to be 
C                           decomposed
C       NROW         Int    Declared number of rows of array A
C       NCOL         Int    Declared number of columns of array A and 
C                           rows and columns of array AH
C       MCON         Int    Number of equality constraints (MCON.LE.N)
C                           Internally reduced if equality constraints
C                           are linearly dependent
C       M            Int    Current number of rows of matrix A
C       N            Int    Current number of columns of matrix A
C     * IRANKC       Int    Prescribed maximum pseudo-rank of 
C                           constrained part of matrix A (IRANKC.LE.MCON)
C     * IRANK        Int    Prescribed maximum pseudo-rank of matrix A
C                           (IRANK.LE.N)
C     * COND         Dble   Permitted upper bound for the subcondition
C                           of the least squares part of A, .i.e.
C                           DABS(D(IRANKC+1)/D(IRANK))
C       KRED         Int    Type of operation
C                           >=0  Householder triangularization
C                                (build up pseudo-inverse,if IRANK.LT.N)
C                           < 0  Reduction of pseudo-rank of matrix A, 
C                                skipping Householder triangularization,
C                                 build-up new pseudo-inverse
C
C*    Output parameters
C     =================
C
C       A(NROW,NCOL)  Dble   Array holding the (M,N)-output consisting
C                            of the transformed matrix in the upper 
C                            right triangle and the performed House-
C                            holder transf. in the lower left triangle.
C     * IRANKC        Int    New pseudo-rank of constrained part of
C                            matrix A, determined so that
C                            DABS(D(1)/D(IRANKC))<1/EPMACH
C     * IRANK         Int    New pseudo-rank of matrix A, determined
C                            so that DABS(D(IRANKC+1)/D(IRANK)) < COND
C       D(IRANK)      Dble   Diagonal elements of upper triangular matr.
C       PIVOT(N)      Int    Index vector storing permutation of columns
C                            due to pivoting
C     * COND          Dble   The sub-condition number belonging to the
C                            least squares part of A.
C                            (in case of rank reduction:
C                             sub-condition number which led to
C                             rank reduction)
C                            COND=0 indicates COND=infinity
C       AH(NCOL,NCOL) Dble   In case of rank-defect used to compute the
C                            pseudo-inverse (currently used will be an
C                            (N,N)-part of this array)
C       V(N)          Dble   V(1) holds on output the sub-condition
C                            number belonging to the constrained part
C                            of A.
C       IERR          Int    Error indicator:
C                            = 0 : DECCON computations are successfull.
C                            =-2 : Numerically negative diagonal element
C                                  encountered during computation of
C                                  pseudo inverse - due to extremely bad
C                                  conditioned Matrix A. DECCON is
C                                  unable to continue rank-reduction.
C
C*    Workspace parameters
C     ====================
C
C       V(N)         Dble   Workspace array
C
C*    Subroutines called: D1MACH
C
C*    Machine constants used
C     ======================
C
C     EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL D1MACH
      INTRINSIC DABS,DSQRT
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION REDUCE
      PARAMETER (REDUCE=0.05D0)
      INTEGER L1
      DOUBLE PRECISION S1
      INTEGER I,II,IRANKH,IRK1,I1,J,JD,JJ,K,K1,LEVEL,MH
      DOUBLE PRECISION DD,D1MACH,H,HMAX,S,SH,SMALLS,T,COND1
C*    Begin
C     --------------------------------------------------------------
C     1 Initialization
      EPMACH = D1MACH(3)
      COND1 = ONE/(EPMACH*10.0D0)
	T=ZERO
      IF(IRANK.GT.N) IRANK = N
      IF(IRANK.GT.M) IRANK = M
C     --------------------------------------------------------------
C     1.1 Special case M=1 and N=1
      IF(M.EQ.1.AND.N.EQ.1)THEN
        PIVOT(1)=1
        D(1)=A(1,1)
        COND = ONE
        RETURN
      ENDIF
      IF(KRED.GE.0)THEN
C       ------------------------------------------------------------
C       1.1 Initialize pivot-array
        DO 11 J=1,N
          PIVOT(J)=J
11      CONTINUE
C       ------------------------------------------------------------
C       2. Constrained Householder triangularization    
        JD = 1
        IRANC1 = IRANKC + 1
        MH = MCON
        IRANKH = IRANKC
        CONDH = COND1
        IF(MH.EQ.0) THEN
          IRANKH = IRANK
          MH = M
          CONDH = COND
        ENDIF
        IRK1 = IRANK
        DO 2  K=1,IRK1
2000      LEVEL = 1
          IF(K.NE.N)THEN
            K1 = K+1
C           DO (Until)
20          CONTINUE
              IF(JD.NE.0)THEN
                DO 201 J=K,N
                  S = ZERO
                  DO 2011 L1=K,MH
                    S = S+A(L1,J)**2
2011              CONTINUE
                  D(J)=S
201             CONTINUE
              ENDIF
C             ------------------------------------------------------
C             2.1 Column pivoting
              S1 = D(K)
              JJ = K
              DO 21 L1=K,N
                IF(D(L1).GT.S1) THEN
                  S1=D(L1)
                  JJ = L1
                ENDIF
21            CONTINUE
              H = D(JJ)
              IF(JD.EQ.1) HMAX = H/(COND1*0.1D0)
              JD = 0
              IF(H.LT.HMAX) JD = 1
            IF(.NOT.(H.GE.HMAX)) GOTO 20
C           UNTIL ( expression - negated above)
            IF(JJ.NE.K)THEN
C             ------------------------------------------------------
C             2.2 Column interchange
              I = PIVOT(K)
              PIVOT(K)=PIVOT(JJ)
              PIVOT(JJ)=I
              D(JJ)=D(K)
              DO 221 L1=1,M
                S1=A(L1,JJ)
                A(L1,JJ)=A(L1,K)
                A(L1,K)=S1
221           CONTINUE
            ENDIF
          ENDIF
          H = ZERO
          DO 222 L1=K,MH
            H = H+A(L1,K)**2
222       CONTINUE
          T = DSQRT(H)
C         ----------------------------------------------------------
C         2.3.0 A-priori test on pseudo-rank
          IF  ( K.EQ.1 .OR. K.EQ.IRANC1 )  DD = T/CONDH
          IF(T.LE.DD .OR. K.GT.IRANKH)THEN
C           ------------------------------------------------------
C           2.3.1 Rank reduction
            IRANKH = K-1
            IF  (MH.NE.MCON)  THEN
              IRANK = IRANKH
              IF (IRANKC.EQ.IRANK) THEN
                LEVEL = 4
              ELSE
                LEVEL = 3
              ENDIF
            ELSE
              IRANKC = IRANKH
              IF  (IRANKC.NE.MCON) THEN
                MH = M
                IRANKH = IRANK
                JD = 1
                CONDH = COND
                GOTO 2000
              ENDIF
            ENDIF
          ENDIF
          IF (LEVEL.EQ.1) THEN
C           ------------------------------------------------------
C           2.4 Householder step
            S = A(K,K)
            T = -DSIGN(T,S)
            D(K)=T
            A(K,K)=S-T
            IF(K.NE.N)THEN
              T = ONE/(H-S*T)
              DO 24 J=K1,N
                S = ZERO
                DO 241 L1=K,MH
                  S = S+A(L1,K)*A(L1,J)
241             CONTINUE
                S = S*T
                S1 =-S
                DO 242 L1=K,M
                  A(L1,J) = A(L1,J)+A(L1,K)*S1
242             CONTINUE
                D(J)=D(J)-A(K,J)**2
24            CONTINUE
              IF(K.EQ.IRANKC)THEN
                MH = M
                JD = 1
                CONDH = COND
                IRANKH=IRANK
                IF (K.EQ.IRK1) LEVEL=3
              ENDIF
            ELSE
              LEVEL = 4
            ENDIF
          ENDIF
C       Exit Do 2 If ... 
          IF(LEVEL.GT.1) GOTO  2999
2       CONTINUE
C       ENDDO
2999    CONTINUE
      ELSE
        K = -1
        LEVEL = 3
      ENDIF
C     --------------------------------------------------------------
C     3 Rank-deficient pseudo-inverse
      IF(LEVEL.EQ.3)THEN
        IRK1 = IRANK+1
        DO 3 J=IRK1,N
          DO 31 II=1,IRANK
            I = IRK1-II
            S = A(I,J)
            IF(II.NE.1)THEN
              SH = ZERO
              DO 311 L1=I1,IRANK
                SH=SH+A(I,L1)*V(L1)
311           CONTINUE
              S = S-SH
            ENDIF
            I1 = I
            V(I)=S/D(I)
            AH(I,J)=V(I)
31        CONTINUE
          DO 32 I=IRK1,J
            S = ZERO
            DO 321 L1=1,I-1
              S = S+AH(L1,I)*V(L1)
321         CONTINUE
            IF(I.NE.J)THEN
              V(I)=-S/D(I)
              AH(I,J)=-V(I)
            ENDIF
32        CONTINUE
          IF (S.GT.-ONE) THEN
            D(J)=DSQRT(S+ONE)
          ELSE 
            IERR=-2
            GOTO 999
          ENDIF
3       CONTINUE
      ENDIF
C    --------------------------------------------------------------
C     9 Exit
C     9.1 Subcondition of constrained part of A 
      IF (IRANKC.NE.0) THEN
        V(1) = D(IRANKC)
        IF(V(1).NE.ZERO) V(1) = DABS(D(1)/V(1))
      ELSE
        V(1)=ZERO
      ENDIF
C     9.1 Subcondition of least squares part of A
      IF (K.EQ.IRANK) T = D(IRANK)
      IF (IRANKC+1.LE.IRANK .AND. T.NE.ZERO) THEN
        COND = DABS(D(IRANKC+1)/T)
      ELSE
        COND = ZERO
      ENDIF
      IERR=0
999   RETURN
      END
      SUBROUTINE SOLCON(A,NROW,NCOL,MCON,M,N,X,B,IRANKC,IRANK,D,PIVOT,
     *KRED,AH,V)
C*    Begin Prologue SOLCON
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION X(NCOL),B(NROW),D(NCOL),V(NCOL)
      INTEGER NROW,NCOL,MCON,M,N,IRANKC,IRANK,KRED
      INTEGER PIVOT(NCOL)
C     ------------------------------------------------------------
C    
C*    Summary
C     =======
C
C     Best constrained linear least squares solution of (M,N)-
C     system . First MCON rows comprise MCON equality constraints.
C     To be used in connection with subroutine DECCON
C     References:       See DECCON
C     Related Programs: DECCON
C    
C*    Parameters:
C     ===========
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       A(M,N), NROW, NCOL, M, N, MCON, IRANKC, IRANK,
C       D(N), PIVOT(N), AH(N,N), KRED
C                           See input- respective output-parameters
C                           description of subroutine DECCON
C     * B(M)         Dble   Right-hand side of linear system, if
C                           KRED.GE.0
C                           Right-hand side of upper linear system,
C                           if KRED.LT.0
C
C*    Output parameters
C     =================
C
C       X(N)         Dble   Best LSQ-solution of linear system
C       B(M)         Dble   Right-hand of upper trigular system
C                           (transformed right-hand side of linear
C                            system)
C
C*    Workspace parameters
C     ====================
C
C       V(N)         Dble   Workspace array
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,L2
      INTEGER I,II,I1,IRANC1,IRK1,J,JJ,J1,MH
      DOUBLE PRECISION S,SH
C*    Begin
C     ------------------------------------------------------------
C     1 Solution for pseudo-rank zero
      IF(IRANK.EQ.0)THEN
        DO 11 L1=1,N
          X(L1)=ZERO
11      CONTINUE
        RETURN
      ENDIF
      IF  ((IRANK.LE.IRANKC).AND.(IRANK.NE.N)) THEN
        IRANC1 = IRANKC + 1
        DO 12 L1=IRANC1,N
          V(L1) = ZERO
12      CONTINUE
      ENDIF
      IF(KRED.GE.0.AND.(M.NE.1.OR.N.NE.1))THEN
C       ----------------------------------------------------------
C       2 Constrained householder transformations of right-hand side
        MH = MCON
        IF(IRANKC.EQ.0) MH = M
        DO 21 J=1,IRANK
          S = ZERO
          DO 211 L1=J,MH
            S = S+A(L1,J)*B(L1)
211       CONTINUE
          S = S/(D(J)*A(J,J))
          DO 212 L1=J,M
            B(L1)=B(L1)+A(L1,J)*S
212       CONTINUE
          IF(J.EQ.IRANKC) MH = M
21      CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     3 Solution of upper triangular system
      IRK1 = IRANK+1
      DO 31 II=1,IRANK
        I = IRK1-II
        I1 = I+1
        S = B(I)
        IF(II.NE.1)THEN 
          SH = ZERO
          DO 311 L1=I1,IRANK 
            SH=SH+A(I,L1)*V(L1)
311       CONTINUE
          S = S-SH
        ENDIF
        V(I)=S/D(I)
31    CONTINUE
      IF((IRANK.NE.N).AND.(IRANK.NE.IRANKC))THEN
C       ----------------------------------------------------------
C       3.2 Computation of the best constrained least squares-
C           solution
        DO 321 J=IRK1,N
          S = ZERO
          DO 3211 L1=1,J-1
            S = S+AH(L1,J)*V(L1)
3211      CONTINUE
          V(J)=-S/D(J)
321     CONTINUE
        DO 322 JJ=1,N
          J = N-JJ+1
          S = ZERO
          IF(JJ.NE.1)THEN
            DO 3221 L1=J1,N
              S=S+AH(J,L1)*V(L1)
3221        CONTINUE
          ENDIF
          IF(JJ.NE.1.AND.J.LE.IRANK)THEN
            V(J)=V(J)-S
          ELSE
            J1 = J
            V(J)=-(S+V(J))/D(J)
          ENDIF
322     CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     4 Back-permutation of solution components
      DO 4 L1=1,N
        L2 = PIVOT(L1)
        X(L2) = V(L1)
4     CONTINUE
      RETURN
      END
C*    End package
C
C
C*    Group  Time monitor package
C
C*    Begin Prologue
C     ------------------------------------------------------------
C
C*  Title
C    
C     Monitor - A package for making multiple time measurements and
C               summary statistics
C
C*  Written by        U. Nowak, L. Weimann 
C*  Version           1.0
C*  Revision          January 1991
C*  Latest Change     January 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C*  Copyright     (c) Konrad Zuse Zentrum fuer
C                     Informationstechnik Berlin
C                     Heilbronner Str. 10, D-1000 Berlin 31
C                     phone 0049+30+89604-0, 
C                     telefax 0049+30+89604-125
C*  Contact           Lutz Weimann 
C                     ZIB, Numerical Software Development 
C                     phone: 0049+30+89604-185 ;
C                     e-mail: 
C                     RFC822 notation: weimann@sc.zib-berlin.de
C                     X.400: C=de;A=dbp;P=zib-berlin;OU=sc;S=Weimann
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from aquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C  ---------------------------------------------------------------
C
C*    Summary:
C
C     Monitor is a package for generating time and summary statistics
C     about the execution of multiple program parts of any program.
C     Nested measurements of program parts are possible.
C     ------------------------------------------------------------
C
C*    Usage:
C
C     The usage of Monitor is naturally divided into three phases:
C     1. the initialization and setup phase before the start of
C        the program or subroutines package to be measured;
C     2. the run phase of the program to be measured;
C     3. the final evaluation call.
C
C     The phase 1 must start with exactly one call of the subroutine
C     MONINI, which passes a title string and a logical unit for
C     later statistics output and possible error messages to the
C     package. This call follows a number of calls of the subroutine
C     MONDEF, where each call associates an identification string
C     to a positive integer number, called the measurement index
C     - up to maxtab, where maxtab is a package constant. Multiple
C     measurement indices may be used for measurements of multiple
C     program parts. The index 0 must also be associated with some
C     identification string, and corresponds to all parts of the
C     measured program from the measurement start call till the final
C     evaluation call, which are not associated with specific positive
C     measurement indices. After all necessary MONDEF calls are done,
C     the measurements are started at begin of the program to be
C     measured by a parameterless call of MONSRT.
C     In phase 2, each program part to be measured must be immediately
C     preceeded by a call of the subroutine MONON with the associated 
C     measurement index, and must be immediately followed by a call of
C     the subroutine MONOFF with the same measurement index. Measure-
C     ments of nested program parts are possible, and nesting is allowed
C     up to the number mnest, where mnest is a package constant.
C     Calling MONOFF without a preceeding MONON call with the same 
C     measurement index, or calling one of these subroutines with a
C     measurement index not previously defined by a MONDEF call causes
C     an error stop of the program. 
C     Finally at the end of the program to be measured, the parameter-
C     less call of the subroutine MONEND closes all measurements and
C     prints the summary statistics.
C     As delivered, maxtab has a value 20 and mnest a value 10, but
C     both constants may be increased, if needed, to any possible
C     integer value, by simply changing it's values in the first 
C     parameter statement of the subroutine MONTOR below.
C
C*    Subroutines and their parameters:
C     =================================
C
C     MONINI(CIDENT,LUMON)  : Initialize Monitor
C       CIDENT  char*20  Identification string for the total measurement
C                        ( printed in summary )
C       LUMON   int      The logical unit for printing out the summary
C
C     MONDEF(MESIND,CIDMES) : Define one measurement index
C       MESIND  int      >=1 : measurement index for a specific part
C                        = 0 : measurement index for all remaining parts
C                              (i.e. not belonging to parts with 
C                               index >=1)
C       CIDMES  char*15  Identification string for the part associated
C                        with MESIND ( printed in summary )
C
C     MONSRT                : Start measurements
C       (no parameters)
C
C     MONON(MESIND)         : Start measurement of a specific part
C       MESIND  int      >=1 : measurement index for a specific part
C
C     MONOFF(MESIND)        : Stop measurement of a specific part
C       MESIND  int      >=1 : measurement index for a specific part
C
C     MONEND                : Finish measurements and print summary
C       (no parameters)
C
C
C*    Example:
C       Calling sequence:
C
C       CALL MONINI (' Example',6)
C       CALL MONDEF (0,'Solver')
C       CALL MONDEF (1,'User function')
C       CALL MONDEF (2,'User matrix')
C       CALL MONSRT ()
C       ...
C       program to be measured (part without specific measurement index)
C       ...
C 1     CONTINUE      
C       ...
C       CALL MONON (2)
C       ...  user matrix code ...
C       CALL MONOFF(2)
C       ...
C       program to be measured (part without specific measurement index)
C       ...
C       CALL MONON (1)
C       ...  user function code ...
C       CALL MONOFF(1)
C       ...
C       program to be measured (part without specific measurement index)
C       ...
C       IF (no termination) GOTO 1
C       ...
C       CALL MONEND ()
C     ------------------------------------------------------------
C 
      SUBROUTINE MONTOR
      PARAMETER(MAXTAB=20,MNEST=10)
      CHARACTER*15 NAME(MAXTAB),NAME0
      CHARACTER*20 TEXT 
      CHARACTER*(*) TEXTH 
      CHARACTER*(*) NAMEH   
      REAL SEC(MAXTAB),ASEC(MAXTAB),PC1(MAXTAB),PC2(MAXTAB)
      INTEGER COUNT(MAXTAB),INDACT(MNEST)
      LOGICAL QON(MAXTAB)
      INTEGER IOUNIT
C
      SAVE SEC,COUNT,ASEC,PC1,PC2,INDXO,TIME1,TIME0,MAXIND,NAME
      SAVE SEC0,NAME0,TEXT,MONI,QON,IONCNT,INDACT
C
C
      DATA MONI/6/ , INFO/1/ , IGRAPH/1/
C
      RETURN
C
C     initialize monitor
C
      ENTRY MONINI (TEXTH,IOUNIT)
C
      MONI=IOUNIT
      MAXIND=0
      TEXT=TEXTH
      DO 100 I=1,MAXTAB
        SEC(I)=0.
        ASEC(I)=0.
        COUNT(I)=0
        QON(I)=.FALSE.
100   CONTINUE
      DO 105 I=1,MNEST
        INDACT(I)=0
105   CONTINUE
C
      SEC0=0.
      IONCNT=0
      RETURN
C
C     define one monitor entry
C
      ENTRY MONDEF(INDX,NAMEH)
      IF(INDX.LT.0 .OR. INDX.GT.MAXTAB) GOTO 1190
      IF (INDX.GT.MAXIND) MAXIND=INDX
      IF (INDX.GT.0) THEN
        IF (COUNT(INDX).GT.0) GOTO 1290
      ENDIF
      IF (INDX.EQ.0) THEN
        NAME0 = NAMEH
      ELSE
        NAME(INDX) = NAMEH
      ENDIF
      RETURN
C
C     start monitor measurements
C 
      ENTRY MONSRT()
      CALL SECOND (TIME1)
C
C      if(igraph.gt.0) call gmini(maxind,name0,name)
C
      RETURN
C
C     start one measurement
C
      ENTRY MONON (INDX)
      IF(INDX.GT.MAXIND.OR.INDX.LE.0) GOTO 1010
      IF (QON(INDX)) GOTO 1030
      CALL SECOND(ASEC(INDX))
      QON(INDX)=.TRUE.
      IF (IONCNT.EQ.0) THEN
        SEC0=SEC0+ASEC(INDX)-TIME1
      ELSE
        INDXO=INDACT(IONCNT)
        SEC(INDXO)=SEC(INDXO)+ASEC(INDX)-ASEC(INDXO)
      ENDIF
      IONCNT=IONCNT+1
      INDACT(IONCNT)=INDX
      IF(INFO.GT.1) WRITE(MONI,*) ' enter',NAME(INDX),ASEC(INDX)
C
C      if(igraph.gt.0) call gmon(indx,sec0)
C
      RETURN
C
C     stop one measurement
C
      ENTRY MONOFF (INDX)
      IF(INDX.GT.MAXIND.OR.INDX.LE.0) GOTO 1010
      IF (.NOT. QON(INDX)) GOTO 1040
      CALL SECOND(TIME2)
      QON(INDX)=.FALSE.
      SEC(INDX)=SEC(INDX)+TIME2-ASEC(INDX)
      COUNT(INDX)=COUNT(INDX)+1
      IONCNT=IONCNT-1
      IF (IONCNT.EQ.0) THEN
        TIME1=TIME2
      ELSE
        ASEC(INDACT(IONCNT))=TIME2
      ENDIF
      IF(INFO.GT.1) WRITE(MONI,*) ' exit ',NAME(INDX),TIME2
C
C      if(igraph.gt.0) call gmoff(indx,sec(indx))
C
      RETURN
C
C     terminate monitor and print statistics
C
      ENTRY MONEND
      CALL SECOND (TIME0)
      SEC0=SEC0+TIME0-TIME1
C
      SUM=1.E-10
      DO 200 I=1,MAXIND
      SUM=SUM+SEC(I)
      IF(COUNT(I).LE.0) GOTO 200
      ASEC(I)=SEC(I)/FLOAT(COUNT(I))
200   CONTINUE
      SUM0=SUM+SEC0
C
      DO 250 I=1,MAXIND
      PC1(I)=100.*SEC(I)/SUM0
      PC2(I)=100.*SEC(I)/SUM
250   CONTINUE
      PC10=100.*SEC0/SUM0
      PC20=100.*SEC0/SUM
C
      WRITE(MONI,9500)
      WRITE(MONI,9510)
      WRITE(MONI,9505)
9500  FORMAT(///)
9510  FORMAT(1X,75('#'))
9505  FORMAT(' #',73X,'#')
      WRITE(MONI,9505)
      WRITE(MONI,9512) TEXT
9512  FORMAT(' #   Results from time monitor program for: ',A29,2X,'#')
      WRITE(MONI,9505)
      WRITE(MONI,9514) SUM0,SUM
9514  FORMAT(' #   Total time:',F11.3,5X,'Sum of parts:',F11.3,19X,'#')
      WRITE(MONI,9505)
      WRITE(MONI,9520)
9520  FORMAT(' #   ',2X,'name',12X,'calls',7X,'time',4X,'av-time',
     1       4X,'% total',6X,'% sum   #')
C
      I0=1
      WRITE(MONI,9550) NAME0,I0,SEC0,SEC0,PC10,PC20
9550  FORMAT(' #   ',A15,I8,F11.3,F11.4,F11.2,F11.2,'   #')
C
      DO 300 I=1,MAXIND
      WRITE(MONI,9550) NAME(I),COUNT(I),SEC(I),ASEC(I),PC1(I),PC2(I)
300   CONTINUE
C
C
      WRITE(MONI,9505)
      WRITE(MONI,9510)
      WRITE(MONI,9500)
C
C
C      IF(IGRAPH.GT.0) CALL GMEND
C
      RETURN
C
C  error exits
C
1010  CONTINUE
      WRITE(MONI,9010) INDX
9010  FORMAT(/,' error in subroutine monon or monoff',/,
     $         '   indx out of range    indx=',I4)
      GOTO 1111
C
1020  CONTINUE
      WRITE(MONI,9020) INDX
9020  FORMAT(/,' error in subroutine monoff',/,'   indx out of range',/,
     1         '   indx=',I4)
      GOTO 1111
C
1030  CONTINUE
      WRITE(MONI,9030) INDX
9030  FORMAT(/,' error in subroutine monon',/,
     $         '   measurement is already running for ',
     1         '   indx=',I4)
      GOTO 1111
C
1040  CONTINUE
      WRITE(MONI,9040) INDX
9040  FORMAT(/,' error in subroutine monoff',/,
     $         '   measurement has never been activated for ',
     1         '   indx=',I4)
      GOTO 1111
C
1190  CONTINUE
      WRITE(MONI,9190) MAXTAB,INDX
9190  FORMAT(/,' error in subroutine mondef',/,'   indx gt ',I4,/,
     1         '   indx=',I4)
      GOTO 1111
C
1290  CONTINUE
      WRITE(MONI,9290) INDX
9290  FORMAT(/,' error in subroutine mondef',/,'   indx = ',I4,
     1         '   already in use' )
      GOTO 1111
C
1111  STOP
C
C  end subroutine monitor
C
      END
C
C*    Group  Machine dependent subroutines and functions
C
      SUBROUTINE SECOND(RTIME)
C
C Implement a dummy SECOND() routine
C
      RTIME=RTIME+0.001
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
C     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
C     SIGNIFICANT BYTE IS STORED FIRST.
C
      DATA SMALL(1),SMALL(2) /          0,    1048576 /
      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /
C
C
      IF (I .LT. 1  .OR.  I .GT. 6) GOTO 999
      IF (I .LE. 5 ) THEN
        D1MACH = DMACH(I)
      ELSE IF (I .EQ. 6) THEN
C       D1MACH = DSQRT(DMACH(1)/DMACH(3))
        D1MACH = 4.94D-32
      ENDIF
      RETURN
  999 WRITE(6,1999) I
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
