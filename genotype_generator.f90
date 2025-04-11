!==================================
  MODULE metatrom_v7_param
! common parameters and variables
!==================================
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: sp=SELECTED_REAL_KIND(7)   ! precision up to 7 decimal digits
  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(14)   ! precision up to 14 decimal digits
  INTEGER, PARAMETER :: saved_reps=100   ! number of replicates whose variables are to be stored per metapopulation
  INTEGER, PARAMETER :: tot_HWE_rep=1000    ! number of random generations of genotypic allocations for HW proportions

  INTEGER, DIMENSION(2) :: iseed    ! stores 2 random seeds
  INTEGER :: idum
  REAL :: distribution_shape	! shape parameter for the gamma distribution

  TYPE candidate
    INTEGER :: pop          ! population to which the candidate belongs
  END TYPE candidate

  TYPE neu_locus
    INTEGER :: nalle        ! number of distinct alleles
    INTEGER :: chrom        ! chromosome to which it belongs
    INTEGER :: posit        ! locus position in the chromosome
    REAL(KIND=dp) :: distr  ! parameter for shaping distribution of allele frequencies
    REAL(KIND=dp) :: locFst ! target Fst value for SA
  END TYPE neu_locus

  TYPE sel_locus
    INTEGER :: nalle        ! number of distinct alleles
    INTEGER :: chrom        ! chromosome to which it belongs
    INTEGER :: funct        ! trait number
    INTEGER :: posit        ! locus position in the chromosome
  END TYPE sel_locus

  TYPE expected_additive_value
	REAL(KIND=dp) :: mean						! expected mean for additive value
	REAL(KIND=dp) :: mean_threshold	! variation of mean additive value
	REAL(KIND=dp) :: sd							! expected standard error for additive value
  END TYPE expected_additive_value

  TYPE(candidate), ALLOCATABLE, DIMENSION(:) :: indv        ! stores pedigree
  TYPE(neu_locus), ALLOCATABLE, DIMENSION(:) :: neu_gene    ! stores gene information
  TYPE(sel_locus), ALLOCATABLE, DIMENSION(:) :: sel_gene    ! stores gene information
  TYPE(expected_additive_value), ALLOCATABLE, DIMENSION(:) :: target_additive_value    ! stores additive value information

  INTEGER :: done_replicate
  INTEGER :: counter
  INTEGER :: tot_neu_loci       ! number of neutral loci in the system
  INTEGER :: samp_neu_loci      ! number of neutral loci to be sampled
  INTEGER :: tot_sel_loci       ! number of selected loci in the system
  INTEGER :: num_traits         ! number of selective traits in the system
  INTEGER :: num_pop_tot        ! number of populations in the original metapopulation
  INTEGER :: num_tot_ind        ! total number of individuals across populations
  INTEGER :: num_tot_haplot     ! total number of haplotypes: 2*num_tot_ind
  INTEGER :: max_neu_allele     ! maximum number of neutral alleles across loci
  INTEGER :: max_sel_allele     ! maximum number of selected alleles across loci
  INTEGER :: max_loci_trait     ! maximum number of loci for a trait
  INTEGER :: max_ind_pop        ! maximum number of individuals in a population
  INTEGER :: jlast 				! rank code for last individual
  INTEGER :: num_meta           ! number of generated metapopulations
  INTEGER :: current_meta       ! current metapopulation
  INTEGER :: targetFst_mode     ! Global Fst (0) or Fst especified by locus (1)
  INTEGER :: input_line					! to known line of reading problem
	INTEGER :: progress_counter		! to track progress of the while loops
	character*1 :: tab = char(9)

  INTEGER, ALLOCATABLE, DIMENSION(:) :: trait_sel_loci

  INTEGER, ALLOCATABLE, DIMENSION(:) :: num_sel_loci        ! number of selective loci in the system
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ind_per_pop         ! number of individuals per i population

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: posit_sel_loci    ! genome position per selective locus
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: posit_ind_pop     ! position in pedigree according to each pop
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neu_copies_dist   ! distribution of neutral allelic copies
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sel_copies_dist   ! distribution of selective allelic copies

  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: neu_genome  ! neutral genome per haplotype, locus, individual
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: sel_genome  ! selected genome per haplotype, locus, individual

  REAL(KIND=dp) :: xrand            ! variable for random numbers
  REAL(KIND=dp) :: Ht_Meta          ! average Ht at neutral loci in the metapopulation as if parametric
  REAL(KIND=dp) :: S_Ht_Meta        ! sum average Ht's at neutral loci across metapopulations
  REAL(KIND=dp) :: Sq_Ht_Meta       ! sum of squares of average Ht's at neutral loci across metapopulations
  REAL(KIND=dp) :: Hs_Meta          ! observed Hs at neutral loci in the metapopulation as if parametric
  REAL(KIND=dp) :: S_Hs_Meta        ! sum of observed Hs's across metapopulations
  REAL(KIND=dp) :: Sq_Hs_Meta       ! sum of squares of observed Hs's across metapopulations
  REAL(KIND=dp) :: exp_Fst          ! expected Fst: target of the SA algorithm
  REAL(KIND=dp) :: ini_Fst_Meta     ! current initial Fst in the metapopulation
  REAL(KIND=dp) :: ini_S_Fst_Meta   ! sum of initial Fst's across metapopulations
  REAL(KIND=dp) :: ini_Sq_Fst_Meta  ! sum of squares of initial Fst's across metapopulations
  REAL(KIND=dp) :: obs_Fst_Meta     ! current final Fst in the metapopulation
  REAL(KIND=dp) :: obs_S_Fst_Meta   ! sum of final Fst's across metapopulations
  REAL(KIND=dp) :: obs_Sq_Fst_Meta  ! sum of squares of final Fst's across metapopulations


  REAL(KIND=dp), DIMENSION(2) :: xdum

  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: loci_contrib_VA   ! contributions of selective loci to additive variance

  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_locFst_Meta       ! current Fst per locus in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_S_locFst_Meta     ! sum of Fst's per locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_Sq_locFst_Meta    ! sum of squares of Fst's per locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_locQst_Meta       ! current Qst per locus in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_S_locQst_Meta     ! sum of Qst's per locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_Sq_locQst_Meta    ! sum of squares of Qst's per locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_locHt_Meta        ! current Het per neutral locus in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_S_locHt_Meta      ! sum of Het's per neutral locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_Sq_locHt_Meta     ! sum of squares of Het's per neutral locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: locHtsel_Meta         ! current Ht per selected locus in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_locHtsel_Meta       ! sum of Ht's per selected locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_locHtsel_Meta      ! sum of squares of Ht's per selected locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: locHssel_Meta         ! current Hs per selected locus in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_locHssel_Meta       ! sum of Hs's per selected locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_locHssel_Meta      ! sum of squares of Hs's per selected locus across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: exp_Qst               ! expected Qst: target of the SA algorithm
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ini_Qst_Meta          ! current initial Qst in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_ini_Qst_Meta        ! sum of initial Qst's across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_ini_Qst_Meta       ! sum of squares of initial Qst's across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: obs_Qst_Meta          ! current final Qst in the metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_obs_Qst_Meta        ! sum of final Qst's across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_obs_Qst_Meta       ! sum of squares of initial Qst's across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VGenot_tot            ! genotypic variance per trait
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_VGenot_tot          ! sum of genotypic variances per trait across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_VGenot_tot         ! sum of squares of genotypic variances per trait across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: VEnvir_tot            ! environmental variance per trait
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: locHs_Meta            ! observed Hs per locus as if parametric
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_locHs_Meta          ! sum of observed Hs per locus as if parametric
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_locHs_Meta         ! sum of squares of observed Hs per locus as if parametric
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: locNa_Meta            ! observed Ne per locus in current Meta
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_locNa_Meta          ! sum of observed Ne per locus across Metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_locNa_Meta         ! sum of squares of observed Ne per locus across Metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: locNasel_Meta         ! observed Ne per selected locus in current Meta
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S_locNasel_Meta       ! sum of observed Ne per selected locus across Metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Sq_locNasel_Meta      ! sum of squares of observed Ne per selected locus across Metapopulations


  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: VGenot_pop              ! genotypic variance per trait et pop.
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: S_VGenot_pop            ! sum of genotypic variances per trait & pop across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: Sq_VGenot_pop           ! sum od squares of genotypic variances per trait & pop  across metapopulations
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: genotypic_scale         ! genotypic values per locus
!  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: genotypic_scale         ! genotypic values per locus

  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: ini_link_trait    ! genotypic covariance due to linkage
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: end_link_trait    ! genotypic covariance due to linkage

  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: neu_meta_freq     ! allele frequencies at neutral loci in metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: S_neu_meta_freq   ! sum of allele frequencies at neutral loci in metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: Sq_neu_meta_freq  ! SQ of allele frequencies at neutral loci in metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: sel_meta_freq     ! allele frequencies at selected loci in metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: S_sel_meta_freq   ! sum of allele frequencies at selected loci in metapopulation
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: Sq_sel_meta_freq  ! SQ of allele frequencies at selected loci in metapopulation

  REAL(KIND=dp):: chisq_threshold    ! chi2 significativity threshold 
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: chisq_before_SA   ! CHI squared test HWE across loci and populations before SA
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: chisq_after_SA    ! CHI squared test HWE across loci and populations after SA

  REAL, ALLOCATABLE, DIMENSION(:) :: mean_additive_value, sd_additive_value  ! mean additive value and sqrt(additive variance) realized (ie after computation) 

  CHARACTER(LEN=10) :: ctime    ! stores system time
  CHARACTER(LEN=50) :: heading  ! stores inputs headings

! parameters for tunning SIMULATED ANNEALING ROUTINE FOR Fst
! Tsteps_Fst: maximum # of loops before reaching cooled state
!  INTEGER, PARAMETER :: Tsteps_Fst=1600  ! Initial (100)
  INTEGER, PARAMETER :: Tsteps_Fst=600  ! Initial (100)
! Nover_Fst: maximum # of loops per temperature level
!  INTEGER, PARAMETER :: Nover_Fst=800   ! Initial (50)
  INTEGER, PARAMETER :: Nover_Fst=300   ! Initial (50)
! Nlimit_Fst: # of successful hits per temperature level
  !INTEGER, PARAMETER :: Nlimit_Fst=300  ! Initial (20)
  INTEGER, PARAMETER :: Nlimit_Fst=100  ! Initial (20)
  REAL(KIND=dp), PARAMETER :: Tini_Fst=1.0_dp   ! Initial temperature
  REAL(KIND=dp), PARAMETER :: Tfactr_Fst=0.9_dp ! Cooling rate of temperature change
! Cooling rate: 0.9 (medium), 0.98 (slow), 0.99 (very slow)

! parameters for tunning SIMULATED ANNEALING ROUTINE FOR Qst
! Tsteps_Qst: maximum # of loops before reaching cooled state
!  INTEGER, PARAMETER :: Tsteps_Qst=1600  ! Initial (100)
  INTEGER, PARAMETER :: Tsteps_Qst=600  ! Initial (100)
! Nover_Qst: maximum # of loops per temperature level
!  INTEGER, PARAMETER :: Nover_Qst=800   ! Initial (50)
  INTEGER, PARAMETER :: Nover_Qst=300   ! Initial (50)
! Nlimit_Qst: # of successful hits per temperature level
!  INTEGER, PARAMETER :: Nlimit_Qst=300   ! Initial (20)
  INTEGER, PARAMETER :: Nlimit_Qst=100   ! Initial (20)
  REAL(KIND=dp), PARAMETER :: Tini_Qst=1.0_dp   ! Initial temperature
  REAL(KIND=dp), PARAMETER :: Tfactr_Qst=0.9_dp ! Cooling rate of temperature change

! parameters for tunning SIMULATED ANNEALING ROUTINE FOR Fst by locus
! Tsteps_locFst: maximum # of loops before reaching cooled state
!  INTEGER, PARAMETER :: Tsteps_locFst=400   ! Initial (100)
  INTEGER, PARAMETER :: Tsteps_locFst=200   ! Initial (100)
! Nover_locFst: maximum # of loops per temperature level
!  INTEGER, PARAMETER :: Nover_locFst=200    ! Initial (50)
  INTEGER, PARAMETER :: Nover_locFst=100    ! Initial (50)
! Nlimit_locFst: # of successful hits per temperature level
!  INTEGER, PARAMETER :: Nlimit_locFst=80    ! Initial (20)
  INTEGER, PARAMETER :: Nlimit_locFst=30    ! Initial (20)

  REAL(KIND=dp), PARAMETER :: Tini_locFst=1.0_dp    ! Initial temperature
  REAL(KIND=dp), PARAMETER :: Tfactr_locFst=0.9_dp  ! Cooling rate of temperature change
!==================================
  END MODULE metatrom_v7_param
!==================================

!==================================
MODULE random
!==================================

! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma
!     Chi-squared                     random_chisq
!     Exponential                     random_exponential
!     Weibull                         random_Weibull
!     Beta                            random_beta
!     t                               random_t
!     Multivariate normal             random_mvnorm
!     Generalized inverse Gaussian    random_inv_gauss
!     Poisson                         random_Poisson
!     Binomial                        random_binomial1   *
!                                     random_binomial2   *
!     Negative binomial               random_neg_binomial
!     von Mises                       random_von_Mises
!     Cauchy                          random_Cauchy
!
!  Generate a random ordering of the integers 1 .. N
!                                     random_order
!     Initialize (seed) the uniform random number generator for ANY compiler
!                                     seed_random_number

!     Lognormal - see note below.

!  ** Two functions are provided for the binomial distribution.
!  If the parameter values remain constant, it is recommended that the
!  first function is used (random_binomial1).   If one or both of the
!  parameters change, use the second function (random_binomial2).

! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.

! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.

! The module uses the following functions which are included here:
! bin_prob to calculate a single binomial probability
! lngamma  to calculate the logarithm to base e of the gamma function

! Some of the code is adapted from Dagpunar's book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! In most of Dagpunar's routines, there is a test to see whether the value
! of one or two floating-point parameters has changed since the last call.
! These tests have been replaced by using a logical variable FIRST.
! This should be set to .TRUE. on the first call using new values of the
! parameters, and .FALSE. if the parameter values are the same as for the
! previous call.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lognormal distribution
! If X has a lognormal distribution, then log(X) is normally distributed.
! Here the logarithm is the natural logarithm, that is to base e, sometimes
! denoted as ln.  To generate random variates from this distribution, generate
! a random deviate from the normal distribution with mean and variance equal
! to the mean and variance of the logarithms of X, then take its exponential.

! Relationship between the mean & variance of log(X) and the mean & variance
! of X, when X has a lognormal distribution.
! Let m = mean of log(X), and s^2 = variance of log(X)
! Then
! mean of X     = exp(m + 0.5s^2)
! variance of X = (mean(X))^2.[exp(s^2) - 1]

! In the reverse direction (rarely used)
! variance of log(X) = log[1 + var(X)/(mean(X))^2]
! mean of log(X)     = log(mean(X) - 0.5var(log(X))

! N.B. The above formulae relate to population parameters; they will only be
!      approximate if applied to sample values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Version 1.13, 2 October 2000
! Changes from version 1.01
! 1. The random_order, random_Poisson & random_binomial routines have been
!    replaced with more efficient routines.
! 2. A routine, seed_random_number, has been added to seed the uniform random
!    number generator.   This requires input of the required number of seeds
!    for the particular compiler from a specified I/O unit such as a keyboard.
! 3. Made compatible with Lahey's ELF90.
! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
! 5. INTENT for array f corrected in random_mvnorm.

!     Author: Alan Miller
!     e-mail: amiller @ bigpond.net.au

IMPLICIT NONE
REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      vsmall = TINY(1.0), vlarge = HUGE(1.0)
PRIVATE            :: integral
INTEGER, PARAMETER :: doubleprec = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal



FUNCTION random_gamma(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s > one) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < one) THEN
  fn_val = random_gamma2(s, first)
ELSE
  fn_val = random_exponential()
END IF

RETURN
END FUNCTION random_gamma



FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (first) THEN
  d = s - one/3.
  c = one/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = random_normal()
    v = (one + c*x)**3
    IF (v > zero) EXIT
  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < one - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1



FUNCTION random_gamma2(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

!     Local variables
REAL       :: r, x, w
REAL, SAVE :: a, p, c, uf, vr, d

IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2


FUNCTION random_exponential() RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

REAL  :: fn_val

!     Local variable
REAL  :: r

DO
  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO

fn_val = -LOG(r)
RETURN

END FUNCTION random_exponential


FUNCTION random_binomial1(n, p, first) RESULT(ival)

! FUNCTION GENERATES A RANDOM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.

!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= REAL <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= INTEGER)
!    FIRST = .TRUE. for the first call using the current parameter values
!          = .FALSE. if the values of (n,p) are unchanged from last call

! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN)    :: p
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival

!     Local variables

INTEGER         :: ru, rd
INTEGER, SAVE   :: r0
REAL            :: u, pd, pu
REAL, SAVE      :: odds_ratio, p_r
REAL, PARAMETER :: zero = 0.0, one = 1.0

IF (first) THEN
  r0 = (n+1)*p
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)
END IF

CALL RANDOM_NUMBER(u)
u = u - p_r
IF (u < zero) THEN
  ival = r0
  RETURN
END IF

pu = p_r
ru = r0
pd = p_r
rd = r0
DO
  rd = rd - 1
  IF (rd >= 0) THEN
    pd = pd * (rd+1) / (odds_ratio * (n-rd))
    u = u - pd
    IF (u < zero) THEN
      ival = rd
      RETURN
    END IF
  END IF

  ru = ru + 1
  IF (ru <= n) THEN
    pu = pu * (n-ru+1) * odds_ratio / ru
    u = u - pu
    IF (u < zero) THEN
      ival = ru
      RETURN
    END IF
  END IF
END DO

!     This point should not be reached, but just in case:

ival = r0
RETURN

END FUNCTION random_binomial1



FUNCTION bin_prob(n, p, r) RESULT(fn_val)
!     Calculate a binomial probability

INTEGER, INTENT(IN) :: n, r
REAL, INTENT(IN)    :: p
REAL                :: fn_val

!     Local variable
REAL                :: one = 1.0

fn_val = EXP( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
              + r*LOG(p) + (n-r)*LOG(one - p) )
RETURN

END FUNCTION bin_prob



FUNCTION lngamma(x) RESULT(fn_val)

! Logarithm to base e of the gamma function.
!
! Accurate to about 1.e-14.
! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 28 February 1988

REAL (doubleprec), INTENT(IN) :: x
REAL (doubleprec)             :: fn_val

!       Local variables

REAL (doubleprec) :: a1 = -4.166666666554424D-02, a2 = 2.430554511376954D-03,  &
             a3 = -7.685928044064347D-04, a4 = 5.660478426014386D-04,  &
             temp, arg, product, lnrt2pi = 9.189385332046727D-1,       &
             pi = 3.141592653589793D0
LOGICAL   :: reflect

!       lngamma is not defined if x = 0 or a negative integer.

IF (x > 0.d0) GO TO 10
IF (x /= INT(x)) GO TO 10
fn_val = 0.d0
RETURN

!       If x < 0, use the reflection formula:
!               gamma(x) * gamma(1-x) = pi * cosec(pi.x)

10 reflect = (x < 0.d0)
IF (reflect) THEN
  arg = 1.d0 - x
ELSE
  arg = x
END IF

!       Increase the argument, if necessary, to make it > 10.

product = 1.d0
20 IF (arg <= 10.d0) THEN
  product = product * arg
  arg = arg + 1.d0
  GO TO 20
END IF

!  Use a polynomial approximation to Stirling's formula.
!  N.B. The real Stirling's formula is used here, not the simpler, but less
!       accurate formula given by De Moivre in a letter to Stirling, which
!       is the one usually quoted.

arg = arg - 0.5D0
temp = 1.d0/arg**2
fn_val = lnrt2pi + arg * (LOG(arg) - 1.d0 + &
                  (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - LOG(product)
IF (reflect) THEN
  temp = SIN(pi * x)
  fn_val = LOG(pi/temp) - fn_val
END IF
RETURN
END FUNCTION lngamma




SUBROUTINE integral(a, b, result, dk)

!     Gaussian integration of exp(k.cosx) from a to b.

REAL (doubleprec), INTENT(IN) :: dk
REAL, INTENT(IN)      :: a, b
REAL, INTENT(OUT)     :: result

!     Local variables

REAL (doubleprec)  :: xmid, range, x1, x2,                                    &
  x(3) = (/0.238619186083197_doubleprec, 0.661209386466265_doubleprec, 0.932469514203152_doubleprec/), &
  w(3) = (/0.467913934572691_doubleprec, 0.360761573048139_doubleprec, 0.171324492379170_doubleprec/)
INTEGER    :: i

xmid = (a + b)/2._doubleprec
range = (b - a)/2._doubleprec

result = 0._doubleprec
DO i = 1, 3
  x1 = xmid + x(i)*range
  x2 = xmid - x(i)*range
  result = result + w(i)*(EXP(dk*COS(x1)) + EXP(dk*COS(x2)))
END DO

result = result * range
RETURN
END SUBROUTINE integral



!==================================
END MODULE random
!==================================

!==================================
  MODULE metatrom_v7_stat
! basic statistic functions
!==================================
  SAVE
  INTEGER, PARAMETER, PRIVATE :: dp=SELECTED_REAL_KIND(14)

  CONTAINS
!----------------------------------
  FUNCTION real_mean(values)
! calculates mean of a sample vector
!----------------------------------
  IMPLICIT NONE
  REAL(KIND=dp) :: real_mean
  INTEGER :: n
  REAL(KIND=dp), DIMENSION(:) :: values
  n=SIZE(values)
  IF (n==0) THEN
    PRINT *, 'ERROR!!, attempt to calculate mean of zero elements in',&
    &' FUNCTION real_mean'
  ENDIF
  real_mean=SUM(values)/n
!----------------------------------
  END FUNCTION real_mean
!----------------------------------
!----------------------------------
  FUNCTION real_var(values)
! calculates variance of a sample vector
!----------------------------------
  IMPLICIT NONE
  REAL(KIND=dp) :: real_var
  INTEGER :: n
  REAL(KIND=dp), DIMENSION(:) :: values
  n=SIZE(values)
  IF (n<2) THEN
    PRINT *, 'ERROR!!, attempt to calculate variance of only one element in',&
    &' FUNCTION real_var'
  ENDIF
  IF(SUM(values**2)>((SUM(values)**2)/n))THEN
    real_var=(SUM(values**2)-(SUM(values)**2)/n)/(n-1)
  ELSE
    real_var=0.0d0
  END IF
!----------------------------------
  END FUNCTION real_var
!----------------------------------
!----------------------------------
  FUNCTION real_cov(data1,data2)
! calculates covariance between two
! equal sized sample vectors
!----------------------------------
  IMPLICIT NONE
  REAL(KIND=dp) :: real_cov
  INTEGER :: n
  REAL(KIND=dp), DIMENSION(:) :: data1,data2
  n=SIZE(data1)
  IF (n<2) THEN
    PRINT *,'ERROR!!, attempt to calculate covariance of only one pair in',&
    &' FUNCTION real_cov'
  ELSE IF (SIZE(data1)/=SIZE(data2)) THEN
    PRINT *,'ERROR!!, unequal number of elements for covariance in',&
    &' FUNCTION real_cov'
  ENDIF
  IF(SUM(data1*data2)==(SUM(data1)*SUM(data2))/n)THEN
    real_cov=0.0d0
  ELSE
    real_cov=(SUM(data1*data2)-(SUM(data1)*SUM(data2))/n)/n
  END IF
!----------------------------------
  END FUNCTION real_cov
!----------------------------------
!----------------------------------
  FUNCTION real_cor(data1,data2)
! calculates correlation between two
! equal sized sample vectors
!----------------------------------
  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND=dp) :: real_cor
  REAL(KIND=dp) :: var_1,var_2,covar_12
  REAL(KIND=dp), DIMENSION(:) :: data1,data2
  n=SIZE(data1)
  IF (n<2) THEN
    PRINT *, 'ERROR!!, attempt to calculate correlation of only one pair in',&
    &' FUNCTION real_cor'
  ELSE IF (SIZE(data1)/=SIZE(data2)) THEN
    PRINT *, 'ERROR!!, unequal number of elements for correlation in',&
    &' FUNCTION real_cor'
  ENDIF
  var_1=real_var(data1)
  var_2=real_var(data2)
  covar_12=(SUM(data1*data2)-(SUM(data1)*SUM(data2))/n)/n
  IF(var_1*var_2>0.0d0)THEN
    real_cor=covar_12/DSQRT(var_1*var_2)
  ELSE
    real_cor=0.0d0
  END IF
!----------------------------------
  END FUNCTION real_cor
!----------------------------------
!----------------------------------
  FUNCTION real_reg(data1,data2)
! calculates linear regression of
! data1 on data2: data1 = b*data2 + e
!----------------------------------
  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND=dp) :: real_reg
  REAL(KIND=dp), DIMENSION(:) :: data1,data2
  n=SIZE(data1)
  IF (n<2) THEN
    PRINT *, 'ERROR!!, attempt to calculate regression of only one pair in',&
      &' FUNCTION real_reg'
    STOP
  ELSE IF (SIZE(data1)/=SIZE(data2)) THEN
    PRINT *,'ERROR!!, unequal number of elements for regression in',&
      &' FUNCTION real_reg'
    STOP
  ENDIF
  real_reg=real_cov(data1,data2)/real_var(data2)
!----------------------------------
  END FUNCTION real_reg
!----------------------------------
!----------------------------------
  FUNCTION mean_rep(data1,elements)
! calculates mean from sum and
! # of elements in summation
!----------------------------------
  IMPLICIT NONE
  INTEGER :: elements
  REAL(KIND=dp) :: mean_rep
  REAL(KIND=dp) :: data1
  mean_rep=data1/elements
!----------------------------------
  END FUNCTION mean_rep
!----------------------------------
!----------------------------------
  FUNCTION vare_rep(data1,data2,elements)
! calculates stand.deviation from sum,
! sum of squares, and # of elements
! in summations
!----------------------------------
  IMPLICIT NONE
  INTEGER :: elements
  REAL(KIND=dp) :: vare_rep
  REAL(KIND=dp) :: data1,data2
  IF(elements<=1)THEN
    vare_rep=0.0_dp
  ELSE
    vare_rep=(data1-((data2**2)/elements))/(elements-1)
  END IF
!----------------------------------
  END FUNCTION vare_rep
!----------------------------------
!==================================
  FUNCTION n_rec(row,col,max_dim)
! row must be <= than col
! max_dim is the DIM of the matrix
!==================================
  IMPLICIT NONE
  INTEGER :: n_rec
  INTEGER :: max_dim
  INTEGER :: row,col
  n_rec=(((row-1)*max_dim)-(((row*row-row)/2)))+col
!==================================
  END FUNCTION n_rec
!==================================

!----------------------------------
FUNCTION ran_gamma(s, first) RESULT(fn_val)
!----------------------------------
USE random
IMPLICIT NONE

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS).

! Generates a random gamma deviate for shape parameter s >= 1.

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (s < 1.0) THEN
   WRITE(*, *) 'Shape parameter must be >= 1'
   STOP
END IF

IF (first) THEN
  d = s - 1./3.
  c = 1.0/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = random_normal()
    v = (1.0 + c*x)**3
    IF (v > 0.0) EXIT
  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < 1.0 - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < 0.5*x**2 + d*(1.0 - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
!----------------------------------
END FUNCTION ran_gamma
!----------------------------------

!==================================
  END MODULE metatrom_v7_stat
!==================================

!==================================
  MODULE metatrom_v7_normal_table
!==================================
  USE metatrom_v7_param
  SAVE

  CONTAINS

!----------------------------------
  REAL(KIND=dp) FUNCTION gcef(q,ix)
! calculates normal deviate x from lower tail proportion p
! LSR 20/07/04 (from JAW)
!----------------------------------
  IMPLICIT NONE
  INTEGER, OPTIONAL :: ix
  REAL(KIND=dp) :: q
  REAL(KIND=dp) :: p
  REAL(KIND=dp) :: pp
  REAL(KIND=dp) :: u
  REAL(KIND=dp) :: t
  REAL(KIND=dp) :: x
  REAL(KIND=dp) :: zero=0.0_dp
  REAL(KIND=dp) :: one=1.0_dp
  REAL(KIND=dp) :: half=0.5_dp
  REAL(KIND=dp) :: a1=2.515517_dp
  REAL(KIND=dp) :: a2=0.802853_dp
  REAL(KIND=dp) :: a3=0.010328_dp
  REAL(KIND=dp) :: b1=1.432788_dp
  REAL(KIND=dp) :: b2=0.189269_dp
  REAL(KIND=dp) :: b3=0.001308_dp

  p=REAL(q,kind=dp)
  SELECT CASE((p>=zero).and.(p<=one))
  CASE(.true.)
  IF(PRESENT(ix)) ix=0
  IF(p>=one) THEN
    gcef=7.0_dp                 !   JAW
  ELSE IF(p<1.0e-10_dp) THEN
    gcef=-7.0_dp                !   JAW
  ELSE
! calculates for probabilities gt half by using remainder
    IF(p>half) THEN
      pp=one-p
    ELSE
      pp=p
    END IF
    u=LOG(one/(pp*pp))
    t=DSQRT(u)
    x=(a1+(a2*t)+(a3*u))
    x=x/(one+(b1*t)+(b2*u)+b3*EXP(3.0_dp*LOG(t)))
    IF(p>half) THEN
      gcef=t-x
    ELSE
      gcef=x-t
    END IF
  END IF

  CASE DEFAULT
  IF(PRESENT(ix)) ix=1
    PRINT *, 'GCEF: probability out of bounds; ignore results'
    gcef=zero
  END SELECT
!----------------------------------
  END FUNCTION gcef
!----------------------------------
!----------------------------------
  REAL(KIND=dp) FUNCTION alea(idum)
!----------------------------------
! Long period (>2*10^18) random number generator of L'Ecuyer
! with Bays-Durham shuffle and added safeguards. Returns a
! uniform random deviate between 0.0 and 1.0 (exclusive of the
! endpoints values). Call with idum a negative integer to
! initialize; thereafter, do not alter idum between successive
! deviates in a sequence. rnmx should approximate the largest
! floating value that is less than 1
!
! From NUMERICAL RECIPES IN FORTRAN 77, pp:272
!
    INTEGER idum
    INTEGER im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
    REAL am_alea,eps,rnmx
    PARAMETER (im1=2147483563)
    PARAMETER (im2=2147483399)
    PARAMETER (am_alea=1.0_dp/im1)
    PARAMETER (imm1=im1-1)
    PARAMETER (ia1=40014)
    PARAMETER (ia2=40692)
    PARAMETER (iq1=53668)
    PARAMETER (iq2=52774)
    PARAMETER (ir1=12211)
    PARAMETER (ir2=3791)
    PARAMETER (ntab=32)
    PARAMETER (ndiv=1+imm1/ntab)
    PARAMETER (eps=1.2e-7)
    PARAMETER (rnmx=1.0_dp-eps)
    INTEGER idum2,j_alea,k_alea,iv(ntab),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/,iv/ntab*0/,iy/0/
    IF(idum.LE.0)THEN
      idum=MAX(-idum,1)
      idum2=idum
      DO j_alea=ntab+8,1,-1
        k_alea=idum/iq1
        idum=ia1*(idum-k_alea*iq1)-k_alea*ir1
        IF(idum.LT.0) idum=idum+im1
        IF(j_alea.LE.ntab) iv(j_alea)=idum
      END DO
      iy=iv(1)
    END IF
    k_alea=idum/iq1
    idum=ia1*(idum-k_alea*iq1)-k_alea*ir1
    IF(idum.LT.0) idum=idum+im1
    k_alea=idum2/iq2
    idum2=ia2*(idum2-k_alea*iq2)-k_alea*ir2
    IF(idum2.LT.0) idum2=idum2+im2
    j_alea=1+iy/ndiv
    iy=iv(j_alea)-idum2
    iv(j_alea)=idum
    IF(iy.LT.1) iy=iy+imm1
    alea=MIN(am_alea*iy,rnmx)
!----------------------------------
  END FUNCTION alea
!----------------------------------
!==================================
  END MODULE metatrom_v7_normal_table
!==================================
!==================================
  MODULE metatrom_v7_routines
!==================================
  USE metatrom_v7_param
  USE metatrom_v7_stat
  USE metatrom_v7_normal_table
  USE random
  SAVE

  CONTAINS
!----------------------------------
  SUBROUTINE allelic_freq
! contribution of loci to additive variance,
! allelic frequencies, mean loci effect and
! additive value per indiv
!----------------------------------

INTEGER :: counter
INTEGER :: progress_counter_1
REAL(KIND=dp) :: sum_loci_contrib ! sum of loci contributions
REAL(KIND=dp) :: x_al1, x_al2 ! frequencies of alleles according the gamma distribution (ie before to be rescaled)
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: al1, al2 ! alleles 1 and 2
REAL(KIND=dp) , ALLOCATABLE, DIMENSION(:) ::  loci_contrib_to_VA ! contribution of each locus to the additive variance
REAL(KIND=dp) , ALLOCATABLE, DIMENSION(:) ::  VAl ! mean effect for each locus (=alpha in Falconer), addtive variance for a given locus
REAL(KIND=dp) , ALLOCATABLE, DIMENSION(:) ::  additive_value ! additive value computed according allele frequencies and genotypic value (=a in Falconer). Note that mean effect = genotypic value because there is no dominancy
REAL(KIND=dp) :: ldum

INTEGER, ALLOCATABLE, DIMENSION(:) :: count_sel

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: al1_neutral, al2_neutral ! alleles 1 and 2 for neutral markers

INTEGER :: i,j,k,l,m
REAL(KIND=4), ALLOCATABLE, DIMENSION(:) :: vec_freq
REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: vec_sq_freq

FORALL(i=1:num_traits+1) mean_additive_value(i)=0.0_dp
FORALL(i=1:num_traits+1) sd_additive_value(i)=0.0_dp

! selected loci
  ALLOCATE(count_sel(num_traits))
  FORALL(i=1:num_traits) count_sel(i)=0
  DO i=1,tot_sel_loci
!   determining position of selected loci in the selected genome
    count_sel(sel_gene(i)%funct)=count_sel(sel_gene(i)%funct)+1
    posit_sel_loci(sel_gene(i)%funct,count_sel(sel_gene(i)%funct))=i
  END DO

  DEALLOCATE(count_sel)

! Allocate a pop and a position in this pop to individuals
  jlast=0
  DO i=1,num_pop_tot
    DO j=1,ind_per_pop(i)
      jlast=jlast+1
      indv(jlast)%pop=i
      posit_ind_pop(i,j)=jlast
    END DO
  END DO

! Allocate an allele number to selected loci : all selected loci are bi-allelics

DO i=1, tot_sel_loci
sel_gene(i)%nalle=2
END DO

! Size definition of vectors and tables
ALLOCATE(vec_freq(tot_sel_loci))
ALLOCATE(VAl(tot_sel_loci))
!ALLOCATE(mean_effect_l(tot_sel_loci))
ALLOCATE(loci_contrib_to_VA(tot_sel_loci))
ALLOCATE(al1(tot_sel_loci,num_tot_ind))
ALLOCATE(al2(tot_sel_loci,num_tot_ind))
ALLOCATE(additive_value(num_tot_ind))



DO i=1, num_traits
!write(*,*) 'trait n°', i


	IF (target_additive_value(i)%sd>0) THEN

		! Loop restart until to find an additive value and additive variance very close of the targeted values
		progress_counter_1 = 0
		!DO WHILE (ABS((target_additive_value(i)%mean - mean_additive_value(i))/target_additive_value(i)%mean) &
		!>0.01_dp .OR. ABS((target_additive_value(i)%sd-sd_additive_value(i))/target_additive_value(i)%sd)>0.01_dp)

		DO WHILE (ABS(mean_additive_value(i))>target_additive_value(i)%mean_threshold &
		.OR. ABS((target_additive_value(i)%sd-sd_additive_value(i))/target_additive_value(i)%sd)>0.01_dp)

			progress_counter_1 = progress_counter_1+1
			IF (MODULO(progress_counter_1, 5000)==0) WRITE(*,*) '[allelic_freq] search additive variance and value close to target, &
					iteration: ', progress_counter_1

			! initialization of intermediate parameters
			FORALL(i=1:tot_sel_loci) vec_freq(i)=0.0_dp
			FORALL(i=1:tot_sel_loci) VAl(i)=0.0_dp
			!FORALL(i=1:tot_sel_loci) mean_effect_l(i)=0.0_dp
			FORALL(i=1:tot_sel_loci) loci_contrib_to_VA(i)=0.0_dp
			FORALL(i=1:tot_sel_loci, j=1:num_tot_ind) al1(i,j)=0.0_dp
			FORALL(i=1:tot_sel_loci, j=1:num_tot_ind) al2(i,j)=0.0_dp
			FORALL(j=1:num_tot_ind) additive_value(j)=0.0_dp

			IF(i==1) THEN
				counter=0
			ELSE
				counter = trait_sel_loci(i-1)
			END IF
			!write(*,*) 'initial counter', counter
			! Initialization of the gamma distribution
			vec_freq = ran_gamma(distribution_shape, .TRUE.)

			!-----------------------------------------------
			! Computation of the contributions of loci to VA
			!-----------------------------------------------
			! Sampling a contribution to VA for each locus in a gamma distribution
			sum_loci_contrib=0.0

			DO j=1, tot_sel_loci
				IF (sel_gene(j)%funct==i) THEN
					counter=counter+1
					loci_contrib_to_VA(counter)=ran_gamma(distribution_shape, .FALSE.)
					!Write(*,*) 'loci contribution to VA before rescaling', j, loci_contrib_to_VA(j)
				END IF
			END DO
				!WRITE(*,*) 'contrib loci to VA before rescale :', loci_contrib_to_VA
			DO j=1,counter
				sum_loci_contrib=sum_loci_contrib+loci_contrib_to_VA(j)
			END DO
			!WRITE(*,*) 'sum of contrib :', sum_loci_contrib

		 	!DO j = 1, loci_nb
			! Rescaling contributions to VA in order that sum of relative contributions =1
			!loci_contrib_to_VA(j)=loci_contrib_to_VA(j)/sum_loci_contrib
		 	!WRITE(*,*) 'contrib loci to VA after rescale :', loci_contrib_to_VA(j)
			! Computation of additive variance at each locus
			!VAl(j)=loci_contrib_to_VA(j)*target_sd_additive_value**2
			!WRITE(*,*) 'VAl :', VAl(j)
			!END DO


			! Rescaling contributions to VA in order that sum of relative contributions =1
			loci_contrib_to_VA=loci_contrib_to_VA/sum_loci_contrib
		 	!WRITE(*,*) 'contrib loci to VA after rescale :', loci_contrib_to_VA
			! Computation of additive variance at each locus
			VAl=loci_contrib_to_VA*target_additive_value(i)%sd**2

			!counter = counter - trait_sel_loci(i)
			DO j=counter-trait_sel_loci(i)+1,counter
				!IF(sel_gene(j)%funct==i) THEN
					!counter=counter+1
					loci_contrib_VA(j)=VAl(j)
				!END IF
			END DO
			!write(*,*) 'allelic_freq, loci contrib to Va', loci_contrib_VA
			! TEMPORARY to solve bugs
						!counter = counter - trait_sel_loci(i)
						!DO j=1, tot_sel_loci
							!IF (sel_gene(j)%funct==i) THEN
								!counter=counter+1
								!loci_contrib_VA(j)=target_additive_value(i)%sd**2/trait_sel_loci(i)
							!END IF
						!END DO
			!write(*,*) 'allelic_freq, loci contrib to Va', loci_contrib_VA
			!-----------------------------------------------
			! End computation of the contributions of loci to VA
			!-----------------------------------------------
			!-----------------------------------------------
			! Computation of allelic frequencies for each locus
			!-----------------------------------------------
			x_al1=0.0
			x_al2=0.0
			counter = counter - trait_sel_loci(i)
			DO j=1, tot_sel_loci
				IF (sel_gene(j)%funct==i) THEN
					counter=counter+1
					! Sampling allelic frequencies in a gamma distribution
					x_al1 = ran_gamma(distribution_shape, .FALSE.)
					x_al2 = ran_gamma(distribution_shape, .FALSE.)

					! Rescaling allelic frequencies in order that sum of p+q=1
					vec_freq(counter) = x_al1/(x_al1+x_al2)
					! Computation of a mean effect at each locus
					genotypic_scale(counter,1)=sqrt(VAl(counter)/(2*vec_freq(counter)*(1-vec_freq(counter))))
				END IF
					!WRITE(*,*) 'VAl', j, VAl(j)
					!WRITE(*,*) 'p', j, vec_freq(i,j)
					!WRITE(*,*) 'alpha', j, mean_effect_l(j)

				!	WRITE(31,*) vec_freq(i,j)
			END DO

			! TEMPORARY to solve bugs
			!genotypic_scale=0
						!DO j=1, tot_sel_loci
							!genotypic_scale(j,1)=1
						!END DO
			!write(*,*) 'genotypic_scale', genotypic_scale(:,1)

			!IF (i==1) THEN
			! OPEN(UNIT=3,FILE='mean_effects_1.txt',STATUS='unknown')
			!DO j = 1, num_loci
			!write(3,*) mean_effect_l(j)
			!END DO
			! close(3)
			!ELSE IF(i==2) THEN
			! OPEN(UNIT=4,FILE='mean_effects_2.txt',STATUS='unknown')
			!DO j = 1, num_loci
			!write(4,*) mean_effect_l(j)
			!END DO
			! close(4)
			!END IF
			!-----------------------------------------------
			! End computation of allelic frequencies for each locus
			!-----------------------------------------------

			!-----------------------------------------------
			! Sampling of alleles at each locus and for all individuals
			!-----------------------------------------------

			!IF(i==1) THEN
				!OPEN(UNIT=2,FILE='alleles_caracter1.txt',STATUS='unknown')

			DO k=1,num_tot_ind
				counter = counter - trait_sel_loci(i)
				DO j=1, tot_sel_loci
					IF (sel_gene(j)%funct==i) THEN
						counter=counter+1
						al1(counter,k) = random_binomial1(1,vec_freq(counter), .TRUE.)+1
						sel_genome(1,counter,k) = al1(counter,k)
						al2(counter,k) = random_binomial1(1,vec_freq(counter), .TRUE.)+1
						sel_genome(2,counter,k) = al2(counter,k)
						!WRITE(*,*)i, j, al1(counter,k), al2(counter,k)

						!write(2,*) i, counter, k, vec_freq(counter), al1(counter,k), al2(counter,k) 
					END IF
				END DO
			END DO
			!write(*,*) 'sel_genome', sel_genome(2,:,1:5)
						!close(2)

			!ELSE IF (i==2) THEN
			!	OPEN(UNIT=3,FILE='alleles_caracter2.txt',STATUS='unknown')
			!
			!			DO k=1,num_tot_ind
			!				DO j=1,num_loci
			!
			!				IF (i==1) THEN
			!					al1(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
			!					sel_genome(1,j,k) = al1(j,k)
			!					al2(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
			!					sel_genome(2,j,k) = al2(j,k)
			!					!WRITE(*,*)i, j, al1(i,j), al2(i,j)
			!		
			!				ELSE IF (i==2) THEN 
			!					al1(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
			!					sel_genome(1,j+trait_sel_loci(i-1),k) = al1(j,k)
			!					al2(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
			!					sel_genome(2,j+trait_sel_loci(i-1),k) = al2(j,k)
			!				END IF 
			!
			!				write(3,*) i, j, k, vec_freq(j), al1(j,k), al2(j,k) 
			!				END DO
			!			END DO
			!			close(3)
			!END IF 

			!-----------------------------------------------
			! End sampling of alleles at each locus and for all individuals
			!-----------------------------------------------

			!-----------------------------------------------
			! Computation of additive value for each individual
			!-----------------------------------------------
			!	OPEN(UNIT=5,FILE='verif_VA.txt',STATUS='unknown')
			DO k=1,num_tot_ind
				additive_value(k)=0
				counter = counter - trait_sel_loci(i)
				DO j=1, tot_sel_loci
					IF (sel_gene(j)%funct==i) THEN
						counter=counter+1			
						! For homozygous 11, additive value = additive value - a (Falconer)
						IF (al1(counter,k)+al2(counter,k)==2) THEN
							additive_value(k) = additive_value(k) - genotypic_scale(counter,1)
						! For homozygous 22, additive value = additive value + a (Falconer)
						ELSE IF (al1(counter,k)+al2(counter,k)==4) THEN
							additive_value(k) = additive_value(k) + genotypic_scale(counter,1)
						! For heterozygous, additive value = additive value + 0 (Falconer)
						END IF
					!write(*,*) i, j, k, al1(counter,k)+al2(counter,k), genotypic_scale(counter,1), additive_value(k)
					!write(*,*) i, j, k, al1(counter,k)+al2(counter,k)
					END IF 
				END DO
				!IF(i==1 .AND. k==5) write(*,*) 'mean_effect_l : ',  genotypic_scale(1:20,1)
				additive_value(k)=additive_value(k)+target_additive_value(i)%mean
			END DO
			!write(*,*) 'counter', counter
			! close(5)
			!-----------------------------------------------
			! End computation of additive value for each individual
			!-----------------------------------------------
			!WRITE(*,*) 'additive value', additive_value
			! Computation of mean additive value and sqrt(additive variance)
			mean_additive_value(i)= real_mean(additive_value)
			sd_additive_value(i)= sqrt(real_var(additive_value))
			!write(*,*) 'additive_value', additive_value
		END DO

!			write(*,*) 'mean additive value', mean_additive_value(i)
!			write(*,*) 'sd additive value', sd_additive_value(i)

!write(*,*) 'test fin de boucle while'

	!WRITE(*,*) 'additive value', i, 15, additive_value(15)-target_additive_value(i)%mean

!  ALLOCATE(vec_sq_freq(num_loci))

!  FORALL(j=1:num_loci) vec_sq_freq(j)=0.0_dp

!  DEALLOCATE(vec_sq_freq)



		!WRITE(*,*) 'mean additivity realized : ', character1_mean_additive_value
		!WRITE(*,*) 'sd additivity realized : ', character1_sd_additive_value


	! If target_additive_value(i)%sd = 0
	ELSE

		! Initialization of the gamma distribution
	  	vec_freq = ran_gamma(distribution_shape, .TRUE.)

		!-----------------------------------------------
		! Computation of the contributions of loci to VA
		!-----------------------------------------------
		! Sampling a contribution to VA for each locus in a gamma distribution
		DO j = 1, num_loci
			loci_contrib_to_VA(j)=ran_gamma(distribution_shape, .FALSE.)
		END DO

		!	WRITE(*,*) 'contrib loci à la VA avant rescale :', loci_contrib_to_VA
		sum_loci_contrib=sum(loci_contrib_to_VA)
		!	WRITE(*,*) 'somme des contrib :', sum_loci_contrib
		! Rescaling contributions to VA in order that sum of relative contributions =1
		loci_contrib_to_VA=loci_contrib_to_VA/sum_loci_contrib
	 	! 	WRITE(*,*) 'contrib loci à la VA après rescale :', loci_contrib_to_VA
		! Computation of additive variance at each locus
		VAl=loci_contrib_to_VA*target_additive_value(i)%sd**2
		!	WRITE(*,*) 'VAl :', VAl

		!-----------------------------------------------
		! End computation of the contributions of loci to VA
		!-----------------------------------------------

		!-----------------------------------------------
		! Computation of allelic frequencies for each locus
		!-----------------------------------------------
	  	DO j = 1, num_loci
		! Sampling allelic frequencies in a gamma distribution
			x_al1 = ran_gamma(distribution_shape, .FALSE.)
			x_al2 = ran_gamma(distribution_shape, .FALSE.)

		! Rescaling allelic frequencies in order that sum of p+q=1
			vec_freq(j) = x_al1/(x_al1+x_al2)
		! Computation of a mean effect at each locus
!			mean_effect_l(j) = 0

				IF (i==1) THEN
				genotypic_scale(j,1)=0
				ELSE IF(i==2) THEN
				genotypic_scale(j+trait_sel_loci(1),1)=0
				END IF

			!	WRITE(*,*) j, vec_freq(i,j)
			!	WRITE(31,*) vec_freq(i,j)
		END DO

		!-----------------------------------------------
		! Computation of allelic frequencies for each locus
		!-----------------------------------------------

		!-----------------------------------------------
		! Sampling of alleles at each locus and for all individuals
		!-----------------------------------------------

!	OPEN(UNIT=2,FILE='alleles_caracter1.txt',STATUS='unknown')

		DO k=1,num_tot_ind
			DO j=1,num_loci

			IF (i==1) THEN
				al1(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
				sel_genome(1,j,k) = al1(j,k)
				al2(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
				sel_genome(2,j,k) = al2(j,k)
				!WRITE(*,*)i, j, al1(i,j), al2(i,j)
			ELSE IF (i==2) THEN
				al1(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
				sel_genome(1,j+trait_sel_loci(i-1),k) = al1(j,k)
				al2(j,k) = random_binomial1(1,vec_freq(j), .TRUE.)+1
				sel_genome(2,j+trait_sel_loci(i-1),k) = al2(j,k)
			END IF
		
!			write(2,*) i, j, k, vec_freq(j), al1(j,k), al2(j,k) 
			END DO
		END DO
!		close(2)



		!-----------------------------------------------
		! End sampling of alleles at each locus and for all individuals
		!-----------------------------------------------

		!-----------------------------------------------
		! Attribution of additive value for each individual
		!-----------------------------------------------
		DO k=1,num_tot_ind
			additive_value(k)=target_additive_value(i)%mean
		END DO

		!-----------------------------------------------
		! End attribution of additive value for each individual
		!-----------------------------------------------

		! Computation of mean additive value and sqrt(additive variance)
		mean_additive_value(i)= target_additive_value(i)%mean
		sd_additive_value(i)=0

	
!		WRITE(*,*) 'mean additivity realized : ', mean_additive_value(i)
!		WRITE(*,*) 'sd additivity realized : ', sd_additive_value(i)


!  ALLOCATE(vec_sq_freq(num_loci))

!  DO j=1,num_loci
!    FORALL(j=1:num_loci) vec_sq_freq(j)=0.0_dp
!      vec_sq_freq(j)=vec_freq(j)*vec_freq(j)+(1-vec_freq(j))*(1-vec_freq(j))
!    obs_locHt_Meta(j)=1.0_dp-SUM(vec_sq_freq)
!    locNa_Meta(j)=1.0_dp/SUM(vec_sq_freq)
!  END DO
!  IF(tot_neu_loci>0) Ht_Meta=SUM(obs_locHt_Meta)/tot_neu_loci

!  DEALLOCATE(vec_sq_freq)

	END IF



END DO

!DO j = 1, tot_sel_loci
!	WRITE(*,*) j, genotypic_scale(j,1)
!END DO

! DEALLOCATING MEMORY BEFORE EXIT
DEALLOCATE(vec_freq)
DEALLOCATE(VAl)
!DEALLOCATE(mean_effect_l)
DEALLOCATE(loci_contrib_to_VA)
DEALLOCATE(al1)
DEALLOCATE(al2)
DEALLOCATE(additive_value)

! ***********************************************************************
! ***********************************************************************
! For neutral markers (microsatelites + SNP)
! ***********************************************************************
! ***********************************************************************


ALLOCATE(al1_neutral(num_tot_ind, tot_neu_loci))
ALLOCATE(al2_neutral(num_tot_ind, tot_neu_loci))
ldum=REAL(idum)
al1_neutral=0.0
al2_neutral=0.0

!	OPEN(UNIT=6,FILE='alleles_microsat.txt',STATUS='unknown')
IF(tot_neu_loci>0)THEN
	DO i=1, num_tot_ind
		DO m=1, tot_neu_loci
			CALL random_number(ldum)
			al1_neutral(i,m)=FLOOR(ldum*max_neu_allele)+11
			neu_genome(1,m,i) = al1_neutral(i,m)
			CALL random_number(ldum)
			al2_neutral(i,m)=FLOOR(ldum*max_neu_allele)+11
			neu_genome(2,m,i) = al2_neutral(i,m)
!write(6,*) al1_neutral(i,m), al2_neutral(i,m)
		END DO
	END DO
END IF
!			close(6)


DEALLOCATE(al1_neutral)
DEALLOCATE(al2_neutral)

!write(*,*) 'end allelic freq'
!----------------------------------
  END SUBROUTINE allelic_freq
!----------------------------------

!----------------------------------
  SUBROUTINE meta_Ht_neu
! calculates Ht and Na across metapopulation
! these values remain cte
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j,k, l
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: copies
  INTEGER :: allele
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: vec_sq_freq
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: vec_freq

  ALLOCATE(copies(tot_neu_loci, max_neu_allele))

  ALLOCATE(vec_sq_freq(max_neu_allele))
  ALLOCATE(vec_freq(tot_neu_loci,max_neu_allele))

  FORALL(i=1:tot_neu_loci,j=1:max_neu_allele) vec_freq(i,j)=0.0_dp
  FORALL(i=1:tot_neu_loci,j=1:max_neu_allele) copies(i,j)=0.0_dp

  DO l=1, 2
	DO k=1, num_tot_ind
		DO j=1,tot_neu_loci
			allele = neu_genome(l,j,k)-10
			copies(j, allele) = copies(j, allele)+1
		END DO
	 END DO
  END DO

  DO i=1,tot_neu_loci
    FORALL(i=1:max_neu_allele) vec_sq_freq(i)=0.0_dp
    DO j=1,neu_gene(i)%nalle
      vec_freq(i,j)=REAL(copies(i,j),dp)/num_tot_haplot
      vec_sq_freq(j)=vec_freq(i,j)*vec_freq(i,j)
    END DO
    obs_locHt_Meta(i)=1.0_dp-SUM(vec_sq_freq)
    locNa_Meta(i)=1.0_dp/SUM(vec_sq_freq)
  END DO
  IF(tot_neu_loci>0) Ht_Meta=SUM(obs_locHt_Meta)/tot_neu_loci

!write(*,*) 'ht_meta', Ht_Meta

  DEALLOCATE(vec_sq_freq)
  DEALLOCATE(vec_freq)
  DEALLOCATE(copies)

!----------------------------------
  END SUBROUTINE meta_Ht_neu
!----------------------------------

!----------------------------------
  SUBROUTINE do_neutral_freq
! calculates allele frequencies for neutral loci
!----------------------------------
 INTEGER :: i,j,k
  INTEGER :: allele,pop
  INTEGER, DIMENSION(num_pop_tot,tot_neu_loci,max_neu_allele) :: do_neutral_copies

  FORALL(i=1:num_pop_tot,j=1:tot_neu_loci,k=1:max_neu_allele) do_neutral_copies(i,j,k)=0
  FORALL(i=1:num_pop_tot,j=1:tot_neu_loci,k=1:max_neu_allele) neu_meta_freq(i,j,k)=0.0_dp

  allele=0
  pop=0
! counts # allele copies per allele & per locus & per population
  DO i=1,tot_neu_loci
    DO k=1,num_tot_ind
      pop=indv(k)%pop
      allele=neu_genome(1,i,k)-10
      do_neutral_copies(pop,i,allele)=do_neutral_copies(pop,i,allele)+1
      allele=neu_genome(2,i,k)-10
      do_neutral_copies(pop,i,allele)=do_neutral_copies(pop,i,allele)+1
    END DO
  END DO
!write(*,*) 'copies', do_neutral_copies
! calculates allele frequencies per allele & per locus & per population
  DO i=1,num_pop_tot
    DO j=1,tot_neu_loci
      DO k=1,neu_gene(j)%nalle
        IF(do_neutral_copies(i,j,k)==0)THEN
          neu_meta_freq(i,j,k)=0.0_dp
        ELSE
          neu_meta_freq(i,j,k)=REAL(do_neutral_copies(i,j,k),dp)/SUM(do_neutral_copies(i,j,:))
        END IF
      END DO
    END DO
  END DO
!----------------------------------
  END SUBROUTINE do_neutral_freq
!----------------------------------

!----------------------------------
  SUBROUTINE do_Fst
!----------------------------------
  INTEGER :: i,j,k
  REAL(KIND=dp) :: do_Fst_SQfreq
  REAL(KIND=dp) :: current_Fst
  REAL(KIND=dp), DIMENSION(num_pop_tot) :: do_Fst_Hpop
  REAL(KIND=dp), DIMENSION(num_pop_tot,tot_neu_loci) :: do_Fst_Hloc_pop

  FORALL(i=1:num_pop_tot) do_Fst_Hpop(i)=0.0_dp
  FORALL(i=1:num_pop_tot,j=1:tot_neu_loci) do_Fst_Hloc_pop(i,j)=0.0_dp

! calculates expected heterozygosities per locus & population
  DO i=1,num_pop_tot
    DO j=1,tot_neu_loci
      do_Fst_SQfreq=0.0_dp
      DO k=1,neu_gene(j)%nalle
        do_Fst_SQfreq=do_Fst_SQfreq+neu_meta_freq(i,j,k)*neu_meta_freq(i,j,k)
      END DO
      do_Fst_Hloc_pop(i,j)=1.0_dp-do_Fst_SQfreq
    END DO
    do_Fst_Hpop(i)=SUM(do_Fst_Hloc_pop(i,:))/tot_neu_loci
  END DO

! calculates expected # heterozygotes across populations
! weighted average by population size following Nei (1977)
  Hs_Meta=0.0_dp
  DO i=1,num_pop_tot
    Hs_Meta=Hs_Meta+(REAL(ind_per_pop(i),dp)*do_Fst_Hpop(i))/num_tot_ind
  END DO
! unweighted average
! Hs_Meta=real_mean(do_Fst_Hpop)
  IF(done_replicate==0)THEN
!   calculates Fst
    current_Fst=(Ht_Meta-Hs_Meta)/Ht_Meta
    ini_Fst_Meta=current_Fst
  ELSE IF(done_replicate==1)THEN
    current_Fst=(Ht_Meta-Hs_Meta)/Ht_Meta
    obs_Fst_Meta=current_Fst
!   Fst per locus
    FORALL(i=1:tot_neu_loci) locHs_Meta(i)=0.0_dp
    DO i=1,tot_neu_loci
!     weighted average by population size following Nei (1977)
      DO j=1,num_pop_tot
        locHs_Meta(i)=locHs_Meta(i)+(REAL(ind_per_pop(j),dp)*do_Fst_Hloc_pop(j,i))/num_tot_ind
      END DO
!     unweighted average
!     locHs_Meta(i)=real_mean(do_Fst_Hloc_pop(:,i))
      current_Fst=(obs_locHt_Meta(i)-locHs_Meta(i))/obs_locHt_Meta(i)
      obs_locFst_Meta(i)=current_Fst
    END DO
  END IF
!----------------------------------
  END SUBROUTINE do_Fst
!----------------------------------

!----------------------------------
  SUBROUTINE MonteCarlo_test_HWE
! Monte Carlo routine for Hardy-Weinberg proportions
! generating exact test
!----------------------------------
  INTEGER :: i,j,k
  INTEGER :: allele_1, allele_2,pop,samp_rep
  INTEGER :: rec_allele,max_rec_allele
  INTEGER :: times_larger_chisq
  INTEGER :: num_alleles
  INTEGER, ALLOCATABLE, DIMENSION(:) :: rand_genot_counts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: pop_genotypes
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neu_genot_counts
  REAL(KIND=dp) :: allele_freq_prod
  REAL(KIND=dp) :: rand_chisq
  REAL(KIND=dp) :: sum_rand_chisq
  REAL(KIND=dp) :: sumq_rand_chisq
  REAL(KIND=dp), DIMENSION(num_pop_tot,tot_neu_loci) :: chisq_obs
  REAL(KIND=dp), DIMENSION(3,num_pop_tot,tot_neu_loci) :: chisq_exp
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rand_obs_genot_freq
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: neu_obs_genot_freq
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: neu_exp_genot_freq

  HWE_over_loci: DO i=1,tot_neu_loci
!   calculates number of elements upper triangle + diagonal
!   of matrix alleles * alleles
    num_alleles=neu_gene(i)%nalle
    max_rec_allele=((num_alleles*num_alleles-num_alleles)/2)+num_alleles
    ALLOCATE(neu_genot_counts(num_pop_tot,max_rec_allele))
    FORALL(j=1:num_pop_tot,k=1:max_rec_allele) neu_genot_counts(j,k)=0
    DO k=1,num_tot_ind
!     counts # individuals per genotypic class
      pop=indv(k)%pop
      allele_1=neu_genome(1,i,k)-10	! To avoid confusion between selected and neutral loci we coded neutral loci from 11 to tot_neu_loci+10 --> allele 11 is the first neutral allele
      allele_2=neu_genome(2,i,k)-10	! To avoid confusion between selected and neutral loci we coded neutral loci from 11 to tot_neu_loci+10 --> allele 11 is the first neutral allele
      IF(allele_1<=allele_2) rec_allele=n_rec(allele_1,allele_2,num_alleles)
      IF(allele_2<allele_1) rec_allele=n_rec(allele_2,allele_1,num_alleles)
      neu_genot_counts(pop,rec_allele)=neu_genot_counts(pop,rec_allele)+1
    END DO

!   calculates observed genotypic frequencies
    ALLOCATE(neu_obs_genot_freq(num_pop_tot,max_rec_allele))
    FORALL(j=1:num_pop_tot,k=1:max_rec_allele) neu_obs_genot_freq(j,k)=0.0_dp
    DO j=1,num_pop_tot
      DO k=1,max_rec_allele
        IF(neu_genot_counts(j,k)==0)THEN
          neu_obs_genot_freq(j,k)=0.0_dp
        ELSE
          neu_obs_genot_freq(j,k)=REAL(neu_genot_counts(j,k),dp)/ind_per_pop(j)
        END IF
      END DO
    END DO
    DEALLOCATE(neu_genot_counts)
!   calculates expected genotypic frequencies
    ALLOCATE(neu_exp_genot_freq(num_pop_tot,max_rec_allele))
    FORALL(j=1:num_pop_tot,k=1:max_rec_allele) neu_exp_genot_freq(j,k)=0.0_dp
    allele_freq_prod=0.0_dp
    DO j=1,num_pop_tot
      DO allele_1=1,num_alleles
        DO allele_2=allele_1,num_alleles
          rec_allele=n_rec(allele_1,allele_2,num_alleles)
          IF((neu_meta_freq(j,i,allele_1)==0.0_dp).OR.(neu_meta_freq(j,i,allele_2)==0.0_dp))THEN
            allele_freq_prod=0.0_dp
          ELSE
            allele_freq_prod=neu_meta_freq(j,i,allele_1)*neu_meta_freq(j,i,allele_2)
          END IF
          neu_exp_genot_freq(j,rec_allele)=allele_freq_prod
        END DO
      END DO
    END DO
!   calculates observed Chi squared
    FORALL(j=1:num_pop_tot) chisq_obs(j,i)=0.0_dp
    DO j=1,num_pop_tot
      DO k=1,max_rec_allele
!       following SUM((obs-exp)^2)/exp over genotypic classes
        IF(neu_exp_genot_freq(j,k)>0.0_dp)THEN
          chisq_obs(j,i)=chisq_obs(j,i)+((neu_obs_genot_freq(j,k)-neu_exp_genot_freq(j,k))**2)/neu_exp_genot_freq(j,k)
        END IF
      END DO
    END DO
    DEALLOCATE(neu_obs_genot_freq)
!write(*,*) 'Monte Carlo, neu_obs_genot_freq', neu_obs_genot_freq(1,1)
!write(*,*) 'Monte Carlo, neu_exp_genot_freq', neu_exp_genot_freq(1,1)
!   calculating random allocations of genotypes per population
    HWE_over_pops: DO j=1,num_pop_tot
!     dimension for pop_genotypes
      counter=0
      ALLOCATE(pop_genotypes(ind_per_pop(j)+ind_per_pop(j)))
!     fill the matrix with individuals of that pop
      DO k=1,num_tot_ind
        IF(indv(k)%pop==j)THEN
          counter=counter+1
          pop_genotypes(counter)=neu_genome(1,i,k)
          pop_genotypes(ind_per_pop(j)+counter)=neu_genome(2,i,k)
        END IF
      END DO
      times_larger_chisq=0
      sum_rand_chisq=0.0_dp
      sumq_rand_chisq=0.0_dp
      HWE_sampling_chisq: DO samp_rep=1,tot_HWE_rep
!       allocate alleles anew amongst individuals
        DO k=1,ind_per_pop(j)+ind_per_pop(j)
          xrand=alea(idum)
          rec_allele=FLOOR(xrand*(ind_per_pop(j)+ind_per_pop(j)))+1
          allele_1=pop_genotypes(k)
          pop_genotypes(k)=pop_genotypes(rec_allele)
          pop_genotypes(rec_allele)=allele_1
        END DO
        ALLOCATE(rand_genot_counts(max_rec_allele))
        FORALL(k=1:max_rec_allele) rand_genot_counts(k)=0
        DO k=1,ind_per_pop(j)
!         counts # individuals per genotypic class
          allele_1=pop_genotypes(k)-10	! To avoid confusion between selected and neutral loci we coded neutral loci from 11 to 20 --> allele 11 is the first neutral allele
          allele_2=pop_genotypes(ind_per_pop(j)+k)-10	! To avoid confusion between selected and neutral loci we coded neutral loci from 11 to 20 --> allele 11 is the first neutral allele
          IF(allele_1<=allele_2) rec_allele=n_rec(allele_1,allele_2,num_alleles)
          IF(allele_2<allele_1) rec_allele=n_rec(allele_2,allele_1,num_alleles)
          rand_genot_counts(rec_allele)=rand_genot_counts(rec_allele)+1
        END DO
!       calculates observed genotypic frequencies
        ALLOCATE(rand_obs_genot_freq(max_rec_allele))
        FORALL(k=1:max_rec_allele) rand_obs_genot_freq(k)=0.0_dp
        DO k=1,max_rec_allele
          IF(rand_genot_counts(k)==0)THEN
            rand_obs_genot_freq(k)=0.0_dp
          ELSE
            rand_obs_genot_freq(k)=REAL(rand_genot_counts(k),dp)/ind_per_pop(j)
          END IF
        END DO
        DEALLOCATE(rand_genot_counts)
        rand_chisq=0.0_dp
        DO k=1,max_rec_allele
!         following SUM((obs-exp)^2)/exp over genotypic classes
          IF(neu_exp_genot_freq(j,k)>0.0_dp)THEN
            rand_chisq=rand_chisq+((rand_obs_genot_freq(k)-neu_exp_genot_freq(j,k))**2)/neu_exp_genot_freq(j,k)
          END IF
        END DO
        DEALLOCATE(rand_obs_genot_freq)
        sum_rand_chisq=sum_rand_chisq+rand_chisq
        sumq_rand_chisq=sumq_rand_chisq+rand_chisq*rand_chisq
        IF(rand_chisq<chisq_obs(j,i))THEN
          times_larger_chisq=times_larger_chisq+1
        END IF
      END DO HWE_sampling_chisq
      DEALLOCATE(pop_genotypes)
!     calculates expected Chi squared distribution, it provides:
!     chisq_exp(1,j,i) mean of the distribution of Chi squared
!     chisq_exp(2,j,i) standard deviation of the Chi squared distribution
!     chisq_exp(3,j,i) percentage of times over tot_HWE_rep when rendomized Chi squared values
!                    resulted larger than observed value
      chisq_exp(1,j,i)=sum_rand_chisq/tot_HWE_rep
      chisq_exp(2,j,i)=(sumq_rand_chisq-((sum_rand_chisq**2)/tot_HWE_rep))/(tot_HWE_rep-1)
      IF(times_larger_chisq==0)THEN
        chisq_exp(3,j,i)=0.0_dp
      ELSE
        chisq_exp(3,j,i)=REAL(times_larger_chisq,dp)/tot_HWE_rep
      END IF
    END DO HWE_over_pops
    DEALLOCATE(neu_exp_genot_freq)
  END DO HWE_over_loci
  IF(done_replicate==0)THEN
    DO i=1,tot_neu_loci
      DO j=1,num_pop_tot
        chisq_before_SA(1,j,i)=chisq_obs(j,i)
        chisq_before_SA(2,j,i)=chisq_exp(1,j,i)
        chisq_before_SA(3,j,i)=chisq_exp(2,j,i)
        chisq_before_SA(4,j,i)=chisq_exp(3,j,i)
      END DO
    END DO
  ELSE IF(done_replicate==1)THEN
    DO i=1,tot_neu_loci
      DO j=1,num_pop_tot
        chisq_after_SA(1,j,i)=chisq_obs(j,i)
        chisq_after_SA(2,j,i)=chisq_exp(1,j,i)
        chisq_after_SA(3,j,i)=chisq_exp(2,j,i)
        chisq_after_SA(4,j,i)=chisq_exp(3,j,i)
      END DO
    END DO
  END IF
!----------------------------------
  END SUBROUTINE MonteCarlo_test_HWE
!----------------------------------

!----------------------------------
  SUBROUTINE do_Qst
! allows only biallelic loci & additive effects;
! assuming VarBetween/(VarBetween+2VarWithin)
! following Spitze in Genetics 135: 367-374.
! Individual values are assumed to be
! the true genotypic values
!----------------------------------
  INTEGER :: i_Qst,j_Qst
  REAL(KIND=dp) :: weight_pop_size
  REAL(KIND=dp), DIMENSION(tot_sel_loci,num_tot_ind) :: locus_genot
  REAL(KIND=dp), DIMENSION(num_traits,num_tot_ind) :: indv_genot
  REAL(KIND=dp), DIMENSION(num_traits,num_pop_tot) :: mean_trait_pop
  REAL(KIND=dp), DIMENSION(num_traits) :: mean_trait_tot
  REAL(KIND=dp), DIMENSION(num_traits) :: SSamong
  REAL(KIND=dp), DIMENSION(num_traits) :: SSwithin
  REAL(KIND=dp), DIMENSION(num_traits) :: MSamong
  REAL(KIND=dp), DIMENSION(num_traits) :: MSwithin
  REAL(KIND=dp), DIMENSION(num_traits) :: expMSamong
! variables used if done_replicate == 1
  REAL(KIND=dp) :: sum_tot_genot
  REAL(KIND=dp) :: sumQ_tot_genot
  REAL(KIND=dp), DIMENSION(num_traits) :: sum_Qst
  REAL(KIND=dp), DIMENSION(num_traits) :: sumQ_Qst
  REAL(KIND=dp), DIMENSION(num_pop_tot) :: sum_genot
  REAL(KIND=dp), DIMENSION(num_pop_tot) :: sumQ_genot
  REAL(KIND=dp), DIMENSION(tot_sel_loci,num_pop_tot) :: mean_trait_pop_loc
  REAL(KIND=dp), DIMENSION(tot_sel_loci) :: mean_trait_tot_loc
  REAL(KIND=dp), DIMENSION(tot_sel_loci) :: SSamong_loc
  REAL(KIND=dp), DIMENSION(tot_sel_loci) :: SSwithin_loc
  REAL(KIND=dp), DIMENSION(tot_sel_loci) :: MSamong_loc
  REAL(KIND=dp), DIMENSION(tot_sel_loci) :: MSwithin_loc
  REAL(KIND=dp), DIMENSION(tot_sel_loci) :: expMSamong_loc

  INTEGER, ALLOCATABLE, DIMENSION(:) :: count    

	ALLOCATE(count(num_traits))

! genotypic values per locus, according to the genotypic scale per locus
  FORALL(i_Qst=1:tot_sel_loci,j_Qst=1:num_tot_ind) locus_genot(i_Qst,j_Qst)=0.0_dp
  DO i_Qst=1,tot_sel_loci
    DO j_Qst=1,num_tot_ind
      IF(sel_genome(1,i_Qst,j_Qst)+sel_genome(2,i_Qst,j_Qst)==2)THEN
        locus_genot(i_Qst,j_Qst)=-genotypic_scale(i_Qst,1)
      ELSE IF(sel_genome(1,i_Qst,j_Qst)+sel_genome(2,i_Qst,j_Qst)==3)THEN
        locus_genot(i_Qst,j_Qst)=genotypic_scale(i_Qst,2)
      ELSE IF(sel_genome(1,i_Qst,j_Qst)+sel_genome(2,i_Qst,j_Qst)==4)THEN
        locus_genot(i_Qst,j_Qst)=genotypic_scale(i_Qst,1)
      END IF
    END DO
  END DO

!write(*,*) 'do QST, sel_genome, indiv 1, al1', sel_genome(1,:,1)
!write(*,*) 'do QST, sel_genome, indiv 1, al2', sel_genome(2,:,1)
!write(*,*) 'do QST, sel_genome, indiv 2, al1', sel_genome(1,:,2)
!write(*,*) 'do QST, sel_genome, indiv 2, al2', sel_genome(2,:,2)
!write(*,*) 'genotypic_scale', genotypic_scale
!write(*,*) 'do QST, locus_genot, indiv 1', locus_genot(:,1)
!write(*,*) 'do QST, locus_genot, indiv 2', locus_genot(:,2)

! genotypic values per individual: sum across loci for same trait
  FORALL(i_Qst=1:num_traits,j_Qst=1:num_tot_ind) indv_genot(i_Qst,j_Qst)=0.0_dp
  DO i_Qst=1,num_tot_ind
		count(1)=0
		count(2) = trait_sel_loci(1)
!write(*,*) counter
    DO j_Qst=1,tot_sel_loci
!write(*,*) counter(sel_gene(j_Qst)%funct)
				count(sel_gene(j_Qst)%funct) = count(sel_gene(j_Qst)%funct)+1
!write(*,*) counter
      indv_genot(sel_gene(j_Qst)%funct,i_Qst)=&
      & indv_genot(sel_gene(j_Qst)%funct,i_Qst)+locus_genot(count(sel_gene(j_Qst)%funct),i_Qst)
    END DO
  END DO
!write(*,*) 'do QST, indv_genot, indiv 1', indv_genot(1,1)
!write(*,*) 'do QST, indv_genot, indiv 2', indv_genot(2,1)


DEALLOCATE(count)

! means per population & global means (per trait)
  FORALL(i_Qst=1:num_traits,j_Qst=1:num_pop_tot) mean_trait_pop(i_Qst,j_Qst)=0.0_dp
  FORALL(i_Qst=1:num_traits) mean_trait_tot(i_Qst)=0.0_dp
  DO i_Qst=1,num_traits
    DO j_Qst=1,num_tot_ind
      mean_trait_pop(i_Qst,indv(j_Qst)%pop)=&
      & mean_trait_pop(i_Qst,indv(j_Qst)%pop)+indv_genot(i_Qst,j_Qst)
    END DO
  END DO
  DO i_Qst=1,num_traits
    mean_trait_tot(i_Qst)=SUM(mean_trait_pop(i_Qst,:))/num_tot_ind
    DO j_Qst=1,num_pop_tot
      mean_trait_pop(i_Qst,j_Qst)=mean_trait_pop(i_Qst,j_Qst)/ind_per_pop(j_Qst)
    END DO
  END DO

!write(*,*) 'do QST, num_traits', num_traits
!write(*,*) 'do QST, mean_trait_pop, trait1', mean_trait_pop(1,:)
!write(*,*) 'do QST, mean_trait_pop, trait2', mean_trait_pop(2,:)
!write(*,*) 'do QST, mean_trait_tot', mean_trait_tot
! SSamong & MSamong (per trait)
  FORALL(i_Qst=1:num_traits) SSamong(i_Qst)=0.0_dp
  FORALL(i_Qst=1:num_traits) MSamong(i_Qst)=0.0_dp
  DO i_Qst=1,num_traits
    DO j_Qst=1,num_pop_tot
      SSamong(i_Qst)=SSamong(i_Qst)+&
      & ind_per_pop(j_Qst)*(mean_trait_pop(i_Qst,j_Qst)-mean_trait_tot(i_Qst))**2
    END DO
    MSamong(i_Qst)=SSamong(i_Qst)/(num_pop_tot-1)
  END DO

!write(*,*) 'do QST, MSamong', MSamong

! SSwithin & MSwithin (per trait)
  FORALL(i_Qst=1:num_traits) SSwithin(i_Qst)=0.0_dp
  FORALL(i_Qst=1:num_traits) MSwithin(i_Qst)=0.0_dp
!  DO i_Qst=1,num_traits
!    DO j_Qst=1,num_tot_ind
!indv_genot(i_Qst,j_Qst)=(indv_genot(i_Qst,j_Qst)-mean_trait_pop(i_Qst,indv(j_Qst)%pop))**2
!    END DO
!  END DO
  DO i_Qst=1,num_traits
    DO j_Qst=1,num_tot_ind
!      SSwithin(i_Qst)=SSwithin(i_Qst)+indv_genot(i_Qst,j_Qst)
      SSwithin(i_Qst)=SSwithin(i_Qst)+(indv_genot(i_Qst,j_Qst)-mean_trait_pop(i_Qst,indv(j_Qst)%pop))**2
    END DO
    MSwithin(i_Qst)=SSwithin(i_Qst)/(num_tot_ind-num_pop_tot)
  END DO

!write(*,*) 'do QST, MSwithin', MSwithin
! weighted population size if unequal
  weight_pop_size=0.0_dp
  IF(num_pop_tot>1)THEN
    weight_pop_size=&
    & (1.0_dp/(num_pop_tot-1))*(REAL(num_tot_ind,dp)-(REAL(SUM(ind_per_pop**2),dp)/num_tot_ind))
  ELSE IF(num_pop_tot==1)THEN
    weight_pop_size=REAL(num_tot_ind,dp)
  END IF
!write(*,*) 'do QST, weight_pop_size', weight_pop_size
! expected MSamong
  FORALL(i_Qst=1:num_traits) expMSamong(i_Qst)=0.0_dp
  FORALL(i_Qst=1:num_traits) expMSamong(i_Qst)=(MSamong(i_Qst)-MSwithin(i_Qst))/weight_pop_size

  IF(done_replicate==0)THEN
!   Qst before SA (per trait)
    FORALL(i_Qst=1:num_traits) ini_Qst_Meta(i_Qst)=0.0_dp
    DO i_Qst=1,num_traits
      ini_Qst_Meta(i_Qst)=expMSamong(i_Qst)/(expMSamong(i_Qst)+2.0_dp*MSwithin(i_Qst))
    END DO
  ELSE IF(done_replicate==1)THEN
!   Qst after SA (per trait)
    FORALL(i_Qst=1:num_traits) obs_Qst_Meta(i_Qst)=0.0_dp
    DO i_Qst=1,num_traits
      obs_Qst_Meta(i_Qst)=expMSamong(i_Qst)/(expMSamong(i_Qst)+2.0_dp*MSwithin(i_Qst))
!write(*,*) 'do QST, SSamong', SSamong(i_Qst)
!write(*,*) 'do QST, MSwithin', MSwithin(i_Qst)
!write(*,*) 'do QST, expMSamong', expMSamong(i_Qst)
!write(*,*) 'do Qst obs', obs_Qst_Meta
    END DO
!   means per population & global means (per locus)
    FORALL(i_Qst=1:tot_sel_loci,j_Qst=1:num_pop_tot) mean_trait_pop_loc(i_Qst,j_Qst)=0.0_dp
    FORALL(i_Qst=1:tot_sel_loci) mean_trait_tot_loc(i_Qst)=0.0_dp
    DO i_Qst=1,tot_sel_loci
      DO j_Qst=1,num_tot_ind
        counter=indv(j_Qst)%pop
        mean_trait_pop_loc(i_Qst,counter)=mean_trait_pop_loc(i_Qst,counter)+locus_genot(i_Qst,j_Qst)
      END DO
    END DO
    DO i_Qst=1,tot_sel_loci
      mean_trait_tot_loc(i_Qst)=SUM(mean_trait_pop_loc(i_Qst,:))/num_tot_ind
      DO j_Qst=1,num_pop_tot
        mean_trait_pop_loc(i_Qst,j_Qst)=mean_trait_pop_loc(i_Qst,j_Qst)/ind_per_pop(j_Qst)
      END DO
    END DO

!write(*,*) 'do Qst mean_trait_pop_loc', mean_trait_pop_loc
!   SSamong & MSamong (per locus)
    FORALL(i_Qst=1:tot_sel_loci) SSamong_loc(i_Qst)=0.0_dp
    FORALL(i_Qst=1:tot_sel_loci) MSamong_loc(i_Qst)=0.0_dp
    DO i_Qst=1,tot_sel_loci
      DO j_Qst=1,num_pop_tot
        SSamong_loc(i_Qst)=SSamong_loc(i_Qst)+&
        & ind_per_pop(j_Qst)*(mean_trait_pop_loc(i_Qst,j_Qst)-mean_trait_tot_loc(i_Qst))**2
      END DO
      MSamong_loc(i_Qst)=SSamong_loc(i_Qst)/(num_pop_tot-1)
!write(*,*) 'do Qst MSamong_loc', MSamong_loc
    END DO
!   SSwithin & MSwithin (per locus)
    FORALL(i_Qst=1:tot_sel_loci) SSwithin_loc(i_Qst)=0.0_dp
    FORALL(i_Qst=1:tot_sel_loci) MSwithin_loc(i_Qst)=0.0_dp
    DO i_Qst=1,tot_sel_loci
      DO j_Qst=1,num_tot_ind
        locus_genot(i_Qst,j_Qst)=&
        & (locus_genot(i_Qst,j_Qst)-mean_trait_pop_loc(i_Qst,indv(j_Qst)%pop))**2
      END DO
    END DO
    DO i_Qst=1,tot_sel_loci
      DO j_Qst=1,num_tot_ind
        SSwithin_loc(i_Qst)=SSwithin_loc(i_Qst)+locus_genot(i_Qst,j_Qst)
      END DO
      MSwithin_loc(i_Qst)=SSwithin_loc(i_Qst)/(num_tot_ind-num_pop_tot)
    END DO
!write(*,*) 'do Qst MSwithin_loc', MSwithin_loc
!   expected MSamong (per locus)
    FORALL(i_Qst=1:tot_sel_loci) expMSamong_loc(i_Qst)=0.0_dp
    FORALL(i_Qst=1:tot_sel_loci) expMSamong_loc(i_Qst)=&
    &(MSamong_loc(i_Qst)-MSwithin_loc(i_Qst))/weight_pop_size
!   Qst (per locus)
    FORALL(i_Qst=1:num_traits) sum_Qst(i_Qst)=0.0_dp
    FORALL(i_Qst=1:num_traits) sumQ_Qst(i_Qst)=0.0_dp
    DO i_Qst=1,tot_sel_loci
      obs_locQst_Meta(i_Qst)=expMSamong_loc(i_Qst)/(expMSamong_loc(i_Qst)+2.0_dp*MSwithin_loc(i_Qst))
    END DO
!write(*,*) 'do Qst obs_locQst_Meta', obs_locQst_Meta
!   genotypic variances per trait et per trait*pop
!   environmental variances per trait
    DO i_Qst=1,num_traits
      FORALL(i_Qst=1:num_pop_tot) sum_genot(i_Qst)=0.0_dp
      FORALL(i_Qst=1:num_pop_tot) sumQ_genot(i_Qst)=0.0_dp
      sum_tot_genot=0.0_dp
      sumQ_tot_genot=0.0_dp
      DO j_Qst=1,num_tot_ind
        sum_genot(indv(j_Qst)%pop)=sum_genot(indv(j_Qst)%pop)+indv_genot(i_Qst,j_Qst)
        sumQ_genot(indv(j_Qst)%pop)=sumQ_genot(indv(j_Qst)%pop)+&
        & indv_genot(i_Qst,j_Qst)*indv_genot(i_Qst,j_Qst)
        sum_tot_genot=sum_tot_genot+indv_genot(i_Qst,j_Qst)
        sumQ_tot_genot=sumQ_tot_genot+indv_genot(i_Qst,j_Qst)*indv_genot(i_Qst,j_Qst)
      END DO
      FORALL(j_Qst=1:num_pop_tot) VGenot_pop(j_Qst,i_Qst)=&
      & (sumQ_genot(j_Qst)-((sum_genot(j_Qst)*sum_genot(j_Qst))/ind_per_pop(j_Qst)))/&
      & (ind_per_pop(j_Qst)-1)
      VGenot_tot(i_Qst)=(sumQ_tot_genot-((sum_tot_genot*sum_tot_genot)/num_tot_ind))/(num_tot_ind-1)

!write(*,*) 'VGenot_tot', VGenot_tot(i_Qst)

    END DO
  END IF
!----------------------------------
  END SUBROUTINE do_Qst
!----------------------------------

!----------------------------------
  SUBROUTINE do_linkage
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j,k,m
  INTEGER :: trait_i,trait_j
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rvec_1
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rvec_2
  REAL(KIND=dp), DIMENSION(tot_sel_loci,num_tot_ind) :: locus_genot

! genotypic values per locus, according to the genotypic scale per locus
  FORALL(i=1:tot_sel_loci,j=1:num_tot_ind) locus_genot(i,j)=0.0_dp
  DO i=1,tot_sel_loci
    DO j=1,num_tot_ind
      IF(sel_genome(1,i,j)+sel_genome(2,i,j)==2)THEN
        locus_genot(i,j)=-genotypic_scale(i,1)
      ELSE IF(sel_genome(1,i,j)+sel_genome(2,i,j)==3)THEN
        locus_genot(i,j)=genotypic_scale(i,2)
      ELSE IF(sel_genome(1,i,j)+sel_genome(2,i,j)==4)THEN
        locus_genot(i,j)=genotypic_scale(i,1)
      END IF
    END DO
  END DO

! genotypic covariances between loci per population
  IF(done_replicate==0)THEN
    FORALL(i=1:num_pop_tot,j=1:num_traits,k=1:num_traits) ini_link_trait(i,j,k)=0.0_dp
    DO i=1,tot_sel_loci-1
      DO j=i+1,tot_sel_loci
        trait_i=sel_gene(i)%funct
        trait_j=sel_gene(j)%funct
        DO k=1,num_pop_tot
          ALLOCATE(rvec_1(ind_per_pop(k)),rvec_2(ind_per_pop(k)))
          FORALL(m=1:ind_per_pop(k)) rvec_1(m)=0.0_dp
          FORALL(m=1:ind_per_pop(k)) rvec_2(m)=0.0_dp
          counter=0
          DO m=1,num_tot_ind
            IF(indv(m)%pop==k)THEN
              counter=counter+1
              rvec_1(counter)=locus_genot(i,m)
              rvec_2(counter)=locus_genot(j,m)
            END IF
          END DO
          ini_link_trait(k,trait_i,trait_j)=ini_link_trait(k,trait_i,trait_j)+real_cov(rvec_1,rvec_2)
          ini_link_trait(k,trait_j,trait_i)=ini_link_trait(k,trait_j,trait_i)+real_cov(rvec_2,rvec_1)
          DEALLOCATE(rvec_1,rvec_2)
        END DO
      END DO
    END DO
  ELSE IF(done_replicate==1)THEN
    FORALL(i=1:num_pop_tot,j=1:num_traits,k=1:num_traits) end_link_trait(i,j,k)=0.0_dp
    DO i=1,tot_sel_loci-1
      DO j=i+1,tot_sel_loci
        trait_i=sel_gene(i)%funct
        trait_j=sel_gene(j)%funct
        DO k=1,num_pop_tot
          ALLOCATE(rvec_1(ind_per_pop(k)),rvec_2(ind_per_pop(k)))
          FORALL(m=1:ind_per_pop(k)) rvec_1(m)=0.0_dp
          FORALL(m=1:ind_per_pop(k)) rvec_2(m)=0.0_dp
          counter=0
          DO m=1,num_tot_ind
            IF(indv(m)%pop==k)THEN
              counter=counter+1
              rvec_1(counter)=locus_genot(i,m)
              rvec_2(counter)=locus_genot(j,m)
            END IF
          END DO
          end_link_trait(k,trait_i,trait_j)=end_link_trait(k,trait_i,trait_j)+real_cov(rvec_1,rvec_2)
          end_link_trait(k,trait_j,trait_i)=end_link_trait(k,trait_j,trait_i)+real_cov(rvec_2,rvec_1)
          DEALLOCATE(rvec_1,rvec_2)
        END DO
      END DO
    END DO
  END IF
!----------------------------------
  END SUBROUTINE do_linkage
!----------------------------------

!----------------------------------
  SUBROUTINE SA_Fst
! LSR 20/07/04
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER :: nsucc_sa
  INTEGER :: sampled_locus
  REAL(KIND=dp) :: fbest_sa
  REAL(KIND=dp) :: score_sa
  REAL(KIND=dp) :: delta_sa
  REAL(KIND=dp) :: tempo_sa
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_trans    ! test locus backup for transfer
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: best_sa ! best genome backup for SA
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: test_sa ! test genome backup for SA

  ALLOCATE(test_trans(2,num_tot_ind))
  ALLOCATE(best_sa(2,tot_neu_loci,num_tot_ind))
  ALLOCATE(test_sa(2,tot_neu_loci,num_tot_ind))

  FORALL(i=1:2,j=1:num_tot_ind) test_trans(i,j)=0
  FORALL(i=1:2,j=1:tot_neu_loci,k=1:num_tot_ind) best_sa(i,j,k)=0
  FORALL(i=1:2,j=1:tot_neu_loci,k=1:num_tot_ind) test_sa(i,j,k)=0

! copies neutral loci genome into SA matrices
  FORALL(i=1:2,j=1:tot_neu_loci,k=1:num_tot_ind) best_sa(i,j,k)=neu_genome(i,j,k)

! SIMULATED ANNEALING ALGORITHM
! (Pham & Karaboga, 2000)
    tempo_sa=Tini_Fst
    fbest_sa=Fst_funct(best_sa)

    across_Fst_temp: DO i=1,Tsteps_Fst
      nsucc_sa=0
      within_Fst_temp: DO j=1,Nover_Fst
        test_sa=best_sa
!       picks one locus for Fst_trans
        xrand=alea(idum)
        sampled_locus=FLOOR(xrand*tot_neu_loci)+1
        FORALL(k=1:num_tot_ind) test_trans(1,k)=test_sa(1,sampled_locus,k)
        FORALL(k=1:num_tot_ind) test_trans(2,k)=test_sa(2,sampled_locus,k)
        CALL Fst_trans(test_trans)
        FORALL(k=1:num_tot_ind) test_sa(1,sampled_locus,k)=test_trans(1,k)
        FORALL(k=1:num_tot_ind) test_sa(2,sampled_locus,k)=test_trans(2,k)
!       scores new arrangement
        score_sa=Fst_funct(test_sa)
        delta_sa=score_sa-fbest_sa
        xrand=alea(idum)
        IF((delta_sa<=0.0d0).OR.(xrand<exp(-delta_sa/tempo_sa)))THEN
          nsucc_sa=nsucc_sa+1
          fbest_sa=score_sa
          best_sa=test_sa
        END IF
        IF(nsucc_sa>=Nlimit_Fst) EXIT within_Fst_temp
      END DO within_Fst_temp
      IF(nsucc_sa==0) EXIT across_Fst_temp
      tempo_sa=tempo_sa*Tfactr_Fst
    END DO across_Fst_temp

! copies SA matrices into neutral loci genome
  FORALL(i=1:2,j=1:tot_neu_loci,k=1:num_tot_ind) neu_genome(i,j,k)=best_sa(i,j,k)

  DEALLOCATE(best_sa)
  DEALLOCATE(test_sa)
  DEALLOCATE(test_trans)

  CONTAINS
!----------------------------------
    SUBROUTINE Fst_trans(Fst_trans_mat)
!   swaps alleles in between individuals
!   coming from different populations
!----------------------------------
    INTEGER, DIMENSION(2,num_tot_ind), INTENT(INOUT) :: Fst_trans_mat
    INTEGER :: swap1,swap1_pop,swap1_haplot,ind_pop1
    INTEGER :: swap2,swap2_pop,swap2_haplot,ind_pop2
    INTEGER :: buffer
		INTEGER :: progress_counter_2

!   picks a 1st population
    xrand=alea(idum)
    swap1_pop=FLOOR(xrand*num_pop_tot)+1
    swap2_pop=swap1_pop
    ind_pop1=ind_per_pop(swap1_pop)
!   picks a 2nd population
		progress_counter_2 = 0
    DO WHILE(swap1_pop==swap2_pop)

			progress_counter_2 = progress_counter_2+1
			IF (MODULO(progress_counter_2, 100)==0) WRITE(*,*) '[SA_Fst] swaping alleles to reach target Fst, &
					iteration: ' , progress_counter_2

      xrand=alea(idum)
      swap2_pop=FLOOR(xrand*num_pop_tot)+1
    END DO
    ind_pop2=ind_per_pop(swap2_pop)
!   picks a 1st candidate & haplotype within 1st population
    xrand=alea(idum)
    swap1=posit_ind_pop(swap1_pop,(FLOOR(xrand*ind_pop1)+1))
    xrand=alea(idum)
    swap1_haplot=FLOOR(xrand*2)+1
!   picks a 2nd candidate & haplotype within 2nd population
    xrand=alea(idum)
    swap2=posit_ind_pop(swap2_pop,(FLOOR(xrand*ind_pop2)+1))
    xrand=alea(idum)
    swap2_haplot=FLOOR(xrand*2)+1
!   swaps alleles
    buffer=Fst_trans_mat(swap1_haplot,swap1)
    Fst_trans_mat(swap1_haplot,swap1)=Fst_trans_mat(swap2_haplot,swap2)
    Fst_trans_mat(swap2_haplot,swap2)=buffer
!----------------------------------
    END SUBROUTINE Fst_trans
!----------------------------------
!----------------------------------
    FUNCTION Fst_funct(Fst_funct_mat)
!----------------------------------
    INTEGER :: i_Fst,j_Fst,k_Fst
    INTEGER :: allele,pop
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: Fst_funct_mat
    INTEGER, DIMENSION(num_pop_tot,tot_neu_loci,max_neu_allele) :: Fst_funct_copies
    REAL(KIND=dp) :: Fst_funct
    REAL(KIND=dp) :: Fst_funct_expHet
    REAL(KIND=dp) :: Fst_funct_SQfreq
    REAL(KIND=dp) :: Fst_funct_value
    REAL(KIND=dp), DIMENSION(num_pop_tot) :: Fst_funct_Hpop
    REAL(KIND=dp), DIMENSION(num_pop_tot,tot_neu_loci) :: Fst_funct_Het
    REAL(KIND=dp), DIMENSION(num_pop_tot,tot_neu_loci,max_neu_allele) :: Fst_funct_freqs

!   to save computing time make Fst_matrices global in metatrom_v7_param
!   and allocation in the main body of the program...

    FORALL(i_Fst=1:num_pop_tot,j_Fst=1:tot_neu_loci,k_Fst=1:max_neu_allele)&
    & Fst_funct_copies(i_Fst,j_Fst,k_Fst)=0
    FORALL(i_Fst=1:num_pop_tot) Fst_funct_Hpop(i_Fst)=0.0_dp
    FORALL(i_Fst=1:num_pop_tot,j_Fst=1:tot_neu_loci) Fst_funct_Het(i_Fst,j_Fst)=0.0_dp
    FORALL(i_Fst=1:num_pop_tot,j_Fst=1:tot_neu_loci,k_Fst=1:max_neu_allele)&
    & Fst_funct_freqs(i_Fst,j_Fst,k_Fst)=0.0_dp
    allele=0
    pop=0

!   counts # allele copies per locus & population
    DO i_Fst=1,tot_neu_loci
      DO j_Fst=1,2
        DO k_Fst=1,num_tot_ind
          allele=Fst_funct_mat(j_Fst,i_Fst,k_Fst)-10	! To avoid confusion between selectd and neutral loci we coded neutral loci from 11 to 20 --> allele 11 is the first neutral allele
          pop=indv(k_Fst)%pop
          Fst_funct_copies(pop,i_Fst,allele)=Fst_funct_copies(pop,i_Fst,allele)+1
        END DO
      END DO
    END DO
!   calculates allele frequencies per locus & population
    DO i_Fst=1,num_pop_tot
      DO j_Fst=1,tot_neu_loci
        DO k_Fst=1,max_neu_allele
          Fst_funct_freqs(i_Fst,j_Fst,k_Fst)=&
          & REAL(Fst_funct_copies(i_Fst,j_Fst,k_Fst),dp)/SUM(Fst_funct_copies(i_Fst,j_Fst,:))
        END DO
      END DO
    END DO
!   calculates expected heterozygosities per locus & population
    DO i_Fst=1,num_pop_tot
      DO j_Fst=1,tot_neu_loci
        Fst_funct_SQfreq=0.0_dp
        DO k_Fst=1,max_neu_allele
          Fst_funct_SQfreq=Fst_funct_SQfreq+Fst_funct_freqs(i_Fst,j_Fst,k_Fst)*Fst_funct_freqs(i_Fst,j_Fst,k_Fst)
        END DO
        Fst_funct_Het(i_Fst,j_Fst)=1.0_dp-Fst_funct_SQfreq
      END DO
      Fst_funct_Hpop(i_Fst)=SUM(Fst_funct_Het(i_Fst,:))/tot_neu_loci
    END DO
!   calculates expected # heterozygotes across populations
    Fst_funct_expHet=0.0_dp
!   weighted average by population size following Nei (1977)
    DO i_Fst=1,num_pop_tot
      Fst_funct_expHet=Fst_funct_expHet+(REAL(ind_per_pop(i_Fst),dp)*Fst_funct_Hpop(i_Fst))/num_tot_ind
    END DO
!   unweighted average
!   Fst_funct_expHet=real_mean(Fst_funct_Hpop)
    Fst_funct_value=0.0_dp
    Fst_funct_value=(Ht_Meta-Fst_funct_expHet)/Ht_Meta

!   SOLUTION: the squared distance between observed and expected values
    Fst_funct=0.0_dp
    Fst_funct=(exp_Fst-Fst_funct_value)*(exp_Fst-Fst_funct_value)
!----------------------------------
    END FUNCTION Fst_funct
!----------------------------------
!----------------------------------
  END SUBROUTINE SA_Fst
!----------------------------------

!----------------------------------
  SUBROUTINE SA_loc_Fst
! LSR 10/08/05
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j
  INTEGER :: sa_locus
  INTEGER :: nsucc_sa
  INTEGER :: sampled_locus
  REAL(KIND=dp) :: fbest_sa
  REAL(KIND=dp) :: score_sa
  REAL(KIND=dp) :: delta_sa
  REAL(KIND=dp) :: tempo_sa
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: best_sa   ! best genome backup for SA
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_sa   ! test genome backup for SA

  ALLOCATE(best_sa(2,num_tot_ind))
  ALLOCATE(test_sa(2,num_tot_ind))
  FORALL(i=1:2,j=1:num_tot_ind) best_sa(i,j)=0
  FORALL(i=1:2,j=1:num_tot_ind) test_sa(i,j)=0

  DO sa_locus=1,tot_neu_loci
!   copies neutral loci genome into SA matrices
    FORALL(i=1:2,j=1:num_tot_ind) best_sa(i,j)=neu_genome(i,sa_locus,j)
!   SIMULATED ANNEALING ALGORITHM
!   (Pham & Karaboga, 2000)
      tempo_sa=Tini_locFst
      fbest_sa=Fst_loc_funct(sa_locus,best_sa)
      across_Fst_loc_temp: DO i=1,Tsteps_locFst
        nsucc_sa=0
        within_Fst_loc_temp: DO j=1,Nover_locFst
          test_sa=best_sa
          CALL Fst_loc_trans(test_sa)
!         scores new arrangement
          score_sa=Fst_loc_funct(sa_locus,test_sa)
          delta_sa=score_sa-fbest_sa
          xrand=alea(idum)
          IF((delta_sa<=0.0d0).OR.(xrand<exp(-delta_sa/tempo_sa)))THEN
            nsucc_sa=nsucc_sa+1
            fbest_sa=score_sa
            best_sa=test_sa
          END IF
          IF(nsucc_sa>=Nlimit_locFst) EXIT within_Fst_loc_temp
        END DO within_Fst_loc_temp
        IF(nsucc_sa==0) EXIT across_Fst_loc_temp
        tempo_sa=tempo_sa*Tfactr_locFst
      END DO across_Fst_loc_temp
!   copies SA matrices into neutral loci genome
    FORALL(i=1:2,j=1:num_tot_ind) neu_genome(i,sa_locus,j)=best_sa(i,j)
  END DO ! sa_locus=1,tot_neu_loci

  DEALLOCATE(best_sa)
  DEALLOCATE(test_sa)

  CONTAINS
!----------------------------------
    SUBROUTINE Fst_loc_trans(Fst_trans_mat)
!   swaps alleles in between individuals
!   coming from different populations
!----------------------------------
    INTEGER, DIMENSION(2,num_tot_ind), INTENT(INOUT) :: Fst_trans_mat
    INTEGER :: swap1,swap1_pop,swap1_haplot,ind_pop1
    INTEGER :: swap2,swap2_pop,swap2_haplot,ind_pop2
    INTEGER :: buffer
    INTEGER :: progress_counter_3

!   picks a 1st population
    xrand=alea(idum)
    swap1_pop=FLOOR(xrand*num_pop_tot)+1
    swap2_pop=swap1_pop
    ind_pop1=ind_per_pop(swap1_pop)
!   picks a 2nd population
		progress_counter_3=0
    DO WHILE(swap1_pop==swap2_pop)

			progress_counter_3 = progress_counter_3+1
			IF (MODULO(progress_counter_3, 100)==0) WRITE(*,*) '[SA_loc_Fst] swaping alleles to reach target Fst per locus, &
					iteration: ', progress_counter_3

      xrand=alea(idum)
      swap2_pop=FLOOR(xrand*num_pop_tot)+1
    END DO
    ind_pop2=ind_per_pop(swap2_pop)
!   picks a 1st candidate & haplotype within 1st population
    xrand=alea(idum)
    swap1=posit_ind_pop(swap1_pop,(FLOOR(xrand*ind_pop1)+1))
    xrand=alea(idum)
    swap1_haplot=FLOOR(xrand*2)+1
!   picks a 2nd candidate & haplotype within 2nd population
    xrand=alea(idum)
    swap2=posit_ind_pop(swap2_pop,(FLOOR(xrand*ind_pop2)+1))
    xrand=alea(idum)
    swap2_haplot=FLOOR(xrand*2)+1
!   swaps alleles
    buffer=Fst_trans_mat(swap1_haplot,swap1)
    Fst_trans_mat(swap1_haplot,swap1)=Fst_trans_mat(swap2_haplot,swap2)
    Fst_trans_mat(swap2_haplot,swap2)=buffer
!----------------------------------
    END SUBROUTINE Fst_loc_trans
!----------------------------------
!----------------------------------
    FUNCTION Fst_loc_funct(Fst_funct_loc,Fst_funct_mat)
!----------------------------------
    INTEGER :: i_Fst,j_Fst,k_Fst
    INTEGER :: allele,pop
    INTEGER, INTENT(IN) :: Fst_funct_loc
    INTEGER, DIMENSION(:,:), INTENT(IN) :: Fst_funct_mat
    INTEGER, DIMENSION(num_pop_tot,max_neu_allele) :: Fst_funct_copies
    REAL(KIND=dp) :: Fst_loc_funct
    REAL(KIND=dp) :: Fst_funct_expHet
    REAL(KIND=dp) :: Fst_funct_SQfreq
    REAL(KIND=dp) :: Fst_funct_value
    REAL(KIND=dp), DIMENSION(num_pop_tot) :: Fst_funct_Het
    REAL(KIND=dp), DIMENSION(num_pop_tot,max_neu_allele) :: Fst_funct_freqs

    FORALL(i_Fst=1:num_pop_tot,j_Fst=1:max_neu_allele) Fst_funct_copies(i_Fst,j_Fst)=0
    FORALL(i_Fst=1:num_pop_tot) Fst_funct_Het(i_Fst)=0.0_dp
    FORALL(i_Fst=1:num_pop_tot,j_Fst=1:max_neu_allele) Fst_funct_freqs(i_Fst,j_Fst)=0.0_dp
    allele=0
    pop=0

!   counts # allele copies per locus & population
    DO j_Fst=1,2
      DO k_Fst=1,num_tot_ind
        allele=Fst_funct_mat(j_Fst,k_Fst)
        pop=indv(k_Fst)%pop
        Fst_funct_copies(pop,allele)=Fst_funct_copies(pop,allele)+1
      END DO
    END DO
!   calculates allele frequencies per locus & population
    DO i_Fst=1,num_pop_tot
      DO k_Fst=1,max_neu_allele
        Fst_funct_freqs(i_Fst,k_Fst)=REAL(Fst_funct_copies(i_Fst,k_Fst),dp)/SUM(Fst_funct_copies(i_Fst,:))
      END DO
    END DO
!   calculates expected heterozygosities per locus & population
    DO i_Fst=1,num_pop_tot
      Fst_funct_SQfreq=0.0_dp
      DO k_Fst=1,max_neu_allele
        Fst_funct_SQfreq=Fst_funct_SQfreq+Fst_funct_freqs(i_Fst,k_Fst)*Fst_funct_freqs(i_Fst,k_Fst)
      END DO
      Fst_funct_Het(i_Fst)=1.0_dp-Fst_funct_SQfreq
    END DO
!   calculates expected # heterozygotes across populations
    Fst_funct_expHet=0.0_dp
!   weighted average by population size following Nei (1977)
    DO i_Fst=1,num_pop_tot
      Fst_funct_expHet=Fst_funct_expHet+(REAL(ind_per_pop(i_Fst),dp)*Fst_funct_Het(i_Fst))/num_tot_ind
    END DO
!   unweighted average
    Fst_funct_value=0.0_dp
    Fst_funct_value=(obs_locHt_Meta(Fst_funct_loc)-Fst_funct_expHet)/obs_locHt_Meta(Fst_funct_loc)
!   SOLUTION: the squared distance between observed and expected values
    Fst_loc_funct=0.0_dp
    Fst_loc_funct=(neu_gene(Fst_funct_loc)%locFst-Fst_funct_value)*(neu_gene(Fst_funct_loc)%locFst-Fst_funct_value)
!----------------------------------
    END FUNCTION Fst_loc_funct
!----------------------------------
!----------------------------------
  END SUBROUTINE SA_loc_Fst
!----------------------------------


!----------------------------------
  SUBROUTINE SA_Qst
! LSR 20/07/04
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j,k,Qst_trait
  INTEGER :: nsucc_sa
  INTEGER :: sampled_locus
  REAL(KIND=dp) :: fbest_sa
  REAL(KIND=dp) :: score_sa
  REAL(KIND=dp) :: delta_sa
  REAL(KIND=dp) :: tempo_sa
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_trans    ! test locus backup for transfer
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: best_sa ! best genome backup for SA
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: test_sa ! test genome backup for SA

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: locus_genot
    REAL(KIND=dp), DIMENSION(num_tot_ind) :: indv_genot
    REAL(KIND=dp), DIMENSION(num_pop_tot) :: mean_trait_pop
    REAL(KIND=dp) :: mean_trait_tot
    REAL(KIND=dp) :: SSamong
    REAL(KIND=dp) :: SSwithin
    REAL(KIND=dp) :: MSamong
    REAL(KIND=dp) :: MSwithin
    REAL(KIND=dp) :: expMSamong


  DO Qst_trait=1,num_traits

    ALLOCATE(test_trans(2,num_tot_ind))
    ALLOCATE(best_sa(2,num_sel_loci(Qst_trait),num_tot_ind))
    ALLOCATE(test_sa(2,num_sel_loci(Qst_trait),num_tot_ind))

    FORALL(i=1:2,j=1:num_tot_ind) test_trans(i,j)=0
    FORALL(i=1:2,j=1:num_sel_loci(Qst_trait),k=1:num_tot_ind) best_sa(i,j,k)=0
    FORALL(i=1:2,j=1:num_sel_loci(Qst_trait),k=1:num_tot_ind) test_sa(i,j,k)=0

!write(*,*) 'before translocation sel_genome', sel_genome(1,1,:)
!   copies selected loci genome into SA matrices
    counter=0
    DO i=1,tot_sel_loci
      IF(sel_gene(i)%funct==Qst_trait)THEN
        counter=counter+1
        DO j=1,num_tot_ind
          best_sa(1,counter,j)=sel_genome(1,i,j)
          best_sa(2,counter,j)=sel_genome(2,i,j)
        END DO
      END IF
    END DO

!   SIMULATED ANNEALING ALGORITHM
!   (Pham & Karaboga, 2000)
    tempo_sa=Tini_Qst
    fbest_sa=Qst_funct(best_sa,Qst_trait)
      across_Qst_temp: DO i=1,Tsteps_Qst
        nsucc_sa=0
        within_Qst_temp: DO j=1,Nover_Qst
          test_sa=best_sa
!         picks one locus of (Qst_trait) trait for Qst_trans
          xrand=alea(idum)
          sampled_locus=FLOOR(xrand*num_sel_loci(Qst_trait))+1
          FORALL(k=1:num_tot_ind) test_trans(1,k)=test_sa(1,sampled_locus,k)
          FORALL(k=1:num_tot_ind) test_trans(2,k)=test_sa(2,sampled_locus,k)
          CALL Qst_trans(test_trans)
          FORALL(k=1:num_tot_ind) test_sa(1,sampled_locus,k)=test_trans(1,k)
          FORALL(k=1:num_tot_ind) test_sa(2,sampled_locus,k)=test_trans(2,k)
!         scores new arrangement
          score_sa=Qst_funct(test_sa,Qst_trait)
          delta_sa=score_sa-fbest_sa
          xrand=alea(idum)
          IF((delta_sa<=0.0d0).OR.(xrand<exp(-delta_sa/tempo_sa)))THEN
            nsucc_sa=nsucc_sa+1
            fbest_sa=score_sa
            best_sa=test_sa
          END IF
          IF(nsucc_sa>=Nlimit_Qst) EXIT within_Qst_temp
        END DO within_Qst_temp
        IF(nsucc_sa==0) EXIT across_Qst_temp
        tempo_sa=tempo_sa*Tfactr_Qst
      END DO across_Qst_temp

!   copies SA matrices into selected loci genome
    counter=0
    DO i=1,tot_sel_loci
      IF(sel_gene(i)%funct==Qst_trait)THEN
        counter=counter+1
        DO j=1,num_tot_ind
          sel_genome(1,i,j)=best_sa(1,counter,j)
          sel_genome(2,i,j)=best_sa(2,counter,j)
        END DO
      END IF
    END DO

    DEALLOCATE(best_sa)
    DEALLOCATE(test_sa)
    DEALLOCATE(test_trans)

  END DO

  CONTAINS
!----------------------------------
    SUBROUTINE Qst_trans(Qst_trans_mat)
!   swaps alleles in between individuals
!   coming from different populations
!----------------------------------
    INTEGER, DIMENSION(2,num_tot_ind), INTENT(INOUT) :: Qst_trans_mat
    INTEGER :: swap1,swap1_pop,swap1_haplot,ind_pop1
    INTEGER :: swap2,swap2_pop,swap2_haplot,ind_pop2
    INTEGER :: buffer
    INTEGER :: progress_counter_4

!   picks a 1st population
    xrand=alea(idum)
    swap1_pop=FLOOR(xrand*num_pop_tot)+1
    swap2_pop=swap1_pop
    ind_pop1=ind_per_pop(swap1_pop)
!   picks a 2nd population
		progress_counter_4=0
    DO WHILE(swap1_pop==swap2_pop)

			progress_counter_4 = progress_counter_4+1
			IF (MODULO(progress_counter_4, 100)==0) WRITE(*,*) '[SA_Qst] swaping alleles to reach target Qst, &
					iteration: ', progress_counter_4

      xrand=alea(idum)
      swap2_pop=FLOOR(xrand*num_pop_tot)+1
    END DO
    ind_pop2=ind_per_pop(swap2_pop)
!   picks a 1st candidate & haplotype within 1st population
    xrand=alea(idum)
    swap1=posit_ind_pop(swap1_pop,(FLOOR(xrand*ind_pop1)+1))
    xrand=alea(idum)
    swap1_haplot=FLOOR(xrand*2)+1
!   picks a 2nd candidate & haplotype within 2nd population
    xrand=alea(idum)
    swap2=posit_ind_pop(swap2_pop,(FLOOR(xrand*ind_pop2)+1))
    xrand=alea(idum)
    swap2_haplot=FLOOR(xrand*2)+1
!   swaps alleles
    buffer=Qst_trans_mat(swap1_haplot,swap1)
    Qst_trans_mat(swap1_haplot,swap1)=Qst_trans_mat(swap2_haplot,swap2)
    Qst_trans_mat(swap2_haplot,swap2)=buffer
!----------------------------------
    END SUBROUTINE Qst_trans
!----------------------------------
!----------------------------------
    FUNCTION Qst_funct(Qst_funct_mat,Qst_funct_trait)
!   allows only biallelic loci
!   & additive effects
!
!   Qst_funct_mat(haplotype,selected_loci,individual) stores the selected genome
!   Qst_funct_trait gives the curreny trait
!----------------------------------
    INTEGER :: i_Qst,j_Qst
    INTEGER :: u_bound
    INTEGER :: position
    INTEGER, INTENT(IN) :: Qst_funct_trait
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: Qst_funct_mat  !assumed-shape dummy array
    REAL(KIND=dp) :: Qst_funct
    REAL(KIND=dp) :: Qst_funct_value
    REAL(KIND=dp) :: weight_pop_size
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: locus_genot

    u_bound=UBOUND(Qst_funct_mat,2) ! number of selected loci for current trait
    ALLOCATE(locus_genot(u_bound,num_tot_ind))

!   genotypic values per locus, according to the genotypic scale per locus
    FORALL(i_Qst=1:u_bound,j_Qst=1:num_tot_ind) locus_genot(i_Qst,j_Qst)=0.0_dp
    DO i_Qst=1,u_bound
      position=posit_sel_loci(Qst_funct_trait,i_Qst)
      DO j_Qst=1,num_tot_ind
        IF(Qst_funct_mat(1,i_Qst,j_Qst)+Qst_funct_mat(2,i_Qst,j_Qst)==2)THEN
          locus_genot(i_Qst,j_Qst)=-genotypic_scale(position,1)
        ELSE IF(Qst_funct_mat(1,i_Qst,j_Qst)+Qst_funct_mat(2,i_Qst,j_Qst)==3)THEN
          locus_genot(i_Qst,j_Qst)=genotypic_scale(position,2)
        ELSE IF(Qst_funct_mat(1,i_Qst,j_Qst)+Qst_funct_mat(2,i_Qst,j_Qst)==4)THEN
          locus_genot(i_Qst,j_Qst)=genotypic_scale(position,1)
        END IF
      END DO
    END DO
!   genotypic values per individual: sum across loci for same trait
    FORALL(i_Qst=1:num_tot_ind) indv_genot(i_Qst)=0.0_dp
    DO i_Qst=1,num_tot_ind
      DO j_Qst=1,u_bound
        indv_genot(i_Qst)=indv_genot(i_Qst)+locus_genot(j_Qst,i_Qst)
      END DO
    END DO

!   free memory allocation
    DEALLOCATE(locus_genot)
!   means per population & global means (per trait)
    FORALL(i_Qst=1:num_pop_tot) mean_trait_pop(i_Qst)=0.0_dp
    mean_trait_tot=0.0_dp
    DO i_Qst=1,num_tot_ind
      mean_trait_pop(indv(i_Qst)%pop)=mean_trait_pop(indv(i_Qst)%pop)+indv_genot(i_Qst)
    END DO
    mean_trait_tot=SUM(mean_trait_pop)/num_tot_ind
    DO i_Qst=1,num_pop_tot
      mean_trait_pop(i_Qst)=mean_trait_pop(i_Qst)/ind_per_pop(i_Qst)
    END DO

!   SSamong & MSamong
    SSamong=0.0_dp
    MSamong=0.0_dp
    DO i_Qst=1,num_pop_tot
      SSamong=SSamong+ind_per_pop(i_Qst)*(mean_trait_pop(i_Qst)-mean_trait_tot)**2
    END DO
    MSamong=SSamong/(num_pop_tot-1)
!   SSwithin & MSwithin
    SSwithin=0.0_dp
    MSwithin=0.0_dp
    DO i_Qst=1,num_tot_ind
      indv_genot(i_Qst)=(indv_genot(i_Qst)-mean_trait_pop(indv(i_Qst)%pop))**2
    END DO
    DO i_Qst=1,num_tot_ind
      SSwithin=SSwithin+indv_genot(i_Qst)
    END DO
    MSwithin=SSwithin/(num_tot_ind-num_pop_tot)
!   weighted population size
    weight_pop_size=0.0_dp
    IF(num_pop_tot>1)THEN
      weight_pop_size=&
      & (1.0_dp/(num_pop_tot-1))*(REAL(num_tot_ind,dp)-(REAL(SUM(ind_per_pop**2),dp)/num_tot_ind))
    ELSE IF(num_pop_tot==1)THEN
      weight_pop_size=REAL(num_tot_ind,dp)
    END IF
!   expected MSamong
    expMSamong=0.0_dp
    expMSamong=(MSamong-MSwithin)/weight_pop_size
!   Qst value
!   assuming VarBetween/(VarBetween+2VarWithin)
    Qst_funct_value=0.0_dp
    Qst_funct_value=expMSamong/(expMSamong+2.0_dp*MSwithin)
!   SOLUTION: the squared distance between observed and expected values
    Qst_funct=0.0_dp
    Qst_funct=(exp_Qst(Qst_funct_trait)-Qst_funct_value)*(exp_Qst(Qst_funct_trait)-Qst_funct_value)
!----------------------------------
    END FUNCTION Qst_funct
!----------------------------------
!----------------------------------
  END SUBROUTINE SA_Qst
!----------------------------------

!----------------------------------
  SUBROUTINE meta_GD_sel
! calculates several DG estimates
! in the metapopulation after SA
!----------------------------------
  INTEGER :: i_Qst,j_Qst,k_Qst
  INTEGER :: pop,allele
  INTEGER, DIMENSION(num_pop_tot,tot_sel_loci,max_sel_allele) :: do_Qst_copies
  REAL(KIND=dp) :: do_Qst_SQfreq
  REAL(KIND=dp), DIMENSION(tot_sel_loci,max_sel_allele) :: do_Qst_t_freqs
  REAL(KIND=dp), DIMENSION(num_pop_tot,tot_sel_loci) :: do_Qst_Hloc_pop
  REAL(KIND=dp), DIMENSION(num_pop_tot) :: do_Qst_Hpop

  FORALL(i_Qst=1:num_pop_tot,j_Qst=1:tot_sel_loci,k_Qst=1:max_sel_allele) do_Qst_copies(i_Qst,j_Qst,k_Qst)=0
  FORALL(i_Qst=1:tot_sel_loci,j_Qst=1:max_sel_allele) do_Qst_t_freqs(i_Qst,j_Qst)=0.0_dp
  FORALL(i_Qst=1:num_pop_tot,j_Qst=1:tot_sel_loci,k_Qst=1:max_sel_allele) sel_meta_freq(i_Qst,j_Qst,k_Qst)=0
  FORALL(i_Qst=1:num_pop_tot,j_Qst=1:tot_sel_loci) do_Qst_Hloc_pop(i_Qst,j_Qst)=0.0_dp
  FORALL(i_Qst=1:num_pop_tot) do_Qst_Hpop(i_Qst)=0.0_dp

! counts # allele copies per allele & per locus & per population
! counts # heterozygotes per locus & per population

  DO i_Qst=1,num_tot_ind
    pop=indv(i_Qst)%pop
    DO j_Qst=1,tot_sel_loci
      allele=sel_genome(1,j_Qst,i_Qst)
      do_Qst_copies(pop,j_Qst,allele)=do_Qst_copies(pop,j_Qst,allele)+1
      allele=sel_genome(2,j_Qst,i_Qst)
      do_Qst_copies(pop,j_Qst,allele)=do_Qst_copies(pop,j_Qst,allele)+1
    END DO
  END DO

! calculates allele frequencies per allele & per locus & per population

  DO i_Qst=1,num_pop_tot
    DO j_Qst=1,tot_sel_loci
      DO k_Qst=1,sel_gene(j_Qst)%nalle
!       allele frequencies
        sel_meta_freq(i_Qst,j_Qst,k_Qst)=&
        & REAL(do_Qst_copies(i_Qst,j_Qst,k_Qst),dp)/SUM(do_Qst_copies(i_Qst,j_Qst,:))
      END DO
    END DO
  END DO

! calculates allele frequencies per allele & per locus, across populations
  DO i_Qst=1,tot_sel_loci
    DO j_Qst=1,sel_gene(i_Qst)%nalle
!     allele frequencies
      do_Qst_t_freqs(i_Qst,j_Qst)=REAL(SUM(do_Qst_copies(:,i_Qst,j_Qst)),dp)/SUM(do_Qst_copies(:,i_Qst,:))
    END DO
  END DO

! calculates Ht per locus
  DO i_Qst=1,tot_sel_loci
    do_Qst_SQfreq=0.0_dp
    DO j_Qst=1,sel_gene(i_Qst)%nalle
      do_Qst_SQfreq=do_Qst_SQfreq+(do_Qst_t_freqs(i_Qst,j_Qst)*do_Qst_t_freqs(i_Qst,j_Qst))
    END DO
    locHtsel_Meta(i_Qst)=1.0_dp-do_Qst_SQfreq
    locNasel_Meta(i_Qst)=1.0_dp/do_Qst_SQfreq
  END DO

! calculates Hs_I per population & locus
  DO i_Qst=1,num_pop_tot
    DO j_Qst=1,tot_sel_loci
      do_Qst_SQfreq=0.0_dp
      DO k_Qst=1,sel_gene(j_Qst)%nalle
        do_Qst_SQfreq=do_Qst_SQfreq+(sel_meta_freq(i_Qst,j_Qst,k_Qst)*sel_meta_freq(i_Qst,j_Qst,k_Qst))
      END DO
!     expected heterozygosities per locus & population
      do_Qst_Hloc_pop(i_Qst,j_Qst)=1.0_dp-do_Qst_SQfreq
    END DO
!   expected heterozygosities per population
    do_Qst_Hpop(i_Qst)=SUM(do_Qst_Hloc_pop(i_Qst,:))/tot_sel_loci
  END DO

! calculates Hs per locus
  FORALL(i_Qst=1:tot_sel_loci) locHssel_Meta(i_Qst)=0.0_dp
  DO i_Qst=1,tot_sel_loci
!   weighted average by population size following Nei (1977)
    DO j_Qst=1,num_pop_tot
      locHssel_Meta(i_Qst)=locHssel_Meta(i_Qst)+(do_Qst_Hloc_pop(j_Qst,i_Qst)*REAL(ind_per_pop(j_Qst),dp))/num_tot_ind
    END DO
  END DO

write(*,*) 
!----------------------------------
  END SUBROUTINE meta_GD_sel
!----------------------------------

!----------------------------------
  SUBROUTINE meta_stats
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j,k

  ini_S_Fst_Meta=ini_S_Fst_Meta+ini_Fst_Meta
  ini_Sq_Fst_Meta=ini_Sq_Fst_Meta+(ini_Fst_Meta*ini_Fst_Meta)

  obs_S_Fst_Meta=obs_S_Fst_Meta+obs_Fst_Meta
  obs_Sq_Fst_Meta=obs_Sq_Fst_Meta+(obs_Fst_Meta*obs_Fst_Meta)

  S_Ht_Meta=S_Ht_Meta+Ht_Meta
  Sq_Ht_Meta=Sq_Ht_Meta+(Ht_Meta*Ht_Meta)

  S_Hs_Meta=S_Hs_Meta+Hs_Meta
  Sq_Hs_Meta=Sq_Hs_Meta+(Hs_Meta*Hs_Meta)

  DO i=1,num_traits
    S_ini_Qst_Meta(i)=S_ini_Qst_Meta(i)+ini_Qst_Meta(i)
    Sq_ini_Qst_Meta(i)=Sq_ini_Qst_Meta(i)+(ini_Qst_Meta(i)*ini_Qst_Meta(i))
  END DO

  DO i=1,num_traits
    S_obs_Qst_Meta(i)=S_obs_Qst_Meta(i)+obs_Qst_Meta(i)
    Sq_obs_Qst_Meta(i)=Sq_obs_Qst_Meta(i)+(obs_Qst_Meta(i)*obs_Qst_Meta(i))
  END DO

  DO i=1,num_traits
    S_VGenot_tot(i)=S_VGenot_tot(i)+VGenot_tot(i)
    Sq_VGenot_tot(i)=Sq_VGenot_tot(i)+(VGenot_tot(i)*VGenot_tot(i))
  END DO

  DO i=1,tot_neu_loci
    S_locHs_Meta(i)=S_locHs_Meta(i)+locHs_Meta(i)
    Sq_locHs_Meta(i)=Sq_locHs_Meta(i)+(locHs_Meta(i)*locHs_Meta(i))
  END DO

  DO i=1,tot_neu_loci
    obs_S_locFst_Meta(i)=obs_S_locFst_Meta(i)+obs_locFst_Meta(i)
    obs_Sq_locFst_Meta(i)=obs_Sq_locFst_Meta(i)+(obs_locFst_Meta(i)*obs_locFst_Meta(i))
  END DO

  DO i=1,tot_neu_loci
    S_locNa_Meta(i)=S_locNa_Meta(i)+locNa_Meta(i)
    Sq_locNa_Meta(i)=Sq_locNa_Meta(i)+(locNa_Meta(i)*locNa_Meta(i))
  END DO

  DO i=1,tot_sel_loci
    S_locNasel_Meta(i)=S_locNasel_Meta(i)+locNasel_Meta(i)
    Sq_locNasel_Meta(i)=Sq_locNasel_Meta(i)+(locNasel_Meta(i)*locNasel_Meta(i))
  END DO

  DO i=1,tot_sel_loci
    obs_S_locQst_Meta(i)=obs_S_locQst_Meta(i)+obs_locQst_Meta(i)
    obs_Sq_locQst_Meta(i)=obs_Sq_locQst_Meta(i)+(obs_locQst_Meta(i)*obs_locQst_Meta(i))
  END DO

  DO i=1,tot_neu_loci
    obs_S_locHt_Meta(i)=obs_S_locHt_Meta(i)+obs_locHt_Meta(i)
    obs_Sq_locHt_Meta(i)=obs_Sq_locHt_Meta(i)+(obs_locHt_Meta(i)*obs_locHt_Meta(i))
  END DO

  DO i=1,tot_sel_loci
    S_locHtsel_Meta(i)=S_locHtsel_Meta(i)+locHtsel_Meta(i)
    Sq_locHtsel_Meta(i)=Sq_locHtsel_Meta(i)+(locHtsel_Meta(i)*locHtsel_Meta(i))
  END DO

  DO i=1,tot_sel_loci
    S_locHssel_Meta(i)=S_locHssel_Meta(i)+locHssel_Meta(i)
    Sq_locHssel_Meta(i)=Sq_locHssel_Meta(i)+(locHssel_Meta(i)*locHssel_Meta(i))
  END DO

  DO i=1,num_pop_tot
    DO j=1,num_traits
      S_VGenot_pop(i,j)=S_VGenot_pop(i,j)+VGenot_pop(i,j)
      Sq_VGenot_pop(i,j)=Sq_VGenot_pop(i,j)+(VGenot_pop(i,j)*VGenot_pop(i,j))
    END DO
  END DO

  DO i=1,num_pop_tot
    DO j=1,tot_neu_loci
      DO k=1,neu_gene(j)%nalle
        S_neu_meta_freq(i,j,k)=S_neu_meta_freq(i,j,k)+neu_meta_freq(i,j,k)
        Sq_neu_meta_freq(i,j,k)=Sq_neu_meta_freq(i,j,k)+(neu_meta_freq(i,j,k)*neu_meta_freq(i,j,k))
      END DO
    END DO
  END DO

  DO i=1,num_pop_tot
    DO j=1,tot_sel_loci
      DO k=1,sel_gene(j)%nalle
        S_sel_meta_freq(i,j,k)=S_sel_meta_freq(i,j,k)+sel_meta_freq(i,j,k)
        Sq_sel_meta_freq(i,j,k)=Sq_sel_meta_freq(i,j,k)+(sel_meta_freq(i,j,k)*sel_meta_freq(i,j,k))
      END DO
    END DO
  END DO
!----------------------------------
  END SUBROUTINE meta_stats
!----------------------------------


!----------------------------------
  SUBROUTINE write_output
!----------------------------------
  IMPLICIT NONE
  INTEGER :: i,j,k
	INTEGER :: counter
	REAL(KIND=dp) :: alpha
  REAL(KIND=dp) :: mean,vare
  REAL(KIND=dp) :: mean_mean
  REAL(KIND=dp) :: mean_vare
  REAL(KIND=dp) :: vare_mean
  REAL(KIND=dp) :: vare_vare

  mean=0.0_dp
  vare=0.0_dp

  mean_mean=0.0_dp
  mean_vare=0.0_dp
  vare_mean=0.0_dp
  vare_vare=0.0_dp

  WRITE(21,*)'# Output file of genotype_generator. Can be read by Capsis-Luberon2.'
  WRITE(21,*)'# ', fdate()
  WRITE(21,*)' '

!write(*,100) 'exp_Qst: ', exp_Qst
!100 format(A,F11.2)

  mean=mean_rep(obs_S_Fst_Meta,num_meta)
  WRITE (21,100) 'final_Fst = ', mean
	100 format(A,F11.6)
  !WRITE (21,'(2F11.6)') mean

  WRITE(21,*)' '

  WRITE(21,'(A,I5,A)') 'FloatMatrix name=trait_final_Qst nlin=', num_traits, ' ncol=2'
  WRITE(21, '(A)') '# first column is the trait'
  DO i=1,num_traits
    mean=mean_rep(S_obs_Qst_Meta(i),num_meta)
    WRITE(21,'(I5, A, F11.6)') i, tab, mean
  END DO

  WRITE(21,*)' '

  IF(num_pop_tot==1) THEN
WRITE(21,'(A,I5,A)') 'FloatMatrix name=trait_genotypic_variance nlin=', num_traits, ' ncol=2'
  WRITE(21, '(A)') '# first column is the trait'
  DO i=1,num_traits
    WRITE(21,'(I5, A, F11.4)') i, tab, sd_additive_value(i)*sd_additive_value(i)
  END DO
ELSE
WRITE(21,'(A,I5,A)') 'FloatMatrix name=trait_genotypic_variance nlin=', num_traits, ' ncol=2'
  WRITE(21, '(A)') '# first column is the trait'
  DO i=1,num_traits
    mean=mean_rep(S_VGenot_tot(i),num_meta)
    WRITE(21,'(I5, A, F11.4)') i, tab, mean
  END DO
END IF

  WRITE(21,*)' '

    WRITE(21,'(A,I5,A)') 'FloatMatrix name=trait_disequilibrium_covariance nlin=', num_pop_tot*num_traits*num_traits, ' ncol=4'
    WRITE(21,'(A)') '# population_trait1_trait2'
    DO i=1,num_pop_tot
      DO j=1,num_traits
        DO k=1,num_traits
          WRITE(21,'(I5,A,I5,A,I5,A,F11.6)')i,tab,j,tab,k,tab,end_link_trait(i,j,k)
        END DO
      END DO
    END DO

  WRITE(21,*)' '

  WRITE(21,'(A,I5,A)') 'FloatMatrix name=neutral_locus_Fst nlin=', tot_neu_loci, ' ncol=2'
  WRITE(21, '(A)') '# first column is the position of the locus in the genome (1, 2...)'
	DO i=1,tot_neu_loci
    mean=mean_rep(obs_S_locFst_Meta(i),num_meta)
    WRITE(21,'(I5,A,F11.6)') neu_gene(i)%posit, tab, mean
  END DO 

  WRITE(21,*)' '

	WRITE(21,'(A,I5,A)') 'FloatMatrix name=selected_locus_Qst nlin=', tot_sel_loci, ' ncol=2'
	WRITE(21, '(A)') '# first column is the position of the locus in the genome (1, 2...)'
	DO i=1, num_traits
		DO j=1, tot_sel_loci
			IF (sel_gene(j)%funct==i) THEN
				mean=mean_rep(obs_S_locQst_Meta(j),num_meta)
				WRITE(21,'(I5,A,F11.6)') sel_gene(j)%posit, tab, mean
			END IF
		END DO
	END DO


!  DO i=1,tot_sel_loci
!    mean=mean_rep(obs_S_locQst_Meta(i),num_meta)
!    WRITE(21,'(I5,A,F11.6)')sel_gene(i)%posit-1, tab, mean
!  END DO

  WRITE(21,*)' '

  WRITE(21,'(A,I5,A)') 'FloatMatrix name=selected_locus_genotypic_scale nlin=', tot_sel_loci, ' ncol=2'
  WRITE(21, '(A)') '# first column is the position of the locus in the genome (1, 2...)'

	counter=0
	DO i=1, num_traits
		DO j=1, tot_sel_loci
			IF (sel_gene(j)%funct==i) THEN
				counter=counter+1
				alpha=genotypic_scale(counter,1)
				alpha=alpha*1000
        WRITE(21,'(I5,A,F11.2)') sel_gene(j)%posit,tab, alpha ! posit is the name (1, 2, 3 ...) of the locus whereas posit-1 is the index (0, 1, 2 ...) of the locus
			END IF
		END DO
  END DO

  WRITE(21,*)' '
 
! OPEN(UNIT=3,FILE='allelic_freq_p.txt',STATUS='unknown')

!  DO i=1,num_pop_tot
!    DO j=1,tot_sel_loci
!      DO k=1,sel_gene(j)%nalle
!        write(3, *) mean_rep(S_sel_meta_freq(i,j,k),num_meta)
!      END DO
!    END DO
!  END DO

! close(3)

	counter=0
  WRITE(21,'(A,I5,A)') 'FloatMatrix name=selected_locus_contrib_to_VA nlin=', tot_sel_loci, ' ncol=2'
  WRITE(21, '(A)') '# first column is the position of the locus in the genome (1, 2...)'
	DO i=1, num_traits
		DO j=1, tot_sel_loci
			IF (sel_gene(j)%funct==i) THEN
			 counter=counter+1
       WRITE(21,'(I5,A,F11.6)', advance='no') sel_gene(j)%posit, tab, loci_contrib_VA(counter)
   	write(21, *) '' 
			END IF
		END DO
   END DO

  WRITE(21,*)' '

  WRITE(21,'(A,I5,A,I5)') 'IntMatrix name=neutral_locus_posit nlin=', 1, ' ncol=', tot_neu_loci
	DO j=1, tot_neu_loci
    WRITE(21,'(I5,A)', advance='no') neu_gene(j)%posit, tab
	END DO

   	write(21, *) ' ' 
   	write(21, *) ' ' 

  WRITE(21,'(A,I5,A,I5)') 'IntMatrix name=selected_locus_posit nlin=', 1, ' ncol=', tot_sel_loci
	DO i=1, num_traits
		DO j=1, tot_sel_loci
			IF (sel_gene(j)%funct==i) THEN
       WRITE(21,'(I5,A)', advance='no') sel_gene(j)%posit, tab
			END IF
		END DO
   END DO

   	write(21, *) ' '
   	write(21, *) ' '  


    WRITE(21,'(A,I5,A,I5)') 'IntMatrix name=maternal_selected_genome nlin=', num_tot_ind, ' ncol=', tot_sel_loci
    DO i=1,num_tot_ind
      DO j=1,tot_sel_loci
        WRITE(21,'(I5,A)', advance='no') sel_genome(1, j, i),tab
      END DO
   	write(21, *) '' 
    END DO

  WRITE(21,*)' '

    WRITE(21,'(A,I5,A,I5)') 'IntMatrix name=paternal_selected_genome nlin=', num_tot_ind, ' ncol=', tot_sel_loci
    DO i=1,num_tot_ind
      DO j=1,tot_sel_loci
        WRITE(21,'(I5,A)', advance='no') sel_genome(2, j, i),tab
      END DO
   	write(21, *) '' 
    END DO

  WRITE(21,*)' '

    WRITE(21,'(A,I5,A,I5)') 'IntMatrix name=maternal_neutral_genome nlin=', num_tot_ind, ' ncol=', tot_neu_loci
    DO i=1,num_tot_ind
      DO j=1,tot_neu_loci
        WRITE(21,'(I5,A)', advance='no') neu_genome(1, j, i),tab
      END DO
   	write(21, *) '' 
    END DO

  WRITE(21,*)' '

    WRITE(21,'(A,I5,A,I5)') 'IntMatrix name=paternal_neutral_genome nlin=', num_tot_ind, ' ncol=', tot_neu_loci
    DO i=1,num_tot_ind
      DO j=1,tot_neu_loci
        WRITE(21,'(I5,A)', advance='no') neu_genome(2, j, i),tab
      END DO
   	write(21, *) '' 
    END DO

!----------------------------------
  END SUBROUTINE write_output
!----------------------------------


!----------------------------------
  SUBROUTINE read_real(out)
!----------------------------------

REAL(KIND=dp), INTENT(OUT) :: out
INTEGER  :: my_iostat
CHARACTER (256) :: my_iomsg

input_line = input_line + 1

READ(31, *, iostat=my_iostat, iomsg=my_iomsg) out
IF (my_iostat /= 0) THEN
	WRITE(*,*) 'Problem in input file. See line: ', input_line
	WRITE(*,*) my_iomsg
	STOP 9
END IF
!----------------------------------
  END SUBROUTINE read_real
!----------------------------------

!----------------------------------
  SUBROUTINE read_integer(out)
!----------------------------------

INTEGER, intent(out) :: out
integer  :: my_iostat
character (256) :: my_iomsg

input_line = input_line + 1

READ(31, *, iostat=my_iostat, iomsg=my_iomsg) out
IF (my_iostat /= 0) THEN
	WRITE(*,*) 'Problem in input file. See line: ', input_line
	WRITE(*,*) my_iomsg
	STOP 9
END IF
!----------------------------------
  END SUBROUTINE read_integer
!----------------------------------

!----------------------------------
  SUBROUTINE read_character
!----------------------------------

integer  :: my_iostat
character (256) :: my_iomsg

input_line = input_line + 1

READ(31,'(A)', iostat=my_iostat, iomsg=my_iomsg) heading
IF (my_iostat /= 0) THEN
	WRITE(*,*) 'Problem in input file. See line: ', input_line
	WRITE(*,*) my_iomsg
	STOP 9
END IF
!----------------------------------
  END SUBROUTINE read_character
!----------------------------------

!==================================
  END MODULE metatrom_v7_routines
!==================================

!==================================
  PROGRAM metatrom_v7
!==================================
  USE metatrom_v7_param
  USE metatrom_v7_stat
  USE metatrom_v7_routines
  USE random
  IMPLICIT NONE
  SAVE

  INTEGER :: i,j,k
  INTEGER :: i_samp,j_samp,k_samp
  INTEGER :: restart
	REAL(KIND=dp) :: aux

  integer  :: my_iostat
  character (256) :: my_iomsg

  CALL DATE_AND_TIME(time=ctime)
  PRINT *,' ----------- Genotype generator (metatrom) -----------'
	PRINT *, ctime, 'starting simulation ...'


  OPEN(UNIT=31,FILE='input',STATUS='old',iostat=my_iostat, iomsg=my_iomsg)
IF (my_iostat /= 0) THEN
	WRITE(*,*) 'Cannot open input file.'
	WRITE(*,*) my_iomsg
END IF


  input_line = 0
  !READ(31,'(A)') heading
	!CALL read_character
  !READ(31,'(A)') output
	!CALL read_character
  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) distribution_shape
  CALL read_real(aux)
	distribution_shape=SNGL(aux)
  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) num_tot_ind
	CALL read_integer(num_tot_ind)
  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) tot_neu_loci
	CALL read_integer(tot_neu_loci)
  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) tot_sel_loci
	CALL read_integer(tot_sel_loci)
  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) num_traits
	CALL read_integer(num_traits)

  ALLOCATE(trait_sel_loci(num_traits))
  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,num_traits
    !READ(31,*) trait_sel_loci(i)
		CALL read_integer(trait_sel_loci(i))
  END DO
 ! READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) num_pop_tot
	CALL read_integer(num_pop_tot)
!  IF(num_pop_tot<=1)THEN
!    WRITE(*,*)'ERROR!, num_pop_tot  '
!    WRITE(*,*)'MUST BE > 1  '
!    WRITE(*,*)'ABORT run when input num_pop_tot...'
!    STOP
!  END IF

  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) targetFst_mode
	CALL read_integer(targetFst_mode)
  IF(targetFst_mode==0)THEN
    !READ(31,'(A)') heading
		CALL read_character
    !READ(31,*) exp_Fst
		CALL read_real(exp_Fst)
  END IF
  !READ(31,'(A)') heading
		CALL read_character
  !READ(31,*) num_meta
	CALL read_integer(num_meta)

! target additive value
  ALLOCATE(target_additive_value(num_traits)) ! TYPE

  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,num_traits
    !READ(31,*) target_additive_value(i)%mean
		CALL read_real(target_additive_value(i)%mean_threshold)
  END DO
  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,num_traits
    !READ(31,*) target_additive_value(i)%sd
		CALL read_real(target_additive_value(i)%sd)
  END DO

  !READ(31,'(A)') heading
	CALL read_character
  !READ(31,*) chisq_threshold
  CALL read_real(chisq_threshold)

! neutral loci
  ALLOCATE(neu_gene(tot_neu_loci)) ! TYPE

  max_neu_allele=0
  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,tot_neu_loci
    !READ(31,*) neu_gene(i)%chrom
		CALL read_integer(neu_gene(i)%chrom)
  END DO
  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,tot_neu_loci
    !READ(31,*) neu_gene(i)%posit
		CALL read_integer(neu_gene(i)%posit)
  END DO

  IF(targetFst_mode==1)THEN
    !READ(31,'(A)') heading
		CALL read_character
    DO i=1,tot_neu_loci
      !READ(31,*) neu_gene(i)%locFst
			CALL read_real(neu_gene(i)%locFst)
    END DO
  END IF

  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,tot_neu_loci
    !READ(31,*) neu_gene(i)%nalle
		CALL read_integer(neu_gene(i)%nalle)
    IF(neu_gene(i)%nalle>max_neu_allele) max_neu_allele=neu_gene(i)%nalle
  END DO

! selected loci
  ALLOCATE(sel_gene(tot_sel_loci)) ! TYPE
  ALLOCATE(num_sel_loci(num_traits))
    FORALL(i=1:num_traits) num_sel_loci(i)=0

  max_sel_allele=2
  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,tot_sel_loci
    !READ(31,*) sel_gene(i)%funct
		CALL read_integer(sel_gene(i)%funct)
  END DO

  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,tot_sel_loci
    !READ(31,*) sel_gene(i)%chrom
		CALL read_integer(sel_gene(i)%chrom)
  END DO

  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,tot_sel_loci
    !READ(31,*) sel_gene(i)%posit
		CALL read_integer(sel_gene(i)%posit)
  END DO
  max_loci_trait=0
  DO i=1,num_traits
    num_sel_loci(i)=COUNT(sel_gene(:)%funct==i)
    IF(num_sel_loci(i)>max_loci_trait) max_loci_trait=num_sel_loci(i)
  END DO
  IF(SUM(num_sel_loci)/=tot_sel_loci)THEN
    PRINT *, 'ERROR, SUM(num_sel_loci)/=tot_sel_loci...'
    PRINT *, 'aborted...'
    STOP
  END IF

! ALLOCATE

! storage ChiSq test for HWE before/after SA
  ALLOCATE(chisq_before_SA(4,num_pop_tot,tot_neu_loci))
    FORALL(i=1:4,j=1:num_pop_tot,k=1:tot_neu_loci) chisq_before_SA(i,j,k)=0.0_dp
  ALLOCATE(chisq_after_SA(4,num_pop_tot,tot_neu_loci))
    FORALL(i=1:4,j=1:num_pop_tot,k=1:tot_neu_loci) chisq_after_SA(i,j,k)=0.0_dp

  ALLOCATE(exp_Qst(num_traits))
    FORALL(i=1:num_traits) exp_Qst(i)=0.0_dp
  ALLOCATE(ini_Qst_Meta(num_traits))
!   set to zero in subroutine
  ALLOCATE(S_ini_Qst_Meta(num_traits))
    FORALL(i=1:num_traits) S_ini_Qst_Meta(i)=0.0_dp
  ALLOCATE(Sq_ini_Qst_Meta(num_traits))
    FORALL(i=1:num_traits) Sq_ini_Qst_Meta(i)=0.0_dp
  ALLOCATE(obs_Qst_Meta(num_traits))
!   set to zero in subroutine
  ALLOCATE(S_obs_Qst_Meta(num_traits))
    FORALL(i=1:num_traits) S_obs_Qst_Meta(i)=0.0_dp
  ALLOCATE(Sq_obs_Qst_Meta(num_traits))
    FORALL(i=1:num_traits) Sq_obs_Qst_Meta(i)=0.0_dp

 ! ALLOCATE(target_additive_value(num_traits)%mean)
    FORALL(i=1:num_traits) target_additive_value(i)%mean=0.0_dp

  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,num_traits
    !READ(31,*) exp_Qst(i)
		CALL read_real(exp_Qst(i))
  END DO

  ALLOCATE(ind_per_pop(num_pop_tot))
    FORALL(i=1:num_pop_tot) ind_per_pop(i)=0

  max_ind_pop=0
  !READ(31,'(A)') heading
	CALL read_character
  DO i=1,num_pop_tot
    !READ(31,*) ind_per_pop(i)
		CALL read_integer(ind_per_pop(i))
    IF(ind_per_pop(i)>max_ind_pop) max_ind_pop=ind_per_pop(i)
  END DO

  num_tot_ind=SUM(ind_per_pop)
  num_tot_haplot=2*num_tot_ind

! List of parameters read
write(*,*) 'List of parameters read in input file:'
write(*,*) 'distribution_shape: ', distribution_shape
write(*,*) 'num_tot_ind: ', num_tot_ind
write(*,*) 'tot_neu_loci: ', tot_neu_loci
write(*,*) 'tot_sel_loci: ', tot_sel_loci
write(*,*) 'num_traits: ', num_traits
write(*,*) 'trait_sel_loci: ', trait_sel_loci
write(*,*) 'num_pop_tot: ', num_pop_tot
write(*,*) 'targetFst_mode: ', targetFst_mode
write(*,*) 'exp_Fst: ', exp_Fst
write(*,*) 'num_meta: ', num_meta

write(*,*) 'target_additive_value_mean, size: ', size(target_additive_value), '...'
do i=1, size(target_additive_value)
 write(*,'(F11.6)',advance='no') target_additive_value(i)%mean
end do
write(*,*)

write(*,*) 'threshold_value_for_additive_mean, size: ', size(target_additive_value), '...'
do i=1, size(target_additive_value)
 write(*,'(F11.6)',advance='no') target_additive_value(i)%mean_threshold
end do
write(*,*)

write(*,*) 'target_additive_value_sd, size: ', size(target_additive_value), '...'
do i=1, size(target_additive_value)
 write(*,'(F11.6)',advance='no') target_additive_value(i)%sd
end do
write(*,*)

write(*,*) 'chisq_threshold: ', chisq_threshold

write(*,*) 'neu_gene_chrom, size: ', size(neu_gene), '...'
do i=1, size(neu_gene)
 write(*,'(I5)',advance='no') neu_gene(i)%chrom
end do
write(*,*)

write(*,*) 'neu_gene_posit, size: ', size(neu_gene), '...'
do i=1, size(neu_gene)
 write(*,'(I5)',advance='no') neu_gene(i)%posit
end do
write(*,*)

write(*,*) 'neu_gene_nalle, size: ', size(neu_gene), '...'
do i=1, size(neu_gene)
 write(*,'(I5)',advance='no') neu_gene(i)%nalle
end do
write(*,*)

write(*,*) 'sel_gene_funct, size: ', size(sel_gene), '...'
do i=1, size(sel_gene)
 write(*,'(I5)',advance='no') sel_gene(i)%funct
end do
write(*,*)

write(*,*) 'sel_gene_chrom, size: ', size(sel_gene), '...'
do i=1, size(sel_gene)
 write(*,'(I5)',advance='no') sel_gene(i)%chrom
end do
write(*,*)

write(*,*) 'sel_gene_posit, size: ', size(sel_gene), '...'
do i=1, size(sel_gene)
 write(*,'(I5)',advance='no') sel_gene(i)%posit
end do
write(*,*)

!write(*,100) 'exp_Qst: ', exp_Qst
!100 format(A,F11.2)
write(*,*) 'exp_Qst: ', exp_Qst

write(*,*) 'ind_per_pop: ', ind_per_pop


  ALLOCATE(genotypic_scale(tot_sel_loci,2))
    FORALL(i=1:tot_sel_loci,j=1:2) genotypic_scale(i,j)=0.0_dp


ALLOCATE(loci_contrib_VA(tot_sel_loci))
FORALL(i=1:tot_sel_loci) loci_contrib_VA(i)=0.0_dp

ALLOCATE(mean_additive_value(num_traits+1))		! a mean additive value per trait and another for neutral genome
ALLOCATE(sd_additive_value(num_traits+1))
! set & initialize additional ALLOCATABLE VARIABLES
  ALLOCATE(indv(num_tot_ind))
!   set to zero in subroutine
  ALLOCATE(neu_copies_dist(tot_neu_loci,max_neu_allele))
!   set to zero in subroutine
  ALLOCATE(sel_copies_dist(tot_sel_loci,max_sel_allele))
!   set to zero in subroutine
  ALLOCATE(neu_genome(2,tot_neu_loci,num_tot_ind))
    FORALL(i=1:tot_neu_loci,j=1:num_tot_ind) neu_genome(1,i,j)=0
    FORALL(i=1:tot_neu_loci,j=1:num_tot_ind) neu_genome(2,i,j)=0
  ALLOCATE(sel_genome(2,tot_sel_loci,num_tot_ind))
    FORALL(i=1:tot_sel_loci,j=1:num_tot_ind) sel_genome(1,i,j)=0
    FORALL(i=1:tot_sel_loci,j=1:num_tot_ind) sel_genome(2,i,j)=0
  ALLOCATE(posit_sel_loci(num_traits,max_loci_trait))
    FORALL(i=1:num_traits,j=1:max_loci_trait) posit_sel_loci(i,j)=0
  ALLOCATE(posit_ind_pop(num_pop_tot,max_ind_pop))
    FORALL(i=1:num_pop_tot,j=1:max_ind_pop) posit_ind_pop(i,j)=0
  ALLOCATE(ini_link_trait(num_pop_tot,num_traits,num_traits))
!   set to zero in subroutine
  ALLOCATE(end_link_trait(num_pop_tot,num_traits,num_traits))
!   set to zero in subroutine
  ALLOCATE(VGenot_tot(num_traits))
    FORALL(i=1:num_traits) VGenot_tot(i)=0.0_dp
  ALLOCATE(S_VGenot_tot(num_traits))
    FORALL(i=1:num_traits) S_VGenot_tot(i)=0.0_dp
  ALLOCATE(Sq_VGenot_tot(num_traits))
    FORALL(i=1:num_traits) Sq_VGenot_tot(i)=0.0_dp
  ALLOCATE(VEnvir_tot(num_traits))
    FORALL(i=1:num_traits) VEnvir_tot(i)=0.0_dp
  ALLOCATE(VGenot_pop(num_pop_tot,num_traits))
    FORALL(i=1:num_pop_tot,j=1:num_traits) VGenot_pop(i,j)=0.0_dp
  ALLOCATE(S_VGenot_pop(num_pop_tot,num_traits))
    FORALL(i=1:num_pop_tot,j=1:num_traits) S_VGenot_pop(i,j)=0.0_dp
  ALLOCATE(Sq_VGenot_pop(num_pop_tot,num_traits))
    FORALL(i=1:num_pop_tot,j=1:num_traits) Sq_VGenot_pop(i,j)=0.0_dp
  ALLOCATE(obs_locFst_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) obs_locFst_Meta(i)=0.0_dp
  ALLOCATE(obs_S_locFst_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) obs_S_locFst_Meta(i)=0.0_dp
  ALLOCATE(obs_Sq_locFst_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) obs_Sq_locFst_Meta(i)=0.0_dp
  ALLOCATE(obs_locQst_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) obs_locQst_Meta(i)=0.0_dp
  ALLOCATE(obs_S_locQst_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) obs_S_locQst_Meta(i)=0.0_dp
  ALLOCATE(obs_Sq_locQst_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) obs_Sq_locQst_Meta(i)=0.0_dp
  ALLOCATE(obs_locHt_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) obs_locHt_Meta(i)=0.0_dp
  ALLOCATE(obs_S_locHt_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) obs_S_locHt_Meta(i)=0.0_dp
  ALLOCATE(obs_Sq_locHt_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) obs_Sq_locHt_Meta(i)=0.0_dp
  ALLOCATE(locHtsel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) locHtsel_Meta(i)=0.0_dp
  ALLOCATE(S_locHtsel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) S_locHtsel_Meta(i)=0.0_dp
  ALLOCATE(Sq_locHtsel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) Sq_locHtsel_Meta(i)=0.0_dp
  ALLOCATE(locHssel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) locHssel_Meta(i)=0.0_dp
  ALLOCATE(S_locHssel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) S_locHssel_Meta(i)=0.0_dp
  ALLOCATE(Sq_locHssel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) Sq_locHssel_Meta(i)=0.0_dp
  ALLOCATE(locHs_Meta(tot_neu_loci))
!   set to zero in subroutine
  ALLOCATE(S_locHs_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) S_locHs_Meta(i)=0.0_dp
  ALLOCATE(Sq_locHs_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) Sq_locHs_Meta(i)=0.0_dp
  ALLOCATE(locNa_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) locNa_Meta(i)=0.0_dp
  ALLOCATE(S_locNa_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) S_locNa_Meta(i)=0.0_dp
  ALLOCATE(Sq_locNa_Meta(tot_neu_loci))
    FORALL(i=1:tot_neu_loci) Sq_locNa_Meta(i)=0.0_dp
  ALLOCATE(locNasel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) locNasel_Meta(i)=0.0_dp
  ALLOCATE(S_locNasel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) S_locNasel_Meta(i)=0.0_dp
  ALLOCATE(Sq_locNasel_Meta(tot_sel_loci))
    FORALL(i=1:tot_sel_loci) Sq_locNasel_Meta(i)=0.0_dp
  ALLOCATE(neu_meta_freq(num_pop_tot,tot_neu_loci,max_neu_allele))
    FORALL(i=1:num_pop_tot,j=1:tot_neu_loci,k=1:max_neu_allele) neu_meta_freq(i,j,k)=0.0_dp
  ALLOCATE(S_neu_meta_freq(num_pop_tot,tot_neu_loci,max_neu_allele))
    FORALL(i=1:num_pop_tot,j=1:tot_neu_loci,k=1:max_neu_allele) S_neu_meta_freq(i,j,k)=0.0_dp
  ALLOCATE(Sq_neu_meta_freq(num_pop_tot,tot_neu_loci,max_neu_allele))
    FORALL(i=1:num_pop_tot,j=1:tot_neu_loci,k=1:max_neu_allele) Sq_neu_meta_freq(i,j,k)=0.0_dp
  ALLOCATE(sel_meta_freq(num_pop_tot,tot_sel_loci,max_sel_allele))
    FORALL(i=1:num_pop_tot,j=1:tot_sel_loci,k=1:max_sel_allele) sel_meta_freq(i,j,k)=0.0_dp
  ALLOCATE(S_sel_meta_freq(num_pop_tot,tot_sel_loci,max_sel_allele))
    FORALL(i=1:num_pop_tot,j=1:tot_sel_loci,k=1:max_sel_allele) S_sel_meta_freq(i,j,k)=0.0_dp
  ALLOCATE(Sq_sel_meta_freq(num_pop_tot,tot_sel_loci,max_sel_allele))
    FORALL(i=1:num_pop_tot,j=1:tot_sel_loci,k=1:max_sel_allele) Sq_sel_meta_freq(i,j,k)=0.0_dp



  Ht_Meta=0.0_dp
  S_Ht_Meta=0.0_dp
  Sq_Ht_Meta=0.0_dp
  Hs_Meta=0.0_dp
  S_Hs_Meta=0.0_dp
  Sq_Hs_Meta=0.0_dp
  ini_Fst_Meta=0.0_dp
  ini_S_Fst_Meta=0.0_dp
  ini_Sq_Fst_Meta=0.0_dp
  obs_Fst_Meta=0.0_dp
  obs_S_Fst_Meta=0.0_dp
  obs_Sq_Fst_Meta=0.0_dp



! end set & initialize

  CLOSE(31) ! CLOSE input file input

! be good & set & record a random number seed
!-----------------------------------------------
! old version from compag visual FORTRAN
!  iseed(1)=0
!  iseed(2)=0
!  xdum(1)=0.0_dp
!  xdum(2)=0.0_dp
!  CALL system_clock(iseed(1))
!  iseed(2)=-0.5*iseed(1)
!  CALL random_seed(put=iseed)
!  CALL random_number(xdum)
!  iseed(:)=INT(1.0e8*xdum(:))
!  IF(1.0e8*xdum(1)-REAL(iseed(1),KIND(xdum))<0.5) iseed(1)=-iseed(1)
!  IF(1.0e8*xdum(2)-REAL(iseed(2),KIND(xdum))>0.5) iseed(2)=-iseed(2)
!  CALL random_seed(put=iseed)
!-----------------------------------------------
! alternative version for compatibility GNU gfortran
!
  CALL system_clock(iseed(1))
  iseed(2)=-0.5*iseed(1)
  i=SIZE(iseed)
  CALL random_seed
  CALL random_seed(SIZE=i)
  CALL random_seed(PUT=iseed(1:i))
  CALL random_number(xdum)
  iseed(:)=INT(1.0e8*xdum(:))
  IF(1.0e8*xdum(1)-REAL(iseed(1),KIND(xdum))<0.5) iseed(1)=-iseed(1)
  IF(1.0e8*xdum(2)-REAL(iseed(2),KIND(xdum))>0.5) iseed(2)=-iseed(2)
  CALL random_seed(PUT=iseed(1:i))
!-----------------------------------------------

  IF(iseed(1)<0)THEN
    idum=iseed(1)
  ELSE
    idum=-iseed(1)
  END IF
  IF(iseed(2)<0)THEN
    idum=iseed(2)
  ELSE
    idum=-iseed(2)
  END IF

  OPEN(UNIT=21,FILE='output',STATUS='unknown')
!  WRITE(21,'(A)') 'Random number seeds ... '
!  WRITE(21,'(2I16)') iseed

! end setting seeds


restart=1

IF(num_pop_tot==1) THEN
	CALL allelic_freq
	CALL write_output
ELSE IF (num_pop_tot>1) THEN

	progress_counter=0
	DO WHILE (restart == 1)

		progress_counter = progress_counter+1
		IF (MODULO(progress_counter, 100)==0) WRITE(*,*) '[program] sampling allelic frequencies, iteration: ', progress_counter

		done_replicate=0
		CALL allelic_freq
		IF (tot_neu_loci/=0) THEN
			CALL meta_Ht_neu
			CALL do_neutral_freq
			CALL do_Fst
			CALL MonteCarlo_test_HWE
			HWE : DO i=1,tot_neu_loci
				DO j=1,num_pop_tot
					IF (chisq_before_SA(4,j,i)<= chisq_threshold) THEN 
						restart=1
						write(*,*) 'restart allelic frequences sample'
						EXIT HWE
					ELSE 
						restart=0
								!write(*,*) 'allelic frequences OK'
					END IF
				END DO
			END DO HWE
		ELSE
			restart=0
		END IF

	END DO
	CALL do_Qst

END IF


!DO k=1, num_traits
!	write(*,*) 'VGenot_tot1', VGenot_tot(k)
!END DO

IF (num_pop_tot>1) THEN
	CALL do_linkage

	restart=1

	progress_counter = 0
	DO WHILE (restart == 1)

		progress_counter = progress_counter+1
		IF (MODULO(progress_counter, 100)==0) WRITE(*,*) '[program] switching of alleles, iteration: ', progress_counter

			IF(tot_neu_loci>0) THEN

				IF (num_pop_tot > 1) THEN
					IF(targetFst_mode==0)THEN
						CALL SA_Fst
					ELSE IF(targetFst_mode==1)THEN
						CALL SA_loc_Fst
					END IF
					CALL SA_Qst

					done_replicate=1
					CALL do_neutral_freq
					CALL do_Fst
					CALL MonteCarlo_test_HWE
					HWE2 : DO i=1,tot_neu_loci
						DO j=1,num_pop_tot
							IF (chisq_after_SA(4,j,i)<=chisq_threshold) THEN 
								restart=1
								write(*,*) 'restart switch of alleles'
								EXIT HWE2
							ELSE 
								restart=0
							END IF
						END DO
					END DO HWE2

				CALL do_Qst

				END IF

		ELSE IF (tot_neu_loci==0) THEN

			IF (num_pop_tot > 1) THEN
				CALL SA_Qst
				done_replicate=1
				restart=0
				CALL do_Qst
			END IF

		END IF

	END DO
		 
	!CALL do_Qst

	!DO k=1, num_traits
	!	write(*,*) 'VGenot_tot2', VGenot_tot(k)
	!END DO

	CALL meta_GD_sel
	CALL do_linkage

	CALL meta_stats
	CALL write_output

END IF

 CLOSE(21)

write(*,*) 'END OF METATROM'

! DEALLOCATING MEMORY BEFORE EXIT
  DEALLOCATE(neu_gene)
  DEALLOCATE(sel_gene)
  DEALLOCATE(target_additive_value)
  DEALLOCATE(ind_per_pop)
  DEALLOCATE(trait_sel_loci)
  DEALLOCATE(num_sel_loci)
  DEALLOCATE(indv)
  DEALLOCATE(neu_copies_dist)
  DEALLOCATE(sel_copies_dist)
  DEALLOCATE(neu_genome)
  DEALLOCATE(sel_genome)
  DEALLOCATE(genotypic_scale)
  DEALLOCATE(posit_sel_loci)
  DEALLOCATE(posit_ind_pop)
  DEALLOCATE(exp_Qst)
  DEALLOCATE(ini_Qst_Meta)
  DEALLOCATE(S_ini_Qst_Meta)
  DEALLOCATE(Sq_ini_Qst_Meta)
  DEALLOCATE(obs_Qst_Meta)
  DEALLOCATE(S_obs_Qst_Meta)
  DEALLOCATE(Sq_obs_Qst_Meta)
  DEALLOCATE(ini_link_trait)
  DEALLOCATE(end_link_trait)
  DEALLOCATE(VGenot_tot)
  DEALLOCATE(S_VGenot_tot)
  DEALLOCATE(Sq_VGenot_tot)
  DEALLOCATE(VEnvir_tot)
  DEALLOCATE(VGenot_pop)
  DEALLOCATE(S_VGenot_pop)
  DEALLOCATE(Sq_VGenot_pop)
  DEALLOCATE(obs_locFst_Meta)
  DEALLOCATE(obs_S_locFst_Meta)
  DEALLOCATE(obs_Sq_locFst_Meta)
  DEALLOCATE(obs_locQst_Meta)
  DEALLOCATE(obs_S_locQst_Meta)
  DEALLOCATE(obs_Sq_locQst_Meta)
  DEALLOCATE(obs_locHt_Meta)
  DEALLOCATE(obs_S_locHt_Meta)
  DEALLOCATE(obs_Sq_locHt_Meta)
  DEALLOCATE(locHtsel_Meta)
  DEALLOCATE(S_locHtsel_Meta)
  DEALLOCATE(Sq_locHtsel_Meta)
  DEALLOCATE(locHssel_Meta)
  DEALLOCATE(S_locHssel_Meta)
  DEALLOCATE(Sq_locHssel_Meta)
  DEALLOCATE(locHs_Meta)
  DEALLOCATE(S_locHs_Meta)
  DEALLOCATE(Sq_locHs_Meta)
  DEALLOCATE(locNa_Meta)
  DEALLOCATE(S_locNa_Meta)
  DEALLOCATE(Sq_locNa_Meta)
  DEALLOCATE(locNasel_Meta)
  DEALLOCATE(S_locNasel_Meta)
  DEALLOCATE(Sq_locNasel_Meta)
  DEALLOCATE(neu_meta_freq)
  DEALLOCATE(S_neu_meta_freq)
  DEALLOCATE(Sq_neu_meta_freq)
  DEALLOCATE(sel_meta_freq)
  DEALLOCATE(S_sel_meta_freq)
  DEALLOCATE(Sq_sel_meta_freq)
  DEALLOCATE(chisq_before_SA,chisq_after_SA)
  DEALLOCATE(mean_additive_value)
  DEALLOCATE(sd_additive_value)
  DEALLOCATE(loci_contrib_VA)

!==================================
  END PROGRAM metatrom_v7
!==================================
