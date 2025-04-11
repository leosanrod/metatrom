#include <iostream>
#include <vector>
#include <string>
#include <cmath>

// Precision definitions
constexpr int sp = 7;  // Precision up to 7 decimal digits
constexpr int dp = 14; // Precision up to 14 decimal digits

// Constants
constexpr int saved_reps = 100;  // Number of replicates to store per metapopulation
constexpr int tot_HWE_rep = 1000; // Number of random generations for HW proportions

// Random seed storage
struct RandomSeed {
    int seed1;
    int seed2;
};

// Data structure definitions
struct Candidate {
    int pop; // Population to which the candidate belongs
};

struct NeuLocus {
    int nalle;           // Number of distinct alleles
    int chrom;           // Chromosome
    int posit;           // Locus position in the chromosome
    double distr;        // Parameter for shaping distribution of allele frequencies
    double locFst;       // Target Fst value
};

struct SelLocus {
    int nalle;           // Number of distinct alleles
    int chrom;           // Chromosome
    int funct;           // Trait number
    int posit;           // Locus position in the chromosome
};

struct ExpectedAdditiveValue {
    double mean;          // Expected mean for additive value
    double mean_threshold; // Variation of mean additive value
    double sd;            // Expected standard deviation for additive value
};

// Global variables
RandomSeed iseed;
int idum;
float distribution_shape;

// Pedigree and gene information
std::vector<Candidate> indv;                     // Stores pedigree
std::vector<NeuLocus> neu_gene;                  // Stores gene information
std::vector<SelLocus> sel_gene;                  // Stores gene information
std::vector<ExpectedAdditiveValue> target_additive_value; // Additive value information

// General counters and system-wide parameters
int done_replicate;
int counter;
int tot_neu_loci;        // Number of neutral loci in the system
int samp_neu_loci;       // Number of neutral loci to sample
int tot_sel_loci;        // Number of selected loci in the system
int num_traits;          // Number of selective traits
int num_pop_tot;         // Number of populations in the original metapopulation
int num_tot_ind;         // Total number of individuals across populations
int num_tot_haplot;      // Total number of haplotypes: 2*num_tot_ind
int max_neu_allele;      // Maximum number of neutral alleles across loci
int max_sel_allele;      // Maximum number of selected alleles across loci
int max_loci_trait;      // Maximum number of loci for a trait
int max_ind_pop;         // Maximum number of individuals in a population
int jlast;               // Rank code for the last individual
int num_meta;            // Number of generated metapopulations
int current_meta;        // Current metapopulation
int targetFst_mode;      // Global Fst (0) or Fst specified by locus (1)
int input_line;          // To track the line of reading problems
int progress_counter;    // To track progress in loops

// Other global variables
std::vector<int> trait_sel_loci;
std::vector<int> num_sel_loci;         // Number of selective loci in the system
std::vector<int> ind_per_pop;          // Number of individuals per population

// Data for genetics and genome positioning
std::vector<std::vector<int>> posit_sel_loci;    // Genome position per selective locus
std::vector<std::vector<int>> posit_ind_pop;     // Position in pedigree per population
std::vector<std::vector<int>> neu_copies_dist;   // Distribution of neutral allelic copies
std::vector<std::vector<int>> sel_copies_dist;   // Distribution of selected allelic copies

// Genome data
std::vector<std::vector<std::vector<int>>> neu_genome;    // Neutral genome per haplotype, locus, individual
std::vector<std::vector<std::vector<int>>> sel_genome;    // Selected genome per haplotype, locus, individual

// Variables for random numbers and statistical calculations
double xrand;         // Variable for random numbers
double Ht_Meta;       // Average Ht at neutral loci in the metapopulation
double S_Ht_Meta;     // Sum of average Ht's across metapopulations
double Sq_Ht_Meta;    // Sum of squares of average Ht's across metapopulations
double Hs_Meta;       // Observed Hs at neutral loci in the metapopulation
double S_Hs_Meta;     // Sum of observed Hs's across metapopulations
double Sq_Hs_Meta;    // Sum of squares of observed Hs's across metapopulations
double exp_Fst;       // Expected Fst: Target of the SA algorithm
double ini_Fst_Meta;  // Current initial Fst in the metapopulation
double ini_S_Fst_Meta; // Sum of initial Fst's across metapopulations
double ini_Sq_Fst_Meta; // Sum of squares of initial Fst's across metapopulations
double obs_Fst_Meta;  // Current final Fst in the metapopulation
double obs_S_Fst_Meta; // Sum of final Fst's across metapopulations
double obs_Sq_Fst_Meta; // Sum of squares of final Fst's across metapopulations

// Random numbers and additional statistics
std::array<double, 2> xdum;
std::vector<double> loci_contrib_VA;   // Contributions of loci to additive variance
std::vector<double> obs_locFst_Meta;  // Current Fst per locus in the metapopulation
std::vector<double> obs_S_locFst_Meta; // Sum of Fst's per locus across metapopulations
std::vector<double> obs_Sq_locFst_Meta; // Sum of squares of Fst's per locus across metapopulations
std::vector<double> obs_locQst_Meta;  // Current Qst per locus in the metapopulation
std::vector<double> obs_S_locQst_Meta; // Sum of Qst's per locus across metapopulations
std::vector<double> obs_Sq_locQst_Meta; // Sum of squares of Qst's per locus across metapopulations
std::vector<double> obs_locHt_Meta;   // Current Het per neutral locus in the metapopulation
std::vector<double> obs_S_locHt_Meta; // Sum of Het's per neutral locus across metapopulations
std::vector<double> obs_Sq_locHt_Meta; // Sum of squares of Het's per neutral locus across metapopulations
std::vector<double> locHtsel_Meta;    // Current Ht per selected locus
std::vector<double> S_locHtsel_Meta;  // Sum of Ht's per selected locus across metapopulations
std::vector<double> Sq_locHtsel_Meta; // Sum of squares of Ht's per selected locus across metapopulations
std::vector<double> locHssel_Meta;    // Current Hs per selected locus
std::vector<double> S_locHssel_Meta;  // Sum of Hs's per selected locus across metapopulations
std::vector<double> Sq_locHssel_Meta; // Sum of squares of Hs's per selected locus across metapopulations

std::vector<double> exp_Qst;          // Expected Qst
std::vector<double> ini_Qst_Meta;     // Initial Qst in metapopulation
std::vector<double> S_ini_Qst_Meta;   // Sum of initial Qst's across metapopulations
std::vector<double> Sq_ini_Qst_Meta;  // Sum of squares of initial Qst's across metapopulations
std::vector<double> obs_Qst_Meta;     // Current final Qst in the metapopulation
std::vector<double> S_obs_Qst_Meta;   // Sum of final Qst's across metapopulations
std::vector<double> Sq_obs_Qst_Meta;  // Sum of squares of final Qst's across metapopopulations

std::vector<double> VGenot_tot;       // Genotypic variance per trait
std::vector<double> S_VGenot_tot;     // Sum of genotypic variances across metapopulations
std::vector<double> Sq_VGenot_tot;    // Sum of squares of genotypic variances across metapopulations
std::vector<double> VEnvir_tot;       // Environmental variance per trait