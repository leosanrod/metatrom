#include <vector>
#include <cmath>
#include <array>

// Variables for metapopulation statistics
std::vector<std::vector<double>> VGenot_pop;           // Genotypic variance per trait and population
std::vector<std::vector<double>> S_VGenot_pop;         // Sum of genotypic variances across metapopulations
std::vector<std::vector<double>> Sq_VGenot_pop;        // Sum of squares of genotypic variances across metapopulations
std::vector<std::vector<double>> genotypic_scale;      // Genotypic values per locus

std::vector<std::vector<std::vector<double>>> ini_link_trait; // Genotypic covariance due to linkage
std::vector<std::vector<std::vector<double>>> end_link_trait; // Genotypic covariance due to linkage

// Allele frequency statistics
std::vector<std::vector<std::vector<double>>> neu_meta_freq;   // Allele frequencies at neutral loci in metapopulation
std::vector<std::vector<std::vector<double>>> S_neu_meta_freq; // Sum of allele frequencies at neutral loci
std::vector<std::vector<std::vector<double>>> Sq_neu_meta_freq; // Square of allele frequencies at neutral loci
std::vector<std::vector<std::vector<double>>> sel_meta_freq;   // Allele frequencies at selected loci in metapopulation
std::vector<std::vector<std::vector<double>>> S_sel_meta_freq; // Sum of allele frequencies at selected loci
std::vector<std::vector<std::vector<double>>> Sq_sel_meta_freq; // Square of allele frequencies at selected loci

// Chi-squared test statistics
double chisq_threshold; // Chi-squared significance threshold
std::vector<std::vector<std::vector<double>>> chisq_before_SA; // Chi-squared test HWE across loci/populations before SA
std::vector<std::vector<std::vector<double>>> chisq_after_SA;  // Chi-squared test HWE across loci/populations after SA

// Additive value statistics
std::vector<double> mean_additive_value; // Mean additive value realized
std::vector<double> sd_additive_value;   // Standard deviation (sqrt of additive variance) realized

// Time and input headers
std::string ctime;   // Stores system time
std::string heading; // Stores input headings

// Simulated Annealing Parameters for Fst
constexpr int Tsteps_Fst = 600; // Maximum number of loops before reaching cooled state
constexpr int Nover_Fst = 300;  // Maximum number of loops per temperature level
constexpr int Nlimit_Fst = 100; // Number of successful hits per temperature level
constexpr double Tini_Fst = 1.0; // Initial temperature
constexpr double Tfactr_Fst = 0.9; // Cooling rate of temperature change

// Simulated Annealing Parameters for Qst
constexpr int Tsteps_Qst = 600; // Maximum number of loops before reaching cooled state
constexpr int Nover_Qst = 300;  // Maximum number of loops per temperature level
constexpr int Nlimit_Qst = 100; // Number of successful hits per temperature level
constexpr double Tini_Qst = 1.0; // Initial temperature
constexpr double Tfactr_Qst = 0.9; // Cooling rate of temperature change

// Simulated Annealing Parameters for Fst by Locus
constexpr int Tsteps_locFst = 200; // Maximum number of loops before reaching cooled state
constexpr int Nover_locFst = 100;  // Maximum number of loops per temperature level
constexpr int Nlimit_locFst = 30;  // Number of successful hits per temperature level
constexpr double Tini_locFst = 1.0; // Initial temperature
constexpr double Tfactr_locFst = 0.9; // Cooling rate of temperature change