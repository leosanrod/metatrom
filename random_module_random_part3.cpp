#include <cmath>
#include <random>
#include <stdexcept>

// Constants
const double zero = 0.0;
const double one = 1.0;

// Random number generator
std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

// Helper function to generate uniform random numbers
double generate_uniform() {
    return uniform_dist(generator);
}

// Function: random_exponential
// Generates a random number from an exponential distribution.
double random_exponential() {
    double r;

    do {
        r = generate_uniform();
    } while (r <= zero);

    return -std::log(r);
}

// Function: random_binomial
// Generates a random number from a binomial distribution.
int random_binomial(int n, double p) {
    if (p < zero || p > one) {
        throw std::invalid_argument("Probability must be between 0 and 1.");
    }
    if (n < 0) {
        throw std::invalid_argument("Number of trials must be non-negative.");
    }

    std::binomial_distribution<int> binomial_dist(n, p);
    return binomial_dist(generator);
}

// Function: random_poisson
// Generates a random number from a Poisson distribution.
int random_poisson(double lambda) {
    if (lambda <= zero) {
        throw std::invalid_argument("Lambda must be positive.");
    }

    std::poisson_distribution<int> poisson_dist(lambda);
    return poisson_dist(generator);
}

// Function: random_multinomial
// Generates a set of random numbers from a multinomial distribution.
std::vector<int> random_multinomial(int n, const std::vector<double>& probabilities) {
    if (n < 0) {
        throw std::invalid_argument("Number of trials must be non-negative.");
    }

    double sum_probabilities = 0.0;
    for (double p : probabilities) {
        if (p < zero) {
            throw std::invalid_argument("Probabilities must be non-negative.");
        }
        sum_probabilities += p;
    }
    if (std::abs(sum_probabilities - one) > 1e-6) {
        throw std::invalid_argument("Probabilities must sum to 1.");
    }

    std::vector<int> counts(probabilities.size(), 0);
    double remaining_probability = one;

    for (size_t i = 0; i < probabilities.size() - 1; ++i) {
        std::binomial_distribution<int> binomial_dist(n, probabilities[i] / remaining_probability);
        counts[i] = binomial_dist(generator);
        n -= counts[i];
        remaining_probability -= probabilities[i];
    }

    counts.back() = n;
    return counts;
}

// Function: random_uniform
// Generates a random number from a uniform distribution between a and b.
double random_uniform(double a, double b) {
    if (a >= b) {
        throw std::invalid_argument("Lower bound must be less than upper bound.");
    }

    std::uniform_real_distribution<double> uniform_dist(a, b);
    return uniform_dist(generator);
}

// Function: random_integer
// Generates a random integer between a and b (inclusive).
int random_integer(int a, int b) {
    if (a > b) {
        throw std::invalid_argument("Lower bound must be less than or equal to upper bound.");
    }

    std::uniform_int_distribution<int> uniform_dist(a, b);
    return uniform_dist(generator);
}

// Function: set_random_seed
// Sets the random seed for reproducibility.
void set_random_seed(int seed) {
    generator.seed(seed);
}