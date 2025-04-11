#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <limits>

// Constants
const double zero = 0.0;
const double one = 1.0;
const double half = 0.5;
const double two = 2.0;
const double vsmall = std::numeric_limits<double>::min();
const double vlarge = std::numeric_limits<double>::max();

// Random number generator
std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

// Helper function to generate uniform random numbers
double generate_uniform() {
    return uniform_dist(generator);
}

// Function: random_normal
// Generates a random number with a normal (Gaussian) distribution.
double random_normal() {
    // Adapted from the algorithm in ACM Transactions on Mathematical Software
    // Algorithm constants
    const double s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472;
    const double r1 = 0.27597, r2 = 0.27846;

    double u, v, x, y, q;

    while (true) {
        // Generate uniform random numbers
        u = generate_uniform();
        v = 1.7156 * (generate_uniform() - half);

        // Evaluate the quadratic form
        x = u - s;
        y = std::abs(v) - t;
        q = x * x + y * (a * y - b * x);

        // Accept P if inside inner ellipse
        if (q < r1) break;

        // Reject P if outside outer ellipse
        if (q > r2) continue;

        // Reject P if outside acceptance region
        if (v * v < -4.0 * std::log(u) * u * u) break;
    }

    // Return ratio of P's coordinates as the normal deviate
    return v / u;
}