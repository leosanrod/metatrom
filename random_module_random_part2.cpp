#include <cmath>
#include <stdexcept>
#include <random>

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

// Function: random_gamma
// Generates a random number from a gamma distribution.
double random_gamma(double shape) {
    if (shape <= zero) {
        throw std::invalid_argument("Shape parameter must be positive.");
    }

    if (shape > one) {
        return random_gamma1(shape);
    } else if (shape < one) {
        return random_gamma2(shape);
    } else {
        return random_exponential();
    }
}

// Helper function for gamma distribution when shape > 1
double random_gamma1(double shape) {
    double d = shape - one / 3.0;
    double c = one / std::sqrt(9.0 * d);

    double u, v, x;

    while (true) {
        do {
            x = random_normal();
            v = std::pow(one + c * x, 3);
        } while (v <= zero);

        u = generate_uniform();

        if (u < one - 0.0331 * x * x * x * x) {
            return d * v;
        }

        if (std::log(u) < 0.5 * x * x + d * (one - v + std::log(v))) {
            return d * v;
        }
    }
}

// Helper function for gamma distribution when shape < 1
double random_gamma2(double shape) {
    double a = one - shape;
    double p = a / (a + shape * std::exp(-a));
    double c = one / shape;
    double uf = p * std::pow(std::numeric_limits<double>::min() / a, shape);
    double vr = one - std::numeric_limits<double>::min();
    double d = a * std::log(a);

    double r, x, w;

    while (true) {
        r = generate_uniform();
        if (r >= vr) continue;

        if (r > p) {
            x = a - std::log((one - r) / (one - p));
            w = a * std::log(x) - d;
        } else if (r > uf) {
            x = a * std::pow(r / p, c);
            w = x;
        } else {
            return zero;
        }

        r = generate_uniform();
        if (one - r <= w && r > zero) {
            if (r * (w + one) >= one) continue;
            if (-std::log(r) <= w) continue;
        }

        return x;
    }
}