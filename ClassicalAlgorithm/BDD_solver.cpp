#include "enumeration.hpp"
#include "pruner.hpp"
#include "utility.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace std;

// c[i] = ||b^*_i||^2
// pruning parameter (len:d, Reverse order)

vector<double> pruning_optimize(const vector<vector<int>> &B, const vector<vector<double>> &GSO, const vector<double> &c, int n, int m, double sigma, int beta_prime, const vector<vector<int>> &e_sample, double p0, double t0 = 0.01, double t1 = 0.0001, int s = 20000, int r = -1) {
    if (r == -1) r = n;
    int d = B.size();
    vector<double> b_norm(beta_prime);
    for (int i = 0; i < beta_prime; i++) b_norm[i] = sqrt(c[d - 1 - i]);
    int T = e_sample.size();

    vector<vector<vector<double>>> norms(T, vector<vector<double>>(n, vector<double>(beta_prime, 0.0)));
    for (int t = 0; t < T; t++) {
        for (int h = 0; h < n; h++) {
            vector<int> rot_e = rot(e_sample[t], n, m, h);
            for (int i = 0; i < beta_prime; i++) {
                for (int j = 0; j < d; j++) norms[t][h][i] += GSO[d - 1 - i][j] * rot_e[j];
                norms[t][h][i] /= b_norm[i];
                norms[t][h][i] *= norms[t][h][i];
                if (i > 0) norms[t][h][i] += norms[t][h][i - 1];
            }
        }
    }

    vector<double> R = simulated_annealing(b_norm, beta_prime, sigma, norms, t0, t1, s, p0, r);
    int left = d - beta_prime, right = d;
    vector<double> R_all(d);
    for (int i = 0; i < left; i++) R_all[i] = 1.5 * d * sigma * sigma;
    for (int i = left; i < d; i++) R_all[i] = R[d - 1 - i];

    return R_all;
}

// enumeration algorithm;target vector v, parameter R + Babai's nearest plane algorithm
// output difference between lattice vector L(B) close enough to v and v (output 0 if not find) 

vector<int> sol(const vector<vector<int>> &B, const vector<vector<double>> &GSO, const vector<vector<double>> &mu, const vector<double> &c, vector<int> v, int beta_prime, const vector<double> &R) {
    int d = B.size();
    int left = d - beta_prime, right = d;

    vector<double> mu_v(d, 0);
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) mu_v[i] += GSO[i][j] * v[j];
        mu_v[i] /= c[i];
    }

    vector<int> x(d, 0);
    if (enumerate(mu, c, R, mu_v, x, left, right - 1, 0.0)) {
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) v[j] -= x[i] * B[i][j];
        }
        return v;
    }

    return vector<int>(d, 0);
}

PYBIND11_MODULE(BDD_solver, m) {
    m.doc() = "BDD solver for proporsed algorithm";
    m.def("pruning_optimize", &pruning_optimize);
    m.def("sol", &sol);
}
