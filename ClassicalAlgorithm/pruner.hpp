#pragma once
#include <cmath>
#include <random>
#include <vector>
using namespace std;

const double pi = acos(-1.0);
const double e = exp(1.0);

// estimated number of nodes N'(R) in enumeration tree

double estimated_nodes(const vector<double> &b_norm, const vector<double> &R) {
    int beta_prime = R.size();
    vector<vector<double>> comb(beta_prime + 1, vector<double>(beta_prime + 1, 0.0));
    comb[0][0] = 1.0;
    for (int i = 0; i < beta_prime; i++) {
        for (int j = 0; j <= i; j++) {
            comb[i + 1][j] += comb[i][j];
            comb[i + 1][j + 1] += comb[i][j];
        }
    }
    vector<double> nodes(beta_prime);
    double nodes_sum = 0.0;

    for (int k = 1; k <= beta_prime / 2; k++) {
        vector<double> b(k);
        for (int i = 0; i < k; i++) b[i] = R[2 * i] / R[2 * k - 2];
        vector<double> F(k + 1, 0.0);
        F[k] = 1.0;
        for (int j = 0; j < k; j++) {
            double h = F[k - j];
            for (int i = 0; i <= k - j; i++) {
                F[k - j - i] -= h * comb[k - j][i];
                h *= -b[j];
            }
        }
        double coef = F[0];
        for (int i = 0; i < k; i++) {
            coef *= R[2 * k - 2] * pi;
            coef /= (i + 1) * b_norm[2 * i] * b_norm[2 * i + 1];
        }
        nodes[2 * k - 2] = coef;
        nodes[2 * k - 1] = coef;
        nodes_sum += coef * 2;
    }

    return nodes_sum;
};

// estimated probability and loops (P'(R), E[loop(R)])

pair<double, double> estimated_probability_and_loops(const vector<vector<vector<double>>> &norms, const vector<double> &R, int r) {
    int T = norms.size();
    int beta_prime = R.size();
    int succ = 0, loops = 0;
    for (int t = 0; t < T; t++) {
        for (int h = 0; h < r; h++) {
            bool ok = true;
            loops++;
            for (int i = 1; i < beta_prime; i += 2) {
                if (norms[t][h][i] > R[i]) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                succ++;
                break;
            }
        }
    }
    return make_pair(double(succ) / double(T), double(loops) / double(T));
}

// estimated probability and cost (P'(R), N'(R)*E[loop(R)])

pair<double, double> estimated_probability_and_cost(const vector<double> &b_norm, const vector<vector<vector<double>>> &norms, const vector<double> &R, int r) {
    double nodes = estimated_nodes(b_norm, R);
    auto [prob, loops] = estimated_probability_and_loops(norms, R, r);
    return make_pair(prob, nodes * loops);
}

// minimize cost N'(R)*E[loop(R)] under success probability P'(R)>=p_0 
// R holds in the form of　R ^2

vector<double> simulated_annealing(const vector<double> &b_norm, int beta_prime, double sigma, const vector<vector<vector<double>>> &norms, double t0, double t1, int s, double p0, int r) {
    vector<double> R(beta_prime, 0.0);
    while (estimated_probability_and_loops(norms, R, r).first < p0) {
        for (int i = 0; i < beta_prime / 2; i++) {
            R[2 * i] += 0.01 * 2 * (i + 1) * sigma * sigma;
            R[2 * i + 1] = R[2 * i];
        }
    }
    double score = 1.0;
    double initial_cost = estimated_probability_and_cost(b_norm, norms, R, r).second;

    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_real_distribution<> dist1(-1.0 * sigma * sigma, 1.0 * sigma * sigma), dist2(0.0, 1.0);
    uniform_int_distribution<> dist3(0, beta_prime / 2 - 1);
    double bound = 1.0 * beta_prime * sigma * sigma;

    for (int i = 0; i <= s; i++) {
        double t = pow(t0, 1.0 - double(i) / double(s)) * pow(t1, double(i) / double(s));
        int j = dist3(engine);
        double x = dist1(engine);
        double pert = R[2 * j] + x;
        pert = max(pert, j > 0 ? R[2 * (j - 1)] : R[0]);
        pert = min(pert, j < beta_prime / 2 - 1 ? R[2 * (j + 1)] : bound);

        vector<double> next_R = R;
        next_R[2 * j] = pert;
        next_R[2 * j + 1] = pert;

        auto [prob, cost] = estimated_probability_and_cost(b_norm, norms, next_R, r);
        if (prob >= p0) {
            double next_score = cost / initial_cost;
            double p = exp((score - next_score) / t);
            if (dist2(engine) <= p) {
                R = next_R;
                score = next_score;
            }
        }
    }

    return R;
}
