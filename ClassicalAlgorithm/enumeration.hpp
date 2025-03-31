#pragma once
#include <math.h>
#include <vector>
using namespace std;

// c[i] = ||b^*_i||^2
// v[i] = Σ[j=0,...,d-1] mu_v[j] b^*_j
// i < left: Babai's nearest plane algorithm
// if find:true, else:false

bool enumerate(const vector<vector<double>> &mu, const vector<double> &c, const vector<double> &R, const vector<double> &mu_v, vector<int> &x, int left, int i, double norm = 0.0) {
    if (i == -1) return true;
    int d = mu.size();
    double center = 0.0;
    for (int j = i + 1; j < d; j++) center -= mu[j][i] * x[j];
    center += mu_v[i];
    int mid = round(center);

    {
        double next_norm = norm + (mid - center) * (mid - center) * c[i];
        if (next_norm > R[i]) return false;
        x[i] = mid;
        if (enumerate(mu, c, R, mu_v, x, left, i - 1, next_norm)) return true;
    }

    if (i >= left) {
        for (int t = 1;; t++) {
            int dir = (mid < center ? t : -t);
            {
                double next_norm = norm + (mid + dir - center) * (mid + dir - center) * c[i];
                if (next_norm > R[i]) return false;
                x[i] = mid + dir;
                if ((enumerate(mu, c, R, mu_v, x, left, i - 1, next_norm))) return true;
            }
            {
                double next_norm = norm + (mid - dir - center) * (mid - dir - center) * c[i];
                if (next_norm > R[i]) return false;
                x[i] = mid - dir;
                if ((enumerate(mu, c, R, mu_v, x, left, i - 1, next_norm))) return true;
            }
        }
    }

    return false;
}
