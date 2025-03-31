#pragma once
#include <vector>
using namespace std;

// inner dot

template <class T>
T dot(const vector<T> &a, const vector<T> &b) {
    T ret = 0;
    int d = a.size();
    for (int i = 0; i < d; i++) ret += a[i] * b[i];

    return ret;
}

// rot^h(v)

vector<int> rot(const vector<int> &v, int n, int m, int h) {
    int d = n * m;
    vector<int> rot_v(d);
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < n; l++) {
            if (l + h < n) {
                rot_v[n * i + l + h] = v[l];
            } else {
                rot_v[n * i + l + h - n] = -v[l];
            }
        }
    }

    return rot_v;
}

// GSO[i] = b^*_i
// c[i] = ||b^*_i||^2
// output π_k(v)

vector<double> projection(const vector<vector<double>> &GSO, const vector<double> &c, int k, const vector<int> &v) {
    int d = GSO.size();
    vector<double> ret(d);
    for (int i = 0; i < d; i++) ret[i] = v[i];
    for (int i = 0; i < k; i++) {
        double coef = 0.0;
        for (int j = 0; j < d; j++) coef += GSO[i][j] * v[j];
        coef /= c[i];
        for (int j = 0; j < d; j++) ret[j] -= coef * GSO[i][j];
    }
    return ret;
}

// GSO basis (GSO) and GSO coefficients (mu) of lattice basis B

pair<vector<vector<double>>, vector<vector<double>>> calc_GSO(const vector<vector<int>> &B) {
    int d = B.size();
    vector<vector<double>> GSO(d, vector<double>(d, 0.0)), mu(d, vector<double>(d, 0.0));
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) GSO[i][j] = B[i][j];
    }
    for (int i = 0; i < d; i++) {
        for (int k = 0; k < i; k++) {
            mu[i][k] = dot(GSO[i], GSO[k]) / dot(GSO[k], GSO[k]);
            for (int j = 0; j < d; j++) GSO[i][j] -= mu[i][k] * GSO[k][j];
        }
    }
    return make_pair(GSO, mu);
}
