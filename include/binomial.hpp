#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

const int MAX_VAL = 30;
inline long long BINOM[MAX_VAL][MAX_VAL];

inline void initBinomial() {
    for (int i = 0; i < MAX_VAL; ++i) {
        BINOM[i][0] = 1;
        for (int j = 1; j <= i; ++j) {
            BINOM[i][j] = BINOM[i-1][j-1] + BINOM[i-1][j];
        }
    }
}

#endif