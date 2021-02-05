#ifndef EVOLBIO_ALGOS_H
#define EVOLBIO_ALGOS_H
#include "definitions.h"
#include "structures.h"

using namespace structures;

namespace algos {
    vvld random_distance_matrix(ll n) {
        vvld dist(n, vld(n)), pos(n, vld(n));;
        ll dim = 2*n;
        F(i,n) F(j,dim) pos[i][j] = uniform_distribution_0_1(generator) * 5;
        F(i,n) F(j,n) {
            std::vector<ld> diff(pos[i]);
            F(k,dim) diff[k] -= pos[j][k];
            F(k,dim) diff[k] *= diff[k];
            dist[i][j] = sqrtl(std::accumulate(A(diff), 0.l));
        }
        return dist;
    }

    std::vector<ubtree> gen_trees(ll n) {
        std::vector<ubtree> trees;
        vl choices(n - 2);
        choices[0] = 1;
        F(i,n - 3) choices[i + 1] = 2 + choices[i];
        vl curr(n - 2);
        while (1) {
            trees.emplace_back(n, curr);
            bool inc = true;
            curr[n-3]++;
            for (ll i = n - 3; i > 0; i--)
                if (inc and curr[i] == choices[i]) curr[i - 1]++, curr[i] = 0;
                else inc = false;
            if (curr[0]) break;
        }

        return trees;
    }
} // namespace algos

#endif //EVOLBIO_ALGOS_H
