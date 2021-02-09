#ifndef EVOLBIO_ALGOS_H
#define EVOLBIO_ALGOS_H
#include "definitions.h"
#include "structures.h"

using namespace structures;

namespace algos {

    vvld random_distance_matrix(ll n) {
        vvld dist(n, vld(n)), pos(n, vld(2 * n));;
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

    // clustering algos

    ld UPGMA_dist(ll u, ll v, vvld &dist, std::unordered_map<ll, sl> &clusters) {
        sl c1 = clusters[u], c2 = clusters[v]; // optimization - just pass in the right clusters
        ld tot = 0;
        for (ll x : c1) for (ll y : c2) tot += dist[x][y];
        tot /= c1.size() * c2.size();
        return tot;
    }

    ubtree UPGMA(vvld dist) {
        ll n = dist.size();
        std::set<std::pair<ld, pl>> choices;
        F(i,n) FS(j,i+1,n) choices.insert({dist[i][j], {i, j}});
        graph gr(2 * (n - 1));
        std::unordered_map<ll, sl> clusters;
        F(i,n) clusters[i] = {i};
        sl rem;
        F(i,n) rem.insert(i);
        F(i,n-2) {
            auto ptr = choices.begin();
            for (ll x : {ptr->Y.X, ptr->Y.Y}) {
                gr.join(x, i + n);
                clusters[i + n].insert(A(clusters[x]));
                TR(j, rem) {
                    choices.erase({dist[j][x], {j, x}});
                    choices.erase({dist[x][j], {x, j}});
                }
                rem.erase(x);
            }
            rem.insert(i + n);
            dist.emplace_back();
            F(j,i+n) dist[i+n].push_back(UPGMA_dist(i + n, j, dist, clusters));
            clusters[i+n].insert(i+n); // has to be inserted after UPGMA_dist call
            F(j,i+n) dist[j].push_back(dist[i+n][j]);
            dist[i+n].push_back(0);
            F(j,i+n) if (!clusters[i+n].count(j)) choices.insert({dist[j][i+n], {j, i+n}});
        }
        auto ptr = choices.begin();
        ll u = ptr->second.first, v = ptr->second.second;
        gr.join(u, v);

        return gr.to_ubtree();
    }

    std::vector<ubtree> gen_trees(ll n) {
        std::vector<ubtree> trees;
        vl choices(n - 2);
        choices[0] = 1;
        F(i,n - 3) choices[i + 1] = 2 + choices[i];
        vl curr(n - 2);
        while (true) {
            trees.emplace_back(curr);
            bool inc = true;
            curr[n-3]++;
            for (ll i = n - 3; i > 0; i--)
                if (inc and curr[i] == choices[i]) curr[i - 1]++, curr[i] = 0;
                else inc = false;
            if (curr[0]) break;
        }

        return trees;
    }

    void fitch_dfs(vpl &dp, phylogeny &p, ll u, ll col) {
        for (ll v : p.c[u])
            fitch_dfs(dp, p, v, col);
        if (p.c[u].empty()) { // leaf
            dp[u] = {0, codon_to_mask[p.g[u][col]]};
        } else if (p.c[u].size() == 1) { // root
            ll child = *p.c[u].begin();
            dp[u] = {dp[child].X + !(codon_to_mask[p.g[u][col]] & (ull)dp[child].Y), codon_to_mask[p.g[u][col]]};
        } else { // two children
            ll left = *p.c[u].begin(), right = *(++p.c[u].begin());
            ull left_mask = dp[left].Y, right_mask = dp[right].Y;
            dp[u].X = dp[left].X + dp[right].X + !(right_mask & left_mask);
            if (!(right_mask & left_mask)) dp[u].Y = right_mask | left_mask;
            else dp[u].Y = right_mask & left_mask;
        }
    }

    ll fitch(phylogeny p) {
        ll ans = 0;
        vpl dp(p.m);
        F(alignment_column,p.genome_length) {
            if (alignment_column) std::fill(A(dp), pl{0, 0});
            fitch_dfs(dp, p, 0, alignment_column);
            ans += dp[0].X;
        }
        return ans;
    }

} // namespace algos

#endif //EVOLBIO_ALGOS_H
