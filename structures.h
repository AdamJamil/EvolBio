#ifndef EVOLBIO_STRUCTURES_H
#define EVOLBIO_STRUCTURES_H
#include "definitions.h"

namespace structures {
    class ubtree {
    public:
        // n is number of leaves (including root !), m is total vertices
        ll n, m;
        // adjacency list representation, but only including children
        // first idx MUST represent first x_1
        std::vector<std::unordered_set<ll>> c;
        // every non-root vertex has a parent
        vl par;
        // label[i] is the vertex where the edge above it is edge i
        vl label;

        ubtree(ll leaves, const vl& choice) {
            n = leaves;
            m = 2 * n - 2;
            c.resize(m);
            par.resize(m);
            c[0].insert(1);      // children of root
            par[0] = -1;            // parent of root
            par[1] = 0;             // parent of second leaf
            label.emplace_back(1);  // edge 0 is above vertex 1; this is invariant under adding vertices!

            for (ll j = 3; j <= n; j++) {
                ll p = choice[j - 3];
                ll v = label[p];
                ll u = par[v];
                ll x = j + (n - 3);
                ll w = j - 1;
                label.emplace_back(x); // inner node gets edge k + 1
                label.emplace_back(w); // new leaf gets edge k + 2
                c[u].erase(v);         // delete old edge from u to v
                c[u].insert(x);        // connect u to x
                par[x] = u;
                c[x].insert(v);        // connect x to new children
                c[x].insert(w);
                par[v] = x;
                par[w] = x;
            }
        }

        friend std::ostream& operator<<(std::ostream& cout, const ubtree& obj) {
            std::set<pl> dfs{{0,0}}; // store (depth, val)
            while (dfs.size()) {
                auto ptr = dfs.begin(); // have to do it backwards - can't easily erase .rbegin()
                dfs.erase(ptr);
                ll u = ptr->second, d = ptr->first;
                for (ll v : obj.c[u]) dfs.insert({d - 1, v});
                F(i,-d) cout << "  ";
                cout << "|-" << ((u < obj.n) ? std::to_string(u + 1) : "-") << '\n';
            }
            return cout;
        }
    };

} // namespace structures

#endif //EVOLBIO_STRUCTURES_H
