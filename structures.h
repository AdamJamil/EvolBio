#ifndef EVOLBIO_STRUCTURES_H
#define EVOLBIO_STRUCTURES_H
#include "definitions.h"

namespace structures {
    class ubtree {
    private:
        vl label;
    public:
        // n is number of leaves (including root !), m is total vertices
        ll n{}, m{};
        // adjacency list representation, but only including children
        // first idx MUST represent first x_1
        // after re-rooting, this is relaxed further to the following:
            // in order to print, the root must be a leaf.
            // we can relax even this by changing how print works, so that we print a label instead of the actual idx
        std::vector<std::unordered_set<ll>> c;
        // every non-root vertex has a parent
        vl par;
        // label[i] is the vertex where the edge above it is edge i
            // unless generated from combo constructor, this will be empty!

        ubtree() = default;

        ubtree empty_copy() const {
            ubtree s;
            s.c.resize(this->c.size());
            s.par.resize(this->par.size());
            s.label.resize(this->label.size());
            s.n = this->n;
            s.m = this->m;

            return s;
        }

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

        ubtree reroot_at(ll x) {
            vl fringe{x};
            sl seen{x};
            ubtree s = empty_copy();

            while (!fringe.empty()) {
                vl next;
                for (ll u : fringe) {
                    for (ll v : c[u]) if (!seen.count(v)) {
                            seen.insert(v);
                            next.push_back(v);
                            s.par[v] = u;
                            s.c[u].insert(v);
                        }
                    if (par[u] != -1 and !seen.count(par[u])) {
                        ll v = par[u];
                        seen.insert(v);
                        next.push_back(v);
                        s.par[v] = u;
                        s.c[u].insert(v);
                    }
                }
                fringe = next;
            }

            return s;
        }

        friend std::ostream& operator<<(std::ostream& cout, const ubtree& obj) {
            std::set<pl> dfs{{0,0}}; // store (depth, val)
            while (!dfs.empty()) {
                auto ptr = dfs.begin(); // have to do it backwards - can't easily erase .rbegin()
                dfs.erase(ptr);
                ll u = ptr->second, d = ptr->first;
                for (ll v : obj.c[u]) dfs.insert({d - 1, v});
                F(i,-d) cout << "  ";
                cout << "|-" << ((u < obj.n) ? std::to_string(u + 1) : "-") << '\n';
            }
            return cout;
        }
    }; // class ubtree

    class genome {
    private:
        // stores two consecutive bits for one codon
        // 0 = A
        // 1 = T
        // 2 = C
        // 3 = G
        bool *codon;
        size_t sz;
    public:
        explicit genome(std::string s) {
            sz = s.size();
            codon = (bool *) malloc(s.size() * 2 * sizeof(bool));
            F(i,s.size())
                if (s[i] == 'A')      codon[2*i] = 0, codon[2*i + 1] = 0;
                else if (s[i] == 'T') codon[2*i] = 0, codon[2*i + 1] = 1;
                else if (s[i] == 'C') codon[2*i] = 1, codon[2*i + 1] = 0;
                else if (s[i] == 'G') codon[2*i] = 1, codon[2*i + 1] = 1;
        }

        int operator[](int idx) {
            return codon[2*idx] + 2 * codon[2*idx + 1];
        }
    };

    class phylogeny : public ubtree {
    public:
        genome g;
        explicit phylogeny(genome h) : ubtree(), g(h) {};
    };

} // namespace structures

#endif //EVOLBIO_STRUCTURES_H
