#ifndef EVOLBIO_STRUCTURES_H
#define EVOLBIO_STRUCTURES_H
#include <utility>

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

        explicit ubtree(ll leaves) {
            n = leaves;
            m = 2 * n - 2;
            c.resize(m);
            par.resize(m);
            par[0] = -1;
        }

        [[nodiscard]] ubtree empty_copy() const {
            ubtree s(n);
            s.label.resize(this->label.size());
            return s;
        }

        explicit ubtree(const vl& choice) : ubtree(choice.size() + 2) {
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

    enum codon { A, C, T, G };

    class genome : public std::vector<codon> {
    public:
        genome(const std::string& s)  {
            this->resize(s.size());
            std::unordered_map<char, codon> label{{'A', A}, {'C', C}, {'T', T}, {'G', G}};
            F(i,s.size()) (*this)[i] = label[s[i]];
        }
        genome(const char s[]) : genome((std::string)s) {};
    };

    class phylogeny : public ubtree {
    public:
        std::vector<genome> g;
        size_t genome_length;
        phylogeny(ubtree t, std::vector<genome> h) : ubtree(std::move(t)), g(std::move(h)) {
            assert(g.size() == n);
            genome_length = g[0].size();
            F(i, n) assert(g[i].size() == genome_length); // all genomes must be same length
        };
    };

} // namespace structures

#endif //EVOLBIO_STRUCTURES_H
