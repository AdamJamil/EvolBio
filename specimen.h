#include "definitions.h"

#ifndef EVOLBIO_SPECIMEN_H
#define EVOLBIO_SPECIMEN_H


class specimen  {
public:
    std::array<std::array<ll, 2>, ALLELE_COUNT> gene{};

    static specimen offspring(std::vector<specimen> parents) {
        specimen ret{};
        F(i,ALLELE_COUNT) F(j,2) {
            ll res = coin_flip(generator);
            ret.gene[i][j] = parents[j].gene[i][res];
            if (parents[j].gene[i][res] > 1) {
                P(j << " " << i << " " << res)
                P(parents[j])
                P(parents[j].gene[i][res])
            }
//            assert(ret.gene[i][j] <= 1);
//            assert(ret.gene[i][j] >= 0);
        }
        return ret;
    }

    static specimen rand(std::array<ld, ALLELE_COUNT> frequency) {
        specimen ret{};
        F(i,ALLELE_COUNT) F(j,2) {
            if (i) ret.gene[i][j] = uniform_distribution_0_1(generator) < frequency[i];
            else ret.gene[i][j] = j * (uniform_distribution_0_1(generator) < frequency[i]);
        }
        return ret;
    }

    friend std::ostream& operator<<(std::ostream& cout, const specimen &g)
    {
        cout << "{";
        F(i,ALLELE_COUNT) {
            if (i) cout << ", " << g.gene[i][0] << g.gene[i][1];
            else cout << ((MALE(g)) ? "MALE" : "FEMALE");
        }
        return cout << "}";
    }
};


#endif //EVOLBIO_SPECIMEN_H
