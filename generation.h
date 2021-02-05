#include "specimen.h"
#include "definitions.h"

#ifndef EVOLBIO_GENERATION_H
#define EVOLBIO_GENERATION_H

namespace dynamics {
    class generation {
    public:
        ll n;
        std::vector<specimen> population;

        generation(ll n, std::array<ld, LOCI> frequency) {
            this->n = n;
            assert(!(n % 2));
            F(i, n) population.push_back(specimen::rand(frequency));
        }

        explicit generation(ll n) {
            this->n = n;
            assert(!(n % 2));
        };

        generation randomMating() {
            std::vector<ll> males, females;
            F(i, n)
                if (MALE(population[i])) males.push_back(i);
                else females.push_back(i);

            assert(males.size());
            assert(females.size());

            std::uniform_int_distribution<ll> male_dist(0, males.size() - 1), female_dist(0, females.size() - 1);

            generation next(n);
            F(i, n) {
                specimen male_par = population[males[male_dist(generator)]];
                specimen female_par = population[females[female_dist(generator)]];
                next.population.push_back(specimen::offspring({male_par, female_par}));
            }

            return next;
        }

        std::array<ld, LOCI> frequency_analysis() {
            std::array<ld, LOCI> frequency{};
            F(i, n) F(j, LOCI) F(k, 2) frequency[j] += ((ld) population[i].gene[j][k]) / (n * 2) * (2 - !!j);
            return frequency;
        }

        friend std::ostream &operator<<(std::ostream &cout, const generation &g) {
            cout << "[";
            F(i, g.n) {
                if (i) cout << ", ";
                cout << g.population[i];
            }
            return cout << "]";
        }
    }; // class generation
} // namespace dynamics


#endif //EVOLBIO_GENERATION_H
