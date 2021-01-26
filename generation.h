#include "specimen.h"
#include "definitions.h"

#ifndef EVOLBIO_GENERATION_H
#define EVOLBIO_GENERATION_H


class generation {
public:
    ll n;
    std::vector<specimen> population;

    generation(ll n, std::array<ld, ALLELE_COUNT> frequency) {
        this->n = n;
        assert(!(n%2));
        F(i,n) population.push_back(specimen::rand(frequency));
    }

    explicit generation(ll n) {
        this->n = n;
        assert(!(n%2));
    };

    generation randomMating() {
        std::vector<specimen> males, females;
        F(i,n)
            if (MALE(population[i])) males.push_back(population[i]);
            else females.push_back(population[i]);

        assert(males.size());
        assert(females.size());

        std::uniform_int_distribution<ll> male_dist(0,males.size() - 1), female_dist(0, females.size() - 1);

        generation next(n);
        F(i,n) {
            specimen male_par = males[male_dist(generator)], female_par = females[female_dist(generator)];
            next.population.push_back(specimen::offspring({male_par, female_par}));
        }

        return next;
    }

    std::array<ld, ALLELE_COUNT> frequency_analysis() {
        std::array<ld, ALLELE_COUNT> frequency{};
        F(i,n) F(j,ALLELE_COUNT) F(k,2) frequency[j] += ((ld) population[i].gene[j][k]) / (n * 2) * (2 - !!j);
        return frequency;
    }

    friend std::ostream& operator<<(std::ostream& cout, const generation &g) {
        cout << "[";
        F(i,g.n) {
            if (i) cout << ", ";
            cout << g.population[i];
        }
        return cout << "]";
    }
};


#endif //EVOLBIO_GENERATION_H
