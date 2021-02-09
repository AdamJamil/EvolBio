#include "definitions.h"
#include "generation.h"
#include "structures.h"
#include "algos.h"

int main() {
    sl vals{1};
    vvld dist = {{0,2,5,5}, {2,0,5,5}, {5,5,0,2}, {5,5,2,0}};
    std::cout << algos::UPGMA(dist);
    return 0;
}
