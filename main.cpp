#include "definitions.h"
#include "generation.h"
#include "structures.h"
#include "algos.h"

int main() {
    auto trees = algos::gen_trees(4);
    std::cout << *trees.begin() << std::endl;
    phylogeny p(*trees.begin(), {"ATC", "ATG", "ACC", "ACG"});
    std::cout << algos::fitch(p) << std::endl;
    return 0;
}
