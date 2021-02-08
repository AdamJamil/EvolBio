#include "definitions.h"
#include "generation.h"
#include "structures.h"
#include "algos.h"

int main() {
    ubtree t(5, {0, 2, 2});
    auto trees = algos::gen_trees(8);
    tr(x, trees) std::cout << x << std::endl;
    return 0;
}
