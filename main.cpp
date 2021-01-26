#include "definitions.h"
#include "generation.h"

int main() {
    generation g(10000, {0.5, 0.1, 0.4, 0.3});
    D(g.frequency_analysis())
    g = g.randomMating();
//    std::cout << g << std::endl;
    D(g.frequency_analysis())
    return 0;
}
