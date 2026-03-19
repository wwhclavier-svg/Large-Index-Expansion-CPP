#include "LayerRecursionCore.hpp"
#include "Combinatorics.hpp"   // 已包含在头文件中，但这里再包含一次确保可见

namespace LayerRecursionCore {

void equationVariable(std::vector<std::array<int,4>>& eqnvar,
                      int order, int level,
                      std::vector<int>& seed,
                      int nb, int ne,
                      const std::vector<std::array<int,4>>& indepSet)
{ 
    eqnvar.clear();
    int ncurr = lastNonZero(seed);
    for(int j = max(ncurr,0); j < ne; ++j) { 
        for(int i = 0; i < nb; ++i) {
            eqnvar.push_back({order, level+1, getIndexOffSet(level,seed,1,j),i});
        }
    }
    for(auto &indep : indepSet) { eqnvar.push_back(indep); }
}

} // namespace LayerRecursionCore