#include "layerRecursion.hpp"

// 实现非模板函数
using namespace std;

int lastNonZero(const vector<int>& list) {
    int i = static_cast<int>(list.size()) - 1;
    while (i >= 0 && list[i] == 0) {
        --i;
    }
    return i;
}

void findPartitions(int level, int n, int current_index, vector<int>& buffer, vector<vector<int>>& current_level_results) {
    if (current_index == n) {
        if (level == 0) {
            current_level_results.push_back(buffer);
        }
        return;
    }
    for (int i = level; i >= 0; --i) {
        buffer[current_index] = i;
        findPartitions(level - i, n, current_index + 1, buffer, current_level_results);
    }
}

vector<vector<vector<int>>> seedGenerator(int maxlevel, int n) {
    vector<vector<vector<int>>> all_results(maxlevel + 1);
    vector<int> buffer(n);
    for (int L = 0; L <= maxlevel; ++L) {
        findPartitions(L, n, 0, buffer, all_results[L]);
    }
    return all_results;
}

int getIndexOffSet(int old_level, vector<int> &seed, int k, int i) {
    seed[i] += k;
    int new_level = old_level + k;
    int idx = getIndex(seed, new_level);
    seed[i] -= k;
    return idx;
}

bool notin(const std::vector<int>& list, int a) {
    return std::find(list.begin(), list.end(), a) == list.end();
}

void equationVariable(vector<std::array<int,4>>& eqnvar, int order, int level, vector<int> seed, int nb, int ne, vector<std::array<int,4>>& indepSet) 
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