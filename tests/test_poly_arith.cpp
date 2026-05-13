// Quick unit test for PolyArith.hpp
#include "../include/PolyArith.hpp"
#include <iostream>
#include <cassert>

using namespace PolyArith;

int main() {
    const int64_t mod = 179424673;
    int nVars = 3;

    // Test 1: polyAdd
    {
        auto a = makeMonomial({1,0,0}, 3, mod);  // 3*x1
        auto b = makeMonomial({1,0,0}, 5, mod);  // 5*x1
        auto s = polyAdd(a, b, mod);
        assert(s.size() == 1);
        assert(s[0].coeff == 8);
        assert(s[0].exps == std::vector<int>({1,0,0}));
        std::cout << "PASS: polyAdd\n";
    }

    // Test 2: polyAdd with cancellation
    {
        auto a = makeMonomial({1,0,0}, 5, mod);
        auto b = makeMonomial({1,0,0}, -5, mod);
        auto s = polyAdd(a, b, mod);
        assert(s.empty());
        std::cout << "PASS: polyAdd cancellation\n";
    }

    // Test 3: polyMul
    {
        auto a = makeMonomial({1,0,0}, 2, mod);  // 2*x1
        auto b = makeMonomial({0,1,0}, 3, mod);  // 3*x2
        auto p = polyMul(a, b, mod);
        assert(p.size() == 1);
        assert(p[0].coeff == 6);
        assert(p[0].exps == std::vector<int>({1,1,0}));
        std::cout << "PASS: polyMul\n";
    }

    // Test 4: polySubstitute — replace x1 → 2*x2 + 1
    {
        // poly = 3*x1, replace x1 → 2*x2 + 1 → result = 3*(2*x2+1) = 6*x2 + 3
        auto p = makeMonomial({1,0,0}, 3, mod);
        Polynomial repl;
        repl.push_back({std::vector<int>(nVars, 0), int64_t(1)});     // 1
        repl.push_back({std::vector<int>({0,1,0}), int64_t(2)});       // 2*x2
        auto s = polySubstitute(p, 0, repl, mod);
        assert(s.size() == 2);
        std::cout << "PASS: polySubstitute\n";
    }

    // Test 5: monomialRulesPower
    {
        Polynomial p;
        p.push_back({std::vector<int>({1,0}), int64_t(3)});
        p.push_back({std::vector<int>({0,1}), int64_t(5)});
        p.push_back({std::vector<int>({1,0}), int64_t(2)}); // same as first
        auto rules = monomialRulesPower(p);
        assert((rules[{1,0}] == 5));  // 3+2
        assert((rules[{0,1}] == 5));
        std::cout << "PASS: monomialRulesPower\n";
    }

    // Test 6: polyToSingularString
    {
        auto p = makeMonomial({2,1,0}, 3, mod);
        std::string s = polyToSingularString(p, {"A1","A2","A3"});
        std::cout << "  Singular: " << s << std::endl;
        // Should be something like "3*A1^2*A2"
    }

    // Test 7: canonicalize with mod reduction
    {
        Polynomial p;
        p.push_back({std::vector<int>({1,0}), mod - 3}); // -3 mod p
        p.push_back({std::vector<int>({1,0}), int64_t(3)});
        canonicalize(p, mod);
        assert(p.empty()); // should cancel
        std::cout << "PASS: canonicalize mod\n";
    }

    std::cout << "\nAll PolyArith tests passed!\n";
    return 0;
}
