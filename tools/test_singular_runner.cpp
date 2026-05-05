// test_singular_runner.cpp — Verify SingularRunner C++ infrastructure
#include <iostream>
#include "SingularRunner.hpp"

using namespace std;

int main() {
    bool allPass = true;

    // Test 1: Groebner basis
    cout << "=== Test 1: Groebner basis ===" << endl;
    try {
        vector<string> ideal = {"x1^2-2", "x1*x2-1", "x3"};
        auto gb = SingularRunner::groebnerBasis(ideal, {"x1","x2","x3"}, 179424673);
        cout << "GB: ";
        for (auto& p : gb) cout << p << " ";
        cout << endl << "PASS" << endl;
    } catch (const exception& e) {
        cerr << "FAIL: " << e.what() << endl;
        allPass = false;
    }

    // Test 2: Minimal associated primes
    cout << "\n=== Test 2: minAssGTZ ===" << endl;
    try {
        vector<string> ideal = {"x1^2-2", "x1*x2-1", "x3"};
        auto [primelist, dims] = SingularRunner::minimalAssPrimes(
            ideal, {"x1","x2","x3"}, 179424673, "lp");

        cout << "Found " << primelist.size() << " component(s):";
        for (size_t i = 0; i < primelist.size(); ++i) {
            cout << "\n  dim=" << dims[i] << ": ";
            for (auto& p : primelist[i]) cout << "[" << p << "] ";
        }
        cout << endl << "PASS" << endl;
    } catch (const exception& e) {
        cerr << "FAIL: " << e.what() << endl;
        allPass = false;
    }

    // Test 3: Primary decomposition with lex output
    cout << "\n=== Test 3: primdecGTZE + lex ===" << endl;
    try {
        vector<string> ideal = {"x1^2-2", "x1*x2-1", "x3"};
        auto [primelist, dims] = SingularRunner::primaryDecompLex(
            ideal, {"x1","x2","x3"}, 179424673);

        cout << "Found " << primelist.size() << " component(s):";
        for (size_t i = 0; i < primelist.size(); ++i) {
            cout << "\n  dim=" << dims[i] << ": ";
            for (auto& p : primelist[i]) cout << "[" << p << "] ";
        }
        // Should have 2 zero-dim components
        if (primelist.size() == 2 && dims[0] == 0 && dims[1] == 0) {
            cout << endl << "PASS" << endl;
        } else {
            cout << endl << "FAIL: expected 2 zero-dim components" << endl;
            allPass = false;
        }
    } catch (const exception& e) {
        cerr << "FAIL: " << e.what() << endl;
        allPass = false;
    }

    // Test 4: idealDim
    cout << "\n=== Test 4: dim ===" << endl;
    try {
        vector<string> ideal = {"x1^2-2", "x1*x2-1", "x3"};
        int d = SingularRunner::idealDim(ideal, {"x1","x2","x3"}, 179424673);
        cout << "dim = " << d << " (expected 0)" << endl;
        cout << (d == 0 ? "PASS" : "FAIL") << endl;
        if (d != 0) allPass = false;
    } catch (const exception& e) {
        cerr << "FAIL: " << e.what() << endl;
        allPass = false;
    }

    cout << "\n=== " << (allPass ? "ALL PASS" : "SOME FAILED") << " ===" << endl;
    return allPass ? 0 : 1;
}
