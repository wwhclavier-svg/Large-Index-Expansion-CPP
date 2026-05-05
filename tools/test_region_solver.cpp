// test_region_solver.cpp — Verify RegionSolver Singular pipeline
#include <iostream>
#include "RegionSolver.hpp"
#include "SingularRunner.hpp"

using namespace std;

int main() {
    bool allPass = true;

    // Test: simple nonlinear system over finite field
    // Equations in A/B form (mimicking a simple IBP case)
    // A1 = propagator denominator variable
    // B1 = 1/A1
    //
    // System: A1^2 - 2 = 0, A1*B1 - 1 = 0
    // This is the simplest possible case: one propagator, characteristic 0 sector
    int64_t p = 179424673;
    int ne = 1;

    cout << "=== RegionSolver: Simple nonlinear system ===" << endl;
    cout << "Modulus: " << p << ", ne=" << ne << endl;

    // Build A/B equations
    vector<string> eqs;
    eqs.push_back("A1^2-2");     // A1^2 = 2 → quotient ring Z_p[A1]/(A1^2-2)
    eqs.push_back("A1*B1-1");    // B1 = 1/A1 = A1/2 (in the field)

    // Solve
    try {
        vector<int> sector = {1}; // top sector, active
        auto regions = RegionSolver::solveRegion(eqs, sector, ne, p);

        cout << "\nFound " << regions.size() << " region(s)" << endl;
        for (size_t i = 0; i < regions.size(); ++i) {
            auto& reg = regions[i];
            cout << "Region " << i << ": nb=" << reg.nb
                 << " VarIndep=" << reg.VarIndep.size()
                 << " VarDep=" << reg.VarDep.size() << endl;

            cout << "  Variable degrees: ";
            for (int d : reg.VarDeg) cout << d << " ";
            cout << endl;

            cout << "  Monomial basis (" << reg.nb << "): ";
            for (auto& m : reg.MonomialBasis) cout << "[" << m << "] ";
            cout << endl;

            // Verify: nb should be > 0 for a reasonable system
            if (reg.nb > 0) {
                cout << "  PASS" << endl;
            } else {
                cout << "  FAIL: nb=0" << endl;
                allPass = false;
            }
        }
    } catch (const exception& e) {
        cerr << "FAIL: " << e.what() << endl;
        allPass = false;
    }

    // Test 2: 2-variable system (more interesting)
    cout << "\n=== RegionSolver: 2-variable system ===" << endl;
    ne = 2;
    vector<string> eqs2 = {
        "A1^2-2",      // A1^2 = 2
        "A1*A2-1",     // A2 = 1/A1 = A1/2
        "A1*B1-1",     // B1 = 1/A1
        "A2*B2-1"      // B2 = 1/A2 = A1 (since A2 = A1/2, B2 = 2/A1)
    };

    try {
        vector<int> sector2 = {1, 1};
        auto regions = RegionSolver::solveRegion(eqs2, sector2, ne, p);

        cout << "\nFound " << regions.size() << " region(s)" << endl;
        for (size_t i = 0; i < regions.size(); ++i) {
            auto& reg = regions[i];
            cout << "Region " << i << ": nb=" << reg.nb << endl;
            cout << "  VarIndep:";
            for (auto& v : reg.VarIndep) cout << " " << v;
            cout << "\n  VarDep:";
            for (auto& v : reg.VarDep) cout << " " << v;
            cout << "\n  GB:";
            for (auto& g : reg.MinPoly) cout << " [" << g << "]";
            cout << "\n  MonomialBasis (" << reg.nb << "):";
            for (auto& m : reg.MonomialBasis) cout << " [" << m << "]";
            cout << endl;
        }
        if (!regions.empty()) cout << "PASS" << endl;
    } catch (const exception& e) {
        cerr << "FAIL: " << e.what() << endl;
        allPass = false;
    }

    cout << "\n=== " << (allPass ? "ALL PASS" : "SOME FAILED") << " ===" << endl;
    return allPass ? 0 : 1;
}
