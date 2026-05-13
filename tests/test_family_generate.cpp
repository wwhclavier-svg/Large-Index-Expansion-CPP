// test_family_generate.cpp — Phase 0: Round-trip verification of binary writers
#include <iostream>
#include <string>
#include "firefly/FFInt.hpp"
#include "BinaryIBPWriter.hpp"
#include "BinaryRingWriter.hpp"

using namespace std;
using namespace firefly;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <family_name>" << endl;
        cerr << "  Verifies round-trip: .bin → read raw → write → compare" << endl;
        return 1;
    }

    string family = argv[1];
    string ibpFile  = "data/IBPMat_" + family + ".bin";
    string ringFile = "data/RingData_" + family + ".bin";

    cout << "=== Binary Writer Round-Trip Verification ===" << endl;
    bool allPass = true;

    // Test 1: IBP Matrix round-trip
    cout << "\n--- IBP Matrix: " << ibpFile << " ---" << endl;
    try {
        if (!verifyRoundTrip<FFInt>(ibpFile))
            allPass = false;
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        allPass = false;
    }

    // Test 2: Ring Data round-trip
    cout << "\n--- Ring Data: " << ringFile << " ---" << endl;
    try {
        if (!verifyRingRoundTrip<FFInt>(ringFile))
            allPass = false;
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        allPass = false;
    }

    cout << "\n=== Result: " << (allPass ? "ALL PASS" : "SOME FAILED") << " ===" << endl;
    return allPass ? 0 : 1;
}
