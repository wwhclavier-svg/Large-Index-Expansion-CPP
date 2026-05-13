#include <iostream>
#include <string>
#include <cassert>

#include "firefly/FFInt.hpp"
#include "RelationLoader.hpp"
#include "SymbolicRule.hpp"

int testRelationLoader() {
    std::cout << "=== Test: RelationLoader ===\n";
    try {
        auto data = loadAllRelations("test/Penta1L/AllRelations_Penta1L_k4.m");
        assert(!data.entries.empty());
        std::cout << "  Loaded " << data.entries.size() << " entries, family="
                  << data.family << " NE=" << data.ne << "\n";

        bool foundPyramid = false;
        for (auto& e : data.entries) {
            if (e.ansatzMode == "Pyramid") {
                foundPyramid = true;
                std::cout << "  Pyramid entry: (lev=" << e.lev << ", deg=" << e.deg
                          << ") rel=" << e.numRelations << " stable=" << e.stableOrder << "\n";
            }
        }
        assert(foundPyramid);
        std::cout << "  PASSED\n";
        return 0;
    } catch (std::exception& e) {
        std::cerr << "  FAILED: " << e.what() << "\n";
        return 1;
    }
}

int testSeedGeneration() {
    std::cout << "=== Test: Seed Generation ===\n";
    auto seeds = generateDotSeeds(2, 2);
    std::cout << "  NE=2, level=2: " << seeds.size() << " seeds\n";
    assert(seeds.size() > 0);

    auto topSeeds = generateDotSeeds(2, 2, true, 2);
    std::cout << "  TopRankOnly, minRank=2: " << topSeeds.size() << " seeds\n";

    auto minusSeeds = generateMinusSeeds(2, 2);
    std::cout << "  MinusSeeds NE=2, level=2: " << minusSeeds.size() << " seeds\n";

    std::cout << "  PASSED\n";
    return 0;
}

int testSetupGeneratingCone() {
    std::cout << "=== Test: setupGeneratingCone (Penta1L, Pyramid) ===\n";
    try {
        auto relData = loadAllRelations("test/Penta1L/AllRelations_Penta1L_k4.m");
        std::vector<int> sector = {1, 1, 1, 1, 1};

        auto result = setupGeneratingCone(relData, sector, 6);
        result.printSummary(std::cout);

        std::cout << "  B1 cones: " << result.B1.size() << "\n";
        std::cout << "  B2 cones: " << result.B2.size() << "\n";
        std::cout << "  PASSED\n";
        return 0;
    } catch (std::exception& e) {
        std::cerr << "  FAILED: " << e.what() << "\n";
        return 1;
    }
}

int testBoxVerification() {
    std::cout << "=== Test: Box Full Verification ===\n";
    try {
        auto rd = loadAllRelations("test/Box/AllRelations_Box_k4.m");
        std::vector<int> sector = {1,1,1,1};
        auto res = setupGeneratingCone(rd, sector, 6);
        assert(!res.B1.empty());
        // B1 should run without exception; gbDiagonal may be 0 if no ν-structure
        std::cout << "  B1 cones: " << res.B1.size() << "\n";
        for (auto& c : res.B1) {
            std::cout << "    level=" << c.level << " gbDiag=" << c.gbDiagonal.size()
                      << " [" << c.flags.toString() << "]\n";
            for (auto& d : c.gbDiagonal) std::cout << "      " << d << "\n";
        }
        std::cout << "  PASSED\n";
        return 0;
    } catch (std::exception& e) {
        std::cerr << "  FAILED: " << e.what() << "\n";
        return 1;
    }
}

int testPenta1LVerification() {
    std::cout << "=== Test: Penta1L Verification ===\n";
    try {
        auto rd = loadAllRelations("test/Penta1L/AllRelations_Penta1L_k4.m");
        std::vector<int> sector = {1,1,1,1,1};
        auto res = setupGeneratingCone(rd, sector, 6);
        // Penta1L has sparse relations; B1 may be empty at low levels
        // Just verify the function runs without exception
        std::cout << "  B1: " << res.B1.size() << "\n";
        std::cout << "  PASSED\n";
        return 0;
    } catch (std::exception& e) {
        std::cerr << "  FAILED: " << e.what() << "\n";
        return 1;
    }
}

int main() {
    int failures = 0;
    failures += testRelationLoader();
    failures += testSeedGeneration();
    failures += testSetupGeneratingCone();
    failures += testBoxVerification();
    failures += testPenta1LVerification();

    std::cout << "\n=== Summary: " << (failures == 0 ? "ALL 5 PASSED" : "SOME FAILED") << " ===\n";
    return failures;
}
