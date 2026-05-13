#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <set>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"
#include "SeriesCoefficientIO.hpp"
#include "firefly/FFInt.hpp"
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

// ============================================================
// Section 1: Utility — Print with indentation
// ============================================================
void printHeader(const string& title) {
    cout << "\n" << string(60, '=') << endl;
    cout << "  " << title << endl;
    cout << string(60, '=') << endl;
}

// ============================================================
// Section 2: Reference Data Loading from Mathematica JSON
// ============================================================
struct RefCoefficient {
    int order;
    vector<int> seed;
    int j;       // 1-indexed (Mathematica convention)
    string value; // string representation for exact comparison
};

struct RefRegion {
    int region_index;
    int ne;
    int nb;
    int nibp;
    vector<int> sector;
    // coefficients grouped by order
    vector<vector<RefCoefficient>> coeffs_by_order;
};

struct RefData {
    string family;
    int order;
    int64_t modulus;
    vector<RefRegion> regions;

    bool loaded = false;
};

RefData loadMathematicaReference(const string& filename) {
    RefData ref;
    ifstream in(filename);
    if (!in) {
        cout << "[INFO] No Mathematica reference file found: " << filename << endl;
        cout << "[INFO] Will perform self-consistency checks instead." << endl;
        return ref;
    }

    json j;
    try {
        in >> j;
    } catch (const json::parse_error& e) {
        cerr << "[WARN] Failed to parse reference JSON: " << e.what() << endl;
        return ref;
    }

    ref.family = j.value("family", "unknown");
    ref.order = j.value("order", 0);
    ref.modulus = j.value("modulus", 0);

    for (const auto& rj : j["regions"]) {
        RefRegion reg;
        reg.region_index = rj.value("region_index", 0);
        reg.ne = rj.value("ne", 0);
        reg.nb = 1; // will be inferred from data
        reg.nibp = 0;

        if (rj.contains("sector")) {
            for (const auto& s : rj["sector"])
                reg.sector.push_back(s.get<int>());
        }

        if (rj.contains("coefficients_by_order")) {
            for (const auto& ord_j : rj["coefficients_by_order"]) {
                int ord = ord_j.value("order", 0);
                vector<RefCoefficient> coeffs;

                if (ord_j.contains("coefficients")) {
                    for (const auto& cj : ord_j["coefficients"]) {
                        RefCoefficient rc;
                        rc.order = ord;
                        rc.j = cj.value("j", 1);
                        rc.value = cj.value("value", "0");

                        if (cj.contains("seed")) {
                            for (const auto& s : cj["seed"])
                                rc.seed.push_back(s.get<int>());
                        }

                        reg.nb = max(reg.nb, rc.j);
                        coeffs.push_back(rc);
                    }
                }
                reg.coeffs_by_order.push_back(coeffs);
            }
        }

        ref.regions.push_back(reg);
    }

    ref.loaded = true;
    cout << "[INFO] Loaded Mathematica reference: " << filename << endl;
    cout << "       family=" << ref.family
         << "  order=" << ref.order
         << "  modulus=" << ref.modulus
         << "  regions=" << ref.regions.size() << endl;

    return ref;
}

// ============================================================
// Section 3: Coefficient Extraction from C++ seriesCoefficient
// ============================================================

template<typename T>
struct FlatCoefficient {
    int order;
    vector<int> seed;
    int level;
    int cid;
    int j;       // 0-indexed (C++ convention)
    int i;       // solution component
    T value;
};

template<typename T>
vector<FlatCoefficient<T>> extractCoefficients(
    const seriesCoefficient<T>& C,
    int max_order,
    int incre,
    int ne,
    int nb,
    long long (&BINOM)[MAX_VAL][MAX_VAL])
{
    vector<FlatCoefficient<T>> result;

    for (int k = 0; k <= max_order; ++k) {
        int lmax = incre * k;
        for (int l = 0; l <= lmax; ++l) {
            long long num_seeds = BINOM[l + ne - 1][ne - 1];
            for (int cid = 0; cid < num_seeds; ++cid) {
                vector<int> seed = readIndex(cid, l, ne);
                for (int j = 0; j < nb; ++j) {
                    // Extract particular solution (i=0) and all homogeneous
                    // components (i=1..nimax). For Mathematica comparison
                    // we mainly care about i=0 for the primary branch.
                    for (int i = 0; i <= C.getNimax(); ++i) {
                        T val = C(k, l, cid, j, i);
                        if (val != T(0)) {
                            result.push_back({k, seed, l, cid, j, i, val});
                        }
                    }
                }
            }
        }
    }

    return result;
}

// ============================================================
// Section 4: Comparison Logic
// ============================================================

// Parse string value from Mathematica JSON into a numeric type
template<typename T>
T parseMathematicaValue(const string& s);

template<>
double parseMathematicaValue<double>(const string& s) {
    // Mathematica InputForm: "1", "-3/2", "0.5", etc.
    // Handle rational numbers
    size_t slash = s.find('/');
    if (slash != string::npos) {
        double num = stod(s.substr(0, slash));
        double den = stod(s.substr(slash + 1));
        return num / den;
    }
    return stod(s);
}

template<>
firefly::FFInt parseMathematicaValue<firefly::FFInt>(const string& s) {
    // For finite field, values are integers modulo prime
    // Handle possible negative integers (MMA may output negative ints
    // that need to be mapped to [0, prime-1])
    size_t slash = s.find('/');
    if (slash != string::npos) {
        // Rational in finite field: numerator * inverse(denominator) mod p
        int64_t num = stoll(s.substr(0, slash));
        int64_t den = stoll(s.substr(slash + 1));
        firefly::FFInt fnum(num);
        firefly::FFInt fden(den);
        return fnum / fden;  // FireFly handles modular inverse
    }
    return firefly::FFInt(stoll(s));
}

template<typename T>
struct ComparisonResult {
    int total_compared = 0;
    int total_match = 0;
    int total_mismatch = 0;
    int total_only_in_cpp = 0;
    int total_only_in_mma = 0;

    vector<string> mismatches;  // First few mismatch details

    void addMismatch(const string& detail) {
        if (mismatches.size() < 20) {
            mismatches.push_back(detail);
        }
    }

    void report() const {
        cout << "\n  Comparison summary:" << endl;
        cout << "    Total compared:    " << total_compared << endl;
        cout << "    Matching:          " << total_match << endl;
        cout << "    Mismatching:       " << total_mismatch << endl;
        cout << "    Only in C++:       " << total_only_in_cpp << endl;
        cout << "    Only in MMA ref:   " << total_only_in_mma << endl;

        if (!mismatches.empty()) {
            cout << "\n  First " << mismatches.size() << " mismatches:" << endl;
            for (const auto& m : mismatches) {
                cout << "    " << m << endl;
            }
        }
    }

    bool passed() const {
        return total_mismatch == 0 && total_only_in_cpp == 0 && total_only_in_mma == 0;
    }
};

template<typename T>
ComparisonResult<T> compareWithReference(
    const vector<seriesCoefficient<T>>& cppResults,
    const RefRegion& refRegion,
    int order,
    int incre,
    int ne,
    int nb,
    long long (&BINOM)[MAX_VAL][MAX_VAL])
{
    ComparisonResult<T> result;

    // Build lookup from Mathematica reference: (order, seed_key, j) -> value
    map<tuple<int, string, int>, string> mmaLookup;
    for (const auto& ord_coeffs : refRegion.coeffs_by_order) {
        for (const auto& rc : ord_coeffs) {
            // Build seed key
            string seed_key;
            for (size_t i = 0; i < rc.seed.size(); ++i) {
                if (i > 0) seed_key += ",";
                seed_key += to_string(rc.seed[i]);
            }
            mmaLookup[{rc.order, seed_key, rc.j}] = rc.value;
        }
    }

    // Iterate over C++ coefficients and compare
    for (size_t sol_idx = 0; sol_idx < cppResults.size(); ++sol_idx) {
        const auto& C = cppResults[sol_idx];
        int nimax = C.getNimax();

        for (int k = 0; k <= order; ++k) {
            int lmax = incre * k;
            for (int l = 0; l <= lmax; ++l) {
                long long num_seeds = BINOM[l + ne - 1][ne - 1];
                for (int cid = 0; cid < num_seeds; ++cid) {
                    vector<int> seed = readIndex(cid, l, ne);

                    // Build seed key
                    string seed_key;
                    for (size_t i = 0; i < seed.size(); ++i) {
                        if (i > 0) seed_key += ",";
                        seed_key += to_string(seed[i]);
                    }

                    for (int j = 0; j < nb; ++j) {
                        // Get C++ value (particular solution i=0)
                        T cpp_val = C(k, l, cid, j, 0);

                        auto mma_it = mmaLookup.find({k, seed_key, j + 1}); // j+1: 1-indexed

                        if (mma_it != mmaLookup.end()) {
                            result.total_compared++;
                            T mma_val = parseMathematicaValue<T>(mma_it->second);

                            bool match;
                            if constexpr (is_same_v<T, double>) {
                                match = (abs(cpp_val - mma_val) < 1e-10 ||
                                         abs(cpp_val - mma_val) / max(1.0, abs(mma_val)) < 1e-10);
                            } else {
                                match = (cpp_val == mma_val);
                            }

                            if (match) {
                                result.total_match++;
                            } else {
                                result.total_mismatch++;
                                ostringstream oss;
                                oss << "order=" << k << " seed=[" << seed_key
                                    << "] j=" << (j+1)
                                    << " C++=" << cpp_val << " MMA=" << mma_val;
                                result.addMismatch(oss.str());
                            }
                        } else {
                            // C++ has non-zero but MMA doesn't mention it
                            if (cpp_val != T(0)) {
                                result.total_only_in_cpp++;
                            }
                        }
                    }
                }
            }
        }
    }

    // Check for coefficients only in MMA but not in C++
    for (const auto& [key, mma_str] : mmaLookup) {
        auto [ord, seed_key, j] = key;

        // Parse seed
        vector<int> seed;
        size_t pos = 0, last = 0;
        while ((pos = seed_key.find(',', last)) != string::npos) {
            seed.push_back(stoi(seed_key.substr(last, pos - last)));
            last = pos + 1;
        }
        seed.push_back(stoi(seed_key.substr(last)));

        int l = 0;
        for (int s : seed) l += s;
        int cid = getIndex(seed, l);

        // Check if it exists in any C++ branch
        bool found = false;
        for (const auto& C : cppResults) {
            if (ord <= C.getKmax() && l <= incre * ord &&
                cid < BINOM[l + ne - 1][ne - 1] && (j-1) < nb) {
                T cpp_val = C(ord, l, cid, j-1, 0);
                T mma_val = parseMathematicaValue<T>(mma_str);

                bool match;
                if constexpr (is_same_v<T, double>) {
                    match = (abs(cpp_val - mma_val) < 1e-10);
                } else {
                    match = (cpp_val == mma_val);
                }

                if (match) {
                    found = true;
                    break;
                }
            }
        }

        if (!found) {
            // Check if MMA value is actually zero
            T mma_val = parseMathematicaValue<T>(mma_str);
            if (mma_val != T(0)) {
                result.total_only_in_mma++;
            }
        }
    }

    return result;
}

// ============================================================
// Section 5: Self-Consistency Checks (when no MMA reference)
// ============================================================

template<typename T>
bool checkInitialCondition(const vector<seriesCoefficient<T>>& cppResults,
                           int ne, int nb) {
    printHeader("Self-Check: Initial Condition");
    bool ok = true;

    for (size_t s = 0; s < cppResults.size(); ++s) {
        T val0 = cppResults[s](0, 0, 0, 0, 0);  // C(0,0,0,0,0)

        if (s == 0) {
            // First branch: should be 1
            if (val0 == T(1)) {
                cout << "  [PASS] Branch " << s << ": C(0,0,0,0,0) = " << val0 << endl;
            } else {
                cout << "  [FAIL] Branch " << s << ": C(0,0,0,0,0) = " << val0
                     << " (expected 1)" << endl;
                ok = false;
            }
        } else {
            // Other branches: should have C(0,0,0,0,0) = 0
            // (they come from higher-order branching)
            // Actually, they should have their own starting values
            cout << "  [INFO] Branch " << s << ": C(0,0,0,0,0) = " << val0 << endl;
        }
    }

    return ok;
}

template<typename T>
bool checkNonDecreasingBranches(const vector<seriesCoefficient<T>>& cppResults) {
    printHeader("Self-Check: Branch Count");

    cout << "  Number of solution branches: " << cppResults.size() << endl;
    if (cppResults.size() >= 1) {
        cout << "  [PASS] At least one solution branch produced." << endl;
        return true;
    } else {
        cout << "  [FAIL] No solution branches." << endl;
        return false;
    }
}

template<typename T>
bool checkCoefficientSanity(const vector<seriesCoefficient<T>>& cppResults,
                            int order, int incre, int ne, int nb,
                            long long (&BINOM)[MAX_VAL][MAX_VAL]) {
    printHeader("Self-Check: Coefficient Sanity");
    bool ok = true;
    int nonZeroCount = 0;

    for (size_t s = 0; s < cppResults.size(); ++s) {
        const auto& C = cppResults[s];
        int nimax = C.getNimax();

        for (int k = 0; k <= order; ++k) {
            int lmax = incre * k;
            for (int l = 0; l <= lmax; ++l) {
                long long num_seeds = BINOM[l + ne - 1][ne - 1];
                for (int cid = 0; cid < num_seeds; ++cid) {
                    for (int j = 0; j < nb; ++j) {
                        for (int i = 0; i <= nimax; ++i) {
                            T val = C(k, l, cid, j, i);
                            if (val != T(0)) nonZeroCount++;

                            // Check no NaN for double
                            if constexpr (is_same_v<T, double>) {
                                if (isnan(val)) {
                                    cout << "  [FAIL] NaN at (s=" << s << ",k=" << k
                                         << ",l=" << l << ",cid=" << cid
                                         << ",j=" << j << ",i=" << i << ")" << endl;
                                    ok = false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cout << "  Non-zero coefficients: " << nonZeroCount << endl;
    if (ok) {
        cout << "  [PASS] No NaN values detected." << endl;
    }
    return ok;
}


// ============================================================
// Section 6: Main Test Driver
// ============================================================

template<typename T>
int runComparisonTest(const string& famname, int order, int incre,
                      const string& dataDir, const string& refFile) {

    cout << "\n========== Comparison Test: C++ vs Mathematica ==========" << endl;
    cout << "Family:     " << famname << endl;
    cout << "Order:      " << order << endl;
    cout << "Increment:  " << incre << endl;
    cout << "Data dir:   " << dataDir << endl;
    cout << "Type:       " << (is_same_v<T, double> ? "double" : "FFInt") << endl;

    // ---- Step 1: Load reference data ----
    printHeader("Step 1: Load Mathematica Reference");
    RefData ref = loadMathematicaReference(refFile);

    // ---- Step 2: Load IBP binary matrices ----
    printHeader("Step 2: Load IBP Binary Matrices");
    string ibpFile = dataDir + "/IBPMat_" + famname + ".bin";
    vector<IBPMatrixE<T>> ibpMatrices;

    try {
        ibpMatrices = loadAllIBPMatricesBinary<T>(ibpFile);
        cout << "  Loaded " << ibpMatrices.size() << " IBP matrices." << endl;

        if (!ibpMatrices.empty()) {
            const auto& m = ibpMatrices[0];
            cout << "  First matrix: ne=" << m.ne
                 << "  nb=" << m.nb
                 << "  nibp=" << m.nibp << endl;
        }
    } catch (const exception& e) {
        cerr << "  [ERROR] Failed to load IBP matrices: " << e.what() << endl;
        return 1;
    }

    // ---- Step 3: Run C++ layer recursion ----
    printHeader("Step 3: Run C++ Layer Recursion");

    auto start = chrono::high_resolution_clock::now();
    auto allResults = batchProcessRecursion<T>(ibpMatrices, order, incre);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    cout << "  Computation time: " << fixed << setprecision(3)
         << elapsed.count() << " s" << endl;

    // Print branch counts per matrix
    cout << "  Solution branches per matrix:" << endl;
    for (size_t i = 0; i < allResults.size(); ++i) {
        cout << "    Matrix " << i << ": " << allResults[i].size() << " branches";
        if (!allResults[i].empty()) {
            cout << "  (nb=" << allResults[i][0].getNb()
                 << ", nimax=" << allResults[i][0].getNimax() << ")";
        }
        cout << endl;
    }

    // ---- Step 4: Compare with reference (if available) ----
    bool allPassed = true;

    if (ref.loaded) {
        printHeader("Step 4: Compare with Mathematica Reference");

        for (size_t i = 0; i < min(allResults.size(), ref.regions.size()); ++i) {
            cout << "\n  --- Region " << i << " ---" << endl;

            const auto& cppBranches = allResults[i];
            const auto& mmaRegion = ref.regions[i];

            if (cppBranches.empty()) {
                cout << "  [SKIP] No C++ results for this region." << endl;
                continue;
            }

            int ne = cppBranches[0].getNe();
            int nb = cppBranches[0].getNb();

            auto comp = compareWithReference<T>(
                cppBranches, mmaRegion, order, incre, ne, nb, BINOM);

            comp.report();

            if (!comp.passed()) {
                allPassed = false;
            }
        }

        // Check count mismatch
        if (allResults.size() != ref.regions.size()) {
            cout << "\n  [WARN] Region count mismatch: C++=" << allResults.size()
                 << "  MMA=" << ref.regions.size() << endl;
        }

    } else {
        printHeader("Step 4: Self-Consistency Checks (no MMA reference)");

        bool check1 = checkNonDecreasingBranches<T>(
            allResults.empty() ? vector<seriesCoefficient<T>>{} : allResults[0]);

        bool check2 = true;
        if (!allResults.empty() && !allResults[0].empty()) {
            int ne = allResults[0][0].getNe();
            int nb = allResults[0][0].getNb();
            check2 = checkInitialCondition<T>(allResults[0], ne, nb);
        }

        bool check3 = true;
        if (!allResults.empty()) {
            for (size_t i = 0; i < allResults.size() && check3; ++i) {
                if (!allResults[i].empty()) {
                    int ne = allResults[i][0].getNe();
                    int nb = allResults[i][0].getNb();
                    check3 = checkCoefficientSanity<T>(
                        allResults[i], order, incre, ne, nb, BINOM);
                }
            }
        }

        allPassed = check1 && check2 && check3;
    }

    // ---- Step 5: Print final verdict ----
    printHeader("Final Verdict");
    if (allPassed) {
        cout << "  ALL CHECKS PASSED" << endl;
        return 0;
    } else {
        cout << "  SOME CHECKS FAILED" << endl;
        return 1;
    }
}


// ============================================================
// Section 7: Main Entry Point
// ============================================================

int main(int argc, char* argv[]) {
    // --- Parse arguments ---
    string famname = "bub00";
    int order = 4;
    int incre = 2;
    string dataDir = ".";
    string refFile = "";  // If empty, auto-construct from famname

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--famname" && i + 1 < argc) {
            famname = argv[++i];
        } else if (arg == "--order" && i + 1 < argc) {
            order = stoi(argv[++i]);
        } else if (arg == "--incre" && i + 1 < argc) {
            incre = stoi(argv[++i]);
        } else if (arg == "--datadir" && i + 1 < argc) {
            dataDir = argv[++i];
        } else if (arg == "--ref" && i + 1 < argc) {
            refFile = argv[++i];
        } else if (arg == "--double") {
            // Use double precision (handled below)
        } else if (arg == "--help" || arg == "-h") {
            cout << "Usage: test_compare_MMA [options]" << endl;
            cout << "  --famname NAME   Family name (default: bub00)" << endl;
            cout << "  --order N        Expansion order (default: 4)" << endl;
            cout << "  --incre N        Level increment (default: 2)" << endl;
            cout << "  --datadir DIR    Data directory (default: .)" << endl;
            cout << "  --ref FILE       Mathematica reference JSON file" << endl;
            cout << "  --double         Use double precision instead of FFInt" << endl;
            return 0;
        }
    }

    // Auto-construct reference file path
    if (refFile.empty()) {
        refFile = dataDir + "/RefCoeff_" + famname + ".json";
    }

    // --- Initialize ---
    initBinomial();

    // --- Run test ---
    // Default: use FFInt (finite field). Use --double for double precision.
    bool useDouble = false;
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "--double") {
            useDouble = true;
            break;
        }
    }

    if (useDouble) {
        return runComparisonTest<double>(famname, order, incre, dataDir, refFile);
    } else {
        // Read modulus from binary file to set FFInt prime
        // The loader sets it automatically
        return runComparisonTest<firefly::FFInt>(famname, order, incre, dataDir, refFile);
    }
}
