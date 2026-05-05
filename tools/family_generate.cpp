// family_generate.cpp — CLI: family.json → .bin files
//
// Phase 1: uses stub intermediate data loader (from existing .bin files).
// Future: replace stub with Singular-based Steps 1-2.
//
// Usage:
//   ./family_generate families/SR5m.json --source-dir .
//   ./family_generate families/SR5m.json --output /tmp/out --diff

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/stat.h>

#include "FamilyConfig.hpp"
#include "BinaryIBPWriter.hpp"
#include "BinaryRingWriter.hpp"

using namespace std;
using namespace firefly;

// ==========================================
// Stub loader: reads intermediate data from existing .bin files.
// In future phases, this will be replaced by Steps 1-2 (Singular).
// ==========================================
struct StubIntermediate {
    vector<RegimeSparseData<FFInt>> regimes;
    vector<RingCellRaw<FFInt>>      ringCells;
};

inline bool loadStubIntermediate(const FamilyDef& fam, StubIntermediate& out,
                                  const string& sourceDir) {
    string ibpFile  = sourceDir + "/IBPMat_" + fam.name + ".bin";
    string ringFile = sourceDir + "/RingData_" + fam.name + ".bin";

    out.regimes = readIBPMatrixRaw<FFInt>(ibpFile);

    auto cells = readRingDataRaw<FFInt>(ringFile);
    out.ringCells = std::move(cells);

    int ne = out.regimes.empty() ? 0 : out.regimes[0].ne;

    cout << "Stub loaded " << out.regimes.size() << " regime(s) from " << ibpFile
         << " and " << out.ringCells.size() << " ring cell(s) from " << ringFile
         << " (ne=" << ne << ")" << endl;
    return true;
}

// Byte-by-byte file comparison
inline bool filesIdentical(const string& a, const string& b) {
    ifstream fa(a, ios::binary | ios::ate);
    ifstream fb(b, ios::binary | ios::ate);
    if (!fa || !fb) return false;
    if (fa.tellg() != fb.tellg()) return false;
    fa.seekg(0); fb.seekg(0);
    return equal(istreambuf_iterator<char>(fa), istreambuf_iterator<char>(),
                 istreambuf_iterator<char>(fb));
}

inline bool mkdirP(const string& path) {
    string cmd = "mkdir -p \"" + path + "\"";
    return system(cmd.c_str()) == 0;
}

// ==========================================
// Main
// ==========================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << R"(Usage: family_generate <family.json> [options]

  Reads a family.json config, generates IBPMat_*.bin and RingData_*.bin.

Options:
  --source-dir DIR   Source directory for stub intermediate data (default: .)
  --output DIR       Output directory (default: .)
  --diff             Compare generated .bin files against source-dir originals

  Current phase uses stub intermediate data from existing .bin files.
  Future phases will compute via Singular CAS subprocess.
)";
        return 1;
    }

    string jsonFile  = argv[1];
    string sourceDir = ".";
    string outputDir = ".";
    bool   doDiff    = false;

    for (int i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--source-dir" && i + 1 < argc) {
            sourceDir = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            outputDir = argv[++i];
        } else if (arg == "--diff") {
            doDiff = true;
        } else {
            cerr << "Unknown flag: " << arg << endl;
            return 1;
        }
    }

    // 0. Ensure output directory exists
    if (outputDir != ".") mkdirP(outputDir);

    // 1. Parse family config
    cout << "=== FamilyGenerate ===" << endl;
    FamilyDef fam;
    try {
        fam = parseFamilyConfig(jsonFile);
    } catch (const exception& e) {
        cerr << "Failed to parse " << jsonFile << ": " << e.what() << endl;
        return 1;
    }

    cout << "Family: " << fam.name
         << " (L=" << fam.nLoop() << " E=" << fam.nExt()
         << " N=" << fam.nProp() << ")" << endl;
    cout << "Modulus: " << fam.modulus << endl;

    // Set FFInt modulus
    FFInt::set_new_prime(static_cast<uint64_t>(fam.modulus));

    // 2. Load intermediate data (stub: from existing .bin in source dir)
    StubIntermediate data;
    if (!loadStubIntermediate(fam, data, sourceDir)) {
        cerr << "Failed to load intermediate data." << endl;
        return 1;
    }

    int ne = data.regimes.empty() ? 0 : data.regimes[0].ne;

    // 3. Generate output files
    string ibpOut  = outputDir + "/IBPMat_" + fam.name + ".bin";
    string ringOut = outputDir + "/RingData_" + fam.name + ".bin";

    writeIBPMatrixBinary<FFInt>(ibpOut, data.regimes);
    writeRingDataBinary<FFInt>(ringOut, data.ringCells, ne);

    // 4. Optional diff against source originals
    if (doDiff) {
        string ibpSrc  = sourceDir + "/IBPMat_" + fam.name + ".bin";
        string ringSrc = sourceDir + "/RingData_" + fam.name + ".bin";
        bool ibpOk  = filesIdentical(ibpSrc, ibpOut);
        bool ringOk = filesIdentical(ringSrc, ringOut);

        cout << "\n--- Diff ---" << endl;
        cout << "IBP:  " << (ibpOk  ? "IDENTICAL" : "DIFFER") << endl;
        cout << "Ring: " << (ringOk ? "IDENTICAL" : "DIFFER") << endl;
        cout << "=== " << (ibpOk && ringOk ? "DIFF PASS" : "DIFF FAIL") << " ===" << endl;
        return ibpOk && ringOk ? 0 : 1;
    }

    cout << "\nDone. Output:" << endl;
    cout << "  " << ibpOut << endl;
    cout << "  " << ringOut << endl;
    return 0;
}
