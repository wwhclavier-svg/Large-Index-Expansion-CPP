// diagnose_sr5m.cpp — Check exists flags per regime
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>

using namespace std;

int32_t readBE32(ifstream& in) {
    uint8_t b[4];
    in.read(reinterpret_cast<char*>(b), 4);
    return (int32_t)((uint32_t)b[0] << 24 | (uint32_t)b[1] << 16 |
                     (uint32_t)b[2] << 8  | (uint32_t)b[3]);
}

int64_t readBE64(ifstream& in) {
    uint8_t b[8];
    in.read(reinterpret_cast<char*>(b), 8);
    return (int64_t)((uint64_t)b[0] << 56 | (uint64_t)b[1] << 48 |
                     (uint64_t)b[2] << 40 | (uint64_t)b[3] << 32 |
                     (uint64_t)b[4] << 24 | (uint64_t)b[5] << 16 |
                     (uint64_t)b[6] << 8  | (uint64_t)b[7]);
}

int main() {
    string filename = "IBPMat_SR5m.bin";
    ifstream in(filename, ios::binary);
    if (!in) { cerr << "Cannot open " << filename << endl; return 1; }

    char magic[4];
    in.read(magic, 4);
    int32_t numRegs = readBE32(in);
    cout << "Regimes: " << numRegs << endl;

    vector<string> opOrder = {"M1","N1","K1","F0","F2","K1s","K2s","F2s"};

    for (int r = 0; r < numRegs && r < 5; ++r) {
        int32_t nibp = readBE32(in);
        int32_t ne   = readBE32(in);
        int32_t nb   = readBE32(in);
        int32_t incre = readBE32(in);
        int64_t mod  = readBE64(in);

        cout << "\nR" << r << " (nb=" << nb << " ne=" << ne << " nibp=" << nibp << "):" << endl;

        for (const auto& op : opOrder) {
            uint8_t exists;
            in.read(reinterpret_cast<char*>(&exists), 1);
            cout << "  " << op << ": exists=" << (int)exists;

            if (exists) {
                int32_t dims_len = readBE32(in);
                cout << " dims_len=" << dims_len << " [";
                vector<int32_t> dims(dims_len);
                for (int i = 0; i < dims_len; ++i) {
                    dims[i] = readBE32(in);
                    cout << dims[i] << (i+1 < dims_len ? "," : "");
                }
                cout << "]";

                int32_t rp_len = readBE32(in);
                cout << " rowPtr_len=" << rp_len;
                for (int i = 0; i < rp_len; ++i) readBE32(in); // skip

                int32_t ci_len = readBE32(in);
                cout << " colIdx_len=" << ci_len;
                for (int i = 0; i < ci_len; ++i) readBE32(in); // skip

                int32_t vals_len = readBE32(in);
                cout << " vals_len=" << vals_len;
                if (vals_len == 0) cout << " *** EMPTY TENSOR ***";
                for (int i = 0; i < vals_len; ++i) readBE64(in); // skip
            }
            cout << endl;
        }
    }
    return 0;
}
