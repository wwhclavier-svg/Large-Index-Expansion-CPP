// diagnose_ibp_reference.cpp — read reference .bin and print all values
#include <iostream>
#include "firefly/FFInt.hpp"
#include "IBPMatrixLoader_Binary.hpp"

using namespace std;
using firefly::FFInt;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <IBPMat_*.bin>" << endl;
        return 1;
    }
    FFInt::set_new_prime(179424673);
    auto mats = loadAllIBPMatricesBinary<FFInt>(argv[1], true);

    for (size_t r = 0; r < mats.size(); ++r) {
        auto& m = mats[r];
        cout << "\n=== Regime " << r << " ===" << endl;
        cout << "nibp=" << m.nibp << " ne=" << m.ne << " nb=" << m.nb << endl;

        auto print4D = [&](const char* name,
            const vector<vector<vector<FFInt>>>& data) {
            cout << "\n" << name << ":" << endl;
            for (int i = 0; i < m.nibp; ++i) {
                for (int j = 0; j < m.ne; ++j) {
                    cout << "  [" << i << "][" << j << "]: ";
                    for (int a = 0; a < m.nb; ++a)
                        for (int b = 0; b < m.nb; ++b)
                            cout << data[i][j][a*m.nb + b].n << " ";
                    cout << endl;
                }
            }
        };

        auto print3D = [&](const char* name,
            const vector<vector<FFInt>>& data) {
            cout << "\n" << name << ":" << endl;
            for (int i = 0; i < m.nibp; ++i) {
                cout << "  [" << i << "]: ";
                for (int a = 0; a < m.nb; ++a)
                    for (int b = 0; b < m.nb; ++b)
                        cout << data[i][a*m.nb + b].n << " ";
                cout << endl;
            }
        };

        auto print5D = [&](const char* name,
            const vector<vector<vector<vector<FFInt>>>>& data) {
            cout << "\n" << name << ":" << endl;
            for (int i = 0; i < m.nibp; ++i) {
                for (int j = 0; j < m.ne; ++j) {
                    for (int k = 0; k < m.ne; ++k) {
                        cout << "  [" << i << "][" << j << "][" << k << "]: ";
                        for (int a = 0; a < m.nb; ++a)
                            for (int b = 0; b < m.nb; ++b)
                                cout << data[i][j][k][a*m.nb + b].n << " ";
                        cout << endl;
                    }
                }
            }
        };

        print4D("M1", m.M1);
        print4D("N1", m.N1);
        print4D("K1", m.K1);
        print4D("K1s", m.K1s);
        print4D("K2s", m.K2s);
        print3D("F0", m.F0);
        print5D("F2", m.F2);
        print5D("F2s", m.F2s);
    }
    return 0;
}
