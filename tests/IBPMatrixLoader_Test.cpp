#include "IBPMatrixLoader.hpp"

template<typename T>
void printMatrix(vector<vector<vector<T>>> Mat){
    for(int i = 0; i < Mat.size(); ++i) {
        for(int j = 0; j < Mat[0].size(); ++j) {
            cout << Mat[i][j][0];
            if(j < Mat[0].size()-1) cout << "\t";
        }
        cout << endl;
    }
}

int main() {
    // vector<IBPMatrixE<double>> matlist;
    try {
        // 假设 json 文件名为 input_sparse.json
        auto ibpmatlist = loadAllIBPMatrices<double>("IBPMatAll_SR.json");
        
        // 验证 M1 (3D)
        auto mat = ibpmatlist[0];
        cout << "nibp = " << mat.nibp << "  ne = " <<mat.ne << "  nb = " << mat.nb << endl;
        if (mat.M1.size() > 0) 
            cout << "M1 Dims: " << mat.M1.size() << "x" << mat.M1[0].size() << "x" << mat.M1[0][0].size() << endl;
        
        printMatrix(mat.M1);
            
        // 验证 F2 (4D)
        if (mat.F2.size() > 0)
            cout << "F2 Dims: " << mat.F2.size() << "x" << mat.F2[0].size() << "x" << mat.F2[0][0].size() << "x" << mat.F2[0][0][0].size() << endl;

        for(int i = 0; i < mat.F2.size(); ++i) {
            cout << "F2["<<i<<"] = " <<endl;
            printMatrix(mat.F2[i]);
        }

    } catch (exception& e) {
        cerr << e.what() << endl;
    }
    return 0;
}