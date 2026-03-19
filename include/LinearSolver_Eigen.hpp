#ifndef LINEAR_SOLVER_EIGEN_HPP
#define LINEAR_SOLVER_EIGEN_HPP

#include <Eigen/Dense>

using namespace Eigen;

template <typename T>
void truncate(T &x)
{
    if(abs(x) < 1e-20) { x=static_cast<T>(0.0); }
}

/**
 * 线性方程组解结构体
 * Mext: 第一列是特解 v，后续列是基础解系 M
 * S: 独立变量（自由元）在原始 A 中的列索引
 */
template <typename T>
struct LinearSystemResult {
    bool hasSolution;
    vector<vector<T>> Mext; // 拼接矩阵 (v | M)
    vector<int> S;          // 独立变量索引集合
};

/**
 * 线性求解器
 */
template <typename T>
LinearSystemResult<T> solveLinearSystem_Eigen(const vector<vector<T>>& A_in, const vector<T>& b_in = {}) {
    int rows = A_in.size();
    if (rows == 0) return {false, {}, {}};
    int cols = A_in[0].size();

    Matrix<T, Dynamic, Dynamic> A(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) A(i, j) = A_in[i][j];
    
    Matrix<T, Dynamic, 1> b;
    if (b_in.empty()) {
        b = Matrix<T, Dynamic, 1>::Zero(rows);
    } else {
        // 检查维度匹配
        //if (b_in.size() != rows) { cerr << "Error: Rows not match"; return; }
        b = Map<const Matrix<T, Dynamic, 1>>(b_in.data(), b_in.size());
    }

    // 1. 使用全主元 LU
    FullPivLU<Matrix<T, Dynamic, Dynamic>> lu(A);

    // 2. 关键：手动设置一个稍大的阈值来“杀掉”数值噪声
    // 对于 IBP 这种有理数系数，1e-10 通常比 1e-12 更稳健
    T precisionLimit = 1e-10; 
    lu.setThreshold(precisionLimit);

    int r = lu.rank();
    int numFreeVars = cols - r;

    // 3. 求解特解
    Matrix<T, Dynamic, 1> v = lu.solve(b);

    // 4. 获取通解并进行“强制截断”
    Matrix<T, Dynamic, Dynamic> kernel = lu.kernel();
    
    // --- 强制对齐逻辑 ---
    // 如果 kernel 拿多了，只取前 numFreeVars 列
    // 如果 kernel 拿少了，说明 rank 判定的阈值和 kernel 内部不一致，以 rank 为准
    Matrix<T, Dynamic, Dynamic> finalKernel;
    if (kernel.cols() > numFreeVars) {
        finalKernel = kernel.leftCols(numFreeVars);
    } else if (kernel.cols() < numFreeVars) {
        // 这种情况极少发生，如果发生，通常意味着矩阵几乎全零
        finalKernel = kernel; 
        numFreeVars = kernel.cols();
        r = cols - numFreeVars;
    } else {
        finalKernel = kernel;
    }

    LinearSystemResult<T> res;
    // 验证解的有效性
    if (!(A * v).isApprox(b, 1e-8)) {
        res.hasSolution = false;
        return res;
    }
    res.hasSolution = true;

    // 5. 组装结果
    res.Mext.assign(cols, vector<T>(1 + numFreeVars));
    for (int i = 0; i < cols; ++i) {
        res.Mext[i][0] = v(i);
        for (int j = 0; j < numFreeVars; ++j) {
            res.Mext[i][j + 1] = finalKernel(i, j);
            truncate(res.Mext[i][j + 1]);
        }
    }

    // 6. 独立变量索引
    auto q_indices = lu.permutationQ().indices();
    for (int i = r; i < cols; ++i) {
        res.S.push_back(q_indices[i]);
    }

    return res;
}

template<typename T>
void printResult(const LinearSystemResult<T> res) {
    cout << "solved" <<endl;
    for(int i = 0; i < res.Mext.size(); ++i) {
        for(int  j = 0 ; j < res.Mext[i].size(); ++j) {
            cout << res.Mext[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "indep var: " << endl;
    for(int i = 0; i < res.S.size(); ++i) { 
        cout << res.S[i] << "  ";
    }
    cout << endl;
}

#endif // LINEAR_SOLVER_EIGEN_HPP