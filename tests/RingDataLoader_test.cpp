#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

// 定义数据结构
using SectorVec = std::vector<int>;
using MatrixVec = std::vector<double>; // 扁平化的矩阵
using MatrixList = std::vector<MatrixVec>; // 一个 Ring 里的多个矩阵

struct SimulationData {
    std::vector<SectorVec> sectorlist;
    std::vector<MatrixList> A_list;
    std::vector<MatrixList> Ainv_list;
};

SimulationData LoadSimulationData(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // 1. 读取全局头信息
    int totalCells = 0;
    int numMatrices = 0; // ne
    
    file.read(reinterpret_cast<char*>(&totalCells), sizeof(int));
    file.read(reinterpret_cast<char*>(&numMatrices), sizeof(int));

    std::cout << "Loading " << totalCells << " cells, each with " 
              << numMatrices << " matrices." << std::endl;

    SimulationData data;
    // 预分配外层内存
    data.sectorlist.resize(totalCells);
    data.A_list.resize(totalCells);
    data.Ainv_list.resize(totalCells);

    for (int i = 0; i < totalCells; ++i) {
        // --- 读取 LimitSector ---
        int sectorSize = 0;
        file.read(reinterpret_cast<char*>(&sectorSize), sizeof(int));
        
        data.sectorlist[i].resize(sectorSize);
        file.read(reinterpret_cast<char*>(data.sectorlist[i].data()), 
                  sectorSize * sizeof(int));

        // --- 读取矩阵维度 nb ---
        int nb = 0;
        file.read(reinterpret_cast<char*>(&nb), sizeof(int));
        int matLen = nb * nb; // 扁平化后的长度

        // --- 准备 A_list 和 Ainv_list 的容器 ---
        data.A_list[i].resize(numMatrices);
        data.Ainv_list[i].resize(numMatrices);

        // --- 读取 A 矩阵数据 ---
        for (int m = 0; m < numMatrices; ++m) {
            data.A_list[i][m].resize(matLen);
            // 直接读取 double 数组
            file.read(reinterpret_cast<char*>(data.A_list[i][m].data()), 
                      matLen * sizeof(double));
        }

        // --- 读取 Ainv 矩阵数据 ---
        for (int m = 0; m < numMatrices; ++m) {
            data.Ainv_list[i][m].resize(matLen);
            file.read(reinterpret_cast<char*>(data.Ainv_list[i][m].data()), 
                      matLen * sizeof(double));
        }
    }

    file.close();
    return data;
}

int main() {
    try {
        SimulationData data = LoadSimulationData("IBPRingData_DB.bin");

        // --- 验证读取结果 ---
        if (!data.sectorlist.empty()) {
            // 打印第一个 cell 的 LimitSector
            std::cout << "Cell 0 LimitSector: ";
            for (int v : data.sectorlist[0]) std::cout << v << " ";
            std::cout << std::endl;

            // 验证 A 矩阵 (Cell 0, Matrix 0, Element 0)
            if (!data.A_list[0].empty() && !data.A_list[0][0].empty()) {
                std::cout << "A[0][0] first element: " << data.A_list[0][0][0] << std::endl;
            }
            
            // 验证 Ainv 矩阵
            if (!data.Ainv_list[0].empty() && !data.Ainv_list[0][0].empty()) {
                std::cout << "Ainv[0][0] first element: " << data.Ainv_list[0][0][0] << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}