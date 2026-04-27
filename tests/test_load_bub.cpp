#include <iostream>
#include <string>
#include "IBPMatrixLoader_Binary.hpp"
#include "firefly/FFInt.hpp"

int main(int argc, char* argv[]) {
    try {
        std::string family = "bub00";
        if (argc > 1) {
            family = argv[1];
        }
        
        std::string filename = "IBPMat_" + family + ".bin";
        std::cout << "=== Test Loading " << filename << " ===" << std::endl;
        
        auto ibpmatlist = loadAllIBPMatricesBinary<firefly::FFInt>(filename);
        std::cout << "Loaded " << ibpmatlist.size() << " matrices" << std::endl;
        
        if (!ibpmatlist.empty()) {
            std::cout << "First matrix: ne=" << ibpmatlist[0].ne 
                      << ", nb=" << ibpmatlist[0].nb 
                      << ", nibp=" << ibpmatlist[0].nibp << std::endl;
        }
        
        std::cout << "\n=== Success ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
