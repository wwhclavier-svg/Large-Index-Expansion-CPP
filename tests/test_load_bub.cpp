#include <iostream>
#include "IBPMatrixLoader_Binary.hpp"
#include "firefly/FFInt.hpp"

int main() {
    try {
        std::cout << "=== Test Loading IBPMat_bub.bin ===" << std::endl;
        
        auto ibpmatlist = loadAllIBPMatricesBinary<firefly::FFInt>("IBPMat_bub.bin");
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
