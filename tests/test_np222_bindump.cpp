// test_np222_bindump.cpp — dump IBPMat_NP222.bin 每个 regime 的元数据
#include "../include/IBPMatrixLoader_Binary.hpp"
#include "firefly/FFInt.hpp"
#include <iostream>

int main() {
    auto regs = loadAllIBPMatricesBinary<firefly::FFInt>("data/IBPMat_NP222.bin");
    std::cout << "numRegs=" << regs.size() << std::endl;
    for (size_t r = 0; r < regs.size(); ++r) {
        const auto& m = regs[r];
        std::cout << "reg " << r << ": nibp=" << m.nibp
                  << " ne=" << m.ne << " nb=" << m.nb
                  << " incre=" << m.incre;
        // 打印各算子的顶层尺寸辅助判断
        std::cout << "  M1=" << m.M1.size()
                  << " K1s=" << m.K1s.size()
                  << " K2s=" << m.K2s.size()
                  << " F2s=" << m.F2s.size();
        std::cout << std::endl;
    }
    return 0;
}
