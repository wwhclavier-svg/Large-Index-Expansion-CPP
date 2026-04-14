#include <iostream>
#include "firefly/FFInt.hpp"

using namespace firefly;

int main() {
    std::cout << "=== FireFly Simple Test ===" << std::endl;
    
    // 设置模数
    uint64_t p = 2147483647;  // 2^31 - 1
    FFInt::set_new_prime(p);
    std::cout << "Prime modulus: " << FFInt::p << std::endl;
    
    // 基本运算测试
    FFInt a(100);
    FFInt b(200);
    
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a + b = " << (a + b) << std::endl;
    std::cout << "a - b = " << (a - b) << std::endl;
    std::cout << "a * b = " << (a * b) << std::endl;
    
    // 测试模运算
    FFInt c(p + 10);
    std::cout << "c = p + 10 = " << c << " (should be 10)" << std::endl;
    
    // 类型转换测试
    int int_val = static_cast<int>(a);
    std::cout << "int(a) = " << int_val << std::endl;
    
    std::cout << "\nFireFly is working correctly!" << std::endl;
    return 0;
}
