// 验证 FireFly 是否正确工作的简单测试
#include <iostream>
#include <cstdint>

// 在包含 FireFly 之前定义 uint64_t
#include <cstdint>

#include "firefly/FFInt.hpp"

using namespace firefly;

int main() {
    std::cout << "=== FireFly Verification Test ===" << std::endl;
    
    try {
        // 设置模数
        uint64_t p = 2147483647;  // 2^31 - 1，一个常用质数
        FFInt::set_new_prime(p);
        
        std::cout << "Prime modulus p = " << FFInt::p << std::endl;
        
        // 创建一些 FFInt 值
        FFInt a(100);
        FFInt b(200);
        
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        
        // 测试基本运算
        FFInt sum = a + b;
        FFInt diff = a - b;
        FFInt prod = a * b;
        
        std::cout << "a + b = " << sum << std::endl;
        std::cout << "a - b = " << diff << std::endl;
        std::cout << "a * b = " << prod << std::endl;
        
        // 验证模运算
        FFInt c(p + 10);
        std::cout << "c = p + 10 mod p = " << c << " (expected: 10)" << std::endl;
        
        // 测试除法（模逆）
        FFInt quotient = b / a;
        std::cout << "b / a = " << quotient << std::endl;
        
        // 验证：(b/a) * a 应该等于 b
        FFInt check = quotient * a;
        std::cout << "(b/a) * a = " << check << " (should equal b = " << b << ")" << std::endl;
        
        if (check == b) {
            std::cout << "\n✓ FireFly is working correctly!" << std::endl;
            return 0;
        } else {
            std::cout << "\n✗ FireFly verification failed!" << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
