#ifndef FFINT_HPP
#define FFINT_HPP

#include <cstdint>
#include <iostream>
#include <string>

namespace firefly {

// 最小化的 FFInt 实现，用于编译测试
class FFInt {
public:
    uint64_t n;
    static uint64_t p;  // 声明（定义在 .cpp 文件中）
    
    FFInt() : n(0) {}
    FFInt(uint64_t val) : n(val % p) {}
    
    static void set_new_prime(uint64_t prime) { p = prime; }
    
    FFInt operator+(const FFInt& other) const { return FFInt((n + other.n) % p); }
    FFInt operator-(const FFInt& other) const { return FFInt((n + p - other.n) % p); }
    FFInt operator*(const FFInt& other) const { return FFInt((n * other.n) % p); }
    FFInt operator/(const FFInt& other) const;
    
    FFInt& operator+=(const FFInt& other) { n = (n + other.n) % p; return *this; }
    FFInt& operator-=(const FFInt& other) { n = (n + p - other.n) % p; return *this; }
    FFInt& operator*=(const FFInt& other) { n = (n * other.n) % p; return *this; }
    FFInt& operator/=(const FFInt& other);
    
    bool operator==(const FFInt& other) const { return n == other.n; }
    bool operator!=(const FFInt& other) const { return n != other.n; }
    bool operator<(const FFInt& other) const { return n < other.n; }
    bool operator>(const FFInt& other) const { return n > other.n; }
    bool operator<=(const FFInt& other) const { return n <= other.n; }
    bool operator>=(const FFInt& other) const { return n >= other.n; }
    
    FFInt operator-() const { return FFInt((p - n) % p); }
    
    explicit operator uint64_t() const { return n; }
    explicit operator int() const { return static_cast<int>(n); }
    explicit operator double() const { return static_cast<double>(n); }
    
    friend std::ostream& operator<<(std::ostream& os, const FFInt& val) {
        os << val.n;
        return os;
    }
};

// 内联定义，避免多重定义
inline uint64_t FFInt::p = 2147483647;  // 默认质数

inline FFInt FFInt::operator/(const FFInt& other) const {
    // 简化实现，实际需要模逆
    return FFInt(n / other.n);
}

inline FFInt& FFInt::operator/=(const FFInt& other) {
    n = (n / other.n) % p;
    return *this;
}

} // namespace firefly

#endif // FFINT_HPP
