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

inline uint64_t modPow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}

inline FFInt FFInt::operator/(const FFInt& other) const {
    // Compute modular inverse using Fermat's little theorem: a^(-1) = a^(p-2) mod p
    uint64_t inv = modPow(other.n, p - 2, p);
    return FFInt((n * inv) % p);
}

inline FFInt& FFInt::operator/=(const FFInt& other) {
    uint64_t inv = modPow(other.n, p - 2, p);
    n = (n * inv) % p;
    return *this;
}

} // namespace firefly

#endif // FFINT_HPP
