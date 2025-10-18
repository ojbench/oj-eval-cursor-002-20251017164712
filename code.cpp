#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// ====== Begin Interface (from int2048.h) ======
#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

namespace sjtu {
class int2048 {
private:
  static const unsigned int BASE = 10000U; // 1e4 base per limb
  std::vector<unsigned int> limbs;              // little-endian limbs
  bool negative = false;                        // sign flag, false means non-negative

  // remove leading zero limbs and fix sign if zero
  void normalize();

  // absolute-value helpers (do not touch sign)
  static int compareAbs(const int2048 &a, const int2048 &b);
  static int2048 addAbs(const int2048 &a, const int2048 &b);
  // precondition: |a| >= |b|
  static int2048 subAbs(const int2048 &a, const int2048 &b);
  static int2048 mulAbs(const int2048 &a, const int2048 &b);
  // returns (q, r) with 0 <= r < |b|, both non-negative (absolute division)
  static void divModAbs(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r);

public:
  // 构造函数
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // ===================================
  // Integer1
  // ===================================

  // 读入一个大整数
  void read(const std::string &);
  // 输出储存的大整数，无需换行
  void print();

  // 加上一个大整数
  int2048 &add(const int2048 &);
  // 返回两个大整数之和
  friend int2048 add(int2048, const int2048 &);

  // 减去一个大整数
  int2048 &minus(const int2048 &);
  // 返回两个大整数之差
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif
// ====== End Interface ======

// ====== Begin Implementation (from int2048.cpp) ======
namespace sjtu {

using cd = std::complex<double>;
static void fft(std::vector<cd> &a, bool invert);

// Normalize: remove leading zeros; if zero, make sign non-negative
void int2048::normalize() {
  while (!limbs.empty() && limbs.back() == 0) limbs.pop_back();
  if (limbs.empty()) negative = false;
}

int int2048::compareAbs(const int2048 &a, const int2048 &b) {
  if (a.limbs.size() != b.limbs.size()) return a.limbs.size() < b.limbs.size() ? -1 : 1;
  for (int i = (int)a.limbs.size() - 1; i >= 0; --i) {
    if (a.limbs[i] != b.limbs[i]) return a.limbs[i] < b.limbs[i] ? -1 : 1;
  }
  return 0;
}

int2048 int2048::addAbs(const int2048 &a, const int2048 &b) {
  int2048 res;
  const unsigned int BASE = int2048::BASE;
  const size_t n = a.limbs.size();
  const size_t m = b.limbs.size();
  const size_t L = n > m ? n : m;
  res.limbs.resize(L);
  unsigned int carry = 0;
  for (size_t i = 0; i < L; ++i) {
    unsigned long long sum = carry;
    if (i < n) sum += a.limbs[i];
    if (i < m) sum += b.limbs[i];
    res.limbs[i] = (unsigned int)(sum % BASE);
    carry = (unsigned int)(sum / BASE);
  }
  if (carry) res.limbs.push_back(carry);
  return res;
}

int2048 int2048::subAbs(const int2048 &a, const int2048 &b) {
  // pre: |a| >= |b|
  int2048 res;
  const unsigned int BASE = int2048::BASE;
  const size_t n = a.limbs.size();
  const size_t m = b.limbs.size();
  res.limbs.resize(n);
  long long borrow = 0;
  for (size_t i = 0; i < n; ++i) {
    long long cur = (long long)a.limbs[i] - (i < m ? b.limbs[i] : 0) - borrow;
    if (cur < 0) {
      cur += BASE;
      borrow = 1;
    } else {
      borrow = 0;
    }
    res.limbs[i] = (unsigned int)cur;
  }
  res.normalize();
  return res;
}

static void fft(std::vector<cd> &a, bool invert) {
  const size_t n = a.size();
  for (size_t i = 1, j = 0; i < n; ++i) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = (invert ? -1.0 : 1.0) * (2.0) * 3.14159265358979323846 / (double)len;
    cd wlen = std::polar(1.0, ang);
    for (size_t i = 0; i < n; i += len) {
      cd w(1.0, 0.0);
      for (size_t j = 0; j < (len >> 1); ++j) {
        cd u = a[i + j];
        cd v = a[i + j + (len >> 1)] * w;
        a[i + j] = u + v;
        a[i + j + (len >> 1)] = u - v;
        w *= wlen;
      }
    }
  }
  if (invert) {
    for (size_t i = 0; i < n; ++i) a[i] /= (double)n;
  }
}

int2048 int2048::mulAbs(const int2048 &a, const int2048 &b) {
  int2048 res;
  if (a.limbs.empty() || b.limbs.empty()) return res; // zero
  const unsigned int BASE = int2048::BASE;
  std::vector<cd> fa, fb;
  size_t n = 1;
  size_t need = a.limbs.size() + b.limbs.size();
  while (n < need) n <<= 1;
  fa.resize(n);
  fb.resize(n);
  for (size_t i = 0; i < n; ++i) {
    fa[i] = cd(i < a.limbs.size() ? (double)a.limbs[i] : 0.0, 0.0);
    fb[i] = cd(i < b.limbs.size() ? (double)b.limbs[i] : 0.0, 0.0);
  }
  fft(fa, false);
  fft(fb, false);
  for (size_t i = 0; i < n; ++i) fa[i] *= fb[i];
  fft(fa, true);
  res.limbs.resize(need);
  unsigned long long carry = 0ULL;
  for (size_t i = 0; i < need; ++i) {
    double val = fa[i].real();
    if (val < 0) val = 0; // guard tiny negative due to precision
    unsigned long long t = (unsigned long long)(val + 0.5) + carry;
    res.limbs[i] = (unsigned int)(t % BASE);
    carry = t / BASE;
  }
  while (carry) {
    res.limbs.push_back((unsigned int)(carry % BASE));
    carry /= BASE;
  }
  res.normalize();
  return res;
}

void int2048::divModAbs(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r) {
  q.limbs.clear(); q.negative = false;
  r.limbs.clear(); r.negative = false;
  if (b.limbs.empty()) return; // undefined
  if (a.limbs.empty()) return; // zero
  int cmp = compareAbs(a, b);
  if (cmp < 0) { r = a; return; }
  if (b.limbs.size() == 1) {
    unsigned int d = b.limbs[0];
    q.limbs.resize(a.limbs.size());
    unsigned long long rem = 0ULL;
    for (int i = (int)a.limbs.size() - 1; i >= 0; --i) {
      unsigned long long cur = a.limbs[i] + rem * (unsigned long long)BASE;
      q.limbs[i] = (unsigned int)(cur / d);
      rem = cur % d;
    }
    q.normalize();
    if (rem) r.limbs.push_back((unsigned int)rem);
    r.normalize();
    return;
  }
  const unsigned int BASE = int2048::BASE;
  const size_t n = a.limbs.size();
  const size_t m = b.limbs.size();
  std::vector<unsigned int> u(n + 1, 0);
  for (size_t i = 0; i < n; ++i) u[i] = a.limbs[i];
  std::vector<unsigned int> v = b.limbs;
  q.limbs.assign(n - m + 1, 0);

  unsigned int norm = (unsigned int)(BASE / ((unsigned long long)v[m - 1] + 1ULL));
  if (norm > 1) {
    unsigned long long carry = 0ULL;
    for (size_t i = 0; i < n; ++i) {
      unsigned long long cur = (unsigned long long)u[i] * norm + carry;
      u[i] = (unsigned int)(cur % BASE);
      carry = cur / BASE;
    }
    u[n] = (unsigned int)carry;
    carry = 0ULL;
    for (size_t i = 0; i < m; ++i) {
      unsigned long long cur = (unsigned long long)v[i] * norm + carry;
      v[i] = (unsigned int)(cur % BASE);
      carry = cur / BASE;
    }
  }

  for (int k = (int)(n - m); k >= 0; --k) {
    unsigned long long u2 = (unsigned long long)u[k + m] * BASE + u[k + m - 1];
    unsigned long long qhat = u2 / v[m - 1];
    unsigned long long rhat = u2 % v[m - 1];
    if (qhat >= BASE) { qhat = BASE - 1; rhat += v[m - 1]; }
    while (m >= 2 && qhat * v[m - 2] > rhat * BASE + u[k + m - 2]) {
      --qhat;
      rhat += v[m - 1];
      if (rhat >= BASE) break;
    }

    long long borrow = 0;
    unsigned long long carry = 0ULL;
    for (size_t j = 0; j < m; ++j) {
      unsigned long long prod = (unsigned long long)v[j] * (unsigned long long)qhat + carry;
      long long cur = (long long)u[k + j] - (long long)(prod % BASE) - borrow;
      carry = prod / BASE;
      if (cur < 0) { cur += BASE; borrow = 1; } else { borrow = 0; }
      u[k + j] = (unsigned int)cur;
    }
    long long cur = (long long)u[k + m] - (long long)carry - borrow;
    if (cur < 0) {
      --qhat;
      unsigned long long c = 0ULL;
      for (size_t j = 0; j < m; ++j) {
        unsigned long long sum = (unsigned long long)u[k + j] + (unsigned long long)v[j] + c;
        u[k + j] = (unsigned int)(sum % BASE);
        c = sum / BASE;
      }
      u[k + m] = (unsigned int)((unsigned long long)u[k + m] + c);
    } else {
      u[k + m] = (unsigned int)cur;
    }
    q.limbs[k] = (unsigned int)qhat;
  }

  r.limbs.assign(m, 0);
  for (size_t i = 0; i < m; ++i) r.limbs[i] = u[i];
  if (norm > 1) {
    unsigned long long carry = 0ULL;
    for (int i = (int)m - 1; i >= 0; --i) {
      unsigned long long cur = (unsigned long long)r.limbs[i] + carry * BASE;
      r.limbs[i] = (unsigned int)(cur / norm);
      carry = cur % norm;
    }
  }
  q.normalize();
  r.normalize();
}

int2048::int2048() {}

int2048::int2048(long long val) {
  if (val < 0) { negative = true; val = -val; } else { negative = false; }
  const unsigned int BASE = int2048::BASE;
  while (val) {
    limbs.push_back((unsigned int)(val % BASE));
    val /= BASE;
  }
  normalize();
}

int2048::int2048(const std::string &s) { read(s); }

int2048::int2048(const int2048 &other) { limbs = other.limbs; negative = other.negative; }

void int2048::read(const std::string &s) {
  limbs.clear(); negative = false;
  if (s.empty()) return;
  size_t i = 0; bool neg = false;
  if (s[0] == '+') { i = 1; }
  else if (s[0] == '-') { neg = true; i = 1; }
  while (i < s.size() && s[i] == '0') ++i;
  const unsigned int BASE = int2048::BASE;
  const unsigned int WIDTH = (BASE == 1000U) ? 3U : (BASE == 10000U ? 4U : 9U);
  if (i == s.size()) { negative = false; return; }
  std::vector<unsigned int> parts;
  size_t n = s.size();
  for (size_t j = n; j > i;) {
    size_t start = (j >= WIDTH) ? j - WIDTH : i;
    if (start < i) start = i;
    unsigned int x = 0;
    for (size_t k = start; k < j; ++k) x = x * 10 + (unsigned int)(s[k] - '0');
    parts.push_back(x);
    j = start;
    if (j == i) break;
  }
  limbs = parts;
  negative = neg;
  normalize();
}

void int2048::print() {
  if (limbs.empty()) { std::cout << 0; return; }
  if (negative) std::cout << '-';
  const unsigned int BASE = int2048::BASE;
  const unsigned int WIDTH = (BASE == 1000U) ? 3U : (BASE == 10000U ? 4U : 9U);
  std::cout << limbs.back();
  for (int i = (int)limbs.size() - 2; i >= 0; --i) {
    unsigned int x = limbs[i];
    unsigned int pow10 = 1;
    for (unsigned int t = 1; t < WIDTH; ++t) pow10 *= 10U;
    while (pow10) {
      unsigned int digit = x / pow10; x %= pow10; pow10 /= 10U;
      std::cout << (char)('0' + digit);
    }
  }
}

int2048 &int2048::add(const int2048 &rhs) {
  if (rhs.limbs.empty()) return *this;
  if (limbs.empty()) { *this = rhs; return *this; }
  if (negative == rhs.negative) {
    *this = addAbs(*this, rhs);
    this->negative = negative;
  } else {
    int cmp = compareAbs(*this, rhs);
    if (cmp == 0) {
      limbs.clear(); negative = false;
    } else if (cmp > 0) {
      *this = subAbs(*this, rhs);
    } else {
      *this = subAbs(rhs, *this);
      this->negative = rhs.negative;
    }
  }
  return *this;
}

int2048 add(int2048 a, const int2048 &b) { return a.add(b); }

int2048 &int2048::minus(const int2048 &rhs) {
  if (rhs.limbs.empty()) return *this;
  int2048 tmp = rhs; tmp.negative = !tmp.negative;
  return this->add(tmp);
}

int2048 minus(int2048 a, const int2048 &b) { return a.minus(b); }

int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const {
  int2048 t(*this);
  if (!t.limbs.empty()) t.negative = !t.negative;
  return t;
}

int2048 &int2048::operator=(const int2048 &rhs) {
  if (this == &rhs) return *this;
  limbs = rhs.limbs; negative = rhs.negative; return *this;
}

int2048 &int2048::operator+=(const int2048 &rhs) { return this->add(rhs); }
int2048 operator+(int2048 a, const int2048 &b) { return a.add(b); }

int2048 &int2048::operator-=(const int2048 &rhs) { return this->minus(rhs); }
int2048 operator-(int2048 a, const int2048 &b) { return a.minus(b); }

int2048 &int2048::operator*=(const int2048 &rhs) {
  if (limbs.empty() || rhs.limbs.empty()) { limbs.clear(); negative = false; return *this; }
  bool sign = (negative != rhs.negative);
  int2048 prod = mulAbs(*this, rhs);
  *this = prod; this->negative = sign && !this->limbs.empty();
  return *this;
}
int2048 operator*(int2048 a, const int2048 &b) { a *= b; return a; }

int2048 &int2048::operator/=(const int2048 &rhs) {
  if (rhs.limbs.empty()) return *this;
  if (limbs.empty()) { negative = false; return *this; }
  int2048 aAbs = *this; aAbs.negative = false;
  int2048 bAbs = rhs;  bAbs.negative = false;
  int2048 q, r;
  divModAbs(aAbs, bAbs, q, r);
  bool signsDifferent = (negative != rhs.negative);
  if (!signsDifferent) {
    *this = q; this->negative = negative; this->normalize();
  } else {
    if (r.limbs.empty()) {
      *this = q; this->negative = true; this->normalize();
    } else {
      q = addAbs(q, int2048(1));
      *this = q; this->negative = true; this->normalize();
    }
  }
  return *this;
}
int2048 operator/(int2048 a, const int2048 &b) { a /= b; return a; }

int2048 &int2048::operator%=(const int2048 &rhs) {
  if (rhs.limbs.empty()) return *this;
  if (limbs.empty()) { negative = false; return *this; }
  int2048 q = (*this) / rhs;
  int2048 prod = q * rhs;
  *this = ::sjtu::minus(*this, prod);
  return *this;
}
int2048 operator%(int2048 a, const int2048 &b) { a %= b; return a; }

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s; is >> s; x.read(s); return is;
}
std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.limbs.empty()) { os << 0; return os; }
  if (x.negative) os << '-';
  const unsigned int BASE = int2048::BASE;
  const unsigned int WIDTH = (BASE == 1000U) ? 3U : (BASE == 10000U ? 4U : 9U);
  os << x.limbs.back();
  for (int i = (int)x.limbs.size() - 2; i >= 0; --i) {
    unsigned int xk = x.limbs[i];
    unsigned int pow10 = 1;
    for (unsigned int t = 1; t < WIDTH; ++t) pow10 *= 10U;
    while (pow10) {
      unsigned int digit = xk / pow10; xk %= pow10; pow10 /= 10U;
      os << (char)('0' + digit);
    }
  }
  return os;
}

bool operator==(const int2048 &a, const int2048 &b) {
  if (a.negative != b.negative) return a.limbs.empty() && b.limbs.empty();
  return a.limbs == b.limbs;
}
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) {
  if (a.negative != b.negative) return a.negative && (!a.limbs.empty() || !b.limbs.empty());
  int cmp = int2048::compareAbs(a, b);
  if (!a.negative) return cmp < 0;
  return cmp > 0;
}
bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu
// ====== End Implementation ======
