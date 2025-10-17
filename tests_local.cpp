#include <bits/stdc++.h>
#include "src/include/int2048.h"
using namespace sjtu;

int main(){
  int2048 a("0"), b("0");
  a.read("123456789012345678901234567890");
  b.read("987654321098765432109876543210");
  std::cout << (a + b) << "\n"; // 1111111110111111111011111111100
  std::cout << (b - a) << "\n"; // 864197532086419753208641975320
  std::cout << (a * b) << "\n"; // product
  int2048 c("-10"), d("3");
  std::cout << (c / d) << "\n"; // -4 floor div
  std::cout << (c % d) << "\n"; // 2, since -10 = (-4)*3 + 2
  // stress small
  int2048 x("9999999999999999999999999999999999");
  int2048 y("1");
  std::cout << (x + y) << "\n";
  return 0;
}
