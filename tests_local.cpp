#include <bits/stdc++.h>
#include "src/include/int2048.h"
using namespace sjtu;

static void check(const std::string &expr, const std::string &got, const std::string &exp){
  if(got!=exp){
    std::cerr<<"FAIL "<<expr<<" got="<<got<<" exp="<<exp<<"\n";
    exit(1);
  }
}

static std::string tostr(const int2048 &x){
  std::ostringstream os; os<<x; return os.str();
}

int main(){
  {
    int2048 a("123456789012345678901234567890");
    int2048 b("987654321098765432109876543210");
    check("a+b", tostr(a+b), "1111111110111111111011111111100");
    check("b-a", tostr(b-a), "864197532086419753208641975320");
  }
  {
    int2048 c("-10"), d("3");
    check("c/d", tostr(c/d), "-4");
    check("c%d", tostr(c%d), "2");
  }
  {
    int2048 x("9999999999999999999999999999999999");
    int2048 y("1");
    check("x+y", tostr(x+y), "10000000000000000000000000000000000");
  }
  {
    int2048 a("0"), b("123");
    check("0-123", tostr(a.minus(b)), "-123");
  }
  {
    int2048 a("-1000"), b("1");
    check("-1000+1", tostr(a.add(b)), "-999");
    a.read("-1000");
    check("-1000-1", tostr(a.minus(b)), "-1001");
  }
  {
    int2048 a("1"), b("-1");
    check("1-(-1)", tostr(a.minus(b)), "2");
    a.read("1");
    check("1+(-1)", tostr(a.add(b)), "0");
  }
  std::cout<<"OK\n";
  return 0;
}
