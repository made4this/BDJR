#include "pcmax.h"

#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace PCmax;

inline void trim(std::string& s, const char* t = " \t\n\r\f\v") {
  s.erase(s.find_last_not_of(t) + 1);
  s.erase(0, s.find_first_not_of(t));
}

bool Instance::Read(const char* file) {

  this->map.clear();

  std::ifstream filein(file);
  std::string line;

  if (!filein.good()) return false;

  getline(filein, line);
  trim(line);
  this->m = std::atoi(line.c_str());

  assert(filein.good());
  getline(filein, line);

  while (filein.good()) {
    getline(filein, line);
    trim(line);
    if (line.empty()) continue;
    nat p = std::atoi(line.c_str());
    ++this->map[p];
  }

  return true;
}

std::ostream& PCmax::operator<<(std::ostream& strm, const Instance& I) {
  strm << "Instance[m = " << I.m << ", jobs (processing time, amount): ";

  std::map<nat,nat>::const_iterator it;
  for (it = I.map.begin(); it != I.map.end(); ++it) {
    strm << "(" << it->first << ", " << it->second << ")";
    if (std::next(it) != I.map.end()) strm << ", ";
  }

  return strm << "]";
}

std::ostream& PCmax::operator<<(std::ostream& strm, const Schedule& S) {
  strm << "[";
  nat i = 0;
  for (auto const & conf : S) {
    strm << "(";
    nat j = 0;
    for (nat job : conf) {
      strm << job;
      if (++j < conf.size()) strm << ",";
    }
    strm << ")";
    if (++i < S.size()) strm << ",";
  }
  strm << "]";
  return strm;
}

nat PCmax::LowerBound(const Instance& I) {

  auto map = I.GetMap();
  nat p_max = map.rbegin()->first;
  nat p_m   = 0; // p_m
  nat p_mp1 = 0; // p_{m+1}
  nat m = I.GetM();
  nat n = I.GetN();

  if (n <= m) return p_max;

  nat k = 0;
  nat P = 0;
  for (auto i = map.rbegin(); i != map.rend(); i++) {
    nat p = i->first;
    nat a = i->second;
    P += a*p;
    if (k < m && m <= k+a) p_m = p;
    if (k < m+1 && m+1 <= k+a) p_mp1 = p;
    k += a;
  }

  nat ceiled_avg_m = (P-1) / m + 1;

  return std::max({ceiled_avg_m, p_max, p_m + p_mp1});
}
