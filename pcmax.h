#ifndef PCMAX_H_
#define PCMAX_H_

#include <map>
#include <ostream>
#include <vector>
#include <iostream>

namespace PCmax {

typedef size_t nat; // declare a natural number

class Instance {

 private:
  nat m;
  std::map<nat,nat> map;

  friend std::ostream& operator<<(std::ostream&, const Instance&);

 public:
  Instance() {
    m = 0;
  }

  Instance(nat m, std::map<nat,nat>& jobs) : map(jobs) {
    this->m = m;
  }

  Instance(nat m, nat n, nat a[], nat p[]) {
    this->m = m;
    for (nat i = 0; i < n; ++i) map[p[i]] = a[i];
  }

  /* Read a file of the content form: m \n n \n p1 \n p2 \n p3 ... */
  bool Read(const char* file);

  nat GetM() const {
	  return m;
  }

  void SetM(nat m) {
    this->m = m;
  }

  nat GetN() const {
    nat n = 0;
    for (auto& j : map) n += j.second;
    return n;
  }

  const std::map<nat,nat>& GetMap() const {
    return map;
  }

  std::map<nat,nat>& GetMap() {
	  return map;
  }

  void Clear() {
	  this->m = 0;
	  this->map.clear();
  }
};

std::ostream& operator<<(std::ostream&, const Instance&);

class Schedule : public std::vector<std::vector<nat>> {

 private:
  friend std::ostream& operator<<(std::ostream&, const Schedule&);

 public:
  nat ComputeMakespan() {
    nat C_max = 0;
    for (auto const& conf : *this) {
      nat C = 0;
      for (auto const& load : conf) C += load;
      if (C > C_max) C_max = C;
    }
    return C_max;
  }

  void AddLoad(nat machine, nat load) {
    if (size() <= machine) resize(machine+1);
    at(machine).push_back(load);
  }

  bool IsFeasibleForInstance(const Instance& I) {
    if (size() > I.GetM()) return false;
    std::map<nat,nat> jobs(I.GetMap()); // copy
    for (auto const& conf : *this) {
      for (nat load : conf) {
        if (jobs.find(load) == jobs.end() || jobs[load] == 0) return false;
        --jobs[load];
      }
    }
    for (auto const& j : jobs) if (j.second > 0) return false;
    return true;
  }

};

std::ostream& operator<<(std::ostream&, const Schedule&);

nat LowerBound(const Instance& I);

}

#endif // PCMAX_H_
