#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include "rounding.h"
#include "scheduler.h"
#include "convolution.h"
#include "heuristic.h"
#include "bdjr.h"

using namespace PCmax;

typedef Rounding9<0> R;

inline nat ScheduleHugeJobs(double eps, nat T, std::map<nat,nat>& jobs, Schedule& S) {

  nat huge_machines = 0;

  for (auto i = jobs.rbegin(); i != jobs.rend(); ++i) {
    nat p = i->first;
    nat a = i->second;
    double p1 = p / (double)T; // normalized

    if (p1 < 1-2*eps) break;

    nat remove = a;
    std::map<nat,nat>::reverse_iterator j(i); ++j;
    while (j != jobs.rend() && remove > 0) {
      nat p_ = j->first;
      nat a_ = j->second;

      double p1_ = p_ / (double)T; // normalized
      if (p1_ <= eps) break;

      if (p <= T && p_ <= T-p) {// largest possible job(s)
        nat rem = std::min(a_, remove);
        for (nat k = 0; k < rem; ++k) S.push_back(std::vector<nat>{p,p_});
        huge_machines += rem;
        if (a_ <= remove) {
          jobs.erase(std::next(j).base()); // remove item j
        } else {
          j->second -= rem;
        }
        remove -= rem;
      }
      else ++j;
    }
    huge_machines += remove;
    while (remove-- > 0) {
      // not enough medium jobs - store them alone!
      S.push_back(std::vector<nat>{p});
    }
  }

  return huge_machines;
}

inline bool DualMakespanTask(double eps, const Instance& I, nat T, nat& m_min) {
  nat m = I.GetM();
  std::map<nat,nat> jobs(I.GetMap()); // copy

  nat huge_machines;
  { Schedule S; huge_machines = ScheduleHugeJobs(eps, T, jobs, S); }
  if (huge_machines > m) return false;

  nat m_med = m - huge_machines; // machines for medium jobs

  Vector<R> b;
  b.Reset();

  for (auto i = jobs.rbegin(); i != jobs.rend(); i++) {
  	nat p = i->first;
  	nat a = i->second;
  	double p1 = p / (double)T; // normalized
    //std::cout << p1 << "," << std::flush;
  	if (p1 <= eps) break;
    if (p1 >= 1-2*eps) continue;
    b[R::GetRoundedIndex(p1)] += a;
  }

  m_min = Scheduler<R, FFTConvolution<R>>::MinMachines(b);
  return (m_min <= m_med);
}

inline nat ComputeFirstMakespan(double eps, const Instance& I, nat& m_min) {

  nat l = LowerBound(I);
  nat u = MF::ComputeMakespan(I);
  bool ok = false;
  do {
    nat T = (l + u) / 2;
    nat m;
    if (DualMakespanTask(eps, I, T, m)) {
      ok = true;
      u = T;
      m_min = m;
    } else l = T+1;
  } while (l < u);

  if (!ok) ok = DualMakespanTask(eps, I, u, m_min);
  if (!ok) std::cout << "NOT OK!" << std::endl;

  return u;
}

inline void RoundMediumJobs(double eps, nat T, const std::map<nat,nat>& jobs, Vector<R>& b) {
  for (auto i = jobs.rbegin(); i != jobs.rend(); ++i) {
    nat p = i->first;
    nat a = i->second;
    double p1 = p / (double)T; // normalized

    if (p1 >= 1-2*eps) continue;
    if (p1 <= eps) break;

    b[R::GetRoundedIndex(p1)] += a;
  }
}

inline void UnroundScheduleOfMediumJobs(double eps, nat T, nat huge_machines,
  const std::map<nat,nat>& jobs, std::vector<Vector<R>>& S_, Schedule& S) {
  std::vector<nat> loads(S_.size());
  for (auto i = jobs.rbegin(); i != jobs.rend(); i++) {
    nat p = i->first;
    nat a = i->second;
    double p1 = p / (double)T; // normalized

    if (p1 >= 1-2*eps) continue;
    if (p1 <= eps) break;

    int j = R::GetRoundedIndex(p1);
    for (nat k = 0; k < a; k++) {
      nat u = huge_machines;
      double minload = T;
      nat u_minload = u;
      for (auto& c : S_) {
        nat load = loads[u-huge_machines];
        if (c[j] > 0 && load < minload) {
          u_minload = u;
          minload = load;
        }
        ++u;
      }
      --S_[u_minload-huge_machines][j];
      S.AddLoad(u_minload, p);
      loads[u_minload-huge_machines] += minload;
    }
  }
}

inline void GetSmallJobs(double eps, nat T, const std::map<nat,nat>& jobs, std::map<nat,nat>& small_jobs) {
  for (auto i = jobs.begin(); i != jobs.end(); i++) {
    nat p = i->first;
    nat a = i->second;
    double p1 = p / (double)T; // normalized
    if (p1 > eps) break;
    small_jobs[p] = a;
  }
}

/*
 * The main algorithm.
 */
nat BDJR::ComputeSchedule(const Instance& I, Schedule& S) {

  //double eps = 0.172874755859;
  double eps = 0.1754019165039063;
  nat m = I.GetM();
  std::map<nat,nat> jobs(I.GetMap());// copy

  nat m_min = m;
  nat T = ComputeFirstMakespan(eps, I, m_min);
  std::cout << "First Makespan: T = " << T << std::endl;

  nat huge_machines = ScheduleHugeJobs(eps, T, jobs, S);
  std::cout << "AfterScheduleHugeJobs: S = " << S << std::endl;

  Vector<R> b;
  b.Reset();
  RoundMediumJobs(eps, T, jobs, b);
  std::cout << "AfterRoundMediumJobs: b = " << b << std::endl;

  std::vector<Vector<R>> S_;
  Scheduler<R, FFTConvolution<R>>::ComputeSchedule(b, m_min, S_);
  std::cout << "AfterComputeScheduleBeforeUnround: S_ = "; for (auto& c : S_) std::cout << c << ", "; std::cout << std::endl;

  for (nat i = 0; i < S_.size(); i++) {
    std::vector<nat> u;
    S.push_back(u);
  }

  UnroundScheduleOfMediumJobs(eps, T, huge_machines, jobs, S_, S);
  std::cout << "AfterUnroundBeforeLPT: S = " << S << std::endl;

  std::map<nat,nat> small_jobs;
  GetSmallJobs(eps, T, jobs, small_jobs);
  Instance I_small(m, small_jobs);
  std::cout << I_small << std::endl;

  return LPT::ComputeSchedule(I_small, S);// schedule small jobs on top of S
}
