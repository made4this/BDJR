#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <numeric>

#include "heuristic.h"

using namespace PCmax;

/*
* Returns the makespan of the corresponding LPT schedule.
*/
nat LPT::ComputeMakespan(const Instance& I) {
  nat m = I.GetM();
  auto map = I.GetMap();
  std::vector<nat> mchns(m, 0);
  nat max_load = 0;
  for (auto i = map.rbegin(); i != map.rend(); i++) {
    nat p = i->first;
    nat a = i->second;
    for (nat j = 0; j < a; j++) {
      nat u_min = 0;
      nat min_load = mchns[0];
      for (nat u = 1; u < m; u++) {
        nat load = mchns[u];
        if (load < min_load) {
          u_min = u;
          min_load = load;
        }
      }
      nat new_load = min_load + p;
      if (new_load > max_load) {
        max_load = new_load;
      }
      mchns[u_min] = new_load;
    }
  }
  return max_load;
}

nat LPT::ComputeSchedule(const Instance& I, Schedule& S) {
  nat m = I.GetM();
  auto map = I.GetMap();
  std::vector<nat> mchns(m, 0);
  // find machines loads
  nat i = 0;
  nat max_load = 0;
  for (auto const& conf : S) {
    nat C = 0;
    for (auto const& load : conf) C += load;
    if (C > max_load) max_load = C;
    mchns[i++] = C;
  }
  for (auto i = map.rbegin(); i != map.rend(); i++) {
    nat p = i->first;
    nat a = i->second;
    for (nat j = 0; j < a; j++) {
      nat u_min = 0;
      nat min_load = mchns[0];
      for (nat u = 1; u < m; u++) {
        nat load = mchns[u];
        if (load < min_load) {
          u_min = u;
          min_load = load;
        }
      }
      nat new_load = min_load + p;
      if (new_load > max_load) {
        max_load = new_load;
      }
      mchns[u_min] = new_load;
      S.AddLoad(u_min, p); // schedule job
    }
  }
  return max_load;
}

/*
 * First Fit Decreasing (Dual approach trying makespan T)
 */
static bool FFD(const Instance& I, nat T) {
  nat m = I.GetM();
  auto const & map = I.GetMap();
  std::vector<nat> mchns(m, 0);
  for (auto i = map.rbegin(); i != map.rend(); i++) {
    nat p = i->first;
    nat a = i->second;
    for (nat u = 0; u < m; u++) {
      nat max = (T-mchns[u])/p;
      if (max >= a) {
        mchns[u] += a*p;
        a = 0;
        break;
      } else {
        mchns[u] += max*p;
        a -= max;
      }
    }
    if (a > 0) return false;
  }
  return true;
}

nat MF::ComputeMakespan(const Instance& I) {
  nat l = LowerBound(I);
  nat u = LPT::ComputeMakespan(I);
  std::cout << "Bounds: l = " << l << ", u = " << u << std::endl;

  while (l != u) {
	  nat T = (l + u) / 2;
	  //std::cout << "Testing T = " << T << std::endl;
	  if (FFD(I, T)) u = T; else l = T+1;
  }

  return u;
}

nat MF::ComputeSchedule(const Instance& I, Schedule& S) {
  nat T = MF::ComputeMakespan(I);

  if (!FFD(I, T)) {
    // MF cannot improve upon LPT!
    return LPT::ComputeSchedule(I, S);
  }

  nat m = I.GetM();
  auto const & map = I.GetMap();
  std::vector<nat> mchns(m, 0);
  for (auto i = map.rbegin(); i != map.rend(); i++) {
    nat p = i->first;
    nat a = i->second;
    for (nat u = 0; u < m; u++) {
      nat max = (T-mchns[u])/p;
      if (max >= a) {
        mchns[u] += a*p;
        for (nat l = 0; l < a; l++) S.AddLoad(u, p); // schedule job
        a = 0;
        break;
      } else {
        mchns[u] += max*p;
        for (nat l = 0; l < max; l++) S.AddLoad(u, p); // schedule job
        a -= max;
      }
    }
  }

  return T;
}

nat DJMS::ComputeSchedule(const Instance& I, Schedule& S) {

  Schedule S_a, S_c;
  Instance I_a(I); // copy of I

  nat m = I.GetM();
  nat L_2 = LowerBound(I);
  nat T_c = L_2;
  nat T = LPT::ComputeMakespan(I) + 1;
  nat mac_c = 0;
  nat T_a;

  do {
    // Apply Multifit to schedule the active jobs I_a on the active machines S_a
    S_a.clear();
    T_a = MF::ComputeSchedule(I_a, S_a);
    std::map<nat,nat>& J_a = I_a.GetMap();

    // Store incumbent solution
    if (T > T_a && T > T_c) {
      // TODO seems unneccessary to do this every round
      S.clear();
      S.insert(S.end(), S_a.begin(), S_a.end());
      S.insert(S.end(), S_c.begin(), S_c.end());
      T = (T_a > T_c) ? T_a : T_c;
    }

    // Find the smallest machine load l in S_a that holds l >= L_2
    nat l = T_a;
    for (auto const & conf : S_a) {
      nat load = std::accumulate(conf.begin(), conf.end(), 0);
	    if (load >= L_2 && load < l) l = load;
	  }

    // Close all active machines with load l
    for (auto const & conf : S_a) {
      nat load = std::accumulate(conf.begin(), conf.end(), 0);
      if (load == l) {
        // Close this machine
        for (auto const & j : conf) J_a[j]--;
        S_c.push_back(conf);
        T_c = (load > T_c) ? load : T_c;
        I_a.SetM(I_a.GetM()-1);
        mac_c++;
      }
    }
  } while (T_a > T_c && mac_c < m);

  return T;
}
