#include <iostream>
#include <filesystem>
#include <chrono>

#include "rounding.h"
#include "scheduler.h"
#include "convolution.h"
#include "pcmax.h"
#include "bdjr.h"
#include "heuristic.h"

typedef Rounding9<0> R;

using namespace PCmax;
using namespace std;

inline void ScaleInstance(Instance& I, nat scaleM, nat scaleN) {
  I.SetM(I.GetM() * scaleM);
  auto& jobs = I.GetMap();
  for (auto j = jobs.begin(); j != jobs.end(); j++) {
    j->second *= scaleN;
  }
}

inline nat RunExperiment(const char* name, nat (*ComputeSchedule)(const Instance&, Schedule&),
  const Instance& I, Schedule& S) {

  S.clear();

  auto start = chrono::steady_clock::now();

  nat makespan = ComputeSchedule(I, S);

  auto stop = chrono::steady_clock::now();
  auto ms = chrono::duration_cast<chrono::milliseconds>(stop-start);

  const char* ok = S.IsFeasibleForInstance(I) ? "ok" : "FAIL";

  cout << name << ": " << makespan << " " << ms.count() << " " << ok << S << endl;

  return makespan;
}

int main(int argc, const char* argv[]) {

  nat scaleM = std::atoi(argv[1]);
  nat scaleN = std::atoi(argv[2]);

  for (int i = 3; i < argc; i++) {
    const char* file = argv[i];

    cout << filesystem::path(file).stem().string() << endl;

    Instance I;
    Schedule S;

    if (I.Read(file)) {
      ScaleInstance(I, scaleM, scaleN);
      cout << I << endl;

      RunExperiment("LPT",  &LPT::ComputeSchedule,  I, S);
      RunExperiment("MF",   &MF::ComputeSchedule,   I, S);
      nat djms = RunExperiment("DJMS", &DJMS::ComputeSchedule, I, S);
      nat bdjr = RunExperiment("BDJR", &BDJR::ComputeSchedule, I, S);

      double eps = 0.1754019165039063;
      if (bdjr > (1+eps)*djms) {
        cout << "BAD MAKESPAN" << endl;
      }
    }
  }

  return 0;
}
