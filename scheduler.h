#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include <vector>

#include "vector.h"
#include "rounding.h"
#include "tensor.h"

template<class R, class C>
class Scheduler {
 public:
  static int MinMachines(const Vector<R>& b) {
    if (b.IsZero()) return 0;
    if (b.Volume() <= R::Makespan() - R::Precision()) return 1;
    Tensor<R> t, t_;
    Vector<R> s;
    s.Reset();
    MinMachines(t, t_, b);
    return t.Get(s);
  }

  static void MinMachines(Tensor<R>& t, Tensor<R>& trash, const Vector<R>& anchor);

  static void ComputeSchedule(const Vector<R>& b, int m, std::vector<Vector<R>>& S);
};

#endif
