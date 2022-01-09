#include <iostream>
#include <algorithm>
#include <map>

#include "convolution.h"
#include "scheduler.h"

template<class R, class C>
void Scheduler<R, C>::MinMachines(Tensor<R>& t, Tensor<R>& trash, const Vector<R>& anchor) {
  Vector<R> anchor_;
  for (int i = 0; i < R::Dim(); ++i) {
    anchor_[i] = R::PreviousInterval(i, anchor[i]);
  }
  if (anchor == anchor_) {
    // Repeat until nothing changes
    t.Initialize(anchor_);
    while (1) {
      C::Square(trash, t, anchor, anchor_);
      if (t == trash) break;
      C::Square(t, trash, anchor, anchor_);
      if (t == trash) break;
    }
  } else {
    MinMachines(trash, t, anchor_);
    C::Square(t, trash, anchor, anchor_);
  }
}

template<class R>
std::vector<std::pair<Vector<R>,int>> Backtrace(const Vector<R>& anchor, const Tensor<R>& t, const std::vector<std::pair<Vector<R>,int>>& targets, std::vector<Vector<R>>& S) {

  std::vector<std::pair<Vector<R>,int>> new_targets;

  //#pragma omp parallel for
  for (const auto& target : targets) {
    const Vector<R>& b = target.first;
    int m = target.second;

    std::cout << "target b = " << b  << ", m = " << m  << std::endl;

    Vector<R> b1 = b, b2 = b;
    int m1 = m, m2 = m, bestdiff = m+1;

    //#pragma omp for nowait
    for (std::size_t i = 0; i < Tensor<R>::data_size(); ++i) {
      Vector<R> a1;
      Tensor<R>::FromIndex(a1, i);
      int v1 = t.Get(a1);
      if (v1 < 0 || v1 > m) continue;

      Vector<R> a_1 = a1;
      a_1 += anchor;
      if (a_1.IsZero()) continue;

      Vector<R> a_2 = b;
      a_2 -= a_1; // => a_1 + a_2 = b

      Vector<R> a2 = a_2;
      a2 -= anchor;
      if (!a2.IsWithinDeviation()) continue;

      int v2 = t.Get(a2);
      if (v2 < 0 || v1+v2 != m) continue;

      int diff = std::abs(v1-v2);
      if (diff < bestdiff) {
        b1 = a_1;
        b2 = a_2;
        m1 = v1;
        m2 = v2;
        bestdiff = diff;
        if (bestdiff <= 1) break;
      }
    }

    if (bestdiff == m+1) {
      std::cout << "Failed to resolve " << b << std::endl;
      continue;
    }

    std::cout << "new target b1 = " << b1 << ", m1 = " << m1 << std::endl;
    std::cout << "new target b2 = " << b2 << ", m2 = " << m2 << std::endl;

    std::pair<Vector<R>,int> t1(b1, m1);
    std::pair<Vector<R>,int> t2(b2, m2);

    //#pragma omp critical
    {
      new_targets.push_back(t1);
      new_targets.push_back(t2);
    }
  }

  std::vector<std::pair<Vector<R>,int>> nontrivial_targets;

  for (const auto& target : new_targets) {
    if (target.first.L1Norm() <= R::MaxL1Norm() && target.second <= 1) {
      S.push_back(target.first);
    } else {
      nontrivial_targets.push_back(target);
    }
  }

  return nontrivial_targets;
}

template<class R>
void RemoveReplacementColumns(std::vector<Vector<R>>& S) {
  for (auto& c : S) {
    for (int i = 0; i < R::Dim(); ++i) if (c[i] == -1) {
      for (auto& c_ : S) if (&c_ != &c) {
        if (c_[i] >= 1) {
          c_ += c;
          c.Reset(); // c = zero column
          goto checkNextColumn;
        }
      }
    }
    checkNextColumn:;
  }

  // remove zero columns
  S.erase(
    std::remove_if(
      S.begin(),
      S.end(),
      [](auto const & c) { return c.IsZero(); }
    ),
    S.end()
  );
}

template<class R>
inline void PreviousAnchor(Vector<R>& anchor) {
  for (int i = 0; i < R::Dim(); ++i) {
    anchor[i] = R::PreviousInterval(i, anchor[i]);
  }
}

template<class R, class C>
void Scheduler<R, C>::ComputeSchedule(const Vector<R>& b, int m, std::vector<Vector<R>>& S) {

  if (b.IsZero()) return;
  if (m == 1) S.push_back(b);
  if (m <= 1) return;

  Vector<R> anchor = b;
  std::vector<std::pair<Vector<R>,int>> targets;
  targets.push_back(std::pair<Vector<R>,int>(b, m));

  Vector<R> anchor_;
  anchor_.Reset();

  Tensor<R> t, t_;

  while (!targets.empty()) {

    if (anchor != anchor_) {

      PreviousAnchor(anchor);

      std::cout << "anchor = " << anchor << std::endl;

      t.Reset();
      t_.Reset();
      MinMachines(t, t_, anchor);
      anchor_ = anchor;

      targets = Backtrace(anchor, t, targets, S);
      if (targets.empty()) break;

      PreviousAnchor(anchor);

      std::cout << "anchor = " << anchor << std::endl;
    }

    targets = Backtrace(anchor, t_, targets, S);
  }

  for (auto& c : S) std::cout << c << ", "; std::cout << std::endl;

  RemoveReplacementColumns(S);
}

#define INSTANTIATE_SCHEDULER(R) \
  template class Scheduler<R, NaiveConvolution<R>>; \
  template class Scheduler<R, ParallelNaiveConvolution<R>>; \
  template class Scheduler<R, FFTConvolution<R>>; \
  template class Scheduler<R, DoubleCheckConvolution<R, FFTConvolution<R>, ParallelNaiveConvolution<R>>>;
ROUNDINGS_LIST(INSTANTIATE_SCHEDULER)
