#include <iostream>
#include "tensor.h"

#include "rounding.h"

template<class R>
void Tensor<R>::Initialize(Vector<R>& anchor) {
  Vector<R> s;
  Reset();

  // zero vector
  s.Reset();
  s -= anchor;
  Set(s, 0);

  // non-negative vectors with bounded l1 norm
  for (std::size_t norm=1; norm<=R::MaxL1Norm(); ++norm) {
    std::size_t num_combinations = 1;
    for (std::size_t i=0; i<norm; ++i) {
      num_combinations *= R::Dim();
    }
    for (std::size_t code=0; code<num_combinations; ++code) {
      s.Reset();
      std::size_t remainder = code;
      for (std::size_t i=0; i<norm; ++i) {
        int d = remainder % R::Dim();
        remainder = remainder / R::Dim();
        s[d] = s[d] + 1;
      }
      // This enumerates all non-negative s with bounded l1 norm
      if (s.Volume() <= R::Makespan() - R::Precision()) {
        s -= anchor;
        Set(s, 1);
      }
    }
  }

  // substitutions
  for (int i=0; i<R::Dim(); ++i) {
    for (int j=0; j<R::Dim(); ++j) {
      for (int k=0; k<R::Dim(); ++k) {
        // Substitude i for j and k
        s.Reset();
        s[i] = -1;
        s[j] = 1;
        s[k] = s[k] + 1;
        if (std::abs(s.Volume()) <= R::Precision()) {
          s -= anchor;
          Set(s, 0);
        }
      }
    }
  }
}

#define INSTANTIATE_TENSOR(R) template class Tensor<R>;
ROUNDINGS_LIST(INSTANTIATE_TENSOR)
