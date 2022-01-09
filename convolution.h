#ifndef CONVOLUTION_H_
#define CONVOLUTION_H_

#include "tensor.h"
#include "vector.h"

template<class R>
class NaiveConvolution {
 public:
  static void Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& target_source);

};

template<class R>
class ParallelNaiveConvolution {
 public:
  static void Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& target_source);

};

template<class R>
class FFTConvolution {
 public:
  static void Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& source_anchor);
};

template<class R, class C, class D>
class DoubleCheckConvolution {
 public:
  static void Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& source_anchor);
};


#endif // CONVOLUTION_H_
