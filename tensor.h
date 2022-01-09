#ifndef TENSOR_H_
#define TENSOR_H_

#include <memory>
#include <cstring>
#include <cassert>
#include <atomic>

#include "vector.h"

template<class R>
class Tensor {
 public:

  Tensor() :
    data_(new std::atomic_int[data_size()])
  {
  }

  inline void Reset(int v = -1) {
    for (std::size_t i=0; i<data_size(); ++i) {
      data_[i] = v;
    }
  }

  inline void Set(const Vector<R>& s, int v) {
    data_[ToIndex(s)] = v;
  }
  
  inline int Get(const Vector<R>& s) const {
    return data_[ToIndex(s)];
  }

  void Initialize(Vector<R>& anchor);

  static constexpr std::size_t data_size(int i=R::Dim()) {
    std::size_t size = 1;
    while (i > 0) size *= R::Deviation(i--);
    return size;
  }

  std::atomic_int* data_unsafe() {
    return data_.get();
  }

  static inline std::size_t ToIndex(const Vector<R>& s) {
    assert(s.IsWithinDeviation());
    std::size_t idx = 0;
    for (int i=R::Dim() - 1; i>=0; --i) {
      idx = idx * R::Deviation(i) + s[i];
    }
    assert(idx >= 0);
    assert(idx < data_size());
    return idx;
  }

  static inline void FromIndex(Vector<R>& s, std::size_t idx) {
    assert(idx >= 0);
    assert(idx < data_size());
    for (int i=0; i<R::Dim(); ++i) {
      s[i] = idx % R::Deviation(i);
      idx = idx / R::Deviation(i);
    }
    assert(s.IsWithinDeviation());
  }

  inline bool operator==(const Tensor<R>& rhs) const {
    return !memcmp(data_.get(), rhs.data_.get(), data_size() * sizeof(int));
  }

 private:
  std::unique_ptr<std::atomic_int[]> data_;
};

#endif // TENSOR_H_
