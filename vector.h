#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>
#include <memory>
#include <cstring>

template<class R>
class Vector {
 public:
  Vector() {
#ifndef NDEBUG
    for (int i=0; i<R::Dim(); ++i) {
      data_[i] = 0xDEADBEEF;
    }
#endif // NDEBUG
  }

  inline int& operator[](int i) {
    return data_[i];
  }

  inline int operator[](int i) const {
    return data_[i];
  }

  inline bool operator==(const Vector<R>& rhs) const {
    for (int i=0; i<R::Dim(); ++i) {
      if ((*this)[i] != rhs[i]) return false;
    }
    return true;
  }

  inline bool operator!=(const Vector<R>& rhs) const {
    for (int i=0; i<R::Dim(); ++i) {
      if ((*this)[i] != rhs[i]) return true;
    }
    return false;
  }

  inline bool IsWithinDeviation() const {
    for (int i=0; i<R::Dim(); ++i) {
      if (0 > (*this)[i] || (*this)[i] >= R::Deviation(i)) return false;
    }
    return true;
  }

  inline int L1Norm() const {
    int norm = 0;
    for (int i=0; i<R::Dim(); ++i) {
      norm += std::abs((*this)[i]);
    }
    return norm;
  }

  inline Vector<R>& operator+=(const Vector<R>& rhs) {
    for (int i=0; i<R::Dim(); ++i) {
      (*this)[i] += rhs[i];
    }
    return *this;
  }

  inline Vector<R>& operator-=(const Vector<R>& rhs) {
    for (int i=0; i<R::Dim(); ++i) {
      (*this)[i] -= rhs[i];
    }
    return *this;
  }

  inline double Volume() const {
    double v = 0.0;
    for (int i=0; i<R::Dim(); ++i) {
      v += R::Size(i) * (*this)[i];
    }
    return v;
  }

  inline void Reset() {
    std::memset(data_, 0, sizeof(data_));
  }

  inline bool IsPositive() const {
    for (int i = 0; i < R::Dim(); ++i) {
      if ((*this)[i] < 0) return false;
    }
    return true;
  }

  inline bool IsZero() const {
    for (int i = 0; i < R::Dim(); ++i) {
      if ((*this)[i] != 0) return false;
    }
    return true;
  }

 private:
  int data_[R::Dim()];
};

template<class R>
std::ostream& operator<<(std::ostream& os, const Vector<R>& obj) {
  os << "(";
  for (int i=0; i<R::Dim(); ++i) {
    os << obj[i];
    if (i < R::Dim() - 1) os << ", ";
  }
  os << ")";
  return os;
}

#endif // VECTOR_H_
