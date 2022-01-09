#ifndef ROUNDING_HH_
#define ROUNDING_HH_

#include <cmath>

class Rounding {
 public:

  // TODO: Does this always apply or just for norm 3?
  static constexpr int PreviousInterval(int i, int begin) {
    // if begin is even, then begin / 2 - 1
    // if begin is odd, then (begin-1) / 2 - 1
    return begin % 2 ? (begin - 1) / 2 - 1 : begin / 2 - 1;
  }

};

template<int unused>
class Rounding9 : public Rounding {
 public:
  static constexpr int Dim() { return 9; }

  static constexpr int Deviation(int i) { return 6; }

  static constexpr double Size(int i) {
    return sizes_[i];
  }

  static constexpr int MaxL1Norm() { return 3; }

  static constexpr double Makespan() { return 1.0; }

  static constexpr double Precision() { return 0.00001; }

  static constexpr int GetRoundedIndex(double p) {
    for (int i = 0; i < Dim(); ++i) {
      if (p >= sizes_[i]) return i;
    }
    return Dim();// too small
  }

 private:
  //Best bound for n=9 is eps=0.1754019165039063
  static constexpr double sizes_[] = {
    0.556953353881836,
    0.48461708384337826,
    0.4122990417480468,
    0.3508038330078126,
    0.33470772331794213,
    0.2848099245600205,
    0.24230854192168913,
    0.2061495208740234,
    0.1754019165039063
  };

};

template<int unused>
class Rounding10 : public Rounding {
 public:
  static constexpr int Dim() { return 10; }

  static constexpr int Deviation(int i) { return 6; }

  static constexpr int Size(int i) {
    // SIZES[] = {12, 12 + 2, 12 + 4, 12 + 6, 12 + 9, 24, 24 + 4, 24 + 8, 24 + 12, 24 + 18};
    return i <= 3 ? 12 + 2 * i : (i == 4 ? 12 + 9 : (i <= 8 ? 24 + 4 * i : 24 + 18));
  }

  static constexpr int MaxL1Norm() { return 3; }

  static constexpr int Makespan() { return 72; }

  static constexpr int Precision() { return 0; }

};

template<int K>
class ArithmeticRounding : public Rounding {
 public:
  static constexpr int Dim() { return K; }

  static constexpr int Deviation(int i) { return 6; }

  static constexpr int Size(int i) { return i+1; }

  static constexpr int MaxL1Norm() { return 3; }

  static constexpr int Makespan() { return K+1; }

  static constexpr int Precision() { return 0; }

};

template<class R>
class PaddedRounding : public R {
 public:
  static constexpr int Deviation(int i) {
    return 2 * R::Deviation(i) - 1;
  }

};

template<class R>
class CheapPaddedRounding : public R {
 public:
  static constexpr int Deviation(int i) {
    return R::Deviation(i) * 3 / 2;
  }

};


#define ROUNDINGS_LIST_HELPER(X, R) \
  X(R) \
  X(PaddedRounding<R>) \
  X(CheapPaddedRounding<R>)
#define ROUNDINGS_LIST(X) \
  ROUNDINGS_LIST_HELPER(X, Rounding9<0>) \
  ROUNDINGS_LIST_HELPER(X, Rounding10<0>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<1>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<2>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<3>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<4>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<5>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<6>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<7>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<8>) \
  ROUNDINGS_LIST_HELPER(X, ArithmeticRounding<9>)

#endif
