#include <complex>
#include <fftw3.h>
#include <omp.h>

#include "convolution.h"

#include "rounding.h"

template<class R>
void NaiveConvolution<R>::Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& source_anchor) {
  Vector<R> a;
  Vector<R> b;
  target.Reset();

  for (std::size_t i = 0; i < Tensor<R>::data_size(); ++i) {
    Tensor<R>::FromIndex(a, i);
    int va = source.Get(a);
    if (va == -1) continue;

    for (std::size_t j = i; j < Tensor<R>::data_size(); ++j) {
      Tensor<R>::FromIndex(b, j);
      int vb = source.Get(b);
      if (vb == -1) continue;

      b += source_anchor;
      b += source_anchor;
      b += a;
      b -= target_anchor;

      if (b.IsWithinDeviation()) {
        int vc = target.Get(b);
        if (vc == -1 || vc > va + vb) {
          target.Set(b, va + vb);
        }
      }
    }
  }
}

template<class R>
void ParallelNaiveConvolution<R>::Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& source_anchor) {
  std::cout << "parallel convolution " << Tensor<R>::data_size() << std::endl;
  #pragma omp parallel
  {
    Vector<R> a;
    Vector<R> a2;
    a2.Reset();
    Vector<R> b;
    Vector<R> c;

    #pragma omp for
    for (size_t i = 0; i < Tensor<R>::data_size(); ++i) {
      Tensor<R>::FromIndex(a, i);

      // This optimization cuts the search space down slightly
      a2[R::Dim() - 1] = (a[R::Dim() - 1] + target_anchor[R::Dim() - 1]) / 2 + 1 - source_anchor[R::Dim() - 1];
      size_t first, last;
      if (2 * a2[R::Dim() - 1] < R::Deviation(R::Dim() - 1)) {
       first = 0;
       last = Tensor<R>::ToIndex(a2);
      } else {
       a2[R::Dim() - 1]--;
       first = Tensor<R>::ToIndex(a2);
       last = Tensor<R>::data_size();
      }
      int va = -1;

      for (size_t j = first; j < last; ++j) {
        Tensor<R>::FromIndex(b, j);
        int vb = source.Get(b);
        if (vb == -1) continue;

        c = a;
        c += target_anchor;
        c -= b;
        c -= source_anchor;
        c -= source_anchor;

        if (c.IsWithinDeviation()) {
          int vc = source.Get(c);
          if (vc == -1) continue;
          if (va == -1 || va > vb + vc) {
            va = vb + vc;
          }
        }
      }
      target.Set(a, va);
    }
  }
}

template<class R>
void FFTConvolution<R>::Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& source_anchor) {
  const std::size_t n = Tensor<CheapPaddedRounding<R>>::data_size();
  std::atomic_int m_max, m_min;
  m_min = 9999999;
  m_max = -1;
  #pragma omp parallel for
  for (std::size_t i = 0; i < Tensor<R>::data_size(); ++i) {
    Vector<R> c;
    Tensor<R>::FromIndex(c, i);
    int v = source.Get(c);
    int cur = m_max;
    while (v != -1 && v > cur) {
      if (std::atomic_compare_exchange_weak(&m_max, &cur, v)) {
        break;
      }
    }
    while (v != -1 && v < cur) {
      if (std::atomic_compare_exchange_weak(&m_min, &cur, v)) {
        break;
      }
    }
  }
  target.Reset();

  int combinations = (m_max - m_min + 1) * (m_max - m_min + 2) / 2;

  int deviations[R::Dim()];
  for (std::size_t i = 0; i < R::Dim(); ++i) {
    deviations[i] = CheapPaddedRounding<R>::Deviation(i);
  }

  fftwf_plan plan_forward = fftwf_plan_dft(R::Dim(), deviations, 0, 0, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_plan plan_backward = fftwf_plan_dft(R::Dim(), deviations, 0, 0, FFTW_BACKWARD, FFTW_ESTIMATE);

  #pragma omp parallel num_threads(15)
  {
    Vector<R> a, b;
    a.Reset();
    std::complex<float>* fftw_in = (std::complex<float>*)fftwf_malloc(n * sizeof(std::complex<float>));
    std::complex<float>* fftw_in_ = (std::complex<float>*)fftwf_malloc(n * sizeof(std::complex<float>));
    #pragma omp for
    for (int iter = 0; iter < combinations; ++iter) {
      int remainder = iter;
      int m = m_min;
      int m_ = m_min;
      while (remainder-- > 0) {
        if (m <= m_) {
          m++;
          m_ = m_min;
        }
        else {
          m_++;
        }
      }
      
      std::fill(fftw_in, fftw_in+n, std::complex<float>{});
      std::fill(fftw_in_, fftw_in_+n, std::complex<float>{});

      for (std::size_t i = 0; i < Tensor<R>::data_size(); ++i) {
        Tensor<R>::FromIndex(a, i);
        int v = source.Get(a);
        int idx = Tensor<CheapPaddedRounding<R>>::ToIndex(reinterpret_cast<Vector<CheapPaddedRounding<R>>&>(a));
        if (m == v) {
          fftw_in[idx].real(1.0f / (float)n);
        }
        if (m_ == v) {
          fftw_in_[idx].real(1.0f);
        }
      }

      fftwf_execute_dft(plan_forward, (fftwf_complex*)fftw_in, (fftwf_complex*)fftw_in);
      fftwf_execute_dft(plan_forward, (fftwf_complex*)fftw_in_, (fftwf_complex*)fftw_in_);

      for (std::size_t i = 0; i < Tensor<CheapPaddedRounding<R>>::data_size(); ++i) {
        fftw_in[i] *= fftw_in_[i];
      }

      fftwf_execute_dft(plan_backward, (fftwf_complex*)fftw_in, (fftwf_complex*)fftw_in);

      for (std::size_t i = 0; i < Tensor<R>::data_size(); ++i) {
        Tensor<R>::FromIndex(a, i);
        a += target_anchor;
        a -= source_anchor;
        a -= source_anchor;
        int idx = Tensor<CheapPaddedRounding<R>>::ToIndex(reinterpret_cast<Vector<CheapPaddedRounding<R>>&>(a));
        int v = target.data_unsafe()[i];
        if (fftw_in[idx].real() > 0.5) {
          while (v == -1 || v > m + m_) {
            if (std::atomic_compare_exchange_weak(&target.data_unsafe()[i], &v, m + m_)) {
              break;
            }
          }
        }
      }
    }
    fftwf_free(fftw_in);
    fftwf_free(fftw_in_);
  }
  fftwf_destroy_plan(plan_forward);
  fftwf_destroy_plan(plan_backward);
}

template<class R, class C, class D>
void DoubleCheckConvolution<R, C, D>::Square(Tensor<R>& target, const Tensor<R>& source, const Vector<R>& target_anchor, const Vector<R>& source_anchor) {
  std::hash<std::string> h;
  C::Square(target, source, target_anchor, source_anchor);
  std::string s((char*)target.data_unsafe(), target.data_size() * sizeof(int));

  D::Square(target, source, target_anchor, source_anchor);
  std::string s_((char*)target.data_unsafe(), target.data_size() * sizeof(int));
  std::cout << "Hashes: " << h(s) << ", " << h(s_) << std::endl;
  if (h(s) != h(s_)) {
    std::cerr << "Different hashes" << std::endl;
    exit(1);
  }
}

#define INSTANTIATE_CONVOLUTION(R) \
  template class NaiveConvolution<R>; \
  template class ParallelNaiveConvolution<R>; \
  template class DoubleCheckConvolution<R, FFTConvolution<R>, ParallelNaiveConvolution<R>>; \
  template class FFTConvolution<R>;
ROUNDINGS_LIST(INSTANTIATE_CONVOLUTION)
