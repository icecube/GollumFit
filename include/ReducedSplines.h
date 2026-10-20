#ifndef GOLLUMFIT_REDUCEDSPLINES_H_
#define GOLLUMFIT_REDUCEDSPLINES_H_

#include <cstdint>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>

/**
* @file ReducedSplines.h
* @brief A 3-D photospline table reduced to a 1-D B-spline along its last axis.
*
* The detector-systematics splines are evaluated at (x0, x1, p) where x0 and x1
* are fixed per event (log10 energy, cos zenith) and only p (the nuisance
* parameter) changes between likelihood evaluations. Because the table is a
* tensor-product B-spline,
*
*   value(x0,x1,p) = sum_{i,j,k} c[i,j,k] B_i(x0) B_j(x1) B_k(p)
*                  = sum_k d[k] B_k(p),   d[k] = sum_{i,j} c[i,j,k] B_i(x0) B_j(x1)
*
* so each event only needs its n_k reduced coefficients d[k], and the basis
* functions B_k(p) (order+1 of them are non-zero) are shared by every event in
* one evaluation. The arithmetic is the same as photospline's ndsplineeval,
* only reordered.
*/

namespace gollumfit {

class ReducedAxisSpline {
public:
    static constexpr unsigned int kAxis = 2; ///< axis along which the 1-D spline is kept
    static constexpr unsigned int kMaxOrder = 7;

    /// Basis functions along the kept axis at one parameter value.
    struct Basis {
        bool valid = false;      ///< p inside the knot range (photospline's searchcenters succeeded)
        int first = 0;           ///< index of the first non-zero reduced coefficient
        unsigned int n = 0;      ///< number of non-zero basis functions (order+1)
        float b[kMaxOrder + 1];  ///< B_k(p)
        float db[kMaxOrder + 1]; ///< dB_k/dp
    };

    explicit ReducedAxisSpline(const photospline::splinetable<>& table) : table_(&table) {
        if (table.get_ndim() != 3)
            throw std::runtime_error("ReducedAxisSpline: only 3-dimensional tables are supported");
        for (unsigned int d = 0; d < 3; ++d) {
            order_[d] = table.get_order(d);
            nknots_[d] = table.get_nknots(d);
            naxes_[d] = table.get_ncoeffs(d);
            strides_[d] = table.get_stride(d);
            knots_[d] = table.get_knots(d);
            if (order_[d] > kMaxOrder)
                throw std::runtime_error("ReducedAxisSpline: spline order too large");
        }
        coeffs_ = table.get_coefficients();
    }

    /// Number of reduced coefficients kept per event.
    unsigned int ncoeffs() const { return naxes_[kAxis]; }

    /// Reduced coefficients for a fixed point (x0,x1) of the first two axes.
    /// Returns false, with out[] zero-filled, if the point is outside the
    /// table support (photospline's operator() evaluates to 0 there).
    bool reduce(double x0, double x1, float* out) const {
        const unsigned int nk = naxes_[kAxis];
        std::fill(out, out + nk, 0.f);
        int c0, c1;
        if (!searchcenter(0, x0, c0) || !searchcenter(1, x1, c1))
            return false;
        float b0[kMaxOrder + 1], b1[kMaxOrder + 1];
        photospline::bsplvb_simple(knots_[0], nknots_[0], x0, c0, order_[0] + 1, b0);
        photospline::bsplvb_simple(knots_[1], nknots_[1], x1, c1, order_[1] + 1, b1);
        double acc[256];
        if (nk > 256)
            throw std::runtime_error("ReducedAxisSpline: too many coefficients along the reduced axis");
        std::fill(acc, acc + nk, 0.);
        for (unsigned int i = 0; i <= order_[0]; ++i) {
            for (unsigned int j = 0; j <= order_[1]; ++j) {
                const double w = static_cast<double>(b0[i]) * static_cast<double>(b1[j]);
                const int64_t base = (static_cast<int64_t>(c0) - order_[0] + i) * static_cast<int64_t>(strides_[0])
                                   + (static_cast<int64_t>(c1) - order_[1] + j) * static_cast<int64_t>(strides_[1]);
                const float* col = coeffs_ + base;
                const int64_t s2 = strides_[2];
                for (unsigned int k = 0; k < nk; ++k)
                    acc[k] += w * col[k * s2];
            }
        }
        for (unsigned int k = 0; k < nk; ++k)
            out[k] = static_cast<float>(acc[k]);
        return true;
    }

    /// Basis functions (and derivatives) along the kept axis at parameter value p.
    Basis basis(double p) const {
        Basis bs;
        int c;
        if (!searchcenter(kAxis, p, c))
            return bs;
        bs.valid = true;
        bs.n = order_[kAxis] + 1;
        bs.first = c - static_cast<int>(order_[kAxis]);
        photospline::bsplvb_simple(knots_[kAxis], nknots_[kAxis], p, c, order_[kAxis] + 1, bs.b);
        photospline::bspline_deriv_nonzero(knots_[kAxis], nknots_[kAxis], p, c, order_[kAxis], bs.db);
        return bs;
    }

    /// Spline value from reduced coefficients d and a basis.
    static inline float eval(const float* d, const Basis& bs) {
        float r = 0.f;
        for (unsigned int j = 0; j < bs.n; ++j)
            r += d[bs.first + j] * bs.b[j];
        return r;
    }
    /// d(value)/dp from reduced coefficients d and a basis.
    static inline float evalDeriv(const float* d, const Basis& bs) {
        float r = 0.f;
        for (unsigned int j = 0; j < bs.n; ++j)
            r += d[bs.first + j] * bs.db[j];
        return r;
    }

    const photospline::splinetable<>& table() const { return *table_; }

private:
    const photospline::splinetable<>* table_;
    unsigned int order_[3], nknots_[3], naxes_[3];
    uint64_t strides_[3];
    const double* knots_[3];
    const float* coeffs_;

    /// photospline::splinetable::searchcenters for a single axis.
    bool searchcenter(unsigned int dim, double x, int& center) const {
        const double* knots = knots_[dim];
        const unsigned int nknots = nknots_[dim];
        const int order = order_[dim];
        const int naxes = naxes_[dim];
        if (x <= knots[0] || x > knots[nknots - 1])
            return false;
        if (x < knots[order]) {
            center = order;
            return true;
        } else if (x >= knots[naxes]) {
            center = naxes - 1;
            return true;
        }
        unsigned int lo = order;
        unsigned int hi = nknots - 2;
        do {
            center = (hi + lo) / 2;
            if (x < knots[center])
                hi = center - 1;
            else
                lo = center + 1;
        } while (x < knots[center] || x >= knots[center + 1]);
        if (center == naxes)
            center--;
        return true;
    }
};

} // namespace gollumfit

#endif /* GOLLUMFIT_REDUCEDSPLINES_H_ */
