    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/tensors/Tensor_exp.h

    Copyright (C) 2015

Author: neo <cossu@post.kek.jp>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_MATH_EXP_H
#define GRID_MATH_EXP_H

#define DEFAULT_MAT_EXP 12

namespace Grid {

  ///////////////////////////////////////////////
  // Exponentiate function for scalar, vector, matrix
  ///////////////////////////////////////////////


  template<class vtype> inline iScalar<vtype> Exponentiate(const iScalar<vtype>&r, RealD alpha ,  Integer Nexp = DEFAULT_MAT_EXP)
    {
      iScalar<vtype> ret;
      ret._internal = Exponentiate(r._internal, alpha, Nexp);
      return ret;
    }

template<class vtype, int N> inline iVector<vtype, N> Exponentiate(const iVector<vtype,N>&r, RealD alpha ,  Integer Nexp = DEFAULT_MAT_EXP)
    {
      iVector<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Exponentiate(r._internal[i], alpha, Nexp);
      return ret;
    }



    // Specialisation: Cayley-Hamilton exponential for SU(3)
    template<class vtype, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0>::type * =nullptr>
    inline iMatrix<vtype,3> Exponentiate(const iMatrix<vtype,3> &arg, RealD alpha  , Integer Nexp = DEFAULT_MAT_EXP )
    {
    // for SU(3) 2x faster than the std implementation using Nexp=12
    // notice that it actually computes
    // exp ( input matrix )
    // the i sign is coming from outside
    // input matrix is anti-hermitian NOT hermitian
      typedef iMatrix<vtype,3> mat;
      typedef iScalar<vtype> scalar;
      mat unit(1.0);
      mat temp(unit);
      const Complex one_over_three = 1.0 / 3.0;
      const Complex one_over_two = 1.0 / 2.0;

      scalar c0, c1, tmp, c0max, theta, u, w;
      scalar xi0, u2, w2, cosw;
      scalar fden, h0, h1, h2;
      scalar e2iu, emiu, ixi0, qt;
      scalar f0, f1, f2;
      scalar unity(1.0);

      mat iQ2 = arg*arg*alpha*alpha;
      mat iQ3 = arg*iQ2*alpha;
      // sign in c0 from the conventions on the Ta
      scalar imQ3, reQ2;
      imQ3 = imag( trace(iQ3) );
      reQ2 = real( trace(iQ2) );
      c0 = -imQ3 * one_over_three;
      c1 = -reQ2 * one_over_two;

      // Cayley Hamilton checks to machine precision, tested
      tmp = c1 * one_over_three;
      c0max = 2.0 * pow(tmp, 1.5);

      theta = acos(c0 / c0max) * one_over_three;
      u = sqrt(tmp) * cos(theta);
      w = sqrt(c1) * sin(theta);

      xi0 = sin(w) / w;
      u2 = u * u;
      w2 = w * w;
      cosw = cos(w);

      ixi0 = timesI(xi0);
      emiu = cos(u) - timesI(sin(u));
      e2iu = cos(2.0 * u) + timesI(sin(2.0 * u));

      h0 = e2iu * (u2 - w2) +
           emiu * ((8.0 * u2 * cosw) + (2.0 * u * (3.0 * u2 + w2) * ixi0));
      h1 = e2iu * (2.0 * u) - emiu * ((2.0 * u * cosw) - (3.0 * u2 - w2) * ixi0);
      h2 = e2iu - emiu * (cosw + (3.0 * u) * ixi0);

      fden = unity / (9.0 * u2 - w2);  // reals
      f0 = h0 * fden;
      f1 = h1 * fden;
      f2 = h2 * fden;

      return (f0 * unit + timesMinusI(f1) * arg*alpha - f2 * iQ2);
    }



// General exponential
template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr>
    inline iMatrix<vtype,N> Exponentiate(const iMatrix<vtype,N> &arg, RealD alpha  , Integer Nexp = DEFAULT_MAT_EXP )
    {
    // notice that it actually computes
    // exp ( input matrix )
    // the i sign is coming from outside
    // input matrix is anti-hermitian NOT hermitian
      typedef iMatrix<vtype,N> mat;
      mat unit(1.0);
      mat temp(unit);
      for(int i=Nexp; i>=1;--i){
	      temp *= alpha/RealD(i);
	      temp = unit + temp*arg;
      }
      return temp;

    }

// Exponential of a series
template<class vtype,int N>
    inline iSeries<vtype,N> Exponentiate(const iSeries<vtype,N> &arg, RealD alpha  , Integer Nexp = DEFAULT_MAT_EXP )
    {
      // This way the result will be exact if the constant term of arg is zero (which is the usual case in NSPT)
      // Actually, there is a factor ~2 optimization opportunity here in that case, by avoiding multiplication by zero
      if(Nexp < N-1)
        Nexp = N-1;

      iSeries<vtype,N> unit(1.0);
      iSeries<vtype,N> temp(unit);
      for(int i=Nexp; i>=1;--i){
	      temp *= iScalar<vtype>(alpha/RealD(i));
	      temp = unit + temp*arg;
      }
      return temp;

    }


    template<class vtype> inline iScalar<vtype> Logarithm(const iScalar<vtype>&r)
      {
        iScalar<vtype> ret;
        ret._internal = Logarithm(r._internal);
        return ret;
      }

    template<class vtype, int N> inline iVector<vtype, N> Logarithm(const iVector<vtype,N>&r)
      {
        iVector<vtype, N> ret;
        for (int i = 0; i < N; i++)
          ret._internal[i] = Logarithm(r._internal[i]);
        return ret;
      }

      // NOTE: basic Taylor series which is only valid if arg is close to 1
template<class vtype, int N> inline iSeries<vtype,N> Logarithm(const iSeries<vtype,N> & arg)
{
    iSeries<vtype,N> x = arg;
    x(0) -= vtype(1.0);
    iSeries<vtype,N> temp = x;
    iSeries<vtype,N> xn = x;
    for(int i = 2; i < N; ++i)
    {
        xn = x*xn;
        temp += xn*iScalar<vtype>(1.0/i * (i%2?1:-1));
    }
    return temp;
}

///////////////////////////////////////////////
// Exponentiate "Fast" function for scalar, vector, matrix.
// This version assumes (without check!) the constant term in the series expansion to be trivial
///////////////////////////////////////////////


template<class vtype> inline iScalar<vtype> ExponentiateFast(const iScalar<vtype>&r, RealD alpha)
  {
    iScalar<vtype> ret;
    ret._internal = ExponentiateFast(r._internal, alpha);
    return ret;
  }

template<class vtype, int N> inline iVector<vtype, N> ExponentiateFast(const iVector<vtype,N>&r, RealD alpha)
  {
    iVector<vtype, N> ret;
    for (int i = 0; i < N; i++)
      ret._internal[i] = ExponentiateFast(r._internal[i], alpha);
    return ret;
  }

// NOTE: no overload for iMatrix

/* computes "a = a*b" assuming b(0)=0 and a(0,1,...k-1) = 0 */
template<class vtype, int N>
strong_inline void mulAssignFast(iSeries<vtype,N>& a, const iSeries<vtype,N>& b, int k)
{
    for(int i = N-1; i >= k; --i)
    {
        a._internal[i] = 0;
        for(int j = 1; j <= i-k; ++j)
            a._internal[i] += a._internal[i-j]*b._internal[j];
    }
}

// Exponential of a series
template<class vtype,int N>
inline iSeries<vtype,N> ExponentiateFast(const iSeries<vtype,N> &arg, RealD alpha)
{
    // Assuming arg._internal[0] == 0 saves half the multiplications (for large N)
    // In reality, it seems to save about 35%

    iSeries<vtype,N> xn = arg*alpha;
    iSeries<vtype,N> temp = xn;
    temp._internal[0] = 1.0;

    for(int i = 2; i < N; ++i)
    {
        for(int l = N-1; l >= i; --l)
        {
            mult(&xn._internal[l], &xn._internal[l-1], &arg._internal[1]);
            for(int j = 2; j <= l-i+1; ++j)
                mac(&xn._internal[l], &xn._internal[l-j], &arg._internal[j]);

            xn._internal[l] *= alpha/i;
            temp._internal[l] += xn._internal[l];
        }
    }
    return temp;
}

template<class vtype> inline iScalar<vtype> LogarithmFast(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = LogarithmFast(r._internal);
    return ret;
}

template<class vtype, int N> inline iVector<vtype, N> LogarithmFast(const iVector<vtype,N>&r)
{
    iVector<vtype, N> ret;
    for (int i = 0; i < N; i++)
        ret._internal[i] = LogarithmFast(r._internal[i]);
    return ret;
}

template<class vtype, int N> inline iSeries<vtype,N> LogarithmFast(const iSeries<vtype,N> & arg)
{
    // Assume arg[0] = 1 and use basic Taylor series log(1+x)=x-x^2/2+...
    iSeries<vtype,N> xn = arg;
    xn._internal[0] = 0.0;
    iSeries<vtype,N> temp = xn;

    for(int i = 2; i < N; ++i)
    {
        for(int l = N-1; l >= i; --l)
        {
            mult(&xn._internal[l], &xn._internal[l-1], &arg._internal[1]);
            for(int j = 2; j <= l-i+1; ++j)
                mac(&xn._internal[l], &xn._internal[l-j], &arg._internal[j]);

            temp._internal[l] += xn._internal[l]*RealD(1.0/i * (i%2?1:-1));
        }
    }
    return temp;
}





}
#endif
