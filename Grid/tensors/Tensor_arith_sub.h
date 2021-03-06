    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/tensors/Tensor_arith_sub.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_MATH_ARITH_SUB_H
#define GRID_MATH_ARITH_SUB_H

namespace Grid {


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// SUB         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////


// SUB is simple for now; cannot mix types and straightforward template
// Scalar +/- Scalar
// Vector +/- Vector
// Matrix +/- Matrix
// Matrix /- scalar
template<class vtype,class ltype,class rtype> strong_inline void sub(iScalar<vtype> * __restrict__ ret,
                                                              const iScalar<ltype> * __restrict__ lhs,
                                                              const iScalar<rtype> * __restrict__ rhs)
{
    sub(&ret->_internal,&lhs->_internal,&rhs->_internal);
}

template<class vtype,class ltype,class rtype,int N> strong_inline void sub(iVector<vtype,N> * __restrict__ ret,
                                                                    const iVector<ltype,N> * __restrict__ lhs,
                                                                    const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c=0;c<N;c++){
        ret->_internal[c]=lhs->_internal[c]-rhs->_internal[c];
    }
    return;
}
template<class vtype,class ltype,class rtype, int N> strong_inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iMatrix<ltype,N> * __restrict__ lhs,
                                                                     const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        sub(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> strong_inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iScalar<ltype> * __restrict__ lhs,
                                                                     const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1==c2) {
            sub(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
        } else {
            // Fails -- need unary minus. Catalogue other unops?
            ret->_internal[c1][c2]=zero;
            ret->_internal[c1][c2]=ret->_internal[c1][c2]-rhs->_internal[c1][c2];

        }
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> strong_inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iMatrix<ltype,N> * __restrict__ lhs,
                                                                     const iScalar<rtype> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1==c2)
            sub(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
        else
            ret->_internal[c1][c2]=lhs->_internal[c1][c2];
    }}
    return;
}
template<class vtype,class ltype,class rtype,int N> strong_inline void sub(iSeries<vtype,N> * __restrict__ ret,
                                                                    const iSeries<ltype,N> * __restrict__ lhs,
                                                                    const iSeries<rtype,N> * __restrict__ rhs)
{
    for(int c=0;c<N;c++){
        ret->_internal[c]=lhs->_internal[c]-rhs->_internal[c];
    }
    return;
}

template<class vtype,class ltype,class rtype,int N> strong_inline void sub(iSeries<vtype,N> * __restrict__ ret,
                                                                    const iScalar<ltype> * __restrict__ lhs,
                                                                    const iSeries<rtype,N> * __restrict__ rhs)
{
    ret->_internal[0] = lhs->internal - rhs->_internal[0];
    for(int c=1;c<N;c++){
        ret->_internal[c]= -rhs->_internal[c];
    }
    return;
}
template<class vtype,class ltype,class rtype,int N> strong_inline void sub(iSeries<vtype,N> * __restrict__ ret,
                                                                    const iSeries<ltype,N> * __restrict__ lhs,
                                                                    const iScalar<rtype> * __restrict__ rhs)
{
    ret->_internal[0] = lhs->internal[0] - rhs->_internal;
    for(int c=1;c<N;c++){
        ret->_internal[c]= lhs->_internal[c];
    }
    return;
}



    // - operator for scalar, vector, matrix
template<class ltype,class rtype> strong_inline auto
operator - (const iScalar<ltype>& lhs, const iScalar<rtype>& rhs) -> iScalar<decltype(lhs._internal - rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal-rhs._internal)> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iVector<ltype,N>& lhs,const iVector<rtype,N>& rhs) ->iVector<decltype(lhs._internal[0]-rhs._internal[0]),N>
{
    typedef iVector<decltype(lhs._internal[0]-rhs._internal[0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iMatrix<ltype,N>& lhs,const iMatrix<rtype,N>& rhs) ->iMatrix<decltype(lhs._internal[0][0]-rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]-rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iScalar<ltype>& lhs,const iMatrix<rtype,N>& rhs)->iMatrix<decltype(lhs._internal-rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal-rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iMatrix<ltype,N>& lhs,const iScalar<rtype>& rhs)->iMatrix<decltype(lhs._internal[0][0]-rhs._internal),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]-rhs._internal),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iSeries<ltype,N>& lhs,const iSeries<rtype,N>& rhs) ->iSeries<decltype(lhs._internal[0]-rhs._internal[0]),N>
{
    typedef iSeries<decltype(lhs._internal[0]-rhs._internal[0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iScalar<ltype>& lhs,const iSeries<rtype,N>& rhs) ->iSeries<decltype(lhs._internal-rhs._internal[0]),N>
{
    typedef iSeries<decltype(lhs._internal-rhs._internal[0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
strong_inline auto operator - (const iSeries<ltype,N>& lhs,const iScalar<rtype>& rhs) ->iSeries<decltype(lhs._internal[0]-rhs._internal),N>
{
    typedef iSeries<decltype(lhs._internal[0]-rhs._internal),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}

}

#endif
