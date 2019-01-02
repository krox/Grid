/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 



    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsm.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#include <Grid/qcd/action/fermion/FermionCore.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////
// Default to no assembler implementation
///////////////////////////////////////////////////////////
template<class Impl> void 
WilsonKernels<Impl >::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
{
  assert(0);
}

template<class Impl> void 
WilsonKernels<Impl >::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					     int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
{
  assert(0);
}

template<class Impl> void 
WilsonKernels<Impl >::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
{
  assert(0);
}

template<class Impl> void 
WilsonKernels<Impl >::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					     int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
{
  assert(0);
}

template<class Impl> void 
WilsonKernels<Impl >::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
{
  assert(0);
}

template<class Impl> void 
WilsonKernels<Impl >::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					     int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out)
{
  assert(0);
}

#include <qcd/action/fermion/WilsonKernelsAsmAvx512.h>
#include <qcd/action/fermion/WilsonKernelsAsmQPX.h>

#define INSTANTIATE_ASM(A) \
template void WilsonKernels<A>::AsmDhopSite(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out);\
 \
template void WilsonKernels<A>::AsmDhopSiteDag(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out);\
template void WilsonKernels<A>::AsmDhopSiteInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out);\
 \
template void WilsonKernels<A>::AsmDhopSiteDagInt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out);\
template void WilsonKernels<A>::AsmDhopSiteExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out);\
 \
template void WilsonKernels<A>::AsmDhopSiteDagExt(StencilView &st, DoubledGaugeFieldView &U, SiteHalfSpinor *buf,\
                                  int ss,int ssU,int Ls,int Ns,const FermionFieldView &in, FermionFieldView &out);\

//INSTANTIATE_ASM(WilsonImplF);
//INSTANTIATE_ASM(WilsonImplD);
INSTANTIATE_ASM(GparityWilsonImplF);
INSTANTIATE_ASM(GparityWilsonImplD);
//INSTANTIATE_ASM(ZWilsonImplF);
//INSTANTIATE_ASM(ZWilsonImplD);
//INSTANTIATE_ASM(DomainWallVec5dImplF);
//INSTANTIATE_ASM(DomainWallVec5dImplD);
//INSTANTIATE_ASM(ZDomainWallVec5dImplF);
//INSTANTIATE_ASM(ZDomainWallVec5dImplD);

//INSTANTIATE_ASM(WilsonImplFH);
//INSTANTIATE_ASM(WilsonImplDF);
//INSTANTIATE_ASM(ZWilsonImplFH);
//INSTANTIATE_ASM(ZWilsonImplDF);
INSTANTIATE_ASM(GparityWilsonImplFH);
INSTANTIATE_ASM(GparityWilsonImplDF);
//INSTANTIATE_ASM(DomainWallVec5dImplFH);
//INSTANTIATE_ASM(DomainWallVec5dImplDF);
//INSTANTIATE_ASM(ZDomainWallVec5dImplFH);
//INSTANTIATE_ASM(ZDomainWallVec5dImplDF);

NAMESPACE_END(Grid);

