// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_TEMPLATE_LAYOUT
#define MFEM_TEMPLATE_LAYOUT

#include "../config/tconfig.hpp"
#include "../fem/fespace.hpp"
#include "../general/backends.hpp"

namespace mfem
{

// Layout classes

template <int64_t N1, int64_t S1>
struct OffsetStridedLayout1D;
template <int64_t N1, int64_t S1, int64_t N2, int64_t S2>
struct StridedLayout2D;

template <int64_t N1, int64_t S1>
struct StridedLayout1D
{
   static const int64_t rank = 1;
   static const int64_t dim_1 = N1;
   static const int64_t size = N1;

   MFEM_HOST_DEVICE static inline int64_t ind(int64_t i1)
   {
      return S1*i1;
   }

   template <int64_t M1>
   static OffsetStridedLayout1D<M1,S1> sub(int64_t o1)
   {
      return OffsetStridedLayout1D<M1,S1>(S1*o1);
   }

   // reshape methods

   template <int64_t N1_1, int64_t N1_2>
   static StridedLayout2D<N1_1,S1,N1_2,S1*N1_1> split_1()
   {
      // S1*i1 == S1*(i1_1+N1_1*i1_2)
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1, "invalid dimensions");
      return StridedLayout2D<N1_1,S1,N1_2,S1*N1_1>();
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2>
struct OffsetStridedLayout2D;

template <int64_t N1, int64_t S1>
struct OffsetStridedLayout1D
{
   static const int64_t rank = 1;
   static const int64_t dim_1 = N1;
   static const int64_t size = N1;

   int64_t offset;

   OffsetStridedLayout1D() { }
   OffsetStridedLayout1D(int64_t offset_) : offset(offset_) { }
   MFEM_HOST_DEVICE inline int64_t ind(int64_t i1) const
   {
      return offset+S1*i1;
   }

   template <int64_t M1>
   OffsetStridedLayout1D<M1,S1> sub(int64_t o1) const
   {
      return OffsetStridedLayout1D<M1,S1>(offset+S1*o1);
   }

   // reshape methods

   template <int64_t N1_1, int64_t N1_2>
   OffsetStridedLayout2D<N1_1,S1,N1_2,S1*N1_1> split_1() const
   {
      // S1*i1 == S1*(i1_1+N1_1*i1_2)
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1, "invalid dimensions");
      return OffsetStridedLayout2D<N1_1,S1,N1_2,S1*N1_1>(offset);
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3>
struct StridedLayout3D;
template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3, int64_t N4, int64_t S4>
struct StridedLayout4D;

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2>
struct StridedLayout2D
{
   static const int64_t rank = 2;
   static const int64_t dim_1 = N1;
   static const int64_t dim_2 = N2;
   static const int64_t size = N1*N2;

   MFEM_HOST_DEVICE static inline int64_t ind(int64_t i1, int64_t i2)
   {
      return (S1*i1+S2*i2);
   }
   static OffsetStridedLayout1D<N2,S2> ind1(int64_t i1)
   {
      return OffsetStridedLayout1D<N2,S2>(S1*i1);
   }
   static OffsetStridedLayout1D<N1,S1> ind2(int64_t i2)
   {
      return OffsetStridedLayout1D<N1,S1>(S2*i2);
   }

   template <int64_t M1, int64_t M2>
   static OffsetStridedLayout2D<M1,S1,M2,S2> sub(int64_t o1, int64_t o2)
   {
      return OffsetStridedLayout2D<M1,S1,M2,S2>(S1*o1+S2*o2);
   }

   // reshape methods

   template <int64_t N1_1, int64_t N1_2>
   static StridedLayout3D<N1_1,S1,N1_2,S1*N1_1,N2,S2> split_1()
   {
      // S1*i1+S2*i2 == S1*(i1_1+N1_1*i1_2)+S2*i2
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1, "invalid dimensions");
      return StridedLayout3D<N1_1,S1,N1_2,S1*N1_1,N2,S2>();
   }
   template <int64_t N2_1, int64_t N2_2>
   static StridedLayout3D<N1,S1,N2_1,S2,N2_2,S2*N2_1> split_2()
   {
      // S1*i1+S2*i2 == S1*i1+S2*(i2_1*N2_1*i2_2)
      MFEM_STATIC_ASSERT(N2_1*N2_2 == N2, "invalid dimensions");
      return StridedLayout3D<N1,S1,N2_1,S2,N2_2,S2*N2_1>();
   }
   template <int64_t N1_1, int64_t N1_2, int64_t N2_1, int64_t N2_2>
   static StridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2_1,S2,N2_2,S2*N2_1> split_12()
   {
      // S1*i1+S2*i2 == S1*(i1_1+N1_1*i1_2)+S2*(i2_1+N2_1*i2_2)
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1 && N2_1*N2_2 == N2,
                         "invalid dimensions");
      return StridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2_1,S2,N2_2,S2*N2_1>();
   }
   static StridedLayout1D<N1*N2,(S1<S2)?S1:S2> merge_12()
   {
      // use: (S1*i1+S2*i2) == (S1*(i1+S2/S1*i2))
      //  or  (S1*i1+S2*i2) == (S2*(S1/S2*i1+i2))
      // assuming: S2 == S1*N1 || S1 == S2*N2
      MFEM_STATIC_ASSERT(S2 == S1*N1 || S1 == S2*N2, "invalid reshape");
      return StridedLayout1D<N1*N2,(S1<S2)?S1:S2>();
   }
   static StridedLayout2D<N2,S2,N1,S1> transpose_12()
   {
      return StridedLayout2D<N2,S2,N1,S1>();
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3>
struct OffsetStridedLayout3D;
template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3, int64_t N4, int64_t S4>
struct OffsetStridedLayout4D;

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2>
struct OffsetStridedLayout2D
{
   static const int64_t rank = 2;
   static const int64_t dim_1 = N1;
   static const int64_t dim_2 = N2;
   static const int64_t size = N1*N2;

   int64_t offset;

   OffsetStridedLayout2D() { }
   OffsetStridedLayout2D(int64_t offset_) : offset(offset_) { }
   MFEM_HOST_DEVICE inline int64_t ind(int64_t i1, int64_t i2) const
   {
      return offset+S1*i1+S2*i2;
   }
   OffsetStridedLayout1D<N2,S2> ind1(int64_t i1) const
   {
      return OffsetStridedLayout1D<N2,S2>(offset+S1*i1);
   }
   OffsetStridedLayout1D<N1,S1> ind2(int64_t i2) const
   {
      return OffsetStridedLayout1D<N1,S1>(offset+S2*i2);
   }

   template <int64_t M1, int64_t M2>
   OffsetStridedLayout2D<M1,S1,M2,S2> sub(int64_t o1, int64_t o2) const
   {
      return OffsetStridedLayout2D<M1,S1,M2,S2>(offset+S1*o1+S2*o2);
   }

   // reshape methods

   template <int64_t N1_1, int64_t N1_2>
   OffsetStridedLayout3D<N1_1,S1,N1_2,S1*N1_1,N2,S2> split_1() const
   {
      // S1*i1+S2*i2 == S1*(i1_1+N1_1*i1_2)+S2*i2
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1, "invalid dimensions");
      return OffsetStridedLayout3D<N1_1,S1,N1_2,S1*N1_1,N2,S2>(offset);
   }
   template <int64_t N2_1, int64_t N2_2>
   OffsetStridedLayout3D<N1,S1,N2_1,S2,N2_2,S2*N2_1> split_2() const
   {
      // S1*i1+S2*i2 == S1*i1+S2*(i2_1*N2_1*i2_2)
      MFEM_STATIC_ASSERT(N2_1*N2_2 == N2, "invalid dimensions");
      return OffsetStridedLayout3D<N1,S1,N2_1,S2,N2_2,S2*N2_1>(offset);
   }
   template <int64_t N1_1, int64_t N1_2, int64_t N2_1, int64_t N2_2>
   OffsetStridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2_1,S2,N2_2,S2*N2_1>
   split_12() const
   {
      // S1*i1+S2*i2 == S1*(i1_1+N1_1*i1_2)+S2*(i2_1+N2_1*i2_2)
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1 && N2_1*N2_2 == N2,
                         "invalid dimensions");
      return OffsetStridedLayout4D<
             N1_1,S1,N1_2,S1*N1_1,N2_1,S2,N2_2,S2*N2_1>(offset);
   }
   OffsetStridedLayout1D<N1*N2,(S1<S2)?S1:S2> merge_12() const
   {
      // use: (S1*i1+S2*i2) == (S1*(i1+S2/S1*i2))
      //  or  (S1*i1+S2*i2) == (S2*(S1/S2*i1+i2))
      // assuming: S2 == S1*N1 || S1 == S2*N2
      MFEM_STATIC_ASSERT(S2 == S1*N1 || S1 == S2*N2, "invalid reshape");
      return OffsetStridedLayout1D<N1*N2,(S1<S2)?S1:S2>(offset);
   }
   OffsetStridedLayout2D<N2,S2,N1,S1> transpose_12() const
   {
      return OffsetStridedLayout2D<N2,S2,N1,S1>(offset);
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3>
struct StridedLayout3D
{
   static const int64_t rank = 3;
   static const int64_t dim_1 = N1;
   static const int64_t dim_2 = N2;
   static const int64_t dim_3 = N3;
   static const int64_t size = N1*N2*N3;

   static inline int64_t ind(int64_t i1, int64_t i2, int64_t i3)
   {
      return S1*i1+S2*i2+S3*i3;
   }
   static OffsetStridedLayout2D<N2,S2,N3,S3> ind1(int64_t i1)
   {
      return OffsetStridedLayout2D<N2,S2,N3,S3>(S1*i1);
   }
   static OffsetStridedLayout2D<N1,S1,N3,S3> ind2(int64_t i2)
   {
      return OffsetStridedLayout2D<N1,S1,N3,S3>(S2*i2);
   }
   static OffsetStridedLayout2D<N1,S1,N2,S2> ind3(int64_t i3)
   {
      return OffsetStridedLayout2D<N1,S1,N2,S2>(S3*i3);
   }

   // reshape methods

   static StridedLayout2D<N1*N2,S1,N3,S3> merge_12()
   {
      // use: (S1*i1+S2*i2+S3*i3) == (S1*(i1+S2/S1*i2)+S3*i3)
      // assuming: S2 == S1*N1
      MFEM_STATIC_ASSERT(S2 == S1*N1, "invalid reshape");
      return StridedLayout2D<N1*N2,S1,N3,S3>();
      // alternative:
      // use: (S1*i1+S2*i2+S3*i3) == (S2*(S1/S2*i1+i2)+S3*i3)
      // assuming: S1 == S2*N2
      // result is: StridedLayout2D<N1*N2,S2,N3,S3>
   }
   static StridedLayout2D<N1,S1,N2*N3,S2> merge_23()
   {
      // use: (S1*i1+S2*i2+S3*i3) == (S1*i1+S2*(i2+S3/S2*i3))
      // assuming: S3 == S2*N2
      MFEM_STATIC_ASSERT(S3 == S2*N2, "invalid reshape");
      return StridedLayout2D<N1,S1,N2*N3,S2>();
   }

   template <int64_t N1_1, int64_t N1_2>
   static StridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2,S2,N3,S3> split_1()
   {
      // S1*i1+S2*i2+S3*i3 == S1*(i1_1+N1_1*i1_2)+S2*i2+S3*i3
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1, "invalid dimensions");
      return StridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2,S2,N3,S3>();
   }
   template <int64_t N2_1, int64_t N2_2>
   static StridedLayout4D<N1,S1,N2_1,S2,N2_2,S2*N2_1,N3,S3> split_2()
   {
      // S1*i1+S2*i2+S3*i3 == S1*i1+S2*(i2_1+N2_1*i2_2)+S3*i3
      MFEM_STATIC_ASSERT(N2_1*N2_2 == N2, "invalid dimensions");
      return StridedLayout4D<N1,S1,N2_1,S2,N2_2,S2*N2_1,N3,S3>();
   }
   template <int64_t N3_1, int64_t N3_2>
   static StridedLayout4D<N1,S1,N2,S2,N3_1,S3,N3_2,S3*N3_1> split_3()
   {
      // S1*i1+S2*i2+S3*i3 == S1*i1+S2*i2+S3*(i3_1+N3_1*i3_2)
      MFEM_STATIC_ASSERT(N3_1*N3_2 == N3, "invalid dimensions");
      return StridedLayout4D<N1,S1,N2,S2,N3_1,S3,N3_2,S3*N3_1>();
   }

   static StridedLayout3D<N2,S2,N1,S1,N3,S3> transpose_12()
   {
      return StridedLayout3D<N2,S2,N1,S1,N3,S3>();
   }
   static StridedLayout3D<N3,S3,N2,S2,N1,S1> transpose_13()
   {
      return StridedLayout3D<N3,S3,N2,S2,N1,S1>();
   }
   static StridedLayout3D<N1,S1,N3,S3,N2,S2> transpose_23()
   {
      return StridedLayout3D<N1,S1,N3,S3,N2,S2>();
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3>
struct OffsetStridedLayout3D
{
   static const int64_t rank = 3;
   static const int64_t dim_1 = N1;
   static const int64_t dim_2 = N2;
   static const int64_t dim_3 = N3;
   static const int64_t size = N1*N2*N3;

   int64_t offset;

   OffsetStridedLayout3D() { }
   OffsetStridedLayout3D(int64_t offset_) : offset(offset_) { }
   inline int64_t ind(int64_t i1, int64_t i2, int64_t i3) const
   {
      return offset+S1*i1+S2*i2+S3*i3;
   }
   OffsetStridedLayout2D<N2,S2,N3,S3> ind1(int64_t i1) const
   {
      return OffsetStridedLayout2D<N2,S2,N3,S3>(offset+S1*i1);
   }
   OffsetStridedLayout2D<N1,S1,N3,S3> ind2(int64_t i2) const
   {
      return OffsetStridedLayout2D<N1,S1,N3,S3>(offset+S2*i2);
   }
   OffsetStridedLayout2D<N1,S1,N2,S2> ind3(int64_t i3) const
   {
      return OffsetStridedLayout2D<N1,S1,N2,S2>(offset+S3*i3);
   }

   // reshape methods

   OffsetStridedLayout2D<N1*N2,S1,N3,S3> merge_12() const
   {
      // use: (S1*i1+S2*i2+S3*i3) == (S1*(i1+S2/S1*i2)+S3*i3)
      // assuming: S2 == S1*N1
      MFEM_STATIC_ASSERT(S2 == S1*N1, "invalid reshape");
      return OffsetStridedLayout2D<N1*N2,S1,N3,S3>(offset);
   }
   OffsetStridedLayout2D<N1,S1,N2*N3,S2> merge_23() const
   {
      // use: (S1*i1+S2*i2+S3*i3) == (S1*i1+S2*(i2+S3/S2*i3))
      // assuming: S3 == S2*N2
      MFEM_STATIC_ASSERT(S3 == S2*N2, "invalid reshape");
      return OffsetStridedLayout2D<N1,S1,N2*N3,S2>(offset);
   }

   template <int64_t N1_1, int64_t N1_2>
   OffsetStridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2,S2,N3,S3> split_1() const
   {
      // S1*i1+S2*i2+S3*i3 == S1*(i1_1+N1_1*i1_2)+S2*i2+S3*i3
      MFEM_STATIC_ASSERT(N1_1*N1_2 == N1, "invalid dimensions");
      return OffsetStridedLayout4D<N1_1,S1,N1_2,S1*N1_1,N2,S2,N3,S3>(offset);
   }
   template <int64_t N2_1, int64_t N2_2>
   OffsetStridedLayout4D<N1,S1,N2_1,S2,N2_2,S2*N2_1,N3,S3> split_2() const
   {
      // S1*i1+S2*i2+S3*i3 == S1*i1+S2*(i2_1+N2_1*i2_2)+S3*i3
      MFEM_STATIC_ASSERT(N2_1*N2_2 == N2, "invalid dimensions");
      return OffsetStridedLayout4D<N1,S1,N2_1,S2,N2_2,S2*N2_1,N3,S3>(offset);
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3, int64_t N4, int64_t S4>
struct StridedLayout4D
{
   static const int64_t rank = 4;
   static const int64_t dim_1 = N1;
   static const int64_t dim_2 = N2;
   static const int64_t dim_3 = N3;
   static const int64_t dim_4 = N4;
   static const int64_t size = N1*N2*N3*N4;

   static inline int64_t ind(int64_t i1, int64_t i2, int64_t i3, int64_t i4)
   {
      return S1*i1+S2*i2+S3*i3+S4*i4;
   }
   static OffsetStridedLayout2D<N1,S1,N4,S4> ind23(int64_t i2, int64_t i3)
   {
      return OffsetStridedLayout2D<N1,S1,N4,S4>(S2*i2+S3*i3);
   }
   static OffsetStridedLayout2D<N2,S2,N3,S3> ind14(int64_t i1, int64_t i4)
   {
      return OffsetStridedLayout2D<N2,S2,N3,S3>(S1*i1+S4*i4);
   }
   static OffsetStridedLayout3D<N1,S1,N2,S2,N3,S3> ind4(int64_t i4)
   {
      return OffsetStridedLayout3D<N1,S1,N2,S2,N3,S3>(S4*i4);
   }

   static StridedLayout3D<N1*N2,S1,N3,S3,N4,S4> merge_12()
   {
      // use: (S1*i1+S2*i2+S3*i3+S4*i4) == (S1*(i1+S2/S1*i2)+S3*i3+S4*i4)
      // assuming S2 == S1*N1
      MFEM_STATIC_ASSERT(S2 == S1*N1, "invalid reshape");
      return StridedLayout3D<N1*N2,S1,N3,S3,N4,S4>();
   }
   static StridedLayout3D<N1,S1,N2,S2,N3*N4,S3> merge_34()
   {
      // use: (S1*i1+S2*i2+S3*i3+S4*i4) == (S1*i1+S2*i2+S3*(i3+S4/S3*i4))
      // assuming S4 == S3*N3
      MFEM_STATIC_ASSERT(S4 == S3*N3, "invalid reshape");
      return StridedLayout3D<N1,S1,N2,S2,N3*N4,S3>();
   }
};

template <int64_t N1, int64_t S1, int64_t N2, int64_t S2, int64_t N3, int64_t S3, int64_t N4, int64_t S4>
struct OffsetStridedLayout4D
{
   static const int64_t rank = 4;
   static const int64_t dim_1 = N1;
   static const int64_t dim_2 = N2;
   static const int64_t dim_3 = N3;
   static const int64_t dim_4 = N4;
   static const int64_t size = N1*N2*N3*N4;

   int64_t offset;

   OffsetStridedLayout4D() { }
   OffsetStridedLayout4D(int64_t offset_) : offset(offset_) { }
   inline int64_t ind(int64_t i1, int64_t i2, int64_t i3, int64_t i4) const
   {
      return offset+S1*i1+S2*i2+S3*i3+S4*i4;
   }
};

template <int64_t N1, int64_t N2>
struct ColumnMajorLayout2D
   : public StridedLayout2D<N1,1,N2,N1> { };

template <int64_t N1, int64_t N2, int64_t N3>
struct ColumnMajorLayout3D
   : public StridedLayout3D<N1,1,N2,N1,N3,N1*N2> { };

template <int64_t N1, int64_t N2, int64_t N3, int64_t N4>
struct ColumnMajorLayout4D
   : public StridedLayout4D<N1,1,N2,N1,N3,N1*N2,N4,N1*N2*N3> { };


// Vector layout classes

class DynamicVectorLayout
{
public:
   static const int64_t vec_dim = 0; // 0 - dynamic

protected:
   int64_t scal_stride, comp_stride;
   int64_t num_components;

   void Init(Ordering::Type ordering, int64_t scalar_size, int64_t num_comp)
   {
      num_components = num_comp;
      if (ordering == Ordering::byNODES)
      {
         scal_stride = 1;
         comp_stride = scalar_size;
      }
      else
      {
         scal_stride = num_comp;
         comp_stride = 1;
      }
   }

public:
   DynamicVectorLayout(Ordering::Type ordering, int64_t scalar_size,
                       int64_t num_comp)
   {
      Init(ordering, scalar_size, num_comp);
   }
   DynamicVectorLayout(const FiniteElementSpace &fes)
   {
      Init(fes.GetOrdering(), fes.GetNDofs(), fes.GetVDim());
   }
   // default copy constructor

   int64_t NumComponents() const { return num_components; }

   int64_t ind(int64_t scalar_idx, int64_t comp_idx) const
   {
      return scal_stride * scalar_idx + comp_stride * comp_idx;
   }

   static bool Matches(const FiniteElementSpace &fes)
   {
      return true;
   }
};

// The default value (NumComp = 0) indicates that the number of components is
// dynamic, i.e. it will be specified at run-time.
template <Ordering::Type Ord, int64_t NumComp = 0>
class VectorLayout
{
public:
   static const int64_t vec_dim = NumComp;

protected:
   int64_t num_components, scalar_size;

public:
   VectorLayout(int64_t scalar_size_, int64_t num_comp_ = NumComp)
      : num_components(num_comp_),
        scalar_size(scalar_size_)
   {
      MFEM_ASSERT(NumComp == 0 || num_components == NumComp,
                  "invalid number of components");
   }

   VectorLayout(const FiniteElementSpace &fes)
      : num_components(fes.GetVDim()),
        scalar_size(fes.GetNDofs())
   {
      MFEM_ASSERT(fes.GetOrdering() == Ord, "ordering mismatch");
      MFEM_ASSERT(NumComp == 0 || num_components == NumComp,
                  "invalid number of components");
   }
   // default copy constructor

   int64_t NumComponents() const { return (NumComp ? NumComp : num_components); }

   int64_t ind(int64_t scalar_idx, int64_t comp_idx) const
   {
      if (Ord == Ordering::byNODES)
      {
         return scalar_idx + comp_idx * scalar_size;
      }
      else
      {
         return comp_idx + (NumComp ? NumComp : num_components) * scalar_idx;
      }
   }

   static bool Matches(const FiniteElementSpace &fes)
   {
      return (Ord == fes.GetOrdering() &&
              (NumComp == 0 || NumComp == fes.GetVDim()));
   }
};

class ScalarLayout
{
public:
   static const int64_t vec_dim = 1;

   ScalarLayout() { }

   ScalarLayout(const FiniteElementSpace &fes)
   {
      MFEM_ASSERT(fes.GetVDim() == 1, "invalid number of components");
   }

   int64_t NumComponents() const { return 1; }

   int64_t ind(int64_t scalar_idx, int64_t comp_idx) const { return scalar_idx; }

   static bool Matches(const FiniteElementSpace &fes)
   {
      return (fes.GetVDim() == 1);
   }
};

} // namespace mfem

#endif // MFEM_TEMPLATE_LAYOUT
