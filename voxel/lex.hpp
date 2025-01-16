// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef LEX_HPP
#define LEX_HPP

#include <array>

struct LexIndex
{
   int64_t ndim = 0;
   std::array<int64_t, 3> coords;

   LexIndex() = default;
   LexIndex(int64_t x) : ndim(1), coords({x}) { }
   LexIndex(int64_t x, int64_t y) : ndim(2), coords({x, y}) { }
   LexIndex(int64_t x, int64_t y, int64_t z) : ndim(3), coords({x, y, z}) { }
   LexIndex(const int64_t *xx, int64_t ndim_) : ndim(ndim_)
   {
      std::copy(xx, xx + ndim, coords.begin());
   }

   int64_t operator[](int64_t i) const { return coords[i]; }

   int64_t LinearIndex(const std::vector<int64_t> &n) const
   {
      int64_t shift = 1;
      int64_t idx = 0;
      for (int64_t i = 0; i < ndim; ++i)
      {
         idx += coords[i]*shift;
         shift *= n[i];
      }
      return idx;
   }
};

#endif
