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


#include "quadrilateral.hpp"

namespace mfem
{

Quadrilateral::Quadrilateral( const int64_t *ind, int64_t attr )
   : Element(Geometry::SQUARE)
{
   attribute = attr;
   for (int64_t i=0; i<4; i++)
   {
      indices[i] = ind[i];
   }
}

Quadrilateral::Quadrilateral( int64_t ind1, int64_t ind2, int64_t ind3,
                              int64_t ind4,
                              int64_t attr ) : Element(Geometry::SQUARE)
{
   attribute  = attr;
   indices[0] = ind1;
   indices[1] = ind2;
   indices[2] = ind3;
   indices[3] = ind4;
}

void Quadrilateral::SetVertices(const int64_t *ind)
{
   std::copy(ind, ind + 4, indices);
}

void Quadrilateral::GetVertices(Array<int64_t> &v) const
{
   v.SetSize(4);
   std::copy(indices, indices + 4, v.begin());
}


void Quadrilateral::SetVertices(const Array<int64_t> &v)
{
   MFEM_ASSERT(v.Size() == 4, "!");
   std::copy(v.begin(), v.end(), indices);
}

}
