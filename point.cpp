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


#include "point.hpp"

namespace mfem
{

Point::Point( const int64_t *ind, int64_t attr ) : Element(Geometry::POINT)
{
   attribute = attr;
   indices[0] = ind[0];
}

void Point::GetVertices(Array<int64_t> &v) const
{
   v.SetSize(1);
   v[0] = indices[0];
}

void Point::SetVertices(const Array<int64_t> &v)
{
   MFEM_ASSERT(v.Size() == 1, "!");
   indices[0] = v[0];
}


void Point::SetVertices(const int64_t *ind)
{
   indices[0] = ind[0];
}

}
