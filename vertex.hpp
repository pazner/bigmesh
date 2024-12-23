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

#ifndef MFEM_VERTEX
#define MFEM_VERTEX

#include "globals.hpp"

namespace mfem
{

/// Data type for vertex
class Vertex
{
protected:
   real_t coord[3];

public:
   Vertex() = default;

   // Trivial copy constructor and trivial copy assignment operator

   Vertex(real_t *xx, int64_t dim);
   Vertex(real_t x, real_t y) { coord[0] = x; coord[1] = y; coord[2] = 0.; }
   Vertex(real_t x, real_t y, real_t z)
   { coord[0] = x; coord[1] = y; coord[2] = z; }

   /// Returns pointer to the coordinates of the vertex.
   inline real_t * operator() () const { return (real_t*)coord; }

   /// Returns the i'th coordinate of the vertex.
   inline real_t & operator() (int64_t i) { return coord[i]; }

   /// Returns the i'th coordinate of the vertex.
   inline const real_t & operator() (int64_t i) const { return coord[i]; }

   /// Sets vertex location based on given point p
   void SetCoords(int64_t dim, const real_t *p)
   { for (int64_t i = 0; i < dim; i++) { coord[i] = p[i]; } }

   // Trivial destructor
};

}

#endif
