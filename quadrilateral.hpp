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

#ifndef MFEM_QUADRILATERAL
#define MFEM_QUADRILATERAL

#include "element.hpp"

namespace mfem
{

/// Data type quadrilateral element
class Quadrilateral : public Element
{
protected:
   int indices[4];

public:
   typedef Geometry::Constants<Geometry::SQUARE> geom_t;

   Quadrilateral() : Element(Geometry::SQUARE) {}

   /// Constructs quadrilateral by specifying the indices and the attribute.
   Quadrilateral(const int *ind, int attr = 1);

   /// Constructs quadrilateral by specifying the indices and the attribute.
   Quadrilateral(int ind1, int ind2, int ind3, int ind4, int attr = 1);

   /// Return element's type
   Type GetType() const override { return Element::QUADRILATERAL; }

   /// Get the indices defining the vertices.
   void GetVertices(Array<int> &v) const override;

   /// Set the indices defining the vertices.
   void SetVertices(const Array<int> &v) override;

   /// @note The returned array should NOT be deleted by the caller.
   int * GetVertices () override { return indices; }

   /// Set the vertices according to the given input.
   void SetVertices(const int *ind) override;

   int GetNVertices() const override { return 4; }

   int GetNEdges() const override { return (4); }

   const int *GetEdgeVertices(int ei) const override
   { return geom_t::Edges[ei]; }

   int GetNFaces() const override { return 0; }

   int GetNFaceVertices(int) const override { return 0; }

   const int *GetFaceVertices(int fi) const override { return NULL; }

   Element *Duplicate(Mesh *m) const override
   { return new Quadrilateral(indices, attribute); }

   virtual ~Quadrilateral() = default;
};

}

#endif
