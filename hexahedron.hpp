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

#ifndef MFEM_HEXAHEDRON
#define MFEM_HEXAHEDRON

#include "element.hpp"

namespace mfem
{

/// Data type hexahedron element
class Hexahedron : public Element
{
protected:
   int64_t indices[8];

public:
   typedef Geometry::Constants<Geometry::CUBE> geom_t;

   Hexahedron() : Element(Geometry::CUBE) { }

   /// Constructs hexahedron by specifying the indices and the attribute.
   Hexahedron(const int64_t *ind, int64_t attr = 1);

   /// Constructs hexahedron by specifying the indices and the attribute.
   Hexahedron(int64_t ind1, int64_t ind2, int64_t ind3, int64_t ind4,
              int64_t ind5, int64_t ind6, int64_t ind7, int64_t ind8, int64_t attr = 1);

   /// Return element's type
   Type GetType() const override { return Element::HEXAHEDRON; }

   /// Get the indices defining the vertices.
   void GetVertices(Array<int64_t> &v) const override;

   /// Set the indices defining the vertices.
   void SetVertices(const Array<int64_t> &v) override;

   /// @note The returned array should NOT be deleted by the caller.
   int64_t * GetVertices () override { return indices; }

   /// Set the indices defining the vertices.
   void SetVertices(const int64_t *ind) override;

   int64_t GetNVertices() const override { return 8; }

   int64_t GetNEdges() const override { return 12; }

   const int64_t *GetEdgeVertices(int64_t ei) const override
   { return geom_t::Edges[ei]; }

   int64_t GetNFaces() const override { return 6; }

   int64_t GetNFaceVertices(int64_t) const override { return 4; }

   const int64_t *GetFaceVertices(int64_t fi) const override
   { return geom_t::FaceVert[fi]; }

   Element *Duplicate(Mesh *m) const override
   { return new Hexahedron(indices, attribute); }

   virtual ~Hexahedron() = default;
};

}

#endif
