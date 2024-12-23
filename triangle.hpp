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

#ifndef MFEM_TRIANGLE
#define MFEM_TRIANGLE

#include "element.hpp"

namespace mfem
{

/// Data type triangle element
class Triangle : public Element
{
protected:
   int64_t indices[3];

   unsigned transform;

public:
   typedef Geometry::Constants<Geometry::TRIANGLE> geom_t;

   Triangle() : Element(Geometry::TRIANGLE) { transform = 0; }

   /// Constructs triangle by specifying the indices and the attribute.
   Triangle(const int64_t *ind, int64_t attr = 1);

   /// Constructs triangle by specifying the indices and the attribute.
   Triangle(int64_t ind1, int64_t ind2, int64_t ind3, int64_t attr = 1);

   /// Return element's type.
   Type GetType() const override { return Element::TRIANGLE; }

   /// Return 1 if the element needs refinement in order to get conforming mesh.
   int64_t NeedRefinement(HashTable<Hashed2> &v_to_v) const override;

   /** Reorder the vertices so that the longest edge is from vertex 0
       to vertex 1. If called it should be once from the mesh constructor,
       because the order may be used later for setting the edges. **/
   void MarkEdge(DenseMatrix & pmat);

   static void MarkEdge(int64_t *indices, const DSTable &v_to_v,
                        const int64_t *length);

   /// Mark the longest edge by assuming/changing the order of the vertices.
   void MarkEdge(const DSTable &v_to_v, const int64_t *length) override
   { MarkEdge(indices, v_to_v, length); }

   void ResetTransform(int64_t tr) override { transform = tr; }
   unsigned GetTransform() const override { return transform; }

   /// Add 'tr' to the current chain of coarse-fine transformations.
   void PushTransform(int64_t tr) override
   { transform = (transform << 3) | (tr + 1); }

   /// Calculate point matrix corresponding to a chain of transformations.
   static void GetPointMatrix(unsigned transform, DenseMatrix &pm);

   /// Get the indices defining the vertices.
   void GetVertices(Array<int64_t> &v) const override;

   /// Set the indices defining the vertices.
   void SetVertices(const Array<int64_t> &v) override;

   /// @note The returned array should NOT be deleted by the caller.
   int64_t * GetVertices () override { return indices; }

   /// Set the indices defining the vertices.
   void SetVertices(const int64_t *ind) override;


   int64_t GetNVertices() const override { return 3; }

   int64_t GetNEdges() const override { return (3); }

   const int64_t *GetEdgeVertices(int64_t ei) const override
   { return geom_t::Edges[ei]; }

   int64_t GetNFaces() const override { return 0; }

   int64_t GetNFaceVertices(int64_t) const override { return 0; }

   const int64_t *GetFaceVertices(int64_t fi) const override
   { MFEM_ABORT("not implemented"); return NULL; }

   Element *Duplicate(Mesh *m) const override
   { return new Triangle(indices, attribute); }

   virtual ~Triangle() = default;
};

}

#endif
