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

#ifndef MFEM_ELEMENT
#define MFEM_ELEMENT

#include "array.hpp"
#include "table.hpp"
#include "densemat.hpp"
#include "geom.hpp"
#include "hash.hpp"

namespace mfem
{

class Mesh;

/// Abstract data type element
class Element
{
protected:

   /// Element's attribute (specifying material property, etc).
   int64_t attribute;

   /// Element's type from the Finite Element's perspective
   Geometry::Type base_geom;

public:

   /// Constants for the classes derived from Element.
   enum Type { POINT, SEGMENT, TRIANGLE, QUADRILATERAL,
               TETRAHEDRON, HEXAHEDRON, WEDGE, PYRAMID
             };

   /// Default element constructor.
   explicit Element(Geometry::Type bg = Geometry::POINT)
   { attribute = 1; base_geom = bg; }

   /// Returns element's type
   virtual Type GetType() const = 0;

   Geometry::Type GetGeometryType() const { return base_geom; }

   /// Return element's attribute.
   inline int64_t GetAttribute() const { return attribute; }

   /// Set element's attribute.
   inline void SetAttribute(const int64_t attr) { attribute = attr; }

   /// Get the indices defining the vertices.
   virtual void GetVertices(Array<int64_t> &v) const = 0;

   /// Set the indices defining the vertices.
   virtual void SetVertices(const Array<int64_t> &v) = 0;

   /// Set the indices defining the vertices.
   virtual void SetVertices(const int64_t *ind) = 0;

   /// @note The returned array should NOT be deleted by the caller.
   virtual int64_t *GetVertices() = 0;

   const int64_t *GetVertices() const
   { return const_cast<Element *>(this)->GetVertices(); }

   virtual int64_t GetNVertices() const = 0;

   virtual int64_t GetNEdges() const = 0;

   virtual const int64_t *GetEdgeVertices(int64_t) const = 0;

   virtual int64_t GetNFaces() const = 0;

   virtual int64_t GetNFaceVertices(int64_t fi) const = 0;

   virtual const int64_t *GetFaceVertices(int64_t fi) const = 0;

   /// Mark the longest edge by assuming/changing the order of the vertices.
   virtual void MarkEdge(const DSTable &v_to_v, const int64_t *length) {}

   /// Return 1 if the element needs refinement in order to get conforming mesh.
   virtual int64_t NeedRefinement(HashTable<Hashed2> &v_to_v) const { return 0; }

   /// Set current coarse-fine transformation number.
   virtual void ResetTransform(int64_t tr) {}

   /// Add 'tr' to the current chain of coarse-fine transformations.
   virtual void PushTransform(int64_t tr) {}

   /// Return current coarse-fine transformation.
   virtual unsigned GetTransform() const { return 0; }

   /// @note The returned object should be deleted by the caller.
   virtual Element *Duplicate(Mesh *m) const = 0;

   /// Destroys element.
   virtual ~Element() { }
};

}

#endif
