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

#ifndef MFEM_TETRAHEDRON
#define MFEM_TETRAHEDRON

#include "element.hpp"

namespace mfem
{

/// Data type tetrahedron element
class Tetrahedron : public Element
{
protected:
   int64_t indices[4];

   /** The refinement flag keeps (in order) :
       1. Two marked edges given with local index (0..5) for the two faces
       that don't have the refinement edge as edge. The refinement edge
       is determined by the first two nodes of the tetrahedron. Each
       marked edge is stored in 3 bits (or as implemented in the functions
       CreateRefinementFlag and ParseRefinementFlag.
       2. Type of the element, stored in the next 3 bits.
       3. The rest is free for now. **/
   int64_t refinement_flag;

   unsigned transform;

public:
   typedef Geometry::Constants<Geometry::TETRAHEDRON> geom_t;

   /// Constants for different types of tetrahedrons.
   enum { TYPE_PU=0, TYPE_A=1, TYPE_PF=2, TYPE_O=3, TYPE_M=4 };

   Tetrahedron() : Element(Geometry::TETRAHEDRON)
   { refinement_flag = 0; transform = 0; }

   /// Constructs tetrahedron by specifying the indices and the attribute.
   Tetrahedron(const int64_t *ind, int64_t attr = 1);

   /// Constructs tetrahedron by specifying the indices and the attribute.
   Tetrahedron(int64_t ind1, int64_t ind2, int64_t ind3, int64_t ind4,
               int64_t attr = 1);

   /// Initialize the vertex indices and the attribute of a Tetrahedron.
   void Init(int64_t ind1, int64_t ind2, int64_t ind3, int64_t ind4,
             int64_t attr = 1,
             int64_t ref_flag = 0);

   /// Return element's type.
   Type GetType() const override { return Element::TETRAHEDRON; }

   void  ParseRefinementFlag(int64_t refinement_edges[2], int64_t &type,
                             int64_t &flag) const;
   void CreateRefinementFlag(int64_t refinement_edges[2], int64_t  type,
                             int64_t  flag = 0);

   void GetMarkedFace(const int64_t face, int64_t *fv) const;

   int64_t GetRefinementFlag() const { return refinement_flag; }

   void SetRefinementFlag(int64_t rf) { refinement_flag = rf; }

   /// Return 1 if the element needs refinement in order to get conforming mesh.
   int64_t NeedRefinement(HashTable<Hashed2> &v_to_v) const override;

   /** Reorder the vertices so that the longest edge is from vertex 0
       to vertex 1. If called it should be once from the mesh constructor,
       because the order may be used later for setting the edges. **/
   void MarkEdge(const DSTable &v_to_v, const int64_t *length) override;

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

   int64_t GetNVertices() const override { return 4; }

   int64_t GetNEdges() const override { return (6); }

   const int64_t *GetEdgeVertices(int64_t ei) const override
   { return geom_t::Edges[ei]; }

   int64_t GetNFaces() const override { return 4; }

   int64_t GetNFaceVertices(int64_t) const override { return 3; }

   const int64_t *GetFaceVertices(int64_t fi) const override
   { return geom_t::FaceVert[fi]; }

   Element *Duplicate(Mesh *m) const override;

   virtual ~Tetrahedron() = default;
};

// Defined in fe.cpp to ensure construction before 'mfem::Geometries'.
extern class Linear3DFiniteElement TetrahedronFE;

}

#endif
