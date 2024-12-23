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

#ifndef MFEM_GEOM
#define MFEM_GEOM

#include "densemat.hpp"
#include "intrules.hpp"

namespace mfem
{

/** Types of domains for integration rules and reference finite elements:
    Geometry::POINT    - a point
    Geometry::SEGMENT  - the interval [0,1]
    Geometry::TRIANGLE - triangle with vertices (0,0), (1,0), (0,1)
    Geometry::SQUARE   - the unit square (0,1)x(0,1)
    Geometry::TETRAHEDRON - w/ vert. (0,0,0),(1,0,0),(0,1,0),(0,0,1)
    Geometry::CUBE - the unit cube
    Geometry::PRISM - w/ vert. (0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,1,1)
    Geometry::PYRAMID - w/ vert. (0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1)
*/
class Geometry
{
public:
   enum Type
   {
      INVALID = -1,
      POINT = 0, SEGMENT, TRIANGLE, SQUARE, TETRAHEDRON, CUBE, PRISM, PYRAMID,
      NUM_GEOMETRIES
   };

   static const int64_t NumGeom = NUM_GEOMETRIES;
   static const int64_t MaxDim = 3;
   static const int64_t NumBdrArray[NumGeom];
   static const char *Name[NumGeom];
   static const real_t Volume[NumGeom];
   static const int64_t Dimension[NumGeom];
   static const int64_t DimStart[MaxDim+2]; // including MaxDim+1
   static const int64_t NumVerts[NumGeom];
   static const int64_t NumEdges[NumGeom];
   static const int64_t NumFaces[NumGeom];

   // Structure that holds constants describing the Geometries.
   template <Type Geom> struct Constants;

private:
   IntegrationRule *GeomVert[NumGeom];
   IntegrationPoint GeomCenter[NumGeom];
   DenseMatrix *GeomToPerfGeomJac[NumGeom];
   DenseMatrix *PerfGeomToGeomJac[NumGeom];

public:
   Geometry();
   ~Geometry();

   /** @brief Return an IntegrationRule consisting of all vertices of the given
       Geometry::Type, @a GeomType. */
   const IntegrationRule *GetVertices(int64_t GeomType) const;

   /// Return the center of the given Geometry::Type, @a GeomType.
   const IntegrationPoint &GetCenter(int64_t GeomType) const
   { return GeomCenter[GeomType]; }

   /// Get a random point in the reference element specified by @a GeomType.
   /** This method uses the function rand() for random number generation. */
   static void GetRandomPoint(int64_t GeomType, IntegrationPoint &ip);

   /// Check if the given point is inside the given reference element.
   static bool CheckPoint(int64_t GeomType, const IntegrationPoint &ip);
   /** @brief Check if the given point is inside the given reference element.
       Overload for fuzzy tolerance. */
   static bool CheckPoint(int64_t GeomType, const IntegrationPoint &ip,
                          real_t eps);

   /// Project a point @a end, onto the given Geometry::Type, @a GeomType.
   /** Check if the @a end point is inside the reference element, if not
       overwrite it with the point on the boundary that lies on the line segment
       between @a beg and @a end (@a beg must be inside the element). Return
       true if @a end is inside the element, and false otherwise. */
   static bool ProjectPoint(int64_t GeomType, const IntegrationPoint &beg,
                            IntegrationPoint &end);

   /// Project a point @a ip, onto the given Geometry::Type, @a GeomType.
   /** If @a ip is outside the element, replace it with the point on the
       boundary that is closest to the original @a ip and return false;
       otherwise, return true without changing @a ip. */
   static bool ProjectPoint(int64_t GeomType, IntegrationPoint &ip);

   /// Returns true if the given @a geom is of tensor-product type (i.e. if geom
   /// is a segment, quadrilateral, or hexahedron), returns false otherwise.
   static bool IsTensorProduct(Type geom)
   { return geom == SEGMENT || geom == SQUARE || geom == CUBE; }

   /// Returns the Geometry::Type corresponding to a tensor-product of the
   /// given dimension.
   static Type TensorProductGeometry(int64_t dim)
   {
      switch (dim)
      {
         case 0: return POINT;
         case 1: return SEGMENT;
         case 2: return SQUARE;
         case 3: return CUBE;
         default: MFEM_ABORT("Invalid dimension."); return INVALID;
      }
   }

   /// Return the inverse of the given orientation for the specified geometry type.
   static int64_t GetInverseOrientation(Type geom_type, int64_t orientation);

   /// Return the number of boundary "faces" of a given Geometry::Type.
   int64_t NumBdr(int64_t GeomType) const { return NumBdrArray[GeomType]; }
};

template <> struct
   Geometry::Constants<Geometry::POINT>
{
   static const int64_t Dimension = 0;
   static const int64_t NumVert = 1;

   static const int64_t NumOrient = 1;
   static const int64_t Orient[NumOrient][NumVert];
   static const int64_t InvOrient[NumOrient];
};

template <> struct
   Geometry::Constants<Geometry::SEGMENT>
{
   static const int64_t Dimension = 1;
   static const int64_t NumVert = 2;
   static const int64_t NumEdges = 1;
   static const int64_t Edges[NumEdges][2];

   static const int64_t NumOrient = 2;
   static const int64_t Orient[NumOrient][NumVert];
   static const int64_t InvOrient[NumOrient];
};

template <> struct
   Geometry::Constants<Geometry::TRIANGLE>
{
   static const int64_t Dimension = 2;
   static const int64_t NumVert = 3;
   static const int64_t NumEdges = 3;
   static const int64_t Edges[NumEdges][2];
   // Upper-triangular part of the local vertex-to-vertex graph.
   struct VertToVert
   {
      static const int64_t I[NumVert];
      static const int64_t J[NumEdges][2]; // {end,edge_idx}
   };
   static const int64_t NumFaces = 1;
   static const int64_t FaceVert[NumFaces][NumVert];

   // For a given base tuple v={v0,v1,v2}, the orientation of a permutation
   // u={u0,u1,u2} of v, is an index 'j' such that u[i]=v[Orient[j][i]].
   // The static method Mesh::GetTriOrientation, computes the index 'j' of the
   // permutation that maps the second argument 'test' to the first argument
   // 'base': test[Orient[j][i]]=base[i].
   static const int64_t NumOrient = 6;
   static const int64_t Orient[NumOrient][NumVert];
   // The inverse of orientation 'j' is InvOrient[j].
   static const int64_t InvOrient[NumOrient];
};

template <> struct
   Geometry::Constants<Geometry::SQUARE>
{
   static const int64_t Dimension = 2;
   static const int64_t NumVert = 4;
   static const int64_t NumEdges = 4;
   static const int64_t Edges[NumEdges][2];
   // Upper-triangular part of the local vertex-to-vertex graph.
   struct VertToVert
   {
      static const int64_t I[NumVert];
      static const int64_t J[NumEdges][2]; // {end,edge_idx}
   };
   static const int64_t NumFaces = 1;
   static const int64_t FaceVert[NumFaces][NumVert];

   static const int64_t NumOrient = 8;
   static const int64_t Orient[NumOrient][NumVert];
   static const int64_t InvOrient[NumOrient];
};

template <> struct
   Geometry::Constants<Geometry::TETRAHEDRON>
{
   static const int64_t Dimension = 3;
   static const int64_t NumVert = 4;
   static const int64_t NumEdges = 6;
   static const int64_t Edges[NumEdges][2];
   static const int64_t NumFaces = 4;
   static const int64_t FaceTypes[NumFaces];
   static const int64_t MaxFaceVert = 3;
   static const int64_t FaceVert[NumFaces][MaxFaceVert];
   // Upper-triangular part of the local vertex-to-vertex graph.
   struct VertToVert
   {
      static const int64_t I[NumVert];
      static const int64_t J[NumEdges][2]; // {end,edge_idx}
   };

   static const int64_t NumOrient = 24;
   static const int64_t Orient[NumOrient][NumVert];
   static const int64_t InvOrient[NumOrient];
};

template <> struct
   Geometry::Constants<Geometry::CUBE>
{
   static const int64_t Dimension = 3;
   static const int64_t NumVert = 8;
   static const int64_t NumEdges = 12;
   static const int64_t Edges[NumEdges][2];
   static const int64_t NumFaces = 6;
   static const int64_t FaceTypes[NumFaces];
   static const int64_t MaxFaceVert = 4;
   static const int64_t FaceVert[NumFaces][MaxFaceVert];
   // Upper-triangular part of the local vertex-to-vertex graph.
   struct VertToVert
   {
      static const int64_t I[NumVert];
      static const int64_t J[NumEdges][2]; // {end,edge_idx}
   };
};

template <> struct
   Geometry::Constants<Geometry::PRISM>
{
   static const int64_t Dimension = 3;
   static const int64_t NumVert = 6;
   static const int64_t NumEdges = 9;
   static const int64_t Edges[NumEdges][2];
   static const int64_t NumFaces = 5;
   static const int64_t FaceTypes[NumFaces];
   static const int64_t MaxFaceVert = 4;
   static const int64_t FaceVert[NumFaces][MaxFaceVert];
   // Upper-triangular part of the local vertex-to-vertex graph.
   struct VertToVert
   {
      static const int64_t I[NumVert];
      static const int64_t J[NumEdges][2]; // {end,edge_idx}
   };
};

template <> struct
   Geometry::Constants<Geometry::PYRAMID>
{
   static const int64_t Dimension = 3;
   static const int64_t NumVert = 5;
   static const int64_t NumEdges = 8;
   static const int64_t Edges[NumEdges][2];
   static const int64_t NumFaces = 5;
   static const int64_t FaceTypes[NumFaces];
   static const int64_t MaxFaceVert = 4;
   static const int64_t FaceVert[NumFaces][MaxFaceVert];
   // Upper-triangular part of the local vertex-to-vertex graph.
   struct VertToVert
   {
      static const int64_t I[NumVert];
      static const int64_t J[NumEdges][2]; // {end,edge_idx}
   };
};

// Defined in fe.cpp to ensure construction after 'mfem::TriangleFE' and
// `mfem::TetrahedronFE`.
extern Geometry Geometries;


class RefinedGeometry
{
public:
   int64_t Times, ETimes;
   IntegrationRule RefPts;
   Array<int64_t> RefGeoms, RefEdges;
   int64_t NumBdrEdges; // at the beginning of RefEdges
   int64_t Type;

   RefinedGeometry(int64_t NPts, int64_t NRefG, int64_t NRefE, int64_t NBdrE = 0) :
      RefPts(NPts), RefGeoms(NRefG), RefEdges(NRefE), NumBdrEdges(NBdrE) {}
};

}

#endif
