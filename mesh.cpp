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

// Implementation of data type mesh

#include "mesh.hpp"
#include "point.hpp"
#include "segment.hpp"
#include "quadrilateral.hpp"
#include "triangle.hpp"
#include "tetrahedron.hpp"
#include "hexahedron.hpp"
#include "wedge.hpp"
#include "pyramid.hpp"

#include "sort_pairs.hpp"
#include "binaryio.hpp"
#include "text.hpp"
#include "device.hpp"
// #include "tic_toc.hpp"
#include "kdtree.hpp"
#include "sets.hpp"

#include "zstr.hpp"

#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include <functional>
#include <unordered_map>
#include <unordered_set>

// Include the METIS header, if using version 5. If using METIS 4, the needed
// declarations are inlined below, i.e. no header is needed.
#if defined(MFEM_USE_METIS) && defined(MFEM_USE_METIS_5)
#include "metis.h"
#endif

using namespace std;

namespace mfem
{

//// SFC Ordering //////////////////////////////////////////////////////////////

static int64_t sgn(int64_t x)
{
   return (x < 0) ? -1 : (x > 0) ? 1 : 0;
}

static void HilbertSfc2D(int64_t x, int64_t y, int64_t ax, int64_t ay,
                         int64_t bx, int64_t by,
                         Array<int64_t> &coords)
{
   int64_t w = std::abs(ax + ay);
   int64_t h = std::abs(bx + by);

   int64_t dax = sgn(ax), day = sgn(ay); // unit major direction ("right")
   int64_t dbx = sgn(bx), dby = sgn(by); // unit orthogonal direction ("up")

   if (h == 1) // trivial row fill
   {
      for (int64_t i = 0; i < w; i++, x += dax, y += day)
      {
         coords.Append(x);
         coords.Append(y);
      }
      return;
   }
   if (w == 1) // trivial column fill
   {
      for (int64_t i = 0; i < h; i++, x += dbx, y += dby)
      {
         coords.Append(x);
         coords.Append(y);
      }
      return;
   }

   int64_t ax2 = ax/2, ay2 = ay/2;
   int64_t bx2 = bx/2, by2 = by/2;

   int64_t w2 = std::abs(ax2 + ay2);
   int64_t h2 = std::abs(bx2 + by2);

   if (2*w > 3*h) // long case: split in two parts only
   {
      if ((w2 & 0x1) && (w > 2))
      {
         ax2 += dax, ay2 += day; // prefer even steps
      }

      HilbertSfc2D(x, y, ax2, ay2, bx, by, coords);
      HilbertSfc2D(x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by, coords);
   }
   else // standard case: one step up, one long horizontal step, one step down
   {
      if ((h2 & 0x1) && (h > 2))
      {
         bx2 += dbx, by2 += dby; // prefer even steps
      }

      HilbertSfc2D(x, y, bx2, by2, ax2, ay2, coords);
      HilbertSfc2D(x+bx2, y+by2, ax, ay, bx-bx2, by-by2, coords);
      HilbertSfc2D(x+(ax-dax)+(bx2-dbx), y+(ay-day)+(by2-dby),
                   -bx2, -by2, -(ax-ax2), -(ay-ay2), coords);
   }
}

static void HilbertSfc3D(int64_t x, int64_t y, int64_t z,
                         int64_t ax, int64_t ay, int64_t az,
                         int64_t bx, int64_t by, int64_t bz,
                         int64_t cx, int64_t cy, int64_t cz,
                         Array<int64_t> &coords)
{
   int64_t w = std::abs(ax + ay + az);
   int64_t h = std::abs(bx + by + bz);
   int64_t d = std::abs(cx + cy + cz);

   int64_t dax = sgn(ax), day = sgn(ay), daz = sgn(az); // unit major dir ("right")
   int64_t dbx = sgn(bx), dby = sgn(by),
           dbz = sgn(bz); // unit ortho dir ("forward")
   int64_t dcx = sgn(cx), dcy = sgn(cy), dcz = sgn(cz); // unit ortho dir ("up")

   // trivial row/column fills
   if (h == 1 && d == 1)
   {
      for (int64_t i = 0; i < w; i++, x += dax, y += day, z += daz)
      {
         coords.Append(x);
         coords.Append(y);
         coords.Append(z);
      }
      return;
   }
   if (w == 1 && d == 1)
   {
      for (int64_t i = 0; i < h; i++, x += dbx, y += dby, z += dbz)
      {
         coords.Append(x);
         coords.Append(y);
         coords.Append(z);
      }
      return;
   }
   if (w == 1 && h == 1)
   {
      for (int64_t i = 0; i < d; i++, x += dcx, y += dcy, z += dcz)
      {
         coords.Append(x);
         coords.Append(y);
         coords.Append(z);
      }
      return;
   }

   int64_t ax2 = ax/2, ay2 = ay/2, az2 = az/2;
   int64_t bx2 = bx/2, by2 = by/2, bz2 = bz/2;
   int64_t cx2 = cx/2, cy2 = cy/2, cz2 = cz/2;

   int64_t w2 = std::abs(ax2 + ay2 + az2);
   int64_t h2 = std::abs(bx2 + by2 + bz2);
   int64_t d2 = std::abs(cx2 + cy2 + cz2);

   // prefer even steps
   if ((w2 & 0x1) && (w > 2))
   {
      ax2 += dax, ay2 += day, az2 += daz;
   }
   if ((h2 & 0x1) && (h > 2))
   {
      bx2 += dbx, by2 += dby, bz2 += dbz;
   }
   if ((d2 & 0x1) && (d > 2))
   {
      cx2 += dcx, cy2 += dcy, cz2 += dcz;
   }

   // wide case, split in w only
   if ((2*w > 3*h) && (2*w > 3*d))
   {
      HilbertSfc3D(x, y, z,
                   ax2, ay2, az2,
                   bx, by, bz,
                   cx, cy, cz, coords);

      HilbertSfc3D(x+ax2, y+ay2, z+az2,
                   ax-ax2, ay-ay2, az-az2,
                   bx, by, bz,
                   cx, cy, cz, coords);
   }
   // do not split in d
   else if (3*h > 4*d)
   {
      HilbertSfc3D(x, y, z,
                   bx2, by2, bz2,
                   cx, cy, cz,
                   ax2, ay2, az2, coords);

      HilbertSfc3D(x+bx2, y+by2, z+bz2,
                   ax, ay, az,
                   bx-bx2, by-by2, bz-bz2,
                   cx, cy, cz, coords);

      HilbertSfc3D(x+(ax-dax)+(bx2-dbx),
                   y+(ay-day)+(by2-dby),
                   z+(az-daz)+(bz2-dbz),
                   -bx2, -by2, -bz2,
                   cx, cy, cz,
                   -(ax-ax2), -(ay-ay2), -(az-az2), coords);
   }
   // do not split in h
   else if (3*d > 4*h)
   {
      HilbertSfc3D(x, y, z,
                   cx2, cy2, cz2,
                   ax2, ay2, az2,
                   bx, by, bz, coords);

      HilbertSfc3D(x+cx2, y+cy2, z+cz2,
                   ax, ay, az,
                   bx, by, bz,
                   cx-cx2, cy-cy2, cz-cz2, coords);

      HilbertSfc3D(x+(ax-dax)+(cx2-dcx),
                   y+(ay-day)+(cy2-dcy),
                   z+(az-daz)+(cz2-dcz),
                   -cx2, -cy2, -cz2,
                   -(ax-ax2), -(ay-ay2), -(az-az2),
                   bx, by, bz, coords);
   }
   // regular case, split in all w/h/d
   else
   {
      HilbertSfc3D(x, y, z,
                   bx2, by2, bz2,
                   cx2, cy2, cz2,
                   ax2, ay2, az2, coords);

      HilbertSfc3D(x+bx2, y+by2, z+bz2,
                   cx, cy, cz,
                   ax2, ay2, az2,
                   bx-bx2, by-by2, bz-bz2, coords);

      HilbertSfc3D(x+(bx2-dbx)+(cx-dcx),
                   y+(by2-dby)+(cy-dcy),
                   z+(bz2-dbz)+(cz-dcz),
                   ax, ay, az,
                   -bx2, -by2, -bz2,
                   -(cx-cx2), -(cy-cy2), -(cz-cz2), coords);

      HilbertSfc3D(x+(ax-dax)+bx2+(cx-dcx),
                   y+(ay-day)+by2+(cy-dcy),
                   z+(az-daz)+bz2+(cz-dcz),
                   -cx, -cy, -cz,
                   -(ax-ax2), -(ay-ay2), -(az-az2),
                   bx-bx2, by-by2, bz-bz2, coords);

      HilbertSfc3D(x+(ax-dax)+(bx2-dbx),
                   y+(ay-day)+(by2-dby),
                   z+(az-daz)+(bz2-dbz),
                   -bx2, -by2, -bz2,
                   cx2, cy2, cz2,
                   -(ax-ax2), -(ay-ay2), -(az-az2), coords);
   }
}

void GridSfcOrdering2D(int64_t width, int64_t height, Array<int64_t> &coords)
{
   coords.SetSize(0);
   coords.Reserve(2*width*height);

   if (width >= height)
   {
      HilbertSfc2D(0, 0, width, 0, 0, height, coords);
   }
   else
   {
      HilbertSfc2D(0, 0, 0, height, width, 0, coords);
   }
}

void GridSfcOrdering3D(int64_t width, int64_t height, int64_t depth,
                       Array<int64_t> &coords)
{
   coords.SetSize(0);
   coords.Reserve(3*width*height*depth);

   if (width >= height && width >= depth)
   {
      HilbertSfc3D(0, 0, 0,
                   width, 0, 0,
                   0, height, 0,
                   0, 0, depth, coords);
   }
   else if (height >= width && height >= depth)
   {
      HilbertSfc3D(0, 0, 0,
                   0, height, 0,
                   width, 0, 0,
                   0, 0, depth, coords);
   }
   else // depth >= width && depth >= height
   {
      HilbertSfc3D(0, 0, 0,
                   0, 0, depth,
                   width, 0, 0,
                   0, height, 0, coords);
   }
}

Mesh::FaceInformation Mesh::GetFaceInformation(int64_t f) const
{
   FaceInformation face;
   int64_t e1, e2;
   int64_t inf1, inf2;
   int64_t ncface;
   GetFaceElements(f, &e1, &e2);
   GetFaceInfos(f, &inf1, &inf2, &ncface);
   face.element[0].index = e1;
   face.element[0].location = ElementLocation::Local;
   face.element[0].orientation = inf1%64;
   face.element[0].local_face_id = inf1/64;
   face.element[1].local_face_id = inf2/64;
   face.ncface = ncface;
   face.point_matrix = nullptr;
   // The following figures out face.location, face.conformity,
   // face.element[1].index, and face.element[1].orientation.
   if (f < GetNumFaces()) // Non-ghost face
   {
      if (e2>=0)
      {
         if (ncface==-1)
         {
            face.tag = FaceInfoTag::LocalConforming;
            face.topology = FaceTopology::Conforming;
            face.element[1].location = ElementLocation::Local;
            face.element[0].conformity = ElementConformity::Coincident;
            face.element[1].conformity = ElementConformity::Coincident;
            face.element[1].index = e2;
            face.element[1].orientation = inf2%64;
         }
         else // ncface >= 0
         {
            face.tag = FaceInfoTag::LocalSlaveNonconforming;
            face.topology = FaceTopology::Nonconforming;
            face.element[1].location = ElementLocation::Local;
            face.element[0].conformity = ElementConformity::Coincident;
            face.element[1].conformity = ElementConformity::Superset;
            face.element[1].index = e2;
            MFEM_ASSERT(inf2%64==0, "unexpected slave face orientation.");
            face.element[1].orientation = inf2%64;
            face.point_matrix = nc_faces_info[ncface].PointMatrix;
         }
      }
      else // e2<0
      {
         if (ncface==-1)
         {
            if (inf2<0)
            {
               face.tag = FaceInfoTag::Boundary;
               face.topology = FaceTopology::Boundary;
               face.element[1].location = ElementLocation::NA;
               face.element[0].conformity = ElementConformity::Coincident;
               face.element[1].conformity = ElementConformity::NA;
               face.element[1].index = -1;
               face.element[1].orientation = -1;
            }
            else // inf2 >= 0
            {
               face.tag = FaceInfoTag::SharedConforming;
               face.topology = FaceTopology::Conforming;
               face.element[0].conformity = ElementConformity::Coincident;
               face.element[1].conformity = ElementConformity::Coincident;
               face.element[1].location = ElementLocation::FaceNbr;
               face.element[1].index = -1 - e2;
               face.element[1].orientation = inf2%64;
            }
         }
         else // ncface >= 0
         {
            if (inf2 < 0)
            {
               face.tag = FaceInfoTag::MasterNonconforming;
               face.topology = FaceTopology::Nonconforming;
               face.element[1].location = ElementLocation::NA;
               face.element[0].conformity = ElementConformity::Coincident;
               face.element[1].conformity = ElementConformity::Subset;
               face.element[1].index = -1;
               face.element[1].orientation = -1;
            }
            else
            {
               face.tag = FaceInfoTag::SharedSlaveNonconforming;
               face.topology = FaceTopology::Nonconforming;
               face.element[1].location = ElementLocation::FaceNbr;
               face.element[0].conformity = ElementConformity::Coincident;
               face.element[1].conformity = ElementConformity::Superset;
               face.element[1].index = -1 - e2;
               face.element[1].orientation = inf2%64;
            }
            face.point_matrix = nc_faces_info[ncface].PointMatrix;
         }
      }
   }
   else // Ghost face
   {
      if (e1==-1)
      {
         face.tag = FaceInfoTag::GhostMaster;
         face.topology = FaceTopology::NA;
         face.element[1].location = ElementLocation::NA;
         face.element[0].conformity = ElementConformity::NA;
         face.element[1].conformity = ElementConformity::NA;
         face.element[1].index = -1;
         face.element[1].orientation = -1;
      }
      else
      {
         face.tag = FaceInfoTag::GhostSlave;
         face.topology = FaceTopology::Nonconforming;
         face.element[1].location = ElementLocation::FaceNbr;
         face.element[0].conformity = ElementConformity::Superset;
         face.element[1].conformity = ElementConformity::Coincident;
         face.element[1].index = -1 - e2;
         face.element[1].orientation = inf2%64;
         face.point_matrix = nc_faces_info[ncface].PointMatrix;
      }
   }
   return face;
}

Mesh::FaceInformation::operator Mesh::FaceInfo() const
{
   FaceInfo res {-1, -1, -1, -1, -1};
   switch (tag)
   {
      case FaceInfoTag::LocalConforming:
         res.Elem1No = element[0].index;
         res.Elem2No = element[1].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         res.Elem2Inf = element[1].orientation + element[1].local_face_id*64;
         res.NCFace = ncface;
         break;
      case FaceInfoTag::LocalSlaveNonconforming:
         res.Elem1No = element[0].index;
         res.Elem2No = element[1].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         res.Elem2Inf = element[1].orientation + element[1].local_face_id*64;
         res.NCFace = ncface;
         break;
      case FaceInfoTag::Boundary:
         res.Elem1No = element[0].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         break;
      case FaceInfoTag::SharedConforming:
         res.Elem1No = element[0].index;
         res.Elem2No = -1 - element[1].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         res.Elem2Inf = element[1].orientation + element[1].local_face_id*64;
         break;
      case FaceInfoTag::MasterNonconforming:
         res.Elem1No = element[0].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         break;
      case FaceInfoTag::SharedSlaveNonconforming:
         res.Elem1No = element[0].index;
         res.Elem2No = -1 - element[1].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         res.Elem2Inf = element[1].orientation + element[1].local_face_id*64;
         break;
      case FaceInfoTag::GhostMaster:
         break;
      case FaceInfoTag::GhostSlave:
         res.Elem1No = element[0].index;
         res.Elem2No = -1 - element[1].index;
         res.Elem1Inf = element[0].orientation + element[0].local_face_id*64;
         res.Elem2Inf = element[1].orientation + element[1].local_face_id*64;
         break;
   }
   return res;
}

std::ostream &operator<<(std::ostream &os, const Mesh::FaceInformation& info)
{
   os << "face topology=";
   switch (info.topology)
   {
      case Mesh::FaceTopology::Boundary:
         os << "Boundary";
         break;
      case Mesh::FaceTopology::Conforming:
         os << "Conforming";
         break;
      case Mesh::FaceTopology::Nonconforming:
         os << "Non-conforming";
         break;
      case Mesh::FaceTopology::NA:
         os << "NA";
         break;
   }
   os << '\n';
   os << "element[0].location=";
   switch (info.element[0].location)
   {
      case Mesh::ElementLocation::Local:
         os << "Local";
         break;
      case Mesh::ElementLocation::FaceNbr:
         os << "FaceNbr";
         break;
      case Mesh::ElementLocation::NA:
         os << "NA";
         break;
   }
   os << '\n';
   os << "element[1].location=";
   switch (info.element[1].location)
   {
      case Mesh::ElementLocation::Local:
         os << "Local";
         break;
      case Mesh::ElementLocation::FaceNbr:
         os << "FaceNbr";
         break;
      case Mesh::ElementLocation::NA:
         os << "NA";
         break;
   }
   os << '\n';
   os << "element[0].conformity=";
   switch (info.element[0].conformity)
   {
      case Mesh::ElementConformity::Coincident:
         os << "Coincident";
         break;
      case Mesh::ElementConformity::Superset:
         os << "Superset";
         break;
      case Mesh::ElementConformity::Subset:
         os << "Subset";
         break;
      case Mesh::ElementConformity::NA:
         os << "NA";
         break;
   }
   os << '\n';
   os << "element[1].conformity=";
   switch (info.element[1].conformity)
   {
      case Mesh::ElementConformity::Coincident:
         os << "Coincident";
         break;
      case Mesh::ElementConformity::Superset:
         os << "Superset";
         break;
      case Mesh::ElementConformity::Subset:
         os << "Subset";
         break;
      case Mesh::ElementConformity::NA:
         os << "NA";
         break;
   }
   os << '\n';
   os << "element[0].index=" << info.element[0].index << '\n'
      << "element[1].index=" << info.element[1].index << '\n'
      << "element[0].local_face_id=" << info.element[0].local_face_id << '\n'
      << "element[1].local_face_id=" << info.element[1].local_face_id << '\n'
      << "element[0].orientation=" << info.element[0].orientation << '\n'
      << "element[1].orientation=" << info.element[1].orientation << '\n'
      << "ncface=" << info.ncface << std::endl;
   return os;
}

void Mesh::GetFaceElements(int64_t Face, int64_t *Elem1, int64_t *Elem2) const
{
   *Elem1 = faces_info[Face].Elem1No;
   *Elem2 = faces_info[Face].Elem2No;
}

void Mesh::GetFaceInfos(int64_t Face, int64_t *Inf1, int64_t *Inf2) const
{
   *Inf1 = faces_info[Face].Elem1Inf;
   *Inf2 = faces_info[Face].Elem2Inf;
}

void Mesh::GetFaceInfos(int64_t Face, int64_t *Inf1, int64_t *Inf2,
                        int64_t *NCFace) const
{
   *Inf1   = faces_info[Face].Elem1Inf;
   *Inf2   = faces_info[Face].Elem2Inf;
   *NCFace = faces_info[Face].NCFace;
}

Geometry::Type Mesh::GetFaceGeometry(int64_t Face) const
{
   switch (Dim)
   {
      case 1: return Geometry::POINT;
      case 2: return Geometry::SEGMENT;
      case 3:
         if (Face < NumOfFaces) // local (non-ghost) face
         {
            return faces[Face]->GetGeometryType();
         }
         // ghost face
         const int64_t nc_face_id = faces_info[Face].NCFace;

         MFEM_ASSERT(nc_face_id >= 0, "parent ghost faces are not supported");
         return faces[nc_faces_info[nc_face_id].MasterFace]->GetGeometryType();
   }
   return Geometry::INVALID;
}

Element::Type Mesh::GetFaceElementType(int64_t Face) const
{
   return (Dim == 1) ? Element::POINT : faces[Face]->GetType();
}

Array<int64_t> Mesh::GetFaceToBdrElMap() const
{
   Array<int64_t> face_to_be(Dim == 2 ? NumOfEdges : NumOfFaces);
   face_to_be = -1;
   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      face_to_be[GetBdrElementFaceIndex(i)] = i;
   }
   return face_to_be;
}

void Mesh::Init()
{
   // in order of declaration:
   Dim = spaceDim = 0;
   NumOfVertices = -1;
   NumOfElements = NumOfBdrElements = 0;
   NumOfEdges = NumOfFaces = 0;
   nbInteriorFaces = -1;
   nbBoundaryFaces = -1;
   meshgen = mesh_geoms = 0;
}

void Mesh::InitTables()
{
   el_to_edge =
      el_to_face = el_to_el = bel_to_edge = face_edge = edge_vertex = NULL;
   face_to_elem = NULL;
}

void Mesh::SetEmpty()
{
   Init();
   InitTables();
}

void Mesh::DestroyTables()
{
   delete el_to_edge;
   delete el_to_face;
   delete el_to_el;

   if (Dim == 3)
   {
      delete bel_to_edge;
   }

   delete face_edge;
   delete edge_vertex;

   delete face_to_elem;
   face_to_elem = NULL;
}

void Mesh::DestroyPointers()
{
   for (int64_t i = 0; i < NumOfElements; i++)
   {
      FreeElement(elements[i]);
   }

   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      FreeElement(boundary[i]);
   }

   for (int64_t i = 0; i < faces.Size(); i++)
   {
      FreeElement(faces[i]);
   }

   DestroyTables();
}

void Mesh::Destroy()
{
   DestroyPointers();

   elements.DeleteAll();
   vertices.DeleteAll();
   boundary.DeleteAll();
   faces.DeleteAll();
   faces_info.DeleteAll();
   nc_faces_info.DeleteAll();
   be_to_face.DeleteAll();

   attributes.DeleteAll();
   bdr_attributes.DeleteAll();
}

void Mesh::ResetLazyData()
{
   delete el_to_el;     el_to_el = NULL;
   delete face_edge;    face_edge = NULL;
   delete face_to_elem;    face_to_elem = NULL;
   delete edge_vertex;  edge_vertex = NULL;
   nbInteriorFaces = -1;
   nbBoundaryFaces = -1;
}

void Mesh::SetAttributes()
{
   Array<int64_t> attribs;

   attribs.SetSize(GetNBE());
   for (int64_t i = 0; i < attribs.Size(); i++)
   {
      attribs[i] = GetBdrAttribute(i);
   }
   attribs.Sort();
   attribs.Unique();
   attribs.Copy(bdr_attributes);
   if (bdr_attributes.Size() > 0 && bdr_attributes[0] <= 0)
   {
      MFEM_WARNING("Non-positive attributes on the boundary!");
   }

   attribs.SetSize(GetNE());
   for (int64_t i = 0; i < attribs.Size(); i++)
   {
      attribs[i] = GetAttribute(i);
   }
   attribs.Sort();
   attribs.Unique();
   attribs.Copy(attributes);
   if (attributes.Size() > 0 && attributes[0] <= 0)
   {
      MFEM_WARNING("Non-positive attributes in the domain!");
   }
}

void Mesh::InitMesh(int64_t Dim_, int64_t spaceDim_, int64_t NVert,
                    int64_t NElem, int64_t NBdrElem)
{
   SetEmpty();

   Dim = Dim_;
   spaceDim = spaceDim_;

   NumOfVertices = 0;
   vertices.SetSize(NVert);  // just allocate space for vertices

   NumOfElements = 0;
   elements.SetSize(NElem);  // just allocate space for Element *

   NumOfBdrElements = 0;
   boundary.SetSize(NBdrElem);  // just allocate space for Element *
}

template<typename T>
static void CheckEnlarge(Array<T> &array, int64_t size)
{
   if (size >= array.Size()) { array.SetSize(size + 1); }
}

int64_t Mesh::AddVertex(real_t x, real_t y, real_t z)
{
   CheckEnlarge(vertices, NumOfVertices);
   real_t *v = vertices[NumOfVertices]();
   v[0] = x;
   v[1] = y;
   v[2] = z;
   return NumOfVertices++;
}

int64_t Mesh::AddVertex(const real_t *coords)
{
   CheckEnlarge(vertices, NumOfVertices);
   vertices[NumOfVertices].SetCoords(spaceDim, coords);
   return NumOfVertices++;
}

int64_t Mesh::AddVertex(const Vector &coords)
{
   MFEM_ASSERT(coords.Size() >= spaceDim,
               "invalid 'coords' size: " << coords.Size());
   return AddVertex(coords.GetData());
}

int64_t Mesh::AddVertexAtMeanCenter(const int64_t *vi, int64_t nverts,
                                    int64_t dim)
{
   Vector vii(dim);
   vii = 0.0;
   for (int64_t i = 0; i < nverts; i++)
   {
      real_t *vp = vertices[vi[i]]();
      for (int64_t j = 0; j < dim; j++)
      {
         vii(j) += vp[j];
      }
   }
   vii /= nverts;
   AddVertex(vii);
   return NumOfVertices;
}

int64_t Mesh::AddSegment(int64_t v1, int64_t v2, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Segment(v1, v2, attr);
   return NumOfElements++;
}

int64_t Mesh::AddSegment(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Segment(vi, attr);
   return NumOfElements++;
}

int64_t Mesh::AddTriangle(int64_t v1, int64_t v2, int64_t v3, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Triangle(v1, v2, v3, attr);
   return NumOfElements++;
}

int64_t Mesh::AddTriangle(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Triangle(vi, attr);
   return NumOfElements++;
}

int64_t Mesh::AddQuad(int64_t v1, int64_t v2, int64_t v3, int64_t v4,
                      int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Quadrilateral(v1, v2, v3, v4, attr);
   return NumOfElements++;
}

int64_t Mesh::AddQuad(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Quadrilateral(vi, attr);
   return NumOfElements++;
}

int64_t Mesh::AddTet(int64_t v1, int64_t v2, int64_t v3, int64_t v4,
                     int64_t attr)
{
   int64_t vi[4] = {v1, v2, v3, v4};
   return AddTet(vi, attr);
}

int64_t Mesh::AddTet(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Tetrahedron(vi, attr);
   return NumOfElements++;
}

int64_t Mesh::AddWedge(int64_t v1, int64_t v2, int64_t v3, int64_t v4,
                       int64_t v5, int64_t v6, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Wedge(v1, v2, v3, v4, v5, v6, attr);
   return NumOfElements++;
}

int64_t Mesh::AddWedge(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Wedge(vi, attr);
   return NumOfElements++;
}

int64_t Mesh::AddPyramid(int64_t v1, int64_t v2, int64_t v3, int64_t v4,
                         int64_t v5, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Pyramid(v1, v2, v3, v4, v5, attr);
   return NumOfElements++;
}

int64_t Mesh::AddPyramid(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Pyramid(vi, attr);
   return NumOfElements++;
}

int64_t Mesh::AddHex(int64_t v1, int64_t v2, int64_t v3, int64_t v4, int64_t v5,
                     int64_t v6, int64_t v7, int64_t v8,
                     int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] =
      new Hexahedron(v1, v2, v3, v4, v5, v6, v7, v8, attr);
   return NumOfElements++;
}

int64_t Mesh::AddHex(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = new Hexahedron(vi, attr);
   return NumOfElements++;
}

void Mesh::AddHexAsTets(const int64_t *vi, int64_t attr)
{
   static const int64_t hex_to_tet[6][4] =
   {
      { 0, 1, 2, 6 }, { 0, 5, 1, 6 }, { 0, 4, 5, 6 },
      { 0, 2, 3, 6 }, { 0, 3, 7, 6 }, { 0, 7, 4, 6 }
   };
   int64_t ti[4];

   for (int64_t i = 0; i < 6; i++)
   {
      for (int64_t j = 0; j < 4; j++)
      {
         ti[j] = vi[hex_to_tet[i][j]];
      }
      AddTet(ti, attr);
   }
}

void Mesh::AddHexAsWedges(const int64_t *vi, int64_t attr)
{
   static const int64_t hex_to_wdg[2][6] =
   {
      { 0, 1, 2, 4, 5, 6 }, { 0, 2, 3, 4, 6, 7 }
   };
   int64_t ti[6];

   for (int64_t i = 0; i < 2; i++)
   {
      for (int64_t j = 0; j < 6; j++)
      {
         ti[j] = vi[hex_to_wdg[i][j]];
      }
      AddWedge(ti, attr);
   }
}

void Mesh::AddHexAsPyramids(const int64_t *vi, int64_t attr)
{
   static const int64_t hex_to_pyr[6][5] =
   {
      { 0, 1, 2, 3, 8 }, { 0, 4, 5, 1, 8 }, { 1, 5, 6, 2, 8 },
      { 2, 6, 7, 3, 8 }, { 3, 7, 4, 0, 8 }, { 7, 6, 5, 4, 8 }
   };
   int64_t ti[5];

   for (int64_t i = 0; i < 6; i++)
   {
      for (int64_t j = 0; j < 5; j++)
      {
         ti[j] = vi[hex_to_pyr[i][j]];
      }
      AddPyramid(ti, attr);
   }
}

void Mesh::AddQuadAs4TrisWithPoints(int64_t *vi, int64_t attr)
{
   int64_t num_faces = 4;
   static const int64_t quad_to_tri[4][2] =
   {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}
   };

   int64_t elem_center_index = AddVertexAtMeanCenter(vi, 4, 2) - 1;

   int64_t ti[3];
   ti[2] = elem_center_index;
   for (int64_t i = 0; i < num_faces; i++)
   {
      for (int64_t j = 0; j < 2; j++)
      {
         ti[j] = vi[quad_to_tri[i][j]];
      }
      AddTri(ti, attr);
   }
}

void Mesh::AddQuadAs5QuadsWithPoints(int64_t *vi, int64_t attr)
{
   int64_t num_faces = 4;
   static const int64_t quad_faces[4][2] =
   {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}
   };

   Vector px(4), py(4);
   for (int64_t i = 0; i < 4; i++)
   {
      real_t *vp = vertices[vi[i]]();
      px(i) = vp[0];
      py(i) = vp[1];
   }

   int64_t vnew_index[4];
   real_t vnew[2];
   real_t r = 0.25, s = 0.25;
   vnew[0] = px(0)*(1-r)*(1-s) + px(1)*(r)*(1-s) + px(2)*r*s + px(3)*(1-r)*s;
   vnew[1] = py(0)*(1-r)*(1-s) + py(1)*(r)*(1-s) + py(2)*r*s + py(3)*(1-r)*s;
   AddVertex(vnew);
   vnew_index[0] = NumOfVertices-1;

   r = 0.75, s = 0.25;
   vnew[0] = px(0)*(1-r)*(1-s) + px(1)*(r)*(1-s) + px(2)*r*s + px(3)*(1-r)*s;
   vnew[1] = py(0)*(1-r)*(1-s) + py(1)*(r)*(1-s) + py(2)*r*s + py(3)*(1-r)*s;
   AddVertex(vnew);
   vnew_index[1] = NumOfVertices-1;

   r = 0.75, s = 0.75;
   vnew[0] = px(0)*(1-r)*(1-s) + px(1)*(r)*(1-s) + px(2)*r*s + px(3)*(1-r)*s;
   vnew[1] = py(0)*(1-r)*(1-s) + py(1)*(r)*(1-s) + py(2)*r*s + py(3)*(1-r)*s;
   AddVertex(vnew);
   vnew_index[2] = NumOfVertices-1;

   r = 0.25, s = 0.75;
   vnew[0] = px(0)*(1-r)*(1-s) + px(1)*(r)*(1-s) + px(2)*r*s + px(3)*(1-r)*s;
   vnew[1] = py(0)*(1-r)*(1-s) + py(1)*(r)*(1-s) + py(2)*r*s + py(3)*(1-r)*s;
   AddVertex(vnew);
   vnew_index[3] = NumOfVertices-1;

   static const int64_t quad_faces_new[4][2] =
   {
      { 1, 0}, { 2, 1}, { 3, 2}, { 0, 3}
   };

   int64_t ti[4];
   for (int64_t i = 0; i < num_faces; i++)
   {
      for (int64_t j = 0; j < 2; j++)
      {
         ti[j] = vi[quad_faces[i][j]];
         ti[j+2] = vnew_index[quad_faces_new[i][j]];
      }
      AddQuad(ti, attr);
   }
   AddQuad(vnew_index, attr);
}

void Mesh::AddHexAs24TetsWithPoints(int64_t *vi,
                                    std::map<std::array<int64_t, 4>, int64_t> &hex_face_verts,
                                    int64_t attr)
{
   auto get4arraysorted = [](Array<int64_t> v)
   {
      v.Sort();
      return std::array<int64_t, 4> {v[0], v[1], v[2], v[3]};
   };

   int64_t num_faces = 6;
   static const int64_t hex_to_tet[6][4] =
   {
      { 0, 1, 2, 3 }, { 1, 2, 6, 5 }, { 5, 4, 7, 6},
      { 0, 1, 5, 4 }, { 2, 3, 7, 6 }, { 0,3, 7, 4}
   };

   int64_t elem_center_index = AddVertexAtMeanCenter(vi, 8, 3) - 1;

   Array<int64_t> flist(4);

   // local vertex indices for each of the 4 edges of the face
   static const int64_t tet_face[4][2] =
   {
      {0, 1}, {1, 2}, {3, 2}, {3, 0}
   };

   for (int64_t i = 0; i < num_faces; i++)
   {
      for (int64_t j = 0; j < 4; j++)
      {
         flist[j] = vi[hex_to_tet[i][j]];
      }
      int64_t face_center_index;

      auto t = get4arraysorted(flist);
      auto it = hex_face_verts.find(t);
      if (it == hex_face_verts.end())
      {
         face_center_index = AddVertexAtMeanCenter(flist.GetData(),
                                                   flist.Size(), 3) - 1;
         hex_face_verts.insert({t, face_center_index});
      }
      else
      {
         face_center_index = it->second;
      }
      int64_t fti[4];
      fti[2] = face_center_index;
      fti[3] = elem_center_index;
      for (int64_t j = 0; j < 4; j++)
      {
         for (int64_t k = 0; k < 2; k++)
         {
            fti[k] = flist[tet_face[j][k]];
         }
         AddTet(fti, attr);
      }
   }
}

int64_t Mesh::AddElement(Element *elem)
{
   CheckEnlarge(elements, NumOfElements);
   elements[NumOfElements] = elem;
   return NumOfElements++;
}

int64_t Mesh::AddBdrElement(Element *elem)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = elem;
   return NumOfBdrElements++;
}

void Mesh::AddBdrElements(Array<Element *> &bdr_elems,
                          const Array<int64_t> &new_be_to_face)
{
   boundary.Reserve(boundary.Size() + bdr_elems.Size());
   MFEM_ASSERT(bdr_elems.Size() == new_be_to_face.Size(), "wrong size");
   for (int64_t i = 0; i < bdr_elems.Size(); i++)
   {
      AddBdrElement(bdr_elems[i]);
   }
   be_to_face.Append(new_be_to_face);
}

int64_t Mesh::AddBdrSegment(int64_t v1, int64_t v2, int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Segment(v1, v2, attr);
   return NumOfBdrElements++;
}

int64_t Mesh::AddBdrSegment(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Segment(vi, attr);
   return NumOfBdrElements++;
}

int64_t Mesh::AddBdrTriangle(int64_t v1, int64_t v2, int64_t v3, int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Triangle(v1, v2, v3, attr);
   return NumOfBdrElements++;
}

int64_t Mesh::AddBdrTriangle(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Triangle(vi, attr);
   return NumOfBdrElements++;
}

int64_t Mesh::AddBdrQuad(int64_t v1, int64_t v2, int64_t v3, int64_t v4,
                         int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Quadrilateral(v1, v2, v3, v4, attr);
   return NumOfBdrElements++;
}

int64_t Mesh::AddBdrQuad(const int64_t *vi, int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Quadrilateral(vi, attr);
   return NumOfBdrElements++;
}

void Mesh::AddBdrQuadAsTriangles(const int64_t *vi, int64_t attr)
{
   static const int64_t quad_to_tri[2][3] = { { 0, 1, 2 }, { 0, 2, 3 } };
   int64_t ti[3];

   for (int64_t i = 0; i < 2; i++)
   {
      for (int64_t j = 0; j < 3; j++)
      {
         ti[j] = vi[quad_to_tri[i][j]];
      }
      AddBdrTriangle(ti, attr);
   }
}

int64_t Mesh::AddBdrPoint(int64_t v, int64_t attr)
{
   CheckEnlarge(boundary, NumOfBdrElements);
   boundary[NumOfBdrElements] = new Point(&v, attr);
   return NumOfBdrElements++;
}

void Mesh::GenerateBoundaryElements()
{
   for (auto &b : boundary)
   {
      FreeElement(b);
   }

   if (Dim == 3)
   {
      delete bel_to_edge;
      bel_to_edge = NULL;
   }

   // count the 'NumOfBdrElements'
   NumOfBdrElements = 0;
   for (const auto &fi : faces_info)
   {
      if (fi.Elem2No < 0) { ++NumOfBdrElements; }
   }

   // Add the boundary elements
   boundary.SetSize(NumOfBdrElements);
   be_to_face.SetSize(NumOfBdrElements);
   for (int64_t i = 0, j = 0; i < faces_info.Size(); i++)
   {
      if (faces_info[i].Elem2No < 0)
      {
         boundary[j] = faces[i]->Duplicate(this);
         be_to_face[j++] = i;
      }
   }

   // Note: in 3D, 'bel_to_edge' is destroyed but it's not updated.
}

void Mesh::FinalizeCheck()
{
   MFEM_VERIFY(vertices.Size() == NumOfVertices ||
               vertices.Size() == 0,
               "incorrect number of vertices: preallocated: " << vertices.Size()
               << ", actually added: " << NumOfVertices);
   MFEM_VERIFY(elements.Size() == NumOfElements,
               "incorrect number of elements: preallocated: " << elements.Size()
               << ", actually added: " << NumOfElements);
   MFEM_VERIFY(boundary.Size() == NumOfBdrElements,
               "incorrect number of boundary elements: preallocated: "
               << boundary.Size() << ", actually added: " << NumOfBdrElements);
}

void Mesh::FinalizeTriMesh(int64_t generate_edges, int64_t refine,
                           bool fix_orientation)
{
   FinalizeCheck();

   if (refine)
   {
      MarkTriMeshForRefinement();
   }

   if (generate_edges)
   {
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
      GenerateFaces();
   }
   else
   {
      NumOfEdges = 0;
   }

   NumOfFaces = 0;

   SetAttributes();

   SetMeshGen();
}

void Mesh::FinalizeQuadMesh(int64_t generate_edges, int64_t refine,
                            bool fix_orientation)
{
   FinalizeCheck();

   if (generate_edges)
   {
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
      GenerateFaces();
   }
   else
   {
      NumOfEdges = 0;
   }

   NumOfFaces = 0;

   SetAttributes();

   SetMeshGen();
}

struct HilbertCmp
{
   int64_t coord;
   bool dir;
   const Array<real_t> &points;
   real_t mid;

   HilbertCmp(int64_t coord, bool dir, const Array<real_t> &points, real_t mid)
      : coord(coord), dir(dir), points(points), mid(mid) {}

   bool operator()(int64_t i) const
   {
      return (points[3*i + coord] < mid) != dir;
   }
};

static void HilbertSort2D(int64_t coord1, // major coordinate to sort points by
                          bool dir1,  // sort coord1 ascending/descending?
                          bool dir2,  // sort coord2 ascending/descending?
                          const Array<real_t> &points, int64_t *beg, int64_t *end,
                          real_t xmin, real_t ymin, real_t xmax, real_t ymax)
{
   if (end - beg <= 1) { return; }

   real_t xmid = (xmin + xmax)*0.5;
   real_t ymid = (ymin + ymax)*0.5;

   int64_t coord2 = (coord1 + 1) % 2; // the 'other' coordinate

   // sort (partition) points into four quadrants
   int64_t *p0 = beg, *p4 = end;
   int64_t *p2 = std::partition(p0, p4, HilbertCmp(coord1,  dir1, points, xmid));
   int64_t *p1 = std::partition(p0, p2, HilbertCmp(coord2,  dir2, points, ymid));
   int64_t *p3 = std::partition(p2, p4, HilbertCmp(coord2, !dir2, points, ymid));

   if (p1 != p4)
   {
      HilbertSort2D(coord2, dir2, dir1, points, p0, p1,
                    ymin, xmin, ymid, xmid);
   }
   if (p1 != p0 || p2 != p4)
   {
      HilbertSort2D(coord1, dir1, dir2, points, p1, p2,
                    xmin, ymid, xmid, ymax);
   }
   if (p2 != p0 || p3 != p4)
   {
      HilbertSort2D(coord1, dir1, dir2, points, p2, p3,
                    xmid, ymid, xmax, ymax);
   }
   if (p3 != p0)
   {
      HilbertSort2D(coord2, !dir2, !dir1, points, p3, p4,
                    ymid, xmax, ymin, xmid);
   }
}

static void HilbertSort3D(int64_t coord1, bool dir1, bool dir2, bool dir3,
                          const Array<real_t> &points, int64_t *beg, int64_t *end,
                          real_t xmin, real_t ymin, real_t zmin,
                          real_t xmax, real_t ymax, real_t zmax)
{
   if (end - beg <= 1) { return; }

   real_t xmid = (xmin + xmax)*0.5;
   real_t ymid = (ymin + ymax)*0.5;
   real_t zmid = (zmin + zmax)*0.5;

   int64_t coord2 = (coord1 + 1) % 3;
   int64_t coord3 = (coord1 + 2) % 3;

   // sort (partition) points into eight octants
   int64_t *p0 = beg, *p8 = end;
   int64_t *p4 = std::partition(p0, p8, HilbertCmp(coord1,  dir1, points, xmid));
   int64_t *p2 = std::partition(p0, p4, HilbertCmp(coord2,  dir2, points, ymid));
   int64_t *p6 = std::partition(p4, p8, HilbertCmp(coord2, !dir2, points, ymid));
   int64_t *p1 = std::partition(p0, p2, HilbertCmp(coord3,  dir3, points, zmid));
   int64_t *p3 = std::partition(p2, p4, HilbertCmp(coord3, !dir3, points, zmid));
   int64_t *p5 = std::partition(p4, p6, HilbertCmp(coord3,  dir3, points, zmid));
   int64_t *p7 = std::partition(p6, p8, HilbertCmp(coord3, !dir3, points, zmid));

   if (p1 != p8)
   {
      HilbertSort3D(coord3, dir3, dir1, dir2, points, p0, p1,
                    zmin, xmin, ymin, zmid, xmid, ymid);
   }
   if (p1 != p0 || p2 != p8)
   {
      HilbertSort3D(coord2, dir2, dir3, dir1, points, p1, p2,
                    ymin, zmid, xmin, ymid, zmax, xmid);
   }
   if (p2 != p0 || p3 != p8)
   {
      HilbertSort3D(coord2, dir2, dir3, dir1, points, p2, p3,
                    ymid, zmid, xmin, ymax, zmax, xmid);
   }
   if (p3 != p0 || p4 != p8)
   {
      HilbertSort3D(coord1, dir1, !dir2, !dir3, points, p3, p4,
                    xmin, ymax, zmid, xmid, ymid, zmin);
   }
   if (p4 != p0 || p5 != p8)
   {
      HilbertSort3D(coord1, dir1, !dir2, !dir3, points, p4, p5,
                    xmid, ymax, zmid, xmax, ymid, zmin);
   }
   if (p5 != p0 || p6 != p8)
   {
      HilbertSort3D(coord2, !dir2, dir3, !dir1, points, p5, p6,
                    ymax, zmid, xmax, ymid, zmax, xmid);
   }
   if (p6 != p0 || p7 != p8)
   {
      HilbertSort3D(coord2, !dir2, dir3, !dir1, points, p6, p7,
                    ymid, zmid, xmax, ymin, zmax, xmid);
   }
   if (p7 != p0)
   {
      HilbertSort3D(coord3, !dir3, !dir1, dir2, points, p7, p8,
                    zmid, xmax, ymin, zmin, xmid, ymid);
   }
}

void Mesh::ReorderElements(const Array<int64_t> &ordering,
                           bool reorder_vertices)
{
   MFEM_VERIFY(ordering.Size() == GetNE(), "invalid reordering array.")

   // Data members that need to be updated:

   // - elements   - reorder of the pointers and the vertex ids if reordering
   //                the vertices
   // - vertices   - if reordering the vertices
   // - boundary   - update the vertex ids, if reordering the vertices
   // - faces      - regenerate
   // - faces_info - regenerate

   // Deleted by DeleteTables():
   // - el_to_edge  - rebuild in 2D and 3D only
   // - el_to_face  - rebuild in 3D only
   // - bel_to_edge - rebuild in 3D only
   // - el_to_el    - no need to rebuild
   // - face_edge   - no need to rebuild
   // - edge_vertex - no need to rebuild
   // - geom_factors - no need to rebuild

   // - be_to_face

   // - Nodes

   // Save the locations of the Nodes so we can rebuild them later
   Array<Vector*> old_elem_node_vals;

   // Get the newly ordered elements
   Array<Element *> new_elements(GetNE());
   for (int64_t old_elid = 0; old_elid < ordering.Size(); ++old_elid)
   {
      int64_t new_elid = ordering[old_elid];
      new_elements[new_elid] = elements[old_elid];
   }
   mfem::Swap(elements, new_elements);
   new_elements.DeleteAll();

   if (reorder_vertices)
   {
      // Get the new vertex ordering permutation vectors and fill the new
      // vertices
      Array<int64_t> vertex_ordering(GetNV());
      vertex_ordering = -1;
      Array<Vertex> new_vertices(GetNV());
      int64_t new_vertex_ind = 0;
      for (int64_t new_elid = 0; new_elid < GetNE(); ++new_elid)
      {
         int64_t *elem_vert = elements[new_elid]->GetVertices();
         int64_t nv = elements[new_elid]->GetNVertices();
         for (int64_t vi = 0; vi < nv; ++vi)
         {
            int64_t old_vertex_ind = elem_vert[vi];
            if (vertex_ordering[old_vertex_ind] == -1)
            {
               vertex_ordering[old_vertex_ind] = new_vertex_ind;
               new_vertices[new_vertex_ind] = vertices[old_vertex_ind];
               new_vertex_ind++;
            }
         }
      }
      mfem::Swap(vertices, new_vertices);
      new_vertices.DeleteAll();

      // Replace the vertex ids in the elements with the reordered vertex
      // numbers
      for (int64_t new_elid = 0; new_elid < GetNE(); ++new_elid)
      {
         int64_t *elem_vert = elements[new_elid]->GetVertices();
         int64_t nv = elements[new_elid]->GetNVertices();
         for (int64_t vi = 0; vi < nv; ++vi)
         {
            elem_vert[vi] = vertex_ordering[elem_vert[vi]];
         }
      }

      // Replace the vertex ids in the boundary with reordered vertex numbers
      for (int64_t belid = 0; belid < GetNBE(); ++belid)
      {
         int64_t *be_vert = boundary[belid]->GetVertices();
         int64_t nv = boundary[belid]->GetNVertices();
         for (int64_t vi = 0; vi < nv; ++vi)
         {
            be_vert[vi] = vertex_ordering[be_vert[vi]];
         }
      }
   }

   // Destroy tables that need to be rebuild
   DeleteTables();

   if (Dim > 1)
   {
      // generate el_to_edge, be_to_face (2D), bel_to_edge (3D)
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   }
   if (Dim > 2)
   {
      // generate el_to_face, be_to_face
      GetElementToFaceTable();
   }
   // Update faces and faces_info
   GenerateFaces();
}


void Mesh::MarkForRefinement()
{
   if (meshgen & 1)
   {
      if (Dim == 2)
      {
         MarkTriMeshForRefinement();
      }
      else if (Dim == 3)
      {
         DSTable v_to_v(NumOfVertices);
         GetVertexToVertexTable(v_to_v);
         MarkTetMeshForRefinement(v_to_v);
      }
   }
}

void Mesh::MarkTriMeshForRefinement()
{
   // Mark the longest triangle edge by rotating the indices so that
   // vertex 0 - vertex 1 is the longest edge in the triangle.
   DenseMatrix pmat;
   for (int64_t i = 0; i < NumOfElements; i++)
   {
      if (elements[i]->GetType() == Element::TRIANGLE)
      {
         GetPointMatrix(i, pmat);
         static_cast<Triangle*>(elements[i])->MarkEdge(pmat);
      }
   }
}

void Mesh::GetEdgeOrdering(const DSTable &v_to_v, Array<int64_t> &order)
{
   NumOfEdges = v_to_v.NumberOfEntries();
   order.SetSize(NumOfEdges);
   Array<Pair<real_t, int64_t> > length_idx(NumOfEdges);

   for (int64_t i = 0; i < NumOfVertices; i++)
   {
      for (DSTable::RowIterator it(v_to_v, i); !it; ++it)
      {
         int64_t j = it.Index();
         length_idx[j].one = GetLength(i, it.Column());
         length_idx[j].two = j;
      }
   }

   // Sort by increasing edge-length.
   length_idx.Sort();

   for (int64_t i = 0; i < NumOfEdges; i++)
   {
      order[length_idx[i].two] = i;
   }
}

void Mesh::MarkTetMeshForRefinement(const DSTable &v_to_v)
{
   // Mark the longest tetrahedral edge by rotating the indices so that
   // vertex 0 - vertex 1 is the longest edge in the element.
   Array<int64_t> order;
   GetEdgeOrdering(v_to_v, order);

   for (int64_t i = 0; i < NumOfElements; i++)
   {
      if (elements[i]->GetType() == Element::TETRAHEDRON)
      {
         elements[i]->MarkEdge(v_to_v, order);
      }
   }
   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      if (boundary[i]->GetType() == Element::TRIANGLE)
      {
         boundary[i]->MarkEdge(v_to_v, order);
      }
   }
}

void Mesh::FinalizeTetMesh(int64_t generate_edges, int64_t refine,
                           bool fix_orientation)
{
   FinalizeCheck();

   if (!HasBoundaryElements())
   {
      GetElementToFaceTable();
      GenerateFaces();
      GenerateBoundaryElements();
   }

   if (refine)
   {
      DSTable v_to_v(NumOfVertices);
      GetVertexToVertexTable(v_to_v);
      MarkTetMeshForRefinement(v_to_v);
   }

   GetElementToFaceTable();
   GenerateFaces();

   if (generate_edges == 1)
   {
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   }
   else
   {
      el_to_edge = NULL;  // Not really necessary -- InitTables was called
      bel_to_edge = NULL;
      NumOfEdges = 0;
   }

   SetAttributes();

   SetMeshGen();
}

void Mesh::FinalizeWedgeMesh(int64_t generate_edges, int64_t refine,
                             bool fix_orientation)
{
   FinalizeCheck();

   if (!HasBoundaryElements())
   {
      GetElementToFaceTable();
      GenerateFaces();
      GenerateBoundaryElements();
   }

   GetElementToFaceTable();
   GenerateFaces();

   if (generate_edges == 1)
   {
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   }
   else
   {
      el_to_edge = NULL;  // Not really necessary -- InitTables was called
      bel_to_edge = NULL;
      NumOfEdges = 0;
   }

   SetAttributes();

   SetMeshGen();
}

void Mesh::FinalizeHexMesh(int64_t generate_edges, int64_t refine,
                           bool fix_orientation)
{
   FinalizeCheck();

   GetElementToFaceTable();
   GenerateFaces();

   if (!HasBoundaryElements())
   {
      GenerateBoundaryElements();
   }

   if (generate_edges)
   {
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   }
   else
   {
      NumOfEdges = 0;
   }

   SetAttributes();

   SetMeshGen();
}

void Mesh::FinalizeMesh(int64_t refine, bool fix_orientation)
{
   FinalizeTopology();

   Finalize(refine, fix_orientation);
}

void Mesh::FinalizeTopology(bool generate_bdr)
{
   // Requirements: the following should be defined:
   //   1) Dim
   //   2) NumOfElements, elements
   //   3) NumOfBdrElements, boundary
   //   4) NumOfVertices
   // Optional:
   //   2) ncmesh may be defined
   //   3) el_to_edge may be allocated (it will be re-computed)

   FinalizeCheck();
   bool generate_edges = true;

   if (spaceDim == 0) { spaceDim = Dim; }

   // set the mesh type: 'meshgen', ...
   SetMeshGen();

   // generate the faces
   if (Dim > 2)
   {
      GetElementToFaceTable();
      GenerateFaces();
      if (!HasBoundaryElements() && generate_bdr)
      {
         GenerateBoundaryElements();
         GetElementToFaceTable(); // update be_to_face
      }
   }
   else
   {
      NumOfFaces = 0;
   }

   // generate edges if requested
   if (Dim > 1 && generate_edges)
   {
      // el_to_edge may already be allocated (P2 VTK meshes)
      if (!el_to_edge) { el_to_edge = new Table; }
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
      if (Dim == 2)
      {
         GenerateFaces(); // 'Faces' in 2D refers to the edges
         if (!HasBoundaryElements() && generate_bdr)
         {
            GenerateBoundaryElements();
         }
      }
   }
   else
   {
      NumOfEdges = 0;
   }

   if (Dim == 1)
   {
      GenerateFaces();
      if (!HasBoundaryElements() && generate_bdr)
      {
         // be_to_face will be set inside GenerateBoundaryElements
         GenerateBoundaryElements();
      }
      else
      {
         be_to_face.SetSize(NumOfBdrElements);
         for (int64_t i = 0; i < NumOfBdrElements; ++i)
         {
            be_to_face[i] = boundary[i]->GetVertices()[0];
         }
      }
   }

   // generate the arrays 'attributes' and 'bdr_attributes'
   SetAttributes();
}

void Mesh::Finalize(bool refine, bool fix_orientation)
{
   // Requirements:
   //  1) FinalizeTopology() or equivalent was called
   //  2) if (Nodes == NULL), vertices must be defined
   //  3) if (Nodes != NULL), Nodes must be defined

   const bool check_orientation = true; // for regular elements, not boundary
   const bool may_change_topology =
      ( refine && (Dim > 1 && (meshgen & 1)) ) ||
      ( check_orientation && fix_orientation &&
        (Dim == 2 || (Dim == 3 && (meshgen & 1))) );

   if (refine)
   {
      MarkForRefinement();   // may change topology!
   }

   if (may_change_topology)
   {
      FinalizeTopology(); // Re-computes some data unnecessarily.

      // TODO: maybe introduce Mesh::NODE_REORDER operation and FESpace::
      // NodeReorderMatrix and do Nodes->Update() instead of DoNodeReorder?
   }

#ifdef MFEM_DEBUG
   // For non-orientable surfaces/manifolds, the check below will fail, so we
   // only perform it when Dim == spaceDim.
   if (Dim >= 2 && Dim == spaceDim)
   {
      const int64_t num_faces = GetNumFaces();
      for (int64_t i = 0; i < num_faces; i++)
      {
         MFEM_VERIFY(faces_info[i].Elem2No < 0 ||
                     faces_info[i].Elem2Inf%2 != 0, "Invalid mesh topology."
                     " Interior face with incompatible orientations.");
      }
   }
#endif
}

void Mesh::Make3D(int64_t nx, int64_t ny, int64_t nz, Element::Type type,
                  real_t sx, real_t sy, real_t sz, bool sfc_ordering)
{
   int64_t x, y, z;

   int64_t NVert, NElem, NBdrElem;

   NVert = (nx+1) * (ny+1) * (nz+1);
   NElem = nx * ny * nz;
   NBdrElem = 2*(nx*ny+nx*nz+ny*nz);
   if (type == Element::TETRAHEDRON)
   {
      NElem *= 6;
      NBdrElem *= 2;
   }
   else if (type == Element::WEDGE)
   {
      NElem *= 2;
      NBdrElem += 2*nx*ny;
   }
   else if (type == Element::PYRAMID)
   {
      NElem *= 6;
      NVert += nx * ny * nz;
   }

   InitMesh(3, 3, NVert, NElem, NBdrElem);

   real_t coord[3];
   int64_t ind[9];

   // Sets vertices and the corresponding coordinates
   for (z = 0; z <= nz; z++)
   {
      coord[2] = ((real_t) z / nz) * sz;
      for (y = 0; y <= ny; y++)
      {
         coord[1] = ((real_t) y / ny) * sy;
         for (x = 0; x <= nx; x++)
         {
            coord[0] = ((real_t) x / nx) * sx;
            AddVertex(coord);
         }
      }
   }
   if (type == Element::PYRAMID)
   {
      for (z = 0; z < nz; z++)
      {
         coord[2] = (((real_t) z + 0.5) / nz) * sz;
         for (y = 0; y < ny; y++)
         {
            coord[1] = (((real_t) y + 0.5) / ny) * sy;
            for (x = 0; x < nx; x++)
            {
               coord[0] = (((real_t) x + 0.5) / nx) * sx;
               AddVertex(coord);
            }
         }
      }
   }

#define VTX(XC, YC, ZC) ((XC)+((YC)+(ZC)*(ny+1))*(nx+1))
#define VTXP(XC, YC, ZC) ((nx+1)*(ny+1)*(nz+1)+(XC)+((YC)+(ZC)*ny)*nx)

   // Sets elements and the corresponding indices of vertices
   if (sfc_ordering && type == Element::HEXAHEDRON)
   {
      Array<int64_t> sfc;
      GridSfcOrdering3D(nx, ny, nz, sfc);
      MFEM_VERIFY(sfc.Size() == 3*nx*ny*nz, "");

      for (int64_t k = 0; k < nx*ny*nz; k++)
      {
         x = sfc[3*k + 0];
         y = sfc[3*k + 1];
         z = sfc[3*k + 2];

         // *INDENT-OFF*
         ind[0] = VTX(x  , y  , z  );
         ind[1] = VTX(x+1, y  , z  );
         ind[2] = VTX(x+1, y+1, z  );
         ind[3] = VTX(x  , y+1, z  );
         ind[4] = VTX(x  , y  , z+1);
         ind[5] = VTX(x+1, y  , z+1);
         ind[6] = VTX(x+1, y+1, z+1);
         ind[7] = VTX(x  , y+1, z+1);
         // *INDENT-ON*

         AddHex(ind, 1);
      }
   }
   else
   {
      for (z = 0; z < nz; z++)
      {
         for (y = 0; y < ny; y++)
         {
            for (x = 0; x < nx; x++)
            {
               // *INDENT-OFF*
               ind[0] = VTX(x  , y  , z  );
               ind[1] = VTX(x+1, y  , z  );
               ind[2] = VTX(x+1, y+1, z  );
               ind[3] = VTX(x  , y+1, z  );
               ind[4] = VTX(x  , y  , z+1);
               ind[5] = VTX(x+1, y  , z+1);
               ind[6] = VTX(x+1, y+1, z+1);
               ind[7] = VTX(  x, y+1, z+1);
               // *INDENT-ON*
               if (type == Element::TETRAHEDRON)
               {
                  AddHexAsTets(ind, 1);
               }
               else if (type == Element::WEDGE)
               {
                  AddHexAsWedges(ind, 1);
               }
               else if (type == Element::PYRAMID)
               {
                  ind[8] = VTXP(x, y, z);
                  AddHexAsPyramids(ind, 1);
               }
               else
               {
                  AddHex(ind, 1);
               }
            }
         }
      }
   }

   // Sets boundary elements and the corresponding indices of vertices
   // bottom, bdr. attribute 1
   for (y = 0; y < ny; y++)
   {
      for (x = 0; x < nx; x++)
      {
         // *INDENT-OFF*
         ind[0] = VTX(x  , y  , 0);
         ind[1] = VTX(x  , y+1, 0);
         ind[2] = VTX(x+1, y+1, 0);
         ind[3] = VTX(x+1, y  , 0);
         // *INDENT-ON*
         if (type == Element::TETRAHEDRON)
         {
            AddBdrQuadAsTriangles(ind, 1);
         }
         else if (type == Element::WEDGE)
         {
            AddBdrQuadAsTriangles(ind, 1);
         }
         else
         {
            AddBdrQuad(ind, 1);
         }
      }
   }
   // top, bdr. attribute 6
   for (y = 0; y < ny; y++)
   {
      for (x = 0; x < nx; x++)
      {
         // *INDENT-OFF*
         ind[0] = VTX(x  , y  , nz);
         ind[1] = VTX(x+1, y  , nz);
         ind[2] = VTX(x+1, y+1, nz);
         ind[3] = VTX(x  , y+1, nz);
         // *INDENT-ON*
         if (type == Element::TETRAHEDRON)
         {
            AddBdrQuadAsTriangles(ind, 6);
         }
         else if (type == Element::WEDGE)
         {
            AddBdrQuadAsTriangles(ind, 6);
         }
         else
         {
            AddBdrQuad(ind, 6);
         }
      }
   }
   // left, bdr. attribute 5
   for (z = 0; z < nz; z++)
   {
      for (y = 0; y < ny; y++)
      {
         // *INDENT-OFF*
         ind[0] = VTX(0  , y  , z  );
         ind[1] = VTX(0  , y  , z+1);
         ind[2] = VTX(0  , y+1, z+1);
         ind[3] = VTX(0  , y+1, z  );
         // *INDENT-ON*
         if (type == Element::TETRAHEDRON)
         {
            AddBdrQuadAsTriangles(ind, 5);
         }
         else
         {
            AddBdrQuad(ind, 5);
         }
      }
   }
   // right, bdr. attribute 3
   for (z = 0; z < nz; z++)
   {
      for (y = 0; y < ny; y++)
      {
         // *INDENT-OFF*
         ind[0] = VTX(nx, y  , z  );
         ind[1] = VTX(nx, y+1, z  );
         ind[2] = VTX(nx, y+1, z+1);
         ind[3] = VTX(nx, y  , z+1);
         // *INDENT-ON*
         if (type == Element::TETRAHEDRON)
         {
            AddBdrQuadAsTriangles(ind, 3);
         }
         else
         {
            AddBdrQuad(ind, 3);
         }
      }
   }
   // front, bdr. attribute 2
   for (x = 0; x < nx; x++)
   {
      for (z = 0; z < nz; z++)
      {
         // *INDENT-OFF*
         ind[0] = VTX(x  , 0, z  );
         ind[1] = VTX(x+1, 0, z  );
         ind[2] = VTX(x+1, 0, z+1);
         ind[3] = VTX(x  , 0, z+1);
         // *INDENT-ON*
         if (type == Element::TETRAHEDRON)
         {
            AddBdrQuadAsTriangles(ind, 2);
         }
         else
         {
            AddBdrQuad(ind, 2);
         }
      }
   }
   // back, bdr. attribute 4
   for (x = 0; x < nx; x++)
   {
      for (z = 0; z < nz; z++)
      {
         // *INDENT-OFF*
         ind[0] = VTX(x  , ny, z  );
         ind[1] = VTX(x  , ny, z+1);
         ind[2] = VTX(x+1, ny, z+1);
         ind[3] = VTX(x+1, ny, z  );
         // *INDENT-ON*
         if (type == Element::TETRAHEDRON)
         {
            AddBdrQuadAsTriangles(ind, 4);
         }
         else
         {
            AddBdrQuad(ind, 4);
         }
      }
   }

#undef VTX

#if 0
   ofstream test_stream("debug.mesh");
   Print(test_stream);
   test_stream.close();
#endif

   FinalizeTopology();

   // Finalize(...) can be called after this method, if needed
}


void Mesh::Make2D4TrisFromQuad(int64_t nx, int64_t ny, real_t sx, real_t sy)
{
   SetEmpty();

   Dim = 2;
   spaceDim = 2;

   NumOfVertices = (nx+1) * (ny+1);
   NumOfElements = nx * ny * 4;
   NumOfBdrElements =  (2 * nx + 2 * ny);
   vertices.SetSize(NumOfVertices);
   elements.SetSize(NumOfElements);
   boundary.SetSize(NumOfBdrElements);
   NumOfElements = 0;

   int64_t ind[4];

   // Sets vertices and the corresponding coordinates
   int64_t k = 0;
   for (real_t j = 0; j < ny+1; j++)
   {
      real_t cy = (j / ny) * sy;
      for (real_t i = 0; i < nx+1; i++)
      {
         real_t cx = (i / nx) * sx;
         vertices[k](0) = cx;
         vertices[k](1) = cy;
         k++;
      }
   }

   for (int64_t y = 0; y < ny; y++)
   {
      for (int64_t x = 0; x < nx; x++)
      {
         ind[0] = x + y*(nx+1);
         ind[1] = x + 1 +y*(nx+1);
         ind[2] = x + 1 + (y+1)*(nx+1);
         ind[3] = x + (y+1)*(nx+1);
         AddQuadAs4TrisWithPoints(ind, 1);
      }
   }

   int64_t m = (nx+1)*ny;
   for (int64_t i = 0; i < nx; i++)
   {
      boundary[i] = new Segment(i, i+1, 1);
      boundary[nx+i] = new Segment(m+i+1, m+i, 3);
   }
   m = nx+1;
   for (int64_t j = 0; j < ny; j++)
   {
      boundary[2*nx+j] = new Segment((j+1)*m, j*m, 4);
      boundary[2*nx+ny+j] = new Segment(j*m+nx, (j+1)*m+nx, 2);
   }

   SetMeshGen();

   el_to_edge = new Table;
   NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   GenerateFaces();

   NumOfFaces = 0;

   attributes.Append(1);
   bdr_attributes.Append(1); bdr_attributes.Append(2);
   bdr_attributes.Append(3); bdr_attributes.Append(4);

   FinalizeTopology();
}

void Mesh::Make2D5QuadsFromQuad(int64_t nx, int64_t ny,
                                real_t sx, real_t sy)
{
   SetEmpty();

   Dim = 2;
   spaceDim = 2;

   NumOfElements = nx * ny * 5;
   NumOfVertices = (nx+1) * (ny+1); //it will be enlarged later on
   NumOfBdrElements =  (2 * nx + 2 * ny);
   vertices.SetSize(NumOfVertices);
   elements.SetSize(NumOfElements);
   boundary.SetSize(NumOfBdrElements);
   NumOfElements = 0;

   int64_t ind[4];

   // Sets vertices and the corresponding coordinates
   int64_t k = 0;
   for (real_t j = 0; j < ny+1; j++)
   {
      real_t cy = (j / ny) * sy;
      for (real_t i = 0; i < nx+1; i++)
      {
         real_t cx = (i / nx) * sx;
         vertices[k](0) = cx;
         vertices[k](1) = cy;
         k++;
      }
   }

   for (int64_t y = 0; y < ny; y++)
   {
      for (int64_t x = 0; x < nx; x++)
      {
         ind[0] = x + y*(nx+1);
         ind[1] = x + 1 +y*(nx+1);
         ind[2] = x + 1 + (y+1)*(nx+1);
         ind[3] = x + (y+1)*(nx+1);
         AddQuadAs5QuadsWithPoints(ind, 1);
      }
   }

   int64_t m = (nx+1)*ny;
   for (int64_t i = 0; i < nx; i++)
   {
      boundary[i] = new Segment(i, i+1, 1);
      boundary[nx+i] = new Segment(m+i+1, m+i, 3);
   }
   m = nx+1;
   for (int64_t j = 0; j < ny; j++)
   {
      boundary[2*nx+j] = new Segment((j+1)*m, j*m, 4);
      boundary[2*nx+ny+j] = new Segment(j*m+nx, (j+1)*m+nx, 2);
   }

   SetMeshGen();

   el_to_edge = new Table;
   NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   GenerateFaces();

   NumOfFaces = 0;

   attributes.Append(1);
   bdr_attributes.Append(1); bdr_attributes.Append(2);
   bdr_attributes.Append(3); bdr_attributes.Append(4);

   FinalizeTopology();
}

void Mesh::Make3D24TetsFromHex(int64_t nx, int64_t ny, int64_t nz,
                               real_t sx, real_t sy, real_t sz)
{
   const int64_t NVert = (nx+1) * (ny+1) * (nz+1);
   const int64_t NElem = nx * ny * nz * 24;
   const int64_t NBdrElem = 2*(nx*ny+nx*nz+ny*nz)*4;

   InitMesh(3, 3, NVert, NElem, NBdrElem);

   real_t coord[3];

   // Sets vertices and the corresponding coordinates
   for (real_t z = 0; z <= nz; z++)
   {
      coord[2] = ( z / nz) * sz;
      for (real_t y = 0; y <= ny; y++)
      {
         coord[1] = (y / ny) * sy;
         for (real_t x = 0; x <= nx; x++)
         {
            coord[0] = (x / nx) * sx;
            AddVertex(coord);
         }
      }
   }

   std::map<std::array<int64_t, 4>, int64_t> hex_face_verts;
   auto VertexIndex = [nx, ny](int64_t xc, int64_t yc, int64_t zc)
   {
      return xc + (yc + zc*(ny+1))*(nx+1);
   };

   int64_t ind[9];
   for (int64_t z = 0; z < nz; z++)
   {
      for (int64_t y = 0; y < ny; y++)
      {
         for (int64_t x = 0; x < nx; x++)
         {
            // *INDENT-OFF*
            ind[0] = VertexIndex(x  , y  , z  );
            ind[1] = VertexIndex(x+1, y  , z  );
            ind[2] = VertexIndex(x+1, y+1, z  );
            ind[3] = VertexIndex(x  , y+1, z  );
            ind[4] = VertexIndex(x  , y  , z+1);
            ind[5] = VertexIndex(x+1, y  , z+1);
            ind[6] = VertexIndex(x+1, y+1, z+1);
            ind[7] = VertexIndex(  x, y+1, z+1);
            AddHexAs24TetsWithPoints(ind, hex_face_verts, 1);
         }
      }
   }

   hex_face_verts.clear();

   // Done adding Tets
   // Now figure out elements that are on the boundary
   GetElementToFaceTable(false);
   GenerateFaces();

   // Map to count number of tets sharing a face
   std::map<std::array<int64_t, 3>, int64_t> tet_face_count;
   // Map from tet face defined by three vertices to the local face number
   std::map<std::array<int64_t, 3>, int64_t> face_count_map;

   auto get3array = [](Array<int64_t> v)
   {
       v.Sort();
       return std::array<int64_t, 3>{v[0], v[1], v[2]};
   };

   Array<int64_t> el_faces;
   Array<int64_t> ori;
   Array<int64_t> vertidxs;
   for (int64_t i = 0; i < el_to_face->Size(); i++)
   {
       el_to_face->GetRow(i, el_faces);
       for (int64_t j = 0; j < el_faces.Size(); j++)
       {
           GetFaceVertices(el_faces[j], vertidxs);
           auto t = get3array(vertidxs);
           auto it = tet_face_count.find(t);
           if (it == tet_face_count.end())  //edge does not already exist
           {
               tet_face_count.insert({t, 1});
               face_count_map.insert({t, el_faces[j]});
           }
           else
           {
               it->second++; // increase edge count value by 1.
           }
       }
   }

   for (const auto &edge : tet_face_count)
   {
       if (edge.second == 1)  //if this only appears once, it is a boundary edge
       {
           int64_t facenum = (face_count_map.find(edge.first))->second;
           GetFaceVertices(facenum, vertidxs);
           AddBdrTriangle(vertidxs, 1);
       }
   }

#if 0
   ofstream test_stream("debug.mesh");
   Print(test_stream);
   test_stream.close();
#endif

   FinalizeTopology();
   // Finalize(...) can be called after this method, if needed
}

void Mesh::Make2D(int64_t nx, int64_t ny, Element::Type type,
                  real_t sx, real_t sy,
                  bool generate_edges, bool sfc_ordering)
{
   int64_t i, j, k;

   SetEmpty();

   Dim = spaceDim = 2;

   // Creates quadrilateral mesh
   if (type == Element::QUADRILATERAL)
   {
      NumOfVertices = (nx+1) * (ny+1);
      NumOfElements = nx * ny;
      NumOfBdrElements = 2 * nx + 2 * ny;

      vertices.SetSize(NumOfVertices);
      elements.SetSize(NumOfElements);
      boundary.SetSize(NumOfBdrElements);

      real_t cx, cy;
      int64_t ind[4];

      // Sets vertices and the corresponding coordinates
      k = 0;
      for (j = 0; j < ny+1; j++)
      {
         cy = ((real_t) j / ny) * sy;
         for (i = 0; i < nx+1; i++)
         {
            cx = ((real_t) i / nx) * sx;
            vertices[k](0) = cx;
            vertices[k](1) = cy;
            k++;
         }
      }

      // Sets elements and the corresponding indices of vertices
      if (sfc_ordering)
      {
         Array<int64_t> sfc;
         GridSfcOrdering2D(nx, ny, sfc);
         MFEM_VERIFY(sfc.Size() == 2*nx*ny, "");

         for (k = 0; k < nx*ny; k++)
         {
            i = sfc[2*k + 0];
            j = sfc[2*k + 1];
            ind[0] = i + j*(nx+1);
            ind[1] = i + 1 +j*(nx+1);
            ind[2] = i + 1 + (j+1)*(nx+1);
            ind[3] = i + (j+1)*(nx+1);
            elements[k] = new Quadrilateral(ind);
         }
      }
      else
      {
         k = 0;
         for (j = 0; j < ny; j++)
         {
            for (i = 0; i < nx; i++)
            {
               ind[0] = i + j*(nx+1);
               ind[1] = i + 1 +j*(nx+1);
               ind[2] = i + 1 + (j+1)*(nx+1);
               ind[3] = i + (j+1)*(nx+1);
               elements[k] = new Quadrilateral(ind);
               k++;
            }
         }
      }

      // Sets boundary elements and the corresponding indices of vertices
      int64_t m = (nx+1)*ny;
      for (i = 0; i < nx; i++)
      {
         boundary[i] = new Segment(i, i+1, 1);
         boundary[nx+i] = new Segment(m+i+1, m+i, 3);
      }
      m = nx+1;
      for (j = 0; j < ny; j++)
      {
         boundary[2*nx+j] = new Segment((j+1)*m, j*m, 4);
         boundary[2*nx+ny+j] = new Segment(j*m+nx, (j+1)*m+nx, 2);
      }
   }
   // Creates triangular mesh
   else if (type == Element::TRIANGLE)
   {
      NumOfVertices = (nx+1) * (ny+1);
      NumOfElements = 2 * nx * ny;
      NumOfBdrElements = 2 * nx + 2 * ny;

      vertices.SetSize(NumOfVertices);
      elements.SetSize(NumOfElements);
      boundary.SetSize(NumOfBdrElements);

      real_t cx, cy;
      int64_t ind[3];

      // Sets vertices and the corresponding coordinates
      k = 0;
      for (j = 0; j < ny+1; j++)
      {
         cy = ((real_t) j / ny) * sy;
         for (i = 0; i < nx+1; i++)
         {
            cx = ((real_t) i / nx) * sx;
            vertices[k](0) = cx;
            vertices[k](1) = cy;
            k++;
         }
      }

      // Sets the elements and the corresponding indices of vertices
      k = 0;
      for (j = 0; j < ny; j++)
      {
         for (i = 0; i < nx; i++)
         {
            ind[0] = i + j*(nx+1);
            ind[1] = i + 1 + (j+1)*(nx+1);
            ind[2] = i + (j+1)*(nx+1);
            elements[k] = new Triangle(ind);
            k++;
            ind[1] = i + 1 + j*(nx+1);
            ind[2] = i + 1 + (j+1)*(nx+1);
            elements[k] = new Triangle(ind);
            k++;
         }
      }

      // Sets boundary elements and the corresponding indices of vertices
      int64_t m = (nx+1)*ny;
      for (i = 0; i < nx; i++)
      {
         boundary[i] = new Segment(i, i+1, 1);
         boundary[nx+i] = new Segment(m+i+1, m+i, 3);
      }
      m = nx+1;
      for (j = 0; j < ny; j++)
      {
         boundary[2*nx+j] = new Segment((j+1)*m, j*m, 4);
         boundary[2*nx+ny+j] = new Segment(j*m+nx, (j+1)*m+nx, 2);
      }

      // MarkTriMeshForRefinement(); // done in Finalize(...)
   }
   else
   {
      MFEM_ABORT("Unsupported element type.");
   }

   SetMeshGen();

   if (generate_edges == 1)
   {
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
      GenerateFaces();
   }
   else
   {
      NumOfEdges = 0;
   }

   NumOfFaces = 0;

   attributes.Append(1);
   bdr_attributes.Append(1); bdr_attributes.Append(2);
   bdr_attributes.Append(3); bdr_attributes.Append(4);

   // Finalize(...) can be called after this method, if needed
}

void Mesh::Make1D(int64_t n, real_t sx)
{
   int64_t j, ind[1];

   SetEmpty();

   Dim = 1;
   spaceDim = 1;

   NumOfVertices = n + 1;
   NumOfElements = n;
   NumOfBdrElements = 2;
   vertices.SetSize(NumOfVertices);
   elements.SetSize(NumOfElements);
   boundary.SetSize(NumOfBdrElements);

   // Sets vertices and the corresponding coordinates
   for (j = 0; j < n+1; j++)
   {
      vertices[j](0) = ((real_t) j / n) * sx;
   }

   // Sets elements and the corresponding indices of vertices
   for (j = 0; j < n; j++)
   {
      elements[j] = new Segment(j, j+1, 1);
   }

   // Sets the boundary elements
   ind[0] = 0;
   boundary[0] = new Point(ind, 1);
   ind[0] = n;
   boundary[1] = new Point(ind, 2);

   NumOfEdges = 0;
   NumOfFaces = 0;

   SetMeshGen();
   GenerateFaces();

   // Set be_to_face
   be_to_face.SetSize(2);
   be_to_face[0] = 0;
   be_to_face[1] = n;

   attributes.Append(1);
   bdr_attributes.Append(1); bdr_attributes.Append(2);
}

Mesh::Mesh(const Mesh &mesh, bool copy_nodes)
  : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
{
   Dim = mesh.Dim;
   spaceDim = mesh.spaceDim;

   NumOfVertices = mesh.NumOfVertices;
   NumOfElements = mesh.NumOfElements;
   NumOfBdrElements = mesh.NumOfBdrElements;
   NumOfEdges = mesh.NumOfEdges;
   NumOfFaces = mesh.NumOfFaces;
   nbInteriorFaces = mesh.nbInteriorFaces;
   nbBoundaryFaces = mesh.nbBoundaryFaces;

   meshgen = mesh.meshgen;
   mesh_geoms = mesh.mesh_geoms;

   // Duplicate the elements
   elements.SetSize(NumOfElements);
   for (int64_t i = 0; i < NumOfElements; i++)
   {
      elements[i] = mesh.elements[i]->Duplicate(this);
   }

   // Copy the vertices
   mesh.vertices.Copy(vertices);

   // Duplicate the boundary
   boundary.SetSize(NumOfBdrElements);
   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      boundary[i] = mesh.boundary[i]->Duplicate(this);
   }

   // Copy the element-to-face Table, el_to_face
   el_to_face = (mesh.el_to_face) ? new Table(*mesh.el_to_face) : NULL;

   // Copy the boundary-to-face Array, be_to_face.
   mesh.be_to_face.Copy(be_to_face);

   // Copy the element-to-edge Table, el_to_edge
   el_to_edge = (mesh.el_to_edge) ? new Table(*mesh.el_to_edge) : NULL;

   // Copy the boundary-to-edge Table, bel_to_edge (3D)
   bel_to_edge = (mesh.bel_to_edge) ? new Table(*mesh.bel_to_edge) : NULL;

   // Duplicate the faces and faces_info.
   faces.SetSize(mesh.faces.Size());
   for (int64_t i = 0; i < faces.Size(); i++)
   {
      Element *face = mesh.faces[i]; // in 1D the faces are NULL
      faces[i] = (face) ? face->Duplicate(this) : NULL;
   }
   mesh.faces_info.Copy(faces_info);
   mesh.nc_faces_info.Copy(nc_faces_info);

   // Do NOT copy the element-to-element Table, el_to_el
   el_to_el = NULL;

   // Do NOT copy the face-to-edge Table, face_edge
   face_edge = NULL;
   face_to_elem = NULL;

   // Copy the edge-to-vertex Table, edge_vertex
   edge_vertex = (mesh.edge_vertex) ? new Table(*mesh.edge_vertex) : NULL;

   // Copy the attributes and bdr_attributes
   mesh.attributes.Copy(attributes);
   mesh.bdr_attributes.Copy(bdr_attributes);

   // Copy attribute and bdr_attribute names
   mesh.attribute_sets.Copy(attribute_sets);
   mesh.bdr_attribute_sets.Copy(bdr_attribute_sets);
}

Mesh::Mesh(Mesh &&mesh) : Mesh()
{
   Swap(mesh, true);
}

Mesh& Mesh::operator=(Mesh &&mesh)
{
   Swap(mesh, true);
   return *this;
}

Mesh Mesh::LoadFromFile(const std::string &filename, int64_t generate_edges,
                        int64_t refine, bool fix_orientation)
{
   Mesh mesh;
   named_ifgzstream imesh(filename);
   if (!imesh) { MFEM_ABORT("Mesh file not found: " << filename << '\n'); }
   else { mesh.Load(imesh, generate_edges, refine, fix_orientation); }
   return mesh;
}

Mesh Mesh::MakeCartesian1D(int64_t n, real_t sx)
{
   Mesh mesh;
   mesh.Make1D(n, sx);
   // mesh.Finalize(); not needed in this case
   return mesh;
}

Mesh Mesh::MakeCartesian2D(
   int64_t nx, int64_t ny, Element::Type type, bool generate_edges,
   real_t sx, real_t sy, bool sfc_ordering)
{
   Mesh mesh;
   mesh.Make2D(nx, ny, type, sx, sy, generate_edges, sfc_ordering);
   mesh.Finalize(true); // refine = true
   return mesh;
}

Mesh Mesh::MakeCartesian3D(
   int64_t nx, int64_t ny, int64_t nz, Element::Type type,
   real_t sx, real_t sy, real_t sz, bool sfc_ordering)
{
   Mesh mesh;
   mesh.Make3D(nx, ny, nz, type, sx, sy, sz, sfc_ordering);
   mesh.Finalize(true); // refine = true
   return mesh;
}

Mesh Mesh::MakeCartesian3DWith24TetsPerHex(int64_t nx, int64_t ny, int64_t nz,
                              real_t sx, real_t sy, real_t sz)
{
   Mesh mesh;
   mesh.Make3D24TetsFromHex(nx, ny, nz, sx, sy, sz);
   mesh.Finalize(false, false);
   return mesh;
}

Mesh Mesh::MakeCartesian2DWith4TrisPerQuad(int64_t nx, int64_t ny,
                                           real_t sx, real_t sy)
{
   Mesh mesh;
   mesh.Make2D4TrisFromQuad(nx, ny, sx, sy);
   mesh.Finalize(false, false);
   return mesh;
}

Mesh Mesh::MakeCartesian2DWith5QuadsPerQuad(int64_t nx, int64_t ny,
                                            real_t sx, real_t sy)
{
   Mesh mesh;
   mesh.Make2D5QuadsFromQuad(nx, ny, sx, sy);
   mesh.Finalize(false, false);
   return mesh;
}

Mesh::Mesh(const std::string &filename, int64_t generate_edges, int64_t refine,
           bool fix_orientation)
 : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
{
   // Initialization as in the default constructor
   SetEmpty();

   named_ifgzstream imesh(filename);
   if (!imesh)
   {
   }
   else
   {
      Load(imesh, generate_edges, refine, fix_orientation);
   }
}

Mesh::Mesh(std::istream &input, int64_t generate_edges, int64_t refine,
           bool fix_orientation)
 : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
{
   SetEmpty();
   Load(input, generate_edges, refine, fix_orientation);
}

Mesh::Mesh(real_t *vertices_, int64_t num_vertices,
           int64_t *element_indices, Geometry::Type element_type,
           int64_t *element_attributes, int64_t num_elements,
           int64_t *boundary_indices, Geometry::Type boundary_type,
           int64_t *boundary_attributes, int64_t num_boundary_elements,
           int64_t dimension, int64_t space_dimension)
 : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
{
   if (space_dimension == -1)
   {
      space_dimension = dimension;
   }

   InitMesh(dimension, space_dimension, /*num_vertices*/ 0, num_elements,
            num_boundary_elements);

   int64_t element_index_stride = Geometry::NumVerts[element_type];
   int64_t boundary_index_stride = num_boundary_elements > 0 ?
                               Geometry::NumVerts[boundary_type] : 0;

   // assuming Vertex is POD
   vertices.MakeRef(reinterpret_cast<Vertex*>(vertices_), num_vertices);
   NumOfVertices = num_vertices;

   for (int64_t i = 0; i < num_elements; i++)
   {
      elements[i] = NewElement(element_type);
      elements[i]->SetVertices(element_indices + i * element_index_stride);
      elements[i]->SetAttribute(element_attributes[i]);
   }
   NumOfElements = num_elements;

   for (int64_t i = 0; i < num_boundary_elements; i++)
   {
      boundary[i] = NewElement(boundary_type);
      boundary[i]->SetVertices(boundary_indices + i * boundary_index_stride);
      boundary[i]->SetAttribute(boundary_attributes[i]);
   }
   NumOfBdrElements = num_boundary_elements;

   FinalizeTopology();
}

Element *Mesh::NewElement(int64_t geom)
{
   switch (geom)
   {
      case Geometry::POINT:     return (new Point);
      case Geometry::SEGMENT:   return (new Segment);
      case Geometry::TRIANGLE:  return (new Triangle);
      case Geometry::SQUARE:    return (new Quadrilateral);
      case Geometry::TETRAHEDRON: return (new Tetrahedron);
      case Geometry::CUBE:      return (new Hexahedron);
      case Geometry::PRISM:     return (new Wedge);
      case Geometry::PYRAMID:   return (new Pyramid);
      default:
         MFEM_ABORT("invalid Geometry::Type, geom = " << geom);
   }

   return NULL;
}

Element *Mesh::ReadElementWithoutAttr(std::istream &input)
{
   int64_t geom, nv, *v;
   Element *el;

   input >> geom;
   el = NewElement(geom);
   MFEM_VERIFY(el, "Unsupported element type: " << geom);
   nv = el->GetNVertices();
   v  = el->GetVertices();
   for (int64_t i = 0; i < nv; i++)
   {
      input >> v[i];
   }

   return el;
}

void Mesh::PrintElementWithoutAttr(const Element *el, std::ostream &os)
{
   os << el->GetGeometryType();
   const int64_t nv = el->GetNVertices();
   const int64_t *v = el->GetVertices();
   for (int64_t j = 0; j < nv; j++)
   {
      os << ' ' << v[j];
   }
   os << '\n';
}

Element *Mesh::ReadElement(std::istream &input)
{
   int64_t attr;
   Element *el;

   input >> attr;
   el = ReadElementWithoutAttr(input);
   el->SetAttribute(attr);

   return el;
}

void Mesh::PrintElement(const Element *el, std::ostream &os)
{
   os << el->GetAttribute() << ' ';
   PrintElementWithoutAttr(el, os);
}

void Mesh::SetMeshGen()
{
   meshgen = mesh_geoms = 0;
   for (int64_t i = 0; i < NumOfElements; i++)
   {
      const Element::Type type = GetElement(i)->GetType();
      switch (type)
      {
         case Element::TETRAHEDRON:
            mesh_geoms |= (1 << Geometry::TETRAHEDRON);
         case Element::TRIANGLE:
            mesh_geoms |= (1 << Geometry::TRIANGLE);
         case Element::SEGMENT:
            mesh_geoms |= (1 << Geometry::SEGMENT);
         case Element::POINT:
            mesh_geoms |= (1 << Geometry::POINT);
            meshgen |= 1;
            break;

         case Element::HEXAHEDRON:
            mesh_geoms |= (1 << Geometry::CUBE);
         case Element::QUADRILATERAL:
            mesh_geoms |= (1 << Geometry::SQUARE);
            mesh_geoms |= (1 << Geometry::SEGMENT);
            mesh_geoms |= (1 << Geometry::POINT);
            meshgen |= 2;
            break;

         case Element::WEDGE:
            mesh_geoms |= (1 << Geometry::PRISM);
            mesh_geoms |= (1 << Geometry::SQUARE);
            mesh_geoms |= (1 << Geometry::TRIANGLE);
            mesh_geoms |= (1 << Geometry::SEGMENT);
            mesh_geoms |= (1 << Geometry::POINT);
            meshgen |= 4;
            break;

         case Element::PYRAMID:
            mesh_geoms |= (1 << Geometry::PYRAMID);
            mesh_geoms |= (1 << Geometry::SQUARE);
            mesh_geoms |= (1 << Geometry::TRIANGLE);
            mesh_geoms |= (1 << Geometry::SEGMENT);
            mesh_geoms |= (1 << Geometry::POINT);
            meshgen |= 8;
            break;

         default:
            MFEM_ABORT("invalid element type: " << type);
            break;
      }
   }
}

void Mesh::Loader(std::istream &input, int64_t generate_edges,
                  std::string parse_tag)
{
   int64_t curved = 0, read_gf = 1;
   bool finalize_topo = true;

   if (!input)
   {
      MFEM_ABORT("Input stream is not open");
   }

   Clear();

   string mesh_type;
   input >> ws;
   getline(input, mesh_type);
   filter_dos(mesh_type);

   // MFEM's conforming mesh formats
   int64_t mfem_version = 0;
   if (mesh_type == "MFEM mesh v1.0") { mfem_version = 10; } // serial
   else if (mesh_type == "MFEM mesh v1.2") { mfem_version = 12; } // parallel
   else if (mesh_type == "MFEM mesh v1.3") { mfem_version = 13; } // attr sets

   // MFEM nonconforming mesh format
   // (NOTE: previous v1.1 is now under this branch for backward compatibility)
   int64_t mfem_nc_version = 0;
   if (mesh_type == "MFEM NC mesh v1.0") { mfem_nc_version = 10; }
   else if (mesh_type == "MFEM NC mesh v1.1") { mfem_nc_version = 11; }
   else if (mesh_type == "MFEM mesh v1.1") { mfem_nc_version = 1 /*legacy*/; }

   if (mfem_version)
   {
      // Formats mfem_v12 and newer have a tag indicating the end of the mesh
      // section in the stream. A user provided parse tag can also be provided
      // via the arguments. For example, if this is called from parallel mesh
      // object, it can indicate to read until parallel mesh section begins.
      if (mfem_version >= 12 && parse_tag.empty())
      {
         parse_tag = "mfem_mesh_end";
      }
      ReadMFEMMesh(input, mfem_version, curved);
   }
   else if (mfem_nc_version)
   {
      MFEM_ABORT("");
   }
   else if (mesh_type == "linemesh") // 1D mesh
   {
      ReadLineMesh(input);
   }
   else if (mesh_type == "areamesh2" || mesh_type == "curved_areamesh2")
   {
      if (mesh_type == "curved_areamesh2")
      {
         curved = 1;
      }
      ReadNetgen2DMesh(input, curved);
   }
   else if (mesh_type == "NETGEN" || mesh_type == "NETGEN_Neutral_Format")
   {
      ReadNetgen3DMesh(input);
   }
   else if (mesh_type == "TrueGrid")
   {
      ReadTrueGridMesh(input);
   }
   else if (mesh_type.rfind("# vtk DataFile Version") == 0)
   {
      int64_t major_vtk_version = mesh_type[mesh_type.length()-3] - '0';
      // int64_t minor_vtk_version = mesh_type[mesh_type.length()-1] - '0';
      MFEM_VERIFY(major_vtk_version >= 2 && major_vtk_version <= 4,
                  "Unsupported VTK format");
      ReadVTKMesh(input, curved, read_gf, finalize_topo);
   }
   else if (mesh_type.rfind("<VTKFile ") == 0 || mesh_type.rfind("<?xml") == 0)
   {
      ReadXML_VTKMesh(input, curved, read_gf, finalize_topo, mesh_type);
   }
   else if (mesh_type == "MFEM NURBS mesh v1.0")
   {
      // ReadNURBSMesh(input, curved, read_gf);
      MFEM_ABORT("");
   }
   else if (mesh_type == "MFEM NURBS mesh v1.1")
   {
   //   ReadNURBSMesh(input, curved, read_gf, true);
      MFEM_ABORT("");
   }
   else if (mesh_type == "MFEM INLINE mesh v1.0")
   {
      ReadInlineMesh(input, generate_edges);
      return; // done with inline mesh construction
   }
   else if (mesh_type == "$MeshFormat") // Gmsh
   {
      ReadGmshMesh(input, curved, read_gf);
   }
   else if
   ((mesh_type.size() > 2 &&
     mesh_type[0] == 'C' && mesh_type[1] == 'D' && mesh_type[2] == 'F') ||
    (mesh_type.size() > 3 &&
     mesh_type[1] == 'H' && mesh_type[2] == 'D' && mesh_type[3] == 'F'))
   {
      named_ifgzstream *mesh_input = dynamic_cast<named_ifgzstream *>(&input);
      if (mesh_input)
      {
#ifdef MFEM_USE_NETCDF
         ReadCubit(mesh_input->filename, curved, read_gf);
#else
         MFEM_ABORT("NetCDF support requires configuration with"
                    " MFEM_USE_NETCDF=YES");
         return;
#endif
      }
      else
      {
         MFEM_ABORT("Can not determine Cubit mesh filename!"
                    " Use mfem::named_ifgzstream for input.");
         return;
      }
   }
   else
   {
      MFEM_ABORT("Unknown input mesh format: " << mesh_type);
      return;
   }

   // at this point the following should be defined:
   //  1) Dim
   //  2) NumOfElements, elements
   //  3) NumOfBdrElements, boundary
   //  4) NumOfVertices, with allocated space in vertices
   //  5) curved
   //  5a) if curved == 0, vertices must be defined
   //  5b) if curved != 0 and read_gf != 0,
   //         'input' must point to a GridFunction
   //  5c) if curved != 0 and read_gf == 0,
   //         vertices and Nodes must be defined
   // optional:
   //  1) el_to_edge may be allocated (as in the case of P2 VTK meshes)
   //  2) ncmesh may be allocated

   // FinalizeTopology() will:
   // - assume that generate_edges is true
   // - assume that refine is false
   // - does not check the orientation of regular and boundary elements
   if (finalize_topo)
   {
      // don't generate any boundary elements, especially in parallel
      bool generate_bdr = false;

      FinalizeTopology(generate_bdr);
   }

   // If a parse tag was supplied, keep reading the stream until the tag is
   // encountered.
   if (mfem_version >= 12)
   {
      string line;
      do
      {
         skip_comment_lines(input, '#');
         MFEM_VERIFY(input.good(), "Required mesh-end tag not found");
         getline(input, line);
         filter_dos(line);
         // mfem v1.2 may not have parse_tag in it, e.g. if trying to read a
         // serial mfem v1.2 mesh as parallel with "mfem_serial_mesh_end" as
         // parse_tag. That's why, regardless of parse_tag, we stop reading if
         // we find "mfem_mesh_end" which is required by mfem v1.2 format.
         if (line == "mfem_mesh_end") { break; }
      }
      while (line != parse_tag);
   }
   else if (mfem_nc_version >= 10)
   {
      string ident;
      skip_comment_lines(input, '#');
      input >> ident;
      MFEM_VERIFY(ident == "mfem_mesh_end",
                  "invalid mesh: end of file tag not found");
   }

   // Finalize(...) should be called after this, if needed.
}

Mesh::Mesh(Mesh *mesh_array[], int64_t num_pieces)
 : attribute_sets(attributes), bdr_attribute_sets(bdr_attributes)
{
   int64_t      i, j, ie, ib, iv, *v, nv;
   Element *el;
   Mesh    *m;

   SetEmpty();

   Dim = mesh_array[0]->Dimension();
   spaceDim = mesh_array[0]->SpaceDimension();

   NumOfElements    = 0;
   NumOfBdrElements = 0;
   NumOfVertices    = 0;
   for (i = 0; i < num_pieces; i++)
   {
      m = mesh_array[i];
      NumOfElements    += m->GetNE();
      NumOfBdrElements += m->GetNBE();
      NumOfVertices    += m->GetNV();
   }
   elements.SetSize(NumOfElements);
   boundary.SetSize(NumOfBdrElements);
   vertices.SetSize(NumOfVertices);
   ie = ib = iv = 0;
   for (i = 0; i < num_pieces; i++)
   {
      m = mesh_array[i];
      // copy the elements
      for (j = 0; j < m->GetNE(); j++)
      {
         el = m->GetElement(j)->Duplicate(this);
         v  = el->GetVertices();
         nv = el->GetNVertices();
         for (int64_t k = 0; k < nv; k++)
         {
            v[k] += iv;
         }
         elements[ie++] = el;
      }
      // copy the boundary elements
      for (j = 0; j < m->GetNBE(); j++)
      {
         el = m->GetBdrElement(j)->Duplicate(this);
         v  = el->GetVertices();
         nv = el->GetNVertices();
         for (int64_t k = 0; k < nv; k++)
         {
            v[k] += iv;
         }
         boundary[ib++] = el;
      }
      // copy the vertices
      for (j = 0; j < m->GetNV(); j++)
      {
         vertices[iv++].SetCoords(m->SpaceDimension(), m->GetVertex(j));
      }
   }

   FinalizeTopology();
}

int64_t Mesh::GetNumFaces() const
{
   switch (Dim)
   {
      case 1: return GetNV();
      case 2: return GetNEdges();
      case 3: return GetNFaces();
   }
   return 0;
}

int64_t Mesh::GetNumFacesWithGhost() const
{
   return faces_info.Size();
}

int64_t Mesh::GetNFbyType(FaceType type) const
{
   const bool isInt = type==FaceType::Interior;
   int64_t &nf = isInt ? nbInteriorFaces : nbBoundaryFaces;
   if (nf<0)
   {
      nf = 0;
      for (int64_t f = 0; f < GetNumFacesWithGhost(); ++f)
      {
         FaceInformation face = GetFaceInformation(f);
         if (face.IsOfFaceType(type))
         {
            if (face.IsNonconformingCoarse())
            {
               // We don't count nonconforming coarse faces.
               continue;
            }
            nf++;
         }
      }
   }
   return nf;
}

#if (!defined(MFEM_USE_MPI) || defined(MFEM_DEBUG))
static const char *fixed_or_not[] = { "fixed", "NOT FIXED" };
#endif

int64_t Mesh::GetTriOrientation(const int64_t *base, const int64_t *test)
{
   // Static method.
   // This function computes the index 'j' of the permutation that transforms
   // test into base: test[tri_orientation[j][i]]=base[i].
   // tri_orientation = Geometry::Constants<Geometry::TRIANGLE>::Orient
   int64_t orient;

   if (test[0] == base[0])
      if (test[1] == base[1])
      {
         orient = 0;   //  (0, 1, 2)
      }
      else
      {
         orient = 5;   //  (0, 2, 1)
      }
   else if (test[0] == base[1])
      if (test[1] == base[0])
      {
         orient = 1;   //  (1, 0, 2)
      }
      else
      {
         orient = 2;   //  (1, 2, 0)
      }
   else // test[0] == base[2]
      if (test[1] == base[0])
      {
         orient = 4;   //  (2, 0, 1)
      }
      else
      {
         orient = 3;   //  (2, 1, 0)
      }

#ifdef MFEM_DEBUG
   const int64_t *aor = tri_t::Orient[orient];
   for (int64_t j = 0; j < 3; j++)
      if (test[aor[j]] != base[j])
      {
         mfem::err << "Mesh::GetTriOrientation(...)" << endl;
         mfem::err << " base = [";
         for (int64_t k = 0; k < 3; k++)
         {
            mfem::err << " " << base[k];
         }
         mfem::err << " ]\n test = [";
         for (int64_t k = 0; k < 3; k++)
         {
            mfem::err << " " << test[k];
         }
         mfem::err << " ]" << endl;
         mfem_error();
      }
#endif

   return orient;
}

int64_t Mesh::ComposeTriOrientations(int64_t ori_a_b, int64_t ori_b_c)
{
   // Static method.
   // Given three, possibly different, configurations of triangular face
   // vertices: va, vb, and vc.  This function returns the relative orientation
   // GetTriOrientation(va, vc) by composing previously computed orientations
   // ori_a_b = GetTriOrientation(va, vb) and
   // ori_b_c = GetTriOrientation(vb, vc) without accessing the vertices.

   const int64_t oo[6][6] =
   {
      {0, 1, 2, 3, 4, 5},
      {1, 0, 5, 4, 3, 2},
      {2, 3, 4, 5, 0, 1},
      {3, 2, 1, 0, 5, 4},
      {4, 5, 0, 1, 2, 3},
      {5, 4, 3, 2, 1, 0}
   };

   int64_t ori_a_c = oo[ori_a_b][ori_b_c];
   return ori_a_c;
}

int64_t Mesh::InvertTriOrientation(int64_t ori)
{
   const int64_t inv_ori[6] = {0, 1, 4, 3, 2, 5};
   return inv_ori[ori];
}

int64_t Mesh::GetQuadOrientation(const int64_t *base, const int64_t *test)
{
   int64_t i;

   for (i = 0; i < 4; i++)
      if (test[i] == base[0])
      {
         break;
      }

#ifdef MFEM_DEBUG
   int64_t orient;
   if (test[(i+1)%4] == base[1])
   {
      orient = 2*i;
   }
   else
   {
      orient = 2*i+1;
   }
   const int64_t *aor = quad_t::Orient[orient];
   for (int64_t j = 0; j < 4; j++)
      if (test[aor[j]] != base[j])
      {
         mfem::err << "Mesh::GetQuadOrientation(...)" << endl;
         mfem::err << " base = [";
         for (int64_t k = 0; k < 4; k++)
         {
            mfem::err << " " << base[k];
         }
         mfem::err << " ]\n test = [";
         for (int64_t k = 0; k < 4; k++)
         {
            mfem::err << " " << test[k];
         }
         mfem::err << " ]" << endl;
         mfem_error();
      }
#endif

   if (test[(i+1)%4] == base[1])
   {
      return 2*i;
   }

   return 2*i+1;
}

int64_t Mesh::ComposeQuadOrientations(int64_t ori_a_b, int64_t ori_b_c)
{
   // Static method.
   // Given three, possibly different, configurations of quadrilateral face
   // vertices: va, vb, and vc.  This function returns the relative orientation
   // GetQuadOrientation(va, vc) by composing previously computed orientations
   // ori_a_b = GetQuadOrientation(va, vb) and
   // ori_b_c = GetQuadOrientation(vb, vc) without accessing the vertices.

   const int64_t oo[8][8] =
   {
      {0, 1, 2, 3, 4, 5, 6, 7},
      {1, 0, 3, 2, 5, 4, 7, 6},
      {2, 7, 4, 1, 6, 3, 0, 5},
      {3, 6, 5, 0, 7, 2, 1, 4},
      {4, 5, 6, 7, 0, 1, 2, 3},
      {5, 4, 7, 6, 1, 0, 3, 2},
      {6, 3, 0, 5, 2, 7, 4, 1},
      {7, 2, 1, 4, 3, 6, 5, 0}
   };

   int64_t ori_a_c = oo[ori_a_b][ori_b_c];
   return ori_a_c;
}

int64_t Mesh::InvertQuadOrientation(int64_t ori)
{
   const int64_t inv_ori[8] = {0, 1, 6, 3, 4, 5, 2, 7};
   return inv_ori[ori];
}

int64_t Mesh::GetTetOrientation(const int64_t *base, const int64_t *test)
{
   // Static method.
   // This function computes the index 'j' of the permutation that transforms
   // test into base: test[tet_orientation[j][i]]=base[i].
   // tet_orientation = Geometry::Constants<Geometry::TETRAHEDRON>::Orient
   int64_t orient;

   if (test[0] == base[0])
      if (test[1] == base[1])
         if (test[2] == base[2])
         {
            orient = 0;   //  (0, 1, 2, 3)
         }
         else
         {
            orient = 1;   //  (0, 1, 3, 2)
         }
      else if (test[2] == base[1])
         if (test[3] == base[2])
         {
            orient = 2;   //  (0, 2, 3, 1)
         }
         else
         {
            orient = 3;   //  (0, 2, 1, 3)
         }
      else // test[3] == base[1]
         if (test[1] == base[2])
         {
            orient = 4;   //  (0, 3, 1, 2)
         }
         else
         {
            orient = 5;   //  (0, 3, 2, 1)
         }
   else if (test[1] == base[0])
      if (test[2] == base[1])
         if (test[0] == base[2])
         {
            orient = 6;   //  (1, 2, 0, 3)
         }
         else
         {
            orient = 7;   //  (1, 2, 3, 0)
         }
      else if (test[3] == base[1])
         if (test[2] == base[2])
         {
            orient = 8;   //  (1, 3, 2, 0)
         }
         else
         {
            orient = 9;   //  (1, 3, 0, 2)
         }
      else // test[0] == base[1]
         if (test[3] == base[2])
         {
            orient = 10;   //  (1, 0, 3, 2)
         }
         else
         {
            orient = 11;   //  (1, 0, 2, 3)
         }
   else if (test[2] == base[0])
      if (test[3] == base[1])
         if (test[0] == base[2])
         {
            orient = 12;   //  (2, 3, 0, 1)
         }
         else
         {
            orient = 13;   //  (2, 3, 1, 0)
         }
      else if (test[0] == base[1])
         if (test[1] == base[2])
         {
            orient = 14;   //  (2, 0, 1, 3)
         }
         else
         {
            orient = 15;   //  (2, 0, 3, 1)
         }
      else // test[1] == base[1]
         if (test[3] == base[2])
         {
            orient = 16;   //  (2, 1, 3, 0)
         }
         else
         {
            orient = 17;   //  (2, 1, 0, 3)
         }
   else // (test[3] == base[0])
      if (test[0] == base[1])
         if (test[2] == base[2])
         {
            orient = 18;   //  (3, 0, 2, 1)
         }
         else
         {
            orient = 19;   //  (3, 0, 1, 2)
         }
      else if (test[1] == base[1])
         if (test[0] == base[2])
         {
            orient = 20;   //  (3, 1, 0, 2)
         }
         else
         {
            orient = 21;   //  (3, 1, 2, 0)
         }
      else // test[2] == base[1]
         if (test[1] == base[2])
         {
            orient = 22;   //  (3, 2, 1, 0)
         }
         else
         {
            orient = 23;   //  (3, 2, 0, 1)
         }

#ifdef MFEM_DEBUG
   const int64_t *aor = tet_t::Orient[orient];
   for (int64_t j = 0; j < 4; j++)
      if (test[aor[j]] != base[j])
      {
         mfem_error("Mesh::GetTetOrientation(...)");
      }
#endif

   return orient;
}

int64_t Mesh::GetNumGeometries(int64_t dim) const
{
   MFEM_ASSERT(0 <= dim && dim <= Dim, "invalid dim: " << dim);
   int64_t num_geoms = 0;
   for (int64_t g = Geometry::DimStart[dim]; g < Geometry::DimStart[dim+1]; g++)
   {
      if (HasGeometry(Geometry::Type(g))) { num_geoms++; }
   }
   return num_geoms;
}

void Mesh::GetGeometries(int64_t dim, Array<Geometry::Type> &el_geoms) const
{
   MFEM_ASSERT(0 <= dim && dim <= Dim, "invalid dim: " << dim);
   el_geoms.SetSize(0);
   for (int64_t g = Geometry::DimStart[dim]; g < Geometry::DimStart[dim+1]; g++)
   {
      if (HasGeometry(Geometry::Type(g)))
      {
         el_geoms.Append(Geometry::Type(g));
      }
   }
}

void Mesh::GetElementEdges(int64_t i, Array<int64_t> &edges, Array<int64_t> &cor) const
{
   if (el_to_edge)
   {
      el_to_edge->GetRow(i, edges);
   }
   else
   {
      mfem_error("Mesh::GetElementEdges(...) element to edge table "
                 "is not generated.");
   }

   const int64_t *v = elements[i]->GetVertices();
   const int64_t ne = elements[i]->GetNEdges();
   cor.SetSize(ne);
   for (int64_t j = 0; j < ne; j++)
   {
      const int64_t *e = elements[i]->GetEdgeVertices(j);
      cor[j] = (v[e[0]] < v[e[1]]) ? (1) : (-1);
   }
}

void Mesh::GetBdrElementEdges(int64_t i, Array<int64_t> &edges, Array<int64_t> &cor) const
{
   if (Dim == 2)
   {
      edges.SetSize(1);
      cor.SetSize(1);
      edges[0] = be_to_face[i];
      const int64_t *v = boundary[i]->GetVertices();
      cor[0] = (v[0] < v[1]) ? (1) : (-1);
   }
   else if (Dim == 3)
   {
      if (bel_to_edge)
      {
         bel_to_edge->GetRow(i, edges);
      }
      else
      {
         mfem_error("Mesh::GetBdrElementEdges(...)");
      }

      const int64_t *v = boundary[i]->GetVertices();
      const int64_t ne = boundary[i]->GetNEdges();
      cor.SetSize(ne);
      for (int64_t j = 0; j < ne; j++)
      {
         const int64_t *e = boundary[i]->GetEdgeVertices(j);
         cor[j] = (v[e[0]] < v[e[1]]) ? (1) : (-1);
      }
   }
}

void Mesh::GetFaceEdges(int64_t i, Array<int64_t> &edges, Array<int64_t> &o) const
{
   if (Dim == 2)
   {
      edges.SetSize(1);
      edges[0] = i;
      o.SetSize(1);
      const int64_t *v = faces[i]->GetVertices();
      o[0] = (v[0] < v[1]) ? (1) : (-1);
   }

   if (Dim != 3)
   {
      return;
   }

   GetFaceEdgeTable(); // generate face_edge Table (if not generated)

   face_edge->GetRow(i, edges);

   const int64_t *v = faces[i]->GetVertices();
   const int64_t ne = faces[i]->GetNEdges();
   o.SetSize(ne);
   for (int64_t j = 0; j < ne; j++)
   {
      const int64_t *e = faces[i]->GetEdgeVertices(j);
      o[j] = (v[e[0]] < v[e[1]]) ? (1) : (-1);
   }
}

void Mesh::GetEdgeVertices(int64_t i, Array<int64_t> &vert) const
{
   // the two vertices are sorted: vert[0] < vert[1]
   // this is consistent with the global edge orientation
   // generate edge_vertex Table (if not generated)
   if (!edge_vertex) { GetEdgeVertexTable(); }
   edge_vertex->GetRow(i, vert);
}

Table *Mesh::GetFaceEdgeTable() const
{
   if (face_edge)
   {
      return face_edge;
   }

   if (Dim != 3)
   {
      return NULL;
   }

#ifdef MFEM_DEBUG
   if (faces.Size() != NumOfFaces)
   {
      mfem_error("Mesh::GetFaceEdgeTable : faces were not generated!");
   }
#endif

   DSTable v_to_v(NumOfVertices);
   GetVertexToVertexTable(v_to_v);

   face_edge = new Table;
   GetElementArrayEdgeTable(faces, v_to_v, *face_edge);

   return (face_edge);
}

Table *Mesh::GetEdgeVertexTable() const
{
   if (edge_vertex)
   {
      return edge_vertex;
   }

   DSTable v_to_v(NumOfVertices);
   GetVertexToVertexTable(v_to_v);

   int64_t nedges = v_to_v.NumberOfEntries();
   edge_vertex = new Table(nedges, 2);
   for (int64_t i = 0; i < NumOfVertices; i++)
   {
      for (DSTable::RowIterator it(v_to_v, i); !it; ++it)
      {
         int64_t j = it.Index();
         edge_vertex->Push(j, i);
         edge_vertex->Push(j, it.Column());
      }
   }
   edge_vertex->Finalize();

   return edge_vertex;
}

Table *Mesh::GetVertexToElementTable()
{
   Table *vert_elem = new Table;

   vert_elem->MakeI(NumOfVertices);

   for (int64_t i = 0; i < NumOfElements; i++)
   {
      const int64_t nv = elements[i]->GetNVertices();
      const int64_t *v = elements[i]->GetVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         vert_elem->AddAColumnInRow(v[j]);
      }
   }

   vert_elem->MakeJ();

   for (int64_t i = 0; i < NumOfElements; i++)
   {
      const int64_t nv = elements[i]->GetNVertices();
      const int64_t *v = elements[i]->GetVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         vert_elem->AddConnection(v[j], i);
      }
   }

   vert_elem->ShiftUpI();

   return vert_elem;
}

Table *Mesh::GetVertexToBdrElementTable()
{
   Table *vert_bdr_elem = new Table;

   vert_bdr_elem->MakeI(NumOfVertices);

   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      const int64_t nv = boundary[i]->GetNVertices();
      const int64_t *v = boundary[i]->GetVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         vert_bdr_elem->AddAColumnInRow(v[j]);
      }
   }

   vert_bdr_elem->MakeJ();

   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      const int64_t nv = boundary[i]->GetNVertices();
      const int64_t *v = boundary[i]->GetVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         vert_bdr_elem->AddConnection(v[j], i);
      }
   }

   vert_bdr_elem->ShiftUpI();

   return vert_bdr_elem;
}

Table *Mesh::GetFaceToElementTable() const
{
   Table *face_elem = new Table;

   face_elem->MakeI(faces_info.Size());

   for (int64_t i = 0; i < faces_info.Size(); i++)
   {
      if (faces_info[i].Elem2No >= 0)
      {
         face_elem->AddColumnsInRow(i, 2);
      }
      else
      {
         face_elem->AddAColumnInRow(i);
      }
   }

   face_elem->MakeJ();

   for (int64_t i = 0; i < faces_info.Size(); i++)
   {
      face_elem->AddConnection(i, faces_info[i].Elem1No);
      if (faces_info[i].Elem2No >= 0)
      {
         face_elem->AddConnection(i, faces_info[i].Elem2No);
      }
   }

   face_elem->ShiftUpI();

   return face_elem;
}

void Mesh::GetElementFaces(int64_t i, Array<int64_t> &el_faces, Array<int64_t> &ori) const
{
   MFEM_VERIFY(el_to_face != NULL, "el_to_face not generated");

   el_to_face->GetRow(i, el_faces);

   int64_t n = el_faces.Size();
   ori.SetSize(n);

   for (int64_t j = 0; j < n; j++)
   {
      if (faces_info[el_faces[j]].Elem1No == i)
      {
         ori[j] = faces_info[el_faces[j]].Elem1Inf % 64;
      }
      else
      {
         MFEM_ASSERT(faces_info[el_faces[j]].Elem2No == i, "internal error");
         ori[j] = faces_info[el_faces[j]].Elem2Inf % 64;
      }
   }
}

Array<int64_t> Mesh::FindFaceNeighbors(const int64_t elem) const
{
   if (face_to_elem == NULL)
   {
      face_to_elem = GetFaceToElementTable();
   }

   Array<int64_t> elem_faces;
   Array<int64_t> ori;
   GetElementFaces(elem, elem_faces, ori);

   Array<int64_t> nghb;
   for (auto f : elem_faces)
   {
      Array<int64_t> row;
      face_to_elem->GetRow(f, row);
      for (auto r : row)
      {
         nghb.Append(r);
      }
   }

   nghb.Sort();
   nghb.Unique();

   return nghb;
}

void Mesh::GetBdrElementFace(int64_t i, int64_t *f, int64_t *o) const
{
   *f = GetBdrElementFaceIndex(i);

   const int64_t *fv = (Dim > 1) ? faces[*f]->GetVertices() : NULL;
   const int64_t *bv = boundary[i]->GetVertices();

   // find the orientation of the bdr. elem. w.r.t.
   // the corresponding face element (that's the base)
   switch (GetBdrElementGeometry(i))
   {
      case Geometry::POINT:    *o = 0; break;
      case Geometry::SEGMENT:  *o = (fv[0] == bv[0]) ? 0 : 1; break;
      case Geometry::TRIANGLE: *o = GetTriOrientation(fv, bv); break;
      case Geometry::SQUARE:   *o = GetQuadOrientation(fv, bv); break;
      default: MFEM_ABORT("invalid geometry");
   }
}

void Mesh::GetBdrElementAdjacentElement(int64_t bdr_el, int64_t &el, int64_t &info) const
{
   int64_t fid = GetBdrElementFaceIndex(bdr_el);

   const FaceInfo &fi = faces_info[fid];
   MFEM_ASSERT(fi.Elem1Inf % 64 == 0, "internal error"); // orientation == 0

   const int64_t *fv = (Dim > 1) ? faces[fid]->GetVertices() : NULL;
   const int64_t *bv = boundary[bdr_el]->GetVertices();
   int64_t ori;
   switch (GetBdrElementGeometry(bdr_el))
   {
      case Geometry::POINT:    ori = 0; break;
      case Geometry::SEGMENT:  ori = (fv[0] == bv[0]) ? 0 : 1; break;
      case Geometry::TRIANGLE: ori = GetTriOrientation(fv, bv); break;
      case Geometry::SQUARE:   ori = GetQuadOrientation(fv, bv); break;
      default: MFEM_ABORT("boundary element type not implemented"); ori = 0;
   }
   el   = fi.Elem1No;
   info = fi.Elem1Inf + ori;
}

void Mesh::SetAttribute(int64_t i, int64_t attr)
{
  elements[i]->SetAttribute(attr);
}

Element::Type Mesh::GetElementType(int64_t i) const
{
   return elements[i]->GetType();
}

Element::Type Mesh::GetBdrElementType(int64_t i) const
{
   return boundary[i]->GetType();
}

void Mesh::GetPointMatrix(int64_t i, DenseMatrix &pointmat) const
{
   int64_t k, j, nv;
   const int64_t *v;

   v  = elements[i]->GetVertices();
   nv = elements[i]->GetNVertices();

   pointmat.SetSize(spaceDim, nv);
   for (k = 0; k < spaceDim; k++)
   {
      for (j = 0; j < nv; j++)
      {
         pointmat(k, j) = vertices[v[j]](k);
      }
   }
}

void Mesh::GetBdrPointMatrix(int64_t i,DenseMatrix &pointmat) const
{
   int64_t k, j, nv;
   const int64_t *v;

   v  = boundary[i]->GetVertices();
   nv = boundary[i]->GetNVertices();

   pointmat.SetSize(spaceDim, nv);
   for (k = 0; k < spaceDim; k++)
      for (j = 0; j < nv; j++)
      {
         pointmat(k, j) = vertices[v[j]](k);
      }
}

real_t Mesh::GetLength(int64_t i, int64_t j) const
{
   const real_t *vi = vertices[i]();
   const real_t *vj = vertices[j]();
   real_t length = 0.;

   for (int64_t k = 0; k < spaceDim; k++)
   {
      length += (vi[k]-vj[k])*(vi[k]-vj[k]);
   }

   return sqrt(length);
}

// static method
void Mesh::GetElementArrayEdgeTable(const Array<Element*> &elem_array,
                                    const DSTable &v_to_v, Table &el_to_edge)
{
   el_to_edge.MakeI(elem_array.Size());
   for (int64_t i = 0; i < elem_array.Size(); i++)
   {
      el_to_edge.AddColumnsInRow(i, elem_array[i]->GetNEdges());
   }
   el_to_edge.MakeJ();
   for (int64_t i = 0; i < elem_array.Size(); i++)
   {
      const int64_t *v = elem_array[i]->GetVertices();
      const int64_t ne = elem_array[i]->GetNEdges();
      for (int64_t j = 0; j < ne; j++)
      {
         const int64_t *e = elem_array[i]->GetEdgeVertices(j);
         el_to_edge.AddConnection(i, v_to_v(v[e[0]], v[e[1]]));
      }
   }
   el_to_edge.ShiftUpI();
}

void Mesh::GetVertexToVertexTable(DSTable &v_to_v) const
{
   if (edge_vertex)
   {
      for (int64_t i = 0; i < edge_vertex->Size(); i++)
      {
         const int64_t *v = edge_vertex->GetRow(i);
         v_to_v.Push(v[0], v[1]);
      }
   }
   else
   {
      for (int64_t i = 0; i < NumOfElements; i++)
      {
         const int64_t *v = elements[i]->GetVertices();
         const int64_t ne = elements[i]->GetNEdges();
         for (int64_t j = 0; j < ne; j++)
         {
            const int64_t *e = elements[i]->GetEdgeVertices(j);
            v_to_v.Push(v[e[0]], v[e[1]]);
         }
      }
   }
}

int64_t Mesh::GetElementToEdgeTable(Table &e_to_f)
{
   int64_t i, NumberOfEdges;

   DSTable v_to_v(NumOfVertices);
   GetVertexToVertexTable(v_to_v);

   NumberOfEdges = v_to_v.NumberOfEntries();

   // Fill the element to edge table
   GetElementArrayEdgeTable(elements, v_to_v, e_to_f);

   if (Dim == 2)
   {
      // Initialize the indices for the boundary elements.
      be_to_face.SetSize(NumOfBdrElements);
      for (i = 0; i < NumOfBdrElements; i++)
      {
         const int64_t *v = boundary[i]->GetVertices();
         be_to_face[i] = v_to_v(v[0], v[1]);
      }
   }
   else if (Dim == 3)
   {
      if (bel_to_edge == NULL)
      {
         bel_to_edge = new Table;
      }
      GetElementArrayEdgeTable(boundary, v_to_v, *bel_to_edge);
   }
   else
   {
      mfem_error("1D GetElementToEdgeTable is not yet implemented.");
   }

   // Return the number of edges
   return NumberOfEdges;
}

const Table & Mesh::ElementToElementTable()
{
   if (el_to_el)
   {
      return *el_to_el;
   }

   // Note that, for ParNCMeshes, faces_info will contain also the ghost faces
   MFEM_ASSERT(faces_info.Size() >= GetNumFaces(), "faces were not generated!");

   Array<Connection> conn;
   conn.Reserve(2*faces_info.Size());

   for (int64_t i = 0; i < faces_info.Size(); i++)
   {
      const FaceInfo &fi = faces_info[i];
      if (fi.Elem2No >= 0)
      {
         conn.Append(Connection(fi.Elem1No, fi.Elem2No));
         conn.Append(Connection(fi.Elem2No, fi.Elem1No));
      }
      else if (fi.Elem2Inf >= 0)
      {
         int64_t nbr_elem_idx = NumOfElements - 1 - fi.Elem2No;
         conn.Append(Connection(fi.Elem1No, nbr_elem_idx));
         conn.Append(Connection(nbr_elem_idx, fi.Elem1No));
      }
   }

   conn.Sort();
   conn.Unique();
   el_to_el = new Table(NumOfElements, conn);

   return *el_to_el;
}

const Table & Mesh::ElementToFaceTable() const
{
   if (el_to_face == NULL)
   {
      mfem_error("Mesh::ElementToFaceTable()");
   }
   return *el_to_face;
}

const Table & Mesh::ElementToEdgeTable() const
{
   if (el_to_edge == NULL)
   {
      mfem_error("Mesh::ElementToEdgeTable()");
   }
   return *el_to_edge;
}

void Mesh::AddPointFaceElement(int64_t lf, int64_t gf, int64_t el)
{
   if (faces[gf] == NULL)  // this will be elem1
   {
      faces[gf] = new Point(&gf);
      faces_info[gf].Elem1No  = el;
      faces_info[gf].Elem1Inf = 64 * lf; // face lf with orientation 0
      faces_info[gf].Elem2No  = -1; // in case there's no other side
      faces_info[gf].Elem2Inf = -1; // face is not shared
   }
   else  //  this will be elem2
   {
      /* WARNING: Without the following check the mesh faces_info data structure
         may contain unreliable data. Normally, the order in which elements are
         processed could swap which elements appear as Elem1No and Elem2No. In
         branched meshes, where more than two elements can meet at a given node,
         the indices stored in Elem1No and Elem2No will be the first and last,
         respectively, elements found which touch a given node. This can lead to
         inconsistencies in any algorithms which rely on this data structure. To
         properly support branched meshes this data structure should be extended
         to support multiple elements per face. */
      /*
      MFEM_VERIFY(faces_info[gf].Elem2No < 0, "Invalid mesh topology. "
                  "Interior point found connecting 1D elements "
                  << faces_info[gf].Elem1No << ", " << faces_info[gf].Elem2No
                  << " and " << el << ".");
      */
      faces_info[gf].Elem2No  = el;
      faces_info[gf].Elem2Inf = 64 * lf + 1;
   }
}

void Mesh::AddSegmentFaceElement(int64_t lf, int64_t gf, int64_t el, int64_t v0, int64_t v1)
{
   if (faces[gf] == NULL)  // this will be elem1
   {
      faces[gf] = new Segment(v0, v1);
      faces_info[gf].Elem1No  = el;
      faces_info[gf].Elem1Inf = 64 * lf; // face lf with orientation 0
      faces_info[gf].Elem2No  = -1; // in case there's no other side
      faces_info[gf].Elem2Inf = -1; // face is not shared
   }
   else  //  this will be elem2
   {
      MFEM_VERIFY(faces_info[gf].Elem2No < 0, "Invalid mesh topology.  "
                  "Interior edge found between 2D elements "
                  << faces_info[gf].Elem1No << ", " << faces_info[gf].Elem2No
                  << " and " << el << ".");
      int64_t *v = faces[gf]->GetVertices();
      faces_info[gf].Elem2No  = el;
      if (v[1] == v0 && v[0] == v1)
      {
         faces_info[gf].Elem2Inf = 64 * lf + 1;
      }
      else if (v[0] == v0 && v[1] == v1)
      {
         // Temporarily allow even edge orientations: see the remark in
         // AddTriangleFaceElement().
         // Also, in a non-orientable surface mesh, the orientation will be even
         // for edges that connect elements with opposite orientations.
         faces_info[gf].Elem2Inf = 64 * lf;
      }
      else
      {
         MFEM_ABORT("internal error");
      }
   }
}

void Mesh::AddTriangleFaceElement(int64_t lf, int64_t gf, int64_t el,
                                  int64_t v0, int64_t v1, int64_t v2)
{
   if (faces[gf] == NULL)  // this will be elem1
   {
      faces[gf] = new Triangle(v0, v1, v2);
      faces_info[gf].Elem1No  = el;
      faces_info[gf].Elem1Inf = 64 * lf; // face lf with orientation 0
      faces_info[gf].Elem2No  = -1; // in case there's no other side
      faces_info[gf].Elem2Inf = -1; // face is not shared
   }
   else  //  this will be elem2
   {
      MFEM_VERIFY(faces_info[gf].Elem2No < 0, "Invalid mesh topology.  "
                  "Interior triangular face found connecting elements "
                  << faces_info[gf].Elem1No << ", " << faces_info[gf].Elem2No
                  << " and " << el << ".");
      int64_t orientation, vv[3] = { v0, v1, v2 };
      orientation = GetTriOrientation(faces[gf]->GetVertices(), vv);
      // In a valid mesh, we should have (orientation % 2 != 0), however, if
      // one of the adjacent elements has wrong orientation, both face
      // orientations can be even, until the element orientations are fixed.
      // MFEM_ASSERT(orientation % 2 != 0, "");
      faces_info[gf].Elem2No  = el;
      faces_info[gf].Elem2Inf = 64 * lf + orientation;
   }
}

void Mesh::AddQuadFaceElement(int64_t lf, int64_t gf, int64_t el,
                              int64_t v0, int64_t v1, int64_t v2, int64_t v3)
{
   if (faces_info[gf].Elem1No < 0)  // this will be elem1
   {
      faces[gf] = new Quadrilateral(v0, v1, v2, v3);
      faces_info[gf].Elem1No  = el;
      faces_info[gf].Elem1Inf = 64 * lf; // face lf with orientation 0
      faces_info[gf].Elem2No  = -1; // in case there's no other side
      faces_info[gf].Elem2Inf = -1; // face is not shared
   }
   else  //  this will be elem2
   {
      MFEM_VERIFY(faces_info[gf].Elem2No < 0, "Invalid mesh topology.  "
                  "Interior quadrilateral face found connecting elements "
                  << faces_info[gf].Elem1No << ", " << faces_info[gf].Elem2No
                  << " and " << el << ".");
      int64_t vv[4] = { v0, v1, v2, v3 };
      int64_t oo = GetQuadOrientation(faces[gf]->GetVertices(), vv);
      // Temporarily allow even face orientations: see the remark in
      // AddTriangleFaceElement().
      // MFEM_ASSERT(oo % 2 != 0, "");
      faces_info[gf].Elem2No  = el;
      faces_info[gf].Elem2Inf = 64 * lf + oo;
   }
}

void Mesh::GenerateFaces()
{
   int64_t nfaces = GetNumFaces();
   for (auto &f : faces)
   {
      FreeElement(f);
   }

   // (re)generate the interior faces and the info for them
   faces.SetSize(nfaces);
   faces_info.SetSize(nfaces);
   for (int64_t i = 0; i < nfaces; ++i)
   {
      faces[i] = NULL;
      faces_info[i].Elem1No = -1;
      faces_info[i].NCFace = -1;
   }

   Array<int64_t> v;
   for (int64_t i = 0; i < NumOfElements; ++i)
   {
      elements[i]->GetVertices(v);
      if (Dim == 1)
      {
         AddPointFaceElement(0, v[0], i);
         AddPointFaceElement(1, v[1], i);
      }
      else if (Dim == 2)
      {
         const int64_t * const ef = el_to_edge->GetRow(i);
         const int64_t ne = elements[i]->GetNEdges();
         for (int64_t j = 0; j < ne; j++)
         {
            const int64_t *e = elements[i]->GetEdgeVertices(j);
            AddSegmentFaceElement(j, ef[j], i, v[e[0]], v[e[1]]);
         }
      }
      else
      {
         const int64_t * const ef = el_to_face->GetRow(i);
         switch (GetElementType(i))
         {
            case Element::TETRAHEDRON:
            {
               for (int64_t j = 0; j < 4; j++)
               {
                  const int64_t *fv = tet_t::FaceVert[j];
                  AddTriangleFaceElement(j, ef[j], i,
                                         v[fv[0]], v[fv[1]], v[fv[2]]);
               }
               break;
            }
            case Element::WEDGE:
            {
               for (int64_t j = 0; j < 2; j++)
               {
                  const int64_t *fv = pri_t::FaceVert[j];
                  AddTriangleFaceElement(j, ef[j], i,
                                         v[fv[0]], v[fv[1]], v[fv[2]]);
               }
               for (int64_t j = 2; j < 5; j++)
               {
                  const int64_t *fv = pri_t::FaceVert[j];
                  AddQuadFaceElement(j, ef[j], i,
                                     v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
               }
               break;
            }
            case Element::PYRAMID:
            {
               for (int64_t j = 0; j < 1; j++)
               {
                  const int64_t *fv = pyr_t::FaceVert[j];
                  AddQuadFaceElement(j, ef[j], i,
                                     v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
               }
               for (int64_t j = 1; j < 5; j++)
               {
                  const int64_t *fv = pyr_t::FaceVert[j];
                  AddTriangleFaceElement(j, ef[j], i,
                                         v[fv[0]], v[fv[1]], v[fv[2]]);
               }
               break;
            }
            case Element::HEXAHEDRON:
            {
               for (int64_t j = 0; j < 6; j++)
               {
                  const int64_t *fv = hex_t::FaceVert[j];
                  AddQuadFaceElement(j, ef[j], i,
                                     v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
               }
               break;
            }
            default:
               MFEM_ABORT("Unexpected type of Element.");
         }
      }
   }
}

STable3D *Mesh::GetFacesTable()
{
   STable3D *faces_tbl = new STable3D(NumOfVertices);
   for (int64_t i = 0; i < NumOfElements; i++)
   {
      const int64_t *v = elements[i]->GetVertices();
      switch (GetElementType(i))
      {
         case Element::TETRAHEDRON:
         {
            for (int64_t j = 0; j < 4; j++)
            {
               const int64_t *fv = tet_t::FaceVert[j];
               faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]);
            }
            break;
         }
         case Element::PYRAMID:
         {
            for (int64_t j = 0; j < 1; j++)
            {
               const int64_t *fv = pyr_t::FaceVert[j];
               faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
            }
            for (int64_t j = 1; j < 5; j++)
            {
               const int64_t *fv = pyr_t::FaceVert[j];
               faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]);
            }
            break;
         }
         case Element::WEDGE:
         {
            for (int64_t j = 0; j < 2; j++)
            {
               const int64_t *fv = pri_t::FaceVert[j];
               faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]);
            }
            for (int64_t j = 2; j < 5; j++)
            {
               const int64_t *fv = pri_t::FaceVert[j];
               faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
            }
            break;
         }
         case Element::HEXAHEDRON:
         {
            // find the face by the vertices with the smallest 3 numbers
            // z = 0, y = 0, x = 1, y = 1, x = 0, z = 1
            for (int64_t j = 0; j < 6; j++)
            {
               const int64_t *fv = hex_t::FaceVert[j];
               faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]);
            }
            break;
         }
         default:
            MFEM_ABORT("Unexpected type of Element: " << GetElementType(i));
      }
   }
   return faces_tbl;
}

STable3D *Mesh::GetElementToFaceTable(int64_t ret_ftbl)
{
   Array<int64_t> v;
   STable3D *faces_tbl;

   if (el_to_face != NULL)
   {
      delete el_to_face;
   }
   el_to_face = new Table(NumOfElements, 6);  // must be 6 for hexahedra
   faces_tbl = new STable3D(NumOfVertices);
   for (int64_t i = 0; i < NumOfElements; i++)
   {
      elements[i]->GetVertices(v);
      switch (GetElementType(i))
      {
         case Element::TETRAHEDRON:
         {
            for (int64_t j = 0; j < 4; j++)
            {
               const int64_t *fv = tet_t::FaceVert[j];
               el_to_face->Push(
                  i, faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]));
            }
            break;
         }
         case Element::WEDGE:
         {
            for (int64_t j = 0; j < 2; j++)
            {
               const int64_t *fv = pri_t::FaceVert[j];
               el_to_face->Push(
                  i, faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]));
            }
            for (int64_t j = 2; j < 5; j++)
            {
               const int64_t *fv = pri_t::FaceVert[j];
               el_to_face->Push(
                  i, faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]));
            }
            break;
         }
         case Element::PYRAMID:
         {
            for (int64_t j = 0; j < 1; j++)
            {
               const int64_t *fv = pyr_t::FaceVert[j];
               el_to_face->Push(
                  i, faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]));
            }
            for (int64_t j = 1; j < 5; j++)
            {
               const int64_t *fv = pyr_t::FaceVert[j];
               el_to_face->Push(
                  i, faces_tbl->Push(v[fv[0]], v[fv[1]], v[fv[2]]));
            }
            break;
         }
         case Element::HEXAHEDRON:
         {
            // find the face by the vertices with the smallest 3 numbers
            // z = 0, y = 0, x = 1, y = 1, x = 0, z = 1
            for (int64_t j = 0; j < 6; j++)
            {
               const int64_t *fv = hex_t::FaceVert[j];
               el_to_face->Push(
                  i, faces_tbl->Push4(v[fv[0]], v[fv[1]], v[fv[2]], v[fv[3]]));
            }
            break;
         }
         default:
            MFEM_ABORT("Unexpected type of Element.");
      }
   }
   el_to_face->Finalize();
   NumOfFaces = faces_tbl->NumberOfElements();
   be_to_face.SetSize(NumOfBdrElements);

   for (int64_t i = 0; i < NumOfBdrElements; i++)
   {
      boundary[i]->GetVertices(v);
      switch (GetBdrElementType(i))
      {
         case Element::TRIANGLE:
         {
            be_to_face[i] = (*faces_tbl)(v[0], v[1], v[2]);
            break;
         }
         case Element::QUADRILATERAL:
         {
            be_to_face[i] = (*faces_tbl)(v[0], v[1], v[2], v[3]);
            break;
         }
         default:
            MFEM_ABORT("Unexpected type of boundary Element.");
      }
   }

   if (ret_ftbl)
   {
      return faces_tbl;
   }
   delete faces_tbl;
   return NULL;
}

// shift cyclically 3 integers so that the smallest is first
static inline
void Rotate3(int64_t &a, int64_t &b, int64_t &c)
{
   if (a < b)
   {
      if (a > c)
      {
         ShiftRight(a, b, c);
      }
   }
   else
   {
      if (b < c)
      {
         ShiftRight(c, b, a);
      }
      else
      {
         ShiftRight(a, b, c);
      }
   }
}

void FindPartitioningComponents(Table &elem_elem,
                                const Array<int64_t> &partitioning,
                                Array<int64_t> &component,
                                Array<int64_t> &num_comp);

int64_t *Mesh::GeneratePartitioning(int64_t nparts, int64_t part_method)
{
#ifdef MFEM_USE_METIS

   int64_t print_messages = 1;

   int64_t i, *partitioning;

   ElementToElementTable();

   partitioning = new int64_t[NumOfElements];

   if (nparts == 1)
   {
      for (i = 0; i < NumOfElements; i++)
      {
         partitioning[i] = 0;
      }
   }
   else if (NumOfElements <= nparts)
   {
      for (i = 0; i < NumOfElements; i++)
      {
         partitioning[i] = i;
      }
   }
   else
   {
      idx_t *I, *J, n;
#ifndef MFEM_USE_METIS_5
      idx_t wgtflag = 0;
      idx_t numflag = 0;
      idx_t options[5];
#else
      idx_t ncon = 1;
      idx_t errflag;
      idx_t options[40];
#endif
      idx_t edgecut;

      // In case METIS have been compiled with 64bit indices
      bool freedata = false;
      idx_t mparts = (idx_t) nparts;
      idx_t *mpartitioning;

      n = NumOfElements;
      if (sizeof(idx_t) == sizeof(int64_t))
      {
         I = (idx_t*) el_to_el->GetI();
         J = (idx_t*) el_to_el->GetJ();
         mpartitioning = (idx_t*) partitioning;
      }
      else
      {
         int64_t *iI = el_to_el->GetI();
         int64_t *iJ = el_to_el->GetJ();
         int64_t m = iI[n];
         I = new idx_t[n+1];
         J = new idx_t[m];
         for (int64_t k = 0; k < n+1; k++) { I[k] = iI[k]; }
         for (int64_t k = 0; k < m; k++) { J[k] = iJ[k]; }
         mpartitioning = new idx_t[n];
         freedata = true;
      }
#ifndef MFEM_USE_METIS_5
      options[0] = 0;
#else
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_CONTIG] = 1; // set METIS_OPTION_CONTIG
      // If the mesh is disconnected, disable METIS_OPTION_CONTIG.
      {
         Array<int64_t> part(partitioning, NumOfElements);
         part = 0; // single part for the whole mesh
         Array<int64_t> component; // size will be set to num. elem.
         Array<int64_t> num_comp;  // size will be set to num. parts (1)
         FindPartitioningComponents(*el_to_el, part, component, num_comp);
         if (num_comp[0] > 1) { options[METIS_OPTION_CONTIG] = 0; }
      }
#endif

      // Sort the neighbor lists
      if (part_method >= 0 && part_method <= 2)
      {
         for (i = 0; i < n; i++)
         {
            // Sort in increasing order.
            // std::sort(J+I[i], J+I[i+1]);

            // Sort in decreasing order, as in previous versions of MFEM.
            std::sort(J+I[i], J+I[i+1], std::greater<idx_t>());
         }
      }

      // This function should be used to partition a graph into a small
      // number of partitions (less than 8).
      if (part_method == 0 || part_method == 3)
      {
#ifndef MFEM_USE_METIS_5
         METIS_PartGraphRecursive(&n,
                                  I,
                                  J,
                                  NULL,
                                  NULL,
                                  &wgtflag,
                                  &numflag,
                                  &mparts,
                                  options,
                                  &edgecut,
                                  mpartitioning);
#else
         errflag = METIS_PartGraphRecursive(&n,
                                            &ncon,
                                            I,
                                            J,
                                            NULL,
                                            NULL,
                                            NULL,
                                            &mparts,
                                            NULL,
                                            NULL,
                                            options,
                                            &edgecut,
                                            mpartitioning);
         if (errflag != 1)
         {
            mfem_error("Mesh::GeneratePartitioning: "
                       " error in METIS_PartGraphRecursive!");
         }
#endif
      }

      // This function should be used to partition a graph into a large
      // number of partitions (greater than 8).
      if (part_method == 1 || part_method == 4)
      {
#ifndef MFEM_USE_METIS_5
         METIS_PartGraphKway(&n,
                             I,
                             J,
                             NULL,
                             NULL,
                             &wgtflag,
                             &numflag,
                             &mparts,
                             options,
                             &edgecut,
                             mpartitioning);
#else
         errflag = METIS_PartGraphKway(&n,
                                       &ncon,
                                       I,
                                       J,
                                       NULL,
                                       NULL,
                                       NULL,
                                       &mparts,
                                       NULL,
                                       NULL,
                                       options,
                                       &edgecut,
                                       mpartitioning);
         if (errflag != 1)
         {
            mfem_error("Mesh::GeneratePartitioning: "
                       " error in METIS_PartGraphKway!");
         }
#endif
      }

      // The objective of this partitioning is to minimize the total
      // communication volume
      if (part_method == 2 || part_method == 5)
      {
#ifndef MFEM_USE_METIS_5
         METIS_PartGraphVKway(&n,
                              I,
                              J,
                              NULL,
                              NULL,
                              &wgtflag,
                              &numflag,
                              &mparts,
                              options,
                              &edgecut,
                              mpartitioning);
#else
         options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
         errflag = METIS_PartGraphKway(&n,
                                       &ncon,
                                       I,
                                       J,
                                       NULL,
                                       NULL,
                                       NULL,
                                       &mparts,
                                       NULL,
                                       NULL,
                                       options,
                                       &edgecut,
                                       mpartitioning);
         if (errflag != 1)
         {
            mfem_error("Mesh::GeneratePartitioning: "
                       " error in METIS_PartGraphKway!");
         }
#endif
      }

#ifdef MFEM_DEBUG
      if (print_messages)
      {
         mfem::out << "Mesh::GeneratePartitioning(...): edgecut = "
                   << edgecut << endl;
      }
#endif
      nparts = (int64_t) mparts;
      if (mpartitioning != (idx_t*)partitioning)
      {
         for (int64_t k = 0; k<NumOfElements; k++)
         {
            partitioning[k] = mpartitioning[k];
         }
      }
      if (freedata)
      {
         delete[] I;
         delete[] J;
         delete[] mpartitioning;
      }
   }

   delete el_to_el;
   el_to_el = NULL;

   // Check for empty partitionings (a "feature" in METIS)
   if (nparts > 1 && NumOfElements > nparts)
   {
      Array< Pair<int64_t,int64_t> > psize(nparts);
      int64_t empty_parts;

      // Count how many elements are in each partition, and store the result in
      // psize, where psize[i].one is the number of elements, and psize[i].two
      // is partition index. Keep track of the number of empty parts.
      auto count_partition_elements = [&]()
      {
         for (i = 0; i < nparts; i++)
         {
            psize[i].one = 0;
            psize[i].two = i;
         }

         for (i = 0; i < NumOfElements; i++)
         {
            psize[partitioning[i]].one++;
         }

         empty_parts = 0;
         for (i = 0; i < nparts; i++)
         {
            if (psize[i].one == 0) { empty_parts++; }
         }
      };

      count_partition_elements();

      // This code just split the largest partitionings in two.
      // Do we need to replace it with something better?
      while (empty_parts)
      {
         if (print_messages)
         {
            mfem::err << "Mesh::GeneratePartitioning(...): METIS returned "
                      << empty_parts << " empty parts!"
                      << " Applying a simple fix ..." << endl;
         }

         SortPairs<int64_t,int64_t>(psize, nparts);

         for (i = nparts-1; i > nparts-1-empty_parts; i--)
         {
            psize[i].one /= 2;
         }

         for (int64_t j = 0; j < NumOfElements; j++)
         {
            for (i = nparts-1; i > nparts-1-empty_parts; i--)
            {
               if (psize[i].one == 0 || partitioning[j] != psize[i].two)
               {
                  continue;
               }
               else
               {
                  partitioning[j] = psize[nparts-1-i].two;
                  psize[i].one--;
               }
            }
         }

         // Check for empty partitionings again
         count_partition_elements();
      }
   }

   return partitioning;

#else
   MFEM_ABORT("MFEM was compiled without Metis.");
#endif
}

/* required: 0 <= partitioning[i] < num_part */
void FindPartitioningComponents(Table &elem_elem,
                                const Array<int64_t> &partitioning,
                                Array<int64_t> &component,
                                Array<int64_t> &num_comp)
{
   int64_t i, j, k;
   int64_t num_elem, *i_elem_elem, *j_elem_elem;

   num_elem    = elem_elem.Size();
   i_elem_elem = elem_elem.GetI();
   j_elem_elem = elem_elem.GetJ();

   component.SetSize(num_elem);

   Array<int64_t> elem_stack(num_elem);
   int64_t stack_p, stack_top_p, elem;
   int64_t num_part;

   num_part = -1;
   for (i = 0; i < num_elem; i++)
   {
      if (partitioning[i] > num_part)
      {
         num_part = partitioning[i];
      }
      component[i] = -1;
   }
   num_part++;

   num_comp.SetSize(num_part);
   for (i = 0; i < num_part; i++)
   {
      num_comp[i] = 0;
   }

   stack_p = 0;
   stack_top_p = 0;  // points to the first unused element in the stack
   for (elem = 0; elem < num_elem; elem++)
   {
      if (component[elem] >= 0)
      {
         continue;
      }

      component[elem] = num_comp[partitioning[elem]]++;

      elem_stack[stack_top_p++] = elem;

      for ( ; stack_p < stack_top_p; stack_p++)
      {
         i = elem_stack[stack_p];
         for (j = i_elem_elem[i]; j < i_elem_elem[i+1]; j++)
         {
            k = j_elem_elem[j];
            if (partitioning[k] == partitioning[i])
            {
               if (component[k] < 0)
               {
                  component[k] = component[i];
                  elem_stack[stack_top_p++] = k;
               }
               else if (component[k] != component[i])
               {
                  mfem_error("FindPartitioningComponents");
               }
            }
         }
      }
   }
}

void Mesh::CheckPartitioning(int64_t *partitioning_)
{
   int64_t i, n_empty, n_mcomp;
   Array<int64_t> component, num_comp;
   const Array<int64_t> partitioning(partitioning_, GetNE());

   ElementToElementTable();

   FindPartitioningComponents(*el_to_el, partitioning, component, num_comp);

   n_empty = n_mcomp = 0;
   for (i = 0; i < num_comp.Size(); i++)
      if (num_comp[i] == 0)
      {
         n_empty++;
      }
      else if (num_comp[i] > 1)
      {
         n_mcomp++;
      }

   if (n_empty > 0)
   {
      mfem::out << "Mesh::CheckPartitioning(...) :\n"
                << "The following subdomains are empty :\n";
      for (i = 0; i < num_comp.Size(); i++)
         if (num_comp[i] == 0)
         {
            mfem::out << ' ' << i;
         }
      mfem::out << endl;
   }
   if (n_mcomp > 0)
   {
      mfem::out << "Mesh::CheckPartitioning(...) :\n"
                << "The following subdomains are NOT connected :\n";
      for (i = 0; i < num_comp.Size(); i++)
         if (num_comp[i] > 1)
         {
            mfem::out << ' ' << i;
         }
      mfem::out << endl;
   }
   if (n_empty == 0 && n_mcomp == 0)
      mfem::out << "Mesh::CheckPartitioning(...) : "
                "All subdomains are connected." << endl;

   if (el_to_el)
   {
      delete el_to_el;
   }
   el_to_el = NULL;
}

// compute the coefficients of the polynomial in t:
//   c(0)+c(1)*t+...+c(d)*t^d = det(A+t*B)
// where A, B are (d x d), d=2,3
void DetOfLinComb(const DenseMatrix &A, const DenseMatrix &B, Vector &c)
{
   const real_t *a = A.Data();
   const real_t *b = B.Data();

   c.SetSize(A.Width()+1);
   switch (A.Width())
   {
      case 2:
      {
         // det(A+t*B) = |a0 a2|   / |a0 b2| + |b0 a2| \       |b0 b2|
         //              |a1 a3| + \ |a1 b3|   |b1 a3| / * t + |b1 b3| * t^2
         c(0) = a[0]*a[3]-a[1]*a[2];
         c(1) = a[0]*b[3]-a[1]*b[2]+b[0]*a[3]-b[1]*a[2];
         c(2) = b[0]*b[3]-b[1]*b[2];
      }
      break;

      case 3:
      {
         /*              |a0 a3 a6|
          * det(A+t*B) = |a1 a4 a7| +
          *              |a2 a5 a8|

          *     /  |b0 a3 a6|   |a0 b3 a6|   |a0 a3 b6| \
          *   + |  |b1 a4 a7| + |a1 b4 a7| + |a1 a4 b7| | * t +
          *     \  |b2 a5 a8|   |a2 b5 a8|   |a2 a5 b8| /

          *     /  |a0 b3 b6|   |b0 a3 b6|   |b0 b3 a6| \
          *   + |  |a1 b4 b7| + |b1 a4 b7| + |b1 b4 a7| | * t^2 +
          *     \  |a2 b5 b8|   |b2 a5 b8|   |b2 b5 a8| /

          *     |b0 b3 b6|
          *   + |b1 b4 b7| * t^3
          *     |b2 b5 b8|       */
         c(0) = (a[0] * (a[4] * a[8] - a[5] * a[7]) +
                 a[1] * (a[5] * a[6] - a[3] * a[8]) +
                 a[2] * (a[3] * a[7] - a[4] * a[6]));

         c(1) = (b[0] * (a[4] * a[8] - a[5] * a[7]) +
                 b[1] * (a[5] * a[6] - a[3] * a[8]) +
                 b[2] * (a[3] * a[7] - a[4] * a[6]) +

                 a[0] * (b[4] * a[8] - b[5] * a[7]) +
                 a[1] * (b[5] * a[6] - b[3] * a[8]) +
                 a[2] * (b[3] * a[7] - b[4] * a[6]) +

                 a[0] * (a[4] * b[8] - a[5] * b[7]) +
                 a[1] * (a[5] * b[6] - a[3] * b[8]) +
                 a[2] * (a[3] * b[7] - a[4] * b[6]));

         c(2) = (a[0] * (b[4] * b[8] - b[5] * b[7]) +
                 a[1] * (b[5] * b[6] - b[3] * b[8]) +
                 a[2] * (b[3] * b[7] - b[4] * b[6]) +

                 b[0] * (a[4] * b[8] - a[5] * b[7]) +
                 b[1] * (a[5] * b[6] - a[3] * b[8]) +
                 b[2] * (a[3] * b[7] - a[4] * b[6]) +

                 b[0] * (b[4] * a[8] - b[5] * a[7]) +
                 b[1] * (b[5] * a[6] - b[3] * a[8]) +
                 b[2] * (b[3] * a[7] - b[4] * a[6]));

         c(3) = (b[0] * (b[4] * b[8] - b[5] * b[7]) +
                 b[1] * (b[5] * b[6] - b[3] * b[8]) +
                 b[2] * (b[3] * b[7] - b[4] * b[6]));
      }
      break;

      default:
         mfem_error("DetOfLinComb(...)");
   }
}

// compute the real roots of
//   z(0)+z(1)*x+...+z(d)*x^d = 0,  d=2,3;
// the roots are returned in x, sorted in increasing order;
// it is assumed that x is at least of size d;
// return the number of roots counting multiplicity;
// return -1 if all z(i) are 0.
int64_t FindRoots(const Vector &z, Vector &x)
{
   int64_t d = z.Size()-1;
   if (d > 3 || d < 0)
   {
      mfem_error("FindRoots(...)");
   }

   while (z(d) == 0.0)
   {
      if (d == 0)
      {
         return (-1);
      }
      d--;
   }
   switch (d)
   {
      case 0:
      {
         return 0;
      }

      case 1:
      {
         x(0) = -z(0)/z(1);
         return 1;
      }

      case 2:
      {
         real_t a = z(2), b = z(1), c = z(0);
         real_t D = b*b-4*a*c;
         if (D < 0.0)
         {
            return 0;
         }
         if (D == 0.0)
         {
            x(0) = x(1) = -0.5 * b / a;
            return 2; // root with multiplicity 2
         }
         if (b == 0.0)
         {
            x(0) = -(x(1) = fabs(0.5 * sqrt(D) / a));
            return 2;
         }
         else
         {
            real_t t;
            if (b > 0.0)
            {
               t = -0.5 * (b + sqrt(D));
            }
            else
            {
               t = -0.5 * (b - sqrt(D));
            }
            x(0) = t / a;
            x(1) = c / t;
            if (x(0) > x(1))
            {
               Swap<real_t>(x(0), x(1));
            }
            return 2;
         }
      }

      case 3:
      {
         real_t a = z(2)/z(3), b = z(1)/z(3), c = z(0)/z(3);

         // find the real roots of x^3 + a x^2 + b x + c = 0
         real_t Q = (a * a - 3 * b) / 9;
         real_t R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
         real_t Q3 = Q * Q * Q;
         real_t R2 = R * R;

         if (R2 == Q3)
         {
            if (Q == 0)
            {
               x(0) = x(1) = x(2) = - a / 3;
            }
            else
            {
               real_t sqrtQ = sqrt(Q);

               if (R > 0)
               {
                  x(0) = -2 * sqrtQ - a / 3;
                  x(1) = x(2) = sqrtQ - a / 3;
               }
               else
               {
                  x(0) = x(1) = - sqrtQ - a / 3;
                  x(2) = 2 * sqrtQ - a / 3;
               }
            }
            return 3;
         }
         else if (R2 < Q3)
         {
            real_t theta = acos(R / sqrt(Q3));
            real_t A = -2 * sqrt(Q);
            real_t x0, x1, x2;
            x0 = A * cos(theta / 3) - a / 3;
            x1 = A * cos((theta + 2.0 * M_PI) / 3) - a / 3;
            x2 = A * cos((theta - 2.0 * M_PI) / 3) - a / 3;

            /* Sort x0, x1, x2 */
            if (x0 > x1)
            {
               Swap<real_t>(x0, x1);
            }
            if (x1 > x2)
            {
               Swap<real_t>(x1, x2);
               if (x0 > x1)
               {
                  Swap<real_t>(x0, x1);
               }
            }
            x(0) = x0;
            x(1) = x1;
            x(2) = x2;
            return 3;
         }
         else
         {
            real_t A;
            if (R >= 0.0)
            {
               A = -pow(sqrt(R2 - Q3) + R, 1.0/3.0);
            }
            else
            {
               A =  pow(sqrt(R2 - Q3) - R, 1.0/3.0);
            }
            x(0) = A + Q / A - a / 3;
            return 1;
         }
      }
   }
   return 0;
}

void FindTMax(Vector &c, Vector &x, real_t &tmax,
              const real_t factor, const int64_t Dim)
{
   const real_t c0 = c(0);
   c(0) = c0 * (1.0 - pow(factor, -Dim));
   int64_t nr = FindRoots(c, x);
   for (int64_t j = 0; j < nr; j++)
   {
      if (x(j) > tmax)
      {
         break;
      }
      if (x(j) >= 0.0)
      {
         tmax = x(j);
         break;
      }
   }
   c(0) = c0 * (1.0 - pow(factor, Dim));
   nr = FindRoots(c, x);
   for (int64_t j = 0; j < nr; j++)
   {
      if (x(j) > tmax)
      {
         break;
      }
      if (x(j) >= 0.0)
      {
         tmax = x(j);
         break;
      }
   }
}

void Mesh::Swap(Mesh& other, bool non_geometry)
{
   mfem::Swap(Dim, other.Dim);
   mfem::Swap(spaceDim, other.spaceDim);

   mfem::Swap(NumOfVertices, other.NumOfVertices);
   mfem::Swap(NumOfElements, other.NumOfElements);
   mfem::Swap(NumOfBdrElements, other.NumOfBdrElements);
   mfem::Swap(NumOfEdges, other.NumOfEdges);
   mfem::Swap(NumOfFaces, other.NumOfFaces);

   mfem::Swap(meshgen, other.meshgen);
   mfem::Swap(mesh_geoms, other.mesh_geoms);

   mfem::Swap(elements, other.elements);
   mfem::Swap(vertices, other.vertices);
   mfem::Swap(boundary, other.boundary);
   mfem::Swap(faces, other.faces);
   mfem::Swap(faces_info, other.faces_info);
   mfem::Swap(nc_faces_info, other.nc_faces_info);
   mfem::Swap(nbInteriorFaces, other.nbInteriorFaces);
   mfem::Swap(nbBoundaryFaces, other.nbBoundaryFaces);

   mfem::Swap(el_to_edge, other.el_to_edge);
   mfem::Swap(el_to_face, other.el_to_face);
   mfem::Swap(el_to_el, other.el_to_el);
   mfem::Swap(bel_to_edge, other.bel_to_edge);
   mfem::Swap(be_to_face, other.be_to_face);
   mfem::Swap(face_edge, other.face_edge);
   mfem::Swap(face_to_elem, other.face_to_elem);
   mfem::Swap(edge_vertex, other.edge_vertex);

   mfem::Swap(attributes, other.attributes);
   mfem::Swap(bdr_attributes, other.bdr_attributes);
}

void Mesh::GetElementData(const Array<Element*> &elem_array, int64_t geom,
                          Array<int64_t> &elem_vtx, Array<int64_t> &attr) const
{
   // protected method
   const int64_t nv = Geometry::NumVerts[geom];
   int64_t num_elems = 0;
   for (int64_t i = 0; i < elem_array.Size(); i++)
   {
      if (elem_array[i]->GetGeometryType() == geom)
      {
         num_elems++;
      }
   }
   elem_vtx.SetSize(nv*num_elems);
   attr.SetSize(num_elems);
   elem_vtx.SetSize(0);
   attr.SetSize(0);
   for (int64_t i = 0; i < elem_array.Size(); i++)
   {
      Element *el = elem_array[i];
      if (el->GetGeometryType() != geom) { continue; }

      Array<int64_t> loc_vtx(el->GetVertices(), nv);
      elem_vtx.Append(loc_vtx);
      attr.Append(el->GetAttribute());
   }
}

void Mesh::Printer(std::ostream &os, std::string section_delimiter,
                   const std::string &comments) const
{
   int64_t i, j;

   // serial/parallel conforming mesh format
   const bool set_names = attribute_sets.SetsExist() ||
     bdr_attribute_sets.SetsExist();
   os << (!set_names && section_delimiter.empty()
          ? "MFEM mesh v1.0\n" :
	  (!set_names ? "MFEM mesh v1.2\n" : "MFEM mesh v1.3\n"));

   if (set_names && section_delimiter.empty())
   {
      section_delimiter = "mfem_mesh_end";
   }

   // optional
   if (!comments.empty()) { os << '\n' << comments << '\n'; }

   os <<
      "\n#\n# MFEM Geometry Types (see fem/geom.hpp):\n#\n"
      "# POINT       = 0\n"
      "# SEGMENT     = 1\n"
      "# TRIANGLE    = 2\n"
      "# SQUARE      = 3\n"
      "# TETRAHEDRON = 4\n"
      "# CUBE        = 5\n"
      "# PRISM       = 6\n"
      "# PYRAMID     = 7\n"
      "#\n";

   os << "\ndimension\n" << Dim;

   os << "\n\nelements\n" << NumOfElements << '\n';
   for (i = 0; i < NumOfElements; i++)
   {
      PrintElement(elements[i], os);
   }

   if (set_names)
   {
     os << "\nattribute_sets\n";
     attribute_sets.Print(os);
   }

   os << "\nboundary\n" << NumOfBdrElements << '\n';
   for (i = 0; i < NumOfBdrElements; i++)
   {
      PrintElement(boundary[i], os);
   }

   if (set_names)
   {
     os << "\nbdr_attribute_sets\n";
     bdr_attribute_sets.Print(os);
   }

   os << "\nvertices\n" << NumOfVertices << '\n';
   os << spaceDim << '\n';
   for (i = 0; i < NumOfVertices; i++)
   {
      os << vertices[i](0);
      for (j = 1; j < spaceDim; j++)
      {
         os << ' ' << vertices[i](j);
      }
      os << '\n';
   }
   os.flush();

   if (!section_delimiter.empty())
   {
      os << '\n'
         << section_delimiter << endl; // only with formats v1.2 and above
   }
}

void Mesh::Save(const std::string &fname, int64_t precision) const
{
   ofstream ofs(fname);
   ofs.precision(precision);
   Print(ofs);
}

void Mesh::RemoveUnusedVertices()
{
   Array<int64_t> v2v(GetNV());
   v2v = -1;
   for (int64_t i = 0; i < GetNE(); i++)
   {
      Element *el = GetElement(i);
      int64_t nv = el->GetNVertices();
      int64_t *v = el->GetVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         v2v[v[j]] = 0;
      }
   }
   for (int64_t i = 0; i < GetNBE(); i++)
   {
      Element *el = GetBdrElement(i);
      int64_t *v = el->GetVertices();
      int64_t nv = el->GetNVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         v2v[v[j]] = 0;
      }
   }
   int64_t num_vert = 0;
   for (int64_t i = 0; i < v2v.Size(); i++)
   {
      if (v2v[i] == 0)
      {
         vertices[num_vert] = vertices[i];
         v2v[i] = num_vert++;
      }
   }

   if (num_vert == v2v.Size()) { return; }

   vertices.SetSize(num_vert);
   NumOfVertices = num_vert;
   for (int64_t i = 0; i < GetNE(); i++)
   {
      Element *el = GetElement(i);
      int64_t *v = el->GetVertices();
      int64_t nv = el->GetNVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         v[j] = v2v[v[j]];
      }
   }
   for (int64_t i = 0; i < GetNBE(); i++)
   {
      Element *el = GetBdrElement(i);
      int64_t *v = el->GetVertices();
      int64_t nv = el->GetNVertices();
      for (int64_t j = 0; j < nv; j++)
      {
         v[j] = v2v[v[j]];
      }
   }
   DeleteTables();
   if (Dim > 1)
   {
      // generate el_to_edge, be_to_face (2D), bel_to_edge (3D)
      el_to_edge = new Table;
      NumOfEdges = GetElementToEdgeTable(*el_to_edge);
   }
   if (Dim > 2)
   {
      // generate el_to_face, be_to_face
      GetElementToFaceTable();
   }
   // Update faces and faces_info
   GenerateFaces();
}

void Mesh::RemoveInternalBoundaries()
{
   int64_t num_bdr_elem = 0;
   int64_t new_bel_to_edge_nnz = 0;
   for (int64_t i = 0; i < GetNBE(); i++)
   {
      if (FaceIsInterior(GetBdrElementFaceIndex(i)))
      {
         FreeElement(boundary[i]);
      }
      else
      {
         num_bdr_elem++;
         if (Dim == 3)
         {
            new_bel_to_edge_nnz += bel_to_edge->RowSize(i);
         }
      }
   }

   if (num_bdr_elem == GetNBE()) { return; }

   Array<Element *> new_boundary(num_bdr_elem);
   Array<int64_t> new_be_to_face;
   Table *new_bel_to_edge = NULL;
   new_boundary.SetSize(0);
   new_be_to_face.Reserve(num_bdr_elem);
   if (Dim == 3)
   {
      new_bel_to_edge = new Table;
      new_bel_to_edge->SetDims(num_bdr_elem, new_bel_to_edge_nnz);
   }
   for (int64_t i = 0; i < GetNBE(); i++)
   {
      if (!FaceIsInterior(GetBdrElementFaceIndex(i)))
      {
         new_boundary.Append(boundary[i]);
         int64_t row = new_be_to_face.Size();
         new_be_to_face.Append(be_to_face[i]);
         if (Dim == 3)
         {
            int64_t *e = bel_to_edge->GetRow(i);
            int64_t ne = bel_to_edge->RowSize(i);
            int64_t *new_e = new_bel_to_edge->GetRow(row);
            for (int64_t j = 0; j < ne; j++)
            {
               new_e[j] = e[j];
            }
            new_bel_to_edge->GetI()[row+1] = new_bel_to_edge->GetI()[row] + ne;
         }
      }
   }

   NumOfBdrElements = new_boundary.Size();
   mfem::Swap(boundary, new_boundary);

   mfem::Swap(be_to_face, new_be_to_face);

   if (Dim == 3)
   {
      delete bel_to_edge;
      bel_to_edge = new_bel_to_edge;
   }

   Array<int64_t> attribs(num_bdr_elem);
   for (int64_t i = 0; i < attribs.Size(); i++)
   {
      attribs[i] = GetBdrAttribute(i);
   }
   attribs.Sort();
   attribs.Unique();
   bdr_attributes.DeleteAll();
   attribs.Copy(bdr_attributes);
}

void Mesh::FreeElement(Element *E)
{
   delete E;
}

std::ostream &operator<<(std::ostream &os, const Mesh &mesh)
{
   mesh.Print(os);
   return os;
}

MeshPart::EntityHelper::EntityHelper(
   int64_t dim_, const Array<int64_t> (&entity_to_vertex_)[Geometry::NumGeom])
   : dim(dim_),
     entity_to_vertex(entity_to_vertex_)
{
   int64_t geom_offset = 0;
   for (int64_t g = Geometry::DimStart[dim]; g < Geometry::DimStart[dim+1]; g++)
   {
      geom_offsets[g] = geom_offset;
      geom_offset += entity_to_vertex[g].Size()/Geometry::NumVerts[g];
   }
   geom_offsets[Geometry::DimStart[dim+1]] = geom_offset;
   num_entities = geom_offset;
}

MeshPart::Entity MeshPart::EntityHelper::FindEntity(int64_t bytype_entity_id)
{
   // Find the 'geom' that corresponds to 'bytype_entity_id'
   int64_t geom = Geometry::DimStart[dim];
   while (geom_offsets[geom+1] <= bytype_entity_id) { geom++; }
   MFEM_ASSERT(geom < Geometry::NumGeom, "internal error");
   MFEM_ASSERT(Geometry::Dimension[geom] == dim, "internal error");
   const int64_t nv = Geometry::NumVerts[geom];
   const int64_t geom_elem_id = bytype_entity_id - geom_offsets[geom];
   const int64_t *v = &entity_to_vertex[geom][nv*geom_elem_id];
   return { geom, nv, v };
}

void MeshPart::Print(std::ostream &os) const
{
   os << "MFEM mesh v1.2\n";

   // optional
   os <<
      "\n#\n# MFEM Geometry Types (see mesh/geom.hpp):\n#\n"
      "# POINT       = 0\n"
      "# SEGMENT     = 1\n"
      "# TRIANGLE    = 2\n"
      "# SQUARE      = 3\n"
      "# TETRAHEDRON = 4\n"
      "# CUBE        = 5\n"
      "# PRISM       = 6\n"
      "# PYRAMID     = 7\n"
      "#\n";

   const int64_t dim = dimension;
   os << "\ndimension\n" << dim;

   os << "\n\nelements\n" << num_elements << '\n';
   {
      const bool have_element_map = (element_map.Size() == num_elements);
      MFEM_ASSERT(have_element_map || element_map.Size() == 0,
                  "invalid MeshPart state");
      EntityHelper elem_helper(dim, entity_to_vertex);
      MFEM_ASSERT(elem_helper.num_entities == num_elements,
                  "invalid MeshPart state");
      for (int64_t nat_elem_id = 0; nat_elem_id < num_elements; nat_elem_id++)
      {
         const int64_t bytype_elem_id = have_element_map ?
                                    element_map[nat_elem_id] : nat_elem_id;
         const Entity ent = elem_helper.FindEntity(bytype_elem_id);
         // Print the element
         os << attributes[nat_elem_id] << ' ' << ent.geom;
         for (int64_t i = 0; i < ent.num_verts; i++)
         {
            os << ' ' << ent.verts[i];
         }
         os << '\n';
      }
   }

   os << "\nboundary\n" << num_bdr_elements << '\n';
   {
      const bool have_boundary_map = (boundary_map.Size() == num_bdr_elements);
      MFEM_ASSERT(have_boundary_map || boundary_map.Size() == 0,
                  "invalid MeshPart state");
      EntityHelper bdr_helper(dim-1, entity_to_vertex);
      MFEM_ASSERT(bdr_helper.num_entities == num_bdr_elements,
                  "invalid MeshPart state");
      for (int64_t nat_bdr_id = 0; nat_bdr_id < num_bdr_elements; nat_bdr_id++)
      {
         const int64_t bytype_bdr_id = have_boundary_map ?
                                   boundary_map[nat_bdr_id] : nat_bdr_id;
         const Entity ent = bdr_helper.FindEntity(bytype_bdr_id);
         // Print the boundary element
         os << bdr_attributes[nat_bdr_id] << ' ' << ent.geom;
         for (int64_t i = 0; i < ent.num_verts; i++)
         {
            os << ' ' << ent.verts[i];
         }
         os << '\n';
      }
   }

   os << "\nvertices\n" << num_vertices << '\n';
   const int64_t sdim = space_dimension;
   os << sdim << '\n';
   for (int64_t i = 0; i < num_vertices; i++)
   {
      os << vertex_coordinates[i*sdim];
      for (int64_t d = 1; d < sdim; d++)
      {
         os << ' ' << vertex_coordinates[i*sdim+d];
      }
      os << '\n';
   }

   os << "\nmfem_serial_mesh_end\n";

   // Start: GroupTopology::Save
   const int64_t num_groups = my_groups.Size();
   os << "\ncommunication_groups\n";
   os << "number_of_groups " << num_groups << "\n\n";

   os << "# number of entities in each group, followed by ranks in group\n";
   for (int64_t group_id = 0; group_id < num_groups; ++group_id)
   {
      const int64_t group_size = my_groups.RowSize(group_id);
      const int64_t *group_ptr = my_groups.GetRow(group_id);
      os << group_size;
      for (int64_t group_member_index = 0; group_member_index < group_size;
           ++group_member_index)
      {
         os << ' ' << group_ptr[group_member_index];
      }
      os << '\n';
   }
   // End: GroupTopology::Save

   const Table &g2v  = group_shared_entity_to_vertex[Geometry::POINT];
   const Table &g2ev = group_shared_entity_to_vertex[Geometry::SEGMENT];
   const Table &g2tv = group_shared_entity_to_vertex[Geometry::TRIANGLE];
   const Table &g2qv = group_shared_entity_to_vertex[Geometry::SQUARE];

   MFEM_VERIFY(g2v.RowSize(0) == 0, "internal erroor");
   os << "\ntotal_shared_vertices " << g2v.Size_of_connections() << '\n';
   if (dimension >= 2)
   {
      MFEM_VERIFY(g2ev.RowSize(0) == 0, "internal erroor");
      os << "total_shared_edges " << g2ev.Size_of_connections()/2 << '\n';
   }
   if (dimension >= 3)
   {
      MFEM_VERIFY(g2tv.RowSize(0) == 0, "internal erroor");
      MFEM_VERIFY(g2qv.RowSize(0) == 0, "internal erroor");
      const int64_t total_shared_faces =
         g2tv.Size_of_connections()/3 + g2qv.Size_of_connections()/4;
      os << "total_shared_faces " << total_shared_faces << '\n';
   }
   os << "\n# group 0 has no shared entities\n";
   for (int64_t gr = 1; gr < num_groups; gr++)
   {
      {
         const int64_t  nv = g2v.RowSize(gr);
         const int64_t *sv = g2v.GetRow(gr);
         os << "\n# group " << gr << "\nshared_vertices " << nv << '\n';
         for (int64_t i = 0; i < nv; i++)
         {
            os << sv[i] << '\n';
         }
      }
      if (dimension >= 2)
      {
         const int64_t  ne = g2ev.RowSize(gr)/2;
         const int64_t *se = g2ev.GetRow(gr);
         os << "\nshared_edges " << ne << '\n';
         for (int64_t i = 0; i < ne; i++)
         {
            const int64_t *v = se + 2*i;
            os << v[0] << ' ' << v[1] << '\n';
         }
      }
      if (dimension >= 3)
      {
         const int64_t  nt = g2tv.RowSize(gr)/3;
         const int64_t *st = g2tv.GetRow(gr);
         const int64_t  nq = g2qv.RowSize(gr)/4;
         const int64_t *sq = g2qv.GetRow(gr);
         os << "\nshared_faces " << nt+nq << '\n';
         for (int64_t i = 0; i < nt; i++)
         {
            os << Geometry::TRIANGLE;
            const int64_t *v = st + 3*i;
            for (int64_t j = 0; j < 3; j++) { os << ' ' << v[j]; }
            os << '\n';
         }
         for (int64_t i = 0; i < nq; i++)
         {
            os << Geometry::SQUARE;
            const int64_t *v = sq + 4*i;
            for (int64_t j = 0; j < 4; j++) { os << ' ' << v[j]; }
            os << '\n';
         }
      }
   }

   // Write out section end tag for mesh.
   os << "\nmfem_mesh_end" << endl;
}

Mesh &MeshPart::GetMesh()
{
   if (mesh) { return *mesh; }

   mesh.reset(new Mesh(dimension,
                       num_vertices,
                       num_elements,
                       num_bdr_elements,
                       space_dimension));

   // Add elements
   {
      const bool have_element_map = (element_map.Size() == num_elements);
      MFEM_ASSERT(have_element_map || element_map.Size() == 0,
                  "invalid MeshPart state");
      EntityHelper elem_helper(dimension, entity_to_vertex);
      MFEM_ASSERT(elem_helper.num_entities == num_elements,
                  "invalid MeshPart state");
      const bool have_tet_refine_flags = (tet_refine_flags.Size() > 0);
      for (int64_t nat_elem_id = 0; nat_elem_id < num_elements; nat_elem_id++)
      {
         const int64_t bytype_elem_id = have_element_map ?
                                    element_map[nat_elem_id] : nat_elem_id;
         const Entity ent = elem_helper.FindEntity(bytype_elem_id);
         Element *el = mesh->NewElement(ent.geom);
         el->SetVertices(ent.verts);
         el->SetAttribute(attributes[nat_elem_id]);
         if (ent.geom == Geometry::TETRAHEDRON && have_tet_refine_flags)
         {
            constexpr int64_t geom_tet = Geometry::TETRAHEDRON;
            const int64_t tet_id = (ent.verts - entity_to_vertex[geom_tet])/4;
            const int64_t ref_flag = tet_refine_flags[tet_id];
            static_cast<Tetrahedron*>(el)->SetRefinementFlag(ref_flag);
         }
         mesh->AddElement(el);
      }
   }

   // Add boundary elements
   {
      const bool have_boundary_map = (boundary_map.Size() == num_bdr_elements);
      MFEM_ASSERT(have_boundary_map || boundary_map.Size() == 0,
                  "invalid MeshPart state");
      EntityHelper bdr_helper(dimension-1, entity_to_vertex);
      MFEM_ASSERT(bdr_helper.num_entities == num_bdr_elements,
                  "invalid MeshPart state");
      for (int64_t nat_bdr_id = 0; nat_bdr_id < num_bdr_elements; nat_bdr_id++)
      {
         const int64_t bytype_bdr_id = have_boundary_map ?
                                   boundary_map[nat_bdr_id] : nat_bdr_id;
         const Entity ent = bdr_helper.FindEntity(bytype_bdr_id);
         Element *bdr = mesh->NewElement(ent.geom);
         bdr->SetVertices(ent.verts);
         bdr->SetAttribute(bdr_attributes[nat_bdr_id]);
         mesh->AddBdrElement(bdr);
      }
   }

   // Add vertices
   if (vertex_coordinates.Size() == space_dimension*num_vertices)
   {
      MFEM_ASSERT(!nodes, "invalid MeshPart state");
      for (int64_t vert_id = 0; vert_id < num_vertices; vert_id++)
      {
         mesh->AddVertex(vertex_coordinates + space_dimension*vert_id);
      }
   }
   else
   {
      MFEM_ASSERT(vertex_coordinates.Size() == 0, "invalid MeshPart state");
      for (int64_t vert_id = 0; vert_id < num_vertices; vert_id++)
      {
         mesh->AddVertex(0., 0., 0.);
      }
      // 'mesh.Nodes' cannot be set here -- they can be set later, if needed
   }

   mesh->FinalizeTopology(/* generate_bdr: */ false);

   return *mesh;
}


MeshPartitioner::MeshPartitioner(Mesh &mesh_,
                                 int64_t num_parts_,
                                 const int64_t *partitioning_,
                                 int64_t part_method)
   : mesh(mesh_)
{
   if (partitioning_)
   {
      partitioning.MakeRef(const_cast<int64_t *>(partitioning_), mesh.GetNE(),
                           false);
   }
   else
   {
      // Mesh::GeneratePartitioning always uses new[] to allocate the,
      // partitioning, so we need to tell the memory manager to free it with
      // delete[] (even if a different host memory type has been selected).
      constexpr MemoryType mt = MemoryType::HOST;
      partitioning.MakeRef(mesh.GeneratePartitioning(num_parts_, part_method),
                           mesh.GetNE(), mt, true);
   }

   Transpose(partitioning, part_to_element, num_parts_);
   // Note: the element ids in each row of 'part_to_element' are sorted.

   const int64_t dim = mesh.Dimension();
   if (dim >= 2)
   {
      Transpose(mesh.ElementToEdgeTable(), edge_to_element, mesh.GetNEdges());
   }

   Array<int64_t> boundary_to_part(mesh.GetNBE());
   // Same logic as in ParMesh::BuildLocalBoundary
   if (dim >= 3)
   {
      for (int64_t i = 0; i < boundary_to_part.Size(); i++)
      {
         int64_t face, o, el1, el2;
         mesh.GetBdrElementFace(i, &face, &o);
         mesh.GetFaceElements(face, &el1, &el2);
         boundary_to_part[i] =
            partitioning[(o % 2 == 0 || el2 < 0) ? el1 : el2];
      }
   }
   else if (dim == 2)
   {
      for (int64_t i = 0; i < boundary_to_part.Size(); i++)
      {
         int64_t edge = mesh.GetBdrElementFaceIndex(i);
         int64_t el1 = edge_to_element.GetRow(edge)[0];
         boundary_to_part[i] = partitioning[el1];
      }
   }
   else if (dim == 1)
   {
      for (int64_t i = 0; i < boundary_to_part.Size(); i++)
      {
         int64_t vert = mesh.GetBdrElementFaceIndex(i);
         int64_t el1, el2;
         mesh.GetFaceElements(vert, &el1, &el2);
         boundary_to_part[i] = partitioning[el1];
      }
   }
   Transpose(boundary_to_part, part_to_boundary, num_parts_);
   // Note: the boundary element ids in each row of 'part_to_boundary' are
   // sorted.
   boundary_to_part.DeleteAll();

   Table *vert_element = mesh.GetVertexToElementTable(); // we must delete this
   vertex_to_element.Swap(*vert_element);
   delete vert_element;
}

void MeshPartitioner::ExtractPart(int64_t part_id, MeshPart &mesh_part) const
{
   const int64_t num_parts = part_to_element.Size();

   MFEM_VERIFY(0 <= part_id && part_id < num_parts,
               "invalid part_id = " << part_id
               << ", num_parts = " << num_parts);

   const int64_t dim = mesh.Dimension();
   const int64_t sdim = mesh.SpaceDimension();
   const int64_t num_elems = part_to_element.RowSize(part_id);
   const int64_t *elem_list = part_to_element.GetRow(part_id); // sorted
   const int64_t num_bdr_elems = part_to_boundary.RowSize(part_id);
   const int64_t *bdr_elem_list = part_to_boundary.GetRow(part_id); // sorted

   // Initialize 'mesh_part'
   mesh_part.dimension = dim;
   mesh_part.space_dimension = sdim;
   mesh_part.num_vertices = 0;
   mesh_part.num_elements = num_elems;
   mesh_part.num_bdr_elements = num_bdr_elems;
   for (int64_t g = 0; g < Geometry::NumGeom; g++)
   {
      mesh_part.entity_to_vertex[g].SetSize(0); // can reuse Array allocation
   }
   mesh_part.tet_refine_flags.SetSize(0);
   mesh_part.element_map.SetSize(0); // 0 or 'num_elements', if needed
   mesh_part.boundary_map.SetSize(0); // 0 or 'num_bdr_elements', if needed
   mesh_part.attributes.SetSize(num_elems);
   mesh_part.bdr_attributes.SetSize(num_bdr_elems);
   mesh_part.vertex_coordinates.SetSize(0);

   mesh_part.num_parts = num_parts;
   mesh_part.my_part_id = part_id;
   mesh_part.my_groups.Clear();
   for (int64_t g = 0; g < Geometry::NumGeom; g++)
   {
      mesh_part.group_shared_entity_to_vertex[g].Clear();
   }
   mesh_part.mesh.reset(nullptr);

   // Initialize:
   // - 'mesh_part.entity_to_vertex' for the elements (boundary elements are
   //   set later); vertex ids are global at this point - they will be mapped to
   //   local ids later
   // - 'mesh_part.attributes'
   // - 'mesh_part.tet_refine_flags' if needed
   int64_t geom_marker = 0, num_geom = 0;
   for (int64_t i = 0; i < num_elems; i++)
   {
      const Element *elem = mesh.GetElement(elem_list[i]);
      const int64_t geom = elem->GetGeometryType();
      const int64_t nv = Geometry::NumVerts[geom];
      const int64_t *v = elem->GetVertices();
      MFEM_VERIFY(numeric_limits<int64_t>::max() - nv >=
                  mesh_part.entity_to_vertex[geom].Size(),
                  "overflow in 'entity_to_vertex[geom]', geom: "
                  << Geometry::Name[geom]);
      mesh_part.entity_to_vertex[geom].Append(v, nv);
      mesh_part.attributes[i] = elem->GetAttribute();
      if (geom == Geometry::TETRAHEDRON)
      {
         // Create 'mesh_part.tet_refine_flags' but only if we find at least one
         // non-zero flag in a tetrahedron.
         const Tetrahedron *tet = static_cast<const Tetrahedron*>(elem);
         const int64_t ref_flag = tet->GetRefinementFlag();
         if (mesh_part.tet_refine_flags.Size() == 0)
         {
            if (ref_flag)
            {
               // This is the first time we encounter non-zero 'ref_flag'
               const int64_t num_tets = mesh_part.entity_to_vertex[geom].Size()/nv;
               mesh_part.tet_refine_flags.SetSize(num_tets, 0);
               mesh_part.tet_refine_flags.Last() = ref_flag;
            }
         }
         else
         {
            mesh_part.tet_refine_flags.Append(ref_flag);
         }
      }
      if ((geom_marker & (1 << geom)) == 0)
      {
         geom_marker |= (1 << geom);
         num_geom++;
      }
   }
   MFEM_ASSERT(mesh_part.tet_refine_flags.Size() == 0 ||
               mesh_part.tet_refine_flags.Size() ==
               mesh_part.entity_to_vertex[Geometry::TETRAHEDRON].Size()/4,
               "internal error");
   // Initialize 'mesh_part.element_map' if needed
   if (num_geom > 1)
   {
      int64_t offsets[Geometry::NumGeom];
      int64_t offset = 0;
      for (int64_t g = Geometry::DimStart[dim]; g < Geometry::DimStart[dim+1]; g++)
      {
         offsets[g] = offset;
         offset += mesh_part.entity_to_vertex[g].Size()/Geometry::NumVerts[g];
      }
      mesh_part.element_map.SetSize(num_elems);
      for (int64_t i = 0; i < num_elems; i++)
      {
         const int64_t geom = mesh.GetElementGeometry(elem_list[i]);
         mesh_part.element_map[i] = offsets[geom]++;
      }
   }

   // Initialize:
   // - 'mesh_part.entity_to_vertex' for the boundary elements; vertex ids are
   //   global at this point - they will be mapped to local ids later
   // - 'mesh_part.bdr_attributes'
   geom_marker = 0; num_geom = 0;
   for (int64_t i = 0; i < num_bdr_elems; i++)
   {
      const Element *bdr_elem = mesh.GetBdrElement(bdr_elem_list[i]);
      const int64_t geom = bdr_elem->GetGeometryType();
      const int64_t nv = Geometry::NumVerts[geom];
      const int64_t *v = bdr_elem->GetVertices();
      MFEM_VERIFY(numeric_limits<int64_t>::max() - nv >=
                  mesh_part.entity_to_vertex[geom].Size(),
                  "overflow in 'entity_to_vertex[geom]', geom: "
                  << Geometry::Name[geom]);
      mesh_part.entity_to_vertex[geom].Append(v, nv);
      mesh_part.bdr_attributes[i] = bdr_elem->GetAttribute();
      if ((geom_marker & (1 << geom)) == 0)
      {
         geom_marker |= (1 << geom);
         num_geom++;
      }
   }
   // Initialize 'mesh_part.boundary_map' if needed
   if (num_geom > 1)
   {
      int64_t offsets[Geometry::NumGeom];
      int64_t offset = 0;
      for (int64_t g = Geometry::DimStart[dim-1]; g < Geometry::DimStart[dim]; g++)
      {
         offsets[g] = offset;
         offset += mesh_part.entity_to_vertex[g].Size()/Geometry::NumVerts[g];
      }
      mesh_part.boundary_map.SetSize(num_bdr_elems);
      for (int64_t i = 0; i < num_bdr_elems; i++)
      {
         const int64_t geom = mesh.GetBdrElementGeometry(bdr_elem_list[i]);
         mesh_part.boundary_map[i] = offsets[geom]++;
      }
   }

   // Create the vertex id map, 'vertex_loc_to_glob', which maps local ids to
   // global ones; the map is sorted, preserving the global ordering.
   Array<int64_t> vertex_loc_to_glob;
   {
      std::unordered_set<int64_t> vertex_set;
      for (int64_t i = 0; i < num_elems; i++)
      {
         const Element *elem = mesh.GetElement(elem_list[i]);
         const int64_t geom = elem->GetGeometryType();
         const int64_t nv = Geometry::NumVerts[geom];
         const int64_t *v = elem->GetVertices();
         vertex_set.insert(v, v + nv);
      }
      vertex_loc_to_glob.SetSize(static_cast<int64_t>(vertex_set.size()));
      std::copy(vertex_set.begin(), vertex_set.end(), // src
                vertex_loc_to_glob.begin());          // dest
   }
   vertex_loc_to_glob.Sort();

   // Initialize 'mesh_part.num_vertices'
   mesh_part.num_vertices = vertex_loc_to_glob.Size();

   // Update the vertex ids in the arrays 'mesh_part.entity_to_vertex' from
   // global to local.
   for (int64_t g = 0; g < Geometry::NumGeom; g++)
   {
      Array<int64_t> &vert_array = mesh_part.entity_to_vertex[g];
      for (int64_t i = 0; i < vert_array.Size(); i++)
      {
         const int64_t glob_id = vert_array[i];
         const int64_t loc_id = vertex_loc_to_glob.FindSorted(glob_id);
         MFEM_ASSERT(loc_id >= 0, "internal error: global vertex id not found");
         vert_array[i] = loc_id;
      }
   }

   // Initialize one of 'mesh_part.vertex_coordinates' or 'mesh_part.nodes'
   MFEM_VERIFY(numeric_limits<int64_t>::max()/sdim >= vertex_loc_to_glob.Size(),
               "overflow in 'vertex_coordinates', num_vertices = "
               << vertex_loc_to_glob.Size() << ", sdim = " << sdim);
   mesh_part.vertex_coordinates.SetSize(sdim*vertex_loc_to_glob.Size());
   for (int64_t i = 0; i < vertex_loc_to_glob.Size(); i++)
   {
      const real_t *coord = mesh.GetVertex(vertex_loc_to_glob[i]);
      for (int64_t d = 0; d < sdim; d++)
      {
         mesh_part.vertex_coordinates[i*sdim+d] = coord[d];
      }
   }

   // Begin constructing the "neighbor" groups, i.e. the groups that contain
   // 'part_id'.
   ListOfIntegerSets groups;
   {
      // the first group is the local one
      IntegerSet group;
      group.Recreate(1, &part_id);
      groups.Insert(group);
   }

   // 'shared_faces' : shared face id -> (global_face_id, group_id)
   // Note: 'shared_faces' will be sorted by 'global_face_id'.
   Array<Pair<int64_t,int64_t>> shared_faces;

   // Add "neighbor" groups defined by faces
   // Construct 'shared_faces'.
   if (dim >= 3)
   {
      std::unordered_set<int64_t> face_set;
      // Construct 'face_set'
      const Table &elem_to_face = mesh.ElementToFaceTable();
      for (int64_t loc_elem_id = 0; loc_elem_id < num_elems; loc_elem_id++)
      {
         const int64_t glob_elem_id = elem_list[loc_elem_id];
         const int64_t nfaces = elem_to_face.RowSize(glob_elem_id);
         const int64_t *faces = elem_to_face.GetRow(glob_elem_id);
         face_set.insert(faces, faces + nfaces);
      }
      // Construct 'shared_faces'; add "neighbor" groups defined by faces.
      IntegerSet group;
      for (int64_t glob_face_id : face_set)
      {
         int64_t el[2];
         mesh.GetFaceElements(glob_face_id, &el[0], &el[1]);
         if (el[1] < 0) { continue; }
         el[0] = partitioning[el[0]];
         el[1] = partitioning[el[1]];
         MFEM_ASSERT(el[0] == part_id || el[1] == part_id, "internal error");
         if (el[0] != part_id || el[1] != part_id)
         {
            group.Recreate(2, el);
            const int64_t group_id = groups.Insert(group);
            shared_faces.Append(Pair<int64_t,int64_t>(glob_face_id, group_id));
         }
      }
      shared_faces.Sort(); // sort the shared faces by 'glob_face_id'
   }

   // 'shared_edges' : shared edge id -> (global_edge_id, group_id)
   // Note: 'shared_edges' will be sorted by 'global_edge_id'.
   Array<Pair<int64_t,int64_t>> shared_edges;

   // Add "neighbor" groups defined by edges.
   // Construct 'shared_edges'.
   if (dim >= 2)
   {
      std::unordered_set<int64_t> edge_set;
      // Construct 'edge_set'
      const Table &elem_to_edge = mesh.ElementToEdgeTable();
      for (int64_t loc_elem_id = 0; loc_elem_id < num_elems; loc_elem_id++)
      {
         const int64_t glob_elem_id = elem_list[loc_elem_id];
         const int64_t nedges = elem_to_edge.RowSize(glob_elem_id);
         const int64_t *edges = elem_to_edge.GetRow(glob_elem_id);
         edge_set.insert(edges, edges + nedges);
      }
      // Construct 'shared_edges'; add "neighbor" groups defined by edges.
      IntegerSet group;
      for (int64_t glob_edge_id : edge_set)
      {
         const int64_t nelem = edge_to_element.RowSize(glob_edge_id);
         const int64_t *elem = edge_to_element.GetRow(glob_edge_id);
         Array<int64_t> &gr = group; // reference to the 'group' internal Array
         gr.SetSize(nelem);
         for (int64_t j = 0; j < nelem; j++)
         {
            gr[j] = partitioning[elem[j]];
         }
         gr.Sort();
         gr.Unique();
         MFEM_ASSERT(gr.FindSorted(part_id) >= 0, "internal error");
         if (group.Size() > 1)
         {
            const int64_t group_id = groups.Insert(group);
            shared_edges.Append(Pair<int64_t,int64_t>(glob_edge_id, group_id));
         }
      }
      shared_edges.Sort(); // sort the shared edges by 'glob_edge_id'
   }

   // 'shared_verts' : shared vertex id -> (global_vertex_id, group_id)
   // Note: 'shared_verts' will be sorted by 'global_vertex_id'.
   Array<Pair<int64_t,int64_t>> shared_verts;

   // Add "neighbor" groups defined by vertices.
   // Construct 'shared_verts'.
   {
      IntegerSet group;
      for (int64_t i = 0; i < vertex_loc_to_glob.Size(); i++)
      {
         // 'vertex_to_element' maps global vertex ids to global element ids
         const int64_t glob_vertex_id = vertex_loc_to_glob[i];
         const int64_t nelem = vertex_to_element.RowSize(glob_vertex_id);
         const int64_t *elem = vertex_to_element.GetRow(glob_vertex_id);
         Array<int64_t> &gr = group; // reference to the 'group' internal Array
         gr.SetSize(nelem);
         for (int64_t j = 0; j < nelem; j++)
         {
            gr[j] = partitioning[elem[j]];
         }
         gr.Sort();
         gr.Unique();
         MFEM_ASSERT(gr.FindSorted(part_id) >= 0, "internal error");
         if (group.Size() > 1)
         {
            const int64_t group_id = groups.Insert(group);
            shared_verts.Append(Pair<int64_t,int64_t>(glob_vertex_id, group_id));
         }
      }
   }

   // Done constructing the "neighbor" groups in 'groups'.
   const int64_t num_groups = groups.Size();

   // Define 'mesh_part.my_groups'
   groups.AsTable(mesh_part.my_groups);

   // Construct 'mesh_part.group_shared_entity_to_vertex[Geometry::POINT]'
   Table &group__shared_vertex_to_vertex =
      mesh_part.group_shared_entity_to_vertex[Geometry::POINT];
   group__shared_vertex_to_vertex.MakeI(num_groups);
   for (int64_t sv = 0; sv < shared_verts.Size(); sv++)
   {
      const int64_t group_id = shared_verts[sv].two;
      group__shared_vertex_to_vertex.AddAColumnInRow(group_id);
   }
   group__shared_vertex_to_vertex.MakeJ();
   for (int64_t sv = 0; sv < shared_verts.Size(); sv++)
   {
      const int64_t glob_vertex_id = shared_verts[sv].one;
      const int64_t group_id       = shared_verts[sv].two;
      const int64_t loc_vertex_id = vertex_loc_to_glob.FindSorted(glob_vertex_id);
      MFEM_ASSERT(loc_vertex_id >= 0, "internal error");
      group__shared_vertex_to_vertex.AddConnection(group_id, loc_vertex_id);
   }
   group__shared_vertex_to_vertex.ShiftUpI();

   // Construct 'mesh_part.group_shared_entity_to_vertex[Geometry::SEGMENT]'
   if (dim >= 2)
   {
      Table &group__shared_edge_to_vertex =
         mesh_part.group_shared_entity_to_vertex[Geometry::SEGMENT];
      group__shared_edge_to_vertex.MakeI(num_groups);
      for (int64_t se = 0; se < shared_edges.Size(); se++)
      {
         const int64_t group_id = shared_edges[se].two;
         group__shared_edge_to_vertex.AddColumnsInRow(group_id, 2);
      }
      group__shared_edge_to_vertex.MakeJ();
      const Table &edge_to_vertex = *mesh.GetEdgeVertexTable();
      for (int64_t se = 0; se < shared_edges.Size(); se++)
      {
         const int64_t glob_edge_id = shared_edges[se].one;
         const int64_t group_id     = shared_edges[se].two;
         const int64_t *v = edge_to_vertex.GetRow(glob_edge_id);
         for (int64_t i = 0; i < 2; i++)
         {
            const int64_t loc_vertex_id = vertex_loc_to_glob.FindSorted(v[i]);
            MFEM_ASSERT(loc_vertex_id >= 0, "internal error");
            group__shared_edge_to_vertex.AddConnection(group_id, loc_vertex_id);
         }
      }
      group__shared_edge_to_vertex.ShiftUpI();
   }

   // Construct 'mesh_part.group_shared_entity_to_vertex[Geometry::TRIANGLE]'
   // and 'mesh_part.group_shared_entity_to_vertex[Geometry::SQUARE]'.
   if (dim >= 3)
   {
      Table &group__shared_tria_to_vertex =
         mesh_part.group_shared_entity_to_vertex[Geometry::TRIANGLE];
      Table &group__shared_quad_to_vertex =
         mesh_part.group_shared_entity_to_vertex[Geometry::SQUARE];
      Array<int64_t> vertex_ids;
      group__shared_tria_to_vertex.MakeI(num_groups);
      group__shared_quad_to_vertex.MakeI(num_groups);
      for (int64_t sf = 0; sf < shared_faces.Size(); sf++)
      {
         const int64_t glob_face_id = shared_faces[sf].one;
         const int64_t group_id     = shared_faces[sf].two;
         const int64_t geom         = mesh.GetFaceGeometry(glob_face_id);
         mesh_part.group_shared_entity_to_vertex[geom].
         AddColumnsInRow(group_id, Geometry::NumVerts[geom]);
      }
      group__shared_tria_to_vertex.MakeJ();
      group__shared_quad_to_vertex.MakeJ();
      for (int64_t sf = 0; sf < shared_faces.Size(); sf++)
      {
         const int64_t glob_face_id = shared_faces[sf].one;
         const int64_t group_id     = shared_faces[sf].two;
         const int64_t geom         = mesh.GetFaceGeometry(glob_face_id);
         mesh.GetFaceVertices(glob_face_id, vertex_ids);
         // See also ParMesh::BuildSharedFaceElems.
         if (geom == Geometry::TRIANGLE)
         {
            int64_t glob_el_id[2];
            mesh.GetFaceElements(glob_face_id, &glob_el_id[0], &glob_el_id[1]);
            int64_t side = 0;
            const Element *el = mesh.GetElement(glob_el_id[0]);
            const Tetrahedron *tet = nullptr;
            if (el->GetGeometryType() == Geometry::TETRAHEDRON)
            {
               tet = static_cast<const Tetrahedron*>(el);
            }
            else
            {
               side = 1;
               el = mesh.GetElement(glob_el_id[1]);
               if (el->GetGeometryType() == Geometry::TETRAHEDRON)
               {
                  tet = static_cast<const Tetrahedron*>(el);
               }
            }
            if (tet && tet->GetRefinementFlag())
            {
               // mark the shared face for refinement by reorienting
               // it according to the refinement flag in the tetrahedron
               // to which this shared face belongs to.
               int64_t info[2];
               mesh.GetFaceInfos(glob_face_id, &info[0], &info[1]);
               tet->GetMarkedFace(info[side]/64, &vertex_ids[0]);
            }
         }
         for (int64_t i = 0; i < vertex_ids.Size(); i++)
         {
            const int64_t glob_id = vertex_ids[i];
            const int64_t loc_id = vertex_loc_to_glob.FindSorted(glob_id);
            MFEM_ASSERT(loc_id >= 0, "internal error");
            vertex_ids[i] = loc_id;
         }
         mesh_part.group_shared_entity_to_vertex[geom].
         AddConnections(group_id, vertex_ids, vertex_ids.Size());
      }
      group__shared_tria_to_vertex.ShiftUpI();
      group__shared_quad_to_vertex.ShiftUpI();
   }
}

}
