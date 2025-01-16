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

#include "geom.hpp"
#include "wedge.hpp"
#include "pyramid.hpp"

namespace mfem
{

const char *Geometry::Name[NumGeom] =
{
   "Point", "Segment", "Triangle", "Square", "Tetrahedron", "Cube", "Prism",
   "Pyramid"
};

const real_t Geometry::Volume[NumGeom] =
{ 1.0, 1.0, 0.5, 1.0, 1./6, 1.0, 0.5, 1./3 };

Geometry::Geometry()
{
   // Vertices for Geometry::POINT
   GeomVert[0] =  new IntegrationRule(1);
   GeomVert[0]->IntPoint(0).x = 0.0;

   // Vertices for Geometry::SEGMENT
   GeomVert[1] = new IntegrationRule(2);

   GeomVert[1]->IntPoint(0).x = 0.0;
   GeomVert[1]->IntPoint(1).x = 1.0;

   // Vertices for Geometry::TRIANGLE
   GeomVert[2] = new IntegrationRule(3);

   GeomVert[2]->IntPoint(0).x = 0.0;
   GeomVert[2]->IntPoint(0).y = 0.0;

   GeomVert[2]->IntPoint(1).x = 1.0;
   GeomVert[2]->IntPoint(1).y = 0.0;

   GeomVert[2]->IntPoint(2).x = 0.0;
   GeomVert[2]->IntPoint(2).y = 1.0;

   // Vertices for Geometry::SQUARE
   GeomVert[3] = new IntegrationRule(4);

   GeomVert[3]->IntPoint(0).x = 0.0;
   GeomVert[3]->IntPoint(0).y = 0.0;

   GeomVert[3]->IntPoint(1).x = 1.0;
   GeomVert[3]->IntPoint(1).y = 0.0;

   GeomVert[3]->IntPoint(2).x = 1.0;
   GeomVert[3]->IntPoint(2).y = 1.0;

   GeomVert[3]->IntPoint(3).x = 0.0;
   GeomVert[3]->IntPoint(3).y = 1.0;

   // Vertices for Geometry::TETRAHEDRON
   GeomVert[4] = new IntegrationRule(4);
   GeomVert[4]->IntPoint(0).x = 0.0;
   GeomVert[4]->IntPoint(0).y = 0.0;
   GeomVert[4]->IntPoint(0).z = 0.0;

   GeomVert[4]->IntPoint(1).x = 1.0;
   GeomVert[4]->IntPoint(1).y = 0.0;
   GeomVert[4]->IntPoint(1).z = 0.0;

   GeomVert[4]->IntPoint(2).x = 0.0;
   GeomVert[4]->IntPoint(2).y = 1.0;
   GeomVert[4]->IntPoint(2).z = 0.0;

   GeomVert[4]->IntPoint(3).x = 0.0;
   GeomVert[4]->IntPoint(3).y = 0.0;
   GeomVert[4]->IntPoint(3).z = 1.0;

   // Vertices for Geometry::CUBE
   GeomVert[5] = new IntegrationRule(8);

   GeomVert[5]->IntPoint(0).x = 0.0;
   GeomVert[5]->IntPoint(0).y = 0.0;
   GeomVert[5]->IntPoint(0).z = 0.0;

   GeomVert[5]->IntPoint(1).x = 1.0;
   GeomVert[5]->IntPoint(1).y = 0.0;
   GeomVert[5]->IntPoint(1).z = 0.0;

   GeomVert[5]->IntPoint(2).x = 1.0;
   GeomVert[5]->IntPoint(2).y = 1.0;
   GeomVert[5]->IntPoint(2).z = 0.0;

   GeomVert[5]->IntPoint(3).x = 0.0;
   GeomVert[5]->IntPoint(3).y = 1.0;
   GeomVert[5]->IntPoint(3).z = 0.0;

   GeomVert[5]->IntPoint(4).x = 0.0;
   GeomVert[5]->IntPoint(4).y = 0.0;
   GeomVert[5]->IntPoint(4).z = 1.0;

   GeomVert[5]->IntPoint(5).x = 1.0;
   GeomVert[5]->IntPoint(5).y = 0.0;
   GeomVert[5]->IntPoint(5).z = 1.0;

   GeomVert[5]->IntPoint(6).x = 1.0;
   GeomVert[5]->IntPoint(6).y = 1.0;
   GeomVert[5]->IntPoint(6).z = 1.0;

   GeomVert[5]->IntPoint(7).x = 0.0;
   GeomVert[5]->IntPoint(7).y = 1.0;
   GeomVert[5]->IntPoint(7).z = 1.0;

   // Vertices for Geometry::PRISM
   GeomVert[6] = new IntegrationRule(6);
   GeomVert[6]->IntPoint(0).x = 0.0;
   GeomVert[6]->IntPoint(0).y = 0.0;
   GeomVert[6]->IntPoint(0).z = 0.0;

   GeomVert[6]->IntPoint(1).x = 1.0;
   GeomVert[6]->IntPoint(1).y = 0.0;
   GeomVert[6]->IntPoint(1).z = 0.0;

   GeomVert[6]->IntPoint(2).x = 0.0;
   GeomVert[6]->IntPoint(2).y = 1.0;
   GeomVert[6]->IntPoint(2).z = 0.0;

   GeomVert[6]->IntPoint(3).x = 0.0;
   GeomVert[6]->IntPoint(3).y = 0.0;
   GeomVert[6]->IntPoint(3).z = 1.0;

   GeomVert[6]->IntPoint(4).x = 1.0;
   GeomVert[6]->IntPoint(4).y = 0.0;
   GeomVert[6]->IntPoint(4).z = 1.0;

   GeomVert[6]->IntPoint(5).x = 0.0;
   GeomVert[6]->IntPoint(5).y = 1.0;
   GeomVert[6]->IntPoint(5).z = 1.0;

   // Vertices for Geometry::PYRAMID
   GeomVert[7] = new IntegrationRule(5);
   GeomVert[7]->IntPoint(0).x = 0.0;
   GeomVert[7]->IntPoint(0).y = 0.0;
   GeomVert[7]->IntPoint(0).z = 0.0;

   GeomVert[7]->IntPoint(1).x = 1.0;
   GeomVert[7]->IntPoint(1).y = 0.0;
   GeomVert[7]->IntPoint(1).z = 0.0;

   GeomVert[7]->IntPoint(2).x = 1.0;
   GeomVert[7]->IntPoint(2).y = 1.0;
   GeomVert[7]->IntPoint(2).z = 0.0;

   GeomVert[7]->IntPoint(3).x = 0.0;
   GeomVert[7]->IntPoint(3).y = 1.0;
   GeomVert[7]->IntPoint(3).z = 0.0;

   GeomVert[7]->IntPoint(4).x = 0.0;
   GeomVert[7]->IntPoint(4).y = 0.0;
   GeomVert[7]->IntPoint(4).z = 1.0;

   GeomCenter[POINT].x = 0.0;
   GeomCenter[POINT].y = 0.0;
   GeomCenter[POINT].z = 0.0;

   GeomCenter[SEGMENT].x = 0.5;
   GeomCenter[SEGMENT].y = 0.0;
   GeomCenter[SEGMENT].z = 0.0;

   GeomCenter[TRIANGLE].x = 1.0 / 3.0;
   GeomCenter[TRIANGLE].y = 1.0 / 3.0;
   GeomCenter[TRIANGLE].z = 0.0;

   GeomCenter[SQUARE].x = 0.5;
   GeomCenter[SQUARE].y = 0.5;
   GeomCenter[SQUARE].z = 0.0;

   GeomCenter[TETRAHEDRON].x = 0.25;
   GeomCenter[TETRAHEDRON].y = 0.25;
   GeomCenter[TETRAHEDRON].z = 0.25;

   GeomCenter[CUBE].x = 0.5;
   GeomCenter[CUBE].y = 0.5;
   GeomCenter[CUBE].z = 0.5;

   GeomCenter[PRISM].x = 1.0 / 3.0;
   GeomCenter[PRISM].y = 1.0 / 3.0;
   GeomCenter[PRISM].z = 0.5;

   GeomCenter[PYRAMID].x = 0.375;
   GeomCenter[PYRAMID].y = 0.375;
   GeomCenter[PYRAMID].z = 0.25;
}

template <Geometry::Type GEOM>
int64_t GetInverseOrientation_(int64_t orientation)
{
   using geom_t = Geometry::Constants<GEOM>;
   MFEM_ASSERT(0 <= orientation && orientation < geom_t::NumOrient,
               "Invalid orientation");
   return geom_t::InvOrient[orientation];
}

int64_t Geometry::GetInverseOrientation(Type geom_type, int64_t orientation)
{
   switch (geom_type)
   {
      case Geometry::POINT:
         return GetInverseOrientation_<Geometry::POINT>(orientation);
      case Geometry::SEGMENT:
         return GetInverseOrientation_<Geometry::SEGMENT>(orientation);
      case Geometry::TRIANGLE:
         return GetInverseOrientation_<Geometry::TRIANGLE>(orientation);
      case Geometry::SQUARE:
         return GetInverseOrientation_<Geometry::SQUARE>(orientation);
      case Geometry::TETRAHEDRON:
         return GetInverseOrientation_<Geometry::TETRAHEDRON>(orientation);
      default:
         MFEM_ABORT("Geometry type does not have inverse orientations");
   }
}

Geometry::~Geometry()
{
   for (int64_t i = 0; i < NumGeom; i++)
   {
      delete PerfGeomToGeomJac[i];
      delete GeomToPerfGeomJac[i];
      delete GeomVert[i];
   }
}

const IntegrationRule *Geometry::GetVertices(int64_t GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return GeomVert[0];
      case Geometry::SEGMENT:     return GeomVert[1];
      case Geometry::TRIANGLE:    return GeomVert[2];
      case Geometry::SQUARE:      return GeomVert[3];
      case Geometry::TETRAHEDRON: return GeomVert[4];
      case Geometry::CUBE:        return GeomVert[5];
      case Geometry::PRISM:       return GeomVert[6];
      case Geometry::PYRAMID:     return GeomVert[7];
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         mfem_error("Geometry::GetVertices(...)");
   }
   // make some compilers happy.
   return GeomVert[0];
}

// static method
void Geometry::GetRandomPoint(int64_t GeomType, IntegrationPoint &ip)
{
   switch (GeomType)
   {
      case Geometry::POINT:
         ip.x = 0.0;
         break;
      case Geometry::SEGMENT:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         break;
      case Geometry::TRIANGLE:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         ip.y = real_t(rand()) / real_t(RAND_MAX);
         if (ip.x + ip.y > 1.0)
         {
            ip.x = 1.0 - ip.x;
            ip.y = 1.0 - ip.y;
         }
         break;
      case Geometry::SQUARE:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         ip.y = real_t(rand()) / real_t(RAND_MAX);
         break;
      case Geometry::TETRAHEDRON:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         ip.y = real_t(rand()) / real_t(RAND_MAX);
         ip.z = real_t(rand()) / real_t(RAND_MAX);
         // map to the triangular prism obtained by extruding the reference
         // triangle in z direction
         if (ip.x + ip.y > 1.0)
         {
            ip.x = 1.0 - ip.x;
            ip.y = 1.0 - ip.y;
         }
         // split the prism into 3 parts: 1 is the reference tet, and the
         // other two tets (as given below) are mapped to the reference tet
         if (ip.x + ip.z > 1.0)
         {
            // tet with vertices: (0,0,1),(1,0,1),(0,1,1),(1,0,0)
            ip.x = ip.x + ip.z - 1.0;
            // ip.y = ip.y;
            ip.z = 1.0 - ip.z;
            // mapped to: (0,0,0),(1,0,0),(0,1,0),(0,0,1)
         }
         else if (ip.x + ip.y + ip.z > 1.0)
         {
            // tet with vertices: (0,1,1),(0,1,0),(0,0,1),(1,0,0)
            real_t x = ip.x;
            ip.x = 1.0 - x - ip.z;
            ip.y = 1.0 - x - ip.y;
            ip.z = x;
            // mapped to: (0,0,0),(1,0,0),(0,1,0),(0,0,1)
         }
         break;
      case Geometry::CUBE:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         ip.y = real_t(rand()) / real_t(RAND_MAX);
         ip.z = real_t(rand()) / real_t(RAND_MAX);
         break;
      case Geometry::PRISM:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         ip.y = real_t(rand()) / real_t(RAND_MAX);
         ip.z = real_t(rand()) / real_t(RAND_MAX);
         if (ip.x + ip.y > 1.0)
         {
            ip.x = 1.0 - ip.x;
            ip.y = 1.0 - ip.y;
         }
         break;
      case Geometry::PYRAMID:
         ip.x = real_t(rand()) / real_t(RAND_MAX);
         ip.y = real_t(rand()) / real_t(RAND_MAX);
         ip.z = real_t(rand()) / real_t(RAND_MAX);
         if (ip.x + ip.z > 1.0 && ip.y < ip.x)
         {
            real_t x = ip.x;
            ip.x = ip.y;
            ip.y = 1.0 - ip.z;
            ip.z = 1.0 - x;
         }
         else if (ip.y + ip.z > 1.0)
         {
            real_t z = ip.z;
            ip.z = 1.0 - ip.y;
            ip.y = ip.x;
            ip.x = 1.0 - z;
         }
         break;
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }
}


namespace internal
{

// Fuzzy equality operator with absolute tolerance eps.
inline bool NearlyEqual(real_t x, real_t y, real_t eps)
{
   return std::abs(x-y) <= eps;
}

// Fuzzy greater than comparison operator with absolute tolerance eps.
// Returns true when x is greater than y by at least eps.
inline bool FuzzyGT(real_t x, real_t y, real_t eps)
{
   return (x > y + eps);
}

// Fuzzy less than comparison operator with absolute tolerance eps.
// Returns true when x is less than y by at least eps.
inline bool FuzzyLT(real_t x, real_t y, real_t eps)
{
   return (x < y - eps);
}

}

// static method
bool Geometry::CheckPoint(int64_t GeomType, const IntegrationPoint &ip)
{
   switch (GeomType)
   {
      case Geometry::POINT:
         if (ip.x != 0.0) { return false; }
         break;
      case Geometry::SEGMENT:
         if (ip.x < 0.0 || ip.x > 1.0) { return false; }
         break;
      case Geometry::TRIANGLE:
         if (ip.x < 0.0 || ip.y < 0.0 || ip.x+ip.y > 1.0) { return false; }
         break;
      case Geometry::SQUARE:
         if (ip.x < 0.0 || ip.x > 1.0 || ip.y < 0.0 || ip.y > 1.0)
         { return false; }
         break;
      case Geometry::TETRAHEDRON:
         if (ip.x < 0.0 || ip.y < 0.0 || ip.z < 0.0 ||
             ip.x+ip.y+ip.z > 1.0) { return false; }
         break;
      case Geometry::CUBE:
         if (ip.x < 0.0 || ip.x > 1.0 || ip.y < 0.0 || ip.y > 1.0 ||
             ip.z < 0.0 || ip.z > 1.0) { return false; }
         break;
      case Geometry::PRISM:
         if (ip.x < 0.0 || ip.y < 0.0 || ip.x+ip.y > 1.0 ||
             ip.z < 0.0 || ip.z > 1.0) { return false; }
         break;
      case Geometry::PYRAMID:
         if (ip.x < 0.0 || ip.y < 0.0 || ip.x+ip.z > 1.0 || ip.y+ip.z > 1.0 ||
             ip.z < 0.0 || ip.z > 1.0) { return false; }
         break;
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }
   return true;
}

// static method
bool Geometry::CheckPoint(int64_t GeomType, const IntegrationPoint &ip,
                          real_t eps)
{
   switch (GeomType)
   {
      case Geometry::POINT:
         if (! internal::NearlyEqual(ip.x, 0.0, eps))
         {
            return false;
         }
         break;
      case Geometry::SEGMENT:
         if ( internal::FuzzyLT(ip.x, 0.0, eps)
              || internal::FuzzyGT(ip.x, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::TRIANGLE:
         if ( internal::FuzzyLT(ip.x, 0.0, eps)
              || internal::FuzzyLT(ip.y, 0.0, eps)
              || internal::FuzzyGT(ip.x+ip.y, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::SQUARE:
         if ( internal::FuzzyLT(ip.x, 0.0, eps)
              || internal::FuzzyGT(ip.x, 1.0, eps)
              || internal::FuzzyLT(ip.y, 0.0, eps)
              || internal::FuzzyGT(ip.y, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::TETRAHEDRON:
         if ( internal::FuzzyLT(ip.x, 0.0, eps)
              || internal::FuzzyLT(ip.y, 0.0, eps)
              || internal::FuzzyLT(ip.z, 0.0, eps)
              || internal::FuzzyGT(ip.x+ip.y+ip.z, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::CUBE:
         if ( internal::FuzzyLT(ip.x, 0.0, eps)
              || internal::FuzzyGT(ip.x, 1.0, eps)
              || internal::FuzzyLT(ip.y, 0.0, eps)
              || internal::FuzzyGT(ip.y, 1.0, eps)
              || internal::FuzzyLT(ip.z, 0.0, eps)
              || internal::FuzzyGT(ip.z, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::PRISM:
         if ( internal::FuzzyLT(ip.x, 0.0, eps)
              || internal::FuzzyLT(ip.y, 0.0, eps)
              || internal::FuzzyGT(ip.x+ip.y, 1.0, eps)
              || internal::FuzzyLT(ip.z, 0.0, eps)
              || internal::FuzzyGT(ip.z, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::PYRAMID:
         if (internal::FuzzyLT(ip.x, 0.0, eps)
             || internal::FuzzyLT(ip.y, 0.0, eps)
             || internal::FuzzyGT(ip.x+ip.z, 1.0, eps)
             || internal::FuzzyGT(ip.y+ip.z, 1.0, eps)
             || internal::FuzzyLT(ip.z, 0.0, eps)
             || internal::FuzzyGT(ip.z, 1.0, eps) )
         {
            return false;
         }
         break;
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }
   return true;
}


namespace internal
{

template <int64_t N, int64_t dim>
inline bool IntersectSegment(real_t lbeg[N], real_t lend[N],
                             IntegrationPoint &end)
{
   real_t t = 1.0;
   bool out = false;
   for (int64_t i = 0; i < N; i++)
   {
      lbeg[i] = std::max(lbeg[i], (real_t) 0.0); // remove round-off
      if (lend[i] < 0.0)
      {
         out = true;
         t = std::min(t, lbeg[i]/(lbeg[i]-lend[i]));
      }
   }
   if (out)
   {
      if (dim >= 1) { end.x = t*lend[0] + (1.0-t)*lbeg[0]; }
      if (dim >= 2) { end.y = t*lend[1] + (1.0-t)*lbeg[1]; }
      if (dim >= 3) { end.z = t*lend[2] + (1.0-t)*lbeg[2]; }
      return false;
   }
   return true;
}

inline bool ProjectTriangle(real_t &x, real_t &y)
{
   if (x < 0.0)
   {
      x = 0.0;
      if (y < 0.0)      { y = 0.0; }
      else if (y > 1.0) { y = 1.0; }
      return false;
   }
   if (y < 0.0)
   {
      if (x > 1.0) { x = 1.0; }
      y = 0.0;
      return false;
   }
   const real_t l3 = 1.0-x-y;
   if (l3 < 0.0)
   {
      if (y - x > 1.0)       { x = 0.0; y = 1.0; }
      else if (y - x < -1.0) { x = 1.0; y = 0.0; }
      else                   { x += l3/2; y += l3/2; }
      return false;
   }
   return true;
}

}

// static method
bool Geometry::ProjectPoint(int64_t GeomType, const IntegrationPoint &beg,
                            IntegrationPoint &end)
{
   constexpr real_t fone = 1.0;

   switch (GeomType)
   {
      case Geometry::POINT:
      {
         if (end.x != 0.0) { end.x = 0.0; return false; }
         break;
      }
      case Geometry::SEGMENT:
      {
         if (end.x < 0.0) { end.x = 0.0; return false; }
         if (end.x > 1.0) { end.x = 1.0; return false; }
         break;
      }
      case Geometry::TRIANGLE:
      {
         real_t lend[3] = { end.x, end.y, fone-end.x-end.y };
         real_t lbeg[3] = { beg.x, beg.y, fone-beg.x-beg.y };
         return internal::IntersectSegment<3,2>(lbeg, lend, end);
      }
      case Geometry::SQUARE:
      {
         real_t lend[4] = { end.x, end.y, fone-end.x, fone-end.y };
         real_t lbeg[4] = { beg.x, beg.y, fone-beg.x, fone-beg.y };
         return internal::IntersectSegment<4,2>(lbeg, lend, end);
      }
      case Geometry::TETRAHEDRON:
      {
         real_t lend[4] = { end.x, end.y, end.z, fone-end.x-end.y-end.z };
         real_t lbeg[4] = { beg.x, beg.y, beg.z, fone-beg.x-beg.y-beg.z };
         return internal::IntersectSegment<4,3>(lbeg, lend, end);
      }
      case Geometry::CUBE:
      {
         real_t lend[6] = { end.x, end.y, end.z,
                            fone-end.x, fone-end.y, fone-end.z
                          };
         real_t lbeg[6] = { beg.x, beg.y, beg.z,
                            fone-beg.x, fone-beg.y, fone-beg.z
                          };
         return internal::IntersectSegment<6,3>(lbeg, lend, end);
      }
      case Geometry::PRISM:
      {
         real_t lend[5] = { end.x, end.y, end.z, fone-end.x-end.y, fone-end.z };
         real_t lbeg[5] = { beg.x, beg.y, beg.z, fone-beg.x-beg.y, fone-beg.z };
         return internal::IntersectSegment<5,3>(lbeg, lend, end);
      }
      case Geometry::PYRAMID:
      {
         real_t lend[6] = { end.x, end.y, end.z,
                            fone-end.x-end.z, fone-end.y-end.z, fone-end.z
                          };
         real_t lbeg[6] = { beg.x, beg.y, beg.z,
                            fone-beg.x-beg.z, fone-beg.y-beg.z, fone-beg.z
                          };
         return internal::IntersectSegment<6,3>(lbeg, lend, end);
      }
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }
   return true;
}

// static method
bool Geometry::ProjectPoint(int64_t GeomType, IntegrationPoint &ip)
{
   // If ip is outside the element, replace it with the point on the boundary
   // that is closest to the original ip and return false; otherwise, return
   // true without changing ip.

   switch (GeomType)
   {
      case SEGMENT:
      {
         if (ip.x < 0.0)      { ip.x = 0.0; return false; }
         else if (ip.x > 1.0) { ip.x = 1.0; return false; }
         return true;
      }

      case TRIANGLE:
      {
         return internal::ProjectTriangle(ip.x, ip.y);
      }

      case SQUARE:
      {
         bool in_x, in_y;
         if (ip.x < 0.0)      { in_x = false; ip.x = 0.0; }
         else if (ip.x > 1.0) { in_x = false; ip.x = 1.0; }
         else                 { in_x = true; }
         if (ip.y < 0.0)      { in_y = false; ip.y = 0.0; }
         else if (ip.y > 1.0) { in_y = false; ip.y = 1.0; }
         else                 { in_y = true; }
         return in_x && in_y;
      }

      case TETRAHEDRON:
      {
         if (ip.z < 0.0)
         {
            ip.z = 0.0;
            internal::ProjectTriangle(ip.x, ip.y);
            return false;
         }
         if (ip.y < 0.0)
         {
            ip.y = 0.0;
            internal::ProjectTriangle(ip.x, ip.z);
            return false;
         }
         if (ip.x < 0.0)
         {
            ip.x = 0.0;
            internal::ProjectTriangle(ip.y, ip.z);
            return false;
         }
         const real_t l4 = 1.0-ip.x-ip.y-ip.z;
         if (l4 < 0.0)
         {
            const real_t l4_3 = l4/3;
            ip.x += l4_3;
            ip.y += l4_3;
            internal::ProjectTriangle(ip.x, ip.y);
            ip.z = 1.0-ip.x-ip.y;
            return false;
         }
         return true;
      }

      case CUBE:
      {
         bool in_x, in_y, in_z;
         if (ip.x < 0.0)      { in_x = false; ip.x = 0.0; }
         else if (ip.x > 1.0) { in_x = false; ip.x = 1.0; }
         else                 { in_x = true; }
         if (ip.y < 0.0)      { in_y = false; ip.y = 0.0; }
         else if (ip.y > 1.0) { in_y = false; ip.y = 1.0; }
         else                 { in_y = true; }
         if (ip.z < 0.0)      { in_z = false; ip.z = 0.0; }
         else if (ip.z > 1.0) { in_z = false; ip.z = 1.0; }
         else                 { in_z = true; }
         return in_x && in_y && in_z;
      }

      case PRISM:
      {
         bool in_tri, in_z;
         in_tri = internal::ProjectTriangle(ip.x, ip.y);
         if (ip.z < 0.0)      { in_z = false; ip.z = 0.0; }
         else if (ip.z > 1.0) { in_z = false; ip.z = 1.0; }
         else                 { in_z = true; }
         return in_tri && in_z;
      }

      case PYRAMID:
      {
         if (ip.x < 0.0)
         {
            ip.x = 0.0;
            internal::ProjectTriangle(ip.y, ip.z);
            return false;
         }
         if (ip.y < 0.0)
         {
            ip.y = 0.0;
            internal::ProjectTriangle(ip.x, ip.z);
            return false;
         }
         if (ip.z < 0.0)
         {
            ip.z = 0.0;
            if (ip.x > 1.0) { ip.x = 1.0; }
            if (ip.y > 1.0) { ip.y = 1.0; }
            return false;
         }
         if (ip.x >= ip.y)
         {
            bool in_y = true;
            bool in_tri = internal::ProjectTriangle(ip.x, ip.z);
            if (ip.y > ip.z) { in_y = false; ip.y = ip.z; }
            return in_tri && in_y;
         }
         else
         {
            bool in_x = true;
            bool in_tri = internal::ProjectTriangle(ip.y, ip.z);
            if (ip.x > ip.z) { in_x = false; ip.x = ip.z; }
            return in_tri && in_x;
         }
      }

      case Geometry::POINT:
         MFEM_ABORT("Reference element type is not supported!");
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }
   return true;
}

const int64_t Geometry::NumBdrArray[NumGeom] = { 0, 2, 3, 4, 4, 6, 5, 5 };
const int64_t Geometry::Dimension[NumGeom] = { 0, 1, 2, 2, 3, 3, 3, 3 };
const int64_t Geometry::DimStart[MaxDim+2] =
{ POINT, SEGMENT, TRIANGLE, TETRAHEDRON, NUM_GEOMETRIES };
const int64_t Geometry::NumVerts[NumGeom] = { 1, 2, 3, 4, 4, 8, 6, 5 };
const int64_t Geometry::NumEdges[NumGeom] = { 0, 1, 3, 4, 6, 12, 9, 8 };
const int64_t Geometry::NumFaces[NumGeom] = { 0, 0, 1, 1, 4, 6, 5, 5 };

const int64_t Geometry::
Constants<Geometry::POINT>::Orient[1][1] = {{0}};
const int64_t Geometry::
Constants<Geometry::POINT>::InvOrient[1] = {0};

const int64_t Geometry::
Constants<Geometry::SEGMENT>::Edges[1][2] = { {0, 1} };
const int64_t Geometry::
Constants<Geometry::SEGMENT>::Orient[2][2] = { {0, 1}, {1, 0} };
const int64_t Geometry::
Constants<Geometry::SEGMENT>::InvOrient[2] = { 0, 1 };

const int64_t Geometry::
Constants<Geometry::TRIANGLE>::Edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
const int64_t Geometry::
Constants<Geometry::TRIANGLE>::VertToVert::I[3] = {0, 2, 3};
const int64_t Geometry::
Constants<Geometry::TRIANGLE>::VertToVert::J[3][2] = {{1, 0}, {2, -3}, {2, 1}};
const int64_t Geometry::
Constants<Geometry::TRIANGLE>::FaceVert[1][3] = {{0, 1, 2}};
const int64_t Geometry::
Constants<Geometry::TRIANGLE>::Orient[6][3] =
{
   {0, 1, 2}, {1, 0, 2}, {2, 0, 1},
   {2, 1, 0}, {1, 2, 0}, {0, 2, 1}
};
const int64_t Geometry::
Constants<Geometry::TRIANGLE>::InvOrient[6] = { 0, 1, 4, 3, 2, 5 };

const int64_t Geometry::
Constants<Geometry::SQUARE>::Edges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
const int64_t Geometry::
Constants<Geometry::SQUARE>::VertToVert::I[4] = {0, 2, 3, 4};
const int64_t Geometry::
Constants<Geometry::SQUARE>::VertToVert::J[4][2] =
{{1, 0}, {3, -4}, {2, 1}, {3, 2}};
const int64_t Geometry::
Constants<Geometry::SQUARE>::FaceVert[1][4] = {{0, 1, 2, 3}};
const int64_t Geometry::
Constants<Geometry::SQUARE>::Orient[8][4] =
{
   {0, 1, 2, 3}, {0, 3, 2, 1}, {1, 2, 3, 0}, {1, 0, 3, 2},
   {2, 3, 0, 1}, {2, 1, 0, 3}, {3, 0, 1, 2}, {3, 2, 1, 0}
};
const int64_t Geometry::
Constants<Geometry::SQUARE>::InvOrient[8] = { 0, 1, 6, 3, 4, 5, 2, 7 };

const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::Edges[6][2] =
{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::FaceTypes[4] =
{
   Geometry::TRIANGLE, Geometry::TRIANGLE,
   Geometry::TRIANGLE, Geometry::TRIANGLE
};
const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::FaceVert[4][3] =
{{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::VertToVert::I[4] = {0, 3, 5, 6};
const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::VertToVert::J[6][2] =
{
   {1, 0}, {2, 1}, {3, 2}, // 0,1:0   0,2:1   0,3:2
   {2, 3}, {3, 4},         // 1,2:3   1,3:4
   {3, 5}                  // 2,3:5
};
const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::Orient[24][4] =
{
   {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {0, 2, 1, 3},
   {0, 3, 1, 2}, {0, 3, 2, 1},
   {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 2, 0}, {1, 3, 0, 2},
   {1, 0, 3, 2}, {1, 0, 2, 3},
   {2, 3, 0, 1}, {2, 3, 1, 0}, {2, 0, 1, 3}, {2, 0, 3, 1},
   {2, 1, 3, 0}, {2, 1, 0, 3},
   {3, 0, 2, 1}, {3, 0, 1, 2}, {3, 1, 0, 2}, {3, 1, 2, 0},
   {3, 2, 1, 0}, {3, 2, 0, 1}
};
const int64_t Geometry::
Constants<Geometry::TETRAHEDRON>::InvOrient[24] =
{
   0,   1,  4,  3,  2,  5,
   14, 19, 18, 15, 10, 11,
   12, 23,  6,  9, 20, 17,
   8,   7, 16, 21, 22, 13
};

const int64_t Geometry::
Constants<Geometry::CUBE>::Edges[12][2] =
{
   {0, 1}, {1, 2}, {3, 2}, {0, 3}, {4, 5}, {5, 6},
   {7, 6}, {4, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}
};
const int64_t Geometry::
Constants<Geometry::CUBE>::FaceTypes[6] =
{
   Geometry::SQUARE, Geometry::SQUARE, Geometry::SQUARE,
   Geometry::SQUARE, Geometry::SQUARE, Geometry::SQUARE
};
const int64_t Geometry::
Constants<Geometry::CUBE>::FaceVert[6][4] =
{
   {3, 2, 1, 0}, {0, 1, 5, 4}, {1, 2, 6, 5},
   {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}
};
const int64_t Geometry::
Constants<Geometry::CUBE>::VertToVert::I[8] = {0, 3, 5, 7, 8, 10, 11, 12};
const int64_t Geometry::
Constants<Geometry::CUBE>::VertToVert::J[12][2] =
{
   {1, 0}, {3, 3}, {4, 8}, // 0,1:0   0,3:3   0,4:8
   {2, 1}, {5, 9},         // 1,2:1   1,5:9
   {3,-3}, {6,10},         // 2,3:-3  2,6:10
   {7,11},                 // 3,7:11
   {5, 4}, {7, 7},         // 4,5:4   4,7:7
   {6, 5},                 // 5,6:5
   {7,-7}                  // 6,7:-7
};

const int64_t Geometry::
Constants<Geometry::PRISM>::Edges[9][2] =
{{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
const int64_t Geometry::
Constants<Geometry::PRISM>::FaceTypes[5] =
{
   Geometry::TRIANGLE, Geometry::TRIANGLE,
   Geometry::SQUARE, Geometry::SQUARE, Geometry::SQUARE
};
const int64_t Geometry::
Constants<Geometry::PRISM>::FaceVert[5][4] =
{{0, 2, 1, -1}, {3, 4, 5, -1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
const int64_t Geometry::
Constants<Geometry::PRISM>::VertToVert::I[6] = {0, 3, 5, 6, 8, 9};
const int64_t Geometry::
Constants<Geometry::PRISM>::VertToVert::J[9][2] =
{
   {1, 0}, {2, -3}, {3, 6}, // 0,1:0   0,2:-3  0,3:6
   {2, 1}, {4, 7},          // 1,2:1   1,4:7
   {5, 8},                  // 2,5:8
   {4, 3}, {5, -6},         // 3,4:3   3,5:-6
   {5, 4}                   // 4,5:4
};

const int64_t Geometry::
Constants<Geometry::PYRAMID>::Edges[8][2] =
{{0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
const int64_t Geometry::
Constants<Geometry::PYRAMID>::FaceTypes[5] =
{
   Geometry::SQUARE,
   Geometry::TRIANGLE, Geometry::TRIANGLE,
   Geometry::TRIANGLE, Geometry::TRIANGLE
};
const int64_t Geometry::
Constants<Geometry::PYRAMID>::FaceVert[5][4] =
{{3, 2, 1, 0}, {0, 1, 4, -1}, {1, 2, 4, -1}, {2, 3, 4, -1}, {3, 0, 4, -1}};
const int64_t Geometry::
Constants<Geometry::PYRAMID>::VertToVert::I[5] = {0, 3, 5, 7, 8};
const int64_t Geometry::
Constants<Geometry::PYRAMID>::VertToVert::J[8][2] =
{
   {1, 0}, {3, 3}, {4, 4}, // 0,1:0   0,3:3  0,4:4
   {2, 1}, {4, 5},         // 1,2:1   1,4:5
   {3,-3}, {4, 6},         // 2,3:-3  2,4:6
   {4, 7}                  // 3,4:7
};

Geometry Geometries;

}
