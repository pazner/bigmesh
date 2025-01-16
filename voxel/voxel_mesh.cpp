#include "voxel_mesh.hpp"
#include "geom.hpp"

namespace mfem
{

int64_t FindComponents(const Table &elem_elem, Array<int64_t> &component)
{
   const int64_t num_elem = elem_elem.Size();
   const int64_t *i_elem_elem = elem_elem.GetI();
   const int64_t *j_elem_elem = elem_elem.GetJ();

   component.SetSize(num_elem);
   component = -1;

   Array<int64_t> elem_stack(num_elem);

   int64_t num_comp = 0;

   int64_t stack_p = 0;
   int64_t stack_top_p = 0; // points to the first unused element in the stack
   for (int64_t elem = 0; elem < num_elem; elem++)
   {
      if (component[elem] >= 0) { continue; }

      component[elem] = num_comp;
      ++num_comp;

      elem_stack[stack_top_p++] = elem;

      for ( ; stack_p < stack_top_p; stack_p++)
      {
         const int64_t i = elem_stack[stack_p];
         for (int64_t j = i_elem_elem[i]; j < i_elem_elem[i+1]; j++)
         {
            const int64_t k = j_elem_elem[j];
            if (component[k] < 0)
            {
               component[k] = component[i];
               elem_stack[stack_top_p++] = k;
            }
            else if (component[k] != component[i])
            {
               MFEM_ABORT("");
            }
         }
      }
   }
   return num_comp;
}

VoxelMesh::VoxelMesh(const std::string &filename) : Mesh(filename)
{
   MFEM_VERIFY(GetNE() > 0, "Empty mesh not supported.");

   // h = GetElementSize(0);
   {
      Array<int64_t> verts;
      GetElementVertices(0, verts);
      h = 0.0;
      for (int d = 0; d < spaceDim; ++d)
      {
         h += std::abs(vertices[verts[1]](d) - vertices[verts[0]](d));
      }
   }
   // SetCurvature(0);

   Array<int64_t> components;
   const int64_t n_components = FindComponents(ElementToElementTable(),
                                               components);
   std::cout << "Found " << n_components << " connected components.\n";
   if (n_components > 1)
   {
      Array<int64_t> component_count(n_components);
      component_count = 0;
      for (int64_t c : components) { ++component_count[c]; }

      const auto max_it = std::max_element(component_count.begin(),
                                           component_count.end());
      const int64_t c = std::distance(component_count.begin(), max_it);

      Array<Element*> new_elements(component_count[c]);
      int64_t j = 0;
      for (int64_t i = 0; i < GetNE(); ++i)
      {
         if (components[i] == c)
         {
            new_elements[j] = elements[i];
            ++j;
         }
         else
         {
            FreeElement(elements[i]);
         }
      }
      mfem::Swap(elements, new_elements);
      NumOfElements = j;

      // TODO: retain boundary attributes
      boundary.DeleteAll();
      NumOfBdrElements = 0;

      DeleteTables();
      RemoveUnusedVertices();
      FinalizeMesh();
   }

   Vector mmin, mmax;
   GetBoundingBox(mmin, mmax, 0);

   n.resize(Dim);
   for (int64_t i = 0; i < Dim; ++i) { n[i] = std::ceil((mmax[i] - mmin[i])/h); }

   for (int i = 0; i < NumOfVertices; i++)
   {
      for (int d = 0; d < spaceDim; d++)
      {
         vertices[i](d) -= mmin[d];
      }
   }


   Vector center;
   std::vector<int64_t> center_idx(Dim);
   idx2lex.resize(NumOfElements);
   for (int64_t i = 0; i < NumOfElements; ++i)
   {
      GetElementCenter(i, center);
      for (int64_t d = 0; d < Dim; ++d)
      {
         center_idx[d] = std::floor(center[d] / h);
      }

      LexIndex lex(center_idx.data(), Dim);
      const int64_t lin_idx = lex.LinearIndex(n);

      idx2lex[i] = lex;
      MFEM_VERIFY(lex2idx.find(lin_idx) == lex2idx.end(), "");
      lex2idx[lin_idx] = i;
   }
}

VoxelMesh::VoxelMesh(double h_, const std::vector<int64_t> &n_)
   : Mesh(n_.size(), 0, 0), h(h_), n(n_) { }

VoxelMesh VoxelMesh::Coarsen() const
{
   double new_h = 2.0*h;
   std::vector<int64_t> new_n(Dim);
   for (int64_t i = 0; i < Dim; ++i)
   {
      new_n[i] = static_cast<int64_t>(std::ceil(0.5*n[i]));
   }
   VoxelMesh coarsened_mesh(new_h, new_n);

   std::vector<int64_t> new_n_vert(Dim);
   for (int64_t i = 0; i < Dim; ++i) { new_n_vert[i] = new_n[i] + 1; }

   std::unordered_map<int64_t,int64_t> lex2vert;

   auto coords_to_lex = [&](const std::array<double,3> &coords)
   {
      LexIndex idx;
      idx.ndim = Dim;
      for (int64_t d = 0; d < Dim; ++d)
      {
         idx.coords[d] = static_cast<int64_t>(std::round(coords[d] / h));
      }
      return idx;
   };

   auto maybe_add_vertex = [&](const LexIndex &lex,
                               const std::array<double,3> &coords)
   {
      const int64_t lin_idx = lex.LinearIndex(new_n_vert);
      if (lex2vert.find(lin_idx) == lex2vert.end())
      {
         lex2vert[lin_idx] = coarsened_mesh.NumOfVertices;
         coarsened_mesh.AddVertex(coords.data());
      }
   };

   std::array<double,3> v;

   // Add every even-indexed vertex and all of its neighbors
   for (int64_t iv = 0; iv < NumOfVertices; ++iv)
   {
      // Populate v with the coordinates of vertex iv
      {
         const double *v_iv = vertices[iv]();
         std::copy(v_iv, v_iv + Dim, v.begin());
      }

      // Retain only the even-numbered vertices
      LexIndex lex = coords_to_lex(v);
      bool even = true;
      for (int64_t d = 0; d < Dim; ++d)
      {
         if (lex.coords[d] % 2 != 0)
         {
            even = false;
            break;
         }
      }
      if (!even) { continue; }

      // Add the vertex to the coarse mesh (if it doesn't exist already)
      for (int64_t d = 0; d < Dim; ++d) { lex.coords[d] /= 2; }
      maybe_add_vertex(lex, v);

      // Add all of its neighbors to the coarse mesh (if they don't exist)
      const int64_t ngrid = std::pow(3, Dim);

      std::array<int64_t,3> shift;
      for (int64_t i = 0; i < ngrid; ++i)
      {
         int64_t j = i;
         bool in_bounds = true;
         for (int64_t d = 0; d < Dim; ++d)
         {
            shift[d] = (j % 3) - 1;
            j /= 3;

            if (lex.coords[d] + shift[d] < 0 ||
                lex.coords[d] + shift[d] >= new_n_vert[d])
            {
               in_bounds = false; break;
            }
         }

         if (!in_bounds) { continue; }

         for (int64_t d = 0; d < Dim; ++d)
         {
            v[d] += shift[d]*new_h;
            lex.coords[d] += shift[d];
         }

         // If this vertex hasn't been added yet, add it.
         maybe_add_vertex(lex, v);

         // Reset (shift back)
         for (int64_t d = 0; d < Dim; ++d)
         {
            v[d] -= shift[d]*new_h;
            lex.coords[d] -= shift[d];
         }
      }
   }

   // Get the vertex integration rule for the geometry
   const IntegrationRule &ir = *Geometries.GetVertices(GetElementGeometry(0));

   // Having added all the (potential) coarse mesh vertices, we now add the
   // coarse mesh elements
   for (int64_t ie = 0; ie < NumOfElements; ++ie)
   {
      LexIndex lex = idx2lex[ie];
      // Figure out which macro element we're in
      for (int64_t d = 0; d < Dim; ++d) { lex.coords[d] /= 2; }

      const int64_t lin_idx = lex.LinearIndex(new_n);
      if (coarsened_mesh.lex2idx.find(lin_idx) == coarsened_mesh.lex2idx.end())
      {
         const int64_t attr = GetAttribute(ie);
         coarsened_mesh.lex2idx[lin_idx] = coarsened_mesh.NumOfElements;
         coarsened_mesh.idx2lex.push_back(lex);

         std::vector<int64_t> el_vert(ir.Size());
         for (int64_t iv = 0; iv < ir.Size(); ++iv)
         {
            std::array<double,3> ip;
            ir[iv].Get(ip.data(), Dim);
            for (int64_t d = 0; d < Dim; ++d) { lex.coords[d] += ip[d]; }
            const int64_t v_lin_idx = lex.LinearIndex(new_n_vert);
            el_vert[iv] = lex2vert.at(v_lin_idx);
            // Reset
            for (int64_t d = 0; d < Dim; ++d) { lex.coords[d] -= ip[d]; }
         }

         if (Dim == 1) { MFEM_ABORT("To be implemented"); }
         else if (Dim == 2) { coarsened_mesh.AddQuad(el_vert.data(), attr); }
         else if (Dim == 3) { coarsened_mesh.AddHex(el_vert.data(), attr); }
         else { MFEM_ABORT("Unsupported dimension."); }
      }
   }

   coarsened_mesh.RemoveUnusedVertices();
   coarsened_mesh.FinalizeMesh();

   const int64_t ngrid = pow(2,
                             Dim); // number of fine elements per coarse elements
   // Inherit boundary attributes from parents
   std::unordered_map<int64_t,std::unordered_map<int64_t,int64_t>> el_to_bdr_attr;
   for (int64_t ib = 0; ib < GetNBE(); ++ib)
   {
      const int64_t attr = GetBdrAttribute(ib);
      int64_t ie, info;
      GetBdrElementAdjacentElement(ib, ie, info);
      const int64_t local_side = info/64; // decode info
      el_to_bdr_attr[ie][local_side] = attr;
   }

   for (int64_t ib = 0; ib < coarsened_mesh.GetNBE(); ++ib)
   {
      int64_t attr = -1;
      int64_t ie, info;
      coarsened_mesh.GetBdrElementAdjacentElement(ib, ie, info);
      const int64_t local_side = info/64;
      LexIndex lex = coarsened_mesh.GetLexicographicIndex(ie);
      LexIndex shifted_lex = lex;
      // Convert from coarse lexicographic index to fine index
      for (int64_t d = 0; d < Dim; ++d)
      {
         lex.coords[d] *= 2;
      }
      for (int64_t i = 0; i < ngrid; ++i)
      {
         int64_t j = i;
         std::array<int64_t,3> shift;
         for (int64_t d = 0; d < Dim; ++d)
         {
            shift[d] = j % 2;
            j /= 2;

            shifted_lex.coords[d] = lex.coords[d] + shift[d];
         }
         const int64_t fine_idx = GetElementIndex(shifted_lex);
         const auto result_1 = el_to_bdr_attr.find(fine_idx);

         if (result_1 != el_to_bdr_attr.end())
         {
            const auto &local_attr_map = result_1->second;
            const auto result_2 = local_attr_map.find(local_side);
            if (result_2 != local_attr_map.end())
            {
               if (attr >= 0)
               {
                  attr = std::min(attr, result_2->second);
               }
               else
               {
                  attr = result_2->second;
               }
            }
         }
      }
      MFEM_VERIFY(attr >= 0, "");
      coarsened_mesh.SetBdrAttribute(ib, attr);
   }

   return coarsened_mesh;
}

void GetVoxelParents(const VoxelMesh &coarse_mesh, const VoxelMesh &fine_mesh,
                     Array<ParentIndex> &parents, Array<int64_t> &parent_offsets)
{
   const int64_t dim = coarse_mesh.Dimension();
   const int64_t coarse_ne = coarse_mesh.GetNE();
   const int64_t ngrid = pow(2,
                             dim); // number of fine elements per coarse elements

   parent_offsets.SetSize(coarse_ne + 1);
   parents.SetSize(ngrid*coarse_ne);

   int64_t offset = 0;
   for (int64_t ie = 0; ie < coarse_ne; ++ie)
   {
      parent_offsets[ie] = offset;
      LexIndex lex = coarse_mesh.GetLexicographicIndex(ie);

      // Convert from coarse lexicographic index to fine index
      for (int64_t d = 0; d < dim; ++d)
      {
         lex.coords[d] *= 2;
      }

      for (int64_t i = 0; i < ngrid; ++i)
      {
         int64_t j = i;
         std::array<int64_t,3> shift;
         for (int64_t d = 0; d < dim; ++d)
         {
            shift[d] = j % 2;
            j /= 2;

            lex.coords[d] += shift[d];
         }

         const int64_t fine_idx = fine_mesh.GetElementIndex(lex);
         if (fine_idx >= 0)
         {
            parents[offset] = {fine_idx, i};
            ++offset;
         }

         // Reset
         for (int64_t d = 0; d < dim; ++d)
         {
            lex.coords[d] -= shift[d];
         }
      }
   }
   parent_offsets.Last() = offset;
   parents.SetSize(offset);
}

}
