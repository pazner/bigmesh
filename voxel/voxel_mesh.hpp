#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "mesh.hpp"
#include "lex.hpp"

namespace mfem
{

class VoxelMesh : public Mesh
{
protected:
   double h;
   std::vector<int64_t> n;
   std::unordered_map<int64_t,int64_t> lex2idx;
   std::vector<LexIndex> idx2lex;

   VoxelMesh(double h_, const std::vector<int64_t> &n_);

public:
   VoxelMesh(const std::string &filename);

   VoxelMesh Coarsen() const;

   int64_t GetElementIndex(const LexIndex &lex) const
   {
      const auto result = lex2idx.find(lex.LinearIndex(n));
      if (result != lex2idx.end()) { return result->second; }
      return -1; // invalid index: element not found
   }
   LexIndex GetLexicographicIndex(int64_t idx) const { return idx2lex[idx]; }

   const std::vector<int64_t> &GetVoxelBounds() const { return n; }
};

struct ParentIndex
{
   int64_t element_index;
   int64_t pmat_index;
};

void GetVoxelParents(const VoxelMesh &coarse_mesh, const VoxelMesh &fine_mesh,
                     Array<ParentIndex> &parents, Array<int64_t> &parent_offsets);

}

#endif
