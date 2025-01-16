#ifndef PAR_MG_HPP
#define PAR_MG_HPP

#include "voxel_mesh.hpp"

namespace mfem
{

struct CoarseToFineCommunication
{
   struct CoarseToFineIndex
   {
      int64_t coarse_element_index;
      int64_t pmat_index;
   };
   int64_t rank; // rank owning the fine elements
   std::vector<CoarseToFineIndex> coarse_to_fine;

   CoarseToFineCommunication() = default;
   CoarseToFineCommunication(int64_t rank_) : rank(rank_) { }
};

struct FineToCoarseCommunication
{
   struct FineToCoarseIndex
   {
      int64_t fine_element_index;
      int64_t pmat_index;
   };
   int64_t rank; // rank owning the coarse elements
   std::vector<FineToCoarseIndex> fine_to_coarse;

   FineToCoarseCommunication() = default;
   FineToCoarseCommunication(int64_t rank_) : rank(rank_) { }
};

struct ParVoxelMapping
{
   Array<ParentIndex> local_parents;
   Array<int64_t> local_parent_offsets;

   std::vector<CoarseToFineCommunication> coarse_to_fine;
   std::vector<FineToCoarseCommunication> fine_to_coarse;
};

std::vector<ParVoxelMapping> CreateParVoxelMappings(
   const int64_t nranks,
   const int64_t dim,
   const Array<ParentIndex> &parents,
   const Array<int64_t> &parent_offsets,
   const Array<int64_t> &fine_partitioning,
   const Array<int64_t> &coarse_partitioning);

}

#endif
