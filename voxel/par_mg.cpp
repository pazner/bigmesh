// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "par_mg.hpp"
// #include "voxel_integ.hpp"
// #include "../../fem/picojson.h"

namespace mfem
{

std::vector<ParVoxelMapping> CreateParVoxelMappings(
   const int64_t nranks,
   const int64_t dim,
   const Array<ParentIndex> &parents,
   const Array<int64_t> &parent_offsets,
   const Array<int64_t> &fine_partitioning,
   const Array<int64_t> &coarse_partitioning)
{
   std::vector<ParVoxelMapping> mappings(nranks);

   std::vector<std::vector<int64_t>> local_coarse_elements(nranks);
   std::vector<std::unordered_map<int64_t,int64_t>> global_to_local_coarse(nranks);
   for (int64_t i = 0; i < coarse_partitioning.Size(); ++i)
   {
      const int64_t rank = coarse_partitioning[i];
      global_to_local_coarse[rank][i] = local_coarse_elements[rank].size();
      local_coarse_elements[rank].push_back(i);
   }

   std::vector<std::vector<int64_t>> local_fine_elements(nranks);
   std::vector<std::unordered_map<int64_t,int64_t>> global_to_local_fine(nranks);
   for (int64_t i = 0; i < fine_partitioning.Size(); ++i)
   {
      const int64_t rank = fine_partitioning[i];
      global_to_local_fine[rank][i] = local_fine_elements[rank].size();
      local_fine_elements[rank].push_back(i);
   }

   const int64_t ngrid = pow(2, dim);

   for (int64_t i = 0; i < nranks; ++i)
   {
      const int64_t n_local_coarse = local_coarse_elements[i].size();
      mappings[i].local_parent_offsets.SetSize(n_local_coarse + 1);
      mappings[i].local_parents.SetSize(ngrid*n_local_coarse);
   }

   std::vector<int64_t> offsets(nranks, 0);
   std::vector<std::unordered_map<int64_t,int64_t>> c2f_ranks(nranks);
   std::vector<std::unordered_map<int64_t,int64_t>> f2c_ranks(nranks);
   for (int64_t i = 0; i < coarse_partitioning.Size(); ++i)
   {
      const int64_t coarse_rank = coarse_partitioning[i];
      const int64_t local_coarse = global_to_local_coarse[coarse_rank].at(i);

      mappings[coarse_rank].local_parent_offsets[local_coarse] = offsets[coarse_rank];

      for (int64_t j = parent_offsets[i]; j < parent_offsets[i+1]; ++j)
      {
         const auto &p = parents[j];
         const int64_t fine_rank = fine_partitioning[p.element_index];
         const int64_t local_fine = global_to_local_fine[fine_rank].at(p.element_index);

         if (coarse_rank == fine_rank)
         {
            const int64_t rank = coarse_rank; // same as fine_rank
            mappings[coarse_rank].local_parents[offsets[rank]] = {local_fine, p.pmat_index};
            ++offsets[rank];
         }
         else
         {
            // Create the coarse to fine mapping
            {
               const auto result = c2f_ranks[coarse_rank].find(fine_rank);
               CoarseToFineCommunication *c2f = nullptr;
               if (result != c2f_ranks[coarse_rank].end())
               {
                  c2f = &mappings[coarse_rank].coarse_to_fine[result->second];
               }
               else
               {
                  c2f_ranks[coarse_rank][fine_rank] = mappings[coarse_rank].coarse_to_fine.size();
                  mappings[coarse_rank].coarse_to_fine.emplace_back(fine_rank);
                  c2f = &mappings[coarse_rank].coarse_to_fine.back();
               }
               c2f->coarse_to_fine.push_back({local_coarse, p.pmat_index});
            }

            // Create the fine to coarse mapping
            {
               const auto result = f2c_ranks[fine_rank].find(coarse_rank);
               FineToCoarseCommunication *f2c = nullptr;
               if (result != f2c_ranks[fine_rank].end())
               {
                  f2c = &mappings[fine_rank].fine_to_coarse[result->second];
               }
               else
               {
                  f2c_ranks[fine_rank][coarse_rank] = mappings[fine_rank].fine_to_coarse.size();
                  mappings[fine_rank].fine_to_coarse.emplace_back(coarse_rank);
                  f2c = &mappings[fine_rank].fine_to_coarse.back();
               }
               f2c->fine_to_coarse.push_back({local_fine, p.pmat_index});
            }
         }
      }
   }

   for (int64_t i = 0; i < nranks; ++i)
   {
      mappings[i].local_parent_offsets.Last() = offsets[i];
      mappings[i].local_parents.SetSize(offsets[i]);
   }

   return mappings;
}

}
