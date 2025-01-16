#include "voxel_mesh.hpp"
#include "optparser.hpp"
#include "tic_toc.hpp"
#include "par_mg.hpp"
#include "picojson.h"

#include <sys/stat.h>  // mkdir

using namespace std;
using namespace mfem;

void SavePartitionedMesh(const string &prefix,
                         Mesh &mesh,
                         int64_t np,
                         int64_t *partitioning)
{
   MeshPartitioner partitioner(mesh, np, partitioning);
   MeshPart mesh_part;
   for (int64_t i = 0; i < np; i++)
   {
      partitioner.ExtractPart(i, mesh_part);

      ofstream f(MakeParFilename(prefix + ".mesh.", i));
      f.precision(16);
      mesh_part.Print(f);
   }
}

void SaveParVoxelMappings(const string &prefix,
                          const vector<ParVoxelMapping> &mappings)
{
   const int64_t np = mappings.size();
   for (int64_t i = 0; i < np; ++i)
   {
      ofstream f(MakeParFilename(prefix + ".mapping.", i));

      f << mappings[i].local_parents.Size() << '\n';
      for (const auto &p : mappings[i].local_parents)
      {
         f << p.element_index << '\n' << p.pmat_index << '\n';
      }
      f << '\n';

      f << mappings[i].local_parent_offsets.Size() << '\n';
      for (const int64_t x : mappings[i].local_parent_offsets) { f << x << '\n'; }
      f << '\n';

      f << mappings[i].coarse_to_fine.size() << '\n';
      for (const auto &c2f : mappings[i].coarse_to_fine)
      {
         f << c2f.rank << '\n';
         f << c2f.coarse_to_fine.size() << '\n';
         for (const auto &x : c2f.coarse_to_fine)
         {
            f << x.coarse_element_index << '\n' << x.pmat_index << '\n';
         }
      }
      f << '\n';

      f << mappings[i].fine_to_coarse.size() << '\n';
      for (const auto &f2c : mappings[i].fine_to_coarse)
      {
         f << f2c.rank << '\n';
         f << f2c.fine_to_coarse.size() << '\n';
         for (const auto &x : f2c.fine_to_coarse)
         {
            f << x.fine_element_index << '\n' << x.pmat_index << '\n';
         }
      }
   }
}

int main(int argc, char *argv[])
{
   string mesh_file =
      "/Users/pazner/Documents/portland_state/10_research/13_meshes/bone_72k.mesh";
   string dir = "Voxel";
   int np = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&np, "-np", "--n-partitions", "Number of mesh partitions.");
   args.AddOption(&dir, "-d", "--dir", "Data directory.");
   args.ParseCheck();

   int err_flag = mkdir(dir.c_str(), 0777);
   err_flag = (err_flag && (errno != EEXIST)) ? 1 : 0;
   MFEM_VERIFY(err_flag == 0, "Could not create directory " << dir)

   cout << "\n";
   tic_toc.Restart();
   cout << "Reading fine mesh... " << flush;
   // Read the (fine) mesh from a file
   unique_ptr<VoxelMesh> mesh(new VoxelMesh(mesh_file));
   cout << "Done. " << tic_toc.RealTime() << endl;
   const int64_t dim = mesh->Dimension();

   // Partition the fine mesh
   tic_toc.Restart();
   cout << "Partitioning... " << flush;
   Array<int64_t> partitioning(mesh->GeneratePartitioning(np), mesh->GetNE(),
                               true);
   cout << "Done. " << tic_toc.RealTime() << endl;

   cout << "\n";
   cout << "Level      Elements     Save Mesh      Coarsening      Hierarchy      Mapping"
        << '\n'
        << string(77, '=')
        << endl;

   cout << left << setprecision(5) << fixed;

   int64_t level = 0;
   while (true)
   {
      cout << setw(11) << level << flush;
      cout << setw(13) << mesh->GetNE() << flush;

      tic_toc.Restart();
      const string level_str = dir + "/level_" + to_string(level);
      SavePartitionedMesh(level_str, *mesh, np, partitioning);
      cout << setw(15) << tic_toc.RealTime() << flush;

      const vector<int64_t> &bounds = mesh->GetVoxelBounds();
      if (!all_of(bounds.begin(), bounds.end(), [](int64_t x) { return x >= 4; }))
      {
         break;
      }

      tic_toc.Restart();
      unique_ptr<VoxelMesh> new_mesh(new VoxelMesh(mesh->Coarsen()));
      cout << setw(16) << tic_toc.RealTime() << flush;

      // Get hierarchy information for the new level
      tic_toc.Restart();
      Array<ParentIndex> parents;
      Array<int64_t> parent_offsets;
      GetVoxelParents(*new_mesh, *mesh, parents, parent_offsets);
      // Coarsen the partitioning
      Array<int64_t> new_partitioning(new_mesh->GetNE());
      for (int64_t i = 0; i < new_mesh->GetNE(); ++i)
      {
         const int64_t parent_index = parents[parent_offsets[i]].element_index;
         new_partitioning[i] = partitioning[parent_index];
      }
      cout << setw(15) << tic_toc.RealTime() << flush;

      // Create and save the parallel mappings
      tic_toc.Restart();
      auto mappings = CreateParVoxelMappings(
                         np, dim, parents, parent_offsets, partitioning,
                         new_partitioning);
      SaveParVoxelMappings(level_str, mappings);
      cout << setw(16) << tic_toc.RealTime() << endl;

      swap(new_mesh, mesh);
      Swap(new_partitioning, partitioning);

      ++level;
   }

   cout << endl;

   {
      picojson::object info;
      info["np"] = picojson::value(double(np));
      info["nlevels"] = picojson::value(double(level + 1));
      ofstream f(dir + "/info.json");
      f << picojson::value(info) << '\n';
   }

   return 0;
}
