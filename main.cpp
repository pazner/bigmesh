#include "mesh.hpp"
#include "optparser.hpp"
#include "tic_toc.hpp"

#include <sys/stat.h>  // mkdir
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

void SavePartitionedMesh(const string &prefix, Mesh &mesh, int np,
                         int64_t *partitioning)
{
   MeshPartitioner partitioner(mesh, np, partitioning);
   MeshPart mesh_part;
   for (int i = 0; i < np; i++)
   {
      cout << i << endl;
      partitioner.ExtractPart(i, mesh_part);

      ofstream f(MakeParFilename(prefix + ".mesh.", i));
      f.precision(16);
      mesh_part.Print(f);
   }
}

int main(int argc, char *argv[])
{
   string mesh_file = "";
   string dir = "PartitionedMesh";
   int np = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&np, "-np", "--n-partitions", "Number of mesh partitions.");
   args.AddOption(&dir, "-d", "--dir", "Data directory.");
   args.ParseCheck();

   MFEM_VERIFY(!mesh_file.empty(), "Must provide input mesh file.");

   int err_flag = mkdir(dir.c_str(), 0777);
   err_flag = (err_flag && (errno != EEXIST)) ? 1 : 0;
   MFEM_VERIFY(err_flag == 0, "Could not create directory " << dir)

   cout << "\n";
   tic_toc.Restart();
   cout << "Reading fine mesh... " << flush;
   unique_ptr<Mesh> mesh(new Mesh(mesh_file));
   cout << "Done. " << tic_toc.RealTime() << endl;

   // Partition the fine mesh
   tic_toc.Restart();
   cout << "Partitioning... " << flush;
   Array<int64_t> partitioning(mesh->GeneratePartitioning(np), mesh->GetNE(),
                               true);
   cout << "Done. " << tic_toc.RealTime() << endl;

   cout << "Saving..." << endl;
   SavePartitionedMesh(dir + "/partitioned", *mesh, np, partitioning);
   cout << "Done." << endl;

   return 0;
}
