
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   const char *mesh_file = "data/mesh_100k.mesh";
   // const char *mesh_file = "mesh_300k.mesh";

   Mpi::Init();
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // const char *mesh_file = "mesh.mesh";

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   Mesh mesh(mesh_file);


   char vishost[] = "localhost";
   int  visport   = 19916;
   socketstream mesh_sock(vishost, visport);
   mesh_sock.precision(8);
   mesh_sock << "mesh\n" << mesh << flush;


   std::filebuf fb;
   //   fb.open ("sol.gf",std::ios::in);
   fb.open ("data/eps_r_100k.gf",std::ios::in);
   std::istream is(&fb);
   GridFunction * gf = new GridFunction(&mesh, is);
   fb.close();

   FiniteElementSpace * vfes = gf->FESpace();
   int vdim = vfes->GetVDim();
   const FiniteElementCollection * fec = vfes->FEColl();
   FiniteElementSpace * fes = new FiniteElementSpace(&mesh, fec);
   Array<GridFunction *> gfs(vdim);
   Array<ParGridFunction *> pgfs(vdim);
   Array<socketstream *> ss(vdim);
   Array<ostringstream *> oss(vdim);
   for (int i = 0; i<vdim; i++)
   {
      oss[i] = new ostringstream;
   }
   *oss[0] << " ϵ_xx ";
   *oss[1] << " ϵ_xy";
   *oss[2] << " ϵ_xz ";
   *oss[3] << " ϵ_yx ";
   *oss[4] << " ϵ_yy ";
   *oss[5] << " ϵ_yz ";
   *oss[6] << " ϵ_zx ";
   *oss[7] << " ϵ_zy ";
   *oss[8] << " ϵ_zz ";
   double *data = gf->GetData();

   string keys;
   keys = "keys iXXXXXXXXXXXXXXXXXac\n";

   for (int i = vdim-1; i>-1; i--)
   {
      gfs[i] = new GridFunction(fes,&data[i*fes->GetVSize()]);
      // ss[i] = new socketstream(vishost,visport);
      // ss[i]->precision(8);
      // *ss[i] << "solution\n" << mesh << *gfs[i] << keys
      //  << "window_title '" << oss[i]->str() << "'" << flush;

   }


   int * partitioning = mesh.GeneratePartitioning(num_procs);
   ParMesh pmesh(MPI_COMM_WORLD, mesh);

   for (int i = vdim-1; i>-1; i--)
   {
      pgfs[i] = new ParGridFunction(&pmesh,gfs[i],partitioning);
      ss[i] = new socketstream(vishost,visport);
      ss[i]->precision(8);
      *ss[i] << "parallel " << num_procs << " " << myid << "\n";
      *ss[i] << "solution\n" << pmesh << *pgfs[i] << keys
             << "window_title '" << oss[i]->str() << "'" << flush;

   }





   // socketstream sol_sock(vishost, visport);
   // sol_sock.precision(8);
   // sol_sock << "solution\n" << mesh << *gfs[0] << flush;

   return 0;
}
