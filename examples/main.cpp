#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

inline void ApplyBlkMult(const mfem::DenseTensor &Mat,
                         const mfem::Vector &x,
                         mfem::Vector &y)
{

   const int ndof = Mat.SizeI();
   const int NE = Mat.SizeK();

   auto X  = mfem::Reshape(x.Read(), ndof, NE);
   auto Y  = mfem::Reshape(y.Write(), ndof, NE);
   auto Me = mfem::Reshape(Mat.Read(), ndof, ndof, NE);

   MFEM_FORALL(tid, ndof * NE,
   {
      const int c = tid % ndof;
      const int e = tid / ndof;

      //for (int c=0; c<ndof; ++c)
      {
         double dot = 0;
         for (int r = 0; r < ndof; ++r)
         {
            dot += Me(r, c, e) * X(r, e);
         }
         Y(c, e) = dot;
      }
   });
}

int main(int argc, char *argv[])
{

   std::cout<<"mfem cuda/hip block solver "<<std::endl;

   const char *mesh_file = "../data/periodic-hexagon.mesh";
   int order = 2;
   const char *device_config = "cuda";

   Device device(device_config);
   device.Print();

   // 2. Read the mesh from the given mesh file. We can handle geometrically
   //    periodic meshes in this code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 5. Define the discontinuous DG finite element space of the given
   //    polynomial order on the refined mesh.
   DG_FECollection fec(order, dim, BasisType::GaussLobatto);
   FiniteElementSpace fes(&mesh, &fec);

   cout << "Number of unknowns: " << fes.GetVSize() << endl;

   const int NE     = fes.GetMesh()->GetNE();
   const int ndofs  = fes.GetNDofs()/NE;

   std::cout<<"number of dofs per elem = "<<ndofs<<std::endl;

   MassIntegrator *mass_Int = new MassIntegrator;

   Vector elemMass(ndofs * ndofs * NE);

   //Populate element mass matrix
   mass_Int->AssembleEA(fes, const_cast<mfem::Vector &>(elemMass), false);

   Vector x(ndofs * NE);
   Vector x_ref(ndofs * NE);
   Vector y(ndofs * NE); y.Randomize();

   x_ref = y;
   x = y;

   DenseTensor massMats(ndofs, ndofs, NE);

   const double *h_elemMass = elemMass.HostRead();
   for (int i=0; i<elemMass.Size(); ++i)
   {
      massMats.HostWrite()[i] = h_elemMass[i];
   }

   LibBatchMult libMult(massMats);

   libMult.Mult(y, x);

   ApplyBlkMult(massMats, y, x_ref);

   x_ref -= x;

   if (x_ref.Norml2() > 1e-10)
   {
      std::cout<<"mult error too high! error = "<< x_ref.Norml2()<<std::endl;
   }
   else
   {
      std::cout<<"mult solutions agree! "<<std::endl;
   }



   //Reset
   x_ref = y;
   x = y;

   Array<int> P(ndofs*ndofs*NE);

   DenseTensor BatchLU = massMats;
   BatchLUFactor(BatchLU, P);
   BatchLUSolve(BatchLU, P, x_ref);

   LibBatchSolve libSolver(massMats, LibBatchMode::INVERSE);

   libSolver.Setup();

   libSolver.Mult(y, x);

   //x.Print(mfem::out,9);
   /*
   x_ref.Print();
   std::cout<<"----"<<std::endl;
   */

   x_ref -= x;


   double error = x_ref.Norml2();
   if (error > 1e-10)
   {
      std::cout<<"error too high! error = "<< error<<std::endl;
   }
   else
   {
      std::cout<<"solutions agree! "<<std::endl;
   }

   //x.Print();
   //std::cout<<"-----"<<std::endl;
   //x_ref.Print();


   delete mass_Int;
}
