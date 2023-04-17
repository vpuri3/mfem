//
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // init MPI, Hypre
    Mpi::Init(argc, argv);
    Hypre::Init();

    // Read in mesh
    const char *mesh_file = "../data/square-disc.mesh";
    Mesh serial_mesh(mesh_file);

    ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
    serial_mesh.Clear();

    mesh.UniformRefinement();

    // Create FEM space
    int order = 1;
    H1_FECollection fec(order, mesh.Dimension());
    ParFiniteElementSpace fespace(&mesh, &fec);

    HYPRE_BigInt num_dofs = fespace.GlobalTrueVSize();

    if (Mpi::Root())
    {
        cout << "Number of degrees of freedom: " << num_dofs << endl;
    }

    // Get boundary nodes
    Array<int> boundary_dofs;
    fespace.GetBoundaryTrueDofs(boundary_dofs);

    ParGridFunction x(&fespace);

    // Assign dirichlet BC at boundary points
    x = 0;

    // Right hand side
    ConstantCoefficient one(1.0);
    ParLinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    // Left hand side
    ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator);
    a.Assemble();

    // Linear System
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(boundary_dofs, x, b, A, X, B);

    // preconditioner
    HypreBoomerAMG M(A);

    if (Mpi::Root())
    {
        cout << __LINE__ << endl;
    }
    // preconditioner
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    cg.SetPreconditioner(M);
    cg.SetOperator(A);

    // solve
    cg.Mult(B, X);

    // recover solution as grid function
    a.RecoverFEMSolution(X, b, x);
    //x.Save("sol");
    //mesh.Save("mesh");

    // save paraview file
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection("Ex0p", &mesh);
    pd->SetPrefixPath("results");
    pd->RegisterField("sol", &x);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    delete pd;

    return 0;
}
