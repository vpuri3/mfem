//
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // init MPI, HYPRE
    Mpi::Init();
    int nprocs = Mpi::WorldSize();
    int pid = Mpi::WorldRank();
    Hypre::Init();

    const char *device_config = "cpu";
    Device device(device_config);
    if (pid == 0)
    {
        device.Print();
    }

    // Read in mesh
    const char *mesh_file = "../data/star.mesh";
    Mesh serial_mesh(mesh_file, 1, 1);
    int dim = serial_mesh.Dimension();

    {
        int ref_levels = 2;
        for (int l = 0; l < ref_levels; l++)
        {
            serial_mesh.UniformRefinement();
        }
    }

    ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
    serial_mesh.Clear();

    // Create FEM space
    int order = 7;
    H1_FECollection fec(order, dim);
    ParFiniteElementSpace fespace(&mesh, &fec);

    HYPRE_BigInt num_dofs = fespace.GlobalTrueVSize();

    if (pid == 0)
    {
        cout << "Number of degrees of freedom: " << num_dofs << endl;
    }

    // set up and apply boundary condition
    Array<int> bc_ess_dofs;
    if (mesh.bdr_attributes.Size())
    {
        Array<int> bdr_ess(mesh.bdr_attributes.Max());
        bdr_ess = 1;
        fespace.GetEssentialTrueDofs(bdr_ess, bc_ess_dofs);
    }

    // Right hand side form
    ParLinearForm b(&fespace);
    ConstantCoefficient one(1.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    // Solution function
    ParGridFunction x(&fespace);
    x = 0.0;

    // Left hand side form
    bool partial_assembly = true;
    bool static_cond = false;

    ParBilinearForm a(&fespace);
    if (partial_assembly)
    {
        a.SetAssemblyLevel(AssemblyLevel::PARTIAL);
    } else
    {
        a.SetAssemblyLevel(AssemblyLevel::FULL);
    }

    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    if (static_cond)
    {
        a.EnableStaticCondensation();
    }

    a.Assemble();

    // linear system
    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(bc_ess_dofs, x, b, A, X, B);

    // Preconditioner
    bool algebraic_ceed = false;

    Solver *prec = NULL;

    if (partial_assembly)
    {
        if (UsesTensorBasis(fespace))
        {
            if (algebraic_ceed)
            {
                prec = new ceed::AlgebraicSolver(a, bc_ess_dofs);
            } else
            {
                prec = new OperatorJacobiSmoother(a, bc_ess_dofs);
            }
        }
    } else
    {
        prec = new HypreBoomerAMG;
    }

    // Solver
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    cg.SetPreconditioner(*prec);
    cg.SetOperator(*A);

    // Solve
    cg.Mult(B, X);

    // recover solution
    a.RecoverFEMSolution(X, b, x);

    // save paraview file
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection("Ex1p", &mesh);
    pd->SetPrefixPath("results");
    pd->RegisterField("sol", &x);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    // end
    delete prec;

    return 0;
}
