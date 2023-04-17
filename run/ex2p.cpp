//
#include <iostream>
#include <fstream>
#include <mfem.hpp>

// Solve
//
// -div(sigma(u)) = 0
// sigma(u) = lambda * div(u) * I + mu * (grad * u + u * grad)
//
// BC
// Attribute 1 - u = 0
// Attribute 2 - sigma(u) * n = f
// Rest:         sigma(u) * n = 0
//
//

using namespace std;
using namespace mfem;

int main(int argc, char* argv[])
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

    // read mesh
    const char *mesh_file = "../data/beam-tri.mesh";
    Mesh serial_mesh(mesh_file, 1, 1);
    int dim = serial_mesh.Dimension();

    if (serial_mesh.attributes.Max() < 2 ||
        serial_mesh.bdr_attributes.Max() < 2)
    {
        if (pid == 0)
        {
            cerr << "\nInput mesh should have at least two materials and "
                 << "two boundary attributes\n" << endl;
        }
        return 3;
    }

    {
        int ref_levels = 0;
        for (int i=0; i < ref_levels; i++)
        {
            serial_mesh.UniformRefinement();
        }
    }

    ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
    serial_mesh.Clear();

    {
        int par_ref_levels = 1;
        for (int l = 0; l < par_ref_levels; l++)
        {
            mesh.UniformRefinement();
        }
    }

    // create FEM space
    int order = 8;
    bool reorder_space = false;

    FiniteElementCollection *fec;
    ParFiniteElementSpace *fespace;

    fec = new H1_FECollection(order, dim);
    if (reorder_space)
    {
        fespace = new ParFiniteElementSpace(&mesh, fec, dim, Ordering::byNODES);
    }
    else
    {
        fespace = new ParFiniteElementSpace(&mesh, fec, dim, Ordering::byVDIM);
    }

    HYPRE_BigInt num_dofs = (*fespace).GlobalTrueVSize();

    if (pid == 0)
    {
        cout << "Number of Finite Element DOFs: " << num_dofs << endl;
    }

    // set up and apply boundary condition
    // apply BC at boundary attribute i where bc_dir[i] != 0
    Array<int> bc_dir(mesh.bdr_attributes.Max());
    bc_dir = 0;
    bc_dir[0] = 1;

    // get boundary degrees of freedom corresponding to boundary attributes
    // i where bc_dir[i] != 0
    Array<int> bc_dir_dofs;
    fespace->GetEssentialTrueDofs(bc_dir, bc_dir_dofs);

    // set up boundary forcing function

    VectorArrayCoefficient f(dim);
    for (int i=0; i<dim-1; i++)
    {
        // set X component to 0
        f.Set(i, new ConstantCoefficient(0.0));
    }
    {
        // define forcing function at each boundary attribute
        Vector pull_force(mesh.bdr_attributes.Max());
        pull_force[0] = 0.0;
        pull_force[1] = -1.0e-2;

        f.Set(dim-1, new PWConstCoefficient(pull_force));
    }

    // Right hand side form
    ParLinearForm b(fespace);
    b.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
    b.Assemble();

    // Solution function
    ParGridFunction x(fespace);
    x = 0.0;

    // Left hand side form
    bool partial_assembly = false;
    bool static_cond = false;

    if (pid == 0)
    {
        cout << "mesh attrbutes #: " << mesh.attributes.Max() << endl;
        cout << "bdry attrbutes #: " << mesh.bdr_attributes.Max() << endl;
    }

    // define lame coefficients at each boundary attribute
    Vector lambda (mesh.attributes.Max());
    lambda[1] = 1.0;
    lambda[0] = lambda(1) * 50;
    PWConstCoefficient lambda_func(lambda);

    Vector mu(mesh.attributes.Max());
    mu[1] = 1.0;
    mu[0] = mu(1) * 50;
    PWConstCoefficient mu_func(mu);

    ParBilinearForm a(fespace);
    a.AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

    if (static_cond)
    {
        a.EnableStaticCondensation();
    }
    a.Assemble();

    // Linear system
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(bc_dir_dofs, x, b, A, X, B);

    if (pid == 0){
        cout << "Assembed linear system of size: " << A.GetGlobalNumRows() << endl;
    }

    // Preconditioner
    bool algebraic_ceed = false;
    HypreBoomerAMG amg(A);
    amg.SetSystemsOptions(dim, reorder_space);
    //amg.SetElasticityOptions(fespace);

    // Solver
    HyprePCG pcg(A);
    pcg.SetTol(1e-8);
    pcg.SetMaxIter(500);
    pcg.SetPrintLevel(2);
    pcg.SetPreconditioner(amg);

    // Solve
    pcg.Mult(B, X);

    // recover solution
    a.RecoverFEMSolution(X, b, x);

    // save paraview file
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection("Ex2p", &mesh);
    pd->SetPrefixPath("results");
    pd->RegisterField("sol", &x);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    // cleanup
    //delete pcg;
    //delete amg;
    //delete a;
    //delete b;
    //delete fespace;
    //delete fec;
    mesh.Clear();

    return 0;
}
