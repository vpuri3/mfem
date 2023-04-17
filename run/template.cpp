//
#include <iostream>
#include <fstream>
#include <mfem.hpp>

// Solve

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
    Mesh *serial_mesh = new Mesh(mesh_file);
    int dim = serial_mesh->Dimension();

    ParMesh *mesh = new ParMesh(MPI_COMM_WORLD, *serial_mesh);
    serial_mesh->Clear();

    // create FEM space
    int order = 4;
    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(mesh, fec);

    // set up and apply boundary condition

    // Right hand side form

    // Solution function

    // Left hand side form
    bool partial_assembly = false;
    bool static_cond = false;

    // Linear system

    // Preconditioner
    bool algebraic_ceed = false;

    // Solver

    // Solve

    // recover solution

    // save paraview file
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection("Ex1p", mesh);
    pd->SetPrefixPath("results");
    pd->RegisterField("sol", &x);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    // cleanup

    return 0;
}
