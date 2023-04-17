//
#include "config/_config.hpp"
#include "fem/pbilinearform.hpp"
#include "linalg/hypre.hpp"
#include "linalg/operator.hpp"
#include <iostream>
#include <fstream>
#include <mfem.hpp>

// Solve
//
// du/dt + v \cdot grad(u) = 0
//
// with periodic boundary condition, and given initial condition.
//
// We use a discontinuous galerkin method (FaceIntegrators),
// prepackaged time-integrators (explicit, implicit)

using namespace std;
using namespace mfem;

// choice for prolbme setup (velocit, IC, inflow bc)
int problem;

void vel_func(const Vector & x, Vector &v);
double u0_func(const Vector &x);
double inflow_func(const Vector &x);

// mesh bounding box
Vector bb_min, bb_max;

enum class PrecType : int
{
    ILU = 0,
    AIR = 1
};

class DGSolver : public Solver
{
    private:
        HypreParMatrix &M, &K;
        HypreParMatrix *A; // A = M - dt * K
        SparseMatrix Mdiag;
        GMRESSolver linsolve;
        Solver *prec;
        double dt;
    public:
        DGSolver(HypreParMatrix &M_,
            HypreParMatrix &K_,
            const FiniteElementSpace &fespace,
            PrecType prec_type) :

            M(M_), K(K_), A(NULL), linsolve(M.GetComm()), dt(-1.0)
        {
            int block_size = fespace.GetFE(0)->GetDof();
            if (prec_type == PrecType::ILU)
            {
                prec = new BlockILU(block_size,
                    BlockILU::Reordering::MINIMUM_DISCARDED_FILL);
            }
            else if (prec_type == PrecType::AIR)
            {
                MFEM_ABORT("need MFEM_HYPRE_VERSION >= 21800 for AIR prec.");
            }

            linsolve.iterative_mode = false;
            linsolve.SetRelTol(1e-9);
            linsolve.SetAbsTol(0.0);
            linsolve.SetMaxIter(100);
            linsolve.SetPrintLevel(0);
            linsolve.SetPreconditioner(*prec);

            M.GetDiag(Mdiag);
        }

        void setDT(double dt_)
        {
            if (dt_ != dt)
            {
                dt = dt_;
                delete A;

                A = Add(-dt, K, 0.0, K); // A = -dt * K
                SparseMatrix Adiag;
                A->GetDiag(Adiag);
                Adiag.Add(1.0, Mdiag);

                // this will also call SetOperator on preconditioner
                linsolve.SetOperator(*A);
            }
        }

        void SetOperator(const Operator &op)
        {
            linsolve.SetOperator(op);
        }

        virtual void Mult(const Vector &x, Vector &y) const
        {
            linsolve.Mult(x, y);
        }

        ~DGSolver()
        {
            delete prec;
            delete A;
        }
};

// Time dependent operator for ODE RHS
//
// du/dt = - v * grad(u)
//
// equivalent to
//
// M du/dt = K u + b
//
// This can be written as a general ODE
//
// du/dt = M^{-1} * (K * u + b)
//
class Dudt : public TimeDependentOperator
{
    private:
        OperatorHandle M, K;
        const Vector &b;
        Solver *precM;
        CGSolver solM;
        DGSolver solDG;

        mutable Vector z;

    public:
        Dudt(ParBilinearForm &M_, ParBilinearForm &K_,)
}


int main(int argc, char* argv[])
{
    return 0;
}

/*
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
    const char *mesh_file = "../data/periodic-square.mesh";
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
    pd = new ParaViewDataCollection("Ex9p", mesh);
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
*/
