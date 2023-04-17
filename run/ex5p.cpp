//
#include <iostream>
#include <fstream>
#include <mfem.hpp>

// Solve 2D mixed Darcy problem
//
// k * u + grad p = f
// - div u        = g
//
// with prescribed natural (dirichlet) BC: -p = <given pressure>.
// Use RT space for velocity, L2 space for pressure. We also use
// BlockOperator for solving
//

using namespace std;
using namespace mfem;

void u_exact(const Vector & x, Vector & u);
double pres_exact(const Vector & x);
void forc(const Vector & x, Vector & f);
double divr(const Vector & x);
double f_natr(const Vector & x);

int main(int argc, char *argv[])
{
    StopWatch watch;

    // init MPI, Hypre
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

    // mesh
    const char *mesh_file = "../data/star.mesh";
    Mesh serial_mesh(mesh_file, 1, 1);
    int dim = serial_mesh.Dimension();

    ParMesh *mesh = new ParMesh(MPI_COMM_WORLD, serial_mesh);
    serial_mesh.Clear();

    int ref_levels = 1;
    for (int i=0; i<ref_levels;i++)
    {
        mesh->UniformRefinement();
    }

    // create FE space
    int order = 8;
    // TODO - try with H1/L2 (Taylor-Hood) elements?
    FiniteElementCollection *hdiv_coll = new RT_FECollection(order, dim); //vector
    FiniteElementCollection *l2_coll = new L2_FECollection(order, dim);

    ParFiniteElementSpace *R_space = new ParFiniteElementSpace(mesh, hdiv_coll);
    ParFiniteElementSpace *W_space = new ParFiniteElementSpace(mesh, l2_coll);

    HYPRE_BigInt dimR = R_space->GlobalTrueVSize();
    HYPRE_BigInt dimW = W_space->GlobalTrueVSize();

    if (pid == 0)
    {
        cout << "#========================================#" << endl;
        cout << "dim: "    << dim  << endl;
        cout << "dim(R): " << dimR << endl;
        cout << "dim(W): " << dimW << endl;
        cout << "#========================================#" << endl;
    }

    // define coefficients, analytical solution, and RHS of PDE

    // k * u + grad p = f
    // - div u        = g

    ConstantCoefficient k(1.0);
    VectorFunctionCoefficient f(dim, forc);
    FunctionCoefficient f_n(f_natr);
    FunctionCoefficient g(divr);

    // Define problem block structure
    // computed block offsets locally on each core

    // used for Vector based on dof
    // eg ParGridFunction, ParLinearForm
    Array<int> block_offsets(3); // # variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = R_space->GetVSize();
    block_offsets[2] = W_space->GetVSize();
    block_offsets.PartialSum(); // fill array with cumilative sum of entries

    // used for vector based on trueDof
    // eg HypreParVector for RHS, sol. to linear system
    Array<int> block_trueOffsets(3); // # variables + 1
    block_trueOffsets[0] = 0;
    block_trueOffsets[1] = R_space->TrueVSize();
    block_trueOffsets[2] = W_space->TrueVSize();
    block_trueOffsets.PartialSum();

    // form right-hand-side
    MemoryType mt = device.GetMemoryType();
    BlockVector x(block_offsets, mt), rhs(block_offsets, mt);
    BlockVector trueX(block_trueOffsets, mt), trueRhs(block_trueOffsets, mt);

    ParLinearForm *fform = new ParLinearForm();
    fform->Update(R_space, rhs.GetBlock(0), 0);
    fform->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
    fform->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(f_n));
    fform->Assemble();
    fform->SyncAliasMemory(rhs);
    fform->ParallelAssemble(trueRhs.GetBlock(0));
    trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);

    ParLinearForm *gform = new ParLinearForm();
    gform->Update(W_space, rhs.GetBlock(1), 0);
    gform->AddDomainIntegrator(new DomainLFIntegrator(g));
    gform->Assemble();
    gform->SyncAliasMemory(rhs);
    gform->ParallelAssemble(trueRhs.GetBlock(1));
    trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);

    // assemble left-hand-side
    //
    // Darcy Operator D = [ M  B^T ]
    //                    [ B  0   ]
    //
    // where
    //
    // M = \int_\Omega k u      \cdot v dx, u, v \in RT_space
    // B = \int_\Omega \div( u) \cdot q dx, u \in RT_space, q \in W_space

    bool partial_assembly = false;

    partial_assembly = true; // TODO - fix segfault

    ParBilinearForm *Mform = new ParBilinearForm(R_space);
    ParMixedBilinearForm *Bform = new ParMixedBilinearForm(R_space, W_space);

    if (partial_assembly)
    {
        Mform->SetAssemblyLevel(AssemblyLevel::PARTIAL);
        Bform->SetAssemblyLevel(AssemblyLevel::PARTIAL);
    }

    Mform->AddDomainIntegrator(new VectorFEMassIntegrator(k));
    Mform->Assemble(); // local (to element) assembly

    Bform->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
    Bform->Assemble(); // local (to assemble) assembly

    if (!partial_assembly)
    {
        Mform->Finalize();
        Bform->Finalize();
    }

    // assemble darcy op
    BlockOperator *darcyOp = new BlockOperator(block_trueOffsets);

    TransposeOperator *Bt = NULL;
    OperatorPtr opM, opB;  // for partial asssembly
    HypreParMatrix *M, *B; // for full assembly

    OperatorPtr opI;

    Array<int> empty_tdof_list; // empty list of DOFs for essential BC

    if (partial_assembly)
    {
        Mform->FormSystemMatrix(empty_tdof_list, opM);
        Bform->FormRectangularSystemMatrix(empty_tdof_list, empty_tdof_list, opB);
        Bt = new TransposeOperator(opB.Ptr());

        //OperatorPtr opI = new IdentityOperator(R_space->GetTrueVSize());

        darcyOp->SetBlock(0, 0, opM.Ptr());
        darcyOp->SetBlock(0, 1, Bt, -1.0);
        darcyOp->SetBlock(1, 0, opB.Ptr(), -1.0);
    }
    else
    {
        M = Mform->ParallelAssemble(); // eleminate inter-element dofs
        B = Bform->ParallelAssemble();
        (*B) *= -1;
        Bt = new TransposeOperator(B);

        darcyOp->SetBlock(0, 0, M);
        darcyOp->SetBlock(0, 1, Bt);
        darcyOp->SetBlock(1, 0, B);
    }

    // Preconditioner - TODO
    //
    // P = [ diag(M)           0       ]
    //     [    0     B diag(M)^-1 B^T ]
    //
    // use symmetric Gauss-Seidel to approximate the inverse of the
    // pressure Schur Complement

    /*
    HypreParMatrix *MinvBt = NULL;
    HypreParVector *Mdiag = NULL;
    HypreParMatrix *S = NULL;

    Solver *invM, *invS;

    if (partial_assembly)
    {
        Vector Md(R_space->GetTrueVSize()); // why is this local?
        Vector Mdi(R_space->GetTrueVSize());

        Mform->AssembleDiagonal(Md);

        invM = new OperatorJacobiSmoother(_, empty_tdof_list);
        invS = new OperatorJacobiSmoother(_, empty_tdof_list);
    }
    else
    {
        Mdiag = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),
                                   M->GetRowStarts());
        M->GetDiag(*Mdiag);

        MinvBt;
        S = ParMult(B, MinvBt);

        invM->iterative_mode = false;
        invS->iterative_mode = false;
    }

    BlockDiagonalPreconditioner *darcyPr = new BlockDiagonalPreconditioner(
            block_trueOffsets);

    darcyPr->SetDiagonalBlock(0, invM);
    darcyPr->SetDiagonalBlock(1, invS);
    */

    // Solver
    MINRESSolver solver(MPI_COMM_WORLD);
    solver.SetAbsTol(1e-10);
    solver.SetRelTol(1e-6);
    solver.SetMaxIter(5000);
    solver.SetOperator(*darcyOp);
    //solver.SetPreconditioner(*darcyPr);
    solver.SetPrintLevel(2);

    // Solve
    watch.Clear();
    watch.Start();

    solver.Mult(trueRhs, trueX);

    watch.Stop();

    if (pid == 0)
    {
        if (solver.GetConverged())
        {
            cout << "Solver converged in " << solver.GetNumIterations()
                 << " iterations with residual norm " << solver.GetFinalNorm() << endl;
        }
        else
        {
            cout << "Soler did not converge. Final residual norm "
                 << solver.GetFinalNorm() << endl;
        }

        cout << "Time taken: " << watch.RealTime() << endl;
    }

    // Analysis
    ParGridFunction *u = new ParGridFunction(R_space, x.GetBlock(0), 0);
    ParGridFunction *p = new ParGridFunction(W_space, x.GetBlock(1), 0);

    u->Distribute(&(trueX.GetBlock(0)));
    p->Distribute(&(trueX.GetBlock(1)));

    int quadrature_order = max(2, order * 2 + 1);
    const IntegrationRule *int_rules[Geometry::NumGeom];
    for (int i=0;i<Geometry::NumGeom;i++)
    {
        int_rules[i] = &(IntRules.Get(i, quadrature_order));
    }

    if(pid == 0)
    {
        cout << "Geometry::NumGeom " << Geometry::NumGeom << endl;
    }

    VectorFunctionCoefficient u_ex(dim, u_exact);
    FunctionCoefficient p_ex(pres_exact);

    double err_u = u->ComputeL2Error(u_ex, int_rules);
    double norm_u = ComputeGlobalLpNorm(2, u_ex, *mesh, int_rules);

    double err_p = p->ComputeL2Error(p_ex, int_rules);
    double norm_p = ComputeGlobalLpNorm(2, p_ex, *mesh, int_rules);

    if (pid == 0)
    {
        cout << "|| u_h - u_ex || / || u_ex || = " << err_u / norm_u << "\n";
        cout << "|| p_h - p_ex || / || p_ex || = " << err_p / norm_p << "\n";
    }

    // save paraview file
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection("Ex5p", mesh);
    pd->SetPrefixPath("results");
    pd->RegisterField("velocity", u);
    pd->RegisterField("pressure", p);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    // delete stuff

    return 0;
}

// TODO - do these double allocations count?? can we just use value by reference?

void u_exact(const Vector & x, Vector & u)
{
    double _x(x(0));
    double _y(x(1));
    double _z(0.0);
    if (x.Size() == 3)
    {
        _z = x(2);
    }

    u(0) = -exp(_x) * sin(_y) * cos(_z);
    u(1) = -exp(_x) * cos(_y) * cos(_z);

    if (x.Size() == 3)
    {
        u(2) = exp(_x) * sin(_y) * sin(_z);
    }

    return;
}

double pres_exact(const Vector & x)
{
    double _x(x(0));
    double _y(x(1));
    double _z(0.0);
    if (x.Size() == 3)
    {
        _z = x(2);
    }

    return exp(_x) * sin(_y) * cos(_z);
}

void forc(const Vector & x, Vector & f)
{
    f = 0.0;
}

double divr(const Vector & x)
{
    if (x.Size() == 3)
    {
        return - pres_exact(x);
    }
    else
    {
        return 0.0;
    }
}

double f_natr(const Vector & x)
{
    return - pres_exact(x);
}

