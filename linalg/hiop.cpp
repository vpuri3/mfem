// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include "../config/config.hpp"
#include "hiop.hpp"


#ifdef MFEM_USE_HIOP
#include <iostream>

#include "hiopAlgFilterIPM.hpp"

using namespace hiop;

namespace mfem
{

bool HiopOptimizationProblem::get_prob_sizes(long long &n, long long &m)
{
   n = ntdofs_glob;
   m = problem.GetNumConstraints();

   return true;
}

bool HiopOptimizationProblem::get_starting_point(const long long &n, double *x0)
{
   MFEM_ASSERT(x_start != NULL && ntdofs_loc == x_start->Size(),
               "Starting point is not set properly.");

   memcpy(x0, x_start->GetData(), ntdofs_loc * sizeof(double));

   return true;
}

bool HiopOptimizationProblem::get_vars_info(const long long &n,
                                            double *xlow, double *xupp,
                                            NonlinearityType *type)
{
   MFEM_ASSERT(n == ntdofs_glob, "Global input mismatch.");
   MFEM_ASSERT(problem.x_lo && problem.x_hi,
               "Solution bounds are not set!");

   std::memcpy(xlow, problem.x_lo->GetData(), ntdofs_loc * sizeof(double));
   std::memcpy(xupp, problem.x_hi->GetData(), ntdofs_loc * sizeof(double));

   return true;
}

bool HiopOptimizationProblem::get_cons_info(const long long &m,
                                            double *clow, double *cupp,
                                            NonlinearityType *type)
{
   MFEM_ASSERT(m == m_total, "Global constraint size mismatch.");

   int csize = 0;
   if (problem.C)
   {
      csize = problem.c_e->Size();
      std::memcpy(clow, problem.c_e->GetData(), csize * sizeof(double));
      std::memcpy(cupp, problem.c_e->GetData(), csize * sizeof(double));
   }
   if (problem.D)
   {
      const int dsize = problem.d_lo->Size();
      std::memcpy(clow + csize, problem.d_lo->GetData(), dsize *sizeof(double));
      std::memcpy(cupp + csize, problem.d_hi->GetData(), dsize *sizeof(double));
   }

   return true;
}

bool HiopOptimizationProblem::eval_f(const long long &n, const double *x,
                                     bool new_x, double &obj_value)
{
   MFEM_ASSERT(n == ntdofs_glob, "Global input mismatch.");

   if (new_x) { constr_info_is_current = false; }

   Vector x_vec(ntdofs_loc);
   x_vec = x;
   obj_value = problem.CalcObjective(x_vec);

   return true;
}

bool HiopOptimizationProblem::eval_grad_f(const long long &n, const double *x,
                                          bool new_x, double *gradf)
{
   MFEM_ASSERT(n == ntdofs_glob, "Global input mismatch.");

   if (new_x) { constr_info_is_current = false; }

   Vector x_vec(ntdofs_loc), gradf_vec(ntdofs_loc);
   x_vec = x;
   problem.CalcObjectiveGrad(x_vec, gradf_vec);
   std::memcpy(gradf, gradf_vec.GetData(), ntdofs_loc * sizeof(double));

   return true;
}

bool HiopOptimizationProblem::eval_cons(const long long &n, const long long &m,
                                        const long long &num_cons,
                                        const long long *idx_cons,
                                        const double *x, bool new_x,
                                        double *cons)
{
   MFEM_ASSERT(n == ntdofs_glob, "Global input mismatch.");
   MFEM_ASSERT(m == m_total, "Constraint size mismatch.");
   MFEM_ASSERT(num_cons <= m, "num_cons should be at most m = " << m);

   if (num_cons == 0) { return true; }

   if (new_x) { constr_info_is_current = false; }
   Vector x_vec(ntdofs_loc);
   x_vec = x;
   UpdateConstrValsGrads(x_vec);

   for (int c = 0; c < num_cons; c++)
   {
      MFEM_ASSERT(idx_cons[c] < m_total, "Constraint index is out of bounds.");
      cons[c] = constr_vals(idx_cons[c]);
   }

   return true;
}

bool HiopOptimizationProblem::eval_Jac_cons(const long long &n,
                                            const long long &m,
                                            const long long &num_cons,
                                            const long long *idx_cons,
                                            const double *x, bool new_x,
                                            double **Jac)
{
   MFEM_ASSERT(n == ntdofs_glob, "Global input mismatch.");
   MFEM_ASSERT(m == m_total, "Constraint size mismatch.");
   MFEM_ASSERT(num_cons <= m, "num_cons should be at most m = " << m);

   if (num_cons == 0) { return true; }

   if (new_x) { constr_info_is_current = false; }
   Vector x_vec(ntdofs_loc);
   x_vec = x;
   UpdateConstrValsGrads(x_vec);

   for (int c = 0; c < num_cons; c++)
   {
      MFEM_ASSERT(idx_cons[c] < m_total, "Constraint index is out of bounds.");
      for (int j = 0; j < ntdofs_loc; j++)
      {
         Jac[c][j] = constr_grads(idx_cons[c], j);
      }
   }

   return true;
}

bool HiopOptimizationProblem::get_vecdistrib_info(long long global_n,
                                                  long long *cols)
{
#ifdef MFEM_USE_MPI
   int nranks;
   MPI_Comm_size(comm_, &nranks);

   long long *sizes = new long long[nranks];
   MPI_Allgather(&ntdofs_loc, 1, MPI_LONG_LONG_INT, sizes, 1,
                 MPI_LONG_LONG_INT, comm_);
   cols[0] = 0;
   for (int r = 1; r <= nranks; r++)
   {
      cols[r] = sizes[r-1] + cols[r-1];
   }

   delete [] sizes;
   return true;
#else
   // Returning false means that Hiop runs in non-distributed mode.
   return false;
#endif
}

void HiopOptimizationProblem::UpdateConstrValsGrads(const Vector x)
{
   if (constr_info_is_current) { return; }

   // If needed (e.g. for CG spaces), communication should be handled by the
   // operators' Mult() and GetGradient() methods.

   int cheight = 0;
   if (problem.C)
   {
      cheight = problem.C->Height();

      // Values of C.
      Vector vals_C(constr_vals.GetData(), cheight);
      problem.C->Mult(x, vals_C);

      // Gradients C.
      const Operator &oper_C = problem.C->GetGradient(x);
      const DenseMatrix *grad_C = dynamic_cast<const DenseMatrix *>(&oper_C);
      MFEM_VERIFY(grad_C, "Hiop expects DenseMatrices as operator gradients.");
      MFEM_ASSERT(grad_C->Height() == cheight && grad_C->Width() == ntdofs_loc,
                  "Incorrect dimensions of the C constraint gradient.");
      for (int i = 0; i < cheight; i++)
      {
         for (int j = 0; j < ntdofs_loc; j++)
         {
            constr_grads(i, j) = (*grad_C)(i, j);
         }
      }
   }

   if (problem.D)
   {
      const int dheight = problem.D->Height();

      // Values of D.
      Vector vals_D(constr_vals.GetData() + cheight, dheight);
      problem.D->Mult(x, vals_D);

      // Gradients of D.
      const Operator &oper_D = problem.D->GetGradient(x);
      const DenseMatrix *grad_D = dynamic_cast<const DenseMatrix *>(&oper_D);
      MFEM_VERIFY(grad_D, "Hiop expects DenseMatrices as operator gradients.");
      MFEM_ASSERT(grad_D->Height() == dheight && grad_D->Width() == ntdofs_loc,
                  "Incorrect dimensions of the D constraint gradient.");
      for (int i = 0; i < dheight; i++)
      {
         for (int j = 0; j < ntdofs_loc; j++)
         {
            constr_grads(i + cheight, j) = (*grad_D)(i, j);
         }
      }
   }

   constr_info_is_current = true;
}

HiopNlpOptimizer::HiopNlpOptimizer() : OptimizationSolver(), hiop_problem(NULL)
{ 
#ifdef MFEM_USE_MPI
  //in case a serial driver in parallel MFEM build calls HiOp
  comm_ = MPI_COMM_WORLD;
  int initialized, nret = MPI_Initialized(&initialized); 
  MFEM_ASSERT(MPI_SUCCESS==nret, "failure in calling MPI_Initialized");
  if(!initialized) {
    nret = MPI_Init(NULL, NULL);
    MFEM_ASSERT(MPI_SUCCESS==nret, "failure in calling MPI_Init");
  }
#endif
} 

#ifdef MFEM_USE_MPI
HiopNlpOptimizer::HiopNlpOptimizer(MPI_Comm _comm) 
  : OptimizationSolver(_comm),
    hiop_problem(NULL),
    comm_(_comm) { }
#endif

HiopNlpOptimizer::~HiopNlpOptimizer()
{
   delete hiop_problem;
}

void HiopNlpOptimizer::SetOptimizationProblem(const OptimizationProblem &prob)
{
   problem = &prob;
   height = width = problem->input_size;

   if (hiop_problem) { delete hiop_problem; }

#ifdef MFEM_USE_MPI
   hiop_problem = new HiopOptimizationProblem(comm_, *problem);
#else
   hiop_problem = new HiopOptimizationProblem(*problem);
#endif
}

void HiopNlpOptimizer::Mult(const Vector &xt, Vector &x) const
{
   MFEM_ASSERT(hiop_problem != NULL,
               "Unspecified OptimizationProblem that must be solved.");

   hiop_problem->setStartingPoint(xt);

   hiop::hiopNlpDenseConstraints hiopInstance(*hiop_problem);

   hiopInstance.options->SetNumericValue("rel_tolerance", rel_tol);
   hiopInstance.options->SetNumericValue("tolerance", abs_tol);
   hiopInstance.options->SetIntegerValue("max_iter", max_iter);

   hiopInstance.options->SetStringValue("fixed_var", "relax");
   hiopInstance.options->SetNumericValue("fixed_var_tolerance", 1e-20);
   hiopInstance.options->SetNumericValue("fixed_var_perturb", 1e-9);

   // 0: no output; 3: not too much
   hiopInstance.options->SetIntegerValue("verbosity_level", print_level);

   // Use the IPM solver.
   hiop::hiopAlgFilterIPM solver(&hiopInstance);
   const hiop::hiopSolveStatus status = solver.run();
   final_norm = solver.getObjective();
   final_iter = solver.getNumIterations();

   if (status != hiop::Solve_Success && status != hiop::Solve_Success_RelTol)
   {
      converged = false;
      MFEM_WARNING("HIOP returned with a non-success status: " << status);
   }
   else { converged = true; }

   // Copy the final solution in x.
   solver.getSolution(x.GetData());
}

} // mfem namespace
#endif // MFEM_USE_HIOP
