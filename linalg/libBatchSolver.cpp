// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "libBatchSolver.hpp"

namespace mfem
{

void LibBatchSolve::AssignMatrices(DenseTensor &MatrixBatch_)
{
   LUMatrixBatch = MatrixBatch_;
}

void LibBatchSolve::ComputeLU()
{

   MFEM_SUB_cu_or_hip(blasStatus_t)
   status = MFEM_SUB_cu_or_hip(blasDgetrfBatched)(MFEM_SUB_Cuda_or_Hip(
                                                     BLAS::Handle()),
                                                  mat_size,
                                                  lu_ptr_array.ReadWrite(),
                                                  mat_size,
                                                  P.Write(),
                                                  info_array.Write(),
                                                  num_mats);

   MFEM_VERIFY(status == MFEM_SUB_CU_or_HIP(BLAS_STATUS_SUCCESS),
               "Failed at blasDgetrfBatched");
}

void LibBatchSolve::ComputeInverse()
{

   ComputeLU();

   inv_ptr_array.SetSize(num_mats);

   InvMatrixBatch.SetSize(LUMatrixBatch.SizeI(), LUMatrixBatch.SizeJ(),
                          LUMatrixBatch.SizeK());

   double *inv_ptr_base = InvMatrixBatch.Write();

   for (int i=0; i<num_mats; ++i)
   {
      inv_ptr_array.HostWrite()[i] = inv_ptr_base + i * mat_size * mat_size;
   }

   //Invert matrices
   MFEM_SUB_cu_or_hip(blasStatus_t)
   status = MFEM_SUB_cu_or_hip(blasDgetriBatched)(MFEM_SUB_Cuda_or_Hip(
                                                     BLAS::Handle)(),
                                                  mat_size,
                                                  lu_ptr_array.ReadWrite(),
                                                  mat_size,
                                                  P.ReadWrite(),
                                                  inv_ptr_array.ReadWrite(),
                                                  mat_size,
                                                  info_array.Write(),
                                                  num_mats);

   MFEM_VERIFY(status == MFEM_SUB_CU_or_HIP(BLAS_STATUS_SUCCESS),
               "Failed at blasDgetriBatched");
}

void LibBatchSolve::Setup()
{
   mat_size = LUMatrixBatch.SizeI();
   num_mats = LUMatrixBatch.SizeK();

   std::cout<<"mat_size = "<<mat_size<<" "<<num_mats<<std::endl;

   P.SetSize(mat_size * num_mats);
   lu_ptr_array.SetSize(num_mats);

   info_array.SetSize(num_mats);

   double *lu_ptr_base = LUMatrixBatch.ReadWrite();

   for (int i=0; i<num_mats; ++i)
   {
      lu_ptr_array.HostWrite()[i] = lu_ptr_base + i * mat_size * mat_size;
   }

   if (mode == LibBatchMode::LU) { return ComputeLU(); }
   if (mode == LibBatchMode::INVERSE) { return ComputeInverse(); }
}

void LibBatchSolve::SolveLU(const Vector &b, Vector &x)
{

   x = b;

#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)

   //Need to figure out how this should look...

   /*
   vector_array.SetSize(num_mats);

   double *x_ptr_base = x.ReadWrite();

   double alpha = 1.0;

   for(int i=0; i<num_mats; ++i) {
     vector_array.HostWrite()[i] = x_ptr_base + i * mat_size;
   }

   MFEM_SUB_cu_or_hip(blasStatus_t) status =
     MFEM_SUB_cu_or_hip(blasDtrsmBatched) (MFEM_SUB_Cuda_or_Hip(BLAS::Handle)(),
                                      MFEM_SUB_CU_or_HIP(BLAS_SIDE_LEFT),
                                      MFEM_SUB_CU_or_HIP(BLAS_FILL_MODE_UPPER),
                                      MFEM_SUB_CU_or_HIP(BLAS_OP_N),
                                      MFEM_SUB_CU_or_HIP(BLAS_DIAG_NON_UNIT),
                                      mat_size,
                                      1,
                                      &alpha,
                                      lu_ptr_array.Read(),
                                      mat_size,
                                      vector_array.ReadWrite(),
                                      mat_size);

   MFEM_SUB_cu_or_hip(blasStatus_t) status =
     MFEM_SUB_cu_or_hip(blasDtrsmBatched) (MFEM_SUB_Cuda_or_Hip(BLAS::Handle)(),
                                      MFEM_SUB_CU_or_HIP(BLAS_SIDE_LEFT),
                                      MFEM_SUB_CU_or_HIP(BLAS_FILL_MODE_LOWER),
                                      MFEM_SUB_CU_or_HIP(BLAS_OP_N),
                                      MFEM_SUB_CU_or_HIP(BLAS_DIAG_UNIT),
                                      mat_size,
                                      1,
                                      &alpha,
                                      lu_ptr_array.Read(),
                                      mat_size,
                                      vector_array.ReadWrite(),
                                      mat_size);
   */

   /*
   MFEM_SUB_cu_or_hip(blasStatus_t) status =
     MFEM_SUB_cu_or_hip(blasDgetriBatched)(MFEM_SUB_Cuda_or_Hip(BLAS::Handle)(),
                                           mat_size,
                                           lu_ptr_array.Read(),
                                           mat_size,
                                           P.Write(),
                                           vector_array.ReadWrite(),
                                           mat_size,
                                           info_array.Write(),
                                           num_mats);
   */
#endif
}

void LibBatchSolve::ApplyInverse(const Vector &b, Vector &x)
{

   std::cout<<"Applying inverse"<<std::endl;
   const double alpha = 1.0, beta = 0.0;
   MFEM_SUB_cu_or_hip(blasStatus_t)
   status = MFEM_SUB_cu_or_hip(blasDgemmStridedBatched)(MFEM_SUB_Cuda_or_Hip(
                                                           BLAS::Handle)(),
                                                        MFEM_SUB_CU_or_HIP(BLAS_OP_N),
                                                        MFEM_SUB_CU_or_HIP(BLAS_OP_N),
                                                        mat_size,
                                                        1,
                                                        mat_size,
                                                        &alpha,
                                                        InvMatrixBatch.Read(),
                                                        mat_size,
                                                        mat_size * mat_size,
                                                        b.Read(),
                                                        mat_size,
                                                        mat_size,
                                                        &beta,
                                                        x.Write(),
                                                        mat_size,
                                                        mat_size,
                                                        num_mats);

   MFEM_VERIFY(status == MFEM_SUB_CU_or_HIP(BLAS_STATUS_SUCCESS),
               "blasDgemmStridedBatched");

}

void LibBatchSolve::Mult(const Vector &b, Vector &x)
{
   if (mode == LibBatchMode::LU) { return SolveLU(b, x); }
   if (mode == LibBatchMode::INVERSE) { return ApplyInverse(b, x); }
}

/*
LibBatchMult::LibBatchMult(DenseTensor &MatrixBatch_) :
MatrixBatch(MatrixBatch_)
{
mat_size = MatrixBatch.SizeI();
num_mats = MatrixBatch.SizeK();

mat_ptr_array.SetSize(num_mats);
info_array.SetSize(num_mats);

double *mat_ptr_base = MatrixBatch.ReadWrite();

for(int i=0; i<num_mats; ++i)
{
  mat_ptr_array.HostWrite()[i] = mat_ptr_base + i * mat_size * mat_size;
}

}
*/

void LibBatchMult::Mult(const Vector &b, Vector &x)
{

   const double alpha = 1.0, beta = 0.0;
   MFEM_SUB_cu_or_hip(blasStatus_t)
   status = MFEM_SUB_cu_or_hip(blasDgemmStridedBatched)(MFEM_SUB_Cuda_or_Hip(
                                                           BLAS::Handle)(),
                                                        MFEM_SUB_CU_or_HIP(BLAS_OP_N),
                                                        MFEM_SUB_CU_or_HIP(BLAS_OP_N),
                                                        mat_size,
                                                        1,
                                                        mat_size,
                                                        &alpha,
                                                        MatrixBatch.Read(),
                                                        mat_size,
                                                        mat_size * mat_size,
                                                        b.Read(),
                                                        mat_size,
                                                        mat_size,
                                                        &beta,
                                                        x.Write(),
                                                        mat_size,
                                                        mat_size,
                                                        num_mats);

   MFEM_VERIFY(status == MFEM_SUB_CU_or_HIP(BLAS_STATUS_SUCCESS),
               "blasDgemmStridedBatched");

}


} //mfem namespace
