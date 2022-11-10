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

#ifndef MFEM_LIB_BATCH_SOLVER
#define MFEM_LIB_BATCH_SOLVER

#include "../config/config.hpp"
#include "../general/globals.hpp"
#include "matrix.hpp"
#include "densemat.hpp"

#include <cublas.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>

#if defined(MFEM_USE_CUDA)
#define MFEM_ENABLE_CUDA_OR_HIP
// These macros prefixes the argument with `cuda` or `hip`, e.g. `cudaStream_t` or `hipStream_t`
#define MFEM_SUB_cuda_or_hip(stub) cuda##stub
#define MFEM_SUB_cu_or_hip(stub) cu##stub
#define MFEM_SUB_Cuda_or_Hip(stub) Cuda##stub
#define MFEM_SUB_CUDA_or_HIP(stub) CUDA##stub
#define MFEM_SUB_CU_or_HIP(stub) CU##stub

#endif  // MFEM_USE_CUDA

#if defined(MFEM_USE_HIP)
#define MFEM_ENABLE_CUDA_OR_HIP
// These macros prefixes the argument with `cuda` or `hip`, e.g. `cudaStream_t` or `hipStream_t`
#define MFEM_SUB_cu_or_hip(stub) hip##stub
#define MFEM_SUB_cuda_or_hip(stub) hip##stub
#define MFEM_SUB_Cuda_or_Hip(stub) Hip##stub
#define MFEM_SUB_CUDA_or_HIP(stub) HIP##stub
#define MFEM_SUB_CU_or_HIP(stub) HIP##stub
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#endif  // MFEM_USE_HIP

namespace mfem
{

class MFEM_SUB_Cuda_or_Hip(BLAS)
{
protected:
   MFEM_SUB_cu_or_hip(blasHandle_t) handle = nullptr;

   MFEM_SUB_Cuda_or_Hip(BLAS)()
   {
      MFEM_SUB_cu_or_hip(blasStatus_t) status = MFEM_SUB_cu_or_hip(blasCreate)(
                                                   &handle);
      MFEM_VERIFY(status == MFEM_SUB_CU_or_HIP(BLAS_STATUS_SUCCESS),
                  "Cannot initialize BLAS.");
   }

   ~MFEM_SUB_Cuda_or_Hip(BLAS)() { MFEM_SUB_cu_or_hip(blasDestroy)(handle); }

   static MFEM_SUB_Cuda_or_Hip(BLAS) & Instance()
   {
      static MFEM_SUB_Cuda_or_Hip(BLAS) instance;
      return instance;
   }

public:
   static MFEM_SUB_cu_or_hip(blasHandle_t) Handle() { return Instance().handle; }
};

enum struct LibBatchMode : int {LU, INVERSE};

/**
@Brief use vendor libraries to perfom LU factorization and solves
**/
class LibBatchSolve
{

private:
   LibBatchMode mode;

   DenseTensor LUMatrixBatch, InvMatrixBatch;

   Array<int> P, info_array;

   Array<double *> lu_ptr_array, inv_ptr_array, vector_array;

   void ComputeLU();

   void ComputeInverse();

   void SolveLU(const Vector &x, Vector &b);

   void ApplyInverse(const Vector &x, Vector &b);

   int num_mats, mat_size;

public:
   LibBatchSolve() = delete;

   LibBatchSolve(LibBatchMode mode_) : mode(mode_) {};

   LibBatchSolve(DenseTensor &MatrixBatch_, LibBatchMode mode_)
      : LUMatrixBatch(MatrixBatch_), mode(mode_) {};

   void AssignMatrices(DenseTensor &MatrixBatch_);

   //Compute LU or Inverse
   void Setup();

   //Solve linear system Ax = b
   void Mult(const Vector &b, Vector &x);

};

class LibBatchMult
{

private:

   DenseTensor &MatrixBatch;

   int mat_size, num_mats;

public:
   LibBatchMult() = delete;

   LibBatchMult(DenseTensor &MatrixBatch_) :
      mat_size(MatrixBatch.SizeI()),
      num_mats(MatrixBatch.SizeK()),
      MatrixBatch(MatrixBatch_) {};

   //Mult linear system Ax = b
   void Mult(const Vector &b, Vector &x);

};



} // namespace mfem

#endif
