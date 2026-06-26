using System.Runtime.InteropServices;

namespace ReservoirSimulator
{
    // --- 1. THE NATIVE INTEROP WRAPPER ---
    public static class MklNative
    {
        private const string MklRtDll = "mkl_rt";

        [DllImport(MklRtDll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "pardiso")]
        public static unsafe extern void Pardiso(
            IntPtr* pt,             // Internal solver memory pointers (64 elements)
            int* maxfct,            // Maximum number of factorizations
            int* mnum,              // Matrix identity index (usually 1)
            int* mtype,             // Matrix type (11 = real unsymmetric)
            int* phase,             // Process step (11=Analysis, 22=Factorization, 33=Solve)
            int* n,                 // System dimension size
            double* a,              // CSR matrix values
            int* ia,                // CSR row tracking array
            int* ja,                // CSR column index array
            int* perm,              // Permutation vector array
            int* nrhs,              // Right-hand sides quantity
            int* iparm,             // Control parameters config (64 elements)
            int* msglvl,            // Output detail setting (0=quiet)
            double* b,              // RHS system data input
            double* x,              // Output system destination target
            int* error              // Error feedback code target
        );
    }
    public static unsafe class MklSparseSolver
    {
        public static double[] Solve(double[] values, int[] colIndices, int[] rowOffsets, double[] b)
        {
            int n = b.Length;
            double[] x = new double[n];

            // Allocate PARDISO internal memory handle and control parameters
            IntPtr* pt = stackalloc IntPtr[64];
            for (int i = 0; i < 64; i++) pt[i] = IntPtr.Zero;

            int* iparm = stackalloc int[64];
            for (int i = 0; i < 64; i++) iparm[i] = 0;

            // --- Step 2a: Configure PARDISO Parameter Options ---
            iparm[0] = 1;         // 1 = User-defined parameters (Do not use defaults)
            iparm[1] = 2;         // 2 = Nested dissection from METIS fill-in reduction
            iparm[9] = 13;        // Pivot perturbation 10^(-13)
            iparm[10] = 1;        // Enable scaling vectors
            iparm[12] = 1;        // Enable matching non-symmetric permutation
            iparm[17] = -1;       // Output: Number of non-zeros in factors
            iparm[18] = -1;       // Output: Mflops of factorization
            iparm[34] = 1;        // **CRITICAL**: 0-based indexing for C# (CSR arrays use 0, 1, 2...)

            // --- Step 2b: Execution Settings ---
            int maxfct = 1;       // Max number of factorizations
            int mnum = 1;         // Matrix number to track
            int mtype = 11;       // 11 = Real structurally unsymmetric matrix
            int nrhs = 1;         // Number of right-hand sides
            int msglvl = 0;       // 0 = Suppress stdout terminal messages
            int error = 0;

            // Dummy permutation array required by signature
            int[] perm = new int[n];

            // Pin managed arrays explicitly to secure safe memory addresses for the native boundary
            fixed (double* pValues = values, pB = b, pX = x)
            fixed (int* pColIndices = colIndices, pRowOffsets = rowOffsets, pPerm = perm)
            {
                // Phase 11: Symbolic Analysis and Reordering
                int phase = 11;
                MklNative.Pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pValues, pRowOffsets, pColIndices, pPerm, &nrhs, iparm, &msglvl, pB, pX, &error);
                if (error != 0) throw new Exception($"PARDISO Analysis Phase (11) failed with error code: {error}");

                // Phase 22: Numerical Factorization
                phase = 22;
                MklNative.Pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pValues, pRowOffsets, pColIndices, pPerm, &nrhs, iparm, &msglvl, pB, pX, &error);
                if (error != 0) throw new Exception($"PARDISO Factorization Phase (22) failed with error code: {error}");

                // Phase 33: Forward and Backward Substitution (The actual solve)
                phase = 33;
                MklNative.Pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pValues, pRowOffsets, pColIndices, pPerm, &nrhs, iparm, &msglvl, pB, pX, &error);
                if (error != 0) throw new Exception($"PARDISO Solve Phase (33) failed with error code: {error}");

                // Phase -1: Release all internal memory allocations
                phase = -1;
                MklNative.Pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pValues, pRowOffsets, pColIndices, pPerm, &nrhs, iparm, &msglvl, pB, pX, &error);
            }

            return x;
        }
    }
}
