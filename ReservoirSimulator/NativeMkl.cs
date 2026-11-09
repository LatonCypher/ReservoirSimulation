using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace ReservoirSimulator
{
    public static unsafe class NativeMkl
    {
        private const string MklLib = "mkl_rt";

        [DllImport(MklLib, CallingConvention = CallingConvention.Cdecl, EntryPoint = "pardiso")]
        public static extern void pardiso(
            IntPtr* pt,         // Internal solver memory pointer array (size 64)
            int* maxfct,        // Maximum number of numerical factorizations (usually 1)
            int* mnum,          // Which factorization to use (usually 1)
            int* mtype,         // Matrix type (11 = real and unsymmetric)
            int* phase,         // Solver execution phase (11, 22, 33, 23, -1)
            int* n,             // Number of equations / dimension
            double* a,          // Matrix values
            int* ia,            // Row pointer array (1-based for PARDISO!)
            int* ja,            // Column index array (1-based for PARDISO!)
            int* perm,          // Permutation vector (can be null if iparm[4] = 0)
            int* nrhs,          // Number of right-hand sides (1)
            int* iparm,         // Parameter array (size 64)
            int* msglvl,        // Message level (0 = quiet, 1 = print stats)
            double* b,          // Right-hand side vector
            double* x,          // Solution vector
            int* error          // Error flag
        );
    }
}
