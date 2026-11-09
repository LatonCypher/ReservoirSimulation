using MKLNET;
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

namespace ReservoirSimulator
{
    public sealed unsafe class MklPardisoSolver : IDisposable
    {
        private readonly int _n;
        private readonly int[] _ia; // 1-based row_ptr
        private readonly int[] _ja; // 1-based col_index

        private readonly IntPtr[] _pt = new IntPtr[64];
        private readonly int[] _iparm = new int[64];

        private int _mtype = 11; // 11 = Real and unsymmetric matrix
        private int _maxfct = 1;
        private int _mnum = 1;
        private int _msglvl = 0; // 0 for production, 1 for debug diagnostic output
        private int _nrhs = 1;

        public MklPardisoSolver(int n, int[] row_ptr_0based, int[] col_index_0based)
        {
            _n = n;

            // 1. Shift CSR arrays to 1-based indexing
            _ia = new int[row_ptr_0based.Length];
            _ja = new int[col_index_0based.Length];

            for (int i = 0; i < row_ptr_0based.Length; i++)
                _ia[i] = row_ptr_0based[i] + 1;

            for (int j = 0; j < col_index_0based.Length; j++)
                _ja[j] = col_index_0based[j] + 1;

            // 2. Clear internal solver memory pointer array completely
            Array.Clear(_pt, 0, _pt.Length);

            // =====================================================================
            // CALL VALIDATION HERE (Before Native PARDISO touches the graph)
            // =====================================================================
            ValidateCsrStructure();

            // 3. Setup default PARDISO parameters
            _iparm[0] = 1;  // Use non-default values
            _iparm[1] = 2;  // Nested Dissection (METIS)
                            // ... remaining iparm setup ...

            // 4. Phase 11: Symbolic analysis & reordering
            int phase = 11;
            int error = 0;

            fixed (IntPtr* pPt = _pt)
            fixed (int* pIparm = _iparm)
            fixed (int* pIa = _ia)
            fixed (int* pJa = _ja)
            {
                int maxfct = _maxfct, mnum = _mnum, mtype = _mtype, dim = _n, nrhs = _nrhs, msglvl = _msglvl;

                NativeMkl.pardiso(
                    pPt, &maxfct, &mnum, &mtype, &phase, &dim,
                    null, pIa, pJa, null, &nrhs, pIparm, &msglvl, null, null, &error);
            }

            if (error != 0)
                throw new InvalidOperationException($"PARDISO Phase 11 failed with error code: {error}");
        }

        private void ValidateCsrStructure()
        {
            if (_ia == null || _ja == null)
                throw new ArgumentNullException("ia or ja is null");

            if (_ia.Length != _n + 1)
                throw new InvalidOperationException($"_ia must have length N + 1. Expected {_n + 1}, got {_ia.Length}");

            if (_ia[0] != 1)
                throw new InvalidOperationException($"1-based CSR requires _ia[0] == 1. Got {_ia[0]}");

            if (_ia[_n] != _ja.Length + 1)
                throw new InvalidOperationException($"_ia[N] must equal NNZ + 1. Expected {_ja.Length + 1}, got {_ia[_n]}");

            for (int i = 0; i < _n; i++)
            {
                int rowStart = _ia[i] - 1;
                int rowEnd = _ia[i + 1] - 1;

                if (rowStart > rowEnd)
                    throw new InvalidOperationException($"Row {i} has negative length (start: {rowStart}, end: {rowEnd})");

                if (rowStart == rowEnd)
                    throw new InvalidOperationException($"Row {i} is completely empty (no column entries)!");

                for (int k = rowStart; k < rowEnd; k++)
                {
                    int col = _ja[k];
                    if (col < 1 || col > _n)
                        throw new InvalidOperationException($"Row {i} contains invalid column index {col} at index {k}. Must be 1 <= col <= {_n}");

                    if (k < rowEnd - 1 && col >= _ja[k + 1])
                        throw new InvalidOperationException($"Row {i} has duplicate or unsorted column indices: {col} and {_ja[k + 1]}");
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Solve(double[] values, double[] rhs, double[] solution)
        {
            // Phase 23: Numerical Factorization (Phase 22) + Forward/Backward Solve (Phase 33)
            int phase = 23;
            int error = 0;

            fixed (IntPtr* pPt = _pt)
            fixed (int* pIparm = _iparm)
            fixed (int* pIa = _ia)
            fixed (int* pJa = _ja)
            fixed (double* pA = values)
            fixed (double* pB = rhs)
            fixed (double* pX = solution)
            {
                int maxfct = _maxfct, mnum = _mnum, mtype = _mtype, dim = _n, nrhs = _nrhs, msglvl = _msglvl;

                NativeMkl.pardiso(
                    pPt, &maxfct, &mnum, &mtype, &phase, &dim,
                    pA, pIa, pJa, null, &nrhs, pIparm, &msglvl, pB, pX, &error);
            }

            if (error != 0)
                throw new InvalidOperationException($"PARDISO Solve failed with error code: {error}");
        }

        public void Dispose()
        {
            // Phase -1: Release all internal solver memory
            int phase = -1;
            int error = 0;

            fixed (IntPtr* pPt = _pt)
            fixed (int* pIparm = _iparm)
            {
                int maxfct = _maxfct, mnum = _mnum, mtype = _mtype, dim = _n, nrhs = _nrhs, msglvl = 0;

                NativeMkl.pardiso(
                    pPt, &maxfct, &mnum, &mtype, &phase, &dim,
                    null, null, null, null, &nrhs, pIparm, &msglvl, null, null, &error);
            }
        }
    }
}
