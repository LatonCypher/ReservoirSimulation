using SepalSolver;
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace ReservoirSimulator
{
    public static class ReservoirJacobian
    {
        public static (int[] col_index, int[] row_ptr, double[] values) Create<T>(
            int nc,          // 2 for Vector4d, 3 for Vector6d, 4 for Vector8d
            IReadOnlyList<Grid<T>> Grids,
            IReadOnlyList<Well<T>> Wells) where T : struct, IVectorND<T>
        {
            int numWells = Wells.Count;
            int numActiveGrids = Grids.Count(grid => grid.isActive);
            int totalDofs = nc * numActiveGrids + 2 * numWells;

            // 1. Exact capacity estimation
            int totalGridBlocks = numActiveGrids; // Diagonal blocks
            foreach (var mgrid in Grids)
            {
                if (!mgrid.isActive) continue; // Corrected: skip inactive cells

                foreach (var face in mgrid.Faces)
                {
                    foreach (int c in face.GridIndices)
                    {
                        if (Grids[c].isActive) totalGridBlocks++;
                    }
                }
            }

            int totalActivePerfs = 0;
            foreach (var well in Wells)
            {
                foreach (int c in well.PerfNatIndex)
                {
                    if (Grids[c].isActive) totalActivePerfs++;
                }
            }

            // Exact NNZ capacity:
            // - Grid blocks: totalGridBlocks * (nc * nc)
            // - Well internal blocks: numWells * 4 (2x2 closure block)
            // - Perforation blocks: totalActivePerfs * (2 * nc [Grid-Well] + 2 * nc [Well-Grid]) = totalActivePerfs * (4 * nc)
            int capacity = (totalGridBlocks * nc * nc)
                         + (numWells * 4)
                         + (totalActivePerfs * 4 * nc);

            List<int> I = new(capacity);
            List<int> J = new(capacity);

            // 2. Exact AddBlock mechanics (MOL style)
            void AddBlock(int i, int j)
            {
                int I0 = i * nc;
                int J0 = j * nc;

                for (int ii = 0; ii < nc; ii++)
                {
                    for (int jj = 0; jj < nc; jj++)
                    {
                        I.Add(I0 + ii);
                        J.Add(J0 + jj);
                    }
                }
            }

            // Coupling: Grid cell equation depends on Well variables (Nc x 2)
            void AddGridWellBlock(int gridActiveIdx, int wellIdx)
            {
                int I0 = gridActiveIdx * nc;
                int J0 = nc * numActiveGrids + 2 * wellIdx;

                for (int ii = 0; ii < nc; ii++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        I.Add(I0 + ii);
                        J.Add(J0 + jj);
                    }
                }
            }

            // Coupling: Well rate equation depends on Grid variables (2 x Nc)
            void AddWellGridBlock(int wellIdx, int gridActiveIdx)
            {
                int I0 = nc * numActiveGrids + 2 * wellIdx;
                int J0 = gridActiveIdx * nc;

                for (int ii = 0; ii < 2; ii++)
                {
                    for (int jj = 0; jj < nc; jj++)
                    {
                        I.Add(I0 + ii);
                        J.Add(J0 + jj);
                    }
                }
            }

            // Well internal closure block (2 x 2)
            void AddWellInternalBlock(int wellIdx)
            {
                int I0 = nc * numActiveGrids + 2 * wellIdx;

                for (int ii = 0; ii < 2; ii++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        I.Add(I0 + ii);
                        J.Add(I0 + jj);
                    }
                }
            }

            // ---------------------------------------------------------
            // 3. Grid Cells and Face Neighbor Connections
            // ---------------------------------------------------------
            foreach (var mgrid in Grids)
            {
                if (!mgrid.isActive) continue;
                int i = mgrid.ActiveIndex;

                // Diagonal accumulation block (i, i)
                AddBlock(i, i);

                // Connected face neighbor blocks (i, j)
                foreach (var face in mgrid.Faces)
                {
                    foreach (int c in face.GridIndices)
                    {
                        var ngrid = Grids[c];
                        if (!ngrid.isActive) continue;
                        int j = ngrid.ActiveIndex;
                        AddBlock(i, j);
                    }
                }
            }

            // ---------------------------------------------------------
            // 4. Wells Internal Closure and Perforations
            // ---------------------------------------------------------
            for (int w = 0; w < numWells; w++)
            {
                var well = Wells[w];

                // Well internal block (2 x 2)
                AddWellInternalBlock(w);

                // Perforation couplings
                foreach (int m in well.PerfNatIndex)
                {
                    var mgrid = Grids[m];
                    if (!mgrid.isActive) continue;
                    int i = mgrid.ActiveIndex;

                    AddGridWellBlock(i, w);
                    AddWellGridBlock(w, i);
                }
            }

            // Match your SepalSolver constructor
            double[] ones = SepalSolver.Math.Ones(I.Count);
            var jac = new SparseMatrix(I, J, ones, totalDofs, totalDofs);
            

            string filePath = "C:\\Users\\lateef.a.kareem\\Documents\\AbimbolaPeter\\output.m";

            // Open the file stream and redirect standard output
            using (var fileStream = new FileStream(filePath, FileMode.Create, FileAccess.Write))
            using (var writer = new StreamWriter(fileStream))
            {
                // AutoFlush ensures text writes to disk immediately without waiting for buffer fills
                writer.AutoFlush = true;
                Console.SetOut(writer);
                Console.WriteLine($"I = [{string.Join(", ", I)}];");
                Console.WriteLine($"J = [{string.Join(", ", J)}];");
                Console.WriteLine($"V = [{string.Join(", ", ones)}];");
                Console.WriteLine("I = I+1; J = J + 1;");
                Console.WriteLine($"s = sparse(I, J, V, {totalDofs}, {totalDofs});");
                Console.WriteLine("spy(s);");
            }
            var crs = jac.CompressByRow();
            return ([..crs.col_index], [..crs.row_start], [..crs.values]);
        }


        public readonly struct CellOffsets
        {
            // Flattened indices into values[] for diagonal block (i, i) of size nc x nc
            // Index mapping: [eq * nc + var]
            public readonly int[] Diag;

            public CellOffsets(int nc) => Diag = new int[nc * nc];
        }

        public readonly struct FaceOffsets
        {
            public readonly int NeighborActiveIndex;
            // Flattened indices into values[]:
            // Block (i, i) sensitivities (slots 0 .. nc-1): [eq * nc + var]
            public readonly int[] Self;
            // Block (i, j) sensitivities (slots nc .. 2*nc-1): [eq * nc + var]
            public readonly int[] Neighbor;

            public FaceOffsets(int nc, int neighborActiveIdx)
            {
                NeighborActiveIndex = neighborActiveIdx;
                Self = new int[nc * nc];
                Neighbor = new int[nc * nc];
            }
        }

        public readonly struct WellOffsets
        {
            public readonly int RowP_Pwell;
            public readonly int RowP_Qwell;
            public readonly int RowQ_Pwell;
            public readonly int RowQ_Qwell;

            // Per perforation:
            // [perfIdx]: grid cell equation sensitivities to Pwell: size nc
            public readonly int[][] Rm_dPwell;
            // [perfIdx]: well rate equation sensitivities to grid cell variables: size nc
            public readonly int[][] RowQ_dGrid;

            public WellOffsets(int rpP, int rpQ, int rqP, int rqQ, int[][] rmDPwell, int[][] rqDGrid)
            {
                RowP_Pwell = rpP; RowP_Qwell = rpQ;
                RowQ_Pwell = rqP; RowQ_Qwell = rqQ;
                Rm_dPwell = rmDPwell;
                RowQ_dGrid = rqDGrid;
            }
        }


        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int FindValueIndex(int[] row_ptr, int[] col_index, int row, int col)
        {
            int start = row_ptr[row];
            int length = row_ptr[row + 1] - start;
            int idx = Array.BinarySearch(col_index, start, length, col);
            return idx;
        }




        public static (CellOffsets[] cellOffsets, FaceOffsets[][] faceOffsets, WellOffsets[] wellOffsets)
    PrecomputeOffsets<T>(
        int nc,
        int[] row_ptr,
        int[] col_index,
        IReadOnlyList<Grid<T>> Grids,
        IReadOnlyList<Well<T>> Wells) where T : struct, IVectorND<T>
        {
            int numWells = Wells.Count;
            int numActiveGrids = Grids.Count(grid => grid.isActive);
            var cellOffsets = new CellOffsets[numActiveGrids];
            var faceOffsets = new FaceOffsets[Grids.Count][];
            var wellOffsets = new WellOffsets[numWells];

            // 1. Grid Diagonal & Face Offsets
            for (int m = 0; m < Grids.Count; m++)
            {
                var mgrid = Grids[m];
                if (!mgrid.isActive) continue;

                int i = mgrid.ActiveIndex;
                int baseRow_i = nc * i;
                int baseCol_i = nc * i;

                // Diagonal Block (i, i)
                var cOff = new CellOffsets(nc);
                for (int eq = 0; eq < nc; eq++)
                {
                    for (int v = 0; v < nc; v++)
                    {
                        cOff.Diag[eq * nc + v] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, baseCol_i + v);
                    }
                }
                cellOffsets[i] = cOff;

                // Faces
                int faceCount = mgrid.Faces.Length;
                faceOffsets[m] = new FaceOffsets[faceCount];

                for (int f = 0; f < faceCount; f++)
                {
                    var face = mgrid.Faces[f];
                    if (face.GridIndices.Count == 0) continue;

                    int n = face.GridIndices[0];
                    var ngrid = Grids[n];
                    if (!ngrid.isActive) continue;

                    int j = ngrid.ActiveIndex;
                    int baseCol_j = nc * j;

                    var fOff = new FaceOffsets(nc, j);
                    for (int eq = 0; eq < nc; eq++)
                    {
                        for (int v = 0; v < nc; v++)
                        {
                            // Block (i, i)
                            fOff.Self[eq * nc + v] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, baseCol_i + v);
                            // Block (i, j)
                            fOff.Neighbor[eq * nc + v] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, baseCol_j + v);
                        }
                    }
                    faceOffsets[m][f] = fOff;
                }
            }

            // 2. Well Offsets
            for (int w = 0; w < numWells; w++)
            {
                var well = Wells[w];
                int row_p = nc * numActiveGrids + 2 * w + 0;
                int row_q = nc * numActiveGrids + 2 * w + 1;
                int col_pwell = row_p;
                int col_qwell = row_q;

                int rpP = FindValueIndex(row_ptr, col_index, row_p, col_pwell);
                int rpQ = FindValueIndex(row_ptr, col_index, row_p, col_qwell);
                int rqP = FindValueIndex(row_ptr, col_index, row_q, col_pwell);
                int rqQ = FindValueIndex(row_ptr, col_index, row_q, col_qwell);

                int perfCount = well.PerfNatIndex.Count;
                var rmDPwell = new int[perfCount][];
                var rqDGrid = new int[perfCount][];

                for (int p = 0; p < perfCount; p++)
                {
                    int m = well.PerfNatIndex[p];
                    var mgrid = Grids[m];
                    if (!mgrid.isActive) continue;

                    int i = mgrid.ActiveIndex;
                    int baseRow_i = nc * i;
                    int baseCol_i = nc * i;

                    rmDPwell[p] = new int[nc];
                    rqDGrid[p] = new int[nc];

                    for (int eq = 0; eq < nc; eq++)
                    {
                        rmDPwell[p][eq] = FindValueIndex(row_ptr, col_index, baseRow_i + eq, col_pwell);
                        rqDGrid[p][eq] = FindValueIndex(row_ptr, col_index, row_q, baseCol_i + eq);
                    }
                }

                wellOffsets[w] = new WellOffsets(rpP, rpQ, rqP, rqQ, rmDPwell, rqDGrid);
            }

            return (cellOffsets, faceOffsets, wellOffsets);
        }
    }
}