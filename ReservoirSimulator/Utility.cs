using SepalSolver;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    internal class Utility
    {
        /// <summary>
        /// Calculates the volume of a single tetrahedron using SIMD Vector3d registers.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double CalculateTetrahedronVolume(Vector3d v0, Vector3d v1, Vector3d v2, Vector3d v3)
        {
            // Vectorize edge profiles
            Vector3d edge1 = v1 - v0;
            Vector3d edge2 = v2 - v0;
            Vector3d edge3 = v3 - v0;

            // Triple scalar product: edge1 . (edge2 x edge3)
            Vector3d crossProduct = Vector3d.Cross(edge2, edge3);
            double tripleProduct = Vector3d.Dot(edge1, crossProduct);

            // Volume is exactly 1/6th of the parallelpiped volume
            return Abs(tripleProduct) / 6.0;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double ComputeHexahedronVolume<T>(Grid<T> grid) where T : struct, IVectorND<T>
        {
            // Sanity Check: Ensure we have exactly 8 corners
            if (grid.Corners.Length != 8) throw new ArgumentException("Hexahedron must have exactly 8 corners.");

            // Sum the volumes of the 12 constituent tetrahedrons
            double totalVolume = 0;
            foreach(var face in grid.Faces)
            {
                totalVolume += CalculateTetrahedronVolume(grid.Center, face.Corners[0], face.Corners[1], face.Corners[2]);
                totalVolume += CalculateTetrahedronVolume(grid.Center, face.Corners[0], face.Corners[2], face.Corners[3]);
            }

            // Guard against degenerate/pinched-out dead cells
            if (double.IsNaN(totalVolume) || totalVolume < 1e-6)
            {
                return 0.0;
            }

            return totalVolume;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double ComputePerpendicularDistance(Vector3d cellCenter, Vector3d p0, Vector3d p1, Vector3d p2)
        {
            // Step 1: Calculate two edge vectors from the face corners
            Vector3d v1 = p1 - p0;
            Vector3d v2 = p2 - p0;

            // Step 2 & 3: Compute the face normal and normalize it to a unit vector
            Vector3d normalLarge = Vector3d.Cross(v1, v2);
            double magnitude = normalLarge.Magnitude();

            // Guard against degenerate/zero-area faces
            if (magnitude < 1e-12) return 0.0;

            Vector3d unitNormal = (1.0 / magnitude) * normalLarge;

            // Step 4: Create a vector pointing from a face corner to the block center
            Vector3d vCenter = cellCenter - p0;

            // Step 5: Dot product projects the center vector onto the perpendicular normal track
            double distance = Vector3d.Dot(vCenter, unitNormal);

            // Return absolute distance (direction relative to the normal does not matter for transmissibility)
            return Abs(distance);
        }

        public static void ComputeFaceTransmissibility<T>(Grid<T>[] Grids) where T : struct, IVectorND<T>
        {
            // 1. Gather all faces from all active grid cells into a single flat list
            List<Face> Faces = [];
            foreach (var grid in Grids)
                Faces.AddRange(grid.Faces);

            // 2. Instantiate and build the AABB Tree (Runs in O(N log N))
            AABBTree spatialTree = new();
            spatialTree.Build(Faces);

            // 3. Query the tree for each face to resolve connectivity
            foreach (var currentFace in Faces)
            {
                List<Face> overlappingCandidates = [];

                // Pass the root of the built tree and the current target bounding box
                AABBTree.FindOverlappingFaces(spatialTree.Root, currentFace.BoundingBox, overlappingCandidates);

                // 4. Narrow Phase: Filter candidates and calculate true physical connections
                foreach (var candidate in overlappingCandidates)
                {
                    // Rule A: Skip if the face belongs to the exact same cell
                    if (candidate.ActiveGridIndex == currentFace.ActiveGridIndex)
                        continue;

                    // Rule B: Enforce index ordering constraint (current < candidate)
                    // This guarantees each cell-face pair interaction is evaluated exactly once
                    if (currentFace.ActiveGridIndex >= candidate.ActiveGridIndex)
                        continue;

                    // Rule C: Run your computational geometry suite checks
                    if (AreFacesCoplanar(currentFace, candidate))
                    {
                        // Project both sets of 3D vertices onto a 2D plane and calculate the shared area
                        double sharedArea = Calculate2DIntersectionArea(currentFace, candidate);

                        if (sharedArea > 1e-6) // Ensure it is a valid, physical contact area
                        {
                            // Register the connectivity in your global matrix/transmissibility graph
                            RegisterGridConnection(currentFace, candidate, sharedArea);
                        }
                    }
                }
            }

            foreach (var grid in Grids)
            {
                if (!grid.isActive) continue;
                foreach (var face in grid.Faces)
                    grid.indexOfConnectingGrids.AddRange(face.GridIndices);
            }
        }

        static double Harmmean(double x1, double x2)
        {
            // Guard against zero or negative permeability values
            if (x1 <= 0 || x2 <= 0) return 0.0;
            if (double.IsInfinity(x1)) return x2;
            if(double.IsInfinity(x2)) return x1;
            return 2.0 * x1 * x2 / (x1 + x2);
        }

        static bool AreFacesCoplanar(Face f1, Face f2)
        {
            // Geometric tolerance for double precision alignment (adjustable based on grid skew/warp)
            const double Tolerance = 1e-5;

            // --- Step 1: Calculate Normal Vector for Face 1 via Diagonal Cross Product ---
            // Diagonals for a 4-corner quad are strictly (2 - 0) and (3 - 1)
            Vector3d diag1A = f1.Corners[2] - f1.Corners[0];
            Vector3d diag1B = f1.Corners[3] - f1.Corners[1];

            // Cross Product: diag1A x diag1B
            Vector3d normal1 = new(
                diag1A.Y * diag1B.Z - diag1A.Z * diag1B.Y,
                diag1A.Z * diag1B.X - diag1A.X * diag1B.Z,
                diag1A.X * diag1B.Y - diag1A.Y * diag1B.X
            );

            // Check for degenerate/zero-area face
            double lenSq1 = normal1.X * normal1.X + normal1.Y * normal1.Y + normal1.Z * normal1.Z;
            if (lenSq1 < 1e-12) return false;

            double invLen1 = 1.0 / Sqrt(lenSq1);
            Vector3d unitNormal1 = new(normal1.X * invLen1, normal1.Y * invLen1, normal1.Z * invLen1);

            // --- Step 2: Calculate Normal Vector for Face 2 ---
            Vector3d diag2A = f2.Corners[2] - f2.Corners[0];
            Vector3d diag2B = f2.Corners[3] - f2.Corners[1];

            // Cross Product: diag2A x diag2B
            Vector3d normal2 = new(
                diag2A.Y * diag2B.Z - diag2A.Z * diag2B.Y,
                diag2A.Z * diag2B.X - diag2A.X * diag2B.Z,
                diag2A.X * diag2B.Y - diag2A.Y * diag2B.X
            );

            double lenSq2 = normal2.X * normal2.X + normal2.Y * normal2.Y + normal2.Z * normal2.Z;
            if (lenSq2 < 1e-12) return false;

            double invLen2 = 1.0 / Sqrt(lenSq2);
            Vector3d unitNormal2 = new(normal2.X * invLen2, normal2.Y * invLen2, normal2.Z * invLen2);

            // --- Step 3: Test for Parallel Normal Vectors ---
            // Dot product of the two unit normals
            double normalDot = Abs(unitNormal1.X * unitNormal2.X +
                                   unitNormal1.Y * unitNormal2.Y +
                                   unitNormal1.Z * unitNormal2.Z);

            if (normalDot < (1.0 - Tolerance))
            {
                return false; // Normals are not parallel; faces are not facing the same direction
            }

            // --- Step 4: Test for Same Plane Alignment (Distance Check) ---
            // Displacement vector from the first corner of Face 1 to the first corner of Face 2
            Vector3d displacement = f2.Corners[0] - f1.Corners[0];

            // Dot product between displacement and the face normal
            double planeDistanceDot = Abs(displacement.X * unitNormal1.X +
                                          displacement.Y * unitNormal1.Y +
                                          displacement.Z * unitNormal1.Z);

            // If the distance between the infinite planes is within tolerance, they are coplanar
            return planeDistanceDot <= Tolerance;
        }

        static double Calculate2DIntersectionArea(Face f1, Face f2)
        {
            // --- Step 1: Determine the Best Projection Plane ---
            // Calculate an unnormalized normal vector using Face 1's diagonals
            Vector3d diagA = f1.Corners[2] - f1.Corners[0];
            Vector3d diagB = f1.Corners[3] - f1.Corners[1];
            Vector3d normal = new Vector3d(
                diagA.Y * diagB.Z - diagA.Z * diagB.Y,
                diagA.Z * diagB.X - diagA.X * diagB.Z,
                diagA.X * diagB.Y - diagA.Y * diagB.X
            );

            double nx = Abs(normal.X);
            double ny = Abs(normal.Y);
            double nz = Abs(normal.Z);
            double len = Sqrt(nx * nx + ny * ny + nz * nz);
            if (len < 1e-12) return 0.0; // Degenerate face

            // Identify which axis is dominant to drop it (maximizes 2D area/precision)
            int dropAxis = 2; // Default drop Z (work in XY)
            if (nx > ny && nx > nz) dropAxis = 0; // Drop X (work in YZ)
            else if (ny > nx && ny > nz) dropAxis = 1; // Drop Y (work in XZ)

            // Project 3D vertices to 2D
            Vector2d[] poly1 = ProjectTo2D(f1.Corners, dropAxis);
            Vector2d[] poly2 = ProjectTo2D(f2.Corners, dropAxis);

            // --- Step 2: Clip Polygons Using Sutherland-Hodgman ---
            List<Vector2d> clippedPolygon = ClipPolygons(poly1, poly2);
            if (clippedPolygon.Count < 3) return 0.0; // No meaningful overlap

            // --- Step 3: Compute 2D Area and Scale Back to 3D Space ---
            double area2D = Compute2DArea(clippedPolygon);

            // Correction factor for the dropped projection tilt: 3D_Area = 2D_Area / cos(theta)
            // cos(theta) is equivalent to the dominant component divided by the total normal length
            double cosTheta = (dropAxis == 0 ? nx : (dropAxis == 1 ? ny : nz)) / len;

            return area2D / cosTheta;
        }
        static void RegisterGridConnection(Face f1, Face f2, double area)
        {
            double x1 = f1.Permeability / f1.PerpendicularDistance2Center,
                   x2 = f2.Permeability / f2.PerpendicularDistance2Center;
            double T = area * Harmmean(x1, x2); 
            f1.GridIndices.Add(f2.ActiveGridIndex); f1.Transmissibility.Add(T);
            f2.GridIndices.Add(f1.ActiveGridIndex); f2.Transmissibility.Add(T);
        }

        

        static Vector2d[] ProjectTo2D(Vector3d[] vertices, int dropAxis)
        {
            Vector2d[] result = new Vector2d[4];
            for (int i = 0; i < 4; i++)
            {
                if (dropAxis == 0) result[i] = new Vector2d(vertices[i].Y, vertices[i].Z);
                else if (dropAxis == 1) result[i] = new Vector2d(vertices[i].X, vertices[i].Z);
                else result[i] = new Vector2d(vertices[i].X, vertices[i].Y);
            }
            return result;
        }

        static List<Vector2d> ClipPolygons(Vector2d[] subjectPoly, Vector2d[] clipPoly)
        {
            List<Vector2d> outputList = new List<Vector2d>(subjectPoly);

            // Loop through each edge of the clipping polygon
            for (int i = 0; i < clipPoly.Length; i++)
            {
                Vector2d edgeStart = clipPoly[i];
                Vector2d edgeEnd = clipPoly[(i + 1) % clipPoly.Length];

                List<Vector2d> inputList = outputList;
                outputList = new List<Vector2d>();

                if (inputList.Count == 0) break;

                Vector2d s = inputList[inputList.Count - 1];

                foreach (Vector2d e in inputList)
                {
                    if (IsInsideEdge(e, edgeStart, edgeEnd))
                    {
                        if (!IsInsideEdge(s, edgeStart, edgeEnd))
                        {
                            outputList.Add(ComputeIntersection(s, e, edgeStart, edgeEnd));
                        }
                        outputList.Add(e);
                    }
                    else if (IsInsideEdge(s, edgeStart, edgeEnd))
                    {
                        outputList.Add(ComputeIntersection(s, e, edgeStart, edgeEnd));
                    }
                    s = e;
                }
            }
            return outputList;
        }

        static bool IsInsideEdge(Vector2d p, Vector2d edgeStart, Vector2d edgeEnd)
        {
            // Assumes clipPoly vertices are passed in anti-clockwise sequence
            return (edgeEnd.X - edgeStart.X) * (p.Y - edgeStart.Y) - (edgeEnd.Y - edgeStart.Y) * (p.X - edgeStart.X) >= -1e-10;
        }

        static Vector2d ComputeIntersection(Vector2d s, Vector2d e, Vector2d edgeStart, Vector2d edgeEnd)
        {
            double dcX = edgeStart.X - edgeEnd.X;
            double dcY = edgeStart.Y - edgeEnd.Y;
            double dpX = s.X - e.X;
            double dpY = s.Y - e.Y;

            double n1 = edgeStart.X * edgeEnd.Y - edgeStart.Y * edgeEnd.X;
            double n2 = s.X * e.Y - s.Y * e.X;
            double denominator = dcX * dpY - dcY * dpX;

            if (Abs(denominator) < 1e-12) return s; // Parallel line guard

            return new Vector2d(
                (n1 * dpX - dcX * n2) / denominator,
                (n1 * dpY - dcY * n2) / denominator
            );
        }

        static double Compute2DArea(List<Vector2d> polygon)
        {
            // Standard Shoelace Formula
            double area = 0.0;
            int j = polygon.Count - 1;
            for (int i = 0; i < polygon.Count; i++)
            {
                area += (polygon[j].X + polygon[i].X) * (polygon[j].Y - polygon[i].Y);
                j = i;
            }
            return Abs(area * 0.5);
        }

        public static List<PressureDepthPair> IntegrateLegRK4(double zStart, double pStart, double zEnd, int N, Func<double, double> getGradient)
        {
            var table = new List<PressureDepthPair>();
            double currentZ = zStart;
            double currentP = pStart;
            double h = (zEnd - zStart) / N;
            table.Add(new(currentZ, currentP));
            for (int n = 0; n < N; n++)
            {
                double k1 = getGradient(currentP);
                double k2 = getGradient(currentP + 0.5 * h * k1);
                double k3 = getGradient(currentP + 0.5 * h * k2);
                double k4 = getGradient(currentP + h * k3);
                currentP += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                currentZ += h;
                table.Add(new(currentZ, currentP));
            }
            return table;
        }

        /// <summary>
        /// Linearly interpolates ANY target column from ANY independent search column within a table.
        /// </summary>
        /// <param name="rows">The raw list of table rows from your deck parser.</param>
        /// <param name="targetValue">The known value you are querying with.</param>
        /// <param name="searchSelector">Lambda matching the independent variable column (e.g., r => r.Pc).</param>
        /// <param name="targetSelector">Lambda matching the dependent target column to retrieve (e.g., r => r.Sw).</param>
        public static double Interpolate<T>(
            List<T> rows,
            double targetValue,
            Func<T, double> searchSelector,
            Func<T, double> targetSelector)
        {
            if (rows == null)
                throw new ArgumentNullException(nameof(rows), "Table data cannot be null.");

            int count = rows.Count;
            if (count == 0)
                throw new ArgumentException("Table data cannot be empty.", nameof(rows));

            if (count == 1) return targetSelector(rows[0]);

            // CRITICAL SECURITY FIX: Order the dataset by the actual independent variable we are querying.
            // This guarantees the monotonicity required for an O(log N) binary search.
            var optimizedRows = rows.OrderBy(searchSelector).ToList();

            // Since we explicitly ordered by searchSelector, it is GUARANTEED to be ascending!
            double firstSearchVal = searchSelector(optimizedRows[0]);
            double lastSearchVal = searchSelector(optimizedRows[count - 1]);

            // 3. BOUNDARY PROTECTION: Clamp overshoots cleanly to the table edge values
            if (targetValue <= firstSearchVal) return targetSelector(optimizedRows[0]);
            if (targetValue >= lastSearchVal) return targetSelector(optimizedRows[count - 1]);

            // 4. BINARY SEARCH: Fast interval bracketing on a guaranteed ascending array
            int low = 0;
            int high = count - 1;
            int index = 0;

            while (low <= high)
            {
                int mid = low + (high - low) / 2;
                double midVal = searchSelector(optimizedRows[mid]);

                if (midVal == targetValue)
                {
                    return targetSelector(optimizedRows[mid]);
                }

                if (midVal < targetValue)
                {
                    index = mid; // Lower bound bracket candidate
                    low = mid + 1;
                }
                else
                {
                    high = mid - 1;
                }
            }

            // 5. BLEND THE BRACKETS
            var rowLeft = optimizedRows[index];
            var rowRight = optimizedRows[index + 1];

            double xLeft = searchSelector(rowLeft);
            double xRight = searchSelector(rowRight);

            double yLeft = targetSelector(rowLeft);
            double yRight = targetSelector(rowRight);

            double fraction = (targetValue - xLeft) / (xRight - xLeft);

            return yLeft + fraction * (yRight - yLeft);
        }

        public static void SetAcitveGrids<T>(Grid<T>[] Grids, List<int> Actnum, Func<ADiff<T>, ADiff<T>> Er) where T : struct, IVectorND<T>
        {
            int activeindex = 0;
            for (int m = 0; m < Grids.Length; m++)
            {
                var grid = Grids[m];
                grid.isActive = Actnum[m] == 1;
                if(grid.isActive)
                {
                    grid.Er = Er;
                    grid.MakeFaces();
                    grid.ActiveIndex = activeindex; 
                    foreach(var face in grid.Faces) 
                        face.ActiveGridIndex = activeindex;
                    activeindex++;
                }
            }
        }

        public static void SetPVDO<T>(Grid<T>[] Grids, List<int> Pvtnum,  PVDOAnalysis[] pvdo, double γo0) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Bo = Po => pvdo[Pvtnum[m]].Result.EvaluateBo(Po);
                    Grids[m].μo = Po => pvdo[Pvtnum[m]].Result.EvaluateMuo(Po);
                    Grids[m].γo = Po => γo0/pvdo[Pvtnum[m]].Result.EvaluateBo(Po);
                }
            }
        }

        public static void SetPVDG<T>(Grid<T>[] Grids, List<int> Pvtnum, PVDGAnalysis[] pvdg, double γg0) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Bg = Pg => pvdg[Pvtnum[m]].Result.EvaluateBg(Pg);
                    Grids[m].μg = Pg => pvdg[Pvtnum[m]].Result.EvaluateMug(Pg);
                    Grids[m].γg = Pg => γg0 / pvdg[Pvtnum[m]].Result.EvaluateBg(Pg);
                }
            }
        }

        public static void SetPVTO<T>(Grid<T>[] Grids, List<int> Pvtnum, PVTOAnalysis[] pvto, double γo0, double γg0) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Rs = Po => pvto[Pvtnum[m]].Result.EvaluateRs(Po);
                    Grids[m].Bo = Rs => pvto[Pvtnum[m]].Result.EvaluateBo(Rs);
                    Grids[m].μo = Rs => pvto[Pvtnum[m]].Result.EvaluateMuo(Rs);
                    Grids[m].γo = Rs => (γo0 + 0.178*Rs* γg0) / pvto[Pvtnum[m]].Result.EvaluateBo(Rs);
                }
            }
        }

        public static void SetPVTG<T>(Grid<T>[] Grids, List<int> Pvtnum, PVTGAnalysis[] pvtg, double γo0, double γg0) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Rv = Pg => pvtg[Pvtnum[m]].Result.EvaluateRv(Pg);
                    Grids[m].Bg = Rv => pvtg[Pvtnum[m]].Result.EvaluateBg(Rv);
                    Grids[m].μg = Rv => pvtg[Pvtnum[m]].Result.EvaluateMug(Rv);
                    Grids[m].γg = Rv => (γg0 + 5.615 * Rv * γo0) / pvtg[Pvtnum[m]].Result.EvaluateBg(Rv);
                }
            }
        }

        public static void SetPVTW<T>(Grid<T>[] Grids, List<int> Pvtnum, PVTWAnalysis[] pvtw, double γw0) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Bw = Pw => pvtw[Pvtnum[m]].Result.EvaluateBw(Pw);
                    Grids[m].μw = Pw => pvtw[Pvtnum[m]].Result.EvaluateMuw(Pw);
                    Grids[m].γw = Pw => γw0 / pvtw[Pvtnum[m]].Result.EvaluateBw(Pw);
                }
            }
        }

        public static void SetSWOF<T>(Grid<T>[] Grids, List<int> Satnum, SWOFAnalysis[] swof) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Krw = Sw => swof[Satnum[m]].Result.EvaluateKrw(Sw);
                    Grids[m].Krow = Sw => swof[Satnum[m]].Result.EvaluateKrow(Sw);
                    Grids[m].Pcow = Sw => swof[Satnum[m]].Result.EvaluatePcow(Sw);
                }
            }
        }

        public static void SetSGOF<T>(Grid<T>[] Grids, List<int> Satnum, SGOFAnalysis[] sgof) where T : struct, IVectorND<T>
        {
            for (int m = 0; m < Grids.Length; m++)
            {
                if (Grids[m].isActive)
                {
                    Grids[m].Krg = Sg => sgof[Satnum[m]].Result.EvaluateKrg(Sg);
                    Grids[m].Krog = Sg => sgof[Satnum[m]].Result.EvaluateKrog(Sg);
                    Grids[m].Pcgo = Sg => sgof[Satnum[m]].Result.EvaluatePcgo(Sg);
                }
            }
        }

        public static Grid<T>[] ProcessGrid<T>(int nx, int ny, int nz, double[] coord, double[] zcorn, 
            double[] Φ, double[] permx, double[] permy, double[] permz) where T : struct, IVectorND<T>
        {
            int nxny = nx * ny, ngrids = nxny * nz, nx2 = nx * 2, ny2 = ny * 2;
            Grid<T>[] Grids = new Grid<T>[ngrids];
            for (int m = 0; m < ngrids; m++)
            {
                int k = m / nx * ny, rem = m % nx * ny, j = rem / nx, i = rem % nx;
                Grids[m] = new Grid<T>
                {
                    IJKIndex = (i, j, k),
                    NTG = 1.0,
                    Porosity = Φ[m],
                    Corners = new Vector3d[8]
                };
                Grids[m].Faces[0].Permeability = permz[m];
                Grids[m].Faces[1].Permeability = permy[m];
                Grids[m].Faces[2].Permeability = permx[m];
                Grids[m].Faces[3].Permeability = permx[m];
                Grids[m].Faces[4].Permeability = permy[m];
                Grids[m].Faces[5].Permeability = permz[m];

                // Total number of coordinate lines in X and Y directions
                int numLinesX = nx + 1;

                // Loop through top (0-3) and bottom (4-7) layers
                for (int c = 0; c < 8; c++)
                {
                    // Determine local binary topological offsets
                    int ix = (c == 1 || c == 2 || c == 5 || c == 6) ? 1 : 0;
                    int iy = (c == 2 || c == 3 || c == 6 || c == 7) ? 1 : 0;
                    int iz = (c >= 4) ? 1 : 0;

                    // 1. Get the correct Coordinate Line index for this corner
                    int lineI = i + ix, lineJ = j + iy;
                    int lineStride = lineI + lineJ * numLinesX;

                    // Extract Top and Bottom anchor points from the raw COORD array
                    // Layout of COORD: [xtop, ytop, ztop, xbot, ybot, zbot] for each line
                    int coordIdx = lineStride * 6;
                    Vector3d pTop = new(coord[coordIdx + 0], coord[coordIdx + 1], coord[coordIdx + 2]);
                    Vector3d pBot = new(coord[coordIdx + 3], coord[coordIdx + 4], coord[coordIdx + 5]);

                    // 2. Fetch the explicit target Z-coordinate from ZCORN array
                    // Layout of ZCORN uses independent 2nx * 2ny * 2nz index structure
                    int zIdx_X = i * 2 + ix, zIdx_Y = j * 2 + iy, zIdx_Z = k * 2 + iz;

                    int zcornIdx = zIdx_X + (zIdx_Y * nx2) + (zIdx_Z * nx2 * ny2);
                    double targetZ = zcorn[zcornIdx];

                    // 3. Linear interpolation along the line vector using register math
                    double fraction = (targetZ - pTop.Z) / (pBot.Z - pTop.Z);

                    // Check for dead/pinch-out cells where top and bottom meet
                    if (double.IsNaN(fraction) || double.IsInfinity(fraction)) fraction = 0.0;

                    // X and Y are calculated seamlessly via your custom Vector3d scalar multiplication operator
                    Vector3d xyPlaneInterpolated = pTop + fraction * (pBot - pTop);

                    // Save the exact 3D location (Overriding the Z lane with the explicit target depth)
                    Grids[m].Corners[c] = new Vector3d(xyPlaneInterpolated.X, xyPlaneInterpolated.Y, targetZ);
                }
            }
            return Grids;
        }
        public static Grid<T>[] ProcessGrid<T>(int nx, int ny, int nz, double[] dx, double[] dy, double[] dz,
            double[] zTop, double[] Φ, double[] permx, double[] permy, double[] permz) where T : struct, IVectorND<T>
        {
            int nxny = nx * ny, ngrids = nxny * nz;
            Grid<T>[] Grids = new Grid<T>[ngrids];
            for (int m = 0; m < ngrids; m++)
            {
                int k = m / nxny, rem = m % nxny,
                    j = rem / nx, i = rem % nx;
                var K = Enumerable.Range(0, k);
                var J = Enumerable.Range(0, j);
                var I = Enumerable.Range(0, i);
                Grids[m] = new Grid<T>
                {
                    IJKIndex = (i, j, k),
                    NTG = 1.0,
                    Porosity = Φ[m],
                    Corners = new Vector3d[8]
                };
                Grids[m].Faces[0].Permeability = permz[m];
                Grids[m].Faces[1].Permeability = permy[m];
                Grids[m].Faces[2].Permeability = permx[m];
                Grids[m].Faces[3].Permeability = permx[m];
                Grids[m].Faces[4].Permeability = permy[m];
                Grids[m].Faces[5].Permeability = permz[m];
                double zlo = zTop[j * nx + i] + K.Select(k => dz[i + j * nx + k * nxny]).Sum();
                double zhi = zlo - dz[m];
                double xhi = I.Select(i => dx[i + j * nx + k * nxny]).Sum();
                double xlo = xhi - dx[m];
                double yhi = J.Select(j => dy[i + j * nx + k * nxny]).Sum();
                double ylo = yhi - dy[m];

                // Extract the 8 structural corner points for a Hexahedron matching VTK specifications
                // Bottom plane points (Z-Upper)
                Grids[m].Corners[0] = new Vector3d(xlo, ylo, zhi);     // Point 0
                Grids[m].Corners[1] = new Vector3d(xhi, ylo, zhi);     // Point 1
                Grids[m].Corners[2] = new Vector3d(xhi, yhi, zhi);     // Point 2
                Grids[m].Corners[3] = new Vector3d(xlo, yhi, zhi);     // Point 3
                // Top plane points (Z-Lower)
                Grids[m].Corners[4] = new Vector3d(xlo, ylo, zlo);     // Point 4
                Grids[m].Corners[5] = new Vector3d(xhi, ylo, zlo);     // Point 5
                Grids[m].Corners[6] = new Vector3d(xhi, yhi, zlo);     // Point 6
                Grids[m].Corners[7] = new Vector3d(xlo, yhi, zlo);     // Point 7
            }
            return Grids;
        }
    }
}
