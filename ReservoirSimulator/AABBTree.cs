using System;
using System.Collections.Generic;
using SepalSolver;

namespace ReservoirSimulator
{
    public class AABBNode
    {
        public AABB BoundingBox { get; set; }
        public AABBNode Left { get; set; }
        public AABBNode Right { get; set; }
        public Face Face { get; set; }

        public bool IsLeaf => Face != null;
    }

    public readonly struct AABB
    {
        public readonly Vector3d Min;
        public readonly Vector3d Max;

        public static AABB Empty => new(
            new Vector3d(double.MaxValue, double.MaxValue, double.MaxValue),
            new Vector3d(double.MinValue, double.MinValue, double.MinValue)
        );

        public AABB(Vector3d min, Vector3d max)
        {
            Min = min;
            Max = max;
        }

        public AABB(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
        {
            Min = new Vector3d(minX, minY, minZ);
            Max = new Vector3d(maxX, maxY, maxZ);
        }

        public static AABB CreateMerged(in AABB b1, in AABB b2)
        {
            return new AABB(
                Vector3d.Min(b1.Min, b2.Min),
                Vector3d.Max(b1.Max, b2.Max)
            );
        }

        public Vector3d Centroid => (Min + Max) * 0.5;
        public Vector3d Size => Max - Min;

        public bool Intersects(in AABB target)
        {
            return (Min.X <= target.Max.X && Max.X >= target.Min.X) &&
                   (Min.Y <= target.Max.Y && Max.Y >= target.Min.Y) &&
                   (Min.Z <= target.Max.Z && Max.Z >= target.Min.Z);
        }
    }

    public class AABBTree
    {
        public AABBNode Root { get; private set; }

        public void Build(List<Face> faces)
        {
            if (faces == null || faces.Count == 0) return;

            // Single contiguous array to partition in-place (zero subsequent allocations)
            Face[] faceArray = faces.ToArray();
            Root = BuildRecursive(faceArray, 0, faceArray.Length);
        }

        private static AABBNode BuildRecursive(Face[] faces, int start, int count)
        {
            var node = new AABBNode();

            // 1. Leaf base case: single element
            if (count == 1)
            {
                node.Face = faces[start];
                node.BoundingBox = faces[start].BoundingBox;
                return node;
            }

            // 2. Leaf base case: 2 elements (trivial split without QuickSelect)
            if (count == 2)
            {
                node.Left = new AABBNode { Face = faces[start], BoundingBox = faces[start].BoundingBox };
                node.Right = new AABBNode { Face = faces[start + 1], BoundingBox = faces[start + 1].BoundingBox };
                node.BoundingBox = AABB.CreateMerged(node.Left.BoundingBox, node.Right.BoundingBox);
                return node;
            }

            // 3. Compute centroid bounds
            Vector3d cMin = new(double.MaxValue, double.MaxValue, double.MaxValue);
            Vector3d cMax = new(double.MinValue, double.MinValue, double.MinValue);

            for (int i = start; i < start + count; i++)
            {
                Vector3d c = faces[i].BoundingBox.Centroid;
                cMin = Vector3d.Min(cMin, c);
                cMax = Vector3d.Max(cMax, c);
            }

            Vector3d cSpan = cMax - cMin;

            int mid = start + count / 2;

            // 4. CRITICAL FIX: If all centroids are spatially identical (cSpan near 0 on all axes),
            // skip spatial sorting and do a clean 50/50 index split to prevent infinite recursion!
            if (cSpan.X < 1e-9 && cSpan.Y < 1e-9 && cSpan.Z < 1e-9)
            {
                // Already mid = start + count / 2
            }
            else
            {
                // Pick split axis with largest extent
                int axis = 0;
                if (cSpan.Y > cSpan.X && cSpan.Y > cSpan.Z) axis = 1;
                else if (cSpan.Z > cSpan.X && cSpan.Z > cSpan.Y) axis = 2;

                // In-place QuickSelect
                QuickSelect(faces, start, start + count - 1, mid, axis);
            }

            int leftCount = mid - start;
            int rightCount = count - leftCount;

            // Safety guard: guarantee strictly smaller subproblems
            if (leftCount <= 0 || rightCount <= 0)
            {
                leftCount = count / 2;
                mid = start + leftCount;
                rightCount = count - leftCount;
            }

            // 5. Recursively build child nodes
            node.Left = BuildRecursive(faces, start, leftCount);
            node.Right = BuildRecursive(faces, mid, rightCount);

            node.BoundingBox = AABB.CreateMerged(node.Left.BoundingBox, node.Right.BoundingBox);
            return node;
        }

        // Robust Hoare QuickSelect that never hangs on duplicate keys
        private static void QuickSelect(Face[] faces, int left, int right, int k, int axis)
        {
            while (left < right)
            {
                int pivotIdx = Partition(faces, left, right, axis);
                if (pivotIdx == k) return;
                if (k < pivotIdx)
                    right = pivotIdx - 1;
                else
                    left = pivotIdx + 1;
            }
        }

        private static int Partition(Face[] faces, int left, int right, int axis)
        {
            int mid = left + (right - left) / 2;
            double pivotVal = GetCentroidAxis(faces[mid].BoundingBox, axis);

            // Move pivot to right
            (faces[mid], faces[right]) = (faces[right], faces[mid]);

            int i = left;
            for (int j = left; j < right; j++)
            {
                // Use <= instead of < to avoid lopsided clustering on duplicates
                if (GetCentroidAxis(faces[j].BoundingBox, axis) <= pivotVal)
                {
                    (faces[i], faces[j]) = (faces[j], faces[i]);
                    i++;
                }
            }

            (faces[i], faces[right]) = (faces[right], faces[i]);
            return i;
        }

        private static double GetCentroidAxis(in AABB box, int axis)
        {
            Vector3d c = box.Centroid;
            return axis == 0 ? c.X : (axis == 1 ? c.Y : c.Z);
        }

        // Fast query: iterative traversal with a stack to prevent call-stack overhead
        public static void FindOverlappingFaces(AABBNode root, in AABB targetBox, List<Face> results)
        {
            if (root == null) return;

            // Fixed-size traversal stack (max depth of balanced tree < 64 for billions of faces)
            var stack = new Stack<AABBNode>(32);
            stack.Push(root);

            while (stack.Count > 0)
            {
                var current = stack.Pop();

                if (!current.BoundingBox.Intersects(targetBox))
                    continue;

                if (current.IsLeaf)
                {
                    results.Add(current.Face);
                    continue;
                }

                if (current.Right != null) stack.Push(current.Right);
                if (current.Left != null) stack.Push(current.Left);
            }
        }
    }
}