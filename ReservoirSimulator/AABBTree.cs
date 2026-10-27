using System;
using System.Collections.Generic;

namespace ReservoirSimulator
{
    public class AABBNode
    {
        public AABB BoundingBox { get; set; }
        public AABBNode Left { get; set; }
        public AABBNode Right { get; set; }
        public Face Face { get; set; } // Uses your existing Face class

        public bool IsLeaf => Face != null;
    }

    public struct AABB
    {
        public Vector3d Min { get; set; }
        public Vector3d Max { get; set; }

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

        public static AABB CreateMerged(AABB b1, AABB b2)
        {
            return new AABB(
                Vector3d.Min(b1.Min, b2.Min),
                Vector3d.Max(b1.Max, b2.Max)
            );
        }
    }

    public class AABBTree
    {
        public AABBNode Root { get; private set; }

        public void Build(List<Face> faces)
        {
            if (faces == null || faces.Count == 0) return;
            Root = BuildRecursive(faces);
        }

        private static AABBNode BuildRecursive(List<Face> faces)
        {
            AABBNode node = new AABBNode();

            if (faces.Count == 1)
            {
                node.Face = faces[0];
                node.BoundingBox = faces[0].BoundingBox;
                return node;
            }

            node.BoundingBox = EncapsulateAll(faces);

            Vector3d size = node.BoundingBox.Max - node.BoundingBox.Min;
            int axis = 0; // 0 = X, 1 = Y, 2 = Z
            if (size.Y > size.X) axis = 1;
            if (size.Z > size.Y && size.Z > size.X) axis = 2;

            faces.Sort((a, b) =>
            {
                double centroidA = GetCentroidComponent(a.BoundingBox, axis);
                double centroidB = GetCentroidComponent(b.BoundingBox, axis);
                return centroidA.CompareTo(centroidB);
            });

            int mid = faces.Count / 2;
            List<Face> leftFaces = faces.GetRange(0, mid);
            List<Face> rightFaces = faces.GetRange(mid, faces.Count - mid);

            node.Left = BuildRecursive(leftFaces);
            node.Right = BuildRecursive(rightFaces);

            return node;
        }

        private static AABB EncapsulateAll(List<Face> faces)
        {
            Vector3d min = faces[0].BoundingBox.Min;
            Vector3d max = faces[0].BoundingBox.Max;

            for (int i = 1; i < faces.Count; i++)
            {
                min = Vector3d.Min(min, faces[i].BoundingBox.Min);
                max = Vector3d.Max(max, faces[i].BoundingBox.Max);
            }
            return new AABB(min, max);
        }

        private static double GetCentroidComponent(AABB box, int axis)
        {
            Vector3d centroid = (box.Min + box.Max) * 0.5;
            return axis == 0 ? centroid.X : (axis == 1 ? centroid.Y : centroid.Z);
        }

        public static void FindOverlappingFaces(AABBNode node, AABB targetBox, List<Face> results)
        {
            if (node == null) return;

            if (!Intersects(node.BoundingBox, targetBox)) return;

            if (node.IsLeaf)
            {
                results.Add(node.Face);
                return;
            }

            FindOverlappingFaces(node.Left, targetBox, results);
            FindOverlappingFaces(node.Right, targetBox, results);
        }

        private static bool Intersects(AABB a, AABB b)
        {
            return (a.Min.X <= b.Max.X && a.Max.X >= b.Min.X) &&
                   (a.Min.Y <= b.Max.Y && a.Max.Y >= b.Min.Y) &&
                   (a.Min.Z <= b.Max.Z && a.Max.Z >= b.Min.Z);
        }
    }
}