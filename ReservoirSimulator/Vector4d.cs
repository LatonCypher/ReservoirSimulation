using System;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

namespace ReservoirSimulator
{
    public struct Vector4d
    {
        // The internal 256-bit hardware register containing 4 doubles
        internal Vector256<double> _value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector4d(double x, double y, double z, double w)
        {
            _value = Vector256.Create(x, y, z, w);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal Vector4d(Vector256<double> value)
        {
            _value = value;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector4d SetElement(int index, double value)
        {
            return index switch
            {
                0 => new Vector4d(_value.WithElement(0, value)),
                1 => new Vector4d(_value.WithElement(1, value)),
                2 => new Vector4d(_value.WithElement(2, value)),
                3 => new Vector4d(_value.WithElement(3, value)),
                _ => throw new ArgumentOutOfRangeException(nameof(index))
            };
        }

        // Static properties for rapid initialization
        public static Vector4d Zero
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => new Vector4d(Vector256<double>.Zero);
        }

        // Addition
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector4d operator +(Vector4d a, Vector4d b)
        {
            return new Vector4d(a._value + b._value);
        }

        // Subtraction
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector4d operator -(Vector4d a, Vector4d b)
        {
            return new Vector4d(a._value - b._value);
        }

        // Vector * Vector (Hadamard / Component-wise multiplication)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector4d operator *(Vector4d a, Vector4d b)
        {
            return new Vector4d(a._value * b._value);
        }

        // Scalar * Vector (Crucial for Automatic Differentiation Chain Rule)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector4d operator *(double scalar, Vector4d vector)
        {
            var scalarVector = Vector256.Create(scalar);
            return new Vector4d(scalarVector * vector._value);
        }

        // Scalar * Vector (Crucial for Automatic Differentiation Chain Rule)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector4d operator *(Vector4d vector, double scalar)
        {
            var scalarVector = Vector256.Create(scalar);
            return new Vector4d(scalarVector * vector._value);
        }
    }
}
