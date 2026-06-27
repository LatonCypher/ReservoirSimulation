using System;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace ReservoirSimulator
{
    internal struct AutoDiff
    {
        public double Value;
        public Vector4d Derivative;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public AutoDiff(double value, Vector4d derivative)
        {
            Value = value;
            Derivative = derivative;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator AutoDiff(double a) => new AutoDiff(a, Vector4d.Zero);

        // Addition (u + v) -> d(u+v) = du + dv
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AutoDiff operator +(AutoDiff a, AutoDiff b)
            => new (a.Value + b.Value, a.Derivative + b.Derivative);

        // Subtraction (u - v) -> d(u-v) = du - dv
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AutoDiff operator -(AutoDiff a, AutoDiff b)
            => new (a.Value - b.Value, a.Derivative - b.Derivative);

        // Multiplication (u * v) -> d(uv) = v*du + u*dv
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AutoDiff operator *(AutoDiff a, AutoDiff b)
            => new (a.Value * b.Value, (b.Value * a.Derivative) + (a.Value * b.Derivative));

        // Division (u / v) -> d(u/v) = (v*du - u*dv) / v^2
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AutoDiff operator /(AutoDiff a, AutoDiff b)
        {
            double invB = 1.0 / b.Value;
            double value = a.Value * invB;
            Vector4d derivative = (b.Value * a.Derivative - a.Value * b.Derivative) * (invB * invB);
            return new (value, derivative);
        }
    }
}
