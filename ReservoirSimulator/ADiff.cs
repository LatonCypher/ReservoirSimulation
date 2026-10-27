using System;
using System.Runtime.CompilerServices;
using SepalSolver;

namespace ReservoirSimulator
{
    [method: MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly struct ADiff<T>(double value, T derivative) where T : struct, IVectorND<T>
    {
        public readonly double Value = value;
        public readonly T Derivative = derivative;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator ADiff<T>(double a) => new(a, T.Zero);

        // --- Unary & Negation ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(ADiff<T> a)
            => new(-a.Value, -a.Derivative);

        // --- Addition ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator +(ADiff<T> a, ADiff<T> b)
            => new(a.Value + b.Value, a.Derivative + b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator +(double a, ADiff<T> b)
            => new(a + b.Value, b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator +(ADiff<T> a, double b)
            => new(a.Value + b, a.Derivative);

        // --- Subtraction ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(ADiff<T> a, ADiff<T> b)
            => new(a.Value - b.Value, a.Derivative - b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(double a, ADiff<T> b)
            => new(a - b.Value, -b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(ADiff<T> a, double b)
            => new(a.Value - b, a.Derivative);

        // --- Multiplication ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator *(ADiff<T> a, ADiff<T> b)
            => new(a.Value * b.Value, (b.Value * a.Derivative) + (a.Value * b.Derivative));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator *(double a, ADiff<T> b)
            => new(a * b.Value, a * b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator *(ADiff<T> a, double b)
            => new(a.Value * b, a.Derivative * b);

        // --- Division ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator /(ADiff<T> a, ADiff<T> b)
        {
            double invB = 1.0 / b.Value;
            double val = a.Value * invB;
            T deriv = (a.Derivative - (val * b.Derivative)) * invB;
            return new(val, deriv);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator /(double a, ADiff<T> b)
        {
            double invB = 1.0 / b.Value;
            double val = a * invB;
            T deriv = (-val * invB) * b.Derivative;
            return new(val, deriv);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator /(ADiff<T> a, double b)
        {
            double invB = 1.0 / b;
            return new(a.Value * invB, a.Derivative * invB);
        }

        // --- Relational Operators ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator >(ADiff<T> a, ADiff<T> b) => a.Value > b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator <(ADiff<T> a, ADiff<T> b) => a.Value < b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator <=(ADiff<T> a, ADiff<T> b) => a.Value <= b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator >=(ADiff<T> a, ADiff<T> b) => a.Value >= b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator >(ADiff<T> a, double b) => a.Value > b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator <(ADiff<T> a, double b) => a.Value < b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator <=(ADiff<T> a, double b) => a.Value <= b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator >=(ADiff<T> a, double b) => a.Value >= b;
    }
}