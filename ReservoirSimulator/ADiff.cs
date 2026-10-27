using System.Runtime.CompilerServices;
using SepalSolver;

namespace ReservoirSimulator
{
    [method: MethodImpl(MethodImplOptions.AggressiveInlining)]
    public struct ADiff<T>(double value, T derivative) where T : struct, IVectorND<T>
    {
        public double Value = value;
        public T Derivative = derivative;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator ADiff<T>(double a) => new(a, T.Zero);

        // Negation (-u) -> d(-u) = -du
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(ADiff<T> a)
            => new(-a.Value, -a.Derivative);

        // Addition (u + v) -> d(u+v) = du + dv
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator +(ADiff<T> a, ADiff<T> b)
            => new(a.Value + b.Value, a.Derivative + b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator +(double a, ADiff<T> b)
            => new(a + b.Value, b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator +(ADiff<T> a, double b)
            => new(a.Value + b, a.Derivative);

        // Subtraction (u - v) -> d(u-v) = du - dv
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(ADiff<T> a, ADiff<T> b)
            => new(a.Value - b.Value, a.Derivative - b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(double a, ADiff<T> b)
            => new(a - b.Value, -b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator -(ADiff<T> a, double b)
            => new(a.Value - b, a.Derivative);

        // Multiplication (u * v) -> d(uv) = v*du + u*dv
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator *(ADiff<T> a, ADiff<T> b)
            => new(a.Value * b.Value, (b.Value * a.Derivative) + (a.Value * b.Derivative));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator *(double a, ADiff<T> b)
            => new(a * b.Value, a * b.Derivative);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator *(ADiff<T> a, double b)
            => new(a.Value * b, a.Derivative * b);

        // Division (u / v) -> d(u/v) = (v*du - u*dv) / v^2
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator /(ADiff<T> a, ADiff<T> b)
        {
            double invB = 1.0 / b.Value;
            double value = a.Value * invB;
            T derivative = (b.Value * a.Derivative - a.Value * b.Derivative) * (invB * invB);
            return new(value, derivative);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator /(double a, ADiff<T> b)
        {
            double invB = 1.0 / b.Value;
            double value = a * invB;
            T derivative = (-a * invB * invB) * b.Derivative;
            return new(value, derivative);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> operator /(ADiff<T> a, double b)
        {
            double invB = 1.0 / b;
            double value = a.Value * invB;
            T derivative = a.Derivative * invB;
            return new(value, derivative);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator >(ADiff<T> a, ADiff<T> b) 
            => a.Value > b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator <(ADiff<T> a, ADiff<T> b)
            => a.Value < b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator <=(ADiff<T> a, ADiff<T> b)
            => a.Value <= b.Value;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator >=(ADiff<T> a, ADiff<T> b)
            => a.Value >= b.Value;
    }
}