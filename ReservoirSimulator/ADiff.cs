using System;
using System.Collections.Generic;

internal class ADiff
{
    public static int capacity;
    public double Value;
    public Dictionary<int, double> Derivatives;

    /// <summary>
    /// Resets the instance state so it can be reused like a fresh scratchpad.
    /// </summary>
    public ADiff()
    {
        Value = 0.0;
        Derivatives = new(capacity);
    }
    public ADiff(double value)
    {
        Value = value;
        Derivatives = new(capacity);
    }
    public ADiff(double value, int index)
    {
        Value = value;
        Derivatives = new(capacity)
        { [index] = 1.0 };
    }
    public ADiff Clear()
    {
        Value = 0.0;
        Derivatives.Clear();
        return this;
    }
    public ADiff CopyFrom(ADiff source)
    {
        Value = source.Value;
        Derivatives.Clear();
        foreach (var kvp in source.Derivatives)
        {
            Derivatives[kvp.Key] = kvp.Value;
        }
        return this;
    }
    public static implicit operator ADiff(double a) => new(a);
    public ADiff AddInPlace(ADiff other)
    {
        if(other.Value == 0)
            other.Value = 0;
        Value += other.Value;
        // d(u + v) = du + dv
        foreach (var kvp in other.Derivatives)
        {
            if (Derivatives.TryGetValue(kvp.Key, out double existingDerivative))
                Derivatives[kvp.Key] = existingDerivative + kvp.Value;
            else
                Derivatives[kvp.Key] = kvp.Value;
        }
        return this;
    }
    public ADiff AddInPlace(double other)
    {
        Value += other;
        return this;
    }
    public ADiff AddProductInPlace(ADiff source, double scalar)
    {
        Value += source.Value * scalar;
        //d(u + c*v) = du + c*dv
        foreach (var kvp in source.Derivatives)
        {
            if (Derivatives.TryGetValue(kvp.Key, out double du))
                Derivatives[kvp.Key] = du + (scalar * kvp.Value);
            else
                Derivatives[kvp.Key] = (scalar * kvp.Value);
        }
        return this;
    }
    public ADiff SubtractInPlace(ADiff other)
    {
        Value -= other.Value;
        // d(u - v) = du - dv
        foreach (var kvp in other.Derivatives)
        {
            if (Derivatives.TryGetValue(kvp.Key, out double du))
                Derivatives[kvp.Key] = du - kvp.Value;
            else
                Derivatives[kvp.Key] = -kvp.Value;
        }
        return this;
    }
    public ADiff SubtractInPlace(double other)
    {
        Value -= other;
        return this;
    }
    public ADiff SubtractProductInPlace(ADiff source, double scalar)
    {
        Value -= source.Value * scalar;
        //d(u - c*v) = du - c*dv
        foreach (var kvp in source.Derivatives)
        {
            if (Derivatives.TryGetValue(kvp.Key, out double du))
                Derivatives[kvp.Key] = du - (scalar * kvp.Value);
            else
                Derivatives[kvp.Key] = -(scalar * kvp.Value);        }

        return this;
    }
    public ADiff MultiplyInPlace(ADiff other)
    {
        double u = Value, v = other.Value;

        Value = u * v;
        // d(u * v) = v*du + u*dv
        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= v;

        // Second, add the contribution from 'other' (u * dv)
        foreach (var kvp in other.Derivatives)
        {
            double udv = u * kvp.Value;
            if (Derivatives.TryGetValue(kvp.Key, out double existingvdu))
                Derivatives[kvp.Key] = existingvdu + udv;
            else
                Derivatives[kvp.Key] = udv;
        }
        return this;
    }

    public ADiff MultiplyInPlace(double other)
    {
        Value *= other;
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= other;
        return this;
    }
    public ADiff DivideInPlace(ADiff other)
    {
        double u = Value, v = other.Value;
        if (v == 0.0)
            throw new DivideByZeroException("Cannot divide by an ADiff with a value of zero.");

        Value = u / v;
        double vSquared = v * v;

        foreach (var key in Derivatives.Keys)
            Derivatives[key] /= v;

        foreach (var kvp in other.Derivatives)
        {
            double udvOverVSquared = (u * kvp.Value) / vSquared;
            if (Derivatives.TryGetValue(kvp.Key, out double existingTerm))
                Derivatives[kvp.Key] = existingTerm - udvOverVSquared;
            else
                Derivatives[kvp.Key] = -udvOverVSquared;
        }
        return this;
    }
    public ADiff DivideInPlace(double other)
    {
        if (other == 0.0)
            throw new DivideByZeroException("Cannot divide by an ADiff with a value of zero.");

        Value /= other;

        foreach (var key in Derivatives.Keys)
            Derivatives[key] /= other;
        return this;
    }
    public ADiff NegateInPlace()
    {
        Value = -Value;
        foreach (var key in Derivatives.Keys)
            Derivatives[key] = -Derivatives[key];
        return this;
    }   
    public ADiff ReciprocalInPlace()
    {
        double u = Value;
        if (u == 0.0) throw new DivideByZeroException();
        Value = 1.0 / u;
        double negativeReciprocalSquared = -1.0 / (u * u);

        foreach (var key in Derivatives.Keys)
            Derivatives[key] *= negativeReciprocalSquared;
        
        return this;
    }
    public ADiff LogInPlace()
    {
        double u = Value;
        if (u <= 0.0) 
            throw new ArgumentException("Cannot take logarithm of a non-positive value.");
        Value = Math.Log(u);
        double reciprocal = 1.0 / u;

        foreach (var key in Derivatives.Keys)
            Derivatives[key] *= reciprocal;
    
        return this;
    }
    public ADiff SqrtInPlace()
    {
        double u = Value;
        if (u < 0.0)
            throw new ArgumentException("Cannot take square root of an ADiff with a negative value.");
        Value = Math.Sqrt(u);
        // d(sqrt(u)) = (1/(2*sqrt(u))) * du
        double reciprocalOfTwoSqrt = 0.5 / Value; // This is equivalent to 1/(2*sqrt(u))
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= reciprocalOfTwoSqrt;
        return this;
    }
    public ADiff CbrtInPlace()
    {
        double u = Value;
        Value = Math.Cbrt(u);
        // d(cbrt(u)) = (1/(3*cbrt(u)^2)) * du
        double reciprocalOfThreeCbrtSquared = 1.0 / (3.0 * Value * Value);
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= reciprocalOfThreeCbrtSquared;
        return this;
    }
    public ADiff AsinInPlace()
    {
        double u = Value;
        if (u < -1.0 || u > 1.0)
            throw new ArgumentException("Input to arcsin must be in the range [-1, 1].");
        Value = Math.Asin(u);
        // d(asin(u)) = 1/sqrt(1-u^2) * du
        double reciprocalOfSqrt = 1.0 / Math.Sqrt(1.0 - (u * u));
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= reciprocalOfSqrt;
        return this;
    }
    public ADiff SinInPlace()
    {
        double u = Value;
        var (sinU, cosU) = Math.SinCos(u);
        Value = sinU;
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= cosU;
        return this;
    }
    public ADiff AcosInPlace()
    {
        double u = Value;
        if (u < -1.0 || u > 1.0)
            throw new ArgumentException("Input to arccos must be in the range [-1, 1].");
        Value = Math.Acos(u);
        // d(acos(u)) = -1/sqrt(1-u^2) * du
        double negativeReciprocalOfSqrt = -1.0 / Math.Sqrt(1.0 - (u * u));
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= negativeReciprocalOfSqrt;
        return this;
    }
    public ADiff CosInPlace()
    {
        double u = Value;
        var (sinU, cosU) = Math.SinCos(u);
        Value = cosU;
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= -sinU;
        return this;
    }
    public ADiff AtanInPlace()
    {
        double u = Value;
        Value = Math.Atan(u);
        // d(atan(u)) = 1/(1+u^2) * du
        double reciprocalOfOnePlusUSquared = 1.0 / (1.0 + (u * u));
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= reciprocalOfOnePlusUSquared;
        return this;
    }
    public ADiff TanInPlace()
    {
        double u = Value;
        double tanU = Math.Tan(u);
        Value = tanU;
        // d(tan(u)) = sec^2(u) * du = (1 + tan^2(u)) * du
        double secSquared = 1.0 + (tanU * tanU);
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= secSquared;
        return this;
    }
    public ADiff AsinhInPlace()
    {
        double u = Value;
        Value = Math.Asinh(u);
        double f = 1.0 / Math.Sqrt(u * u + 1.0);
        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= f;
        return this;
    }
    public ADiff SinhInPlace()
    {
        double u = Value, coshU = Math.Cosh(u);
        Value  = Math.Sinh(u);

        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= coshU;
        return this;
    }
    public ADiff AcoshInPlace()
    {
        double u = Value;
        Value = Math.Acosh(u);
        double f = 1.0 / Math.Sqrt(u * u - 1.0);
        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= f;
        return this;
    }
    public ADiff CoshInPlace()
    {
        double u = Value, sinhU = Math.Sinh(u);
        Value  = Math.Cosh(u);

        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= sinhU;
        return this;
    }
    public ADiff AtanhInPlace()
    {
        double u = Value;
        Value = Math.Atanh(u);
        double f = 1.0 / (1.0 - u * u);
        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= f;
        return this;
    }
    public ADiff TanhInPlace()
    {
        double u = Value, sech2U = 1 - u * u;
        Value  = Math.Cosh(u);
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= sech2U;
        return this;
    }
    public ADiff ExpInPlace()
    {
        // Cache the old value because we need it for the derivative calculation
        double u = Value;

        // Update the value: e^u
        Value = Math.Exp(u);

        // We can modify values directly by iterating over keys safely 
        // because we aren't adding or removing keys, just changing their values.

        double scale = Value;
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= scale;

        return this;
    }
    public ADiff PowInPlace(double power)
    {
        // Handle edge case for power of 0 (any number to power of 0 is 1, derivatives become 0)
        if (power == 0.0)
        {
            Value = 1.0;
            Derivatives.Clear();
            return this;
        }

        double u = Value;

        // Update the value: u^power
        Value = Math.Pow(u, power);

        // Chain Rule: d(u^p) = p * u^(p-1) * du
        // Optimization: p * u^(p-1) is exactly equal to: (p * u^p) / u = (power * this.Value) / u
        double scale = (power * Value) / u;

        // Scale all existing derivatives by the chain rule factor
        
        foreach (int key in Derivatives.Keys)
            Derivatives[key] *= scale;
        return this;
    }

    // Operators overloads for syntactic sugar
    public static ADiff operator +(ADiff a, ADiff b) => new ADiff().AddInPlace(a).AddInPlace(b);
    public static ADiff operator -(ADiff a, ADiff b) => new ADiff().AddInPlace(a).SubtractInPlace(b);
    public static ADiff operator *(ADiff a, ADiff b) => new ADiff().AddInPlace(a).MultiplyInPlace(b);
    public static ADiff operator /(ADiff a, ADiff b) => new ADiff().AddInPlace(a).DivideInPlace(b);
    public static ADiff operator -(ADiff a) => new ADiff().AddInPlace(a).NegateInPlace();
    public static ADiff operator +(ADiff a, double b) => new ADiff().AddInPlace(a).AddInPlace(b);
    public static ADiff operator -(ADiff a, double b) => new ADiff().AddInPlace(a).SubtractInPlace(b);
    public static ADiff operator *(ADiff a, double b) => new ADiff().AddInPlace(a).MultiplyInPlace(b); 
    public static ADiff operator *(double b, ADiff a) => new ADiff(b).MultiplyInPlace(a);
    public static ADiff operator /(ADiff a, double b) => new ADiff().AddInPlace(a).DivideInPlace(b);
    public static ADiff operator /(double b, ADiff a) => new ADiff(b).DivideInPlace(a);
    public static bool operator >(ADiff a, ADiff b) => a.Value > b.Value;
    public static bool operator <(ADiff a, ADiff b) => a.Value < b.Value;
    public static bool operator >=(ADiff a, ADiff b) => a.Value >= b.Value;
    public static bool operator <=(ADiff a, ADiff b) => a.Value <= b.Value;
}