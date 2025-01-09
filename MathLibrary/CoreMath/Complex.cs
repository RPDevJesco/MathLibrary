namespace MathLibrary
{
    /// <summary>
    /// Represents a complex number with real and imaginary parts
    /// </summary>
    public struct Complex
    {
        public double Real { get; }
        public double Imaginary { get; }

        public Complex(double real, double imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        // Implicit conversion from double to Complex
        public static implicit operator Complex(double d) => new(d, 0);

        public static Complex FromPolar(double magnitude, double phase)
        {
            return new Complex(
                magnitude * Math.Cos(phase),
                magnitude * Math.Sin(phase)
            );
        }

        public static Complex FromPolarCoordinates(double magnitude, double phase)
        {
            return FromPolar(magnitude, phase);
        }

        public static Complex operator +(Complex a, Complex b)
            => new(a.Real + b.Real, a.Imaginary + b.Imaginary);

        public static Complex operator -(Complex a, Complex b)
            => new(a.Real - b.Real, a.Imaginary - b.Imaginary);

        public static Complex operator *(Complex a, Complex b)
            => new(
                a.Real * b.Real - a.Imaginary * b.Imaginary,
                a.Real * b.Imaginary + a.Imaginary * b.Real
            );

        public static Complex operator /(Complex a, Complex b)
        {
            double denominator = b.Real * b.Real + b.Imaginary * b.Imaginary;
            if (denominator == 0)
                throw new DivideByZeroException("Division by zero in complex number");
            
            return new Complex(
                (a.Real * b.Real + a.Imaginary * b.Imaginary) / denominator,
                (a.Imaginary * b.Real - a.Real * b.Imaginary) / denominator
            );
        }

        public static Complex operator /(Complex a, double b)
        {
            if (b == 0)
                throw new DivideByZeroException("Division by zero");
            return new Complex(a.Real / b, a.Imaginary / b);
        }

        public static Complex operator +(Complex a, double b)
            => new(a.Real + b, a.Imaginary);

        public static Complex operator +(double a, Complex b)
            => new(a + b.Real, b.Imaginary);

        public static Complex operator -(Complex a, double b)
            => new(a.Real - b, a.Imaginary);

        public static Complex operator -(double a, Complex b)
            => new(a - b.Real, -b.Imaginary);

        public static Complex operator *(Complex a, double b)
            => new(a.Real * b, a.Imaginary * b);

        public static Complex operator *(double a, Complex b)
            => new(b.Real * a, b.Imaginary * a);

        public static double Abs(Complex c)
            => Math.Sqrt(c.Real * c.Real + c.Imaginary * c.Imaginary);

        public double Magnitude
            => Abs(this);

        public double Phase
            => Math.Atan2(Imaginary, Real);

        public Complex Conjugate
            => new(Real, -Imaginary);

        public static readonly Complex Zero = new(0, 0);
        public static readonly Complex One = new(1, 0);
        public static readonly Complex ImaginaryOne = new(0, 1);

        public override string ToString()
            => Imaginary >= 0 
                ? $"{Real} + {Imaginary}i" 
                : $"{Real} - {Math.Abs(Imaginary)}i";
    }
}