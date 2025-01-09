namespace MathLibrary
{
    public static class ComplexMath
    {
        public static Complex Sin(Complex z)
        {
            return new Complex(
                Math.Sin(z.Real) * Math.Cosh(z.Imaginary),
                Math.Cos(z.Real) * Math.Sinh(z.Imaginary)
            );
        }

        public static Complex Cos(Complex z)
        {
            return new Complex(
                Math.Cos(z.Real) * Math.Cosh(z.Imaginary),
                -Math.Sin(z.Real) * Math.Sinh(z.Imaginary)
            );
        }

        public static Complex Tan(Complex z)
        {
            Complex sin = Sin(z);
            Complex cos = Cos(z);
            return sin / cos;
        }

        public static Complex Exp(Complex z)
        {
            double expReal = Math.Exp(z.Real);
            return new Complex(
                expReal * Math.Cos(z.Imaginary),
                expReal * Math.Sin(z.Imaginary)
            );
        }

        public static Complex Log(Complex z)
        {
            return new Complex(
                Math.Log(z.Magnitude),
                z.Phase
            );
        }

        public static Complex Sqrt(Complex z)
        {
            double r = Math.Sqrt(z.Magnitude);
            double theta = z.Phase / 2;
            return new Complex(
                r * Math.Cos(theta),
                r * Math.Sin(theta)
            );
        }

        public static Complex Pow(Complex z, Complex w)
        {
            if (z.Real == 0 && z.Imaginary == 0)
            {
                if (w.Real == 0 && w.Imaginary == 0)
                    return new Complex(1, 0);
                return new Complex(0, 0);
            }

            return Exp(w * Log(z));
        }
    }
}