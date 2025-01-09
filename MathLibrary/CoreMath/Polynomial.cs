namespace MathLibrary
{
 /// <summary>
    /// Represents a polynomial with real coefficients
    /// </summary>
    public class Polynomial
    {
        private readonly double[] _coefficients;  // Ascending order: constant term first
        public int Degree => _coefficients.Length - 1;

        public Polynomial(params double[] coefficients)
        {
            if (coefficients == null || coefficients.Length == 0)
                throw new ArgumentException("Polynomial must have at least one coefficient");
            _coefficients = (double[])coefficients.Clone();
        }

        public double Evaluate(double x)
        {
            // Horner's method for polynomial evaluation
            double result = _coefficients[^1];
            for (int i = _coefficients.Length - 2; i >= 0; i--)
                result = result * x + _coefficients[i];
            return result;
        }

        public Polynomial Derivative()
        {
            if (Degree == 0)
                return new Polynomial(0);

            var derivCoeffs = new double[Degree];
            for (int i = 1; i <= Degree; i++)
                derivCoeffs[i - 1] = i * _coefficients[i];
            return new Polynomial(derivCoeffs);
        }

        public Complex[] FindRoots(double tolerance = 1e-10, int maxIterations = 100)
        {
            // Implement root finding using Durand-Kerner method
            int n = Degree;
            if (n == 0)
                return Array.Empty<Complex>();

            // Initialize with rough guesses for roots
            var roots = new Complex[n];
            for (int i = 0; i < n; i++)
                roots[i] = Complex.FromPolarCoordinates(1, 2 * Math.PI * i / n);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                bool converged = true;
                var newRoots = new Complex[n];

                for (int i = 0; i < n; i++)
                {
                    Complex numerator = EvaluateComplex(roots[i]);
                    Complex denominator = Complex.One;

                    for (int j = 0; j < n; j++)
                    {
                        if (j != i)
                            denominator *= roots[i] - roots[j];
                    }

                    newRoots[i] = roots[i] - numerator / denominator;
                    if (Complex.Abs(newRoots[i] - roots[i]) > tolerance)
                        converged = false;
                }

                roots = newRoots;
                if (converged)
                    break;
            }

            return roots;
        }

        private Complex EvaluateComplex(Complex x)
        {
            Complex result = _coefficients[^1];
            for (int i = _coefficients.Length - 2; i >= 0; i--)
                result = result * x + _coefficients[i];
            return result;
        }
    }
}