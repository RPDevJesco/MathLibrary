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
            int n = Degree;
            if (n == 0)
                return Array.Empty<Complex>();

            // Scale coefficients to improve numerical stability
            double maxCoeff = _coefficients.Max(Math.Abs);
            var scaledCoeffs = _coefficients.Select(c => c / maxCoeff).ToArray();

            // Estimate initial radius using bounds on roots
            double rootBound = 1.0;
            for (int i = 0; i < n; i++)
            {
                rootBound = Math.Max(rootBound, Math.Abs(scaledCoeffs[i] / scaledCoeffs[n]));
            }
            rootBound += 1.0; // Add margin for numerical stability

            // Initialize with better guesses for roots using scaled radius
            var roots = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                double angle = 2 * Math.PI * i / n;
                roots[i] = Complex.FromPolar(rootBound, angle);
            }

            // Use Aberth-Ehrlich method (more stable than Durand-Kerner)
            for (int iter = 0; iter < maxIterations; iter++)
            {
                bool converged = true;
                var newRoots = new Complex[n];

                for (int i = 0; i < n; i++)
                {
                    Complex p = EvaluateComplex(roots[i]);
                    Complex pDeriv = EvaluateComplexDerivative(roots[i]);
                    
                    Complex sum = Complex.Zero;
                    for (int j = 0; j < n; j++)
                    {
                        if (j != i)
                        {
                            sum += 1.0 / (roots[i] - roots[j]);
                        }
                    }

                    Complex correction = p / (pDeriv - p * sum);
                    newRoots[i] = roots[i] - correction;

                    if (Complex.Abs(correction) > tolerance * (1 + Complex.Abs(roots[i])))
                    {
                        converged = false;
                    }
                }

                roots = newRoots;
                if (converged)
                    break;
            }

            // Scale roots back
            for (int i = 0; i < n; i++)
            {
                roots[i] = roots[i] * Math.Pow(maxCoeff, 1.0 / n);
            }

            // Clean up small imaginary parts
            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(roots[i].Imaginary) < tolerance)
                {
                    roots[i] = new Complex(roots[i].Real, 0);
                }
            }

            return roots;
        }

        private Complex EvaluateComplexDerivative(Complex x)
        {
            if (_coefficients.Length <= 1)
                return Complex.Zero;

            Complex result = _coefficients[^1] * (_coefficients.Length - 1);
            for (int i = _coefficients.Length - 2; i >= 1; i--)
                result = result * x + _coefficients[i] * i;
            return result;
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