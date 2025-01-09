namespace MathLibrary
{
    public static class SpecialFunctions
    {
        private const double Epsilon = 1e-15;
        private const int MaxIterations = 1000;

        #region Gamma Function and Related
        public static double Gamma(double x)
        {
            // Handle special cases
            if (double.IsNaN(x)) return double.NaN;
            if (double.IsInfinity(x)) return double.PositiveInfinity;
            if (x == 0) return double.PositiveInfinity;
            if (x < 0 && Math.Floor(x) == x) return double.NaN; // Negative integers

            // For negative numbers, use the reflection formula
            if (x < 0)
            {
                return Math.PI / (Math.Sin(Math.PI * x) * Gamma(-x));
            }

            // Lanczos approximation for x > 0
            double[] p = {
                676.5203681218851,
                -1259.1392167224028,
                771.32342877765313,
                -176.61502916214059,
                12.507343278686905,
                -0.13857109526572012,
                9.9843695780195716e-6,
                1.5056327351493116e-7
            };

            x -= 1;
            double t = x + 7.5;
            double sum = 0.99999999999980993;

            for (int i = 0; i < p.Length; i++)
            {
                sum += p[i] / (x + i + 1);
            }

            return Math.Sqrt(2 * Math.PI) * Math.Pow(t, x + 0.5) * Math.Exp(-t) * sum;
        }

        public static double LogGamma(double x)
        {
            if (x <= 0) throw new ArgumentException("LogGamma undefined for x <= 0");

            double temp = (x - 0.5) * Math.Log(x + 4.5) - (x + 4.5);
            double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                        + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                        + 0.00120858003 / (x + 4) - 0.00000536382 / (x + 5);
            
            return temp + Math.Log(ser * Math.Sqrt(2 * Math.PI));
        }

        public static double FactorialInt(int n)
        {
            if (n < 0) throw new ArgumentException("Factorial undefined for negative numbers");
            if (n > 170) throw new ArgumentException("Factorial too large to compute");

            double result = 1.0;
            for (int i = 2; i <= n; i++)
                result *= i;
            return result;
        }
        #endregion

        #region Bessel Functions
        public static double BesselJ0(double x)
        {
            if (Math.Abs(x) < 8.0)
            {
                // Use series expansion for small x
                double y = x * x;
                return 1.0 + y * (-0.25 + y * (0.015625 + y * (-0.000434028 + y * 
                       (6.78168e-6 + y * (-6.78168e-8 + y * 4.7095e-10)))));
            }
            else
            {
                // Use asymptotic approximation for large x
                double ax = Math.Abs(x);
                double z = 8.0 / ax;
                double y = z * z;
                double xx = ax - 0.785398164;
                
                double p0 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + 
                           y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                double q0 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * 
                           (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));
                
                return Math.Sqrt(2.0 / (Math.PI * ax)) * (p0 * Math.Cos(xx) - z * q0 * Math.Sin(xx));
            }
        }

        public static double BesselJ1(double x)
        {
            if (Math.Abs(x) < 8.0)
            {
                // Use series expansion for small x
                double y = x * x;
                return x * (0.5 + y * (-0.0625 + y * (0.00260417 + y * (-6.51042e-5 + 
                       y * (1.07512e-6 + y * (-1.32692e-8))))));
            }
            else
            {
                // Use asymptotic approximation for large x
                double ax = Math.Abs(x);
                double z = 8.0 / ax;
                double y = z * z;
                double xx = ax - 2.356194491;
                
                double p1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 + 
                           y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                double q1 = 0.04687499995 + y * (-0.2002690873e-3 + y * 
                           (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
                
                double ans = Math.Sqrt(2.0 / (Math.PI * ax)) * (p1 * Math.Cos(xx) - z * q1 * Math.Sin(xx));
                return x < 0.0 ? -ans : ans;
            }
        }

        public static double BesselY0(double x)
        {
            if (x <= 0) throw new ArgumentException("BesselY0 undefined for x <= 0");

            if (x < 8.0)
            {
                // Use series expansion and relation to J0
                double y = x * x;
                double ans1 = -2.957821389 + y * (0.7062834065 + y * (-0.0517656398 + 
                             y * (0.00074348057 - y * 2.0092e-5)));
                double ans2 = Math.Log(x/2.0) * BesselJ0(x);
                return ans1 + ans2;
            }
            else
            {
                // Use asymptotic approximation
                double z = 8.0/x;
                double y = z * z;
                double xx = x - 0.785398164;
                
                double p0 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + 
                           y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                double q0 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * 
                           (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));
                
                return Math.Sqrt(2.0/(Math.PI * x)) * (p0 * Math.Sin(xx) + z * q0 * Math.Cos(xx));
            }
        }
        #endregion

        #region Legendre Functions
        public static double LegendreP(int n, double x)
        {
            if (n < 0) throw new ArgumentException("Order must be non-negative");
            if (Math.Abs(x) > 1) throw new ArgumentException("Argument must be in [-1,1]");

            if (n == 0) return 1;
            if (n == 1) return x;

            // Use recurrence relation
            double pnm1 = 1;    // P_(n-1)
            double pn = x;      // P_n
            double pnp1;        // P_(n+1)

            for (int k = 2; k <= n; k++)
            {
                pnp1 = ((2 * k - 1) * x * pn - (k - 1) * pnm1) / k;
                pnm1 = pn;
                pn = pnp1;
            }

            return pn;
        }

        public static double AssociatedLegendreP(int l, int m, double x)
        {
            if (l < 0) throw new ArgumentException("Degree l must be non-negative");
            if (Math.Abs(m) > l) throw new ArgumentException("Order m must satisfy |m| <= l");
            if (Math.Abs(x) > 1) throw new ArgumentException("Argument must be in [-1,1]");

            // Handle negative m using relation P(l,-m) = (-1)^m * (l-m)!/(l+m)! * P(l,m)
            if (m < 0)
            {
                m = -m;
                double factor = (m % 2 == 0 ? 1 : -1);
                for (int i = l - m + 1; i <= l + m; i++)
                    factor /= i;
                return factor * AssociatedLegendreP(l, m, x);
            }

            // Compute P(m,m)
            double pmm = 1.0;
            if (m > 0)
            {
                double somx2 = Math.Sqrt((1 - x) * (1 + x));
                double fact = 1.0;
                for (int i = 1; i <= m; i++)
                {
                    pmm *= -fact * somx2;
                    fact += 2.0;
                }
            }

            if (l == m) return pmm;

            // Compute P(m+1,m)
            double pmmp1 = x * (2 * m + 1) * pmm;
            if (l == m + 1) return pmmp1;

            // Use recurrence relation
            double pll = 0.0;
            for (int ll = m + 2; ll <= l; ll++)
            {
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }

            return pll;
        }
        #endregion

        #region Error Functions
        public static double Erf(double x)
        {
            if (double.IsNaN(x)) return double.NaN;
            if (double.IsPositiveInfinity(x)) return 1.0;
            if (double.IsNegativeInfinity(x)) return -1.0;

            // Constants
            const double a1 = 0.254829592;
            const double a2 = -0.284496736;
            const double a3 = 1.421413741;
            const double a4 = -1.453152027;
            const double a5 = 1.061405429;
            const double p = 0.3275911;

            // Save the sign of x
            int sign = Math.Sign(x);
            x = Math.Abs(x);

            // A&S formula 7.1.26
            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

            return sign * y;
        }

        public static double Erfc(double x)
        {
            return 1 - Erf(x);
        }
        #endregion

        #region Hypergeometric Functions
        public static double Hypergeometric2F1(double a, double b, double c, double z)
        {
            if (Math.Abs(z) >= 1.0)
                throw new ArgumentException("|z| must be less than 1");

            if (c <= 0)
                throw new ArgumentException("c must be positive");

            if (c - a - b <= 0)
                throw new ArgumentException("c - a - b must be positive for convergence");

            double sum = 1.0;
            double term = 1.0;
            int n = 0;

            while (Math.Abs(term) > Epsilon && n < MaxIterations)
            {
                term *= (a + n) * (b + n) * z / ((c + n) * (n + 1));
                sum += term;
                n++;

                if (double.IsInfinity(sum))
                    throw new OverflowException("Series diverged");
            }

            if (n >= MaxIterations)
                throw new InvalidOperationException("Failed to converge");

            return sum;
        }
        #endregion

        #region Helper Functions
        private static double PochhammerSymbol(double x, int n)
        {
            double result = 1;
            for (int i = 0; i < n; i++)
                result *= (x + i);
            return result;
        }

        private static double BinomialCoefficient(int n, int k)
        {
            if (k < 0 || k > n) return 0;
            if (k == 0 || k == n) return 1;
            k = Math.Min(k, n - k);
            
            double result = 1;
            for (int i = 0; i < k; i++)
            {
                result *= (n - i);
                result /= (i + 1);
            }
            return result;
        }
        #endregion
    }
}