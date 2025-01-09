namespace MathLibrary
{
    /// <summary>
    /// Provides numerical integration methods
    /// </summary>
    public static class NumericalIntegration
    {
        public static double SimpsonsRule(Func<double, double> f, double a, double b, int n)
        {
            if (n % 2 != 0)
                throw new ArgumentException("Number of intervals must be even");

            double h = (b - a) / n;
            double sum = f(a) + f(b);

            for (int i = 1; i < n; i++)
            {
                double x = a + i * h;
                sum += f(x) * (i % 2 == 0 ? 2 : 4);
            }

            return sum * h / 3;
        }

        public static double GaussianQuadrature(Func<double, double> f, double a, double b, int n)
        {
            // Implement n-point Gaussian quadrature
            // For simplicity, we'll implement 5-point Gaussian quadrature
            var points = new double[] { -0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664 };
            var weights = new double[] { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 };

            double sum = 0;
            double middle = (b + a) / 2;
            double length = (b - a) / 2;

            for (int i = 0; i < points.Length; i++)
            {
                double x = middle + length * points[i];
                sum += weights[i] * f(x);
            }

            return sum * length;
        }
    }
}