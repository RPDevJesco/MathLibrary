namespace MathLibrary
{
    /// <summary>
    /// Provides differential equation solvers
    /// </summary>
    public static class DifferentialEquations
    {
        public static Vector RungeKutta4(
            Func<double, Vector, Vector> f,
            Vector y0,
            double t0,
            double tf,
            int steps)
        {
            double h = (tf - t0) / steps;
            Vector y = y0;
            double t = t0;

            for (int i = 0; i < steps; i++)
            {
                Vector k1 = f(t, y);
                Vector k2 = f(t + h/2, y + k1 * (h/2));
                Vector k3 = f(t + h/2, y + k2 * (h/2));
                Vector k4 = f(t + h, y + k3 * h);

                y += (k1 + k2 * 2 + k3 * 2 + k4) * (h/6);
                t += h;
            }

            return y;
        }
    }
}