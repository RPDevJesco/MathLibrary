using MathLibrary;

class Program
{
    static void Main()
    {
        Console.WriteLine("MathLibrary Showcase\n");

        ShowBasicCalculator();
        ShowSpecialFunctions();
        ShowVectorOperations();
        ShowMatrixOperations();
        ShowComplexNumbers();
        ShowNumericalIntegration();
        ShowDifferentialEquations();
        ShowPolynomialOperations();
        ShowQuaternions();
        ShowStatistics();
    }

    static void ShowBasicCalculator()
    {
        Console.WriteLine("=== Basic Calculator Operations ===");
        Calculator calc = new Calculator();
        
        var expressions = new[]
        {
            "8 / 2(2 + 2)",
            "3 + 5 * 2",           // Basic arithmetic
            "(3 + 5) * 2",         // Parentheses
            "2^3",                 // Exponentiation
            "10 / 2 + 6",          // Division and addition
            "sin(pi/2)",           // Trigonometric function with constant
            "log(1000)",           // Logarithm
            "sqrt(16)",            // Square root
            "abs(-10)",            // Absolute value
            "-3 + 5",              // Unary minus
            "4^(1/2)",            // Square root using exponentiation
            "e^2",                 // Natural exponential
            "round(3.14159, 2)"    // Rounding to decimal places
        };

        foreach (var expr in expressions)
        {
            try
            {
                Console.WriteLine($"{expr} = {calc.Evaluate(expr)}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"{expr} -> Error: {ex.Message}");
            }
        }
        Console.WriteLine();
    }

    static void ShowSpecialFunctions()
    {
        Console.WriteLine("=== Special Functions ===");
        Calculator calc = new Calculator();

        var specialFunctions = new[]
        {
            "gamma(5)",            // Gamma function
            "besselj0(1)",         // Bessel function J0
            "besselj1(1)",         // Bessel function J1
            "legendre(2, 0.5)",    // Legendre polynomial
            "erf(1)",              // Error function
            "hypergeometric2f1(0.5, 0.75, 2.0, 0.5)"  // Hypergeometric with convergent parameters
        };

        foreach (var func in specialFunctions)
        {
            try
            {
                Console.WriteLine($"{func} = {calc.Evaluate(func)}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"{func} -> Error: {ex.Message}");
            }
        }
        Console.WriteLine();
    }

    static void ShowVectorOperations()
    {
        Console.WriteLine("=== Vector Operations ===");
        
        // Create some test vectors
        var v1 = new Vector(1, 2, 3);
        var v2 = new Vector(4, 5, 6);

        Console.WriteLine($"v1 = {v1}");
        Console.WriteLine($"v2 = {v2}");
        Console.WriteLine($"v1 + v2 = {v1 + v2}");
        Console.WriteLine($"v1 · v2 = {Vector.DotProduct(v1, v2)}");
        Console.WriteLine($"v1 × v2 = {Vector.CrossProduct(v1, v2)}");
        Console.WriteLine($"|v1| = {v1.Magnitude}");
        Console.WriteLine($"normalize(v1) = {v1.Normalize()}");
        Console.WriteLine($"angle(v1, v2) = {Vector.Angle(v1, v2)} radians");

        // Gram-Schmidt orthogonalization
        var vectors = new[] { new Vector(1, 1, 0), new Vector(1, 0, 1), new Vector(0, 1, 1) };
        var orthogonal = Vector.GramSchmidt(vectors);
        Console.WriteLine("\nGram-Schmidt orthogonalization:");
        for (int i = 0; i < orthogonal.Length; i++)
            Console.WriteLine($"e{i + 1} = {orthogonal[i]}");
        Console.WriteLine();
    }

    static void ShowMatrixOperations()
    {
        Console.WriteLine("=== Matrix Operations ===");

        // Create test matrices
        var A = new Matrix(new double[,] {
            { 1, 2, 3 },
            { 4, 5, 6 },
            { 7, 8, 9 }
        });
        var B = new Matrix(new double[,] {
            { 9, 8, 7 },
            { 6, 5, 4 },
            { 3, 2, 1 }
        });

        Console.WriteLine("Matrix A:");
        PrintMatrix(A);
        Console.WriteLine("\nMatrix B:");
        PrintMatrix(B);

        Console.WriteLine("\nA + B:");
        PrintMatrix(A + B);

        Console.WriteLine("\nA * B:");
        PrintMatrix(A * B);

        Console.WriteLine("\nTranspose of A:");
        PrintMatrix(A.Transpose());

        // Matrix decompositions
        Console.WriteLine("\nLU Decomposition of A:");
        var (L, U, P, exchanges) = MatrixDecomposition.LUDecomposition(A);
        Console.WriteLine("L matrix:");
        PrintMatrix(L);
        Console.WriteLine("U matrix:");
        PrintMatrix(U);

        // QR Decomposition
        Console.WriteLine("\nQR Decomposition of A:");
        var (Q, R) = MatrixDecomposition.QRDecomposition(A);
        Console.WriteLine("Q matrix:");
        PrintMatrix(Q);
        Console.WriteLine("R matrix:");
        PrintMatrix(R);

        Console.WriteLine();
    }

    static void ShowComplexNumbers()
    {
        Console.WriteLine("=== Complex Number Operations ===");
        
        var z1 = new Complex(3, 4);
        var z2 = new Complex(1, 2);

        Console.WriteLine($"z1 = {z1}");
        Console.WriteLine($"z2 = {z2}");
        Console.WriteLine($"z1 + z2 = {z1 + z2}");
        Console.WriteLine($"z1 * z2 = {z1 * z2}");
        Console.WriteLine($"z1 / z2 = {z1 / z2}");
        Console.WriteLine($"|z1| = {Complex.Abs(z1)}");
        Console.WriteLine($"z1 conjugate = {z1.Conjugate}");
        Console.WriteLine($"Phase of z1 = {z1.Phase} radians");

        // Complex functions
        Console.WriteLine($"exp(i*π) = {ComplexMath.Exp(new Complex(0, Math.PI))}");
        Console.WriteLine();
    }

    static void ShowNumericalIntegration()
    {
        Console.WriteLine("=== Numerical Integration ===");

        // Test function: f(x) = x^2
        Func<double, double> f = x => x * x;

        // Integrate from 0 to 1 (exact result should be 1/3)
        double simpson = NumericalIntegration.SimpsonsRule(f, 0, 1, 100);
        double gauss = NumericalIntegration.GaussianQuadrature(f, 0, 1, 5);

        Console.WriteLine($"∫x² dx from 0 to 1:");
        Console.WriteLine($"Simpson's Rule: {simpson:F6}");
        Console.WriteLine($"Gaussian Quadrature: {gauss:F6}");
        Console.WriteLine($"Exact value: {1.0/3:F6}");
        Console.WriteLine();
    }

    static void ShowDifferentialEquations()
    {
        Console.WriteLine("=== Differential Equations ===");

        // Solve dy/dt = -y, y(0) = 1
        // Exact solution: y = e^(-t)
        Vector y0 = new Vector(1.0);
        double t0 = 0, tf = 1;
        int steps = 100;

        Vector yfinal = DifferentialEquations.RungeKutta4(
            (t, y) => new Vector(-y[0]),
            y0, t0, tf, steps);

        Console.WriteLine($"Solution to dy/dt = -y at t=1:");
        Console.WriteLine($"Numerical: {yfinal[0]:F6}");
        Console.WriteLine($"Exact: {Math.Exp(-1):F6}");
        Console.WriteLine();
    }

    static void ShowPolynomialOperations()
    {
        Console.WriteLine("=== Polynomial Operations ===");

        // Create polynomial: x^3 - 2x^2 + x - 2
        var poly = new Polynomial(1, -2, 1, -2);
        
        Console.WriteLine("Polynomial: x³ - 2x² + x - 2");
        Console.WriteLine($"Value at x=2: {poly.Evaluate(2)}");

        Console.WriteLine("\nRoots:");
        var roots = poly.FindRoots();
        foreach (var root in roots)
        {
            Console.WriteLine($"  {root}");
        }

        var derivative = poly.Derivative();
        Console.WriteLine($"\nDerivative at x=2: {derivative.Evaluate(2)}");
        Console.WriteLine();
    }

    static void ShowQuaternions()
    {
        Console.WriteLine("=== Quaternion Operations ===");

        // Create quaternion for 90-degree rotation around Y axis
        var q = Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI/2);
        var v = new Vector(1, 0, 0);  // vector pointing along X axis

        Console.WriteLine("Rotating vector (1,0,0) by 90° around Y axis:");
        var rotated = q.Rotate(v);
        Console.WriteLine($"Result: {rotated}");  // Should be approximately (0,0,-1)

        // Show Euler angle conversion
        var qEuler = Quaternion.FromEulerAngles(Math.PI/4, 0, 0);  // 45° rotation around X
        Console.WriteLine($"\nQuaternion from 45° X rotation: {qEuler}");
        Console.WriteLine();
    }

    static void ShowStatistics()
    {
        Console.WriteLine("=== Statistical Functions ===");

        var data = new double[] { 2, 4, 4, 4, 5, 5, 7, 9 };
        Console.WriteLine($"Data: {string.Join(", ", data)}");
        Console.WriteLine($"Mean: {Statistics.Mean(data):F2}");
        Console.WriteLine($"Variance: {Statistics.Variance(data):F2}");
        Console.WriteLine($"Standard Deviation: {Statistics.StandardDeviation(data):F2}");
        Console.WriteLine($"Skewness: {Statistics.Skewness(data):F2}");
        Console.WriteLine($"Kurtosis: {Statistics.Kurtosis(data):F2}");

        var data2 = new double[] { 2, 3, 5, 7, 11, 13 };
        Console.WriteLine($"\nCorrelation between two datasets:");
        Console.WriteLine($"Correlation: {Statistics.Correlation(data.Take(6), data2):F2}");
        Console.WriteLine();
    }
    
    static void PrintMatrix(Matrix m)
    {
        for (int i = 0; i < m.Rows; i++)
        {
            for (int j = 0; j < m.Columns; j++)
            {
                Console.Write($"{m[i, j],8:F3} ");
            }
            Console.WriteLine();
        }
    }
}