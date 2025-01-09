# MathLibrary

## Overview
MathLibrary is a comprehensive numerical and symbolic computation library implemented in C#. It is designed to cater to a wide range of mathematical operations, including basic arithmetic, advanced calculus, linear algebra, numerical methods, and more. The library provides functionality for developers, engineers, data scientists, and educators to seamlessly integrate complex mathematical computations into their projects.

### Key Features
- **Expression Parsing and Evaluation**: Parse and evaluate mathematical expressions in string format.
- **Vector and Matrix Operations**: Includes linear algebra utilities like dot product, cross product, and matrix decompositions.
- **Special Functions**: Gamma functions, Bessel functions, Legendre polynomials, and hypergeometric functions.
- **Numerical Methods**: Supports numerical integration and solutions to ordinary differential equations (ODEs).
- **Polynomial Manipulation**: Evaluate, derive, and find roots of polynomials.
- **Statistical Analysis**: Compute statistical measures like mean, variance, skewness, kurtosis, and correlation.
- **Complex Numbers**: Operations and transformations for complex numbers.
- **Quaternions**: Enables 3D rotations and quaternion arithmetic.
- **Fast Fourier Transform (FFT)**: Perform frequency domain analysis.

### Use Cases
- Scientific and Engineering Applications
- Data Analysis and Signal Processing
- Educational Tools for Numerical Methods
- Game Development (Quaternions, Vectors, Matrix operations)

---

## API Documentation and Class Design

### Core Components

#### `Calculator`
- **Description**: Evaluates mathematical expressions.
- **Methods**:
    - `double Evaluate(string expression)`: Parses and computes the result of the input expression.

#### `ExpressionParser`
- **Description**: Parses mathematical expressions into an abstract syntax tree (AST).
- **Methods**:
    - `ExpressionNode Parse(string expression)`: Converts an expression string into an AST.
    - `void SetVariable(string name, double value)`: Assigns values to variables in expressions.

#### `ExpressionNode`
- **Description**: Represents nodes in the AST for expression evaluation.
- **Properties**:
    - `double? Value`: Represents constant values.
    - `Complex? ComplexValue`: Represents complex values.
    - `IMathOperation Operation`: Represents the operation (e.g., addition).
- **Methods**:
    - `double Evaluate()`: Evaluates the node for real numbers.
    - `Complex EvaluateComplex()`: Evaluates the node for complex numbers.

#### `IMathOperation`
- **Description**: Interface for mathematical operations.
- **Methods**:
    - `double Execute(double left, double right)`: Executes operation for real numbers.
    - `Complex ExecuteComplex(Complex left, Complex right)`: Executes operation for complex numbers.

#### `MathOperations`
- **Description**: Implements operations like addition, subtraction, multiplication, and more.
- **Classes**:
    - `Addition`, `Subtraction`, `Multiplication`, `Division`, `Exponentiation`, `Modulo`, `NthRoot`.

### Advanced Components

#### `Matrix`
- **Description**: Handles matrix operations.
- **Methods**:
    - `Matrix Transpose()`: Computes the transpose of the matrix.
    - `Matrix operator+(Matrix other)`: Adds two matrices.
    - `Matrix operator*(Matrix other)`: Multiplies two matrices.

#### `Vector`
- **Description**: Implements vector operations.
- **Methods**:
    - `double Magnitude`: Computes vector magnitude.
    - `Vector Normalize()`: Normalizes the vector.
    - `static double DotProduct(Vector a, Vector b)`: Computes the dot product.
    - `static Vector CrossProduct(Vector a, Vector b)`: Computes the cross product.

#### `MatrixDecomposition`
- **Description**: Implements LU, QR, and SVD decompositions.
- **Methods**:
    - `LUDecomposition(Matrix A)`: Returns L, U, and P matrices for LU decomposition.
    - `QRDecomposition(Matrix A)`: Computes Q and R matrices for QR decomposition.
    - `SVDDecomposition(Matrix A)`: Computes singular value decomposition.

#### `Statistics`
- **Description**: Provides statistical functions.
- **Methods**:
    - `double Mean(IEnumerable<double> values)`: Calculates the mean.
    - `double Variance(IEnumerable<double> values)`: Calculates the variance.
    - `double Skewness(IEnumerable<double> values)`: Computes skewness.
    - `double Kurtosis(IEnumerable<double> values)`: Computes kurtosis.
    - `double Correlation(IEnumerable<double> x, IEnumerable<double> y)`: Computes the correlation.

#### `SpecialFunctions`
- **Description**: Implements mathematical special functions.
- **Methods**:
    - `double Gamma(double x)`: Computes the Gamma function.
    - `double BesselJ0(double x)`: Computes the Bessel function J0.
    - `double LegendreP(int n, double x)`: Computes Legendre polynomials.

#### `DifferentialEquations`
- **Description**: Solves ordinary differential equations.
- **Methods**:
    - `RungeKutta4(Func<double, Vector, Vector> func, Vector y0, double t0, double tf, int steps)`: Implements the Runge-Kutta method.

#### `NumericalIntegration`
- **Description**: Provides numerical integration methods.
- **Methods**:
    - `double SimpsonsRule(Func<double, double> f, double a, double b, int n)`: Implements Simpson's rule.
    - `double GaussianQuadrature(Func<double, double> f, double a, double b, int n)`: Implements Gaussian quadrature.

#### `Quaternion`
- **Description**: Handles quaternion operations.
- **Methods**:
    - `Quaternion FromAxisAngle(Vector axis, double angle)`: Creates a quaternion from axis-angle representation.
    - `Quaternion Rotate(Vector vector)`: Rotates a vector using the quaternion.

---

## How to Use
1. **Basic Example**:
```csharp
using MathLibrary;

class Program
{
    static void Main()
    {
        Calculator calc = new Calculator();
        Console.WriteLine(calc.Evaluate("3 + 5 * 2")); // Output: 13
    }
}
```

3. **Advanced Use**: Explore the demos in the `Program.cs` file for detailed use cases.

---
## License
This project is licensed under the MIT License.
