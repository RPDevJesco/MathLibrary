namespace MathLibrary
{
    public class Addition : MathOperationBase
    {
        public override double Execute(double left, double right) => left + right;
        public override Complex ExecuteComplex(Complex left, Complex right) => left + right;
        public override int Precedence => 1;  // Lowest precedence
    }

    public class Subtraction : MathOperationBase
    {
        public override double Execute(double left, double right) => left - right;
        public override Complex ExecuteComplex(Complex left, Complex right) => left - right;
        public override int Precedence => 1;  // Same as addition
    }

    public class Multiplication : MathOperationBase
    {
        public override double Execute(double left, double right) => left * right;
        public override Complex ExecuteComplex(Complex left, Complex right) => left * right;
        public override int Precedence => 2;  // Higher than addition/subtraction
    }

    public class Division : MathOperationBase
    {
        public override double Execute(double left, double right)
        {
            if (right == 0)
            {
                if (left == 0)
                    return double.NaN;  // 0/0 is undefined
                if (left > 0)
                    return double.PositiveInfinity;  // Positive divided by zero
                return double.NegativeInfinity;  // Negative divided by zero
            }
            return left / right;
        }

        public override Complex ExecuteComplex(Complex left, Complex right)
        {
            if (right.Real == 0 && right.Imaginary == 0)
            {
                if (left.Real == 0 && left.Imaginary == 0)
                    return new Complex(double.NaN, double.NaN);  // 0/0 is undefined
            
                // This should technically return a complex infinity point,
                // but since Complex doesn't have that concept, we'll use the real infinity
                return new Complex(double.PositiveInfinity, double.PositiveInfinity);
            }
            return left / right;
        }

        public override int Precedence => 2;  // Same as multiplication
    }

    public class Exponentiation : MathOperationBase
    {
        public override double Execute(double left, double right) => Math.Pow(left, right);
        public override Complex ExecuteComplex(Complex left, Complex right) => ComplexMath.Pow(left, right);
        public override int Precedence => 3;  // Highest precedence
    }

    public class Modulo : MathOperationBase
    {
        public override double Execute(double left, double right)
        {
            if (right == 0)
                throw new DivideByZeroException("Modulo by zero is not allowed.");
            return left % right;
        }

        public override Complex ExecuteComplex(Complex left, Complex right)
        {
            throw new NotSupportedException("Complex modulo operation is not supported.");
        }

        public override int Precedence => 2;  // Same as multiplication/division
    }

    public class NthRoot : MathOperationBase
    {
        public override double Execute(double left, double right)
        {
            if (right == 0)
                throw new ArgumentException("Root index cannot be zero.");
            if (left < 0 && right % 2 == 0)
                throw new ArgumentException("Even root of negative number is not real.");
            return Math.Sign(left) * Math.Pow(Math.Abs(left), 1.0 / right);
        }

        public override Complex ExecuteComplex(Complex left, Complex right)
        {
            if (right.Real == 0 && right.Imaginary == 0)
                throw new ArgumentException("Root index cannot be zero.");
            return ComplexMath.Pow(left, new Complex(1.0 / right.Real, -right.Imaginary / 
                (right.Real * right.Real + right.Imaginary * right.Imaginary)));
        }

        public override int Precedence => 3;  // Same as exponentiation
    }
}