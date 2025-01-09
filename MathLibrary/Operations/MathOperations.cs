namespace MathLibrary
{
    public class Addition : MathOperationBase
    {
        public override double Execute(double left, double right) => left + right;
        public override Complex ExecuteComplex(Complex left, Complex right) => left + right;
        public override int Precedence => 1;
    }

    public class Subtraction : MathOperationBase
    {
        public override double Execute(double left, double right) => left - right;
        public override Complex ExecuteComplex(Complex left, Complex right) => left - right;
        public override int Precedence => 1;
    }

    public class Multiplication : MathOperationBase
    {
        public override double Execute(double left, double right) => left * right;
        public override Complex ExecuteComplex(Complex left, Complex right) => left * right;
        public override int Precedence => 2;
    }

    public class Division : MathOperationBase
    {
        public override double Execute(double left, double right)
        {
            if (right == 0)
                throw new DivideByZeroException("Division by zero is not allowed.");
            return left / right;
        }

        public override Complex ExecuteComplex(Complex left, Complex right)
        {
            if (right.Real == 0 && right.Imaginary == 0)
                throw new DivideByZeroException("Division by zero is not allowed.");
            return left / right;
        }

        public override int Precedence => 2;
    }
    
 public class Exponentiation : MathOperationBase
    {
        public override double Execute(double left, double right) => Math.Pow(left, right);
        
        public override Complex ExecuteComplex(Complex left, Complex right)
        {
            // Special cases for real exponents to avoid unnecessary complex calculations
            if (right.Imaginary == 0)
            {
                // If base is negative and exponent is integer, we can do real calculation
                if (left.Real < 0 && left.Imaginary == 0 && Math.Floor(right.Real) == right.Real)
                {
                    return new Complex(Math.Pow(left.Real, right.Real), 0);
                }
            }
            return ComplexMath.Pow(left, right);
        }

        public override int Precedence => 3;
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

        public override int Precedence => 2;
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
            
            // If dealing with real numbers and odd root, we can use the real implementation
            if (right.Imaginary == 0 && left.Imaginary == 0 && 
                Math.Floor(right.Real) == right.Real && 
                (int)right.Real % 2 == 1)
            {
                return new Complex(Execute(left.Real, right.Real), 0);
            }
            
            // For complex numbers or even roots of negative numbers, use complex power
            return ComplexMath.Pow(left, new Complex(1.0 / right.Real, -right.Imaginary / 
                (right.Real * right.Real + right.Imaginary * right.Imaginary)));
        }

        public override int Precedence => 3;
    }
}