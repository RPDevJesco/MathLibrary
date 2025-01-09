namespace MathLibrary
{
    // Base class for mathematical operations
    public abstract class MathOperationBase : IMathOperation
    {
        public abstract double Execute(double left, double right);
        public virtual Complex ExecuteComplex(Complex left, Complex right)
        {
            // Default implementation converts the result of real operation to complex
            // Operations that need different complex behavior should override this
            return new Complex(Execute(left.Real, right.Real), 0);
        }

        public abstract int Precedence { get; }
    }
}