namespace MathLibrary
{
    public interface IMathOperation
    {
        public abstract double Execute(double left, double right);
        public abstract Complex ExecuteComplex(Complex left, Complex right);
        public abstract int Precedence { get; }
    }
}