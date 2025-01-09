namespace MathLibrary
{
    public class ExpressionNode
    {
        public double? Value { get; set; }
        public Complex? ComplexValue { get; set; }
        public IMathOperation Operation { get; set; }
        public ExpressionNode Left { get; set; }
        public ExpressionNode Right { get; set; }
        public List<ExpressionNode> FunctionArguments { get; set; }
        public bool IsFunction => FunctionArguments != null;
        public bool IsComplex => ComplexValue.HasValue;

        // Constructor for numeric values
        public ExpressionNode(double value)
        {
            Value = value;
            ComplexValue = null;
        }

        // Constructor for complex values
        public ExpressionNode(Complex value)
        {
            ComplexValue = value;
            Value = null;
        }

        // Constructor for binary operations
        public ExpressionNode(IMathOperation operation, ExpressionNode left, ExpressionNode right)
        {
            Operation = operation;
            Left = left;
            Right = right;
            FunctionArguments = null;
        }

        // Constructor for function operations
        public ExpressionNode(FunctionOperation operation)
        {
            Operation = operation;
            FunctionArguments = operation.Arguments;
            Left = null;
            Right = null;
        }

        public Complex EvaluateComplex()
        {
            if (ComplexValue.HasValue)
                return ComplexValue.Value;
            
            if (Value.HasValue)
                return new Complex(Value.Value, 0);

            if (IsFunction)
            {
                var args = FunctionArguments.Select(arg => arg.EvaluateComplex()).ToList();
                return ((FunctionOperation)Operation).ExecuteFunctionComplex(args);
            }

            Complex leftValue = Left.EvaluateComplex();
            Complex rightValue = Right.EvaluateComplex();
            return Operation.ExecuteComplex(leftValue, rightValue);
        }

        public double Evaluate()
        {
            if (Value.HasValue)
                return Value.Value;

            if (ComplexValue.HasValue)
                throw new InvalidOperationException("Cannot evaluate complex number as real");

            if (IsFunction)
            {
                var args = FunctionArguments.Select(arg => arg.Evaluate()).ToList();
                return ((FunctionOperation)Operation).ExecuteFunction(args);
            }

            double leftValue = Left.Evaluate();
            double rightValue = Right.Evaluate();
            return Operation.Execute(leftValue, rightValue);
        }
    }
}