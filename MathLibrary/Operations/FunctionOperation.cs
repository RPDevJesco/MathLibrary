namespace MathLibrary
{
    public class FunctionOperation : IMathOperation
    {
        private readonly string _name;
        private readonly Func<List<double>, double> _function;
        private readonly Func<List<Complex>, Complex>? _complexFunction;
        public List<ExpressionNode> Arguments { get; }

        public FunctionOperation(string name, Func<List<double>, double> function, List<ExpressionNode> arguments)
        {
            _name = name;
            _function = function;
            _complexFunction = GetComplexImplementation(name);
            Arguments = arguments;
            ValidateArgumentCount(arguments);
        }

        private Func<List<Complex>, Complex>? GetComplexImplementation(string name)
        {
            switch (name.ToLower())
            {
                case "sin": return args => ComplexMath.Sin(args[0]);
                case "cos": return args => ComplexMath.Cos(args[0]);
                case "tan": return args => ComplexMath.Tan(args[0]);
                case "exp": return args => ComplexMath.Exp(args[0]);
                case "ln":  return args => ComplexMath.Log(args[0]);
                case "sqrt": return args => ComplexMath.Sqrt(args[0]);
                case "abs": return args => new Complex(Complex.Abs(args[0]), 0);
                case "pow": return args => ComplexMath.Pow(args[0], args[1]);
                default: return null;
            }
        }

        private void ValidateArgumentCount(List<ExpressionNode> arguments)
        {
            var (minArgs, maxArgs) = GetArgumentRange(_name);
            if (arguments.Count < minArgs || arguments.Count > maxArgs)
            {
                if (minArgs == maxArgs)
                    throw new ArgumentException($"Function {_name} expects {minArgs} argument(s), but got {arguments.Count}");
                throw new ArgumentException($"Function {_name} expects between {minArgs} and {maxArgs} arguments, but got {arguments.Count}");
            }
        }

        private (int min, int max) GetArgumentRange(string name)
        {
            return name.ToLower() switch
            {
                // Single argument functions
                "sin" or "cos" or "tan" or "asin" or "acos" or "atan" or
                    "sinh" or "cosh" or "tanh" or "sqrt" or "cbrt" or "ln" or
                    "abs" or "sign" or "ceil" or "floor" or "exp" or "gamma" or
                    "besselj0" or "besselj1" or "erf" or "erfc" => (1, 1),

                // Two argument functions
                "pow" or "atan2" or "legendre" => (2, 2),

                // Variable argument functions
                "log" or "round" => (1, 2),
                "magnitude" => (1, int.MaxValue),
                "dot" => (2, int.MaxValue),
                "cross" => (6, 6),
                "hypergeometric2f1" => (4, 4),

                // Default case
                _ => throw new ArgumentException($"Unknown function: {name}")
            };
        }

        public double Execute(double left, double right)
        {
            var args = Arguments.Select(arg => arg.Evaluate()).ToList();
            return _function(args);
        }

        public Complex ExecuteComplex(Complex left, Complex right)
        {
            if (_complexFunction == null)
                throw new NotSupportedException($"Function {_name} does not support complex numbers");

            var args = Arguments.Select(arg => arg.EvaluateComplex()).ToList();
            return _complexFunction(args);
        }

        public double ExecuteFunction(List<double> args)
        {
            try
            {
                return _function(args);
            }
            catch (Exception ex)
            {
                throw new ArgumentException($"Error evaluating function {_name}: {ex.Message}");
            }
        }

        public Complex ExecuteFunctionComplex(List<Complex> args)
        {
            if (_complexFunction == null)
                throw new NotSupportedException($"Function {_name} does not support complex numbers");

            try
            {
                return _complexFunction(args);
            }
            catch (Exception ex)
            {
                throw new ArgumentException($"Error evaluating complex function {_name}: {ex.Message}");
            }
        }

        public int Precedence => 4;
    }
}