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
            switch (_name.ToLower())
            {
                // Single argument functions
                case "sin":
                case "cos":
                case "tan":
                case "asin":
                case "acos":
                case "atan":
                case "sinh":
                case "cosh":
                case "tanh":
                case "sqrt":
                case "cbrt":
                case "ln":
                case "abs":
                case "sign":
                case "ceil":
                case "floor":
                case "exp":
                case "gamma":
                case "besselj0":
                case "besselj1":
                case "erf":
                case "erfc":
                    if (arguments.Count != 1)
                        throw new ArgumentException($"Function {_name} expects 1 argument, but got {arguments.Count}");
                    break;

                // Two argument functions
                case "pow":
                case "atan2":
                case "legendre":
                case "hypergeometric2f1":
                    if (arguments.Count != 2)
                        throw new ArgumentException($"Function {_name} expects 2 arguments, but got {arguments.Count}");
                    break;

                // Variable argument functions
                case "log":
                case "round":
                    if (arguments.Count < 1 || arguments.Count > 2)
                        throw new ArgumentException($"Function {_name} expects 1 or 2 arguments, but got {arguments.Count}");
                    break;

                case "magnitude":
                    if (arguments.Count < 1)
                        throw new ArgumentException($"Function {_name} expects at least 1 argument, but got {arguments.Count}");
                    break;

                case "dot":
                    if (arguments.Count < 2 || arguments.Count % 2 != 0)
                        throw new ArgumentException($"Function {_name} expects an even number of arguments (at least 2), but got {arguments.Count}");
                    break;

                case "cross":
                    if (arguments.Count != 6)
                        throw new ArgumentException($"Function {_name} expects exactly 6 arguments (two 3D vectors), but got {arguments.Count}");
                    break;
            }
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