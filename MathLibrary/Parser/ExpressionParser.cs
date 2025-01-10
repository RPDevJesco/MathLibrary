namespace MathLibrary
{
    public class ExpressionParser
    {
        private readonly Dictionary<string, IMathOperation> _operations;
        private readonly Dictionary<string, double> _constants;
        private Dictionary<string, double> _variables;
        private readonly Dictionary<string, (Func<List<double>, double>, int, int)> _functions;

        public enum TokenType
        {
            Number,
            Operator,
            UnaryOperator,  // New type for unary operators
            LeftParenthesis,
            RightParenthesis,
            Function,
            Variable,
            Constant,
            Comma
        }

        public class Token
        {
            public TokenType Type { get; }
            public string Value { get; }

            public Token(TokenType type, string value)
            {
                Type = type;
                Value = value;
            }

            public override string ToString() => $"{Type,-15} : {Value,-10}";
        }

        public ExpressionParser()
        {
            _operations = InitializeOperations();
            _constants = InitializeConstants();
            _variables = new Dictionary<string, double>();
            _functions = InitializeFunctions();
        }

        private Dictionary<string, IMathOperation> InitializeOperations()
        {
            return new Dictionary<string, IMathOperation>
            {
                ["+"] = new Addition(),
                ["-"] = new Subtraction(),
                ["*"] = new Multiplication(),
                ["/"] = new Division(),
                ["^"] = new Exponentiation(),
                ["%"] = new Modulo(),
                ["root"] = new NthRoot()
            };
        }

        private Dictionary<string, double> InitializeConstants()
        {
            return new Dictionary<string, double>
            {
                ["pi"] = Math.PI,
                ["π"] = Math.PI,
                ["e"] = Math.E,
                ["phi"] = (1 + Math.Sqrt(5)) / 2,
                ["φ"] = (1 + Math.Sqrt(5)) / 2,
                ["inf"] = double.PositiveInfinity,
                ["∞"] = double.PositiveInfinity
            };
        }

        private Dictionary<string, (Func<List<double>, double>, int, int)> InitializeFunctions()
        {
            return new Dictionary<string, (Func<List<double>, double>, int, int)>
            {
                // Basic mathematical functions
                ["abs"] = (args => Math.Abs(args[0]), 1, 1),
                ["sign"] = (args => Math.Sign(args[0]), 1, 1),
                ["ceil"] = (args => Math.Ceiling(args[0]), 1, 1),
                ["floor"] = (args => Math.Floor(args[0]), 1, 1),
                ["round"] = (args => args.Count > 1 ? Math.Round(args[0], (int)args[1]) : Math.Round(args[0]), 1, 2),

                // Trigonometric functions
                ["sin"] = (args => Math.Sin(args[0]), 1, 1),
                ["cos"] = (args => Math.Cos(args[0]), 1, 1),
                ["tan"] = (args => Math.Tan(args[0]), 1, 1),
                ["asin"] = (args => Math.Asin(args[0]), 1, 1),
                ["acos"] = (args => Math.Acos(args[0]), 1, 1),
                ["atan"] = (args => Math.Atan(args[0]), 1, 1),
                ["atan2"] = (args => Math.Atan2(args[0], args[1]), 2, 2),

                // Hyperbolic functions
                ["sinh"] = (args => Math.Sinh(args[0]), 1, 1),
                ["cosh"] = (args => Math.Cosh(args[0]), 1, 1),
                ["tanh"] = (args => Math.Tanh(args[0]), 1, 1),

                // Exponential and logarithmic functions
                ["exp"] = (args => Math.Exp(args[0]), 1, 1),
                ["log"] = (args => args.Count > 1 ? Math.Log(args[0], args[1]) : Math.Log10(args[0]), 1, 2),
                ["ln"] = (args => Math.Log(args[0]), 1, 1),

                // Power and root functions
                ["sqrt"] = (args => Math.Sqrt(args[0]), 1, 1),
                ["cbrt"] = (args => Math.Pow(args[0], 1.0/3.0), 1, 1),
                ["pow"] = (args => Math.Pow(args[0], args[1]), 2, 2),

                // Special functions
                ["besselj0"] = (args => SpecialFunctions.BesselJ0(args[0]), 1, 1),
                ["besselj1"] = (args => SpecialFunctions.BesselJ1(args[0]), 1, 1),
                ["legendre"] = (args => SpecialFunctions.LegendreP((int)args[0], args[1]), 2, 2),
                ["gamma"] = (args => SpecialFunctions.Gamma(args[0]), 1, 1),
                ["erf"] = (args => SpecialFunctions.Erf(args[0]), 1, 1),
                ["erfc"] = (args => SpecialFunctions.Erfc(args[0]), 1, 1),
                ["hypergeometric2f1"] = (args => SpecialFunctions.Hypergeometric2F1(args[0], args[1], args[2], args[3]), 4, 4),

                // Vector operations
                ["magnitude"] = (args => new Vector(args.ToArray()).Magnitude, 1, int.MaxValue),
                ["dot"] = (args => {
                    if (args.Count % 2 != 0)
                        throw new ArgumentException("Dot product requires an even number of arguments");
                    int mid = args.Count / 2;
                    var v1 = new Vector(args.Take(mid).ToArray());
                    var v2 = new Vector(args.Skip(mid).Take(mid).ToArray());
                    return Vector.DotProduct(v1, v2);
                }, 2, int.MaxValue),
                ["cross"] = (args => {
                    if (args.Count != 6)
                        throw new ArgumentException("Cross product requires exactly 6 arguments (two 3D vectors)");
                    var v1 = new Vector(args.Take(3).ToArray());
                    var v2 = new Vector(args.Skip(3).Take(3).ToArray());
                    var result = Vector.CrossProduct(v1, v2);
                    return result.Magnitude;
                }, 6, 6)
            };
        }

        public void SetVariable(string name, double value)
        {
            _variables[name.ToLower()] = value;
        }

        public void ClearVariables()
        {
            _variables.Clear();
        }
        
        private ExpressionNode CreateFunctionNode(string functionName, List<ExpressionNode> arguments)
        {
            var (func, minArgs, maxArgs) = _functions[functionName];
    
            if (arguments.Count < minArgs || arguments.Count > maxArgs)
            {
                if (minArgs == maxArgs)
                    throw new ArgumentException($"Function {functionName} expects {minArgs} argument(s), but got {arguments.Count}");
                throw new ArgumentException($"Function {functionName} expects between {minArgs} and {maxArgs} arguments, but got {arguments.Count}");
            }

            return new ExpressionNode(new FunctionOperation(functionName, args => func(args), arguments));
        }

        public ExpressionNode Parse(string expression)
        {
            if (string.IsNullOrWhiteSpace(expression))
                throw new ArgumentException("Expression cannot be empty");

            try
            {
                var tokens = TokenizeExpression(expression);
                ValidateTokenSequence(tokens);
                return BuildExpressionTree(tokens);
            }
            catch (Exception ex) when (ex is not ArgumentException)
            {
                throw new ArgumentException($"Invalid expression: {expression}", ex);
            }
        }

        private List<Token> TokenizeExpression(string expression)
        {
            var tokens = new List<Token>();
            var currentToken = new System.Text.StringBuilder();
            bool expectingValue = true;  // True if we're expecting a value or unary operator

            for (int i = 0; i < expression.Length; i++)
            {
                char c = expression[i];

                // Skip whitespace
                if (char.IsWhiteSpace(c))
                    continue;

                // Handle consecutive unary operators
                if ((c == '+' || c == '-') && expectingValue)
                {
                    int negCount = 0;
                    while (i < expression.Length && (expression[i] == '+' || expression[i] == '-'))
                    {
                        if (expression[i] == '-')
                            negCount++;
                        i++;
                    }
                    i--;  // Back up one character since we'll increment in the loop
                    
                    // If odd number of minus signs, add a unary minus
                    if (negCount % 2 == 1)
                        tokens.Add(new Token(TokenType.UnaryOperator, "neg"));
                        
                    continue;
                }

                // Handle numbers
                if (char.IsDigit(c) || c == '.')
                {
                    currentToken.Clear();
                    currentToken.Append(c);
                    i++;

                    // Read the integer/decimal part
                    while (i < expression.Length && 
                           (char.IsDigit(expression[i]) || expression[i] == '.'))
                    {
                        currentToken.Append(expression[i]);
                        i++;
                    }

                    // Handle scientific notation
                    if (i < expression.Length && 
                        (expression[i] == 'e' || expression[i] == 'E'))
                    {
                        currentToken.Append(expression[i]);
                        i++;

                        // Handle optional sign in exponent
                        if (i < expression.Length && 
                            (expression[i] == '+' || expression[i] == '-'))
                        {
                            currentToken.Append(expression[i]);
                            i++;
                        }

                        // Read exponent
                        while (i < expression.Length && char.IsDigit(expression[i]))
                        {
                            currentToken.Append(expression[i]);
                            i++;
                        }
                    }
                    i--;

                    // Try to parse the number
                    if (double.TryParse(currentToken.ToString(), 
                        System.Globalization.NumberStyles.Float, 
                        System.Globalization.CultureInfo.InvariantCulture, 
                        out double value))
                    {
                        tokens.Add(new Token(TokenType.Number, currentToken.ToString()));
                        expectingValue = false;
                    }
                    else
                    {
                        throw new ArgumentException($"Invalid number format: {currentToken}");
                    }
                    continue;
                }

                // Handle operators
                if (_operations.ContainsKey(c.ToString()))
                {
                    if (!expectingValue)  // Binary operator
                    {
                        tokens.Add(new Token(TokenType.Operator, c.ToString()));
                        expectingValue = true;
                    }
                    continue;
                }

                // Handle parentheses
                if (c == '(')
                {
                    if (!expectingValue && i > 0 && 
                        (tokens.Last().Type == TokenType.Number || 
                         tokens.Last().Type == TokenType.Constant ||
                         tokens.Last().Type == TokenType.RightParenthesis))
                    {
                        // Implicit multiplication
                        tokens.Add(new Token(TokenType.Operator, "*"));
                    }
                    tokens.Add(new Token(TokenType.LeftParenthesis, "("));
                    expectingValue = true;
                    continue;
                }
                if (c == ')')
                {
                    tokens.Add(new Token(TokenType.RightParenthesis, ")"));
                    expectingValue = false;
                    continue;
                }

                // Handle comma
                if (c == ',')
                {
                    tokens.Add(new Token(TokenType.Comma, ","));
                    expectingValue = true;
                    continue;
                }

                // Handle functions and constants
                if (char.IsLetter(c) || c == '_')
                {
                    if (!expectingValue)
                    {
                        // Implicit multiplication before identifier
                        tokens.Add(new Token(TokenType.Operator, "*"));
                    }

                    currentToken.Clear();
                    currentToken.Append(c);
                    i++;

                    while (i < expression.Length && 
                           (char.IsLetterOrDigit(expression[i]) || expression[i] == '_'))
                    {
                        currentToken.Append(expression[i]);
                        i++;
                    }
                    i--;

                    string identifier = currentToken.ToString().ToLower();
                    if (_functions.ContainsKey(identifier))
                    {
                        tokens.Add(new Token(TokenType.Function, identifier));
                        expectingValue = false;
                    }
                    else if (_constants.ContainsKey(identifier))
                    {
                        tokens.Add(new Token(TokenType.Constant, identifier));
                        expectingValue = false;
                    }
                    else
                    {
                        throw new ArgumentException($"Unexpected token: {identifier}");
                    }
                    continue;
                }

                throw new ArgumentException($"Invalid character in expression: {c}");
            }

            return tokens;
        }

        private void ValidateTokenSequence(List<Token> tokens)
        {
            if (tokens.Count == 0)
                throw new ArgumentException("Empty expression");

            // Check for invalid token sequences
            for (int i = 0; i < tokens.Count - 1; i++)
            {
                var current = tokens[i];
                var next = tokens[i + 1];

                // Two operators in a row (unless second is unary minus)
                if (current.Type == TokenType.Operator && next.Type == TokenType.Operator)
                {
                    if (!(next.Value == "-" && 
                          (i == 0 || tokens[i - 1].Type == TokenType.LeftParenthesis)))
                        throw new ArgumentException($"Invalid operator sequence: {current.Value}{next.Value}");
                }

                // Two numbers in a row
                if (current.Type == TokenType.Number && next.Type == TokenType.Number)
                    throw new ArgumentException("Invalid number sequence");

                // Missing argument between parentheses
                if (current.Type == TokenType.LeftParenthesis && 
                    next.Type == TokenType.RightParenthesis)
                    throw new ArgumentException("Empty parentheses");
            }

            // Check parentheses matching
            int parenthesesCount = 0;
            foreach (var token in tokens)
            {
                if (token.Type == TokenType.LeftParenthesis)
                    parenthesesCount++;
                else if (token.Type == TokenType.RightParenthesis)
                {
                    parenthesesCount--;
                    if (parenthesesCount < 0)
                        throw new ArgumentException("Mismatched parentheses");
                }
            }

            if (parenthesesCount != 0)
                throw new ArgumentException("Mismatched parentheses");
        }

        private ExpressionNode BuildExpressionTree(List<Token> tokens)
        {
            if (tokens.Count == 0)
                throw new ArgumentException("Empty expression");

            return BuildExpressionTreeRecursive(tokens, 0, tokens.Count - 1);
        }

        private void DebugPrintTokens(List<Token> tokens, int start, int end)
        {
            //Console.WriteLine($"Tokens from {start} to {end}:");
            for (int i = start; i <= end; i++)
            {
                //Console.WriteLine($"Token {i}: Type={tokens[i].Type}, Value={tokens[i].Value}");
            }
        }

        private ExpressionNode BuildExpressionTreeRecursive(List<Token> tokens, int start, int end)
        {
            if (start > end)
                throw new ArgumentException("Invalid expression range");

            // Handle single token
            if (start == end)
            {
                var token = tokens[start];
                switch (token.Type)
                {
                    case TokenType.Number:
                        return new ExpressionNode(double.Parse(token.Value, System.Globalization.CultureInfo.InvariantCulture));
                    case TokenType.Variable:
                        if (!_variables.TryGetValue(token.Value, out double varValue))
                            throw new ArgumentException($"Undefined variable: {token.Value}");
                        return new ExpressionNode(varValue);
                    case TokenType.Constant:
                        return new ExpressionNode(_constants[token.Value]);
                    default:
                        throw new ArgumentException($"Unexpected token: {token.Value}");
                }
            }

            // Find the highest-level operator (lowest precedence) outside any parentheses
            int parenthesesCount = 0;
            int operationIndex = -1;
            int lowestPrecedence = int.MaxValue;

            // Scan from right to left to handle precedence correctly
            for (int i = end; i >= start; i--)
            {
                var token = tokens[i];

                if (token.Type == TokenType.RightParenthesis)
                {
                    parenthesesCount++;
                }
                else if (token.Type == TokenType.LeftParenthesis)
                {
                    parenthesesCount--;
                }
                else if (parenthesesCount == 0 && token.Type == TokenType.Operator)
                {
                    var operation = _operations[token.Value];
                    if (operation.Precedence <= lowestPrecedence)
                    {
                        lowestPrecedence = operation.Precedence;
                        operationIndex = i;
                    }
                }
            }
            
            // Handle unary operator
            if (tokens[start].Type == TokenType.UnaryOperator)
            {
                if (tokens[start].Value == "neg")
                {
                    return new ExpressionNode(
                        _operations["*"],
                        new ExpressionNode(-1),
                        BuildExpressionTreeRecursive(tokens, start + 1, end)
                    );
                }
            }

            // Handle binary operations
            if (operationIndex != -1)
            {
                return new ExpressionNode(
                    _operations[tokens[operationIndex].Value],
                    BuildExpressionTreeRecursive(tokens, start, operationIndex - 1),
                    BuildExpressionTreeRecursive(tokens, operationIndex + 1, end)
                );
            }

            // Handle parenthesized expression
            if (tokens[start].Type == TokenType.LeftParenthesis && tokens[end].Type == TokenType.RightParenthesis)
            {
                return BuildExpressionTreeRecursive(tokens, start + 1, end - 1);
            }

            // Handle function calls
            if (tokens[start].Type == TokenType.Function)
            {
                string functionName = tokens[start].Value;

                // Extract arguments
                var arguments = new List<ExpressionNode>();
                int currentArgStart = start + 2; // Skip function name and opening parenthesis
                parenthesesCount = 1;

                for (int i = currentArgStart; i < end; i++)
                {
                    if (tokens[i].Type == TokenType.LeftParenthesis)
                        parenthesesCount++;
                    else if (tokens[i].Type == TokenType.RightParenthesis)
                        parenthesesCount--;
                    else if (parenthesesCount == 1 && tokens[i].Type == TokenType.Comma)
                    {
                        arguments.Add(BuildExpressionTreeRecursive(tokens, currentArgStart, i - 1));
                        currentArgStart = i + 1;
                    }
                }

                if (currentArgStart < end)
                    arguments.Add(BuildExpressionTreeRecursive(tokens, currentArgStart, end - 1));

                if (!_functions.ContainsKey(functionName))
                    throw new ArgumentException($"Unknown function: {functionName}");

                var (func, minArgs, maxArgs) = _functions[functionName];

                // Validate argument count
                if (arguments.Count < minArgs || arguments.Count > maxArgs)
                {
                    if (minArgs == maxArgs)
                        throw new ArgumentException($"Function {functionName} expects {minArgs} argument(s), but got {arguments.Count}");
                    throw new ArgumentException($"Function {functionName} expects between {minArgs} and {maxArgs} arguments, but got {arguments.Count}");
                }

                return new ExpressionNode(new FunctionOperation(functionName, func, arguments));
            }

            throw new ArgumentException($"Invalid expression");
        }
        
        public List<Token> ParseAndPrintTokens(string expression)
        {
            var tokens = TokenizeExpression(expression);
            Console.WriteLine("Tokens:");
            foreach (var token in tokens)
            {
                Console.WriteLine($"  {token}");
            }
            return tokens;
        }
    }
}