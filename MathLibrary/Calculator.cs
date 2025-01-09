namespace MathLibrary
{
    public class Calculator
    {
        private readonly ExpressionParser _parser;

        public Calculator()
        {
            _parser = new ExpressionParser();
        }

        public double Evaluate(string expression)
        {
            var expressionTree = _parser.Parse(expression);
            return expressionTree.Evaluate();
        }
    }
}