namespace MathLibrary
{
    public class Matrix
    {
        private readonly double[,] _elements;
        public int Rows { get; }
        public int Columns { get; }

        public double this[int row, int col]
        {
            get
            {
                ValidateIndices(row, col);
                return _elements[row, col];
            }
            set
            {
                ValidateIndices(row, col);
                _elements[row, col] = value;
            }
        }

        // Constructor for creating a zero matrix of specified size
        public Matrix(int rows, int columns)
        {
            if (rows <= 0 || columns <= 0)
                throw new ArgumentException("Matrix dimensions must be positive");

            Rows = rows;
            Columns = columns;
            _elements = new double[rows, columns];
        }

        // Constructor from 2D array
        public Matrix(double[,] elements)
        {
            if (elements == null)
                throw new ArgumentNullException(nameof(elements));

            Rows = elements.GetLength(0);
            Columns = elements.GetLength(1);
            _elements = (double[,])elements.Clone();
        }

        // Copy constructor
        public Matrix(Matrix other)
        {
            if (other == null)
                throw new ArgumentNullException(nameof(other));

            Rows = other.Rows;
            Columns = other.Columns;
            _elements = (double[,])other._elements.Clone();
        }

        // Static factory method for identity matrix
        public static Matrix Identity(int size)
        {
            var matrix = new Matrix(size, size);
            for (int i = 0; i < size; i++)
                matrix[i, i] = 1.0;
            return matrix;
        }

        // Helper method for index validation
        private void ValidateIndices(int row, int col)
        {
            if (row < 0 || row >= Rows)
                throw new ArgumentOutOfRangeException(nameof(row));
            if (col < 0 || col >= Columns)
                throw new ArgumentOutOfRangeException(nameof(col));
        }

        public Matrix Transpose()
        {
            var result = new Matrix(Columns, Rows);
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    result[j, i] = _elements[i, j];
            return result;
        }

        // Basic matrix operations
        public static Matrix operator +(Matrix a, Matrix b)
        {
            if (a.Rows != b.Rows || a.Columns != b.Columns)
                throw new ArgumentException("Matrix dimensions must match for addition");

            var result = new Matrix(a.Rows, a.Columns);
            for (int i = 0; i < a.Rows; i++)
                for (int j = 0; j < a.Columns; j++)
                    result[i, j] = a[i, j] + b[i, j];
            return result;
        }

        public static Matrix operator -(Matrix a, Matrix b)
        {
            if (a.Rows != b.Rows || a.Columns != b.Columns)
                throw new ArgumentException("Matrix dimensions must match for subtraction");

            var result = new Matrix(a.Rows, a.Columns);
            for (int i = 0; i < a.Rows; i++)
                for (int j = 0; j < a.Columns; j++)
                    result[i, j] = a[i, j] - b[i, j];
            return result;
        }

        public static Matrix operator *(Matrix a, Matrix b)
        {
            if (a.Columns != b.Rows)
                throw new ArgumentException("Matrix dimensions must match for multiplication");

            var result = new Matrix(a.Rows, b.Columns);
            for (int i = 0; i < a.Rows; i++)
                for (int j = 0; j < b.Columns; j++)
                    for (int k = 0; k < a.Columns; k++)
                        result[i, j] += a[i, k] * b[k, j];
            return result;
        }

        public static Matrix operator *(Matrix a, double scalar)
        {
            var result = new Matrix(a.Rows, a.Columns);
            for (int i = 0; i < a.Rows; i++)
                for (int j = 0; j < a.Columns; j++)
                    result[i, j] = a[i, j] * scalar;
            return result;
        }

        public static Matrix operator /(Matrix a, double scalar)
        {
            if (scalar == 0)
                throw new DivideByZeroException();

            var result = new Matrix(a.Rows, a.Columns);
            for (int i = 0; i < a.Rows; i++)
                for (int j = 0; j < a.Columns; j++)
                    result[i, j] = a[i, j] / scalar;
            return result;
        }
    }
}