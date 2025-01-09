namespace MathLibrary
{
    public class Vector
    {
        private readonly double[] _components;

        public int Dimension => _components.Length;

        public double this[int index]
        {
            get
            {
                if (index < 0 || index >= Dimension)
                    throw new ArgumentOutOfRangeException(nameof(index));
                return _components[index];
            }
            set
            {
                if (index < 0 || index >= Dimension)
                    throw new ArgumentOutOfRangeException(nameof(index));
                _components[index] = value;
            }
        }

        public Vector(params double[] components)
        {
            if (components == null || components.Length == 0)
                throw new ArgumentException("Vector must have at least one component");
            _components = (double[])components.Clone();
        }

        public Vector(int dimension)
        {
            if (dimension <= 0)
                throw new ArgumentException("Vector dimension must be positive");
            _components = new double[dimension];
        }

        public static Vector operator +(Vector a, Vector b)
        {
            ThrowIfDimensionMismatch(a, b);
            return new Vector(
                Enumerable.Range(0, a.Dimension)
                         .Select(i => a[i] + b[i])
                         .ToArray()
            );
        }

        public static Vector operator -(Vector a, Vector b)
        {
            ThrowIfDimensionMismatch(a, b);
            return new Vector(
                Enumerable.Range(0, a.Dimension)
                         .Select(i => a[i] - b[i])
                         .ToArray()
            );
        }

        public static Vector operator *(Vector v, double scalar)
        {
            return new Vector(v._components.Select(x => x * scalar).ToArray());
        }

        public static Vector operator /(Vector v, double scalar)
        {
            if (scalar == 0)
                throw new DivideByZeroException("Cannot divide vector by zero");
            return new Vector(v._components.Select(x => x / scalar).ToArray());
        }

        public static double DotProduct(Vector a, Vector b)
        {
            ThrowIfDimensionMismatch(a, b);
            return Enumerable.Range(0, a.Dimension)
                            .Sum(i => a[i] * b[i]);
        }

        public static Vector CrossProduct(Vector a, Vector b)
        {
            if (a.Dimension != 3 || b.Dimension != 3)
                throw new ArgumentException("Cross product is only defined for 3D vectors");

            return new Vector(
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
            );
        }

        public double Magnitude => Math.Sqrt(_components.Sum(x => x * x));

        public Vector Normalize()
        {
            double magnitude = Magnitude;
            if (magnitude == 0)
                throw new InvalidOperationException("Cannot normalize zero vector");
            return this / magnitude;
        }

        public static Vector Project(Vector a, Vector b)
        {
            ThrowIfDimensionMismatch(a, b);
            double scalar = DotProduct(a, b) / DotProduct(b, b);
            return b * scalar;
        }

        public static double Angle(Vector a, Vector b)
        {
            ThrowIfDimensionMismatch(a, b);
            double dotProduct = DotProduct(a, b);
            double magnitudeProduct = a.Magnitude * b.Magnitude;
            
            if (magnitudeProduct == 0)
                throw new ArgumentException("Cannot compute angle with zero vector");

            // Handle numerical precision issues
            double cosTheta = dotProduct / magnitudeProduct;
            cosTheta = Math.Max(-1.0, Math.Min(1.0, cosTheta));
            return Math.Acos(cosTheta);
        }

        public static Vector[] GramSchmidt(Vector[] vectors)
        {
            if (vectors == null || vectors.Length == 0)
                throw new ArgumentException("Empty vector set");

            int dimension = vectors[0].Dimension;
            foreach (var v in vectors)
                if (v.Dimension != dimension)
                    throw new ArgumentException("All vectors must have same dimension");

            var result = new Vector[vectors.Length];
            result[0] = vectors[0].Normalize();

            for (int i = 1; i < vectors.Length; i++)
            {
                Vector v = vectors[i];
                for (int j = 0; j < i; j++)
                    v -= Project(vectors[i], result[j]);
                
                if (v.Magnitude < 1e-10)
                    throw new ArgumentException("Vectors are linearly dependent");
                    
                result[i] = v.Normalize();
            }

            return result;
        }

        private static void ThrowIfDimensionMismatch(Vector a, Vector b)
        {
            if (a.Dimension != b.Dimension)
                throw new ArgumentException($"Vector dimensions must match: {a.Dimension} != {b.Dimension}");
        }

        public override string ToString()
        {
            return $"[{string.Join(", ", _components)}]";
        }
    }
}