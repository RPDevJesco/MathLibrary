namespace MathLibrary
{
    public static class MatrixDecomposition
    {
        private const double Epsilon = 1e-10;

        #region LU Decomposition
        public static (Matrix L, Matrix U, Matrix P, int exchanges) LUDecomposition(Matrix A)
        {
            if (A.Rows != A.Columns)
                throw new ArgumentException("Matrix must be square for LU decomposition");

            int n = A.Rows;
            Matrix L = new Matrix(n, n);
            Matrix U = new Matrix(A); // Create a copy to work with
            Matrix P = Matrix.Identity(n);
            int exchanges = 0;

            // Initialize L diagonal to 1
            for (int i = 0; i < n; i++)
                L[i, i] = 1.0;

            // Compute matrix norm for scaling
            double matrixNorm = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    matrixNorm = Math.Max(matrixNorm, Math.Abs(U[i, j]));

            double tolerance = matrixNorm * Epsilon * n;

            // Compute row scaling factors
            var rowScale = new double[n];
            for (int i = 0; i < n; i++)
            {
                double maxRowElement = 0.0;
                for (int j = 0; j < n; j++)
                {
                    double absValue = Math.Abs(U[i, j]);
                    if (absValue > maxRowElement)
                        maxRowElement = absValue;
                }
                if (maxRowElement < tolerance)
                    throw new InvalidOperationException("Matrix is singular or nearly singular");
                rowScale[i] = 1.0 / maxRowElement;
            }

            for (int k = 0; k < n - 1; k++)
            {
                // Find pivot using scaled comparison
                int pivotRow = k;
                double maxScaledPivot = Math.Abs(U[k, k]) * rowScale[k];

                for (int i = k + 1; i < n; i++)
                {
                    double scaledValue = Math.Abs(U[i, k]) * rowScale[i];
                    if (scaledValue > maxScaledPivot)
                    {
                        maxScaledPivot = scaledValue;
                        pivotRow = i;
                    }
                }

                // Check for numerical stability with tolerance
                if (maxScaledPivot < tolerance)
                {
                    // Try to recover by finding a non-zero element in this column
                    bool found = false;
                    for (int i = k + 1; i < n && !found; i++)
                    {
                        if (Math.Abs(U[i, k]) >= tolerance)
                        {
                            pivotRow = i;
                            found = true;
                        }
                    }
                    if (!found)
                        throw new InvalidOperationException("Matrix is singular or nearly singular");
                }

                // Swap rows if necessary
                if (pivotRow != k)
                {
                    // Swap row scaling factors
                    (rowScale[k], rowScale[pivotRow]) = (rowScale[pivotRow], rowScale[k]);

                    for (int j = 0; j < n; j++)
                    {
                        // Swap rows in U
                        (U[k, j], U[pivotRow, j]) = (U[pivotRow, j], U[k, j]);
                        // Swap rows in P
                        (P[k, j], P[pivotRow, j]) = (P[pivotRow, j], P[k, j]);
                        // Swap rows in L up to column k
                        if (j < k)
                            (L[k, j], L[pivotRow, j]) = (L[pivotRow, j], L[k, j]);
                    }
                    exchanges++;
                }

                // Store the pivot element
                double pivot = U[k, k];

                // Perform elimination with scaled pivoting
                for (int i = k + 1; i < n; i++)
                {
                    double factor = U[i, k] / pivot;
                    if (Math.Abs(factor) > 1e15) // Check for potential overflow
                        throw new InvalidOperationException("Matrix is poorly conditioned");
                        
                    L[i, k] = factor;

                    // Update row
                    U[i, k] = 0.0; // Explicitly set to zero to avoid roundoff
                    for (int j = k + 1; j < n; j++)
                    {
                        U[i, j] -= factor * U[k, j];
                        if (Math.Abs(U[i, j]) < tolerance)
                            U[i, j] = 0.0; // Clean up small values
                    }
                }
            }

            return (L, U, P, exchanges);
        }
        #endregion

        #region QR Decomposition
        public static (Matrix Q, Matrix R) QRDecomposition(Matrix A)
        {
            int m = A.Rows;
            int n = A.Columns;
            Matrix Q = new Matrix(m, n);
            Matrix R = new Matrix(n, n);

            // Copy A into Q
            for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                Q[i, j] = A[i, j];

            // Modified Gram-Schmidt process
            for (int k = 0; k < n; k++)
            {
                // Compute the norm of column k
                double norm = 0;
                for (int i = 0; i < m; i++)
                    norm += Q[i, k] * Q[i, k];
                norm = Math.Sqrt(norm);

                if (norm > Epsilon)
                {
                    // Store R[k,k] and normalize column k of Q
                    R[k, k] = norm;
                    for (int i = 0; i < m; i++)
                        Q[i, k] /= norm;

                    // Remove projections of column k from remaining columns
                    for (int j = k + 1; j < n; j++)
                    {
                        // Compute projection coefficient
                        double dot = 0;
                        for (int i = 0; i < m; i++)
                            dot += Q[i, k] * Q[i, j];

                        // Store in R
                        R[k, j] = dot;

                        // Subtract projection
                        for (int i = 0; i < m; i++)
                            Q[i, j] -= dot * Q[i, k];
                    }
                }
            }

            // Clean small values
            for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                if (Math.Abs(Q[i, j]) < Epsilon)
                    Q[i, j] = 0;

            for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (Math.Abs(R[i, j]) < Epsilon || i > j)  // Set below-diagonal elements to 0
                    R[i, j] = 0;

            return (Q, R);
        }
        #endregion

        #region SVD Decomposition
        public static (Matrix U, Vector S, Matrix V) SVDDecomposition(Matrix A, int maxIterations = 100)
        {
            int m = A.Rows;
            int n = A.Columns;

            // Initialize matrices
            Matrix U = new Matrix(A);
            Vector S = new Vector(Math.Min(m, n));
            Matrix V = Matrix.Identity(n);

            // Bidiagonalization (Householder)
            Bidiagonalize(U, S, V);

            // Diagonalization using QR algorithm with implicit shifts
            DiagonalizeQR(U, S, V, maxIterations);

            // Sort singular values in descending order
            SortSVD(U, S, V);

            return (U, S, V);
        }

        private static void Bidiagonalize(Matrix U, Vector S, Matrix V)
        {
            int m = U.Rows;
            int n = U.Columns;
            Vector superdiagonal = new Vector(Math.Min(m - 1, n));

            for (int k = 0; k < Math.Min(m - 1, n); k++)
            {
                // Left Householder transformation
                double leftNorm = 0;
                for (int i = k; i < m; i++)
                    leftNorm += U[i, k] * U[i, k];
                leftNorm = Math.Sqrt(leftNorm);

                if (leftNorm > Epsilon)
                {
                    if (U[k, k] < 0)
                        leftNorm = -leftNorm;

                    for (int i = k; i < m; i++)
                        U[i, k] /= leftNorm;
                    U[k, k] += 1;

                    // Apply transformation to remaining columns
                    for (int j = k + 1; j < n; j++)
                    {
                        double s = 0;
                        for (int i = k; i < m; i++)
                            s += U[i, k] * U[i, j];
                        s /= U[k, k];

                        for (int i = k; i < m; i++)
                            U[i, j] -= s * U[i, k];
                    }
                }
                S[k] = leftNorm;

                // Right Householder transformation (if needed)
                if (k < n - 1)
                {
                    double rightNorm = 0;
                    for (int j = k + 1; j < n; j++)
                        rightNorm += U[k, j] * U[k, j];
                    rightNorm = Math.Sqrt(rightNorm);

                    if (rightNorm > Epsilon)
                    {
                        if (U[k, k + 1] < 0)
                            rightNorm = -rightNorm;

                        for (int j = k + 1; j < n; j++)
                            U[k, j] /= rightNorm;
                        U[k, k + 1] += 1;

                        // Apply transformation to remaining rows
                        for (int i = k + 1; i < m; i++)
                        {
                            double s = 0;
                            for (int j = k + 1; j < n; j++)
                                s += U[k, j] * U[i, j];
                            s /= U[k, k + 1];

                            for (int j = k + 1; j < n; j++)
                                U[i, j] -= s * U[k, j];
                        }
                    }
                    superdiagonal[k] = rightNorm;
                }
            }
        }

        private static void DiagonalizeQR(Matrix U, Vector S, Matrix V, int maxIterations)
        {
            int m = U.Rows;
            int n = U.Columns;
            int p = Math.Min(m, n);
            Vector e = new Vector(p);
            
            // Initialize diagonal and superdiagonal
            for (int i = 0; i < p - 1; i++)
                e[i] = S[i + 1];
            e[p - 1] = 0;

            // Main iteration loop
            for (int k = p - 1; k >= 0; k--)
            {
                int iterations = 0;
                while (iterations < maxIterations)
                {
                    // Test for splitting
                    int l;
                    for (l = k; l >= 0; l--)
                    {
                        if (Math.Abs(e[l]) <= Epsilon)
                            break;
                        if (Math.Abs(S[l - 1]) <= Epsilon)
                        {
                            l--;
                            break;
                        }
                    }

                    if (l == k)
                    {
                        // Convergence achieved
                        break;
                    }
                    else
                    {
                        // QR step
                        double c, s;
                        for (int i = l + 1; i <= k; i++)
                        {
                            GivensRotation(S[i - 1], e[i - 1], out c, out s);
                            ApplyGivens(U, V, i - 1, i, c, s);
                            S[i - 1] = c * S[i - 1] + s * e[i - 1];
                            e[i - 1] = -s * S[i - 1] + c * e[i - 1];
                        }
                    }
                    iterations++;
                }

                if (iterations >= maxIterations)
                    throw new InvalidOperationException("SVD failed to converge");
            }
        }

        private static void GivensRotation(double a, double b, out double c, out double s)
        {
            if (Math.Abs(b) < Epsilon)
            {
                c = 1;
                s = 0;
            }
            else if (Math.Abs(a) < Epsilon)
            {
                c = 0;
                s = Math.Sign(b);
            }
            else
            {
                double r = Math.Sqrt(a * a + b * b);
                c = Math.Abs(a) / r;
                s = Math.Sign(a) * b / r;
            }
        }

        private static void ApplyGivens(Matrix U, Matrix V, int i, int j, double c, double s)
        {
            for (int k = 0; k < U.Columns; k++)
            {
                double t1 = U[k, i];
                double t2 = U[k, j];
                U[k, i] = c * t1 + s * t2;
                U[k, j] = -s * t1 + c * t2;
            }

            for (int k = 0; k < V.Columns; k++)
            {
                double t1 = V[k, i];
                double t2 = V[k, j];
                V[k, i] = c * t1 + s * t2;
                V[k, j] = -s * t1 + c * t2;
            }
        }

        private static void SortSVD(Matrix U, Vector S, Matrix V)
        {
            int n = S.Dimension;
            for (int i = 0; i < n - 1; i++)
            {
                int maxIndex = i;
                double maxValue = S[i];

                for (int j = i + 1; j < n; j++)
                {
                    if (S[j] > maxValue)
                    {
                        maxIndex = j;
                        maxValue = S[j];
                    }
                }

                if (maxIndex != i)
                {
                    // Swap singular values
                    S[i] = S[maxIndex];
                    S[maxIndex] = maxValue;

                    // Swap columns in U and V
                    for (int k = 0; k < U.Rows; k++)
                        (U[k, i], U[k, maxIndex]) = (U[k, maxIndex], U[k, i]);
                    for (int k = 0; k < V.Rows; k++)
                        (V[k, i], V[k, maxIndex]) = (V[k, maxIndex], V[k, i]);
                }
            }
        }
        #endregion

        #region Helper Methods
        public static bool IsUpperTriangular(Matrix matrix, double tolerance = Epsilon)
        {
            for (int i = 1; i < matrix.Rows; i++)
                for (int j = 0; j < Math.Min(i, matrix.Columns); j++)
                    if (Math.Abs(matrix[i, j]) > tolerance)
                        return false;
            return true;
        }

        public static bool IsLowerTriangular(Matrix matrix, double tolerance = Epsilon)
        {
            for (int i = 0; i < matrix.Rows; i++)
                for (int j = i + 1; j < matrix.Columns; j++)
                    if (Math.Abs(matrix[i, j]) > tolerance)
                        return false;
            return true;
        }

        public static bool IsOrthogonal(Matrix matrix, double tolerance = Epsilon)
        {
            if (matrix.Rows != matrix.Columns)
                return false;

            Matrix product = matrix * matrix.Transpose();
            Matrix identity = Matrix.Identity(matrix.Rows);

            for (int i = 0; i < matrix.Rows; i++)
                for (int j = 0; j < matrix.Columns; j++)
                    if (Math.Abs(product[i, j] - identity[i, j]) > tolerance)
                        return false;

            return true;
        }
        #endregion
    }
}