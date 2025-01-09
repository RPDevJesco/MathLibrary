namespace MathLibrary
{
        /// <summary>
    /// Implements Fast Fourier Transform algorithms
    /// </summary>
    public static class FFT
    {
        public static Complex[] Transform(Complex[] input)
        {
            int n = input.Length;

            // Check if length is power of 2
            if ((n & (n - 1)) != 0)
                throw new ArgumentException("Input length must be power of 2");

            Complex[] output = (Complex[])input.Clone();
            TransformRadix2(output, false);
            return output;
        }

        public static Complex[] InverseTransform(Complex[] input)
        {
            int n = input.Length;

            // Check if length is power of 2
            if ((n & (n - 1)) != 0)
                throw new ArgumentException("Input length must be power of 2");

            Complex[] output = (Complex[])input.Clone();
            TransformRadix2(output, true);

            // Scale the inverse transform
            for (int i = 0; i < n; i++)
                output[i] /= n;

            return output;
        }
        
        private static int ReverseBits(int x, int bits)
        {
            int result = 0;
            for (int i = 0; i < bits; i++)
            {
                result = (result << 1) | (x & 1);
                x >>= 1;
            }
            return result;
        }

        private static void TransformRadix2(Complex[] data, bool inverse)
        {
            int n = data.Length;
            int levels = (int)Math.Log2(n);

            // Bit reversal
            for (int i = 0; i < n; i++)
            {
                int j = ReverseBits(i, levels);
                if (j > i)
                {
                    var temp = data[i];
                    data[i] = data[j];
                    data[j] = temp;
                }
            }

            // Cooley-Tukey FFT
            for (int size = 2; size <= n; size *= 2)
            {
                double angle = 2 * Math.PI / size * (inverse ? 1 : -1);
                Complex wn = new(Math.Cos(angle), Math.Sin(angle));

                for (int i = 0; i < n; i += size)
                {
                    Complex w = Complex.One;
                    for (int j = 0; j < size / 2; j++)
                    {
                        int evenIndex = i + j;
                        int oddIndex = i + j + size / 2;

                        Complex even = data[evenIndex];
                        Complex odd = data[oddIndex];

                        data[evenIndex] = even + w * odd;
                        data[oddIndex] = even - w * odd;

                        w *= wn;
                    }
                }
            }
        }
    }
}