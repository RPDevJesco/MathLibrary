namespace MathLibrary
{
    /// <summary>
    /// Provides statistical functions for data analysis
    /// </summary>
    public static class Statistics
    {
        public static double Mean(IEnumerable<double> values)
        {
            if (!values.Any())
                throw new ArgumentException("Cannot compute mean of empty sequence");

            double sum = 0;
            int count = 0;
            foreach (var value in values)
            {
                checked
                {
                    sum += value;
                    count++;
                }
            }

            return sum / count;
        }

        public static double Variance(IEnumerable<double> values, bool sample = true)
        {
            if (!values.Any())
                throw new ArgumentException("Cannot compute variance of empty sequence");

            double mean = Mean(values);
            double sumSquaredDiffs = 0;
            int count = 0;

            foreach (var value in values)
            {
                double diff = value - mean;
                sumSquaredDiffs += diff * diff;
                count++;
            }

            return sumSquaredDiffs / (sample ? count - 1 : count);
        }

        public static double StandardDeviation(IEnumerable<double> values, bool sample = true)
        {
            return Math.Sqrt(Variance(values, sample));
        }

        public static double Correlation(IEnumerable<double> x, IEnumerable<double> y)
        {
            var xArray = x.ToArray();
            var yArray = y.ToArray();

            if (xArray.Length != yArray.Length)
                throw new ArgumentException("Sequences must have same length");

            if (!xArray.Any())
                throw new ArgumentException("Cannot compute correlation of empty sequences");

            double meanX = Mean(xArray);
            double meanY = Mean(yArray);
            double stdDevX = StandardDeviation(xArray);
            double stdDevY = StandardDeviation(yArray);

            double covariance = 0;
            for (int i = 0; i < xArray.Length; i++)
                covariance += (xArray[i] - meanX) * (yArray[i] - meanY);

            covariance /= xArray.Length - 1;
            return covariance / (stdDevX * stdDevY);
        }

        public static double Skewness(IEnumerable<double> values)
        {
            var data = values.ToList();
            double mean = Mean(data);
            double std = StandardDeviation(data);

            double sum = 0;
            foreach (var value in data)
            {
                double diff = (value - mean) / std;
                sum += diff * diff * diff;
            }

            return sum / (data.Count - 1);
        }

        public static double Kurtosis(IEnumerable<double> values)
        {
            var data = values.ToList();
            double mean = Mean(data);
            double std = StandardDeviation(data);

            double sum = 0;
            foreach (var value in data)
            {
                double diff = (value - mean) / std;
                sum += diff * diff * diff * diff;
            }

            return sum / (data.Count - 1) - 3; // Excess kurtosis
        }
    }
}