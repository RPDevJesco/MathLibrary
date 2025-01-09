namespace MathLibrary
{
    /// <summary>
    /// Represents a quaternion for 3D rotations and calculations
    /// </summary>
    public class Quaternion
    {
        public double W { get; }
        public double X { get; }
        public double Y { get; }
        public double Z { get; }

        public Quaternion(double w, double x, double y, double z)
        {
            W = w;
            X = x;
            Y = y;
            Z = z;
        }

        public static Quaternion FromAxisAngle(Vector axis, double angle)
        {
            if (axis.Dimension != 3)
                throw new ArgumentException("Axis must be a 3D vector");

            double halfAngle = angle / 2;
            double sin = Math.Sin(halfAngle);
            axis = axis.Normalize();

            return new Quaternion(
                Math.Cos(halfAngle),
                axis[0] * sin,
                axis[1] * sin,
                axis[2] * sin
            );
        }

        public static Quaternion FromEulerAngles(double pitch, double yaw, double roll)
        {
            double cy = Math.Cos(yaw * 0.5);
            double sy = Math.Sin(yaw * 0.5);
            double cp = Math.Cos(pitch * 0.5);
            double sp = Math.Sin(pitch * 0.5);
            double cr = Math.Cos(roll * 0.5);
            double sr = Math.Sin(roll * 0.5);

            return new Quaternion(
                cr * cp * cy + sr * sp * sy,
                sr * cp * cy - cr * sp * sy,
                cr * sp * cy + sr * cp * sy,
                cr * cp * sy - sr * sp * cy
            );
        }

        public static Quaternion operator *(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.W * b.W - a.X * b.X - a.Y * b.Y - a.Z * b.Z,
                a.W * b.X + a.X * b.W + a.Y * b.Z - a.Z * b.Y,
                a.W * b.Y - a.X * b.Z + a.Y * b.W + a.Z * b.X,
                a.W * b.Z + a.X * b.Y - a.Y * b.X + a.Z * b.W
            );
        }

        public Vector Rotate(Vector v)
        {
            if (v.Dimension != 3)
                throw new ArgumentException("Vector must be 3D for rotation");

            Quaternion p = new(0, v[0], v[1], v[2]);
            Quaternion q = this;
            Quaternion qInv = new(q.W, -q.X, -q.Y, -q.Z);
            Quaternion result = q * p * qInv;

            return new Vector(result.X, result.Y, result.Z);
        }
        
        public override string ToString()
        {
            var terms = new List<string>();
    
            // Add real part if non-zero
            if (W != 0 || (X == 0 && Y == 0 && Z == 0))
                terms.Add(W.ToString("F3"));
        
            // Add vector parts if non-zero
            if (X != 0) terms.Add($"{X:F3}i");
            if (Y != 0) terms.Add($"{Y:F3}j");
            if (Z != 0) terms.Add($"{Z:F3}k");
    
            return string.Join(" + ", terms).Replace("+ -", "- ");
        }
    }
}