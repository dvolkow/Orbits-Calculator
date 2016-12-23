/**
 * Created by danya on 22.12.16.
 */
public class Point {
        private double x;
        private double y;
        private double z;

        public double getX() {
                return x;
        }

        public double getY() {
                return y;
        }

        public double getZ() {
                return z;
        }

        Point(final double X, final double Y, final double Z) {
                x = X;
                y = Y;
                z = Z;
        }
}
