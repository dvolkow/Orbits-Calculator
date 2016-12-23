/**
 * Created by danya on 22.12.16.
 */
public class Calc {

        private class Result {

                private double a;
                private double e;
                private double i;
                private double omega;
                private double bOmega;

                Result(final double a, final double e, final double i, final double omega, final double bOmega) {
                        this.a = a;
                        this.e = e;
                        this.i = i;
                        this.omega = omega;
                        this.bOmega = bOmega;
                }

                public double getA() {
                        return a;
                }

                public double getE() {
                        return e;
                }

                public double getI() {
                        return i;
                }

                public double getOmega() {
                        return omega;
                }

                public double getBOmega() {
                        return bOmega;
                }

        }

        private Result result;

        Calc() {
                Coordinats coordinats = new Coordinats();
                Consts consts = new Consts();

                double[] t = new double[3];
                t[1] = coordinats.getT0()[1] - coordinats.getT0()[0];
                t[2] = coordinats.getT0()[2] - coordinats.getT0()[1];
                t[0] = t[1] + t[2];

                double[] lambda = new double[3];
                double[] mu = new double[3];
                double[] nu = new double[3];

                for (int i = 0; i < 3; ++i) {
                        lambda[i] = Math.cos(coordinats.getBody()[i].getRa()) * Math.cos(coordinats.getBody()[i].getDec());
                        mu[i] = Math.cos(coordinats.getBody()[i].getDec()) * Math.sin(coordinats.getBody()[i].getRa());
                        nu[i] = Math.sin(coordinats.getBody()[i].getDec());
                }

                double[] tau = new double[3];
                for (int i = 0; i < 3; ++i) {
                        tau[i] = t[i] * consts.getXi();
                }

                double n01 = tau[2] / tau[0];
                double n02 = tau[1] / tau[0];

                double c1 = 1.0/6.0 * tau[1] * tau[2] * (1 + n01);
                double c2 = 1.0/6.0 * tau[1] * tau[2] * (1 + n02);

                double lambda13 = mu[0] * nu[2] - mu[2] * nu[0];
                double mu13 = nu[0] * lambda[2] - nu[2]*lambda[0];
                double nu13 = lambda[0]*mu[2] - lambda[2]*mu[0];

                double[][] D = new double[3][3];
                D[0] = new double[] {lambda[1], lambda[0], lambda[2]};
                D[1] = new double[] {mu[1], mu[0], mu[2]};
                D[2] = new double[] {nu[1], nu[0], nu[2]};

                double det = detMatrix3_3(D);

                double[] U = new double[3];
                for (int i = 0; i < 3; ++i) {
                        U[i] = coordinats.getSun()[i].getX() * lambda13 +
                               coordinats.getSun()[i].getY() * mu13 +
                               coordinats.getSun()[i].getZ() * nu13;
                }

                double P = (U[1] - n01 * U[0] - n02 * U[2]) / det;
                double Q = (c1 * U[0] + c2 * U[2]) / det;
                double C = -(lambda[1] * coordinats.getSun()[1].getX() + mu[1] * coordinats.getSun()[1].getY() +
                        nu[1] * coordinats.getSun()[1].getZ());

                double R2 = Math.pow(coordinats.getSun()[1].getX(), 2) + Math.pow(coordinats.getSun()[1].getY(), 2) +
                        Math.pow(coordinats.getSun()[1].getZ(), 2);

                double[] ro = new double[4];
                ro[2] = 0;
                ro[0]  = -1;
                double r2 = 0;

                while (Math.abs(ro[2] - ro[0]) > consts.getPrescision()) {
                        ro[0] = ro[2];
                        r2 = Math.sqrt(ro[0] * ro[0] + 2 * C * ro[0] + R2);
                        ro[2] = ro[0] - (ro[0] - P + Q / (r2 * r2 * r2)) / (1 - 3 * Q * (ro[0] + C) / Math.pow(r2, 5));
                }

                double n1 = n01 + c1 / Math.pow(r2, 3);
                double n2 = n02 + c2 / Math.pow(r2, 3);


                if (lambda13 >= mu13 && lambda13 >= nu13) {

                        ro[1] = (ro[2] * (mu[1] * nu[2] - nu[1] * mu[2]) - (nu[2] * coordinats.getSun()[1].getY() -
                                mu[2] * coordinats.getSun()[1].getZ()) + n1 * (nu[2] * coordinats.getSun()[0].getY() -
                                mu[2] * coordinats.getSun()[0].getZ()) + n2 * (nu[2] * coordinats.getSun()[2].getY() -
                                mu[2] * coordinats.getSun()[2].getZ())) / (n1 * lambda13);

                } else if (mu13 >= lambda13 && mu13 >= nu13) {

                        ro[1] = (ro[2] * (nu[1] * lambda[2] - lambda[1] * nu[2]) - (lambda[2] * coordinats.getSun()[1].getZ() -
                                nu[2] * coordinats.getSun()[1].getX()) + n1 * (lambda[2] * coordinats.getSun()[0].getZ() -
                                nu[2] * coordinats.getSun()[0].getX()) + n2 * (lambda[2] * coordinats.getSun()[2].getZ() -
                                nu[2] * coordinats.getSun()[2].getX())) / (n1 * mu13);

                } else if (nu13 >= lambda13 && nu13 >= mu13) {

                        ro[1] = (ro[2] * (lambda[1] * mu[2] - mu[1] * lambda[2]) - (mu[2] * coordinats.getSun()[1].getX() -
                                lambda[2] * coordinats.getSun()[1].getY()) + n1 * (mu[2] * coordinats.getSun()[0].getX() -
                                lambda[2] * coordinats.getSun()[0].getY()) + n2 * (mu[2] * coordinats.getSun()[2].getX() -
                                lambda[2] * coordinats.getSun()[2].getY())) / (n1 * nu13);

                }

                if (lambda[2] >= mu[2] && (lambda[2] >= nu[2])) {

                        ro[3] = (lambda[1] * ro[2] - n1 * lambda[0] * ro[1] - coordinats.getSun()[1].getX() +
                                n1 * coordinats.getSun()[0].getX() + n2 * coordinats.getSun()[2].getX()) / (n2 * lambda[2]);

                } else if (mu[2] >= nu[2] && mu[2] >= lambda[2]) {

                        ro[3] = (mu[1] * ro[2] - n1 * mu[0] * ro[1] - coordinats.getSun()[1].getY() +
                                n1 * coordinats.getSun()[0].getY() + n2 * coordinats.getSun()[2].getY()) / (n2 * mu[2]);

                } else if (nu[2] >= mu[2] && nu[2] >= lambda[2]) {

                        ro[3] = (nu[1] * ro[2] - n1 * nu[0] * ro[1] - coordinats.getSun()[1].getZ() +
                                n1 * coordinats.getSun()[0].getZ() + n2 * coordinats.getSun()[2].getZ()) / (n2 * nu[2]);

                }

                double[] y = new double[3];
                double[] z = new double[3];
                double[] x = new double[3];
                double[] r = new double[3];

                for (int i = 0; i < 3; ++i) {
                        x[i] = lambda[i] * ro[i] - coordinats.getSun()[i].getX();
                        y[i] = mu[i] * ro[i] - coordinats.getSun()[i].getY();
                        z[i] = nu[i] * ro[i] - coordinats.getSun()[i].getZ();
                }

                r[0] = Math.sqrt(Math.pow(x[0], 2) + Math.pow(y[0], 2) + Math.pow(z[0], 2));
                r[1] = r2;
                r[2] = Math.sqrt(Math.pow(x[2], 2) + Math.pow(y[2], 2) + Math.pow(z[2], 2));

                double hi = Math.sqrt(2 * (r[0] * r[2] + x[0] * x[2] + y[0] * y[2] + z[0] * z[2]));

                double alpha = 1 - hi / (r[0] + r[2]);
                double beta = 4 * Math.pow(tau[0], 2) / (3 * Math.pow(r[0] + r[2], 3));
                double etha = 1 + beta * (1 + 2.4 * alpha - 1.1 * beta);

                double sigma = (x[0] * x[2] + y[0] * y[2] + z[0] * z[2]) / Math.pow(r[0], 2);

                double x0 = x[2] - sigma * x[0];
                double y0 = y[2] - sigma * y[0];
                double z0 = z[2] - sigma * z[0];
                double r0 = Math.sqrt(Math.pow(x0, 2) + Math.pow(y0, 2) + Math.pow(z0, 2));
                double d_Nu = Math.atan2(r0 / r[2], sigma * r[0] / r[2]);
                double p = Math.pow(r0 * r[0] * etha / tau[0], 2);

                double q1 = p / r[0] - 1;
                double q2 = p / r[2] - 1;

                double Nu_1 = Math.atan2((q1 * Math.cos(d_Nu) - q2) / Math.sin(d_Nu), q1);
//                double Nu_2 = Nu_1 + d_Nu;
                double e = q1 / Math.cos(Nu_1);
                double a = p / (1 - e * e);

                double P_x = x[0] * Math.cos(Nu_1) / r[0] - x0 * Math.sin(Nu_1) / r0;
                double P_y = y[0] * Math.cos(Nu_1) / r[0] - y0 * Math.sin(Nu_1) / r0;
                double P_z = z[0] * Math.cos(Nu_1) / r[0] - z0 * Math.sin(Nu_1) / r0;

                double Q_x = x[0] * Math.sin(Nu_1) / r[0] + x0 * Math.cos(Nu_1) / r0;
                double Q_y = y[0] * Math.sin(Nu_1) / r[0] + y0 * Math.cos(Nu_1) / r0;
                double Q_z = z[0] * Math.sin(Nu_1) / r[0] + z0 * Math.cos(Nu_1) / r0;

                double omega = Math.atan2(P_z * Math.cos(consts.getEps()) - P_y * Math.sin(consts.getEps()),
                        Q_z * Math.cos(consts.getEps()) - Q_y * Math.sin(consts.getEps()));

                if (omega < 0) {
                        omega += 2 * Math.PI;
                }

                double i = Math.asin((P_z * Math.cos(consts.getEps()) - P_y * Math.sin(consts.getEps())) / Math.sin(omega));
                double OMEGA = Math.atan2((P_y * Math.cos(omega) - Q_y * Math.sin(omega)) /
                        Math.cos(consts.getEps()), P_x * Math.cos(omega) - Q_x * Math.sin(omega));

                if (OMEGA < 0) {
                        OMEGA += 2 * Math.PI;
                }

                this.result = new Result(a, e, i, omega, OMEGA);
        }

        /**
        * Show result of calculation
        */
        public void getResult() {
                System.out.println("a = " + result.getA());
                System.out.println("e = " + result.getE());
                System.out.println("i = " + Math.toDegrees(result.getI()));
                System.out.println("omega = " + Math.toDegrees(result.getOmega()));
                System.out.println("Omega = " + Math.toDegrees(result.getBOmega()));
        }

        /**
        * Calculate determinant of matrix with dimensions 3 : 3
        * @param matrix have dimensions [3][3]
        * @return double determinant
        */
        private double detMatrix3_3(final double[][] matrix) {
                return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) - matrix[0][1] *
                        (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) + matrix[0][2] * (matrix[1][0] *
                        matrix[2][1] - matrix[1][1] * matrix[2][0]);
        }

}
