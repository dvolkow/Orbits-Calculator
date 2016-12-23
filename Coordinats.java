/**
 * Created by danya on 22.12.16.
 */
public class Coordinats {

    /**
     * Struct for equatoral coordinats
     */
    public class Eqats {
        /**
         * Right Ascention in radians
         */
        private double Ra;
        /**
         * Declination in radians
         */
        private double Dec;

        public double getRa() {
            return Ra;
        }

        public double getDec() {
            return Dec;
        }

        /**
         * Setter for Ra. Attention!
         * @param Ra in HOURs!
         */
        public void setRa(final double Ra) {
            this.Ra = Math.toRadians(Ra * 15);
        }

        /**
         * Setter for Dec. Attention!
         * @param Dec in degress!
         */
        public void setDec(final double Dec) {
            this.Dec = Math.toRadians(Dec);
        }

        Eqats(final double Ra, final double Dec) {
            setRa(Ra);
            setDec(Dec);
        }
    }

    /**
     * Setting Ulian Times
     */
    private double[] t0;

    /**
     * Coordinats of Sun
     */
    private Point[] sun;

    /**
     * Coordinats of Body
     */
    private Eqats[] body;

    /**
     * Default constructor
     */
    Coordinats() {
        t0 = new double[3];
        t0[0]=2457640.5;
        t0[1]=2457660.5;
        t0[2]=2457680.5;

        sun = new Point[3];
        sun[0] = new Point(-0.97945028, 0.21542557, 0.09339731);
        sun[1] = new Point(-0.99616354, -0.09686955, -0.04198842);
        sun[2] = new Point(-0.89666121, -0.39778502, -0.17244123);

        body = new Eqats[3];
        body[0] = new Eqats(16.28038, -21.46886);
        body[1] = new Eqats(16.96938, -23.14182);
        body[2] = new Eqats(17.56053, -24.06582);

    }

    public double[] getT0() {
        return t0;
    }

    public Point[] getSun() {
        return sun;
    }

    public Eqats[] getBody() {
        return body;
    }
}
