/**
 * Created by danya on 22.12.16.
 */
public class Consts {
    /**
     * Gravitation's parameter of Sun
     */
    private final double xi;
    private final double eps;
    private final double prescision;

    Consts() {
        xi = 0.017202;
        eps = (23 + (26 + 21.4119 / 60.0) / 60.0) * Math.PI / 180.0;
        prescision = 0.000_000_000_1;
    }

    public double getXi() {
        return xi;
    }

    public double getEps() {
        return eps;
    }

    public double getPrescision() {
        return prescision;
    }
}
