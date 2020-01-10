public class Gamma {

    private static final double P[] = {
             1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091,
            -1.231739572450155, 1.208650973866179e-3, -5.395239384953e-6};

    // gamma function approximation
    public double get(double x) {

        double v = P[0];
        for (int i = 1; i < P.length; ++i) {
            v += P[i] / (x + (double) i);
        }

        v *= (Math.sqrt(2. * Math.PI) / x);
        v *= Math.pow(x + 5.5, x + 0.5);
        v *= Math.exp(-(x + 5.5));

        return v;
    }

    // exp. distribution family parameter (depending on alpha)
    public double getLambda(double alpha) {
        return Math.sqrt(get(1. / alpha) / get(3. / alpha));
    }
}
