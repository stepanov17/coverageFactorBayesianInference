import java.util.ArrayList;

public class CalcK_N {

    // number of measurements
    private final static int n = 4;

    private final double alpha, lambda;

    private final double x[], z[];
    private final int nx, nz;

    private final double p0 = 0.99;


    public CalcK_N(double alpha) {

        this.alpha = alpha;
        lambda = new Gamma().getLambda(alpha);
        System.out.println("// n = " + n + ", alpha = " + alpha + ", lambda = " + lambda + ", p0 = " + p0);

        double h = 0.001;
        final double kh = 1.0005;
        final double R  = 2.e3;

        ArrayList<Double> vals = new ArrayList<>();
        for (double v = h; v < R + 1.; v += h) {
            vals.add(v);
            h *= kh;
        }

        int nv = vals.size();

        x = new double[2 * nv + 1];
        z = new double[2 * nv + 1];

        int j = 0;
        for (int i = nv - 1; i >= 0; --i) {
            x[j] = -vals.get(i);
            z[j] = -vals.get(i);
            ++j;
        }
        x[j] = 0.;
        z[j] = 0.;
        ++j;
        for (int i = 0; i < nv; ++i) {
            x[j] = vals.get(i);
            z[j] = vals.get(i);
            ++j;
        }

        nx = x.length;
        nz = z.length;
    }


    double getL(double x[], double px[]) {

        double I = 1.;
        int npx = px.length;

        boolean found = false;
        double L = Double.NaN;

        double Iprev = 0., Lprev = 0.;

        for (int shift = 1; shift < npx / 2; ++shift) {

            double x1 = x[shift - 1], x2 = x[shift], p1 = px[shift - 1], p2 = px[shift];
            I -= 0.5 * (x2 - x1) * (p1 + p2);

            double tmp = x1;

            x1 = x[npx - shift - 1]; x2 = x[npx - shift]; p1 = px[npx - shift - 1]; p2 = px[npx - shift]; 
            I -= 0.5 * (x2 - x1) * (p1 + p2);

            if (!found && (I <= p0)) {

                L = 0.5 * (x2 - tmp);

                if (Math.abs(L - Lprev) > 0.) {

                    double kL = (L - Lprev) / (I - Iprev);
                    double bL = I - kL * L;
                    L = (p0 - bL) / kL;
                }

                return L;
            }

            Iprev = I;
            Lprev = 0.5 * (x2 - tmp);
        }
        return L;
    }

    private double getK(double gamma) {

        double px[] = new double[nx];

        double kn = 1. / (n - 1.);

        for (int ix = 0; ix < nx; ++ix) {

            double I = 0.;
            for (int iz = 1; iz < nz; ++iz) {

                double z1 = z[iz - 1];
                double z2 = z[iz];

                double tmp1 = kn * (gamma * z1 + x[ix]) * (gamma * z1 + x[ix]) + 1.;
                double tmp2 = kn * (gamma * z2 + x[ix]) * (gamma * z2 + x[ix]) + 1.;

                double f1 = Math.pow(tmp1, -0.5 * n) * Math.exp(-Math.pow(Math.abs(z1 / lambda), alpha));
                double f2 = Math.pow(tmp2, -0.5 * n) * Math.exp(-Math.pow(Math.abs(z2 / lambda), alpha));

                I += 0.5 * (z2 - z1) * (f1 + f2);
            }

            px[ix] = I;
        }

        // norming
        double s = 0;
        for (int i = 1; i < nx; ++i) {
            s += 0.5 * (x[i] - x[i - 1]) * (px[i - 1] + px[i]);
        }
        for (int i = 0; i < nx; ++i) { px[i] /= s; }

        double L = getL(x, px);
        return L / Math.sqrt((n - 1.) / (n - 3.) + gamma * gamma);
    }

    private static void calc(double alpha) {

        CalcK_N c = new CalcK_N(alpha);

        ArrayList<Double> Gamma = new ArrayList<>();
        Gamma.add(1.e-6);

        for (int i = 1; i <= 5; ++i)  { Gamma.add(0.1 * i);        }
        for (int i = 1; i <= 14; ++i) { Gamma.add(0.5 + 0.25 * i); }
        for (int i = 1; i <= 4; ++i)  { Gamma.add(4. + 0.5 * i);   }
        Gamma.add(7.);
        Gamma.add(8.);
        Gamma.add(9.);
        Gamma.add(10.);

        for (int i = 0; i < Gamma.size(); ++i) {
            double gamma = Gamma.get(i);
            double k = c.getK(gamma);
            System.out.printf("%6.2f  %.3f\n", gamma, k);
        }
    }

    public static void main(String args[]) {


        ArrayList<Double> Alpha = new ArrayList<>();
        Alpha.add(1.);
        Alpha.add(2.);
        Alpha.add(5.);
        Alpha.add(100.);

        int i = 1;
        for (Double alpha: Alpha) {
            System.out.println("data_" + i + " = [");
            calc(alpha);
            System.out.println("];\n\n");
            ++i;
        }
    }
}
