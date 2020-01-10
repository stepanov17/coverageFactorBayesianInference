
// coverage probability for cov. factor value K = 2

import java.util.ArrayList;

public class CalcP_U {

    private final static int n = 4; // number of measurements

    private final double alpha, lambda;

    private final double x[], z[];
    private final int nx, nz;


    public CalcP_U(double alpha) {

        this.alpha = alpha;
        lambda = new Gamma().getLambda(alpha);
        System.out.println("// n = " + n + ", alpha = " + alpha + ", lambda = " + lambda);

        double h  = 0.001;
        final double kh = 1.002;
        final double R  = 1.e3;

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

    private double f(double x, double z, double gamma) {

        double tmp = 2. * Math.abs(gamma * z + x) / Math.sqrt(n) + 1.;
        return Math.pow(tmp, -n) * Math.exp(-Math.pow(Math.abs(z / lambda), alpha));
    }

    private double getNormFactor(double gamma) {

        double px[] = new double[nx];

        for (int ix = 0; ix < nx; ++ix) {

            double I = 0.;
            for (int iz = 1; iz < nz; ++iz) {

                double z1 = z[iz - 1];
                double z2 = z[iz];
                I += 0.5 * (z2 - z1) * (f(x[ix], z1, gamma) + f(x[ix], z2, gamma));
            }

            px[ix] = I;
        }

        double I = 0.;
        for (int i = 1; i < nx; ++i) {
            I += 0.5 * (x[i] - x[i - 1]) * (px[i - 1] + px[i]);
        }
        return 1. / I;
    }

    private double Iz(double x, double gamma) {

        double I = 0.;
        for (int iz = 1; iz < nz; ++iz) {

            double z1 = z[iz - 1];
            double z2 = z[iz];
            I += 0.5 * (z2 - z1) * (f(x, z1, gamma) + f(x, z2, gamma));
        }
        return I;
    }

    private double getP(double gamma, double K, double fact) {

        double R = K * Math.sqrt(0.5 * n / ((n - 2.) * (n - 3.)) + gamma * gamma);
        double dx = 0.01 * R;

        double p_prev = Iz(-R, gamma);
        double I = 0.;

        for (double x = -R + dx; x < R + 0.1 * dx; x += dx) {
            double p = Iz(x, gamma);
            I += 0.5 * dx * (p + p_prev);
            p_prev = p;
        }

        return fact * I;
    }


    private static void calc(double alpha) {

        CalcP_U c = new CalcP_U(alpha);

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
            double fact = c.getNormFactor(gamma);
            double p2 = c.getP(gamma, 2., fact);
            System.out.printf("%6.2f  %.3f\n", gamma, p2);
        }
    }

    public static void main(String[] args) {

        System.out.println("// unif");

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
