package blackdermantoy;

import java.lang.reflect.Array;
import java.util.Formatter;

import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

/*
 * Black Derman Toy (BDT) Interest Rate Model 
 * 
 * This routine demonstrates the implementation, calibration and application of this model 
 * to option embedded bond pricing.
 * 
 * BDT model was first published in 1990 as in this paper:
 * 
 * Black, F., Derman, E. and Toy, W. (1990) A one-factor model of interest rates and its 
 * application to treasury bond options, Financial Analysts Journal, 46(1), 33 9.
 * 
 * The constrcutor of the module takes a zero rate term structure and its volatility term strcuture
 * as inputs, generates a binomial interest tree. Based on the bond features, the module therefore
 * using the interest tree to price the bond. An example has been included.  
 * 
 * Flanagan's Scientific Java Package can be found here:
 * http://www.ee.ucl.ac.uk/~mflanaga/java/  
 * 
 * (C) Copyright 2010, Changwei Xiong (axcw@hotmail.com)
 * 
 */


public class BlackDermanToy {
	private final double dt;
	private final double dtrt;
	private final double[] yield;
	private final double[] beta;
	private double[] sigma;
	private double[] rate;
	private double[][][] ADT;
	private Minimisation nm = new Minimisation();
	private int step = 0;

	public BlackDermanToy(double[] yld, double[] yldvol, double terms[], double dt)
			throws Exception {
		if (dt <= 0 || yldvol == null || yld == null || yld.length != yldvol.length
				|| yld.length != terms.length) {
			throw new Exception("Invalid parameters for BDT model");
		}

		int size = (int) (terms[terms.length-1]/dt);
		this.yield = new double[size];
		this.beta = new double[size];
		
		
		/*
		 * use a modified Hermite Cubic Spine routine to interpolate between the key terms
		 * of the term structures
		 */
		HermiteCubicSpline interpyld = new HermiteCubicSpline(terms, yld);
		for (int i = 0; i < yield.length; i++) {
			yield[i] = interpyld.value(dt * (i + 1));
		}

		HermiteCubicSpline interpvol = new HermiteCubicSpline(terms, yldvol);
		for (int i = 0; i < beta.length; i++) {
			beta[i] = interpvol.value(dt * (i + 1));
		}

		this.dt = dt;
		this.dtrt = Math.sqrt(dt);

		this.sigma = fillZeros(new double[yield.length]);
		this.rate = fillZeros(new double[yield.length]);
		this.ADT = fillZeros(new double[3][yield.length + 1][yield.length + 1]);
	}

	/*
	 * Performs rounding for option term, maps it to the nearest tree node
	 */
	private int[] round_term(double[] t) {
		int[] k = new int[t.length];
		for (int i = 0; i < t.length; i++) {
			k[i] = (int) Math.round(t[i] / dt);
		}
		return k;
	}

	/*
	 * utility function for debug, visualize the tree
	 */
	private void print_tree(int depth) {
		Formatter f = new Formatter(System.out);
		for (int i = 0; i < depth; i++) {
			f.format("%d : ", i);
			for (int j = 0; j <= i; j++) {
				f.format("  %1.4e,", d2y(D(i, j), dt));
			}
			System.out.println();
		}
	}

	/*
	 * computes bond option value
	 */
	public double get_bond_price(double[][] cashflows, double[][] options) {
		final double[] ct = cashflows[0];
		final double[] cv = cashflows[1];
		final double[] ot = options[0];
		final double[] ov = options[1];
		final double[] ok = options[2];

		final int[] cT = round_term(ct);
		final int[] oT = round_term(ot);
		final int MAX = cT[cT.length - 1];

		build_tree(MAX);

		final double[] V = fillZeros(new double[MAX + 1]);
		for (int i = MAX, c = cT.length - 1, o = oT.length - 1; i >= 0; i--) {
			if (o >= 0 && oT[o] == i) {
				if (ok[o] > 0) {
					applyCeiling(V, ov[o], oT[o]);
				} else {
					applyFloor(V, ov[o], oT[o]);
				}
				o--;
			}
			if (c >= 0 && cT[c] == i) {
				addValues(V, cv[c], cT[c]);
				c--;
			}
			for (int j = 0; j < i; j++) {
				V[j] = (V[j] + V[j + 1]) * D(i - 1, j) / 2;
			}
		}
		return V[0];
	}

	public void build_tree(int depth) {
		for (int i = 0; i <= depth; i++) {
			compute_ith_period(i);
		}
		print_tree(depth);
	}

	/*
	 * Computes the i-th period state prices the results are stored in ADT tree
	 */
	private void compute_ith_period(int i) {
		if (i == 0) {
			ADT[0][0][0] = 1.0;
		} else if (i == 1) {
			double r = yield[0];
			double s = 0.0;
			double d = D(r, s, 0);
			ADT[0][1][0] = ADT[0][0][0] * 0.5 * d;
			ADT[0][1][1] = ADT[0][0][0] * 0.5 * d;

			ADT[1][1][1] = 1.0;
			ADT[2][1][0] = 1.0;
			rate[0] = r;
			sigma[0] = s;
			return;
		} else if (i <= yield.length) {
			double[] init = { Math.sqrt(rate[i - 2]), Math.sqrt(sigma[i - 2]) };
			try {
				this.step = i;

				nm.setNmax(5000);
				nm.nelderMead(GenObjective(), init, 1e-10);
				nm.getParamValues();

				nm.nelderMead(GenObjective(), init, 1e-10, 500);
				double[] results = nm.getParamValues();
				System.out.println("it=" + i + "\t" + nm.getNiter() + "\t" + nm.getMinimum() + "\t");
				rate[i - 1] = Math.pow(results[0], 2);
				sigma[i - 1] = Math.pow(results[1], 2);
			} catch (Exception e) {
				e.printStackTrace();
			}

			double[][] ad = compute_stateprices(rate[i - 1], sigma[i - 1], i);
			for (int n = 0; n < 3; n++) {
				for (int k = 0; k < i + 1; k++) {
					ADT[n][i][k] = ad[n][k];
				}
			}
		}
	}

	/*
	 * A generic method for initiating multi-level arrays to zero.
	 */
	private <T> T fillZeros(T a) {
		for (int i = 0; i < Array.getLength(a); i++) {// fill zeros
			if (a instanceof double[]) {
				Array.setDouble(a, i, 0.0);
			} else {
				fillZeros(Array.get(a, i));
			}
		}
		return a;
	}

	/*
	 * utility method for adding a scalar to a vactor
	 */
	private double[] addValues(double[] a, double val, int k) {
		for (int i = 0; i <= k; i++) {
			a[i] += val;
		}
		return a;
	}

	/*
	 * Ceiling price for option payoff
	 */
	private double[] applyCeiling(double[] a, double val, int k) {
		for (int i = 0; i <= k; i++) {
			if (a[i] > val) {
				a[i] = val;
			}
		}
		return a;
	}

	/*
	 * Flooring price for option payoff
	 */
	private double[] applyFloor(double[] a, double val, int k) {
		for (int i = 0; i <= k; i++) {
			if (a[i] < val) {
				a[i] = val;
			}
		}
		return a;
	}

	/*
	 * Generate objective function object for optimization
	 */
	private MinimisationFunction GenObjective() {
		return new MinimisationFunction() {
			public double function(double[] point) {
				return compute_residual(point);
			}
		};
	}

	private double D(int i, int j) {
		return D(rate[i], sigma[i], j);
	}

	/*
	 * yield to discount factor conversion using annual compounding
	 */
	private double D(double r, double s, int j) {
		double rj = r * Math.exp(2 * s * j * dtrt);
		return Math.pow(1 + rj, -dt);
		// return Math.exp(-rj * dt);
	}

	/*
	 * discount factor to yield conversion
	 */
	private double d2y(double d, double t) {
		return Math.pow(d, -1.0 / t) - 1.0;
		// return -Math.log(d) / t;
	}

	/*
	 * Computes objective function residual for optimization
	 */
	private double compute_residual(double in[]) {
		int i = step;
		double r = in[0] * in[0];
		double s = in[1] * in[1];
		double[][] ad = compute_stateprices(r, s, i);

		double yi = d2y(sum(ad[0]), dt * i--);
		double yu = d2y(sum(ad[1]), dt * i);
		double yd = d2y(sum(ad[2]), dt * i);

		double dy = yield[i] - yi;
		// double db = beta[i] - Math.abs(yu - yd) / 2 / dtrt;
		double db = beta[i] - ((yd == 0) ? 0 : (Math.log(yu / yd) / 2)) / dtrt;
		return Math.sqrt(dy * dy + db * db);
	}

	private double sum(double[] a) {
		double tmp = 0.0;
		for (int i = 0; i < a.length; i++) {
			tmp += a[i];
		}
		return tmp;
	}

	/*
	 * Compute the stateprice trees as defined in BDT model.
	 */
	private double[][] compute_stateprices(double r, double s, int i) {
		double[][] ad = fillZeros(new double[3][i + 1]);
		for (int m = 0; m <= i; m++) {
			if (m == 0) {
				double d = 0.5 * D(r, s, 0);
				ad[0][m] = ADT[0][i - 1][m] * d;
				ad[2][m] = ADT[2][i - 1][m] * d;
			} else if (m > 0 && m < i) {
				double d_d = 0.5 * D(r, s, m - 1);
				double d_u = 0.5 * D(r, s, m);
				for (int n = 0; n < 3; n++) {
					ad[n][m] = ADT[n][i - 1][m - 1] * d_d + ADT[n][i - 1][m] * d_u;
				}
			} else {// m == i
				double d = 0.5 * D(r, s, m - 1);
				ad[0][m] = ADT[0][i - 1][m - 1] * d;
				ad[1][m] = ADT[1][i - 1][m - 1] * d;
			}
		}
		return ad;
	}

	/**
	 * Several test cases here.
	 */
	public static void main(String[] args) {

		/*
		 * Yield and Yield Volatility term structure 
		 * total 21 terms from overnight, 1 week, ..., up to 50 years
		 */
		double[] yields = { 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425,
				0.045, 0.0475, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0625, 0.065, 0.0675, 0.07 };
		double[] volatilities = { 0.6, 0.57, 0.541, 0.514, 0.488, 0.464, 0.441, 0.419, 0.398,
				0.378, 0.374, 0.341, 0.324, 0.308, 0.292, 0.277, 0.264, 0.250, 0.238, 0.226, 0.215 };
		double[] terms = { 0.0, 0.0192, 0.0385, 0.833, 0.167, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7,
				8, 10, 15, 20, 30, 50 };
		double deltaT = 0.2;

		try {
			BlackDermanToy bdt = new BlackDermanToy(yields, volatilities, terms, deltaT);

			/*
			 * An example call option embedded bond. 3-year annually paid fixed
			 * coupon bond with annual coupon 5.25. The call option is at the
			 * end of 2nd year.
			 */
			double[][] cashflows = { { 1, 2, 3 }, { 5.25, 5.25, 105.25 } };
			// {terms[]/* years*/, cashflows[]}
			double[][] options = { { 2 }, { 100 }, { 1 } };
			// {terms[], strikes[], types[] <- (1 for call and -1 for put)}
			System.out.println("\n bond price = " + bdt.get_bond_price(cashflows, options));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}