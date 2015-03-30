package blackdermantoy;

/*
 * A Modified Hermite Cubic Spline Interpolator
 * 
 * (C) Copyright 2010, Changwei Xiong (axcw@hotmail.com)
 * 
 */

public class HermiteCubicSpline {
	private final double[][] H = { { 2, -3, 0, 1 }, { 1, -2, 1, 0 },
			{ -2, 3, 0, 0 }, { 1, -1, 0, 0 } };
	private double[] x, y, dr;
	private int N;

	public HermiteCubicSpline(double[] x, double[] y) throws Exception {
		if (x.length < 2 || x.length != y.length) {
			throw new Exception("Invalid input data");
		}
		this.x = x;
		this.y = y;
		this.N = this.x.length;

		double[] delta = new double[N - 1];// N-1 by 1
		for (int i = 0; i < N - 1; i++) {
			delta[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
		}

		dr = new double[N];
		dr[0] = 0;
		dr[N - 1] = 0;
		for (int i = 1; i < N - 1; i++) {
			if (delta[i] != 0) {
				dr[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
			} else {
				delta[i] = 0;
			}
		}

		double[] a_ = new double[N - 1];// N-1 by 1
		double[] b_ = new double[N - 1];// N-1 by 1
		for (int i = 0; i < N - 1; i++) {

			if (delta[i] != 0) {
				a_[i] = dr[i] / delta[i];
				b_[i] = dr[i + 1] / delta[i];
			} else {
				a_[i] = 0;
				b_[i] = 0;
			}

			double a = a_[i];
			double b = b_[i];
			if (a + b - 2 > 0
					&& 2 * a + b - 3 > 0
					&& a + 2 * b - 3 > 0
					&& a - (2 * a + b - 3) * (2 * a + b - 3) / (a + b - 2)
							/ 3.0 < 0) {
				double factor = 3.0 / Math.sqrt(a * a + b * b);
				a_[i] = a * factor;
				b_[i] = b * factor;
				dr[i] = a_[i] * delta[i];
				dr[i + 1] = b_[i] * delta[i];
			}
		}
	}

	private int index(double s) {
		int stt = 0;
		int end = N - 1;
		while (end - stt > 1) {
			int mid = (stt + end) / 2;
			if (s >= x[mid]) {
				stt = mid;
			} else {
				end = mid;
			}
		}
		return stt;
	}

	public double value(double s) {
		int k = index(s);
		double[] tv = new double[4];
		double[] rv = new double[4];
		double t = (s - x[k]) / (x[k + 1] - x[k]);
		tv[0] = t * t * t;
		tv[1] = t * t;
		tv[2] = t;
		tv[3] = 1;

		for (int i = 0; i < 4; i++) {
			rv[i] = 0.0;
			for (int j = 0; j < 4; j++) {
				rv[i] += H[i][j] * tv[j];
			}
		}

		double dx = x[k + 1] - x[k];
		return rv[0] * y[k] + rv[2] * y[k + 1] + rv[1] * dx * dr[k] + rv[3]
				* dx * dr[k + 1];
	}

	public static void main(String[] args) throws Exception {
		double[] x = { 0, 0.17, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 10, 15,
				20, 30, 40, 50 };
		double[] y = { 1.0569, 1.3101, 1.3294, 1.4141, 1.5228, 1.612, 1.9846,
				2.335, 2.5441, 2.7584, 2.9554, 3.0824, 3.4954, 4.039, 4.1689,
				4.3703, 4.5968, 4.786 };
		HermiteCubicSpline hcs = new HermiteCubicSpline(x, y);
		for (int i = 0; i <= 100; i++) {
			double xx = 0.5 * i;
			System.out.println(xx + "  " + hcs.value(xx));
			//49.5  4.785121943749999
		}
	}
}