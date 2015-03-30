{\rtf1\ansi\ansicpg1252\cocoartf1265\cocoasubrtf210
{\fonttbl\f0\froman\fcharset0 Times-Roman;\f1\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue233;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720

\f0\fs24 \cf0 Below are two examples that demonstrate the application of {\field{\*\fldinst{HYPERLINK "http://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo"}}{\fldrslt 
\i \cf2 \ul \ulc2 MCMC method}} in model parameter estimation. The first example is a standard Linear Regression model and the second is a simple ARMAX(0,0,1)/GARCH(1,1) model (the same model shown above).\
The Linear Regression model is in the form:\
y = a + b*x1 + c*x2 + e,\'a0 \'a0 \'a0 e ~ Normal(0, s)\
In this model, the parameters to be estimated are 
\i a
\i0 , "b", "c", and "s". {\field{\*\fldinst{HYPERLINK "http://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm"}}{\fldrslt 
\i \cf2 \ul \ulc2 Metropolis-Hastings method}} is used to generate the MCMC sample sequencies. For simplicity, the prior distribution of the parameters are assumed to be flat (uniformly distributed). The proposed density utilizes a Normal distribution, in which case the original Metropolis-Hastings method has been reduced to a Random-Walk Metropolis Hastings method. Please refer to the {\field{\*\fldinst{HYPERLINK "http://www.cs.utah.edu/~cxiong/Files/Codes/MCMC_LinearRegress.m"}}{\fldrslt \cf2 \ul \ulc2 Matlab script}} for the details of the implementation. For a simulated case using "a = 1.0, b = 2.0, c = 3.0, s = 1.2", the method estimates the parameters and generates a posterior distribution for the parameter estimations. The histograms are shown below:\
\pard\pardeftab720

\f1 \cf0 
\f0 \
\pard\pardeftab720
\cf0 The ARMAX(0,0,1)/GARCH(1,1) model, which was previously estimated by MLE method, is in the form:\
y(t) = a + b * x(t) + u(t)\
u(t) ~ Normal(0, h(t))\
h(t) = c + d * u(t-1)^2 + e * h(t-1)\
In this model, the parameters to be estimated are "a", "b", "c", "d" and "e". Random-Walk Metropolis Hastings method is again used to generate the MCMC sample sequencies. The prior distribution of the parameters is also assumed to be flat. Please refer to the {\field{\*\fldinst{HYPERLINK "http://www.cs.utah.edu/~cxiong/Files/Codes/MCMC_ARMAX_GARCH.m"}}{\fldrslt \cf2 \ul \ulc2 Matlab script}} for the details of the implementation. It's slightly different from the case of Linear Regression, the parameters in GARCH model are constrained. To handle this issue, any sampled values that exceed the parameter bounds are assumed with zero probability and therefore discarded. For a simulated case using "a = 0.2, b = 0.9, c = 0.2, d = 0.3, e = 0.4", the method estimates the parameters and generates a posterior distribution for the parameter estimations. Note that estimation variances of GARCH parameters are much larger than those of regression parameters, which indicates that the Garch parameters are highly sensitive to input data. The histograms are shown below:\
\pard\pardeftab720

\f1 \cf0 }