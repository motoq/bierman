Bierman
=======

This project is written as a supplement to Gerald J. Bierman's "Factorization
Methods for Discrete Sequential Estimation".  Selected algorithms are be
illustrated in Matlab (actually developed with Octave).  Further documentation
will be developed as examples are finished.

Estimation of a static state is first illustrated.  Range only trackers are
used to estimate the location of an object within a boxed volume.

driver_01.m
```
Compares estimation of a static state using full batch, sequential batch,
stabilized Kalman, Potter mechanization, and U-D sequential methods.
```

driver_02.m
```
Performs speed comparisons of the stabilized Kalman, Potter
mechanization, and U-D sequential estimation methods.
```

driver_householder.m - detour
```
Compares  A WLS solution of the normal equations with one using
Householder triangularization and QR decomposition.
```

driver_lls_hhcurvefit.m - detour
```
Illustrates sequential estimation via Householder SRIF for a linear
problem.
```

driver_lls_qrcurvefit.m - detour
```
Illustrates sequential estimation via QR SRIF for a linear problem.
```

driver_03.m
```
Estimation of a dynamic state - simple wiffle ball in a room trajectory
```

