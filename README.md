Bierman
=======

This project is written as a supplement to Gerald J. Bierman's "Factorization
Methods for Discrete Sequential Estimation".  Selected algorithms are
illustrated in Matlab (actually developed with Octave).  The emphasis is
on illustrating the use of these algorithms for problems requiring
"linearization" along with the "extended" filtering method.

Estimation of a static state is first illustrated.  Range only trackers
are used to estimate the location of an object within a boxed volume.
Estimation of a dynamics state follows by tracking movement of a wiffle
ball within a room using the same range only trackers.

driver_01.m
```
Compares estimation of a static state using full batch, sequential batch,
stabilized Kalman, Potter mechanization, and U-D sequential methods.
```

driver_02.m
```
Performs speed comparisons of the stabilized Kalman, Potter
mechanization, and U-D sequential estimation methods.  Percent
containment is also output.
```

driver_02b.m
```
Performs speed comparisons of the stabilized the U-D and SRIF
methods.  Percent containment is also outout.
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
Illustrates sequential estimation via a QR SRIF for a linear problem.
```

driver_03.m
```
Estimation of a dynamic state - simple wiffle ball in a room trajectory

This driver script is the first of this series of examples illustrating
filtering of a dynamic scenario.  A simple drag free trajectory is created
along with simulated range observations (using range trackers from previous
examples) that are subject to only measurement noise.  No bias effects
have been added.  The goal is to implement a simple filter (observation model
and dynamic model) that outperforms the WLS method (observation model only).
```

driver_04.m
```
Adds drag to the truth model for driver_03 while leaving it out of the
filter model.  Illustrates divergence in the filter and how to fix it
with the inclusion of process noise.
```

