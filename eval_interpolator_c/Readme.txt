Name: farahi Zubair
Grup: 315CD

%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    1. eval_interpolator_c.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%
this function determines how quickly a polynomial converges Interpolation for ContinuousCase
This function returns N which is the number of nodes when the interpolation method is converges.
If the method does not converge, the function returns N = inf and displays in CommandWinow as"The interrogator does not converge.".
- Example:

>> eval_interpolator_c(2,1e-5)
The interrogator does not converge.

ans =

   Inf

>>

>> eval_interpolator_c(4,1e-5)

ans =

    66


%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    2. eval_interpolator_d.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%
this function determines how quickly a polynomial converges Interpolation for Discreet Case.
This function returns N which is the number of nodes when the interpolation method is converges.
If the method does not converge, the function returns N = inf and displays in CommandWinow as"The interrogator does not converge.".
- Example:

>> eval_interpolator_d(1,1e-5)

ans =

   300

>> eval_interpolator_d(5,1e-5)

ans =

   300

>> 