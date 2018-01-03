clear all
clc

% test the interpolators

% continous case
N=zeros(2,6);
eps=1e-5;
for tip=1:6;
    N(1,tip)=eval_interpolator_c(tip, eps);
    N(2,tip)=eval_interpolator_d(tip, eps);
end
N
    