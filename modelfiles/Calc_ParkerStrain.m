function [ intervect,omega0inter,sigma0inter ] = Calc_ParkerStrain( phisg0_ST, w0_ST, s0_ST )
%CALC_PARKERSTRAIN 
% Interpolate values of the Strain functions (using cubic Splines)
format long e

ninter1 = 50;
ninter2 = 1000;
ninter3 = 50;
inf_bound = log(min(phisg0_ST));
bound_1 = log(0.8601);
bound_2 = log(38.57);
sup_bound = log(max(phisg0_ST));
dinter_1 = (bound_1 - inf_bound)/ninter1;
dinter_2 = (bound_2 - bound_1)/ninter2;
dinter_3 = (sup_bound - bound_2)/ninter3;

for i = 1:ninter1       
    intervectlog = inf_bound + dinter_1*(i-1);    
    intervect1(i,1) = exp(intervectlog);
end

for i = 1:ninter2       
    intervectlog = bound_1 + dinter_2*(i-1);    
    intervect2(i,1) = exp(intervectlog);
end

for i = 1:ninter3+1       
    intervectlog = bound_2 + dinter_3*(i-1);    
    intervect3(i,1) = exp(intervectlog);
end

intervect = [phisg0_ST;intervect1;intervect2;intervect3];
intervect = sort(intervect);

diffInterv = diff(intervect);
zeroDiff = find(diffInterv==0);
minDiff = min(diffInterv(diffInterv>0));
intervect(zeroDiff+1) = intervect(zeroDiff) + minDiff*0.9;
intervect = sort(intervect);

% Interpolacio omega0
omega0inter = spline(phisg0_ST,w0_ST,intervect);

% Interpolacio sigma0
sigma0inter = spline(phisg0_ST,s0_ST,intervect);

end

