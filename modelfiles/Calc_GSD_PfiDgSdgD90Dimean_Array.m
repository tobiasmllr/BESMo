function [Pfi,Dg,SDg,D90,Dimean] = Calc_GSD_PfiDgSdgD90Dimean_Array(Fi_cum,Sizes,gsd_MG,nodes_N)
% [Pfi,Dg,SDg,D90,Dimean]=GSD(Fi_cum,Sizes)
% GSD requires 2 inputs to run:  a cumulative distribution function (Fi_cum)
% and the corresponding sizes for the cumulative distribution in mm.
% It provides 4 outputs: grain size frequencies (Pfi), Fs( percent of sand)
% geometric mean (Dg), geometric standard deviation (SDg) and grain size at
% which 90 % is finer than(D90). 
%All size estimations have its output in m.
% It transforms the sizes from mm to psi units to make easier all
% statistical calculations. Finally the results are transformed to m.

format long g

%keyboard
m=gsd_MG+1;
% Size scale transformation (from mm to psi units)
%psi=NaN(length(Sizes),1);
%for d=1:m
%psi=log2(Sizes(:));
psi=log2(Sizes);
%end


% Obtain surface frequencies (0-1 scale)
Pfi=Fi_cum(:,2:end)-Fi_cum(:,1:end-1);
Pfi(Pfi<0)=0;    % replace all values smaller than 0 with 0

% Dg estimation
% 1. Estimate geometric mean of each size class (Dgi)which is equivalent to
% the arithmetic mean in psi units
%for d=1:length(psi)-1
Dim_psi=(psi(:,1:end-1)+psi(:,2:end))/2;
%end

% Calculate the geometric mean of each size class in mm:
Dimean = (2.^Dim_psi);

% 2. Estimate product between Dim and Pfi
Dim_Pfi=Dim_psi.*Pfi;

% 3. Estimate geometric mean in psi units (Dg_psi)
Dg_psi=sum(Dim_Pfi,2);
% 4. Transform Dg in psi units to mm.
Dg=(2.^Dg_psi);
%
% Estimate geometric standard deviation in psi (SDg_psi)
% For this it will be used the equation provided in ch.3, p.11, eq 3.5
% (Parker,2007)
% 1. Estimate ((Dim_psi-Dg_psi)^2)*Pfi
Aux=NaN(nodes_N, m-1);
for d=1:m-1
    Aux(:,d)=((Dim_psi(:,d)-Dg_psi).^2).*Pfi(:,d);
end
% 2. Estimate SD in psi by sum Aux(:,j)and transform. This is a
% dimensionless parameter so it doesn�t required to be transformed from mm
% to m.
SDg_psi=(sum(Aux,2)).^(1/2);
SDg=(2.^SDg_psi); % grain size standard deviation in mm
% Estimate D90
D90=NaN(nodes_N,1);
for n=1:nodes_N
    for d=1:m-1
        if(0.9 >Fi_cum(n,d))&&(0.9<=Fi_cum(n,d+1)) % Identify limits
            slope=(Fi_cum(n,d+1)-Fi_cum(n,d))./(psi(n,d+1)-psi(n,d)); % line slope
            D90_psi=psi(n,d)+(0.9-Fi_cum(n,d))./slope; % linear equation
            D90(n)=2.^D90_psi;
        end
    end
end
end
