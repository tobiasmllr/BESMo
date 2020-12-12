function [ D50, D90, GSD_Sizes, GSD_Pfi ] = fkt_generateGSD(GS_Classes,GS_D50,GS_Sigma,GSDistributionType)
%FKT_GENERATEGSD Summary of this function goes here
%   Detailed explanation goes here

GS_D50_psi = log(GS_D50)/log(2);
GS_inPsi = [0:0.01:10];
GS_inMM = 2.^GS_inPsi;


% create cumulative probability density function:
switch GSDistributionType
    case 'Normal'
        mu = GS_D50_psi;
        GS_curveProbability = cdf('Normal',GS_inPsi,mu,GS_Sigma);
    case 'LogNormal'
        mu = log(GS_D50_psi);
        GS_curveProbability = cdf('LogNormal',GS_inPsi,mu,GS_Sigma);
    otherwise
        error(strcat('This Grain Size Distribution type is not supported: ',GSDistributionType))
end

D5 = GS_inMM(find(GS_curveProbability > 0.05,1,'first'));
D10 = GS_inMM(find(GS_curveProbability > 0.10,1,'first'));
D16 = GS_inMM(find(GS_curveProbability > 0.16,1,'first'));
D25 = GS_inMM(find(GS_curveProbability > 0.25,1,'first'));
D34 = GS_inMM(find(GS_curveProbability > 0.34,1,'first'));
D50 = GS_D50;
D66 = GS_inMM(find(GS_curveProbability > 0.66,1,'first'));
D75 = GS_inMM(find(GS_curveProbability > 0.75,1,'first'));
D84 = GS_inMM(find(GS_curveProbability > 0.84,1,'first'));
D90 = GS_inMM(find(GS_curveProbability > 0.90,1,'first'));
D95 = GS_inMM(find(GS_curveProbability > 0.95,1,'first'));

% GS_Classes =  1; %: D50
% GS_Classes =  3; %: D16 D50 D84
% GS_Classes =  5; %: D10 D16 D50 D84 D90
% GS_Classes =  7; %: D10 D16 D34 D50 D66 D84 D90
% GS_Classes =  9; %: D5 D10 D16 D34 D50 D66 D84 D90 D95
% GS_Classes = 11; %: D5 D10 D16 D25 D34 GS_D50 D66 D75 D84 D90 D95

switch GS_Classes
    case 1
        min = D50 - 0.01;
        max = D50 + 0.01;
        GSD_Sizes = [min max];
        GSD_Pfi = [0 100];
        D90 = D50 + 0.01;
    case 3
        min = D16 - 0.01;
        max = D84 + 0.01;
        GSD_Sizes = [min D16 D50 D84 max];
        GSD_Pfi = [0 16 50 84 100];
    case 5
        min = D10 - 0.01;
        max = D90 + 0.01;
        GSD_Sizes = [min D10 D16 D50 D84 D90 max];
        GSD_Pfi = [0 10 16 50 84 90 100];
    case 7
        min = D10 - 0.01;
        max = D90 + 0.01;
        GSD_Sizes = [min D10 D16 D34 D50 D66 D84 D90 max];
        GSD_Pfi = [0 10 16 34 50 66 84 90 100];
    case 9
        min = D5 - 0.01;
        max = D95 + 0.01;
        GSD_Sizes = [min D5 D10 D16 D34 D50 D66 D84 D90 D95 max];
        GSD_Pfi = [0 5 10 16 34 50 66 84 90 95 100];
    case 11
        min = D5 - 0.01;
        max = D95 + 0.01;
        GSD_Sizes = [min D5 D10 D16 D25 D34 D50 D66 D75 D84 D90 D95 max];
        GSD_Pfi = [0 5 10 16 25 34 50 66 75 84 90 95 100];
    otherwise
        error(strcmp('wrong number of grain size classes given:', int2str(GS_Classes)))
end
         
end

