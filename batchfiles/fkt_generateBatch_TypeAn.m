function [ RunConfig ] = fkt_generateBatch_TypeAn(RunConfig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long e

%% Take values from normal Type A:
% repurpuse value of P_distanceIncrease for number of runs to keep
checkNRuns = RunConfig.P_distanceIncrease;

% use P_distance Start same as P_distanceIncrease
RunConfig.P_distanceIncrease = RunConfig.P_distanceStart;

% Get values from normal Type A:
RunConfig = fkt_generateBatch_TypeA(RunConfig);
PDistA_before = RunConfig.P_distanceArrayInMin;

%% Calculate ideal PDist for new number of runs
square = @(x) x.^2;

idealPDist =  (PDistA_before(end) - PDistA_before(1))/(checkNRuns - 1);
idealPDistArrayLinear = round(PDistA_before(1):idealPDist:idealPDist*checkNRuns);
idealPDistArraySquare = square(idealPDistArrayLinear./idealPDistArrayLinear(end)) .* idealPDistArrayLinear(end);

Array_minDiff = zeros(1,checkNRuns);
Array_index = zeros(1,checkNRuns);
Array_value = zeros(1,checkNRuns);

% find indices and values of closest PDist
for i=2:checkNRuns-1
    [Array_minDiff(i), Array_index(i)] = min(abs(idealPDistArraySquare(i) - PDistA_before(:)));
    Array_value(i) = PDistA_before(Array_index(i));
end

% set fixed starting and ending PDist
Array_index(1) = 1;
Array_index(end) = length(PDistA_before);



% Identify which values to add:
%PossibleAddition = ~ismember([1:1:length(PDistA_before)],Array_indexUnique);
%PossibleDistance = PDistA_before(PossibleAddition);
% Add pulse distances of experimets, if not already present:
Ind10h = find(PDistA_before==10*60);
Ind20h = find(PDistA_before==20*60);
Ind40h = find(PDistA_before==40*60);

% remove non-unique values:
Array_indexUnique = unique([Array_index Ind10h Ind20h Ind40h]);
Array_indexUnique = sort(Array_indexUnique);


% IF ABOVE NUMBER: delete a middle run
if length(Array_indexUnique) > checkNRuns
    NumberTooMany = length(Array_indexUnique) - checkNRuns;
    GotAlreaDistance = PDistA_before(Array_indexUnique);
    removepossible = find(GotAlreaDistance > (40*60),NumberTooMany,'first');
    Array_indexUnique(removepossible) = [];
    if NumberTooMany > 1
        warning('Had to remove more than one pulse period to match simulation run number. Consider increasing number or decreasing max pulse distance!')
    end
end

% IF BELOW NUMBER: fill duplicates with least-represented PDist (longest distance in indexes)
if length(Array_indexUnique) < checkNRuns
    found_all_values = false;
    counter = 1;
    while ~found_all_values
        
        
        %warning('did not find enough pulse configuration to match requirement of RunConfig.P_distanceIncrease! \n Trying to interpolate!')
        while length(Array_indexUnique) < checkNRuns
            gradientUnique=gradient(Array_indexUnique);
            [maxDistance, ind]= max(gradientUnique);
            IndexToAdd = round(ind + maxDistance / 2);
            Array_indexUnique = sort([Array_indexUnique IndexToAdd]);
        end
        Array_indexUnique = unique(Array_indexUnique);
        if length(Array_indexUnique) == checkNRuns
            found_all_values = true;
            break;
        else
            counter = counter + 1;
        end
        if counter > 20
            % Identify which values to add:
            NumberStillNeeded= checkNRuns-length(Array_indexUnique);
            PossibleAddition = ~ismember([1:1:length(PDistA_before)],Array_indexUnique);
            GotAlreaDistance = PDistA_before(Array_indexUnique);
            PossibleDistance = PDistA_before(PossibleAddition);
            IndexToAdd       = find(PossibleAddition,NumberStillNeeded,'last');
            Array_indexUnique = sort([Array_indexUnique IndexToAdd]);
            %addDistances = PossibleDistance(end-2:end);
            
            if length(Array_indexUnique) < checkNRuns
                error('did not find pulse distances. Stuck in loop... Aborting...\n INCREASE MAX PULSE DISTANCE')
            end
        end
    end
end

% store shorter PDist Array:
RunConfig.P_distanceArrayInMin = PDistA_before(Array_indexUnique);
RunConfig.P_Number = RunConfig.P_Number(Array_indexUnique);
RunConfig.P_magnitudeArr = RunConfig.P_magnitudeArr(Array_indexUnique);
end

