% quantify similarity between the output of S3MB and the patient data.
function J  = cost(pars)
global after1st;
global before2nd;
global stiff
global glucose
global oxygen
global days
global radStart
global radDose
global chemoStart
global chemoDose
global resection_cav
global dat_ploidy
global path2params

%% Fixed parameters
tmp=readtable(path2params,'ReadRowNames',true,'ReadVariableNames',false);
carryingCapacity = table2array(tmp('sigma',1));
MINCELLS2CAREABOUT=0.05*carryingCapacity;


%% Parameters to be optimized
%     a_mig = pars(1:size(dat_ploidy,2));
%     b_mig= pars((size(dat_ploidy,2)+1):(size(dat_ploidy,2)*2));
%% 2D matrix of cell densities + 2D matrix of ploidies
%   [TotalLiveCells,totalDeadCells, ~, ~, ~, ploidies]=S3MB(days, oxygen, glucose, radStart, radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, dat_ploidy(1,:), path2params, 'a_mig',a_mig, 'b_mig', b_mig);

%% Parameters to be optimized
if size(dat_ploidy,2)==1
    alpha_death = pars(1:size(dat_ploidy,2));
    v_max = pars((size(dat_ploidy,2)+1):(size(dat_ploidy,2)*2));
    Oth = pars((size(dat_ploidy,2)*2+1):(size(dat_ploidy,2)*3));
    delta_R = pars((size(dat_ploidy,2)*3+1):(size(dat_ploidy,2)*4));
    %% 2D matrix of cell densities + 2D matrix of ploidies
    [TotalLiveCells,totalDeadCells, ~, ~, ~, ploidies]=S3MB(days, oxygen, glucose, radStart, radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, dat_ploidy(1,:), path2params, 'alpha_death',alpha_death, 'v_max', v_max,'Oth',Oth,'delta_R',delta_R);

else
    alpha_death = pars(1:size(dat_ploidy,2));
    Oth = pars((size(dat_ploidy,2)+1):(size(dat_ploidy,2)*2));
    delta_R = repmat(pars(length(pars)),1,size(dat_ploidy,2));
    % disp(alpha_death)
    % disp(Oth)
    % disp(delta_R)
    %% 2D matrix of cell densities + 2D matrix of ploidies
    [TotalLiveCells,totalDeadCells, ~, ~, ~, ploidies]=S3MB(days, oxygen, glucose, radStart, radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, dat_ploidy(1,:), path2params, 'alpha_death',alpha_death, 'Oth', Oth,'delta_R',delta_R);
end
%% Ploidy distribution on day of second resection:
ploidies= ploidies{days};

%% sum of live + dead tumor cell densities
live_dead = TotalLiveCells{days}+totalDeadCells{days};

%% Ploidy composition averaged over space --> column 1: ploidy; column 2: total cell count
ploidy = cellfun(@(x) sum(sum(x)),ploidies);
ploidy = ploidy /sum(ploidy);

%RMSE tumor density
diff=before2nd./(live_dead+MINCELLS2CAREABOUT);
diff(diff>1) = 1./diff(diff>1);
diff= (1- diff).^2;
ii = find(live_dead + before2nd >=MINCELLS2CAREABOUT ); %% true negatives (regions without tumor correctly predicted as such) don't count, otherwise small tumors will have better fit
diff = diff(ii);

%RMSE ploidy
%     diff2= (ploidy-dat_ploidy(2,:)).^2;
diff2= (1-(ploidy./dat_ploidy(2,:))).^2;

J = sqrt(mean(diff(:)) + mean(diff2));
%     disp(J)

%     ii = find(live_dead + before2nd >=MINCELLS2CAREABOUT ); %% true negatives (regions without tumor correctly predicted as such) don't count, otherwise small tumors will have better fit
%     wsd = ws_distance(before2nd(ii), live_dead(ii));
%     wsd2 = ws_distance(max(before2nd(:))*ploidy, max(before2nd(:))*dat_ploidy(2,:));
%     J = mean([wsd,wsd2]);
end