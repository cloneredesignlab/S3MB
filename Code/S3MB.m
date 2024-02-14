%team Violet Aim 2 Code
%IMO workshop 10
%Authors: Frederika, Giada, Thomas and Parag from Team Violet
%11/1/2022
%last update 12/13/2022
function [TotalLiveCells, DeadCells, TotalGlucose, TotalOxygen, Vasc, N_all, allParameters, Diffusion]= S3MB(days, oxygen, glucose, radStart, radDose, chemoStart, chemoDose, imageStiff, resection_cav, after1st, n_ploidy, path2params, varargin)

p = inputParser;


%% Read in default parameters
allParameters = readtable(path2params,'ReadRowNames',true,'ReadVariableNames',false);

%% Overwrite default parameters with what we want to fit:
for var = allParameters.Properties.RowNames'
    default=allParameters(var,1:length(n_ploidy));
    addParameter(p,var{:}, table2array(default));
end
parse(p,varargin{:});

%% O2 dependent parameters should mirror Glucose dependent parameters because atlases are identical atm
allParameters= p.Results; %% return final parameters to user
% % @TODO: remove once we also have O2 atlas
% allParameters.alpha_gamma = allParameters.beta_gamma;
% allParameters.alpha_death = allParameters.beta_death;
% allParameters.iota_O = allParameters.iota_G;
% allParameters.omega_O = allParameters.omega_G;


%Radiation and chemotherapy conditions

doRad=(radStart~=0);                      %make this 1, if the simulation needs to account for effects of radiation and then change the two items below acccording to radiation dose and schedule
%     doRad=0;
Raddays = radStart:1:radStart+41;
Raddays = Raddays(mod(Raddays,7)~=5);
Raddays = Raddays(mod(Raddays,7)~=6);
gy = radDose/(6*5);
doChemo=(chemoStart~=0);                      %make this 1, if the simulation needs to account for effects of chemotherapy and then change the two items below acccording to chemo dose and schedule
%     doChemo=0;
chemoDays = chemoStart:1:chemoStart+42;
dose = chemoDose; % let's say this is 150 mg/m^2 of body


%grid
%this determines the simulation space resolution, make sure stiffness and perfusion initial maps match this resolution and grid size
n=size(glucose,2);
m=size(glucose,1);

%number of ploidies
%     n_ploidy=2;

%importing perfusion map data into indivudal Glucose, Phosphate and Oxygen initial maps
imageG=glucose;
imageP=imageG;
imageO=oxygen;

% @TODO: reading in age specific Glucose atlas

%not sure what this does, but is used for some plotting stuff ---
%Thomas Veith should know.
heatMat=zeros(m,n);
liveCellPlotter=zeros(1,days);
N1CellPlotter=zeros(1,days);
N2CellPlotter=zeros(1,days);

%define cell types and cell type dependent parameters
for pop=1:length(n_ploidy)
    %these are counts tracked only in space and stored anew each time steps, but tracked individually for each cell type
    N{pop}=zeros(m,n);        % tracks number of cells of type pop at each grid region
    N_old{pop}=zeros(m,n);    % stores number of cells of type pop at time t-1 in each grid region
    dyingcells{pop}=zeros(m,n);       % tracks number of cells dying at time t in each grid region
    dividingcells{pop}=zeros(m,n);    % tracks number of cells dividing at time t in each grid region
    convertcells{pop}=zeros(m,n);     % tracks number of cells converting from and to type pop at time t in each grid region
    cellsmigratingto{pop}=zeros(m,n); % tracks number of cells of type pop migrating into each grid region at time t

    %model parameters
    doublingrate(pop)=allParameters.gamma(pop);  %define population doubling rate /day for each cell type "pop" in the system
    a_death(pop)=allParameters.alpha_death(pop);    % oxygen dependent death constant per cell type
    b_death(pop)=allParameters.beta_death(pop);    % glucose dependent death constant per cell type
    a_gamma(pop)=allParameters.alpha_gamma(pop);    % oxygen dependent growth constant per cell type
    b_gamma(pop)=allParameters.beta_gamma(pop);     % glucose dependent growth constant per cell type
    a_mig(pop)=allParameters.a_mig(pop);    %biphasic migration rate constant 1 per cell type
    b_mig(pop)=allParameters.b_mig(pop);   %biphasic migration rate constant 2 per cell type
    c_r(pop)=allParameters.alpha_Tchemo(pop);    %oxygen dependence chemo driven cell death per cell type

    alpha(pop)=allParameters.alpha(pop);  % alpha in the linear quadratic equation of radiation sensitivity
    beta(pop)=allParameters.beta(pop);    % beta in the linear quadratic equation of radiation sensitivity
    rstar(pop)=allParameters.rstar(pop);    % linear quadratic equation of radiation sensitivity?
    delta_chemo(pop)=allParameters.delta_Tchemo(pop);    %Chemo shrinkage rate (how effective chemo is killing tumor cells)
    delta_r(pop)=allParameters.delta_R(pop);    %Starvation shrinkage rate (how effective lack of resources is killing tumor cells)
    v_max(pop)=allParameters.v_max(pop);     %top migration speed based on cell type (idea is higher ploidy makes it harder for cells to migrate, though that might not be necessarily true
end

%define cell type independent parameters

Stiff0 = allParameters.Stiff0(1); %max tissue stiffness
Oth = allParameters.Oth(1);     %threshold oxygen concentration to switch between oxygen and glucose based death and division
% C_f = 0;   %base cell propensity to increase ploidy
% C_b = 0;       %base cell propensity to decrease ploidy
Pth = allParameters.Pth(1);     %the minimum Phosphate required for the cells to multiply
% Gmax =25;     %some glucose dependent growth constant that limits higher ploidy cell growth
%define cell type independent parameters: resource consumption/generation
iota_O = allParameters.iota_O(1); % excess O2 consumption rate by tumor cells,
iota_G = allParameters.iota_G(1); % excess Glucose consumption rate by tumor cells,
iota_P = allParameters.iota_P(1); % excess PO4 consumption rate by tumor cells,
omega_O = allParameters.omega_O(1); % additional O2 generation rate by recruited vasculature
omega_G = allParameters.omega_G(1);% additional Glucose generation rate by recruited vasculature; @TODO: needs to be age specific
omega_P = allParameters.omega_P(1);% additional PO4 generation rate by recruited vasculature
sigma = allParameters.sigma(1); % Carrying capacity
p0 = allParameters.p0(1); % Max angiogenesis probability
c_d = allParameters.c_d(1); % clearance rate of dead cells via autophagy or other related mechanisms
c_0 = allParameters.c_0(1); % constant: radiation-induced change of stiffness
c_1 = allParameters.c_1(1); % constant: change in stiffness from remodelling by tumor cells  
c_2 = allParameters.c_2(1); % constant: quadratic parameter for optimal tissue stiffness for directed cell migration
c_3 = allParameters.c_3(1); % constant: linear parameter for optimal tissue stiffness for directed cell migration

% define initial stiffness

StiffIC= Stiff0*imageStiff;
Stiff = StiffIC;

avgStiff = mean(Stiff(2:m-1,2:n-1),"all");
stdStiff = std(Stiff(2:m-1,2:n-1),0,"all");

% define initial carrying capacity (this is the max cell mass including dead cells that can be in a grid space (purely an available space thing)

K = sigma * ones(m,n) + sigma*resection_cav;

%define initial G, O, P, DeadCells
G = ones(m,n).*imageG;
O = ones(m,n).*imageO;
P = ones(m,n).*imageP;

%initialize total cells etc per grid space...
for i=1:days
    % these are parameters that are tracked in time and space butccollectively for cell type except for N1 and N2
    TotalLiveCells{i} = zeros(m,n); %this is used in calculations
    Totalcells{i} = zeros(m,n);    %this is used in calculations
    DeadCells{i} = zeros(m,n);    %this is used in calculations

    TotalGlucose{i} = zeros(m,n);   %not used in calculations, only for reporting
    TotalOxygen{i} = zeros(m,n);   %not used in calculations, only for reporting
    Diffusion{i} = zeros(m,n);   %not used in calculations, only for reporting

    Vasc{i}=zeros(m,n);     %this is the new vasculature - beyond what brings in the base nutrients to the brain
    Radintensity{i}=zeros(m,n);      %this tracks radiation intensity when applied
    chemoDeath{i}=zeros(m,n);       %this tracks chemotherapy when applied
end


% %Seed initial tumor cells at random locations
%         xyCoords = [randi([2,m],1,1), randi([2,n],1,1)];
%         while (O(xyCoords(1),xyCoords(2))>(mean(O(2:m-1,2:n-1), 'all')*3.25)) || (O(xyCoords(1),xyCoords(2))<(mean(O(2:m-1,2:n-1), 'all')*0.35))
%             xyCoords = [randi([2,m],1,1),randi([2,n],1,1)];
%         end

% use this to force tumors to start at a fixed location - left
% frontal lobe?
%         xyCoords=[14,14];

for i=1:length(n_ploidy)
    N_=N{i};
    %             N_(xyCoords(1),xyCoords(2)) = (0.5+0.5*rand)*100;    %randomly seeds anywhere from 5k to 10k tumor cells in grid space defined by xyCoords
    N_=after1st.*(n_ploidy(i));
    %             disp(max(max(N_)))

    N{i} = N_;
end

% THIS FORCES N1>=N2 INITIAL CONDITIONS
% N{1}(xyCoords(1),xyCoords(2)) = (0.5+0.5*rand)*10000;
% N{2}(xyCoords(1),xyCoords(2))= (10000)-N{1}(xyCoords(1),xyCoords(2));

% grab initial frequency of N1 for stats at end
%     initFreqs=sum(N{1})/sum(N{1} + N{2});

%this is the main code that updates cell populations in every grid space at each time
N_all={};
for t=1:days
    N_all{t} = N;
    N_old=N;

    %calculate total cells in a voxel
    for ct=1:length(n_ploidy)
        TotalLiveCells{t} = TotalLiveCells{t}+N_old{ct};
    end

    if t==1
        Totalcells{t} = TotalLiveCells{t};
        DeadCells{t}=0;
    else
        DeadCells{t} = DeadCells{t-1}*c_d;      %this clears deadcells at the rate of 0.5/time interval, 0.5 is arbitrary and can be changed
        Totalcells{t} = TotalLiveCells{t}+DeadCells{t-1};
        Vasc{t} = Vasc{t-1};
    end

    TotalGlucose{t} = G;        %use this to report Glucose levels in the brain at each time point
    TotalOxygen{t} = O; 
    % Calculate new stiffness and resources based on total cells in a voxel
    %     O = O-(iota_O*O.*TotalLiveCells{t})+omega_O*Vasc{t};
    %     G = G-(iota_G*G.*TotalLiveCells{t})+omega_G*Vasc{t};   % the key assumption is that the existing and new vasculature sufficiently distributes nutrients within a grid space, but diffusion distance of these nutrients into neighboring grid spaces falls off rapidly
    % %     G = G-(16*iota_G*G.*TotalLiveCells{t}).*(O<=Oth)-(iota_G*G.*TotalLiveCells{t}).*(O>Oth)+omega_G*Vasc{t};
    %     P = P-(iota_P*P.*TotalLiveCells{t})+omega_P*Vasc{t};
    confluence = TotalLiveCells{t}./K;
    O = O-(iota_O*O.*confluence)+omega_O*Vasc{t};
    G = G-(16*iota_G*G.*confluence).*(O<=Oth)-(iota_G*G.*confluence).*(O>Oth)+omega_G*Vasc{t};
    P = P-(iota_P*P.*confluence)+omega_P*Vasc{t};

    %     Stiff = StiffIC.*(1+(Totalcells{t}./K).*(Totalcells{t}>0.05*K));   %assumption: crowding of tumor cells increases local matrix density and consequently stiffness, max brain tissue stiffness can only increase by 100%
    Stiff = StiffIC.*(1+c_1*(Totalcells{t}./K).*(Totalcells{t}>0.05*K)); %assumption: tumor cells increases local matrix density and consequently stiffness, max tumor ECM density can reach 22.5 kPa as per David Odde's Nat. Matls paper
    %Stiff = StiffIC+22.5*(Totalcells{t}./K).*(Totalcells{t}>0.05*K)); %assumption: tumor cells increases local matrix density and consequently stiffness, max tumor ECM density can reach 22.5 kPa as per David Odde's Nat. Matls paper
    %kill cells

    if doRad       %if radiation is given
        if ismember(t,Raddays)      % at specific points in time
            tumorRadiation = TotalLiveCells{t};
            tumorRadiation = tumorRadiation(:,:,1) >= 0;   %set this to give radiation in only regions above specific cell count, keep it zero if radiation is given everywhere
            Vasc{t}(tumorRadiation(:,:,1) >= 1) = 0;     %kill vasculature where radiation is given
            % Stiff(tumorRadiation(:,:,1) >= 1) = Stiff(tumorRadiation(:,:,1) >= 1)*0.2;     %lower local stiffness (need to check how this is affected as cell density drops too post radiation, might be some conflict)
            StiffIC(tumorRadiation(:,:,1) >= 1) = StiffIC(tumorRadiation(:,:,1) >= 1)*(c_0^gy);
            Radintensity{t} = gy*tumorRadiation;     %radiation intensity in each space based on dosage
        end
    end

    if doChemo      %if radiation is given
        if ismember(t,chemoDays)      % at specific points in time
            chemoTherapy = TotalLiveCells{t};
            chemoTherapy = chemoTherapy(:,:,1) >= 0;
            chemoDeath{t} = dose*chemoTherapy;     %radiation intensity in each space based on dosage
        end
    end

    for ct = 1:length(n_ploidy)
        d{ct}=1-(O./(a_death(ct)+O)).*(O>Oth)-(G./(b_death(ct)+G)).*(O<Oth);      %resource limited cell death rate
        radiationDeathRate{ct} = ((O>Oth).*(1-exp(-alpha(ct)*Radintensity{t}-beta(ct)*Radintensity{t}.^2))+...
            (Oth>=O).*(1-exp(-(alpha(ct)/rstar(ct))*Radintensity{t}-(beta(ct)/rstar(ct))*Radintensity{t}.^2))) * doRad;         %radiation based cell death rate (multiplying by doRad might be unnecessary
        chemoDeathRate{ct} = (1-(O./(c_r(ct)+O))) * doChemo;
        dyingcells{ct}=delta_r(ct)*d{ct}.*N_old{ct}+radiationDeathRate{ct}.*N_old{ct}+delta_chemo(ct)*chemoDeathRate{ct}.*N_old{ct};            % overall cell death (has additional constants 0.1 and 0.5 which may need to be integrated with existing constants)
        DeadCells{t} = DeadCells{t}+dyingcells{ct};         %update total dead cells in each grid space
    end

    %multiply cells
    for ct = 1:length(n_ploidy)
        gamma{ct}=((O./(a_gamma(ct)+O)).*(O>Oth)+(G./(b_gamma(ct)+G)).*(O<=Oth)).*(P>Pth);   %resource limited growth rate
        %         if ct==1
        dividingcells{ct}=doublingrate(ct)*(gamma{ct}.*N_old{ct}).*(1-Totalcells{t}./K);    %logistic cell growth

%         disp(max(max(Totalcells{t})))
%         disp(min(min(O)))
        %         else
        %             dividingcells{ct}=doublingrate(ct)*(0.5*tanh(G-Gmax/2)+3/2).*(gamma{ct}.*N_old{ct}).*(1-Totalcells{t}./K);     %logistic cell growth for higher ploidies with a glucose effect
        %         end
    end

    for ct=1:length(n_ploidy)
        N{ct}=N_old{ct}-dyingcells{ct}+dividingcells{ct}+convertcells{ct};     %updated cell population because of death, division and conversion... I think the order of events is important.
%         disp(max(max(dividingcells{ct})))
%         disp(max(max(dyingcells{ct})))
%         disp(['maximum total cells before migration:', num2str(max(max(N{ct})))]);
        if max(N{ct}(:))<1
            N{ct}(:)=0;
        end
    end
    % migrate cells ... this also includes code for vasculature recruitment because this is the only time we run through each individual grid space one by one

    for i=2:m-1        % i, j define current grid space
        for j=2:n-1
            neightotallive(i,j)=0.0001;     %used in vasculature growth probability, need to initialize as close to zero as possible
            neightotaldead(i,j)=0.0001;      %used in vasculature growth probability, need to initialize as close to zero as possible
            for ct=1:length(n_ploidy)
                if N{ct}(i,j)>0
                    a=i+randperm(3,3)-2;     % a,b define neigbhors of i, j in a random order (they also include i and j, hence the if statement 3 lines later)
                    b=j+randperm(3,3)-2;
                    for x=1:length(a)
                        for y=1:length(b)
                            if (a(x)~=i) || (b(y)~=j)
                                %                                 k= v_max(ct) * exp(-a_mig(ct)*G(i,j)) * (1-exp(-b_mig(ct)*G(i,j))) * ...
                                %                                     (Stiff(a(x), b(y))<Stiff(1,1))*(Stiff(a(x), b(y))>500) * (-0.25e-6*(Stiff(a(x),b(y)))^2+1.25e-3*(Stiff(a(x),b(y)))-0.5625);   %migration rate from i,j to a, b

                                k= v_max(ct) * exp(-a_mig(ct)*G(i,j)) * (1-exp(-b_mig(ct)*G(i,j))) * ...
                                    (Stiff(a(x), b(y))<Stiff(1,1))*(Stiff(a(x), b(y))>500) * (c_2*(Stiff(a(x),b(y)))^2+c_3*(Stiff(a(x),b(y))));   %migration rate from i,j to a, b. follows David Odde's 2022 Nat Matl paper.

                                cellsmigratingto{ct}(a(x),b(y))=cellsmigratingto{ct}(a(x),b(y))+k*N{ct}(i,j);     % number of cell migrated to a, b... will get updated for all i, j (all neighbors of a, b)
                                N{ct}(i,j)=N{ct}(i,j)-k*N{ct}(i,j);      % remove cells migrating away from i, j from total cells in i and j... repeat for all a, b
                                neightotallive(i,j) = neightotallive(i,j)+TotalLiveCells{t}(a(x),b(y));       %tracking number of live and cells in neighboring grid spaces to i, j
                                neightotaldead(i,j) = neightotaldead(i,j)+DeadCells{t}(a(x),b(y));            %tracking number of live and cells in neighboring grid spaces to i, j

                                Diffusion{t}(i,j) = Diffusion{t}(i,j)+k;
                            end
                        end
                    end
                    %recruit vasculature
                    if Vasc{t}(i,j) == 0 
                        probVasc = p0*  (Totalcells{t}(i,j)>(0.005*K(i,j)) )...
                            *(neightotallive(i,j)*neightotaldead(i,j)) / (neightotallive(i,j)+neightotaldead(i,j))^2; %probability of recruiting vasculature based on number of live and dead cells in the neighboring grid spaces
                        if rand<probVasc
                            Vasc{t}(i,j) = 1;      % recruit vasculature if the above stochastic condition is satisfied
                        end
                    end
                end
            end
        end
    end
    for ct=1:length(n_ploidy)
        N{ct}=N{ct}+cellsmigratingto{ct};      % add all cells that have migrated into a grid space (tracked by a,b coordinates in the above code) to the total number of cells in that grid space
        cellsmigratingto{ct}=zeros(m,n);       % empty the cellsmigratingto counter for next cycle
%         disp(['maximum total cells after migration:', num2str(max(max(N{ct})))]);
%         disp('.')
    end
end

% N{ct} has all the relevant cell population numbers for the time point t

% some specific counters for plotting and analysis (@Tommy might know these better)
heatMat=(heatMat+TotalLiveCells{days});
liveCellPlotter=(liveCellPlotter+cellfun(@(x) sum(x, 'all'), TotalLiveCells));
% N1CellPlotter=(N1CellPlotter+cellfun(@(x) sum(x, 'all'), N1ploidyCells));
% N2CellPlotter=(N2CellPlotter+cellfun(@(x) sum(x, 'all'), N2ploidycells));


% more analysis and plotting stuff that Tommy and Frederika will know better
%     ageCellPloter{find(age==age)}={liveCellPlotter, N1CellPlotter, N2CellPlotter};

%     T = table(initFreqs,totalCellCounts,totalN1Counts,totalN2Counts,totalDeadCells);
%     colors = ['r'; 'g'; 'b';'c'];
% %     figure()
% %     boxplot(table2array(T(:,2:5)), 1:4)
% %     h = findobj(gca,'Tag','Box');
% %     for j=1:length(h)
% %         patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
% %     end
% %     set(gca,'xticklabel',{'totalCellCount' 'totalN1' 'totalN2' 'totalDead'})
%     writetable(T,strcat('runHolder', outputStatsFilename,'.txt'))
%     title(age)

% figure()
% heatmap(heatMat)

%     figure()
%     plot(1:t,liveCellPlotter, 'LineWidth', 3, 'Color', 'b')
%     hold on
%     plot(1:t,N1CellPlotter, 'LineWidth', 3, 'Color', 'r')
%     hold on
%     plot(1:t,N2CellPlotter, 'LineWidth', 3, 'Color', 'g')
%     hold off
%     xlabel('Days')
%     ylabel('Live Cells')
%     legend('Total', 'N1', 'N2')
%     title(age)


% tHolder{find(ages==age)}=T;
% totalCellsHolder{find(ages==age)}=Totalcells;


% figure()
% heatmap(TotalLiveCells{days})
%
% figure()
% plot(1:t,cellfun(@(x) sum(x, 'all'), TotalLiveCells), 'LineWidth', 3, 'Color', 'b')
% hold on
% plot(1:t,cellfun(@(x) sum(x, 'all'), N1ploidyCells), 'LineWidth', 3, 'Color', 'r')
% hold on
% plot(1:t,cellfun(@(x) sum(x, 'all'), N2ploidycells), 'LineWidth', 3, 'Color', 'g')
% hold off
% xlabel('Days')
% ylabel('Live Cells')
% legend('Total', 'N1', 'N2')
end
