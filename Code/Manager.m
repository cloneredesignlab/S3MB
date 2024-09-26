
% global inputDataFiles;
global stiff
global glucose
global oxygen
global after1st
global before2nd
global days
global radStart
global radDose
global chemoStart
global chemoDose
global patientage
global resection_cav
global dat_ploidy
global path2params
global mgmtStatus
global reduceMRIweight

addpath utils
addpath utils/export_fig/
% addpath ~/Projects/code/MatlabCode/scripts/wassersteinFun
% % cd /Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/code/CaseWestern_GBM/IMOworkshop2022/S3MB
% % addpath ../
% % INDIR='../data/';

GROWELSEWHERE=false;
AGE='_ya';
% patient='002';
patients={ 'TCGA-06-0190', 'TCGA-06-0210','001', 'TCGA-06-0125', 'TCGA-06-0171', '002'};
numstartpts=100;

for patient= patients(:)'
    patient=char(patient);
    reduceMRIweight = false;
    if strcmp(patient,'TCGA-06-0190') || strcmp(patient,'TCGA-06-0210')
        reduceMRIweight = true;
    end
    disp(patient);
    %% Patient independent parameters
    poi={'alpha_death','v_max','Oth','delta_R'};
    path2params='AllParams_transition2real.xlsx';
    stiff= readmatrix(['data/',patient,'_stiffness',AGE,'.txt']);
    % stiff= readmatrix('data/',patient,'_stiffness_oa.txt');
    stiff(stiff(:)==0)=1;                %added this to make the outside values 1
    % do this n number of times to increase contrast on stiffness map
    % for i=1:4
    %     stiff(stiff~=1) = (stiff(stiff~=1) - min(stiff(stiff~=1)))./(max(stiff(stiff~=1))-min(stiff(stiff~=1)));
    % end
    glucose= readmatrix(['data/',patient,'_glucose',AGE,'.txt']);
    oxygen= readmatrix(['data/',patient,'_oxygen.txt']);

    %% Patient specific parameters:
    dat_ploidy=[1;1];
    clinical=readtable(['data/',patient,'_clinical.txt']);
    patientage=table2array(clinical(1,'age'));
    days=table2array(clinical(1,'days_2ndSurgery'));
    chemoStart=table2array(clinical(1,'chemoStart'));
    radStart=table2array(clinical(1,'radStart'));
    radDose=table2array(clinical(1,'radDose'));
    chemoDose=table2array(clinical(1,'chemoDose'));
    mgmtStatus = table2array(clinical(1,'MGMT_methylation_status'));
    % age=table2array(clinical(1,'age'));
    before1st= readmatrix(['data/',patient,'_t1t2_before1stT.txt']);
    after1st= readmatrix(['data/',patient,'_t1t2_after1stT.txt']);
    before2nd = readmatrix(['data/',patient,'_t1t2_before2ndT.txt']);
    resection_cav = readmatrix(['data/',patient,'_cavity_after1st.txt']);
    after1st(resection_cav>0)=0;

    if GROWELSEWHERE
        after1st= fliplr(flip(after1st));
        before1st = fliplr(flip(before1st));
        resection_cav = fliplr(flip(resection_cav));
    end

    %% Assume max intensity before 1st surgery corresponds to carrying capacity
    %% Same value is used to normalize after1stT and before2ndT
    allPars=readtable(path2params,'ReadRowNames',true,'ReadVariableNames',true);
    carryingCapacity = table2array(allPars('sigma',1));
    normFactor = carryingCapacity/max(before1st(:));
    before2nd = before2nd*normFactor;
    after1st = after1st*normFactor;
    before1st = before1st*normFactor;
    subplot(1,2,1); hist(after1st/carryingCapacity);
    subplot(1,2,2); hist(before2nd/carryingCapacity);

    %% MRI signal after 1st surgery is a combination of inflammation and tumor cell infiltration
    %% In an attempt to subtract inflammation, we use the signal before the 1st surgery outside of the resection cavity
    % after1st(resection_cav==0)=before1st(resection_cav==0);
    ii=find(after1st>before1st);
    after1st(ii)= before1st(ii);

    %% make sure no tumor cells outside brain:
    % @TODO -- should be fixed in data preprocessing script, not here
    after1st(stiff>=1)=0;


    %% adjust stiffness, oxygen and glucose at resection cavity:
    glucose(stiff==1)=0;
    oxygen(stiff==1)=0;
    % there could be an if statement to check if there is a resection cavity atlas
    % stiff(resection_cav>0) = 0.5*quantile(stiff(:),0.1); %decrease stiffness whereever there is resection cavity
    stiff = stiff - (stiff*0.1).*resection_cav; %decrease stiffness whereever there is resection cavity
    glucose = glucose - (glucose*0.5).*resection_cav; %decrease perfusion whereever there is resection cavity
    oxygen = oxygen - (oxygen*0.25).*resection_cav; %decrease perfusion whereever there is resection cavity
    % oxygen = oxygen + (max(oxygen(:))-oxygen).*resection_cav; %increase perfusion whereever there is resection cavity
    subplot(1,2,1)
    hist(glucose(resection_cav==0))
    subplot(1,2,2)
    hist(glucose(resection_cav>0))
    % Convert to mg/L (disolved Oxygen)
    disp("Stiffness")
    disp(quantile(stiff(stiff<1),0:0.1:1));
    disp("Glucose")
    disp(quantile(glucose(stiff<1),0:0.1:1));
    disp("Oxygen")
    disp(quantile(oxygen(stiff<1),0:0.1:1));
    a=tiledlayout(1,3);
    overlayAtlasTumor(glucose/3,glucose,a); title('Glucose')
    overlayAtlasTumor(oxygen/10,oxygen,a); title('oxygen')
    overlayAtlasTumor(stiff/3,stiff,a); title('stiff')

    %% Plot tumor before vs after
    a=tiledlayout(1,2);
    overlayAtlasTumor(before2nd/3,after1st,a); title('0 days after 1st surgery', 'FontSize',9)
    overlayAtlasTumor(before2nd/3,before2nd,a); title([num2str(days),' days after 1st surgery'], 'FontSize',9)



    %% Fit for one subpopulation only:
    opts = optimoptions(@fmincon);%,'Algorithm','active-set');
    ms = MultiStart('UseParallel',true);
    ms.Display='iter';
    if ~exist([patient,'_fmincon_params.mat'], 'file')
        A=[];
        b= [];
        Aeq=[];
        beq=[];
        pars=table2array(allPars(poi,{'Min','Max','Median'}));
        lb = pars(:,1)';
        ub = pars(:,2)';
        pars = pars(:,3)';
        %% Global optimization
        problem = createOptimProblem('fmincon','objective',...
            @cost,'x0',pars,'lb',lb,'ub',ub,'options',opts);
        rs = RandomStartPointSet('NumStartPoints',numstartpts);
        points = list(rs,problem);
        [pars,fval,exitflag,output,solutions] = run(ms,problem,CustomStartPointSet(points));
        % gs = GlobalSearch;
        % gs.Display='iter';
        % [pars,fval] = run(gs,problem);
        save([patient,'_fmincon_params.mat'], 'solutions');
    end


    %% Look at best fit
    load([patient,'_fmincon_params.mat'])
    dat_ploidy=[1;1];
    [~,ia]=min([solutions.Fval]);
    % ia=3;
    pars=solutions(ia).X;
    [TotalLiveCells,totalDeadCells,TotalGlucose,TotalOxygen, Vasculature, ploidies_all, allParameters]=S3MB(days, oxygen, glucose, ...
        radStart, radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, dat_ploidy(1,:), mgmtStatus{:}, path2params, ...
        'alpha_death',pars(1:size(dat_ploidy,2)), 'v_max',pars((size(dat_ploidy,2)+1):size(dat_ploidy,2)*2),'Oth',...
        pars((2*size(dat_ploidy,2)+1):size(dat_ploidy,2)*3),'delta_R',pars((3*size(dat_ploidy,2)+1):size(dat_ploidy,2)*4));
    %dat_ploidy=[0.6, 0.4; 0.8, 0.2];
    dat_ploidy = readPloidy(patient);
    % set parameters
    path2params=[patient,'_AllParams.xlsx'];
    allPars2=readtable(path2params,'ReadRowNames',true,'ReadVariableNames',true,'VariableNamingRule','preserve');
    allPars2(:,'Clone.1') = allPars(allPars2.Properties.RowNames,'Median');
    allPars2(:,'Clone.2') = allPars2(:,'Clone.1');
    if any(cellfun(@(x) strcmp(x,'Clone.3'), allPars2.Properties.VariableNames))
        allPars2(:,'Clone.3') = allPars2(:,'Clone.1');
    end
    % for a_mig, a=b_mig --> min, max correspond to first and second SP fits:
    allPars2('b_mig',{'Clone.1','Clone.2'})=allPars('b_mig',{'Min','Max'});
    % allPars2('v_max',1:size(dat_ploidy,2))= array2table(pars(2)*ones(1,size(dat_ploidy,2)));
    allPars2('v_max',1:size(dat_ploidy,2))= array2table(allParameters.v_max);
    allPars2('alpha_death',1:size(dat_ploidy,2))= array2table(allParameters.alpha_death);
    allPars2('Oth',1:size(dat_ploidy,2))= array2table(allParameters.Oth);
    allPars2('delta_R',1:size(dat_ploidy,2))= array2table(allParameters.delta_R);
    writetable(allPars2(union(poi,setdiff(allPars2.Properties.RowNames,poi),'stable'),:),['GBM_modelFit/',patient,'_AllParams.xlsx'],"WriteRowNames",true, 'WriteMode','replacefile');


    %% Fit for both subpopulations
    if size(dat_ploidy,2)>1
        if ~exist([patient,'_fmincon_params_sps.mat'], 'file')
            writetable(allPars2,path2params,"WriteRowNames",true, 'WriteMode','replacefile');
            % fit
            pars=table2cell(array2table(pars([1,3])));
            pars = cellfun(@(x) x*ones(1,size(dat_ploidy,2)), pars, 'UniformOutput', false);
            pars = cell2mat(pars);
            % also add delta
            pars=[pars,table2array(allPars2('delta_R','Clone.1'))];
            lb = 0.5*pars;
            ub = 1.5*pars;
            problem = createOptimProblem('fmincon','objective',...
                @cost,'x0',pars,'lb',lb,'ub',ub,'options',opts);
            rs = RandomStartPointSet('NumStartPoints',numstartpts);
            points = list(rs,problem);
            [pars,fval,exitflag,output,solutions] = run(ms,problem,CustomStartPointSet(points));
            save([patient,'_fmincon_params_sps.mat'], 'solutions');
        end


        %% Look at all good fits: CIs for inferred parameters:
        FVALMAX=0.6;
        o=[];
        for tmp = dir(fullfile(pwd,[patient, '_fmincon_params_sps*.mat']))'
            load(tmp.name)
            o_=arrayfun(@(x) x.X', solutions,'UniformOutput',0);
            o_=cell2mat(o_)';
            ii=find([solutions.Fval]<FVALMAX);
            o=[o;o_(ii,1:4)];
        end
        tbl=table(o');
        tbl.Properties.RowNames={'alpha_{death,1}';'alpha_{death,2}';'O_{th,1}';'O_{th,2}'};
        tbl.Name = categorical(tbl.Properties.RowNames)
        tbl.Clone= {'clone1';'clone2';'clone1';'clone2'};
        [~,p]=ttest(o(:,1),o(:,2),'Tail','left')
        disp(tbl.Name(1:2))
        [~,p]=ttest(o(:,3),o(:,4),'Tail','left')
        disp(tbl.Name(3:4))
        cellfun(@(x) quantile(x,[0.025,0.975]), num2cell(o,1), 'UniformOutput', false)
        disp(tbl.Name)
        close all hidden;
        figure('name',['GBM_modelFit/',patient,'_SPparametersInferred_CIs.png'],'Position', [584 645 229 270])
        boxplot(o, tbl.Properties.RowNames,'colorgroup',tbl.Clone, 'BoxStyle', 'filled'); ylabel('fitted value')
        set(gca,'TickLabelInterpreter', 'tex');
        set(gca,'xticklabel',tbl.Properties.RowNames)
        savefigs()


        %% Look at best fit
        load([patient,'_fmincon_params_sps.mat'])
        % dat_ploidy=[0.6, 0.4; 0.8, 0.2];
        dat_ploidy = readPloidy(patient);
        [~,ia]=min([solutions.Fval]);
        pars=solutions(ia).X;
        [TotalLiveCells,totalDeadCells,TotalGlucose, TotalOxygen, Vasculature, ploidies_all, allParameters]=S3MB(days, oxygen, glucose, radStart, ...
            radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, ...
            dat_ploidy(1,:), mgmtStatus{:}, path2params, 'alpha_death',pars(1:size(dat_ploidy,2)), 'Oth',pars((size(dat_ploidy,2)+1):size(dat_ploidy,2)*2), ...
            'delta_R', repmat(pars(length(pars)),1,size(dat_ploidy,2)));
        save(['GBM_modelFit/',patient,'_ploidies_all.mat'], 'ploidies_all');

        % writestruct(allParameters,'~/Downloads/AllParams.xml')
        allPars2('alpha_death',1:size(dat_ploidy,2))= array2table(allParameters.alpha_death);
        allPars2('Oth',1:size(dat_ploidy,2))= array2table(allParameters.Oth);
        allPars2('delta_R',1:size(dat_ploidy,2))= array2table(allParameters.delta_R);
        writetable(allPars2(union(poi,setdiff(allPars2.Properties.RowNames,poi),'stable'),:),['GBM_modelFit/',patient,'_AllParams.xlsx'],"WriteRowNames",true, 'WriteMode','replacefile');
    end

    %% Correlation between vasculature and tumor cell count across voxels?
    close all
    sf=15;
    ii=round(length(Vasculature)/10);
    v=imresize(Vasculature{ii}, sf/100);
    t=imresize(TotalLiveCells{ii}, sf/100);
    figure('name',['GBM_modelFit/',patient,'_VasculatureVsTCcount.png'])
    set(gcf,'Position',[100 100 700 250])
    subplot(1,2,2);
    imagesc(v); title('vasculature')
    subplot(1,2,1);
    ii=find(v>0.05);
    plot(t(ii),v(ii),'.','MarkerSize',30);
    ylabel(['Average vasculature across ',num2str(round(100^2/sf^2)),' voxels']);
    xlabel(['Average # tumor cells across ',num2str(round(100^2/sf^2)),' voxels'])
    if ~isempty(ii)
        corr(v(ii), t(ii))
    end
    disp(fitlm(v(ii), t(ii)))
    savefigs()


    % Ploidy distribution on day of second resection:
    ploidies = ploidies_all{days};
    % Plots:
    close all
    background=stiff*0.25*quantile(before2nd(:),1);
    background(background==max(background))=nan;
    figure('name',['GBM_modelFit/',patient,'_RecapEvolution_FromFirst2SecondSurgery.png'])
    set(gcf,'Position',[100 100 900 450])
    a=tiledlayout(2,3);
    a.TileSpacing = 'compact';
    a.Padding = 'compact';
    % image(TotalLiveCells{days}/2000)
    overlayAtlasTumor(background,after1st,a); title('0 days after 1st surgery', 'FontSize',9)
    overlayAtlasTumor(background,before2nd,a); title([num2str(days),' days after 1st surgery'], 'FontSize',9)
    overlayAtlasTumor(background,ploidies{1},a); title(['ploidy compartment 1 (on day ',num2str(days),' )'], 'FontSize',9)
    if length(ploidies)>1
        overlayAtlasTumor(background,ploidies{2},a);  title(['ploidy compartment 2 (on day ',num2str(days),')'], 'FontSize',9)
        disp(sum(ploidies{1}(:))/( sum(ploidies{2}(:))+sum(ploidies{1}(:))))
        tmp=ploidies{1}./(ploidies{1}+ploidies{2});
        tmp(ploidies{2}<500 & ploidies{1}<500)=nan;
        overlayAtlasTumor([],tmp,a);  title('ploidy compartment ratio 1/(1+2)', 'FontSize',9)
    end
    % subplot(3,2,5); image(background*100); title('Glucose atlas')
    overlayAtlasTumor(background,TotalLiveCells{days},a);  title(['total live cells (on day ',num2str(days),')'], 'FontSize',9)
    % overlayAtlasTumor(background,totalDeadCells{days},a);  title(['total dead cells (on day ',num2str(days),')'], 'FontSize',9)
    figure('name',['GBM_modelFit/',patient,'_GrowthDynamics.png'])
    set(gcf,'Position',[100   100   370   901])
    a=tiledlayout(6,1);
    a.TileSpacing = 'compact';
    a.Padding = 'compact';
    ax1 = nexttile();
    hold(ax1, 'on')
    plot(cellfun(@(x) nanmedian(x(resection_cav>0)), TotalGlucose),'LineWidth',3,'Color','black'); xlabel('day')
    yyaxis right
    plot(cellfun(@(x) nanmedian(x(resection_cav>0)), TotalOxygen),'LineWidth',3); xlabel('day')
    plot(cellfun(@(x) nanmean(x(resection_cav>0)), Vasculature),'.','LineWidth',3,'Color','blue');
    set(gca, 'FontSize',14)
    legend('Glucose (mMol)','Oxygen (mg/L)','Angiogenesis','Location','NorthEastOutside')
    hold off
    ax1 = nexttile();
    hold(ax1, 'on')
    plot(cellfun(@(x) nanmean(x(resection_cav>0)), TotalLiveCells),'LineWidth',3,'Color','cyan');xlabel('day')
    plot(cellfun(@(x) nanmean(x(:)), totalDeadCells),'LineWidth',3,'Color','magenta');
    set(gca, 'YScale', 'log','FontSize',14)
    legend('Alive (# cells)','Dead (# cells)','Location','NorthEastOutside')
    hold off
    ax1 = nexttile();
    hold(ax1, 'on')
    plot(cellfun(@(x) nanmean(x(resection_cav>0)), TotalLiveCells)); ylabel('Live cells (inside cavity)'); xlabel('day')
    ax1 = nexttile();
    hold(ax1, 'on')
    plot(cellfun(@(x) nanmean(x(resection_cav==0)), TotalLiveCells)); ylabel('Live cells (outside cavity)'); xlabel('day')
    ax1 = nexttile();
    hold(ax1, 'on')
    plot(cellfun(@(x) nanmean(x(:)), totalDeadCells)); ylabel('Dead cells (everywhere)'); xlabel('day')
    % fractionLiveCells_K=max(cellfun(@(x) mean(x(resection_cav>0)), TotalLiveCells))/table2array(allPars2('sigma',1));
    % disp(['Fraction carrying capacity: ',num2str(fractionLiveCells_K)])
    savefigs()
    close all

end

% ii= 1:1:days;
% background=stiff*0.25*quantile(before2nd(:),1);
% viewS3MBsim(length(ii), background,TotalLiveCells(ii),totalDeadCells(ii))
% background=0.00001*stiff*quantile(after1st(resection_cav>0),0.1);
% % viewS3MBsim(length(ii), background,TotalOxygen(ii),TotalOxygen(ii))
% viewS3MBsim(length(ii), background,Vasculature(ii),Vasculature(ii))
% % viewS3MBsim(length(ii), background,TotalGlucose(ii),TotalGlucose(ii))

