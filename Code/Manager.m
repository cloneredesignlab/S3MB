
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

addpath utils
addpath utils/export_fig/
% addpath ~/Projects/code/MatlabCode/scripts/wassersteinFun
% % cd /Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/code/CaseWestern_GBM/IMOworkshop2022/S3MB
% % addpath ../
% % INDIR='../data/';

GROWELSEWHERE=false;
AGE='_ya';

%% Patient independent parameters
path2params='AllParams_transition2real.xlsx';
stiff= readmatrix(['data/stiffness',AGE,'.txt']);
% stiff= readmatrix('data/stiffness_oa.txt');
stiff(stiff(:)==0)=1;                %added this to make the outside values 1
% do this n number of times to increase contrast on stiffness map
% for i=1:4
%     stiff(stiff~=1) = (stiff(stiff~=1) - min(stiff(stiff~=1)))./(max(stiff(stiff~=1))-min(stiff(stiff~=1)));
% end
glucose= readmatrix(['data/glucose',AGE,'.txt']);
oxygen= readmatrix('data/oxygen.txt');

%% Patient specific parameters:
dat_ploidy=[1;1];
clinical=readtable('data/clinical.txt');
patientage=table2array(clinical(1,'age'));
days=table2array(clinical(1,'days_2ndSurgery'));
chemoStart=table2array(clinical(1,'chemoStart'));
radStart=table2array(clinical(1,'radStart'));
radDose=table2array(clinical(1,'radDose'));
chemoDose=table2array(clinical(1,'chemoDose'));
% age=table2array(clinical(1,'age'));
before1st= readmatrix('data/t1t2_before1stT.txt');
after1st= readmatrix('data/t1t2_after1stT.txt');
before2nd = readmatrix('data/t1t2_before2ndT.txt');
resection_cav = readmatrix('data/cavity_after1st.txt');
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
overlayAtlasTumor(before2nd/3,before2nd,a); title('415 days after 1st surgery', 'FontSize',9)



%% Fit for one subpopulation only:
A=[];
b= [];
Aeq=[];
beq=[];
pars=table2array(allPars({'alpha_death','v_max','Oth','delta_R'},{'Min','Max','Median'}));
lb = pars(:,1)';
ub = pars(:,2)';
pars = pars(:,3)';
% pars=table2array(allPars('alpha_death',1:size(dat_ploidy,2)));
% pars=[pars,table2array(allPars('v_max',1:size(dat_ploidy,2)))];
% pars=[pars,table2array(allPars('Oth',1:size(dat_ploidy,2)))];
% pars=[pars,table2array(allPars('delta_R',1:size(dat_ploidy,2)))];
% % pars=[pars,table2array(allPars('iota_G',1:size(dat_ploidy,2)))];
% lb = pars-0.9*pars;
% lb(1:size(dat_ploidy,2)*2)=0.5*lb(1:size(dat_ploidy,2)*2);
% ub = pars*2;
% ub((size(dat_ploidy,2)+1):size(dat_ploidy,2)*2)=2*ub((size(dat_ploidy,2)+1):size(dat_ploidy,2)*2);
%% Global optimization
opts = optimoptions(@fmincon);%,'Algorithm','active-set');
problem = createOptimProblem('fmincon','objective',...
    @cost,'x0',pars,'lb',lb,'ub',ub,'options',opts);
rs = RandomStartPointSet('NumStartPoints',10);
points = list(rs,problem);
ms = MultiStart('UseParallel',true);
ms.Display='iter';
[pars,fval,exitflag,output,solutions] = run(ms,problem,CustomStartPointSet(points));
% gs = GlobalSearch;
% gs.Display='iter';
% [pars,fval] = run(gs,problem);
save('fmincon3_params.mat', 'solutions');


%% Look at all good fits: CIs for inferred parameters:
o=[]; 
o_b=[];
for tmp = dir(fullfile(pwd,'fmincon*_params.mat'))'
    load(tmp.name)
    o_=arrayfun(@(x) x.X', solutions,'UniformOutput',0);
    o_b_=arrayfun(@(x) cell2mat(x.X0)', solutions,'UniformOutput',0)
    o_=cell2mat(o_)';
    o_b_=cell2mat(o_b_)';
    ii=find([solutions.Fval]<=0.55);
    o=[o;o_(ii,:)];
    o_b=[o_b;o_b_];
end
tbl=[table(o','VariableNames',{'posterior'}),table(o_b','VariableNames',{'prior'})];
tbl.Properties.RowNames={'alpha_{death}';'v_{max}';'O_{th}';'delta_{R}'};
close all hidden;
figure('name','GBM_modelFit/parametersInferred_CIs.png'); i = 1
for what = tbl.Properties.RowNames'
    subplot(2,2,i); i = i + 1;
    hold on
    h3 = histogram(tbl.prior(what,:),"FaceColor",'blue','Normalization','probability');
    h3.NumBins = 5;
    h3.BinWidth = range(tbl.prior(what,:))/5;
    h3.FaceAlpha =0.5;
    h2 = histogram(tbl.posterior(what,:),"FaceColor",'red','Normalization','probability');
    h2.NumBins = h3.NumBins;
    h2.BinWidth = h3.BinWidth;
    h2.FaceAlpha =0.5;
    legend({['prior (n=',num2str(size(o_b,1)),')'],['posterior (n=',num2str(size(o,1)),')']},'Location','NorthEastOutside');
    xlabel(what)
    ylabel('probability')
    title(['[',num2str(round(quantile(tbl.posterior(what,:),[0.025,0.975]),2)),']'])
end
savefigs()


%% Look at best fit
load('fmincon_params.mat')
dat_ploidy=[1;1];
[~,ia]=min([solutions.Fval]);
% ia=3;
pars=solutions(ia).X;
[TotalLiveCells,totalDeadCells,TotalGlucose,TotalOxygen, Vasculature, ploidies_all, allParameters]=S3MB(days, oxygen, glucose, radStart, radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, dat_ploidy(1,:), path2params, ...
    'alpha_death',pars(1:size(dat_ploidy,2)), 'v_max',pars((size(dat_ploidy,2)+1):size(dat_ploidy,2)*2),'Oth',pars((2*size(dat_ploidy,2)+1):size(dat_ploidy,2)*3),'delta_R',pars((3*size(dat_ploidy,2)+1):size(dat_ploidy,2)*4));
% writestruct(allParameters,'~/Downloads/AllParams.xml')


%% Fit for both subpopulations
dat_ploidy=[0.6, 0.4; 0.8, 0.2]; %readmatrix(inputDataFiles(3));
% set parameters
path2params="AllParams.xlsx";
allPars2=readtable(path2params,'ReadRowNames',true,'ReadVariableNames',true);
allPars2(:,'Value') = allPars(allPars2.Properties.RowNames,'Median');
allPars2(:,'Value_ifDifferentFor2ndClone_') = allPars2(:,'Value');
allPars2({'a_mig','b_mig'},{'Value','Value_ifDifferentFor2ndClone_'})=allPars({'a_mig','b_mig'},{'Min','Max'});  %% for a_mig, a=b_mig --> min, max correspond to first and second SP fits
allPars2('delta_R',1:size(dat_ploidy,2))= array2table(pars(4)*ones(1,size(dat_ploidy,2)));
allPars2('Oth',1:size(dat_ploidy,2))= array2table(pars(3)*ones(1,size(dat_ploidy,2)));
allPars2('v_max',1:size(dat_ploidy,2))= array2table(pars(2)*ones(1,size(dat_ploidy,2)));
allPars2('alpha_death',1:size(dat_ploidy,2))= array2table(pars(1)*ones(1,size(dat_ploidy,2)));
writetable(allPars2,path2params,"WriteRowNames",true);
% fit
pars=table2cell(array2table(pars([1,3])));
pars = cellfun(@(x) x*ones(1,size(dat_ploidy,2)), pars, 'UniformOutput', false);
pars = cell2mat(pars);
% also add delta
pars=[pars,table2array(allPars2('delta_R','Value'))];
lb = 0.5*pars;
ub = 1.5*pars;
problem = createOptimProblem('fmincon','objective',...
    @cost,'x0',pars,'lb',lb,'ub',ub,'options',opts);
rs = RandomStartPointSet('NumStartPoints',10);
points = list(rs,problem);
[pars,fval,exitflag,output,solutions] = run(ms,problem,CustomStartPointSet(points));
save('fmincon_params_sps.mat', 'solutions');


%% Look at all good fits: CIs for inferred parameters:
FVALMAX=0.6;
o=[];
for tmp = dir(fullfile(pwd,'fmincon_params_sps*.mat'))'
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
figure('name','GBM_modelFit/SPparametersInferred_CIs.png','Position', [584 645 229 270])
boxplot(o, tbl.Properties.RowNames,'colorgroup',tbl.Clone, 'BoxStyle', 'filled'); ylabel('fitted value')
set(gca,'TickLabelInterpreter', 'tex');
 set(gca,'xticklabel',tbl.Properties.RowNames)
savefigs()


%% Look at best fit
load('fmincon_params_sps.mat')
dat_ploidy=[0.6, 0.4; 0.8, 0.2];
[~,ia]=min([solutions.Fval]);
ia=2;
pars=solutions(ia).X;
[TotalLiveCells,totalDeadCells,TotalGlucose, TotalOxygen, Vasculature, ploidies_all, allParameters]=S3MB(days, oxygen, glucose, radStart, radDose, chemoStart, chemoDose, stiff, resection_cav, after1st, ...
    dat_ploidy(1,:), path2params, 'alpha_death',pars(1:size(dat_ploidy,2)), 'Oth',pars((size(dat_ploidy,2)+1):size(dat_ploidy,2)*2),'delta_R', repmat(pars(length(pars)),1,size(dat_ploidy,2)));
% writestruct(allParameters,'~/Downloads/AllParams.xml')
allPars2('alpha_death',1:size(dat_ploidy,2))= array2table(allParameters.alpha_death);
allPars2('Oth',1:size(dat_ploidy,2))= array2table(allParameters.Oth);
allPars2('delta_R',1:size(dat_ploidy,2))= array2table(allParameters.delta_R);
poi={'Oth','alpha_death','v_max','delta_R'};
writetable(allPars2(union(poi,setdiff(allPars2.Properties.RowNames,poi),'stable'),:),"GBM_modelFit/AllParams.xlsx","WriteRowNames",true);


%% Correlation between vasculature and tumor cell count across voxels?
close all
sf=15;
ii=round(length(Vasculature)/10);
v=imresize(Vasculature{ii}, sf/100);
t=imresize(TotalLiveCells{ii}, sf/100);
figure('name','GBM_modelFit/VasculatureVsTCcount.png')
set(gcf,'Position',[100 100 700 250])
subplot(1,2,2);
imagesc(v); title('vasculature')
subplot(1,2,1);
ii=find(v>0.05);
plot(t(ii),v(ii),'.','MarkerSize',30);
ylabel(['Average vasculature across ',num2str(round(100^2/sf^2)),' voxels']);
xlabel(['Average # tumor cells across ',num2str(round(100^2/sf^2)),' voxels'])
corr(v(ii), t(ii))
disp(fitlm(v(ii), t(ii)))
savefigs()

% %% write out del n_i for Stefano
% for t =1:length(TotalLiveCells)
%     dNdy = diff(TotalLiveCells{t},1,2);
%     dNdx = diff(TotalLiveCells{t},1,1);
%     writematrix(dNdx,['~/Downloads/dNdx_T',num2str(t),'.txt'])
%     writematrix(dNdy,['~/Downloads/dNdy_T',num2str(t),'.txt'])
% end


% Ploidy distribution on day of second resection:
ploidies = ploidies_all{days};
% Plots:
close all
background=stiff*0.25*quantile(before2nd(:),1);
background(background==max(background))=nan;
figure('name','GBM_modelFit/RecapEvolution_FromFirst2SecondSurgery.png')
set(gcf,'Position',[100 100 900 450])
a=tiledlayout(2,3);
a.TileSpacing = 'compact';
a.Padding = 'compact';
% image(TotalLiveCells{days}/2000)
overlayAtlasTumor(background,after1st,a); title('0 days after 1st surgery', 'FontSize',9)
overlayAtlasTumor(background,before2nd,a); title('415 days after 1st surgery', 'FontSize',9)
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
figure('name','GBM_modelFit/GrowthDynamics.png')
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
fractionLiveCells_K=max(cellfun(@(x) mean(x(resection_cav>0)), TotalLiveCells))/table2array(allPars2('sigma',1));
disp(['Fraction carrying capacity: ',num2str(fractionLiveCells_K)])



ii= 1:1:days;
background=stiff*0.25*quantile(before2nd(:),1);
viewS3MBsim(length(ii), background,TotalLiveCells(ii),totalDeadCells(ii))
background=0.00001*stiff*quantile(after1st(resection_cav>0),0.1);
% viewS3MBsim(length(ii), background,TotalOxygen(ii),TotalOxygen(ii))
viewS3MBsim(length(ii), background,Vasculature(ii),Vasculature(ii))
% viewS3MBsim(length(ii), background,TotalGlucose(ii),TotalGlucose(ii))
