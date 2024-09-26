function  dat_ploidy=readPloidy(patient, varargin)
dat_ploidy_p= readmatrix(['data/',patient,'_ploidyComp_primary.txt']);
dat_ploidy_r= readmatrix(['data/',patient,'_ploidyComp_recurrent.txt']);
dat_ploidy = [dat_ploidy_p(:,2)';dat_ploidy_r(:,2)'];
dat_ploidy = dat_ploidy./repmat(sum(dat_ploidy,2),1,size(dat_ploidy,2));
if ~isempty(varargin)
    includePloidy = varargin{1};
    if includePloidy
        dat_ploidy = [dat_ploidy_p(:,1)'; dat_ploidy];
    end
end
end