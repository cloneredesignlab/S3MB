function overlayAtlasTumor(stiff,TotalLiveCells, a)

%% correct orientation to clinical standart
stiff=flipud(stiff');
TotalLiveCells=flipud(TotalLiveCells');

if ~exist('a','var')
    a = tiledlayout(1,1);
end

ncol=300;
xfold=4;
q=0;

if ~isempty(stiff)
    TotalLiveCells = TotalLiveCells+ stiff;
    % q=ceil(1*ncol* max(stiff(:))/max(TotalLiveCells(:)));
    q=ceil(1*ncol* 1/xfold);
else
    stiff=zeros(size(TotalLiveCells,1),size(TotalLiveCells,2));
end
ax1 = nexttile();
hold(ax1, 'on')
imshow(stiff)

%     imshow('young_stiffness.jpg','Parent',ax1,'InitialMagnification', 1000);

imshow(TotalLiveCells,'Parent',ax1);


% Set transparency
alpha(.85)
colormap( ax1, [gray(q); jet(ncol-q)] );
colorbar( ax1 );
try
    caxis([0,xfold*max(max(stiff))] )
catch

end


end