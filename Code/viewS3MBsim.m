
function [] = viewS3MBsim(days, stiff,TotalLiveCells,totalDeadCells)

themax=max(cellfun(@(x) quantile(x(:),0.99), TotalLiveCells))+max(cellfun(@(x) quantile(x(:),0.99), totalDeadCells));

a = tiledlayout(1,1);

ax1 = nexttile();
imshow(stiff)
hold(ax1, 'on')
exportgraphics(gcf,'~/Downloads/testAnimated.gif','Append',false);
for t=1:days
    hold(ax1, 'off')
    imshow(stiff)
    hold(ax1, 'on')
    title(a,['Time = ',sprintf('%d',t),' days'])
    I = TotalLiveCells{t}+totalDeadCells{t};
    %     themax=max(I(:));

    %     imshow('young_stiffness.jpg','Parent',ax1,'InitialMagnification', 1000);
    imshow(I,'Parent',ax1);

    % Set transparency
    alpha(.5)
    colormap( ax1, jet );
    colorbar( ax1 );
    caxis([0,3*max(max(stiff))] )
    %         caxis([0,themax])

    pause(0.1)
    exportgraphics(gcf,'~/Downloads/testAnimated.gif','Append',true);



end