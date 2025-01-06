function [stuff] = plot_all_EBC(root, out, include_cross_corr, square_or_polar, include_spikemap, chasing_maps, firing_rate)

    addpath '/Users/pearls/Work/RSC_project/MATLAB_scripts/EgocentricBoundaryCells-master/EgocentricBoundaryCells-master/slanCM'
    
    if chasing_maps
        ncols = 11; % this one has padding
        ncols = 5; % this one does not
        nrows = 1;

        if square_or_polar == 1 % square EBC plot
            % the +pi/2 brings "forwards" to "up"
            [x,y] = meshgrid(size(out.cfull.rm,2):-1:1,1:size(out.cfull.rm,1));
        else % polar EBC plot
            % the +pi/2 brings "forwards" to "up"
            binz = out.cfull.params.thetaBins;
            binz(end+1) = -binz(1);
            [t2, r2] = meshgrid(wrapTo2Pi(binz+pi/2), out.cfull.params.distanceBins(1:end-1));
            [x, y] = pol2cart(t2,r2);
        end

        figHandle = figure(1); 
        figHandle.Position = [70 200 1300 400];
        ax = gca;

        %% chasing full
        subplot(nrows,ncols,5)
        rm = out.cfull.rm; 
        %rm(isnan(rm)) = 0;
        occ = out.cfull.occ;
        rm(occ < 50) = nan; %threshold data at 50bins (8ms) for less noisy maps
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
    
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(rm(:), 99)])
        set(gca,'YDir','Normal')  

        if square_or_polar == 1 
            rectangle('Position',[min(min(x)) min(min(y)) max(max(x)) max(max(y))])
        else
            rectangle('Position',[min(min(y)) min(min(y)) max(max(y))*2 max(max(y))*2],'Curvature',[1 1])
        end
        axis off; axis square
        title('Chase')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
        %colorbar
        %% padding for same size
        %[x, y] = pol2cart(t2,r2);
        for iii=6:11
            break
            subplot(nrows,ncols,iii)
            rm = out.cfull.rm; 
            rm(isnan(rm)) = 0;
            if square_or_polar == 0
                rm(:,end+1) = rm(:,1);
            end
            surface(x,y, rm); shading interp
        
            hold on
            set(gca,'XTick',[],'YTick',[])
            
            %colormap(slanCM('viridis'))
            colormap(parula)
            set(gca, 'YDir','Normal','CLim',[0 prctile(out.cfull.rm(:), 99)])
            set(gca,'YDir','Normal')  
            axis off; axis square
            title('Chase')
            if include_cross_corr
                set(gca,"FontSize",7)
            else
                set(gca,"FontSize",7)
            end
        end

        %% chasing odd
        subplot(nrows,ncols,2)
        rm = out.codd.rm;
        %rm(isnan(rm)) = 0;
        occ = out.codd.occ;
        rm(occ < 50) = nan; %threshold data at 50bins (8ms) for less noisy maps
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
        
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(rm(:), 99)])
        set(gca,'YDir','Normal')

        if square_or_polar == 1 
            rectangle('Position',[min(min(x)) min(min(y)) max(max(x)) max(max(y))])
        else
            rectangle('Position',[min(min(y)) min(min(y)) max(max(y))*2 max(max(y))*2],'Curvature',[1 1])
        end
        axis off; axis square
        title('Chase odd')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
        %colorbar
        %% chasing even
        
        subplot(nrows,ncols,3)
     
        rm = out.ceven.rm;
        %rm(isnan(rm)) = 0;
        occ = out.ceven.occ;
        rm(occ < 50) = nan; %threshold data at 50bins (8ms) for less noisy maps
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
        
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(rm(:), 99)])
        set(gca,'YDir','Normal')  

        if square_or_polar == 1 
            rectangle('Position',[min(min(x)) min(min(y)) max(max(x)) max(max(y))])
        else
            rectangle('Position',[min(min(y)) min(min(y)) max(max(y))*2 max(max(y))*2],'Curvature',[1 1])
        end
        axis off; axis square
        title('Chase even')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
        %colorbar
       
        %% even odd chasing CC
        rmA = out.codd.rm_ns(:,1:end);
        rmB = out.ceven.rm_ns(:,1:end);
        occA = out.codd.occ;
        occB = out.ceven.occ;
        rmA(occA < 50) = nan;
        rmB(occB < 50) = nan;
        %rmA = out.codd.rm_ns(:,1:end-1);
        %rmB = out.ceven.rm_ns(:,1:end-1);
    
        offs_dist = (-size(rmA,1)+2):(size(rmA,1)-2);
        if mod(size(rmA,2),2) == 1
            offs_angle = (-floor(size(rmA,2)/2)):(floor(size(rmA,2)/2));
        else
            offs_angle = (-size(rmA,2)/2):(size(rmA,2)/2-1);
        end
        cc_plot = NaN(length(offs_dist), length(offs_angle));
        
        for i=1:length(offs_angle)
            rotB = circshift(rmB,offs_angle(i),2);
            for j=1:length(offs_dist)
                if offs_dist(j) < 0
                    Bpart = rotB((-offs_dist(j)+1):end,:);
                    Apart = rmA(1:(size(rmA,1)+offs_dist(j)),:);
                end
                if offs_dist(j) > 0
                    Bpart = rotB(1:size(rmB,1)-offs_dist(j),:);
                    Apart = rmA(offs_dist(j)+1:end,:);
                end
                if offs_dist(j) == 0
                    Apart = rmA;
                    Bpart = rotB;
                end
                Apart = Apart(:);
                Bpart = Bpart(:);
                not_nan = not(isnan(Apart)) .* not(isnan(Bpart));
                Apart = Apart(not_nan == 1);
                Bpart = Bpart(not_nan == 1);
                
                if length(Apart) > 5 & length(Bpart) > 5
                    corr = corrcoef(Apart, Bpart);
                    cc_plot(j,i) = corr(1,2);
                end
            end
        end
        
        thresh_cc_plot = cc_plot > prctile(cc_plot,75, "all");
        thresh_cc_plot = imfill(thresh_cc_plot, 'holes');
        [labeled_cc_plot, ~] = bwlabel(thresh_cc_plot, 8);
        props = regionprops(labeled_cc_plot, cc_plot, 'all');
    
        bigEnough = vertcat(props.Area) > 9; % large enough to be a blob, currently 3x3 square bin
        if sum(bigEnough) == 0 % all blobs are too small, pick the largest one
            bigEnough = vertcat(props.Area) == max(vertcat(props.Area));
        end
        allBlobCentroids = vertcat(props.Centroid);
        allBlobCentroids(:,1) = allBlobCentroids(:,1) + offs_angle(1)-1;
        allBlobCentroids(:,2) = allBlobCentroids(:,2) + offs_dist(1)-1;
        allBlobCentroids = allBlobCentroids(bigEnough,:);
        centDists = zeros(size(allBlobCentroids,1),1);
        for k=1:size(allBlobCentroids,1)
            centDists(k) = sqrt((allBlobCentroids(k,1))^2 + (allBlobCentroids(k,2))^2);
        end
        [~,whichBlob] = min(centDists);
        coordsWhichBlob = allBlobCentroids(whichBlob,:);
    
        subplot(nrows,ncols,4);
        [x,y] = meshgrid(offs_angle, offs_dist);
        surface(x, y, (cc_plot)); shading interp
        
        hold on
        rectangle('Position',[min(min(x)) min(min(y)) max(max(x))*2+1 max(max(y))*2+1])
        [~, ind] = max(cc_plot,[],'all');
        [row_max,col_max] = ind2sub(size(cc_plot),ind);
        %disp([row_max, col_max])
        %disp([length(offs_angle), length(offs_dist)])
        %disp([offs_angle(col_max), offs_dist(row_max)])

        scatter3(coordsWhichBlob(1), coordsWhichBlob(2),1,20,'black','filled');
    
        set(gca,'XTick',[],'YTick',[])
        ylim([-50,50])
        %colormap(slanCM('viridis'))
        colormap(parula)
        min_clim = mink((min(cc_plot)),2);
        max_clim = maxk((max(cc_plot)),2);
        set(gca, 'YDir','Normal','CLim',[min_clim(2),max_clim(2)])
        set(gca,'YDir','Normal')
        set(gca,'YDir','reverse')
        axis off
        %axis square
        posi = get(gca, 'Position');
        set(gca,'DataAspectRatio',[length(offs_angle) length(offs_dist) 1], 'Position',[posi(1),posi(2)-0.4,posi(3),posi(4)+0.8])
        %set(gca,'DataAspectRatio',[size(out.OFfull.rm_ns,2) size(out.OFfull.rm_ns,1) 1])
        %title('CC, even/odd', 'Position', [0 13.2])
        %colorbar
        set(gca,"FontSize",7)

        %% if include spikemaps
        ax_ch = subplot(nrows,ncols,include_spikemap);
        hold on;
        plot(root.cfull.x,root.cfull.y,'Color',[.7 .7 .7])
        colormap(ax_ch,hsv)
        xlim([min(root.cfull.x) max(root.cfull.x)]); ylim([min(root.cfull.y) max(root.cfull.y)])
        cx = root.cfull.x(root.cfull.spike == 1);
        cy = root.cfull.y(root.cfull.spike == 1);
        ch = root.cfull.md(root.cfull.spike == 1);
        scatter(cx,cy,7,ch,'filled')
        
        scatter(out.cfull.QP(:,1),out.cfull.QP(:,2),30,'k','filled')
        set(gca,'YDir','Normal')
        %title('Trajectory')
        %colorbar('Ticks', [-pi/2, 0, pi/2], 'TickLabels',{0,90,180})
        axis off
        axis square

        stuff = cc_plot;
    else
        if include_cross_corr
            ncols = 9;
        else
            ncols = 6;
        end
        nrows = 1;
    
        if include_spikemap
            ncols = ncols + 2;
        end
    
        if square_or_polar == 1 % square EBC plot
            % the +pi/2 brings "forwards" to "up"
            [x,y] = meshgrid(1:size(out.OFfull.rm,2),1:size(out.OFfull.rm,1));
        else % polar EBC plot
            % the +pi/2 brings "forwards" to "up"
            binz= out.cfull.params.thetaBins;
    
            binz(end+1) = -binz(1);
            [t2, r2] = meshgrid(wrapTo2Pi(binz+pi/2), out.OFfull.params.distanceBins(1:end-1));
            [x, y] = pol2cart(t2,r2);
        end
        figHandle = figure(1); 
        figHandle.Position = [70 200 1300 400];
        ax = gca;
        %% ratemap circular, OF full
        subplot(nrows,ncols,7+include_spikemap*2);
        rm = out.OFfull.rm;
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
    
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(out.OFfull.rm(:), 99)])
        set(gca,'YDir','Normal')
        axis off; axis square;
        title('OF')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
    
        %% ratemap circular, OF odd
        subplot(nrows,ncols,1+include_spikemap);
        % the +pi/2 brings "forwards" to "up"
        %[t2, r2] = meshgrid(wrapTo2Pi(out.OFodd.params.thetaBins+pi/2), out.OFodd.params.distanceBins(1:end-1));
        %[x, y] = pol2cart(t2,r2);
        rm = out.OFodd.rm;
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
    
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(out.OFodd.rm(:), 99)])
        set(gca,'YDir','Normal')  
        axis off; axis square
        title('OF odd')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
        
    
        %% ratemap circular, OF even
        
        if include_cross_corr
            subplot(nrows,ncols,2+include_spikemap)
        else
            subplot(nrows,ncols,2+include_spikemap);
        end
        % the +pi/2 brings "forwards" to "up"
        %[t2, r2] = meshgrid(wrapTo2Pi(out.OFeven.params.thetaBins+pi/2), out.OFeven.params.distanceBins(1:end-1));
        %[x, y] = pol2cart(t2,r2);
        rm = out.OFeven.rm;
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
    
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(out.OFeven.rm(:), 99)])
        set(gca,'YDir','Normal')  
        axis off; axis square
        title('OF even')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
    
        %% ratemap circular, c full
        
        if include_cross_corr
            subplot(nrows,ncols,8+include_spikemap*2)
        else
            subplot(nrows,ncols,8+include_spikemap*2);
        end
        % the +pi/2 brings "forwards" to "up"
        %[t2, r2] = meshgrid(wrapTo2Pi(out.cfull.params.thetaBins+pi/2), out.cfull.params.distanceBins(1:end-1));
        %[x, y] = pol2cart(t2,r2);
        rm = out.cfull.rm;
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
    
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(out.cfull.rm(:), 99)])
        set(gca,'YDir','Normal')  
        axis off; axis square
        title('Chase')
        % Add firing rate annotation
        % 
        if nargin > 6 && ~isempty(firing_rate)
    annotation_text = sprintf('Firing Rate: %.2f Hz', firing_rate);
    annotation('textbox', [0.15, 0.85, 0.3, 0.05], 'String', annotation_text, ...
               'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold', ...
               'BackgroundColor', 'white', 'FitBoxToText', 'on');
        end
        if include_cross_corr
    set(gca,"FontSize",7)
        else
    set(gca,"FontSize",7)
        end
    
        %% ratemap circular, c odd
        
        if include_cross_corr
            subplot(nrows,ncols,4+include_spikemap*2)
        else
            subplot(nrows,ncols,4+include_spikemap*2);
        end
        % the +pi/2 brings "forwards" to "up"
        %[t2, r2] = meshgrid(wrapTo2Pi(out.codd.params.thetaBins+pi/2), out.codd.params.distanceBins(1:end-1));
        %[x, y] = pol2cart(t2,r2);
        rm = out.codd.rm;
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
    
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(out.codd.rm(:), 99)])
        set(gca,'YDir','Normal')  
        axis off; axis square
        title('Chase odd')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
    
        %% ratemap circular, c even
        
        if include_cross_corr
            subplot(nrows,ncols,5+include_spikemap*2)
        else
            subplot(nrows,ncols,5+include_spikemap*2);
        end
        % the +pi/2 brings "forwards" to "up"
        %[t2, r2] = meshgrid(wrapTo2Pi(out.ceven.params.thetaBins+pi/2), out.ceven.params.distanceBins(1:end-1));
        %[x, y] = pol2cart(t2,r2);
        rm = out.ceven.rm;
        if square_or_polar == 0
            rm(:,end+1) = rm(:,1);
        end
        surface(x,y, rm); shading interp
        
        hold on
        set(gca,'XTick',[],'YTick',[])
        
        %colormap(slanCM('viridis'))
        colormap(parula)
        set(gca, 'YDir','Normal','CLim',[0 prctile(out.ceven.rm(:), 99)])
        set(gca,'YDir','Normal')  
        axis off; axis square
        title('Chase even')
        if include_cross_corr
            set(gca,"FontSize",7)
        else
            set(gca,"FontSize",7)
        end
    
        %% Make the cross corr maps, first of odd/even
        rmA = out.OFodd.rm_ns(:,1:end);
        rmB = out.OFeven.rm_ns(:,1:end);
    
        offs_dist = (-size(rmA,1)+2):(size(rmA,1)-2);
        if mod(size(rmA,2),2) == 1
            offs_angle = (-floor(size(rmA,2)/2)):(floor(size(rmA,2)/2));
        else
            offs_angle = (-size(rmA,2)/2):(size(rmA,2)/2-1);
        end
        cc_plot = NaN(length(offs_dist), length(offs_angle));
        
        for i=1:length(offs_angle)
            rotB = circshift(rmB,offs_angle(i),2);
            for j=1:length(offs_dist)
                if offs_dist(j) < 0
                    Bpart = rotB((-offs_dist(j)+1):end,:);
                    Apart = rmA(1:(size(rmA,1)+offs_dist(j)),:);
                end
                if offs_dist(j) > 0
                    Bpart = rotB(1:size(rmB,1)-offs_dist(j),:);
                    Apart = rmA(offs_dist(j)+1:end,:);
                end
                if offs_dist(j) == 0
                    Apart = rmA;
                    Bpart = rotB;
                end
    
                Apart = Apart(:);
                Bpart = Bpart(:);
                not_nan = not(isnan(Apart)) .* not(isnan(Bpart));
                Apart = Apart(not_nan == 1);
                Bpart = Bpart(not_nan == 1);
                
                if length(Apart > 5) & length(Bpart > 5)
                    corr = corrcoef(Apart, Bpart);
                    cc_plot(j,i) = corr(1,2);
                end
            end
        end
    
        thresh_cc_plot = cc_plot > prctile(cc_plot,75, "all");
        thresh_cc_plot = imfill(thresh_cc_plot, 'holes');
        [labeled_cc_plot, ~] = bwlabel(thresh_cc_plot, 8);
        props = regionprops(labeled_cc_plot, cc_plot, 'all');
    
        bigEnough = vertcat(props.Area) > 9; % large enough to be a blob, currently 3x3 square bin
        allBlobCentroids = vertcat(props.Centroid);
        allBlobCentroids(:,1) = allBlobCentroids(:,1) + offs_angle(1)-1;
        allBlobCentroids(:,2) = allBlobCentroids(:,2) + offs_dist(1)-1;
        allBlobCentroids = allBlobCentroids(bigEnough,:);
        centDists = zeros(size(allBlobCentroids,1),1);
        for k=1:size(allBlobCentroids,1)
            centDists(k) = sqrt((allBlobCentroids(k,1))^2 + (allBlobCentroids(k,2))^2);
        end
        [~,whichBlob] = min(centDists);
        coordsWhichBlob = allBlobCentroids(whichBlob,:);
    
        assignin('base','test',cc_plot)
        subplot(nrows,ncols,3+include_spikemap);
        [x,y] = meshgrid(offs_angle, offs_dist);
        surface(x, y, (cc_plot)); shading interp
        
        hold on
    
        [~, ind] = max(cc_plot,[],'all');
        [row_max,col_max] = ind2sub(size(cc_plot),ind);
        scatter3(coordsWhichBlob(1), coordsWhichBlob(2),1,20,'black','filled');
    
        set(gca,'XTick',[],'YTick',[])
        ylim([-50,50])
        %colormap(slanCM('viridis'))
        colormap(parula)
        min_clim = mink((min(cc_plot)),2);
        max_clim = maxk((max(cc_plot)),2);
        set(gca, 'YDir','Normal','CLim',[min_clim(2),max_clim(2)])
        set(gca,'YDir','Normal') 
        set(gca,'YDir','reverse')
        axis off;
        % axis square;
        posi = get(gca, 'Position');
        set(gca,'DataAspectRatio',[length(offs_angle) length(offs_dist) 1], 'Position',[posi(1),posi(2)-0.4,posi(3),posi(4)+0.8])
        %set(gca,'DataAspectRatio',[size(out.OFfull.rm_ns,2) size(out.OFfull.rm_ns,1) 1])
        %title('CC, even/odd', 'Position', [0 13.2])
        %colorbar
        set(gca,"FontSize",7)
    
        %% then chasing odd/even
        rmA = out.codd.rm_ns(:,1:end);
        rmB = out.ceven.rm_ns(:,1:end);
        %rmA = out.codd.rm_ns(:,1:end-1);
        %rmB = out.ceven.rm_ns(:,1:end-1);
    
        offs_dist = (-size(rmA,1)+2):(size(rmA,1)-2);
        if mod(size(rmA,2),2) == 1
            offs_angle = (-floor(size(rmA,2)/2)):(floor(size(rmA,2)/2));
        else
            offs_angle = (-size(rmA,2)/2):(size(rmA,2)/2-1);
        end
        cc_plot = NaN(length(offs_dist), length(offs_angle));
        
        for i=1:length(offs_angle)
            rotB = circshift(rmB,offs_angle(i),2);
            for j=1:length(offs_dist)
                if offs_dist(j) < 0
                    Bpart = rotB((-offs_dist(j)+1):end,:);
                    Apart = rmA(1:(size(rmA,1)+offs_dist(j)),:);
                end
                if offs_dist(j) > 0
                    Bpart = rotB(1:size(rmB,1)-offs_dist(j),:);
                    Apart = rmA(offs_dist(j)+1:end,:);
                end
                if offs_dist(j) == 0
                    Apart = rmA;
                    Bpart = rotB;
                end
                Apart = Apart(:);
                Bpart = Bpart(:);
                not_nan = not(isnan(Apart)) .* not(isnan(Bpart));
                Apart = Apart(not_nan == 1);
                Bpart = Bpart(not_nan == 1);
                
                if length(Apart > 5) & length(Bpart > 5)
                    corr = corrcoef(Apart, Bpart);
                    cc_plot(j,i) = corr(1,2);
                end
            end
        end
    
        thresh_cc_plot = cc_plot > prctile(cc_plot,75, "all");
        thresh_cc_plot = imfill(thresh_cc_plot, 'holes');
        [labeled_cc_plot, ~] = bwlabel(thresh_cc_plot, 8);
        props = regionprops(labeled_cc_plot, cc_plot, 'all');
    
        bigEnough = vertcat(props.Area) > 9; % large enough to be a blob, currently 3x3 square bin
        allBlobCentroids = vertcat(props.Centroid);
        allBlobCentroids(:,1) = allBlobCentroids(:,1) + offs_angle(1)-1;
        allBlobCentroids(:,2) = allBlobCentroids(:,2) + offs_dist(1)-1;
        allBlobCentroids = allBlobCentroids(bigEnough,:);
        centDists = zeros(size(allBlobCentroids,1),1);
        for k=1:size(allBlobCentroids,1)
            centDists(k) = sqrt((allBlobCentroids(k,1))^2 + (allBlobCentroids(k,2))^2);
        end
        [~,whichBlob] = min(centDists);
        coordsWhichBlob = allBlobCentroids(whichBlob,:);
    
        subplot(nrows,ncols,6+include_spikemap*2);
        [x,y] = meshgrid(offs_angle, offs_dist);
        surface(x, y, (cc_plot)); shading interp
        
        hold on
    
        [~, ind] = max(cc_plot,[],'all');
        [row_max,col_max] = ind2sub(size(cc_plot),ind);
        %disp([row_max, col_max])
        %disp([length(offs_angle), length(offs_dist)])
        %disp([offs_angle(col_max), offs_dist(row_max)])
        scatter3(coordsWhichBlob(1), coordsWhichBlob(2),1,20,'black','filled');
    
        set(gca,'XTick',[],'YTick',[])
        ylim([-50,50])
        %colormap(slanCM('viridis'))
        colormap(parula)
        min_clim = mink((min(cc_plot)),2);
        max_clim = maxk((max(cc_plot)),2);
        set(gca, 'YDir','Normal','CLim',[min_clim(2),max_clim(2)])
        set(gca,'YDir','Normal')
        set(gca,'YDir','reverse')
        axis off
        %axis square
        posi = get(gca, 'Position');
        set(gca,'DataAspectRatio',[length(offs_angle) length(offs_dist) 1], 'Position',[posi(1),posi(2)-0.4,posi(3),posi(4)+0.8])
        %set(gca,'DataAspectRatio',[size(out.OFfull.rm_ns,2) size(out.OFfull.rm_ns,1) 1])
        %title('CC, even/odd', 'Position', [0 13.2])
        %colorbar
        set(gca,"FontSize",7)
    
        %% finaly of/chasing full
        rmA = out.OFfull.rm_ns(:,1:end);
        rmB = out.cfull.rm_ns(:,1:end);
        %rmA = out.OFfull.rm_ns(:,1:end-1);
        %rmB = out.cfull.rm_ns(:,1:end-1);
    
        offs_dist = (-size(rmA,1)+2):(size(rmA,1)-2);
        if mod(size(rmA,2),2) == 1
            offs_angle = (-floor(size(rmA,2)/2)):(floor(size(rmA,2)/2));
        else
            offs_angle = (-size(rmA,2)/2):(size(rmA,2)/2-1);
        end
        cc_plot = NaN(length(offs_dist), length(offs_angle));
        
        for i=1:length(offs_angle)
            rotB = circshift(rmB,offs_angle(i),2);
            for j=1:length(offs_dist)
                if offs_dist(j) < 0
                    Bpart = rotB((-offs_dist(j)+1):end,:);
                    Apart = rmA(1:(size(rmA,1)+offs_dist(j)),:);
                end
                if offs_dist(j) > 0
                    Bpart = rotB(1:size(rmB,1)-offs_dist(j),:);
                    Apart = rmA(offs_dist(j)+1:end,:);
                end
                if offs_dist(j) == 0
                    Apart = rmA;
                    Bpart = rotB;
                end
                Apart = Apart(:);
                Bpart = Bpart(:);
                not_nan = not(isnan(Apart)) .* not(isnan(Bpart));
                Apart = Apart(not_nan == 1);
                Bpart = Bpart(not_nan == 1);
                
                if length(Apart) > 5 & length(Bpart) > 5
                    corr = corrcoef(Apart, Bpart);
                    cc_plot(j,i) = corr(1,2);
                end
            end
        end
    
        thresh_cc_plot = cc_plot > prctile(cc_plot,75, "all");
        thresh_cc_plot = imfill(thresh_cc_plot, 'holes');
        [labeled_cc_plot, ~] = bwlabel(thresh_cc_plot, 8);
        props = regionprops(labeled_cc_plot, cc_plot, 'all');
    
        bigEnough = vertcat(props.Area) > 9; % large enough to be a blob, currently 3x3 square bin
        allBlobCentroids = vertcat(props.Centroid);
        allBlobCentroids(:,1) = allBlobCentroids(:,1) + offs_angle(1)-1;
        allBlobCentroids(:,2) = allBlobCentroids(:,2) + offs_dist(1)-1;
        allBlobCentroids = allBlobCentroids(bigEnough,:);
        centDists = zeros(size(allBlobCentroids,1),1);
        for k=1:size(allBlobCentroids,1)
            centDists(k) = sqrt((allBlobCentroids(k,1))^2 + (allBlobCentroids(k,2))^2);
        end
        [~,whichBlob] = min(centDists);
        coordsWhichBlob = allBlobCentroids(whichBlob,:);
    
        subplot(nrows,ncols,9+include_spikemap*2);
        [x,y] = meshgrid(offs_angle, offs_dist);
        surface(x, y, (cc_plot)); shading interp
        
        hold on
    
        [~, ind] = max(cc_plot,[],'all');
        [row_max,col_max] = ind2sub(size(cc_plot),ind);
        %disp([row_max, col_max])
        %disp([offs_angle(col_max), offs_dist(row_max)])
        scatter3(coordsWhichBlob(1), coordsWhichBlob(2),1,20,'black','filled');
    
        set(gca,'XTick',[],'YTick',[])
        ylim([-50,50])
        %colormap(slanCM('viridis'))
        colormap(parula)
        min_clim = mink((min(cc_plot)),2);
        max_clim = maxk((max(cc_plot)),2);
        set(gca, 'YDir','Normal','CLim',[min_clim(2),max_clim(2)])
        set(gca,'YDir','Normal') 
        set(gca,'YDir','reverse')
        axis off
        %axis square
        posi = get(gca, 'Position');
        set(gca,'DataAspectRatio',[length(offs_angle) length(offs_dist) 1], 'Position',[posi(1),posi(2)-0.4,posi(3),posi(4)+0.8])
        %set(gca,'DataAspectRatio',[size(out.OFfull.rm_ns,2) size(out.OFfull.rm_ns,1) 1])
        %title('CC, OF/chase', 'Position', [0 13.2])
        %colorbar
        set(gca,"FontSize",7)
    
        %% if hd colored firing where plot
        if include_spikemap
            ax_of = subplot(nrows,ncols,include_spikemap);
            hold on;
            plot(root.OFfull.x,root.OFfull.y,'Color',[.7 .7 .7])
            colormap(ax_of,hsv)
            xlim([min(root.OFfull.x) max(root.OFfull.x)]); ylim([min(root.OFfull.y) max(root.OFfull.y)])
            cx = root.OFfull.x(root.OFfull.spike == 1);
            cy = root.OFfull.y(root.OFfull.spike == 1);
            ch = root.OFfull.md(root.OFfull.spike == 1);
            %scatter(cx,cy,15,ch,'filled')
            scatter(cx,cy,7,ch,'filled')
            
            scatter(out.OFfull.QP(:,1),out.OFfull.QP(:,2),30,'k','filled')
            set(gca,'YDir','Normal')
            %title('Trajectory')
            %colorbar('Ticks', [-pi/2, 0, pi/2], 'TickLabels',{0,90,180})
            axis off
            axis square
    
            ax_ch = subplot(nrows,ncols,include_spikemap+3+include_cross_corr);
            hold on;
            plot(root.OFfull.x,root.OFfull.y,'Color',[.7 .7 .7])
            colormap(ax_ch,hsv)
            xlim([min(root.cfull.x) max(root.cfull.x)]); ylim([min(root.cfull.y) max(root.cfull.y)])
            cx = root.cfull.x(root.cfull.spike == 1);
            cy = root.cfull.y(root.cfull.spike == 1);
            ch = root.cfull.md(root.cfull.spike == 1);
            scatter(cx,cy,7,ch,'filled')
            
            scatter(out.cfull.QP(:,1),out.cfull.QP(:,2),30,'k','filled')
            set(gca,'YDir','Normal')
            %title('Trajectory')
            %colorbar('Ticks', [-pi/2, 0, pi/2], 'TickLabels',{0,90,180})
            axis off
            axis square
        end
    
        stuff = [coordsWhichBlob(1), coordsWhichBlob(2)];
    end

    
   






































