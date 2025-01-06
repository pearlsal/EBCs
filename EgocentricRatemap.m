function [out] = EgocentricRatemap(root, animal, varargin)

    % calculates the egocentric ratemap and returns structure containing all of
    % the relevant information for plotting or statistical testing. 
    %
    % Note: In order to simplify distribution we assume that the "root" object
    % contains cleaned and epoched data only. 
    if animal == "Luke"
        %luke chasing2
        edges = [-25.8353  118.1022;
        87.6440  118.1022;
        122.7823   73.7226;
        121.0541  -40.1460;
        85.3399  -67.0073;
        -36.2039  -66.4234;
        -69.0380  -29.6350;
        -64.4297   81.8978];
    end

    if animal == "Arwen"
        %Arwen OF1 box edges
        edges = [ -61.0253   68.2044;
        82.0622   69.6058;
        83.7212  -74.7445;
        -60.1959  -78.4818];
    end
    
    %OF edges
    if animal == "Tauriel"
        edges = [ -59.1129   63.5328;
        77.6613   64.0000;
        77.6613  -75.6788;
        -59.4355  -75.6788];
    end
    
    %chasing edges
    if animal == "Tauriel"
        edges = [ -61.8548   67.2701;
        82.4770   65.4015;
        82.8917  -76.1460;
        -61.8548  -75.6788];
    end

    %% Setup and parse:
    p = inputParser;
    degreesBin = 10;
    distsBin = 5;
    maxDistsBin = 60;
    if isfield(root,'bait_angle') % when doing bait maps, increase max dist to 90cm
        disp("bait rate maps")
        degreesBin = 20;
        distsBin = 4.5;
        maxDistsBin = 90;
    end
    p.addParameter('videoSamp', 1);                    % calculate every X frames of video
    p.addParameter('degSamp', degreesBin);                      % Degree bins
    p.addParameter('distanceBins', 0:distsBin:maxDistsBin);          % How far to look (cm)
    p.addParameter('thetaBins', -180:degreesBin:180-degreesBin);             % Specify angle binsize here, not autom. set to 360
    p.addParameter('boundaryMode', edges);                 % 0-> autolines, 1-> click, mat->useit
    p.addParameter('smoothKernel', [1 1 0])
%     p.addParameter('sampRate', 30)                     % Samples per second
    p.parse(varargin{:});

    fn = fieldnames(p.Results);
    for i = 1:length(fn)
        eval([fn{i} '=' mat2str(p.Results.(fn{i})) ';']);
    end  

    %% Unpack behavioral information
    if strcmp(class(root), 'CMBHOME.Session')
        rx = CMBHOME.Utils.ContinuizeEpochs(root.x);% * root.spatial_scale;
        ry = CMBHOME.Utils.ContinuizeEpochs(root.y);% * root.spatial_scale;
        md = CMBHOME.Utils.ContinuizeEpochs(root.headdir);
        if max(md) > 2*pi
            md = deg2rad(md);
        end
        ts = CMBHOME.Utils.ContinuizeEpochs(root.ts);
        sts = CMBHOME.Utils.ContinuizeEpochs(root.cel_ts);
        spk = histc(sts, ts);
        distanceBins = distanceBins/root.spatial_scale;
    else
        rx = root.x;      % x position in cm or pixels
        ry = root.y;      % y position in cm or pixels
        md = root.md;     % movement (or head) direction, in radians
        ts = root.ts;     % time stamps (seconds)
        spk = root.spike; % binary spike train
        MRL_distr_comp = root.mrl; % If you want to compute MRL_distribution, 1 if yes
        MI_distr_comp = root.mi; % If you want to compute skaggs_MI_rate, 1 if yes
        if isfield(root, 'cc')
            CC_distr_comp = root.cc; % If you want to compute
        else
            CC_distr_comp = 0;
        end
        
        if isfield(root,'bait_angle')
            disp("Doing bait bearing plots")
            bait_bearing = 1;
            bait_dist = root.bait_dist;
            bait_angle = root.bait_angle;
        else
            bait_bearing = 0;
        end
    end

    %% Get structure of the environnment
    if numel(boundaryMode)==1
        if boundaryMode == 0 
            % Headless / auto. Automatically detect the edges. 
            % Only works for a rectangular box with no insertions.

            p = [-1000 -1000];
            d = (rx-p(1)).^2 + (ry-p(2)).^2;
            [~,ind] = min(d);
            ll = [rx(ind) ry(ind)];

            p = [1000 -1000];
            d = (rx-p(1)).^2 + (ry-p(2)).^2;
            [~,ind] = min(d);
            lr = [rx(ind) ry(ind)];

            p = [1000 1000];
            d = (rx-p(1)).^2 + (ry-p(2)).^2;
            [~,ind] = min(d);
            ur = [rx(ind) ry(ind)];

            p = [-1000 1000];
            d = (rx-p(1)).^2 + (ry-p(2)).^2;
            [~,ind] = min(d);
            ul = [rx(ind) ry(ind)];


            QP = [ll;lr;ur;ul];
        elseif boundaryMode == 1
            % Finds edges by selecting corners of the environment
            QP = findEdges(root);
            disp(QP)
        end
    else
        % was fed in manually
        QP = boundaryMode;
    end
    disp("pre stuff done")
    %disp(QP)
    %% Calculate distances
    %[dis, ex, ey] = calcDistance(rx,ry,md, QP, degSamp);
    % Need to change this to always calculate for all the distances,
    % because Andy's code does not average across bins........lol
    [dis, ex, ey] = calcDistance(rx,ry,md, QP, 1);
    disp("distances done")
    meanFilterFunction = @(theBlockStructure) mean(theBlockStructure.data(:));
    dis = blockproc(dis(:,1:end-1), [1, degSamp], meanFilterFunction);
    %disp(dis(1,:))
    %disp(size(dis))
    %pause;
    disp("rebinning distances done")
    %% Calculate raw maps:
    if bait_bearing % bait rate map, not boundary
        thetaBins_ = deg2rad(-180:degreesBin:180);
        thetaBins = deg2rad(thetaBins);
        
        occ = NaN(length(thetaBins), length(distanceBins));
        nspk = occ;
        distanceBins(end+1) = Inf;
        ci = find(spk);
        
        for i = 1:length(thetaBins_)-1
            t = bait_angle>=thetaBins_(i) & bait_angle<thetaBins_(i+1);
            for k = 1:length(distanceBins)-1
                inds = bait_dist>=distanceBins(k) & bait_dist<distanceBins(k+1);
                goods = intersect(find(t), find(inds));
                occ(i,k) = length(goods);
                nspk(i,k) = length(intersect(goods,ci));
            end
        end
        distanceBins = distanceBins(1:end-1);
        %thetaBins = thetaBins(1:end-1);
        
        occ = occ(:,1:end-1);
        nspk = nspk(:,1:end-1);
        occ=occ';
        nspk=nspk';
    else
        thetaBins = deg2rad(thetaBins);
        occ = NaN(length(thetaBins), length(distanceBins));
        nspk = occ;
        distanceBins(end+1) = Inf;
        
        ci = find(spk);
    
        for i = 1:length(thetaBins)
            t = dis(:,i);
            for k = 1:length(distanceBins)-1
                inds = t>=distanceBins(k) & t<distanceBins(k+1);
                occ(i,k) = sum (inds);
                inds = find(inds);
                nspk(i,k) = length(intersect(inds,ci));
            end
        end
        distanceBins = distanceBins(1:end-1);
        %Don't need these after all, rather remove the repeat at the end
        % altogether
        %occ(end+1,:) = occ(1,:);
        %nspk(end+1,:) = nspk(1,:);
        % bring back to original dims
        occ = occ(:,1:end-1); occ=occ';
        nspk = nspk(:,1:end-1); nspk=nspk';
        
    end

    occ_ns = occ;
    nspk_ns = nspk;

    rm_ns = (nspk./occ); % non-smoothed ratemap
    disp("raw maps done")
    %% Smoothing
    %occ(:,end+1) = occ(:,1);
    %nspk(:,end+1) = nspk(:,1);
    %Don't add the -180 elem to the end of 180, I think it's an oversight
    
    occ = [occ occ occ];
    nd = numel(thetaBins);
    %nd = numel(thetaBins+1);
    occ = SmoothMat(occ, smoothKernel(1:2), smoothKernel(3));
    occ = occ(:, nd+1:2*nd);

    nspk = [nspk nspk nspk];
    nspk = SmoothMat(nspk,smoothKernel(1:2),smoothKernel(3));   % Smooth it
    nspk = nspk(:,nd+1:2*nd);                       % bring it back
    disp('num spikes');
    disp(nspk);
    disp('occupancy');
    disp(occ);
    
    %occ = occ(:,1:end-1);
    %nspk = nspk(:,1:end-1);
    rm = (nspk./occ);
    disp("smoothing done")

    %% ADDED ANDY's EBVC TEST-FUNCTIONALITY HERE
    %addpath('/Users/martibma/Documents/MATLAB/Packages/CMBHOME')
    %import CMBHOME.Utils.circ.*
    NFR = (sum(nspk) ./ sum(occ));% * root.fs_video;
    
    MRL = circ_r(thetaBins, NFR, 2*pi/(length(thetaBins)), 2);
    %MRL = circ_r(-pi:2*pi/(length(thetaBins)-1):pi, NFR, 2*pi/(length(thetaBins)-1), 2);
    
    PO = circ_mean(thetaBins, NFR, 2);
    %PO = circ_mean(-pi:2*pi/(length(thetaBins)-1):pi, NFR, 2);
    
    a = thetaBins - PO;
    [~, thetaPeak] = min(abs(a));
    
%   Fit only along preferred orientation
            
%   transform to weibel distribution of distance
    y=rm(:,thetaPeak); y = y-min(y); 
    y = y+0.00000000000001;
    y= y/sum(y);
    x = distanceBins(2:end); x = x(:); y = y(:);
    disp('wtf');
    disp(y);
    disp(x);

    fr = y;
    ft = fittype( 'weibull' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0];
    
%   Fit model to data.
try
    [fitresult, ~] = fit(x, y, ft, opts ); 
%   get weibull peak
    
    disp('wtf2');
    disp(y);
    disp(x);
    disp(ft);
    disp(opts);

    WeibullFitX = 0:0.01:max(distanceBins);
    WeibullFit = fitresult(WeibullFitX);
    [~,PrefDistance] = max(WeibullFit);
    PD = WeibullFitX(PrefDistance);
end

if ~exist('PrefDistance','var')
    PrefDistance = NaN;
    PD = NaN;
end
%     }

    %% Compute MRL distribution - by Martin
    if MRL_distr_comp == 1
        % Calculate raw maps:
        disp("Starting MRL distribution computation")
        rng(42)
        distanceBins_p = distanceBins;
        distanceBins_p(end+1) = Inf;
        MRL_distr = zeros(1,100);
        for j = 1:100
            occ_ = NaN(length(thetaBins), length(distanceBins));
            nspk_ = occ_;
            %distanceBins(end+1) = Inf;
            
            %ci_ = find(circshift(spk,randi(length(spk))));
            %ci_ = find(circshift(spk,randi([1200,length(spk)]))); % this one shifts by at least 30sec with 25ms bins
            ci_ = find(circshift(spk,randi([3750,length(spk)]))); % this one shifts by at least 30sec with 8ms bins

            for i = 1:length(thetaBins)
                t = dis(:,i);
                for k = 1:length(distanceBins_p)-1
                    inds_ = t>=distanceBins_p(k) & t<distanceBins_p(k+1);
                    occ_(i,k) = sum (inds_);
                    inds_ = find(inds_);
                    nspk_(i,k) = length(intersect(inds_,ci_));
                end
            end
            %distanceBins = distanceBins(1:end-1);
    
            % bring back to original dims
            occ_ = occ_(:,1:end-1); occ_=occ_';
            nspk_ = nspk_(:,1:end-1); nspk_=nspk_';
           
            %disp("raw maps done")
            % Smoothing
            
            occ_ = [occ_ occ_ occ_];
            nd = numel(thetaBins);
            occ_ = SmoothMat(occ_, smoothKernel(1:2), smoothKernel(3));
            occ_ = occ_(:, nd+1:2*nd);
    
            nspk_ = [nspk_ nspk_ nspk_];
            nspk_ = SmoothMat(nspk_,smoothKernel(1:2),smoothKernel(3));   % Smooth it
            nspk_ = nspk_(:,nd+1:2*nd);                       % bring it back
        
            %rm = (nspk_./occ_);
            %disp("smoothing done")
            
            NFR_ = (sum(nspk_) ./ sum(occ_));% * root.fs_video;
            MRL_ = circ_r(thetaBins, NFR_, 2*pi/(length(thetaBins)), 2);
            %MRL_ = circ_r(-pi:2*pi/(length(thetaBins)-1):pi, NFR_, 2*pi/(length(thetaBins)-1), 2);
            MRL_distr(j) = MRL_;
            disp("iteration "+j+ "done")
        end
    else
        MRL_distr = 0;
    end

     %% Compute MI distribution
    if MI_distr_comp == 1
        % Calculate raw maps:
        disp("Starting MI distribution computation")
        rng(42)
        distanceBins_p = distanceBins;
        distanceBins_p(end+1) = Inf;
        MI_distr = zeros(1,100);
        for k = 1:100
            thetaBins_ = deg2rad(-180:degreesBin:180);
        
            occ_ = NaN(length(thetaBins), length(distanceBins));
            nspk_ = occ_;
            %distanceBins(end+1) = Inf;
            ci_ = find(circshift(spk,randi([3750,length(spk)]))); % this one shifts by at least 30sec with 8ms bins
            
            for i = 1:length(thetaBins_)-1
                t = bait_angle>=thetaBins_(i) & bait_angle<thetaBins_(i+1);
                for l = 1:length(distanceBins_p)-1
                    inds_ = bait_dist>=distanceBins_p(l) & bait_dist<distanceBins_p(l+1);
                    goods_ = intersect(find(t), find(inds_));
                    occ_(i,l) = length(goods_);
                    nspk_(i,l) = length(intersect(goods_,ci_));
                end
            end
            %distanceBins = distanceBins(1:end-1);
            %thetaBins = thetaBins(1:end-1);
            
            occ_ = occ_(:,1:end-1);
            nspk_ = nspk_(:,1:end-1);
            occ_=occ_';
            nspk_=nspk_';  
            
            occ_ns_ = occ_;
            nspk_ns_ = nspk_;
            % Smoothing
            
            occ_ = [occ_ occ_ occ_];
            nd = numel(thetaBins);
            occ_ = SmoothMat(occ_, smoothKernel(1:2), smoothKernel(3));
            occ_ = occ_(:, nd+1:2*nd);
    
            nspk_ = [nspk_ nspk_ nspk_];
            nspk_ = SmoothMat(nspk_,smoothKernel(1:2),smoothKernel(3));   % Smooth it
            nspk_ = nspk_(:,nd+1:2*nd);                       % bring it back
        
            rm_ns_ = (nspk_ns_./occ_ns_);
            %disp("smoothing done")
            
            overall_mean_fr_ = 0;
            p_occ_ = zeros(length(distanceBins)-1, length(thetaBins));
            
            for i=1:length(distanceBins)-1
                for j=1:length(thetaBins)
                    if occ_ns_(i,j) > 50
                        p_occ_(i,j) = occ_ns_(i,j)/sum(occ_ns_(occ_ns_ > 50), 'all') + 0.;
                        overall_mean_fr_ = overall_mean_fr_ + p_occ_(i,j)*rm_ns_(i,j);
                    end
                end
            end
            MI_ = 0;
            for i=1:length(distanceBins)-1
                for j=1:length(thetaBins)-1
                    if occ_ns_(i,j) > 50
                        if rm_ns_(i,j) > 0.000001
                            MI_ = MI_ + p_occ_(i,j)*rm_ns_(i,j)/overall_mean_fr_*log2(rm_ns_(i,j)/overall_mean_fr_);
                        end
                    end
                end
            end
            MI_distr(k) = MI_;
            disp("iteration "+k+ "done")
        end
    else
        MI_distr = 0;
    end

    % calculate skaggs mutual information rate
    overall_mean_fr = 0;
    p_occ = zeros(length(distanceBins)-1, length(thetaBins));
    for i=1:length(distanceBins)-1
        for j=1:length(thetaBins)
            if occ_ns(i,j) > 50 % threshold at 50
                p_occ(i,j) = occ_ns(i,j)/sum(occ_ns(occ_ns > 50), 'all') + 0.;
                overall_mean_fr = overall_mean_fr + p_occ(i,j)*rm_ns(i,j);
            end
        end
    end
    MI = 0;
    for i=1:length(distanceBins)-1
        for j=1:length(thetaBins)
            if occ_ns(i,j) > 50 % threshold at 50
                if rm_ns(i,j) > 0.000001
                    MI = MI + p_occ(i,j)*rm_ns(i,j)/overall_mean_fr*log2(rm_ns(i,j)/overall_mean_fr);
                end
            end
        end
    end
    
    if CC_distr_comp == 1
        if bait_bearing == 0
            % compute 100 shifted rate maps
            disp("Starting CC shifted ratemap computation")
            rng(42)
            distanceBins_p = distanceBins;
            distanceBins_p(end+1) = Inf;
            shifted_rate_map = zeros(size(rm,1), size(rm,2), 100);
            for j = 1:100
                occ_ = NaN(length(thetaBins), length(distanceBins));
                nspk_ = occ_;
                %distanceBins(end+1) = Inf;
                
                %ci_ = find(circshift(spk,randi(length(spk))));
                %ci_ = find(circshift(spk,randi([1200,length(spk)]))); % this one shifts by at least 30sec with 25ms bins
                ci_ = find(circshift(spk,randi([3750,length(spk)]))); % this one shifts by at least 30sec with 8ms bins
    
                for i = 1:length(thetaBins)
                    t = dis(:,i);
                    for k = 1:length(distanceBins_p)-1
                        inds_ = t>=distanceBins_p(k) & t<distanceBins_p(k+1);
                        occ_(i,k) = sum (inds_);
                        inds_ = find(inds_);
                        nspk_(i,k) = length(intersect(inds_,ci_));
                    end
                end
                %distanceBins = distanceBins(1:end-1);
        
                % bring back to original dims
                occ_ = occ_(:,1:end-1); occ_=occ_';
                nspk_ = nspk_(:,1:end-1); nspk_=nspk_';
                rm_ns_ = (nspk_./occ_); % non-smoothed ratemap
                shifted_rate_map(:,:,j) = rm_ns_;
                disp("iteration "+j+ "done")
            end
        else
            % compute 100 shifted rate maps
            disp("Starting CC shifted ratemap computation")
            rng(42)
            distanceBins_p = distanceBins;
            distanceBins_p(end+1) = Inf;
            shifted_rate_map = zeros(size(rm,1), size(rm,2), 100);
            for k = 1:100
                thetaBins_ = deg2rad(-180:degreesBin:180);
            
                occ_ = NaN(length(thetaBins), length(distanceBins));
                nspk_ = occ_;
                %distanceBins(end+1) = Inf;
                ci_ = find(circshift(spk,randi([3750,length(spk)]))); % this one shifts by at least 30sec with 8ms bins
                
                for i = 1:length(thetaBins_)-1
                    t = bait_angle>=thetaBins_(i) & bait_angle<thetaBins_(i+1);
                    for l = 1:length(distanceBins_p)-1
                        inds_ = bait_dist>=distanceBins_p(l) & bait_dist<distanceBins_p(l+1);
                        goods_ = intersect(find(t), find(inds_));
                        occ_(i,l) = length(goods_);
                        nspk_(i,l) = length(intersect(goods_,ci_));
                    end
                end
                %distanceBins = distanceBins(1:end-1);
                %thetaBins = thetaBins(1:end-1);
                
                occ_ = occ_(:,1:end-1);
                nspk_ = nspk_(:,1:end-1);
                occ_=occ_';
                nspk_=nspk_';  
                rm_ns_ = (nspk_./occ_); % non-smoothed ratemap
                rm_ns_(occ_ < 50) = nan; % filter based on the occurance of original map
                shifted_rate_map(:,:,k) = rm_ns_;
                disp("iteration "+k+ "done")
            end
        end
    else
        shifted_rate_map = 0;
    end

    %% package the output
    out.rm_ns = rm_ns;
    out.occ_ns = occ_ns;
    out.nspk_ns = nspk_ns;
    out.occ = occ;
    out.nspk = nspk;
    out.rm = rm;
    out.QP = QP;
    
    out.params.videoSamp = videoSamp;
    out.params.degSamp = degSamp;
    out.params.distanceBins = distanceBins;
    out.params.smoothKernel = smoothKernel;
    out.params.thetaBins = thetaBins;

    out.MRL = MRL; % mean resultant length
    out.PrefOrient = PO; % preferred orientation
    out.PrefDist = PD; % preferred distance
    out.MRL_dist = MRL_distr; % MRL distribution
    out.MI = MI; % mutual information rate
    out.MI_dist = MI_distr; % MI distribution

    out.CCrm_shift = shifted_rate_map; % 100 shifted ratemaps
    
end

%%
function QP = findEdges(root)
    ifEscape = 0;
    h=figure();

    while ~ifEscape  
        figure(h); 
        clf
        
        set(gca,'YDir','Normal'); %colormap(jet);
        clim=get(gca,'clim');set(gca,'clim',clim/50);
        hold on
        plot(root.x, root.y,'k');
        QP = [];
        
        set(h,'Name','Select Corners of Walls. Esc--> done. **Do not complete!**')

        button = 1;

        while button~=27
            [x,y,button] = ginput(1);

            clf
            
            set(gca,'YDir','Normal'); %colormap(jet);
            clim=get(gca,'clim');set(gca,'clim',clim/50);
            hold on
            plot(root.x, root.y,'k');
            
            if ~isempty(QP)
                plot(QP(:,1),QP(:,2),'r')
                plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')
            end

            if button == 32 %space bar
                QP = [QP; NaN NaN];
            elseif button~=27
                QP = [QP; x y];
            end

            plot(QP(:,1),QP(:,2),'r')
            plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')

        end

        %Ask for verification
        edg = splitter(QP);
        clf;
        set(h,'Name','Verify. 0--> Try again; 1--> Confirm')
        plot(root.x, root.y,'k');
        hold on
        
        for m = 1:numel(edg)
            for n = 1:size(edg{m},1)
                sp = squeeze(edg{m}(n,:,1));
                ep = squeeze(edg{m}(n,:,2));
                plot([sp(1) ep(1)],[sp(2) ep(2)],'ro','MarkerFaceColor','r')
                plot([sp(1) ep(1)],[sp(2) ep(2)],'r')
            end
        end

        % set or repeat
        while button ~=48 && button~=49
            [~,~,button]=ginput(1);
        end
        ifEscape = button==49;

    end

    close(h);
    drawnow();
end

%%
function edg = splitter(QP)
    
    inds = find(isnan(QP(:,1)));
    xs=SplitVec(QP(:,1), @(x) isnan(x));
    ys=SplitVec(QP(:,2), @(x) isnan(x));
    
    % split corners
    for m = 1:size(xs,1)
        QP2{m} = [xs{m} ys{m}];
        QP2{m}(find(isnan(QP2{m}(:,1))),:) = [];
    end
    
    for m = 1:numel(QP2)
        for n = 1:size(QP2{m},1)
            sp = n;ep=n+1;
            if ep>size(QP2{m},1), ep=1;end
            edg{m}(n,:,1) = [QP2{m}(sp,1) QP2{m}(sp,2)];
            edg{m}(n,:,2) = [QP2{m}(ep,1) QP2{m}(ep,2)];
        end
    end

end

%% 
function [dis, ex, ey] = calcDistance(rx,ry,md, QP, degSamp_)
    
    mxd = sqrt((max(rx)-min(rx))^2 + (max(ry)-min(ry))^2);
    degs = deg2rad(-180:degSamp_:180);
        
    edg = splitter(QP);
    edg = cell2mat(edg(:));
    dis = NaN(numel(rx),size(edg,1), numel(degs));
    dir = dis;
    
    for i = 1:size(edg,1)
        x1=edg(i,1,1);x2=edg(i,1,2);
        y1=edg(i,2,1);y2=edg(i,2,2);
        for h = 1:numel(degs)
            mdof=degs(h);
            y3=ry;x3=rx;
            y4=ry+mxd*sin(md+mdof);
            x4=rx+mxd*cos(md+mdof);
            
            % find the intersection analytically
            px1 = (x1.*y2-y1.*x2).*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4);
            px2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
            px  = px1./px2;
            
            py1 = (x1.*y2-y1.*x2).*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4);
            py2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
            py = py1./py2;

            d = sqrt((ry-py).^2 + (rx-px).^2);
            dis(:,i,h) = d;
            
            % need to filter down to the right direction ...
            dir(:,i,h) = wrapToPi(atan2(py-ry,px-rx)-(md+mdof));
            
            % filter by bounding box
            bb = [min([x1 x2]) max([x1 x2]); min([y1 y2]) max([y1 y2])];  
            % |xmin, xmax|
            % |ymin, ymax|
            indexes = ~(px>=bb(1,1) & px<=bb(1,2) & py>=bb(2,1) & py<=bb(2,2));
            dis(indexes,i,h) = NaN;
        end
        
    end
    
    
    dis(dis>mxd) = NaN;
    dis(abs(dir)>pi/4) = NaN;
    
    %output
    dis=squeeze(nanmin(dis,[],2));
    dd=repmat(degs,size(rx,1),1) + repmat(md,1,numel(degs));
    dx=dis.*cos(dd); dy=dis.*sin(dd);
    ey=dy+repmat(ry,1,numel(degs));
    ex=dx+repmat(rx,1,numel(degs));
    
end

%%
function mat = SmoothMat(mat, kernel_size, std)
    %
    % Smooths matrix by convolving with 2d gaussian of size
    % kernel_size=[bins_x bins_y] and standard deviation 'std'
    %
    % if std==0, just returns mat


    if nargin<3
        std=1;
    end

    if std == 0, return; end

    [Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));

    kernel = pdf('Normal', Rgrid, 0, std);
    kernel = kernel./sum(sum(kernel));
    mat = conv2(mat, kernel, 'same');

end

