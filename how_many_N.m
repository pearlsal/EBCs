%% Read data and do the EBC test computations

folder_loc = '/Users/pearls/Work/RSC_project/';
which_animal = "PreciousGrape"; %{'Luke', 'Arwen', 'Tauriel'} "PreciousGrape"; %{'Luke', 'Arwen', 'Tauriel', 'PreciousGrape'}
which_session = "c1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = [ "c1", "c2"];
chase_or_chill = "chase";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's, or filter specific parts of sessions (chasing bouts)
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
EBC_or_EBOC = "EBOC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait

CC_stability = true; % compute shifted ratemaps, for using cross-corrs. to determine stability, not peak in ratemap like Andy does
%%
time_intervals = 0;
tracking_interval = 0.008332827537658849;
if shared_cells
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binnedshareddata', string(binsize*1000),'ms.mat'));
end

if not(shared_cells)
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binneddata', string(binsize*1000),'ms.mat'));
end

if which_session == "c1"
    disp(which_session)
    if chase_or_chill == "chase"
        c_intervals = [1340,5674,6283,9610,19185,23374,37078,50685,60603,61954,81994,89642,95876,97557,106424,109134,109938,112715,118202,120445,126240,136080,140506,146529,147420,152989,160061,172408];
    end
    if chase_or_chill == "chill"
        c_intervals = [1,1340,5674,6283,9610,19185,23374,37078,50685,60603,61954,81994,89642,95876,97557,106424,109134,109938,112715,118202,120445,126240,136080,140506,146529,147420,152989,160061,172408,size(data.spikemat,2)];
    end


    original_tracking_interval = 0.008332881468650438;
    new_tracking_interval = binsize;
    interval_diff = floor(new_tracking_interval/original_tracking_interval);
    if interval_diff == 0
        interval_diff = 1;
    end
    
    interval_bins_c = [];
    for j=1:length(c_intervals)
        if mod(j,2) == 1
            c_intervals(j) = floor(c_intervals(j)/interval_diff);
        end
        if mod(j,2) == 0
            c_intervals(j) = ceil(c_intervals(j)/interval_diff);
            interval_bins_c = cat(2, interval_bins_c, c_intervals(j-1):c_intervals(j));
        end
    end
    
    x_pos = data.binned_pos(interval_bins_c,1) * 100;
    y_pos = data.binned_pos(interval_bins_c,2) * 100;
    head_dir = data.binned_hd(interval_bins_c);
    speed = data.binned_speed(interval_bins_c);
    spikemat = data.spikemat(:,interval_bins_c);
    if EBC_or_EBOC == "EBOC"
        bait_angle = data.binned_rel_ha(interval_bins_c);
        bait_dist = data.binned_rel_dist(interval_bins_c);
    end
    
    session = repmat(which_session, size(data.binned_pos,1), 1);
    behav = repmat(chase_or_chill, size(data.binned_pos,1), 1);
    bin_number = interval_bins_c;
    
    data_length = size(data.spikemat(:, interval_bins_c),2);    % Ensure time_intervals is initialized
    if isempty(time_intervals)
        time_intervals = 0; % Or a logical starting point
    end
    time_intervals = time_intervals(end):tracking_interval:(tracking_interval*size(data.spikemat(:,interval_bins_c),2)+time_intervals(end));
else
    disp(which_session)
    x_pos = data.binned_pos(:,1) * 100;
    y_pos = data.binned_pos(:,2) * 100;
    head_dir = data.binned_hd;
    speed = data.binned_speed;
    spikemat = data.spikemat;
    if EBC_or_EBOC == "EBOC"
        bait_angle = data.binned_rel_ha;
        bait_dist = data.binned_rel_dist;
    end
    

    session = repmat(which_session, size(data.binned_pos,1), 1);
    behav = repmat("OF", size(data.binned_pos,1), 1);
    bin_number = 1:size(data.binned_pos,1);
    
    data_length = size(data.spikemat,2);
    tracking_interval = 0.008332827537658849;
    time_intervals = 0.0042:tracking_interval:tracking_interval*data_length;
    N = size(data.spikemat,1);
end


if not(concat_which(1) == which_session)
    disp("Remember to include the original session in the concat_list")
    exit()
end

if concat_sessions
    for i=2:length(concat_which)
        data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',concat_which(i),'_binnedshareddata', string(binsize*1000),'ms.mat'));
        if different_sessions & not(concat_which(i) == "OF1") & not(concat_which(i) == "OF2")
            if which_animal == "Tauriel"
                if concat_which(i) == "c1"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [1308,4932,17604,23452,32330,34869,39026,55359,65412,70414];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,1308,4932,17604,23452,32330,34869,39026,55359,65412,70414,size(data.spikemat,2)];
                    end
                end
    
                if concat_which(i) == "c2"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [2848,12662,23006,36415,54592,56622,58815,61079,72111,75756];
                        c_intervals = [2848,12662,23006,36415,54592,56622,72111,75756];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,2848,12662,23006,36415,54592,56622,72111,75756,size(data.spikemat,2)];
                        % doesn't chase in the 4th bout, can use that here instead
                    end
                end
    
                if concat_which(i) == "c4"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [1045,24906,41991,50003,53071,85051];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,1045,24906,41991,50003,53071,85051,size(data.spikemat,2)];
                    end
                end
    
                if concat_which(i) == "c5"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [3266,22481,29292,32927,36211,41067,59829,80588];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,3266,22481,29292,32927,36211,41067,59829,80588,size(data.spikemat,2)];
                    end
                end
            end

            if which_animal == "Arwen"
                if concat_which(i) == "c1"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [1610,7583,27805,38704,53618,57788];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,1610,7583,27805,38704,53618,57788,size(data.spikemat,2)];
                    end
                end

                if concat_which(i) == "c2"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [3996,4946,13087,20357,26277,35577,44556,47376,49806,53276,54831,61091,62611,65211];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,3996,4946,13087,20357,26277,35577,44556,47376,49806,53276,54831,61091,62611,65211,size(data.spikemat,2)];
                    end
                end

                if concat_which(i) == "c4"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [2490,9030,18678,23108,28007,29057,34017,36677,39773,43193,53881,60621,73509,77469];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,2490,9030,18678,23108,28007,29057,34017,36677,39773,43193,53881,60621,73509,77469,size(data.spikemat,2)];
                    end
                end
            end
            if which_animal == "PreciousGrape"
                if concat_which(i) == "c1"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [1340,5674,6283,9610,19185,23374,37078,50685,60603,61954,81994,89642,95876,97557,106424,109134,109938,112715,118202,120445,126240,136080,140506,146529,147420,152989,160061,172408];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,1340,5674,6283,9610,19185,23374,37078,50685,60603,61954,81994,89642,95876,97557,106424,109134,109938,112715,118202,120445,126240,136080,140506,146529,147420,152989,160061,172408,size(data.spikemat,2)];
                    end
                end

                if concat_which(i) == "c2"
                    disp(concat_which(i))
                    if chase_or_chill == "chase"
                        c_intervals = [22980,32707,37184,44892,52420,59100,70679,83959,88277,89508,103013,108722,119041,120567,127260,142127,145257,152810,158931,172812,181107,190995,192010,206903,212618,227531,229729,231482,244829,259379,275681,283632,283987,289621];
                    end
                    if chase_or_chill == "chill"
                        c_intervals = [1,22980,32707,37184,44892,52420,59100,70679,83959,88277,89508,103013,108722,119041,120567,127260,142127,145257,152810,158931,172812,181107,190995,192010,206903,212618,227531,229729,231482,244829,259379,275681,283632,283987,289621,size(data.spikemat,2)];
                    end
                end
            end 




            % Ensure c_intervals is defined before this point
            original_tracking_interval = 0.008332881468650438;
            new_tracking_interval = binsize;
            interval_diff = floor(new_tracking_interval/original_tracking_interval);
            if interval_diff == 0
                interval_diff = 1;
            end
            
            interval_bins_c = [];
            for j = 1:length(c_intervals)
                if mod(j, 2) == 1
                  c_intervals(j) = max(1, floor(c_intervals(j) / interval_diff));
                else
                c_intervals(j) = min(size(data.binned_pos, 1), ceil(c_intervals(j) / interval_diff));
                interval_bins_c = [interval_bins_c, c_intervals(j-1):c_intervals(j)];
                end
            end

            x_pos = cat(1, x_pos, data.binned_pos(interval_bins_c,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(interval_bins_c,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd(interval_bins_c));
            speed = cat(2, speed, data.binned_speed(interval_bins_c));
            spikemat = cat(2, spikemat, data.spikemat(:,interval_bins_c));
            if EBC_or_EBOC == "EBOC"
                bait_angle = cat(2, bait_angle, data.binned_rel_ha(interval_bins_c));
                bait_dist = cat(2, bait_dist, data.binned_rel_dist(interval_bins_c));
            end

            data_length = data_length + size(data.spikemat(:, interval_bins_c),2);
            new_time_intervals = time_intervals(end):tracking_interval:(tracking_interval * size(data.spikemat(:, interval_bins_c), 2) + time_intervals(end));
            new_time_intervals = new_time_intervals(new_time_intervals > 0); % Filter to ensure positive values
            time_intervals = cat(2, time_intervals, new_time_intervals);
        end
        
        if concat_which(i) == "OF1" | concat_which(i) == "OF2"
            disp(concat_which(i))
            x_pos = cat(1, x_pos, data.binned_pos(:,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(:,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd);
            speed = cat(2, speed, data.binned_speed);
            spikemat = cat(2, spikemat, data.spikemat);
            
            data_length = data_length + size(data.spikemat,2);
            time_intervals = cat(2, time_intervals, time_intervals(end):tracking_interval:(tracking_interval*size(data.spikemat,2)+time_intervals(end)));
        end
   end
end

% Filter out where there is nan hd
ind = ~isnan(head_dir);
x_pos = x_pos(ind);
y_pos = y_pos(ind);
head_dir = head_dir(ind);
spikemat = spikemat(:,ind);
time_intervals = time_intervals(ind);
speed = speed(ind);
if EBC_or_EBOC == "EBOC"
    bait_angle = bait_angle(ind);
    bait_dist = bait_dist(ind);
end

% Filter out where animal is standing still
if speed_threshold
    ind = speed > 5; % speed thresholding by 5cm/s
    x_pos = x_pos(ind);
    y_pos = y_pos(ind);
    head_dir = head_dir(ind);
    spikemat = spikemat(:,ind);
    time_intervals = time_intervals(ind);
    speed = speed(ind);
    if EBC_or_EBOC == "EBOC"
        bait_angle = bait_angle(ind);
        bait_dist = bait_dist(ind);
    end
end

% if doing bait plots, filter angle and dist nan too
if EBC_or_EBOC == "EBOC"
    ind = ~isnan(bait_angle);
    bait_angle = bait_angle(ind);
    bait_dist = bait_dist(ind);
    x_pos = x_pos(ind);
    y_pos = y_pos(ind);
    head_dir = head_dir(ind);
    spikemat = spikemat(:,ind);
    time_intervals = time_intervals(ind);
    speed = speed(ind);

    ind = ~isnan(bait_dist);
    bait_angle = bait_angle(ind);
    bait_dist = bait_dist(ind);
    x_pos = x_pos(ind);
    y_pos = y_pos(ind);
    head_dir = head_dir(ind);
    spikemat = spikemat(:,ind);
    time_intervals = time_intervals(ind);
    speed = speed(ind);

end
N = size(data.spikemat, 1);
%disp(size(speed,2)/data_length)
%test = [x_pos, y_pos];
%hist3(test, 'CdataMode', 'auto')
%colorbar
%view(2)

% Run test computations and save stuff for each neuron
for i=1:N
    disp("iteration number "+i)
    neuron_index = i;
    spike_vec = spikemat(neuron_index,:);
    
    r = struct;
    r.x = x_pos;
    r.y = y_pos;
    r.md = head_dir.';
    r.spike = spike_vec.';
    r.ts = time_intervals.';
    r.mrl = 1;
    r.mi = 0;
    if EBC_or_EBOC == "EBOC"
        r.bait_angle = bait_angle.';
        r.bait_dist = bait_dist.';
    end
    r.cc = 0;
    if CC_stability
        r.cc = 1;
    end
    out = EgocentricRatemap(r, which_animal);
    
    if speed_threshold
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(neuron_index), 'full.mat');
    end

    if not(speed_threshold)
         if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(neuron_index), 'full.mat');
    end

    save(filename, 'out');
    disp("done for EBC")
    %This part is the even odd 30 sek split, better use this than half the
    %data set, since we have so little
    bins_in_chunk = 3750; %30sek chunks
    chunks = floor(size(x_pos,1)/bins_in_chunk);
    r1 = struct;
    r1.x = [];
    r1.y = [];
    r1.md = [];
    r1.spike = [];
    r1.ts = [];
    r1.mrl = 0;
    r1.mi = 0;
    r1.cc = 0;
    r2 = struct;
    r2.x = [];
    r2.y = [];
    r2.md = [];
    r2.spike = [];
    r2.ts = [];
    r2.mrl = 0;
    r2.mi = 0;
    if EBC_or_EBOC == "EBOC"
        r1.bait_angle = [];
        r1.bait_dist = [];
        r2.bait_angle = [];
        r2.bait_dist = [];
    end
    r2.cc = 0;
    if CC_stability
        r1.cc = 1;
        r2.cc = 1;
    end

    for t=1:chunks-1
        if mod(t,2) == 0
            r1.x = cat(1, r1.x, x_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r1.y = cat(1, r1.y, y_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r1.md = cat(1, r1.md, head_dir(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.spike = cat(1, r1.spike, spike_vec(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.ts = cat(1, r1.ts, time_intervals(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            if EBC_or_EBOC == "EBOC"
                r1.bait_angle = cat(1, r1.bait_angle, bait_angle(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                r1.bait_dist = cat(1, r1.bait_dist, bait_dist(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            end
        end
        if (mod(t,2)) == 1
            r2.x = cat(1, r2.x, x_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r2.y = cat(1, r2.y, y_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r2.md = cat(1, r2.md, head_dir(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.spike = cat(1, r2.spike, spike_vec(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.ts = cat(1, r2.ts, time_intervals(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            if EBC_or_EBOC == "EBOC"
                r2.bait_angle = cat(1, r2.bait_angle, bait_angle(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                r2.bait_dist = cat(1, r2.bait_dist, bait_dist(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            end
        end
    end

    if (mod(t+1,2)) == 0
        r1.x = cat(1, r1.x, x_pos((t+1)*bins_in_chunk:end));
        r1.y = cat(1, r1.y, y_pos((t+1)*bins_in_chunk:end));
        r1.md = cat(1, r1.md, head_dir((t+1)*bins_in_chunk:end).');
        r1.spike = cat(1, r1.spike, spike_vec((t+1)*bins_in_chunk:end).');
        r1.ts = cat(1, r1.ts, time_intervals((t+1)*bins_in_chunk:end).');
        if EBC_or_EBOC == "EBOC"
            r1.bait_angle = cat(1, r1.bait_angle, bait_angle(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.bait_dist = cat(1, r1.bait_dist, bait_dist(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
        end
    end
    if (mod(t+1,2)) == 1
        r2.x = cat(1, r2.x, x_pos((t+1)*bins_in_chunk:end));
        r2.y = cat(1, r2.y, y_pos((t+1)*bins_in_chunk:end));
        r2.md = cat(1, r2.md, head_dir((t+1)*bins_in_chunk:end).');
        r2.spike = cat(1, r2.spike, spike_vec((t+1)*bins_in_chunk:end).');
        r2.ts = cat(1, r2.ts, time_intervals((t+1)*bins_in_chunk:end).');
        if EBC_or_EBOC == "EBOC"
            r2.bait_angle = cat(1, r2.bait_angle, bait_angle(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.bait_dist = cat(1, r2.bait_dist, bait_dist(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
        end
    end
    
    % DEPRECATED: DON'T USE FIRST AND LAST HALF, TOO DIFFERENT
    %half_i = round(length(x_pos)/2);
    %r1.x = x_pos(1:half_i);
    %r1.y = y_pos(1:half_i);
    %r1.md = head_dir(1:half_i).';
    %r1.spike = spike_vec(1:half_i).';
    %r1.ts = time_intervals(1:half_i).';
    %r1.mrl = 0;

    %r2.x = x_pos(half_i:end);
    %r2.y = y_pos(half_i:end);
    %r2.md = head_dir(half_i:end).';
    %r2.spike = spike_vec(half_i:end).';
    %r2.ts = time_intervals(half_i:end).';
    %r2.mrl = 0;

    out = EgocentricRatemap(r1, which_animal);
    if speed_threshold
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(neuron_index), 'firsthalf.mat');
    end

    if not(speed_threshold)
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(neuron_index), 'firsthalf.mat');
    end
    
    save(filename, 'out');

    out = EgocentricRatemap(r2, which_animal);
    if speed_threshold
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(neuron_index), 'lasthalf.mat');
    end

    if not(speed_threshold)
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(neuron_index), 'lasthalf.mat');
    end

    save(filename, 'out');

end
CC_stability = true; % use cross-corrs. to determine stability, not peak in ratemap like Andy does
test_nr2 = 2; %{1,2} 1=prefOrient from weibull fit, 2=center of mass of rate map
test_nr3 = 2; %{1,2,3} 1=prefDist from weibull fit, 2=center of mass of rate map, 3=peakBin of ratemap

if shared_cells
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binnedshareddata', string(binsize*1000),'ms.mat'));
end

if not(shared_cells)
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binneddata', string(binsize*1000),'ms.mat'));
end


N = size(data.spikemat,1);
EBC_array = zeros([N,1]);
N_EBC = 0;
cnt = 0;
tuned_p_vals = zeros([N,1]);
Stab_p_vals = zeros([N,1]);
for i=1:N
    if speed_threshold
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        if EBC_or_EBOC == "EBOC"
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/bait_oriented/neuron', string(i), 'full.mat');
 %           filename = "/Users/pearls/Work/RSC_project/EBC_results/neuron"+i+"full.mat";
            full_session = load(filename).out;
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/bait_oriented/neuron', string(i), 'firsthalf.mat');
            %filename = "/Users/pearls/Work/RSC_project/EBC_results/neuron"+i+"firsthalf.mat";
            firsthalf_session = load(filename).out;
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/bait_oriented/neuron', string(i), 'lasthalf.mat');
            %filename = "/Users/pearls/Work/RSC_project/EBC_results/neuron"+i+"lasthalf.mat";
            lasthalf_session = load(filename).out;
        else
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(i), 'full.mat');
            %filename = "/Users/pearls/Work/RSC_project/EBC_results/neuron"+i+"full.mat";
            full_session = load(filename).out;
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(i), 'firsthalf.mat');
            %filename = "/Users/pearls/Work/RSC_project/EBC_results/neuron"+i+"firsthalf.mat";
            firsthalf_session = load(filename).out;
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(i), 'lasthalf.mat');
            %filename = "/Users/pearls/Work/RSC_project/EBC_results/neuron"+i+"lasthalf.mat";
            lasthalf_session = load(filename).out;
        end
    end

    if not(speed_threshold)
        if concat_sessions
            session_name = "";
            for j=1:length(concat_which)
                session_name = session_name + concat_which(j);
                if not(j==length(concat_which))
                    session_name = session_name + "&";
                end
            end
        end
        if not(concat_sessions)
            session_name = which_session;
        end
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(i), 'full.mat');
        %filename = "EBC_results/neuron"+i+"full.mat";
        full_session = load(filename).out;
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(i), 'firsthalf.mat');
        %filename = "EBC_results/neuron"+i+"firsthalf.mat";
        firsthalf_session = load(filename).out;
        filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(i), 'lasthalf.mat');
        %filename = "EBC_results/neuron"+i+"lasthalf.mat";
        lasthalf_session = load(filename).out;
    end
    
    if EBC_or_EBOC == "EBOC"
        perc95 = maxk(full_session.MI_dist,6);
        perc95_value = perc95(end);
        distr = full_session.MI_dist;
        single_point = full_session.MI;
    else
        perc99 = maxk(full_session.MRL_dist,2);
        perc99 = perc99(2);
        perc95 = maxk(full_session.MRL_dist,6);
        perc95 = perc95(6);
        distr = full_session.MRL_dist;
        single_point = full_session.MRL;
    end

    %histogram(full_session.MRL_dist, 'FaceColor', "#1f77b4", 'FaceAlpha', 1);
    %hold on;
    %line([full_session.MRL,full_session.MRL], ylim, 'LineWidth', 3, 'Color', '#ff7f0e');
    %xlabel('MRL score')
    %ylabel('count')

    
    cdfplot(distr)
    set(findall(gca, 'Type', 'Line'),'LineWidth',5);
    hold on
    %line([full_session.MRL,full_session.MRL], ylim, 'LineWidth', 3, 'Color', '#ff7f0e');
    plot(single_point,1,'.', 'MarkerSize',70, 'Color','#ff7f0e')
    set(gca, 'box', 'off')
    grid off
    ylim([0,1.05])
    if single_point > max(distr)
        xlim([0,single_point+0.02])
    end
    view([90 -90])   
    set(gca,'YTick',[0,0.5,1])
    set(gca,'XTick',[gca().XTick(1),(gca().XTick(end)+gca().XTick(1))/2,gca().XTick(end)])
    ylabel('Cumulative dist. function')
    if EBC_or_EBOC  == "EBOC"
        xlabel('MI')
        filename_pdf = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/MI_distr', string(i), '.pdf');
        filename_png = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/MI_distr', string(i), '.png');
    else
        xlabel('MRL')
        filename_pdf = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/MRL_distr', string(i), '.pdf');
        filename_png = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/MRL_distr', string(i), '.png');
    end
    %pbaspect([2 1 1])
    title([])
    tuned_p_vals(i) = sum(single_point > sort(distr));

    
    %exportgraphics(gcf,filename_pdf, 'Resolution', 600)
    %exportgraphics(gcf, filename_png)
    %break
    close
    %set(gca,'YTick',[])

    c1 = single_point > perc95;
    if not(c1) % failed first criterion
        disp("Neuron number "+i+" failed first criterion with MRL "+full_session.MRL+" against 99th perc. "+perc99)
        cnt = cnt+1;
    end
    if true % passed the first criterion
        disp("Neuron number "+i+" passed first criterion")
        if CC_stability
            which_percentile = 95; % 95percentile of shuffled as criterion
            %CC_stability_mat = zeros(2,100);
            CC_shuffled = zeros(1,100);
            for k = 1:101
                if k == 101
                    rmA = firsthalf_session.rm_ns;
                    rmB = lasthalf_session.rm_ns;
                else
                    rmA = firsthalf_session.CCrm_shift(:,:,k);
                    rmB = lasthalf_session.CCrm_shift(:,:,k);
                end
                
                if k == 101
                    perc = maxk(CC_shuffled,100-which_percentile);
                    perc = perc(100-which_percentile);
                    if EBC_or_EBOC  == "EBOC"
                        rmA(firsthalf_session.occ_ns < 50) = nan;
                        rmB(lasthalf_session.occ_ns < 50) = nan;
                    end
                    corr = corrcoef(rmA, rmB, 'Rows', 'complete');
                    c2 = corr(1,2) >= perc;
                    c3 = c2;

                    cdfplot(CC_shuffled)
                    set(findall(gca, 'Type', 'Line'),'LineWidth',5);
                    hold on
                    %line([full_session.MRL,full_session.MRL], ylim, 'LineWidth', 3, 'Color', '#ff7f0e');
                    plot(corr(1,2),1,'.', 'MarkerSize',70, 'Color','#ff7f0e')
                    set(gca, 'box', 'off')
                    grid off
                    ylim([0,1.05])
                    if corr(1,2) > max(CC_shuffled)
                        xlim([min(CC_shuffled),corr(1,2)+0.04])
                    end
                    view([90 -90])   
                    set(gca,'YTick',[0,0.5,1])
                    set(gca,'XTick',[gca().XTick(1),(gca().XTick(end)+gca().XTick(1))/2,gca().XTick(end)])
                    ylabel('Cumulative dist. function')
                    xlabel('Corr.')
                    title([])
                    Stab_p_vals(i) = sum(corr(1,2) > sort(CC_shuffled)) ;

                    if EBC_or_EBOC  == "EBOC"
                        filename_pdf = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_bait', string(i), '.pdf');
                        filename_png = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_bait', string(i), '.png');
                    else
                        filename_pdf = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr', string(i), '.pdf');
                        filename_png = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr', string(i), '.png');
                    end
                    %exportgraphics(gcf,filename_pdf, 'Resolution', 600)
                    %exportgraphics(gcf, filename_png)
                    close

                else
                    corr = corrcoef(rmA, rmB, 'Rows', 'complete');
                    CC_shuffled(k) = corr(1,2);
                end
            end

        else
            PO1 = firsthalf_session.PrefOrient;
            PO2 = lasthalf_session.PrefOrient;
            d1 = abs(PO1 - PO2);
            if d1 > pi
                d1 = 2*pi - d1;
            end
    
            spmat1 = firsthalf_session.nspk_ns;
            spmatb1 = (spmat1 > prctile(firsthalf_session.nspk_ns,75, "all"));
            [rows1, cols1] = ndgrid(1:size(spmat1, 1), 1:size(spmat1, 2));
            rowcentre1 = sum(rows1(spmatb1) .* spmat1(spmatb1)) / sum(spmat1(spmatb1));
            colcentre1 = sum(cols1(spmatb1) .* spmat1(spmatb1)) / sum(spmat1(spmatb1));
            
            spmat2 = lasthalf_session.nspk_ns;
            spmatb2 = (spmat2 > prctile(lasthalf_session.nspk_ns,75, "all"));
            [rows2, cols2] = ndgrid(1:size(spmat2, 1), 1:size(spmat2, 2));
            rowcentre2 = sum(rows2(spmatb2) .* spmat2(spmatb2)) / sum(spmat2(spmatb2));
            colcentre2 = sum(cols2(spmatb2) .* spmat2(spmatb2)) / sum(spmat2(spmatb2));
    
            spmat = full_session.nspk_ns;
            spmatb = (spmat > prctile(full_session.nspk_ns,75, "all"));
            [rows, cols] = ndgrid(1:size(spmat, 1), 1:size(spmat, 2));
            rowcentre = sum(rows(spmatb) .* spmat(spmatb)) / sum(spmat(spmatb));
            colcentre = sum(cols(spmatb) .* spmat(spmatb)) / sum(spmat(spmatb));
    
            disp(["second, alternative", abs(colcentre1-colcentre2)])
            disp([colcentre, colcentre1, colcentre2])
    
            d2 = abs(colcentre1-colcentre2);
            if d2 > size(spmat,2)/2
                d2 = size(spmat,2) - d2;
            end
            
            if test_nr2 == 1
                c2 = d1 < pi/4;
            else
                c2 = (d2*360/(size(spmat,2)-1)) < 45;
            end
        end
        if c2 % passed second criterion
            disp("Neuron number "+i+" passed second criterion")
            if not(CC_stability)
            
                %d1 = abs(firsthalf_session.PrefDist - full_session.PrefDist);
                %change1 = d1/firsthalf_session.PrefDist;
                %d2 = abs(lasthalf_session.PrefDist - full_session.PrefDist);
                %change2 = d2/lasthalf_session.PrefDist;
                %c3 = change1 < 0.5 & change2 < 0.5;
                
                d1 = abs(firsthalf_session.PrefDist - lasthalf_session.PrefDist);
                disp(["first",d1 < full_session.params.distanceBins(end)*1/6])
                disp([full_session.PrefDist,firsthalf_session.PrefDist,lasthalf_session.PrefDist])         
                
                d2 = abs(rowcentre1-rowcentre2);
                disp(["second", abs(rowcentre1-rowcentre2) < full_session.params.distanceBins(end)*1/6/(full_session.params.distanceBins(2)-full_session.params.distanceBins(1))])
                disp([rowcentre,rowcentre1,rowcentre2])
                
                [M1,I1] = max(firsthalf_session.nspk_ns,[], "all");
                [M2,I2] = max(lasthalf_session.nspk_ns,[], "all");
                [M,I] = max(full_session.nspk_ns,[], "all");
                [I_row1, I_col1] = ind2sub(size(firsthalf_session.nspk_ns),I1);
                [I_row2, I_col2] = ind2sub(size(lasthalf_session.nspk_ns),I2);
                [I_row, I_col] = ind2sub(size(full_session.nspk_ns),I);
                d3 = abs(I_row1 - I_row2);
                c3 = d3 < full_session.params.distanceBins(end)*1/6/(full_session.params.distanceBins(2)-full_session.params.distanceBins(1));
                disp(["second",c3])
                disp([I_row,I_row1,I_row2])
    
                if test_nr3 == 1
                    c3 = d1 <= full_session.params.distanceBins(end)*1/6;
                end
                if test_nr3 == 2
                    c3 = d2 <= full_session.params.distanceBins(end)*1/6/(full_session.params.distanceBins(2)-full_session.params.distanceBins(1));
                else
                    c3 = d3 < full_session.params.distanceBins(end)*1/6/(full_session.params.distanceBins(2)-full_session.params.distanceBins(1));
                end
            end

            if c3 % passed third criterion, is an EBC
                if c1
                    EBC_array(i) = 1;
                    if EBC_or_EBOC == "EBOC"
                        disp("Neuron number "+i+" is an EBOC")
                    else
                        disp("Neuron number "+i+" is an EBC")
                    end
                    N_EBC = N_EBC + 1;

                    
                end
            end
        end
    end
end

if EBC_or_EBOC == "EBOC"
    disp("There are "+N_EBC+" total EBOCs out of "+N+" cells")
else
    disp("There are "+N_EBC+" total EBCs out of "+N+" cells")
end  
disp(N-cnt)