%% Read data and do the EBC test computations

% GOTTA RUN ALL THIS SHIT AGAIN BECAUSE I OVERWRITE THE SPEED IN EVERY TIME
% LIKE AN IDIOT, HERE'S WHICH HAVE BEEN RUN SO FAR: Arwen OF1&OF2, Tauriel
% OF1

% ACTUALLY GOTTA RUN IT ALL OVER AGAIN BECAUSE I HAVE USED SC DATA FOR
% ARWEN AND NOT RSC, SO ALL THE DATA THAT'S THERE IS FROM SC
% AND TAURIEL HAS A MIX OF BOTH, WHICH IS PURE MESS...
% Tauriel: OF1, OF1&c1&c2&c4&c5
% Arwen: OF1&OF2&c1&c2&c4, c1&c2&c4

% Read the .mat file with all the binned rat data
folder_loc = '/Users/martibma/Library/CloudStorage/OneDrive-NTNU/PhD/RSC_project/';
which_animal = "Arwen"; %{'Luke', 'Arwen', 'Tauriel'}
which_session = "c1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["c1", "c2", "c4"];
%concat_which = ["OF1", "OF2", "c1", "c2", "c4"];
chase_or_chill = "chase";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
EBC_or_EBOC = "EBC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait

CC_stability = true; % compute shifted ratemaps, for using cross-corrs. to determine stability, not peak in ratemap like Andy does
%%

if shared_cells
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binnedshareddata', string(binsize*1000),'ms.mat'));
end

if not(shared_cells)
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binneddata', string(binsize*1000),'ms.mat'));
end

if which_session == "c1"
    disp(which_session)
    if chase_or_chill == "chase"
        c_intervals = [1610,7583,27805,38704,53618,57788];
    end
    if chase_or_chill == "chill"
        c_intervals = [1,1610,7583,27805,38704,53618,57788,size(data.spikemat,2)];
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
    
    data_length = size(data.spikemat(:, interval_bins_c),2);
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
            time_intervals = cat(2, time_intervals, time_intervals(end):tracking_interval:(tracking_interval*size(data.spikemat(:,interval_bins_c),2)+time_intervals(end)));
        end
        % THINK ITS OKAY FROM HERE, HAVE EXTRACTED THE SAME AS UNDER SO IT
        % SHOULD WORK NOW I HOPE
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(neuron_index), 'full.mat');
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(neuron_index), 'full.mat');
    end

    save(filename, 'out');
    
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(neuron_index), 'firsthalf.mat');
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(neuron_index), 'firsthalf.mat');
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(neuron_index), 'lasthalf.mat');
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(neuron_index), 'lasthalf.mat');
    end

    save(filename, 'out');

end

%% Do the actual test, loading the data

% Read the .mat file with all the binned rat data
folder_loc = '/Users/martibma/Library/CloudStorage/OneDrive-NTNU/PhD/RSC_project/';
which_animal = "Arwen"; %{'Luke', 'Arwen', 'Tauriel'}
which_session = "OF1"; %Luke: {'session2chasing_solo'}, Arwen: {c1, c2, c4, OF1, OF2}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["OF1", "OF2", "c1", "c2", "c4"];
%concat_which = ["c1", "c2", "c4"];
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
EBC_or_EBOC = "EBC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait


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
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/bait_oriented/neuron', string(i), 'full.mat');
            %filename = "EBC_results/neuron"+i+"full.mat";
            full_session = load(filename).out;
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/bait_oriented/neuron', string(i), 'firsthalf.mat');
            %filename = "EBC_results/neuron"+i+"firsthalf.mat";
            firsthalf_session = load(filename).out;
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/bait_oriented/neuron', string(i), 'lasthalf.mat');
            %filename = "EBC_results/neuron"+i+"lasthalf.mat";
            lasthalf_session = load(filename).out;
        else
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(i), 'full.mat');
            %filename = "EBC_results/neuron"+i+"full.mat";
            full_session = load(filename).out;
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(i), 'firsthalf.mat');
            %filename = "EBC_results/neuron"+i+"firsthalf.mat";
            firsthalf_session = load(filename).out;
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/neuron', string(i), 'lasthalf.mat');
            %filename = "EBC_results/neuron"+i+"lasthalf.mat";
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
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(i), 'full.mat');
        %filename = "EBC_results/neuron"+i+"full.mat";
        full_session = load(filename).out;
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(i), 'firsthalf.mat');
        %filename = "EBC_results/neuron"+i+"firsthalf.mat";
        firsthalf_session = load(filename).out;
        filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/neuron', string(i), 'lasthalf.mat');
        %filename = "EBC_results/neuron"+i+"lasthalf.mat";
        lasthalf_session = load(filename).out;
    end
    
    if EBC_or_EBOC == "EBOC"
        perc95 = maxk(full_session.MI_dist,6);
        perc95 = perc95(6);
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
        filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MI_distr', string(i), '.pdf');
        filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MI_distr', string(i), '.png');
    else
        xlabel('MRL')
        filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MRL_distr', string(i), '.pdf');
        filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MRL_distr', string(i), '.png');
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
                %{
                offs_dist = (-size(rmA,1)+3):(size(rmA,1)-3);
                if mod(size(rmA,2),2) == 1
                    offs_angle = (-floor(size(rmA,2)/2)):(floor(size(rmA,2)/2));
                else
                    offs_angle = (-size(rmA,2)/2):(size(rmA,2)/2-1);
                end
                cc_plot = NaN(length(offs_dist), length(offs_angle));
                
                for ii=1:length(offs_angle)
                    rotB = circshift(rmB,offs_angle(ii),2);
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
                            cc_plot(j,ii) = corr(1,2);
                        end
                    end
                end
                if k == 101
                    dist_perc = maxk(CC_stability_mat(1,:),which_percentile);
                    ang_perc = maxk(CC_stability_mat(2,:),which_percentile);
                    dist_perc = dist_perc(which_percentile);
                    ang_perc = ang_perc(which_percentile);

                    %[~, ind] = max(cc_plot,[],'all');
                    %[row_max,col_max] = ind2sub(size(cc_plot),ind);
                    
                    x = 1 : size(cc_plot, 2); % Columns.
                    y = 1 : size(cc_plot, 1); % Rows.
                    [X, Y] = meshgrid(x, y);
                    meanCC = mean(cc_plot(:));
                    col_max = mean(cc_plot(:) .* X(:)) / meanCC;
                    row_max = mean(cc_plot(:) .* Y(:)) / meanCC;
                    %disp([col_max, row_max])
                    surface(X,Y,cc_plot)
                    set(gca, 'Ydir', 'reverse')

                    %spmat = cc_plot;
                    %spmatb = (spmat > prctile(cc_plot,75, "all"));
                    %[rows, cols] = ndgrid(1:size(spmat, 1), 1:size(spmat, 2));
                    %row_max = sum(rows(spmatb) .* spmat(spmatb)) / sum(spmat(spmatb));
                    %col_max = sum(cols(spmatb) .* spmat(spmatb)) / sum(spmat(spmatb));
                    %disp([col_max, row_max])

                    centre_dist = (ceil(length(offs_dist)/2));
                    centre_ang = (floor(length(offs_angle)/2)+1);
                    dist_dist = abs(diff([row_max, centre_dist]));
                    ang_dist = abs(diff([col_max, centre_ang]));
                    %disp([dist_perc, dist_dist])
                    %disp([ang_perc, ang_dist])
                    c2 = dist_dist <= dist_perc;
                    c3 = ang_dist <= ang_perc;
                else
                    [~, ind] = max(cc_plot,[],'all');
                    [row_max,col_max] = ind2sub(size(cc_plot),ind);
                    
                    %x = 1 : size(cc_plot, 2); % Columns.
                    %y = 1 : size(cc_plot, 1); % Rows.
                    %[X, Y] = meshgrid(x, y);
                    %meanCC = mean(cc_plot(:));
                    %col_max = mean(cc_plot(:) .* X(:)) / meanCC;
                    %row_max = mean(cc_plot(:) .* Y(:)) / meanCC; 

                    %spmat = cc_plot;
                    %spmatb = (spmat > prctile(cc_plot,75, "all"));
                    %[rows, cols] = ndgrid(1:size(spmat, 1), 1:size(spmat, 2));
                    %row_max = sum(rows(spmatb) .* spmat(spmatb)) / sum(spmat(spmatb));
                    %col_max = sum(cols(spmatb) .* spmat(spmatb)) / sum(spmat(spmatb));
                    %disp([col_max, row_max])
                    
                    centre_dist = ceil(length(offs_dist)/2);
                    centre_ang = (floor(length(offs_angle)/2)+1);
                    dist_dist = abs(diff([row_max, centre_dist]));
                    ang_dist = abs(diff([col_max, centre_ang]));
                    CC_stability_mat(:,k) = [dist_dist, ang_dist];
                end
                %}
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
                        filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_bait', string(i), '.pdf');
                        filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_bait', string(i), '.png');
                    else
                        filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr', string(i), '.pdf');
                        filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr', string(i), '.png');
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
%%
figure(1)
h = histogram(1-tuned_p_vals/100, 20, 'FaceColor', "#1f77b4", 'FaceAlpha', 1, 'EdgeColor', 'none');
hold on
xlabel('"p-value" estimate')
ylabel('frequency')
set(gca, 'box', 'off')
set(gca,'XTick',[0,1],'YTick',[max(h.Values)])
histogram(ones(h.Values(1),1)*0.025, 'BinEdges',[0,0.05], 'FaceColor', '#ff7f0e', 'FaceAlpha', 1, 'EdgeColor', 'none')

if EBC_or_EBOC == "EBOC"
    filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MI_distr_p_vals.pdf');
    filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MI_distr_p_vals.png');
else
    filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MRL_distr_p_vals.pdf');
    filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/MRL_distr_p_vals.png');
end
exportgraphics(gcf,filename_pdf, 'Resolution', 600)
exportgraphics(gcf, filename_png)
%close

figure(2)
h = histogram(1-Stab_p_vals/100, 20, 'FaceColor', "#1f77b4", 'FaceAlpha', 1, 'EdgeColor', 'none');
hold on
xlabel('"p-value" estimate')
ylabel('frequency')
set(gca, 'box', 'off')
set(gca,'XTick',[0,1],'YTick',[max(h.Values)])
histogram(ones(h.Values(1),1)*0.025, 'BinEdges',[0,0.05], 'FaceColor', '#ff7f0e', 'FaceAlpha', 1, 'EdgeColor', 'none')

if EBC_or_EBOC == "EBOC"
    filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_p_vals_bait.pdf');
    filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_p_vals_bait.png');
else
    filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_p_vals.pdf');
    filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_distr_p_vals.png');
end

exportgraphics(gcf,filename_pdf, 'Resolution', 600)
exportgraphics(gcf, filename_png)
%close

figure(3)
isEBC = EBC_array == 1;
scatter(1-tuned_p_vals(~isEBC)/100+0.01, 1-Stab_p_vals(~isEBC)/100+0.01, '.', 'SizeData', 1600, 'Color', "#1f77b4")
hold on
%plot(full_session.MRL,1,'.', 'MarkerSize',40, 'Color','#ff7f0e')
scatter(1-tuned_p_vals(isEBC)/100+0.01, 1-Stab_p_vals(isEBC)/100+0.01, '.', 'SizeData', 1600, 'Color','#ff7f0e')
xlim([-0.05,1])
ylim([-0.05,1])
if EBC_or_EBOC == "EBOC"
    xlabel('"p-values" MI')
else
    xlabel('"p-values" MRL')
end
ylabel('"p-values" stability')
set(gca,'XTick',[0,0.5,1],'YTick',[0,0.5,1])
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%xlim([0.009,1])
%ylim([0.009,1])
%set(gca, 'YDir','reverse')
%set(gca, 'XDir','reverse')


figure(3)
new_idea_x = 1-tuned_p_vals/100+0.01;
new_idea_y = 1-Stab_p_vals/100+0.01;
%outer_edge_x = new_idea_x > 0.3;
%outer_edge_y = new_idea_y > 0.3;
%new_idea_x(outer_edge_x) = 0.35;
%new_idea_y(outer_edge_y) = 0.35;

%hist3([new_idea_x, new_idea_y], 'NBins', [20 20], 'CdataMode','auto')
%view(2)
histogram2(new_idea_x,new_idea_y,[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
colorbar

if EBC_or_EBOC == "EBOC"
    filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_MI_heatmap.pdf');
    filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_MI_heatmap.png');
else
    filename_pdf = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_MRL_heatmap.pdf');
    filename_png = strcat('EBC_results/', which_animal, '/', which_channels, '/All/Stab_MRL_heatmap.png');
end
exportgraphics(gcf,filename_pdf, 'Resolution', 600)
exportgraphics(gcf, filename_png)

for elem=1:N
    if EBC_array(elem) == 1
        disp(elem)
        disp(data.cell_names(elem,:))
    end
end

%% Plot neurons that are classified as EBC

% Read the .mat file with all the binned rat data
folder_loc = '/Users/martibma/Library/CloudStorage/OneDrive-NTNU/PhD/RSC_project/';
which_animal = "Arwen"; %{'Luke', 'Arwen', 'Tauriel'}
which_session = "OF1"; %Luke: {'session2chasing_solo'}, Arwen: {c1, c2, c4, OF1, OF2}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["OF1", "OF2", "c1", "c2", "c4"];
chase_or_chill = "chill";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
odd_even_maps = false;
EBC_or_EBOC = "EBC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait

folder_loc = '/Users/martibma/Library/CloudStorage/OneDrive-NTNU/PhD/RSC_project/';
which_animal = "Arwen"; %{'Luke', 'Arwen', 'Tauriel'}
which_session = "c1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["c1", "c2", "c4"];
chase_or_chill = "chase";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
odd_even_maps = true;
EBC_or_EBOC = "EBOC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait


if shared_cells
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binnedshareddata', string(binsize*1000),'ms.mat'));
end

if not(shared_cells)
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binneddata', string(binsize*1000),'ms.mat'));
end

if which_session == "c1"
    disp(which_session)
    if chase_or_chill == "chase"
        c_intervals = [1610,7583,27805,38704,53618,57788];
    end
    if chase_or_chill == "chill"
        c_intervals = [1,1610,7583,27805,38704,53618,57788,size(data.spikemat,2)];
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
    
    tracking_interval = 0.008332827537658849;
    data_length = size(data.spikemat(:, interval_bins_c),2);
    time_intervals = 0.0042:tracking_interval:(tracking_interval*data_length);
end
if which_session ~= "c1"
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
        
            x_pos = cat(1, x_pos, data.binned_pos(interval_bins_c,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(interval_bins_c,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd(interval_bins_c));
            speed = cat(2, speed, data.binned_speed(interval_bins_c));
            spikemat = cat(2, spikemat, data.spikemat(:,interval_bins_c));
            if EBC_or_EBOC == "EBOC"
                bait_angle = cat(2, bait_angle, data.binned_rel_ha(interval_bins_c));
                bait_dist = cat(2, bait_dist, data.binned_rel_dist(interval_bins_c));
            end

            session = cat(1, session, repmat(concat_which(i), size(data.binned_pos,1), 1));
            behav = cat(1,behav, repmat(chase_or_chill, size(data.binned_pos,1), 1));
            bin_number = cat(2, bin_number, interval_bins_c);

            data_length = data_length + size(data.spikemat(:, interval_bins_c),2);
            time_intervals = cat(2, time_intervals, time_intervals(end):tracking_interval:(tracking_interval*size(data.spikemat(:,interval_bins_c),2)+time_intervals(end)));
        end
        % THINK ITS OKAY FROM HERE, HAVE EXTRACTED THE SAME AS UNDER SO IT
        % SHOULD WORK NOW I HOPE
        if concat_which(i) == "OF1" | concat_which(i) == "OF2"
            disp(concat_which(i))
            x_pos = cat(1, x_pos, data.binned_pos(:,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(:,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd);
            speed = cat(2, speed, data.binned_speed);
            spikemat = cat(2, spikemat, data.spikemat);

            session = cat(1, session, repmat(concat_which(i), size(data.binned_pos,1), 1));
            behav = cat(1, behav, repmat("OF", size(data.binned_pos,1), 1));
            bin_number = cat(2, bin_number, 1:size(data.binned_pos,1));
            
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

session = session(ind);
behav = behav(ind);
bin_number = bin_number(ind);

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

    session = session(ind);
    behav = behav(ind);
    bin_number = bin_number(ind);
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

    session = session(ind);
    behav = behav(ind);
    bin_number = bin_number(ind);

    ind = ~isnan(bait_dist);
    bait_angle = bait_angle(ind);
    bait_dist = bait_dist(ind);
    x_pos = x_pos(ind);
    y_pos = y_pos(ind);
    head_dir = head_dir(ind);
    spikemat = spikemat(:,ind);
    time_intervals = time_intervals(ind);
    speed = speed(ind);

    session = session(ind);
    behav = behav(ind);
    bin_number = bin_number(ind);

end

Conc_filtered = struct;
Conc_filtered.x = x_pos;
Conc_filtered.y = y_pos;
Conc_filtered.hd = head_dir;
Conc_filtered.spike = spikemat;
Conc_filtered.speed = speed;
Conc_filtered.bait_angle = bait_angle.';
Conc_filtered.bait_dist = bait_dist.';
Conc_filtered.session = char(session);
Conc_filtered.behav = char(behav);
Conc_filtered.bin_number = bin_number.';

save("Chasing_conc_filtered_new.mat","-struct","Conc_filtered")
%save("OF_conc_filtered_new.mat","-struct","Conc_filtered")

N = size(data.spikemat,1);
%EBC_array = ones(length(EBC_array)); % just set them all to one if you wanna plot all neurons
% Run test computations and save stuff for each neuron
observed_mean_spikes_in_bins = zeros(length(EBC_array), 30);
rr = 1;
cc = 1;
nn = 1;
for i=40
    if true
    %if EBC_array(i) == 1
        disp("iteration number "+i)
        neuron_index = i;
        spike_vec = spikemat(neuron_index,:);

        r = struct;
        r.x = x_pos;
        r.y = y_pos;
        r.md = head_dir.';
        r.spike = spike_vec.';
        r.ts = time_intervals.';
        r.mrl = 0;
        r.mi = 0;
        if EBC_or_EBOC == "EBOC"
            r.bait_angle = bait_angle.';
            r.bait_dist = bait_dist.';
        end
        %{
        bins = linspace(-pi-0.000001, pi+0.0000001, 30 + 1);
        x_grid = 0.5*(bins(1:end-1)+bins(2:end));
        
        for xx=1:30
            timesinbin = (head_dir>bins(xx)).*(head_dir<bins(xx+1));
            if(sum(timesinbin)>0)
                observed_mean_spikes_in_bins(i,xx) = mean(spike_vec(timesinbin>0));
                
            end
        end
        subplot(2,6,nn)
        plot(x_grid, observed_mean_spikes_in_bins(i,:)/0.008)
        nn = nn+1;
        title(i)
        %}
        
        cfull = EgocentricRatemap(r, which_animal);
        %plotEBC(r, cfull, i)

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
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBOC_ratemap_neuron', string(i), '.png');
            else
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBC_ratemap_neuron', string(i), '.png');
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
            if EBC_or_EBOC == "EBOC"
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBOC_ratemap_neuron', string(i), '.png');
            else
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBC_ratemap_neuron', string(i), '.png');
            end
        end
        %saveas(gcf, filename)
        %close
        if odd_even_maps
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
                
            for t=1:chunks-1
                if mod(t,2) == 0
                    r1.x = cat(1, r1.x, x_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
                    r1.y = cat(1, r1.y, y_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
                    r1.md = cat(1, r1.md, head_dir(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r1.spike = cat(1, r1.spike, spike_vec(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r1.ts = cat(1, r1.ts, time_intervals(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r1.bait_angle = cat(1, r1.bait_angle, bait_angle(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r1.bait_dist = cat(1, r1.bait_dist, bait_dist(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                end
                if (mod(t,2)) == 1
                    r2.x = cat(1, r2.x, x_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
                    r2.y = cat(1, r2.y, y_pos(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
                    r2.md = cat(1, r2.md, head_dir(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r2.spike = cat(1, r2.spike, spike_vec(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r2.ts = cat(1, r2.ts, time_intervals(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r2.bait_angle = cat(1, r2.bait_angle, bait_angle(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                    r2.bait_dist = cat(1, r2.bait_dist, bait_dist(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
                end
            end
        
            if (mod(t+1,2)) == 0
                r1.x = cat(1, r1.x, x_pos((t+1)*bins_in_chunk:end));
                r1.y = cat(1, r1.y, y_pos((t+1)*bins_in_chunk:end));
                r1.md = cat(1, r1.md, head_dir((t+1)*bins_in_chunk:end).');
                r1.spike = cat(1, r1.spike, spike_vec((t+1)*bins_in_chunk:end).');
                r1.ts = cat(1, r1.ts, time_intervals((t+1)*bins_in_chunk:end).');
                r1.bait_angle = cat(1, r1.bait_angle, bait_angle((t+1)*bins_in_chunk:end).');
                r1.bait_dist = cat(1, r1.bait_dist, bait_dist((t+1)*bins_in_chunk:end).');
            end
            if (mod(t+1,2)) == 1
                r2.x = cat(1, r2.x, x_pos((t+1)*bins_in_chunk:end));
                r2.y = cat(1, r2.y, y_pos((t+1)*bins_in_chunk:end));
                r2.md = cat(1, r2.md, head_dir((t+1)*bins_in_chunk:end).');
                r2.spike = cat(1, r2.spike, spike_vec((t+1)*bins_in_chunk:end).');
                r2.ts = cat(1, r2.ts, time_intervals((t+1)*bins_in_chunk:end).');
                r2.bait_angle = cat(1, r2.bait_angle, bait_angle((t+1)*bins_in_chunk:end).');
                r2.bait_dist = cat(1, r2.bait_dist, bait_dist((t+1)*bins_in_chunk:end).');
            end
    
            codd = EgocentricRatemap(r1, which_animal);
            %plotEBC(r1, codd, i)    
    
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
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBC_ratemap_neuron', string(i), 'even.png');
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
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBC_ratemap_neuron', string(i), 'even.png');
            end
            
            %saveas(gcf, filename)
            %close
    
            ceven = EgocentricRatemap(r2, which_animal);
            %plotEBC(r2, ceven, i)    
    
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
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBC_ratemap_neuron', string(i), 'odd.png');
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
                filename = strcat('EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBC_ratemap_neuron', string(i), 'odd.png');
            end
            
            %saveas(gcf, filename)
            %close

            root = struct;
            root.cfull = r;
            root.codd = r1;
            root.ceven = r2;

            out = struct;
            out.cfull = cfull;
            out.codd = codd;
            out.ceven = ceven;
            throwaway = plot_all_EBC(root, out, 1, 0, 1, 1);
            %plotEBC(r,cfull,2)
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/All/EBOC_ratemap&CC_neuron', string(i), '_polar.pdf');
            %exportgraphics(gcf,filename, 'Resolution', 600) 
            %saveas(gcf, filename)
            filename = strcat('EBC_results/', which_animal, '/', which_channels, '/All/EBOC_ratemap&CC_neuron', string(i), '_polar.png');
            %exportgraphics(gcf,filename) 
            %saveas(gcf, filename)
            close
        end
    end
    
end

%hist3([bait_angle', bait_dist'], 'NBins', [36 28], 'CdataMode','auto')
%colorbar
%view(2)
%% Plot all 6 rate maps for all neurons (chase, of, even, odd) (and cross corr as well, if you want, and spike maps)

% Read the .mat file with all the binned rat data
folder_loc = '/Users/martibma/Library/CloudStorage/OneDrive-NTNU/PhD/RSC_project/';
which_animal = "Arwen"; %{'Luke', 'Arwen', 'Tauriel'}
which_session = "OF1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["OF1", "OF2", "c1", "c2", "c4"];
chase_or_chill = "chill";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's, or filter specific parts of sessions (chasing bouts)
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
EBC_or_EBOC = "EBC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait

if shared_cells
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binnedshareddata', string(binsize*1000),'ms.mat'));
end

if not(shared_cells)
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binneddata', string(binsize*1000),'ms.mat'));
end

x_pos = data.binned_pos(:,1) * 100;
y_pos = data.binned_pos(:,2) * 100;
head_dir = data.binned_hd;
speed = data.binned_speed;
spikemat = data.spikemat;

session = repmat(which_session, size(data.binned_pos,1), 1);
behav = NaN(size(data.binned_pos,1), 1);
bin_number = 1:size(data.binned_pos,1);

data_length = size(data.spikemat,2);
tracking_interval = 0.008332827537658849;
time_intervals = 0.0042:tracking_interval:tracking_interval*data_length;
N = size(data.spikemat,1);

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

            x_pos = cat(1, x_pos, data.binned_pos(interval_bins_c,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(interval_bins_c,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd(interval_bins_c));
            speed = cat(2, speed, data.binned_speed(interval_bins_c));
            spikemat = cat(2, spikemat, data.spikemat(:,interval_bins_c));
            session = cat(1, session, repmat(concat_which(i), size(data.binned_pos,1), 1));
            behav = cat(1,behav, repmat(chase_or_chill, size(data.binned_pos,1), 1));
            bin_number = cat(2, bin_number, interval_bins_c);

            data_length = data_length + size(data.spikemat(:, interval_bins_c),2);
            time_intervals = cat(2, time_intervals, time_intervals(end):tracking_interval:(tracking_interval*size(data.spikemat(:,interval_bins_c),2)+time_intervals(end)));
        end
        % THINK ITS OKAY FROM HERE, HAVE EXTRACTED THE SAME AS UNDER SO IT
        % SHOULD WORK NOW I HOPE
        if concat_which(i) == "OF1" | concat_which(i) == "OF2"
            disp(concat_which(i))
            x_pos = cat(1, x_pos, data.binned_pos(:,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(:,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd);
            speed = cat(2, speed, data.binned_speed);
            spikemat = cat(2, spikemat, data.spikemat);
            session = cat(1, session, repmat(concat_which(i), size(data.binned_pos,1), 1));
            behav = cat(1, behav, NaN(size(data.binned_pos,1), 1));
            bin_number = cat(2, bin_number, 1:size(data.binned_pos,1));
            
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

session = session(ind);
behav = behav(ind);
bin_number = bin_number(ind);

% Filter out where animal is standing still
if speed_threshold
    ind = speed > 5; % speed thresholding by 5cm/s
    x_pos = x_pos(ind);
    y_pos = y_pos(ind);
    head_dir = head_dir(ind);
    spikemat = spikemat(:,ind);
    time_intervals = time_intervals(ind);
    speed = speed(ind);

    session = session(ind);
    behav = behav(ind);
    bin_number = bin_number(ind);
end

x_pos_OF = x_pos;
y_pos_OF = y_pos;
head_dir_OF = head_dir;
spikemat_OF = spikemat;
time_intervals_OF = time_intervals;
speed_OF = speed;

%OF_conc_filtered = struct;
%OF_conc_filtered.x = x_pos;
%OF_conc_filtered.y = y_pos;
%OF_conc_filtered.hd = head_dir;
%OF_conc_filtered.spike = spikemat;
%OF_conc_filtered.session = session;
%OF_conc_filtered.behav = behav;
%OF_conc_filtered.bin_number = bin_number.';

%save('OF_conc_filtered.mat','OF_conc_filtered')

% Read the .mat file with all the binned rat data
folder_loc = '/Users/martibma/Library/CloudStorage/OneDrive-NTNU/PhD/RSC_project/';
which_animal = "Arwen"; %{'Luke', 'Arwen', 'Tauriel'}
which_session = "c1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["c1", "c2", "c4"];
chase_or_chill = "chase";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}

if shared_cells
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binnedshareddata', string(binsize*1000),'ms.mat'));
end

if not(shared_cells)
    data = load(strcat(folder_loc,'Data/', which_animal,'/',which_channels,'_',which_session,'_binneddata', string(binsize*1000),'ms.mat'));
end

if which_session == "c1"
    disp(which_session)
    if chase_or_chill == "chase"
        c_intervals = [1610,7583,27805,38704,53618,57788];
    end
    if chase_or_chill == "chill"
        c_intervals = [1,1610,7583,27805,38704,53618,57788,size(data.spikemat,2)];
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
    
    data_length = size(data.spikemat(:, interval_bins_c),2);
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
    behav = NaN(size(data.binned_pos,1), 1);
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

            x_pos = cat(1, x_pos, data.binned_pos(interval_bins_c,1) * 100);
            y_pos = cat(1, y_pos, data.binned_pos(interval_bins_c,2) * 100);
            head_dir = cat(2, head_dir, data.binned_hd(interval_bins_c));
            speed = cat(2, speed, data.binned_speed(interval_bins_c));
            spikemat = cat(2, spikemat, data.spikemat(:,interval_bins_c));

            data_length = data_length + size(data.spikemat(:, interval_bins_c),2);
            time_intervals = cat(2, time_intervals, time_intervals(end):tracking_interval:(tracking_interval*size(data.spikemat(:,interval_bins_c),2)+time_intervals(end)));
        end
        % THINK ITS OKAY FROM HERE, HAVE EXTRACTED THE SAME AS UNDER SO IT
        % SHOULD WORK NOW I HOPE
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

% Filter out where animal is standing still
if speed_threshold
    ind = speed > 5; % speed thresholding by 5cm/s
    x_pos = x_pos(ind);
    y_pos = y_pos(ind);
    head_dir = head_dir(ind);
    spikemat = spikemat(:,ind);
    time_intervals = time_intervals(ind);
    speed = speed(ind);
end
 
x_pos_c = x_pos;
y_pos_c = y_pos;
head_dir_c = head_dir;
spikemat_c = spikemat;
time_intervals_c = time_intervals;
speed_c = speed;


N = size(data.spikemat,1);
CC_shift_mat = zeros(2,N);
for i=1:N
    disp("iteration number "+i)
    neuron_index = i;
    spike_vec_OF = spikemat_OF(neuron_index,:);
    root = struct;

    r = struct;
    r.x = x_pos_OF;
    r.y = y_pos_OF;
    r.md = head_dir_OF.';
    r.spike = spike_vec_OF.';
    r.ts = time_intervals_OF.';
    r.mrl = 0;
    r.mi = 0;

    OFfull = EgocentricRatemap(r, which_animal);
    root.OFfull = r;

    %This part is the even odd 30 sek split, better use this than half the
    %data set, since we have so little
    bins_in_chunk = 3750; %30sek chunks
    chunks = floor(size(x_pos_OF,1)/bins_in_chunk);
    r1 = struct;
    r1.x = [];
    r1.y = [];
    r1.md = [];
    r1.spike = [];
    r1.ts = [];
    r1.mrl = 0;
    r1.mi = 0;
    r2 = struct;
    r2.x = [];
    r2.y = [];
    r2.md = [];
    r2.spike = [];
    r2.ts = [];
    r2.mrl = 0;
    r2.mi = 0;

    for t=1:chunks-1
        if mod(t,2) == 0
            r1.x = cat(1, r1.x, x_pos_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r1.y = cat(1, r1.y, y_pos_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r1.md = cat(1, r1.md, head_dir_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.spike = cat(1, r1.spike, spike_vec_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.ts = cat(1, r1.ts, time_intervals_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
        end
        if (mod(t,2)) == 1
            r2.x = cat(1, r2.x, x_pos_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r2.y = cat(1, r2.y, y_pos_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r2.md = cat(1, r2.md, head_dir_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.spike = cat(1, r2.spike, spike_vec_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.ts = cat(1, r2.ts, time_intervals_OF(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
        end
    end

    if (mod(t+1,2)) == 0
        r1.x = cat(1, r1.x, x_pos_OF((t+1)*bins_in_chunk:end));
        r1.y = cat(1, r1.y, y_pos_OF((t+1)*bins_in_chunk:end));
        r1.md = cat(1, r1.md, head_dir_OF((t+1)*bins_in_chunk:end).');
        r1.spike = cat(1, r1.spike, spike_vec_OF((t+1)*bins_in_chunk:end).');
        r1.ts = cat(1, r1.ts, time_intervals_OF((t+1)*bins_in_chunk:end).');
    end
    if (mod(t+1,2)) == 1
        r2.x = cat(1, r2.x, x_pos_OF((t+1)*bins_in_chunk:end));
        r2.y = cat(1, r2.y, y_pos_OF((t+1)*bins_in_chunk:end));
        r2.md = cat(1, r2.md, head_dir_OF((t+1)*bins_in_chunk:end).');
        r2.spike = cat(1, r2.spike, spike_vec_OF((t+1)*bins_in_chunk:end).');
        r2.ts = cat(1, r2.ts, time_intervals_OF((t+1)*bins_in_chunk:end).');
    end

    OFeven = EgocentricRatemap(r1, which_animal);  
    OFodd = EgocentricRatemap(r2, which_animal);
    root.OFeven = r1;
    root.OFodd = r2;

    spike_vec_c = spikemat_c(neuron_index,:);

    r = struct;
    r.x = x_pos_c;
    r.y = y_pos_c;
    r.md = head_dir_c.';
    r.spike = spike_vec_c.';
    r.ts = time_intervals_c.';
    r.mrl = 0;
    r.mi = 0;

    cfull = EgocentricRatemap(r, which_animal);
    root.cfull = r;

    %This part is the even odd 30 sek split, better use this than half the
    %data set, since we have so little
    bins_in_chunk = 3750; %30sek chunks
    chunks = floor(size(x_pos_c,1)/bins_in_chunk);
    r1 = struct;
    r1.x = [];
    r1.y = [];
    r1.md = [];
    r1.spike = [];
    r1.ts = [];
    r1.mrl = 0;
    r1.mi = 0;
    r2 = struct;
    r2.x = [];
    r2.y = [];
    r2.md = [];
    r2.spike = [];
    r2.ts = [];
    r2.mrl = 0;
    r2.mi = 0;

    for t=1:chunks-1
        if mod(t,2) == 0
            r1.x = cat(1, r1.x, x_pos_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r1.y = cat(1, r1.y, y_pos_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r1.md = cat(1, r1.md, head_dir_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.spike = cat(1, r1.spike, spike_vec_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r1.ts = cat(1, r1.ts, time_intervals_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
        end
        if (mod(t,2)) == 1
            r2.x = cat(1, r2.x, x_pos_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r2.y = cat(1, r2.y, y_pos_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1));
            r2.md = cat(1, r2.md, head_dir_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.spike = cat(1, r2.spike, spike_vec_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
            r2.ts = cat(1, r2.ts, time_intervals_c(t*bins_in_chunk:(t+1)*bins_in_chunk-1).');
        end
    end

    if (mod(t+1,2)) == 0
        r1.x = cat(1, r1.x, x_pos_c((t+1)*bins_in_chunk:end));
        r1.y = cat(1, r1.y, y_pos_c((t+1)*bins_in_chunk:end));
        r1.md = cat(1, r1.md, head_dir_c((t+1)*bins_in_chunk:end).');
        r1.spike = cat(1, r1.spike, spike_vec_c((t+1)*bins_in_chunk:end).');
        r1.ts = cat(1, r1.ts, time_intervals_c((t+1)*bins_in_chunk:end).');
    end
    if (mod(t+1,2)) == 1
        r2.x = cat(1, r2.x, x_pos_c((t+1)*bins_in_chunk:end));
        r2.y = cat(1, r2.y, y_pos_c((t+1)*bins_in_chunk:end));
        r2.md = cat(1, r2.md, head_dir_c((t+1)*bins_in_chunk:end).');
        r2.spike = cat(1, r2.spike, spike_vec_c((t+1)*bins_in_chunk:end).');
        r2.ts = cat(1, r2.ts, time_intervals_c((t+1)*bins_in_chunk:end).');
    end

    ceven = EgocentricRatemap(r1, which_animal);  
    codd = EgocentricRatemap(r2, which_animal);
    root.ceven = r1;
    root.codd = r2;
   
    out = struct;
    out.OFfull = OFfull;
    out.OFodd = OFodd;
    out.OFeven = OFeven;
    out.cfull = cfull;
    out.codd = codd;
    out.ceven = ceven;
    CC_shift = plot_all_EBC(root, out, 1, 0, 1, 0);
    CC_shift_mat(:,i) = CC_shift;   
    filename = strcat('EBC_results/', which_animal, '/', which_channels, '/All/EBC_ratemap&CC_neuron', string(i), '_polar.pdf');
    %exportgraphics(gcf,filename, 'Resolution', 600) 
    %saveas(gcf, filename)
    filename = strcat('EBC_results/', which_animal, '/', which_channels, '/All/EBC_ratemap&CC_neuron', string(i), '_polar.png');
    %exportgraphics(gcf,filename) 
    %saveas(gcf, filename)
    
    close
end

%save('OF_conc_filtered.mat','root.OFfull')

%% Shift in CC center
isEBC_of = zeros([N,1]);
isEBC_of([22 25 27 28 29 30 33 36 37 39 48 56]) = 1;
isEBC_ch = zeros([N,1]);
isEBC_ch([5 13 22 24 25 28 29 32 33 35 37 39 41 48 53]) = 1;
sz = 250;
%c1 = "#1F77B4";
c1 = [0.12156862745098039 0.4666666666666667 0.7058823529411765];
%c2 = "#FF7F0E";
c2 = [1 0.4980392156862745 0.054901960784313725];
figure(1)
scatter(CC_shift_mat(1,~isEBC_of), CC_shift_mat(2,~isEBC_of), sz, c1, '+', 'Linewidth', 2)
hold on
scatter(CC_shift_mat(1,~~isEBC_of), CC_shift_mat(2,~~isEBC_of), sz, c2, '+', 'Linewidth', 2)
scatter(CC_shift_mat(1,~isEBC_ch), CC_shift_mat(2,~isEBC_ch), sz, c1, 'x', 'Linewidth', 2)
scatter(CC_shift_mat(1,~~isEBC_ch), CC_shift_mat(2,~~isEBC_ch), sz, c2, 'x', 'Linewidth', 2)
%scatter(CC_shift_mat(1,:),CC_shift_mat(2,:), '.', 'SizeData', 400)
xline(0, '--', 'LineWidth', 2.0)
yline(0, '--', 'LineWidth', 2.0)
xlabel('Shift in angle bin')
ylabel('Shift in distance bin')
xlim([-18,18])
ylim([-12,12])
set(gca,'XTick',[-18,0,18],'YTick',[-12,0,12])
exportgraphics(gcf,strcat('EBC_results/', which_animal, '/', which_channels, '/All/shift_from_OF_to_ch_both.pdf'), 'Resolution', 600) 
exportgraphics(gcf,strcat('EBC_results/', which_animal, '/', which_channels, '/All/shift_from_OF_to_ch_both.png'))

figure(2)
histogram2(CC_shift_mat(1,:),CC_shift_mat(2,:), 6,'FaceColor','flat','ShowEmptyBins','on');
colormap(parula)
xlabel('Shift in angle bin')
ylabel('Shift in distance bin')

%colorbar
%view(2)

%% Wall or center tuning curves looks better?

which_neuron = 39;

spikes = spikemat_OF(which_neuron,:);
spikes = spikes.';
x = x_pos_OF;
y = y_pos_OF;
hd = head_dir_OF.';
x_max = max(x);
x_min = min(x);
y_max = max(y);
y_min = min(y);

l_mat = zeros(length(x), 4);
l_mat(:, 1) = y_max - y;     % -1, north side according to pos
l_mat(:, 2) = x_max - x;     % 0,  east side according to pos 
l_mat(:, 3) = -(y_min - y);  % 1,  south side according to pos
l_mat(:, 4) = -(x_min - x);  % 2,  west side according to pos

[l, which_side] = min(l_mat,[],2);
which_side = which_side - 2;

phi = hd - pi/2*which_side;
%phi = phi/pi*180;
%phi = abs(phi);
phi(phi > pi) = phi(phi > pi) - pi;
phi(phi < -pi) = phi(phi < -pi) + pi;

same_color = zeros(2,2);

figure(1)
subplot(1,2,1)
%histogram2(l,phi);
histogram2(l,phi,[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
same_color(1,:) = clim;
%colorbar


phi_bins = linspace(-pi-0.000001, pi+0.0000001, 20 + 1);
l_bins = linspace(0.000001, 75+0.0000001, 20 + 1);
phi_grid = 0.5*(phi_bins(1:end-1)+phi_bins(2:end));
l_grid = 0.5*(l_bins(1:end-1)+l_bins(2:end));
observed_mean_spikes_in_bins = zeros(20, 20);

for xx=1:20
    for yy=1:20
        timesinbin = (phi>phi_bins(xx)).*(phi<phi_bins(xx+1)).*(l>l_bins(yy)).*(l<l_bins(yy+1));
        if(sum(timesinbin)>50)
            observed_mean_spikes_in_bins(xx,yy) = mean(spikes(timesinbin>0));
        elseif(sum(timesinbin)<=50)
            observed_mean_spikes_in_bins(xx,yy) = nan;
        end

    end
end
    
figure(2)
subplot(1,2,1)
surface(l_grid, phi_grid, observed_mean_spikes_in_bins)

test_root1 = root.OFfull;
test_root1.bait_angle = phi;
test_root1.bait_dist = l;
%test1 = EgocentricRatemap(test_root1, which_animal);
%plotEBC(test_root1, test1,2)


x_center = (x_max + x_min)/2;
centered_x_max = x_max - x_center;
centered_x_min = x_min - x_center;
y_center = (y_max + y_min)/2;
centered_y_max = y_max - y_center;
centered_y_min = y_min - y_center;
l = sqrt((0 - (x-x_center)).^2 + (0 - (y-y_center)).^2);
x1 = cos(hd);
y1 = cos(hd + pi/2);
x2 = 0 - (x-x_center);
y2 = 0 - (y-y_center);
%phi = acos((x1.*x2 + y1.*y2)./l);

dot = x1.*x2 + y1.*y2;     
det = x1.*y2 - y1.*x2;     
phi = atan2(det, dot);

% Testing stuff, seems to work correct now
%testhd = pi/2;
%testx = -70;
%testy = 1;
%testx1 = cos(testhd);
%testy1 = cos(testhd - pi/2);
%testx2 = -testx - 0;
%testy2 = -testy - 0;
%testdot = testx1*testx2 + testy1*testy2;
%testdet = testx1*testy2 - testy1*testx2;
%testangle = atan2(testdet, testdot)*180/pi

figure(1)
subplot(1,2,2)
h=histogram2(l,phi,[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
X_edges = h.XBinEdges;
Y_edges = h.YBinEdges;
%test_values = histogram2(l,phi,[20 20],'XBinLimits',[0, 100],'YBinLimits',[-pi, pi],'DisplayStyle','tile','ShowEmptyBins','on').Values;
same_color(2,:) = clim;

for whatever = 1:2
  subplot (1,2,whatever)
  clim([min(same_color(:,1)),max(same_color(:,2))]) % From the smallest low to the largest high 
end
colorbar

%figure(3)
%histogram2(l,phi)


phi_bins = linspace(-pi-0.000001, pi+0.0000001, 20 + 1);
l_bins = linspace(0.000001, 100+0.0000001, 20 + 1);
phi_grid = 0.5*(phi_bins(1:end-1)+phi_bins(2:end));
l_grid = 0.5*(l_bins(1:end-1)+l_bins(2:end));
observed_mean_spikes_in_bins = zeros(20, 20);

for xx=1:20
    for yy=1:20
        timesinbin = (phi>phi_bins(xx)).*(phi<phi_bins(xx+1)).*(l>l_bins(yy)).*(l<l_bins(yy+1));
        if(sum(timesinbin)>50)
            observed_mean_spikes_in_bins(xx,yy) = mean(spikes(timesinbin>0));
        elseif(sum(timesinbin)<=50)
            observed_mean_spikes_in_bins(xx,yy) = nan;
        end
    end
end
    
figure(2)
subplot(1,2,2)
surface(l_grid, phi_grid, observed_mean_spikes_in_bins)
subplot(1,2,1)
h = imagesc(l_grid, phi_grid, observed_mean_spikes_in_bins)
set(h, 'AlphaData', ~isnan(observed_mean_spikes_in_bins))
set(gca,'YDir','normal')

test_r2 = struct;
test_r2.x = x_pos_OF;
test_r2.y = y_pos_OF;
test_r2.md = head_dir_OF.';
test_r2.spike = spikes;
test_r2.ts = time_intervals_OF.';
test_r2.mrl = 0;

test_r2.bait_angle = phi;
test_r2.bait_dist = l;
%test2 = EgocentricRatemap(test_r2, which_animal);
%plotEBC(test_r2, test2,4)
