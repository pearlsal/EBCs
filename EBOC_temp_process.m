folder_loc = '/Users/pearls/Work/RSC_project/';
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
