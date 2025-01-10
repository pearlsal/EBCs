
% Read the .mat file with all the binned rat data
folder_loc = '/Users/pearls/Work/RSC_project/';
which_animal = "PreciousGrape"; %{'Luke', 'Arwen', 'Tauriel'} "PreciousGrape"; %{'Luke', 'Arwen', 'Tauriel', 'PreciousGrape'}
which_session = "OF1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["OF1", "OF2", "c1", "c2",];
chase_or_chill = "chill";
different_sessions = true; % if you wanna concatenate over both OF and (chilling bouts) c's, or filter specific parts of sessions (chasing bouts)
shared_cells = true; % this one must be true if using concatenated sessions, and should probably be true otherwise as well
which_channels = 'RSC'; % {RSC, SC, ALL}
EBC_or_EBOC = "EBOC"; % EgocentricBoundary-vectorCell or EgocentricBaitOrientedCell, ie. plot bearing to wall or to bait

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
folder_loc = '/Users/pearls/Work/RSC_project/';
which_animal = "PreciousGrape"; %{'Luke', 'Arwen', 'Tauriel'} "PreciousGrape"; %{'Luke', 'Arwen', 'Tauriel', 'PreciousGrape'}
which_session = "c1"; %Luke: {'session2chasing_solo'}, Arwen: {OF1, OF2, c1, c2, c4}, Tauriel: {OF1, c1, c2, c4, c5}
binsize = 0.008; %binsize in s, {0.008, 0.025, 0.050, 0.100}, 0.008 = same as recording, ie. 8ms ish
speed_threshold = true; % filter if the animal isn't moving
concat_sessions = true; % if we want to concatenate OF(any) sessions 
concat_which = ["c1", "c2"];
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

            original_tracking_interval = 0.008332881468650438;
            new_tracking_interval = binsize;
            interval_diff = floor(new_tracking_interval/original_tracking_interval);
            if interval_diff == 0
                interval_diff = 1;
            end
            
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