
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
    %%%%%%%%%%%%%%%%%%%%%%%%
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
    total_spikes = sum(data.spikemat(i,:));
    total_time_seconds = length(data.spikemat(i,:)) * binsize; % Total time in seconds
    firing_rate_hz = total_spikes / total_time_seconds;

    figure('Name', sprintf('Neuron %d', i), 'Position', [70 200 1300 400]);
    CC_shift = plot_all_EBC(root, out, 1, 0, 1, 0, firing_rate_hz);
    
    CC_shift_mat(:,i) = CC_shift;   

    filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/EBC_ratemap&CC_neuron', string(i), '_polar.pdf');
    disp(filename)
    %exportgraphics(gcf,filename, 'Resolution', 600) 
    %saveas(gcf, filename)
    filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/EBC_ratemap&CC_neuron', string(i), '_polar.png');
    %exportgraphics(gcf,filename) 
    saveas(gcf, filename)
    
    close
end
