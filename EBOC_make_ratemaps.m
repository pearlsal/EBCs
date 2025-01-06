
N = size(data.spikemat,1);
%EBC_array = ones(length(EBC_array)); % just set them all to one if you wanna plot all neurons
% Run test computations and save stuff for each neuron
%observed_mean_spikes_in_bins = zeros(length(EBC_array), 30);
rr = 1;
cc = 1;
nn = 1;
for i=1:N
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
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBOC_ratemap_neuron', string(i), '.png');
            else
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBC_ratemap_neuron', string(i), '.png');
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
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBOC_ratemap_neuron', string(i), '.png');
            else
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBC_ratemap_neuron', string(i), '.png');
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
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBC_ratemap_neuron', string(i), 'even.png');
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
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBC_ratemap_neuron', string(i), 'even.png');
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
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/speed_filtered/EBC_ratemap_neuron', string(i), 'odd.png');
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
                filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/', session_name, '/not_speed_filtered/EBC_ratemap_neuron', string(i), 'odd.png');
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
            total_spikes = sum(data.spikemat(i,:));
            total_time_seconds = length(data.spikemat(i,:)) * binsize; % Total time in seconds
            firing_rate_hz = total_spikes / total_time_seconds;

            figure('Name', sprintf('Neuron %d', i), 'Position', [70 200 1300 400]);
            throwaway = plot_all_EBC(root, out, 1, 0, 1, 1, firing_rate_hz);

% Add a new subplot for firing rate
            subplot(1, 7, 7)  % Adjust the subplot layout as needed
            bar(firing_rate_hz)
            title('Firing Rate')
            xlabel('Neuron')
            ylabel('Rate (Hz)')
            set(gca, "FontSize", 7)
            % Add text to display the firing rate value
            text(1, firing_rate_hz, sprintf('%.2f Hz', firing_rate_hz), ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 8)

            % Adjust y-axis limit to accommodate the text
            ylim([0, firing_rate_hz * 1.1])
            % Adjust figure properties to fit the page
            set(gcf, 'PaperPositionMode', 'auto');
            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 13 4]);
            
            % Save the figure
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/EBOC_ratemap&CC_neuron', string(i), '_polar.pdf');
            exportgraphics(gcf, filename, 'Resolution', 600)
            saveas(gcf, filename)
            disp(filename)
            
            filename = strcat('/Users/pearls/Work/RSC_project/EBC_results/', which_animal, '/', which_channels, '/All/EBOC_ratemap&CC_neuron', string(i), '_polar.png');
            print(gcf, filename, '-dpng', '-r300')
            exportgraphics(gcf, filename)
            saveas(gcf, filename)
            
            close
            
            end
    end
    
end 