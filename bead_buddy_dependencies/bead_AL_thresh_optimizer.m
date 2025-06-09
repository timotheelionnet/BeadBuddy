
% my_path = "H:\Finn_DriveH\20240226_PPARG_Staggered_Timecourse\beads\session2_beads\split_channels\C1-session2_beads_1.tif";
% cur_thresh = 30;

% bead_cnt = count_binary_beads_for_pipe(my_path);


disp('~~~~~~~~~')
disp('Running AL threshhold optimizer')
disp(strcat('Approximate bead count is: ', string(n_beads)))

my_path = curImPath;
% cur_thresh = 1;
bead_cnt = n_beads;

loc = perform_detection_on_single_img_once_Finn(my_path, 0, 0, 1,2, 'SD', cur_thresh); %Manual training on first image for ith channel

prev_spot_cnt = height(loc);
%% while loop to optimize spot count

we_stuck = 0 ;% are we not getting any better and stuck in the loop...
discrep = bead_cnt - prev_spot_cnt;
disp('discrepancy is')
disp(strcat(string(discrep), ' beads'))


% are we within 10% of the bead count?
if abs(discrep) > 0.1 *bead_cnt
    disp('optimizing detection thresh')
    opto = true;
    
    
    iter_cnt = 0; % track futile while loop iterations
   
    while opto == true
        
        % if we AL detect too many, then discrep is negative
        % we need to raise thresh
        if discrep < 0

            % remember out thresh values and report
            prev_thresh = cur_thresh;

            cur_thresh = cur_thresh + 0.5;
        

            
            disp(fprintf('Adjusting Thresh from %d to %d', [prev_thresh, cur_thresh]));
            

        elseif discrep > 0

            % remember out thresh values and report
            prev_thresh = cur_thresh;

            cur_thresh = cur_thresh - 0.5;
        

            
            disp(fprintf('Adjusting Thresh from %d to %d\n', [prev_thresh, cur_thresh]));
        end
     
        % run AL again with new thresh
        loc = perform_detection_on_single_img_once_Finn(my_path, 0, 0, 1,2, 'SD', cur_thresh); %Manual training on first image for ith channel
        
        % did we get better?
        cur_spot_cnt = height(loc);
    
        cur_discrep = bead_cnt - cur_spot_cnt;

        if abs(cur_discrep) < abs(discrep)
            disp('Getting closer');
            % 'discrep'
            % cur_discrep
        else
            % somehow got worse, just keep the last thresh
            cur_thresh = prev_thresh;
            break
            

        end
        
        % if we are still within 5 % of the last count, then we will only
        % try two more times
        if (prev_spot_cnt - 0.05* cur_spot_cnt < cur_spot_cnt ) && ( cur_spot_cnt < prev_spot_cnt + 0.05* cur_spot_cnt) 
            iter_cnt = iter_cnt +1;
        end
        
        if iter_cnt >3
            we_stuck = 1;
        end

        if cur_thresh == 0

            we_suck =1;
        end

        if we_stuck ==1
            break
        end
    
        prev_spot_cnt = height(loc);
    end
end