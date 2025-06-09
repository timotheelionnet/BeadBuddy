

% table must ahve 
function nn_res = get_nn_matrix_from_2arrays_with_ref_coords(tab1, tab2, nDims, nnThresh, mutualOnly)

% if mutalOnly is off, then we will only look for nns to tab1

% tab1 = locRef;
% tab2 = locReg;
% nDims = 3;
% nnThresh = 100;
% mutualOnly =1;

% OUTPUT is dr, dx, dy, dz, x1, y1, z1




nn_res = [];


%Row counter of ith NN matrix
myDistMat = compute_distMatrix_from_loc(tab1, tab2, nDims);


% track misses
missed = 0;
            
%Loop through each row (ie ref spot) and find nearest neighbor
for p = 1:height(myDistMat)

    cur_ref_loc = tab1(p, 1:3); % get xyz of ref loc
    
    %Store index of nn as nnIdx
    nnIdx = find(myDistMat(p,:) == min(myDistMat(p,:)));
                
    %Is nn within our nn threshold?
    if myDistMat(p, nnIdx) <= nnThresh

        if mutualOnly == 1
    
            % Is spot p also the NN for spot nnIdx?
            mutualNnIdx = find(myDistMat(:,nnIdx) == min(myDistMat(:,nnIdx)));
    
            if mutualNnIdx ==p
                
                %add it to the nn dataframe: dr, dx, dy, dz 
                
                drxyz  = [myDistMat(p, nnIdx), tab2(nnIdx,1:3) - tab1(p,1:3), cur_ref_loc];

                nn_res = [nn_res ; drxyz ];

            else 
                missed = missed + 1;
              
                
            
            end
        else

            drxyz  = [myDistMat(p, nnIdx), tab2(nnIdx,1:3) - tab1(p,1:3), cur_ref_loc];

            nn_res = [nn_res ; drxyz ];
                
            
        end
    else
        %This ref spot did not have a match in the reg channel under the
        %thresh
        missed = missed +1;
        
    end          
        

       
end
       
    
% Store number of unmatched tab1 spots 
disp(strcat( string(missed/height(tab1)), ' of table1 spots were not matched (', string(missed), ' spots)'));




end