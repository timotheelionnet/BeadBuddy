function [] = make_new_dir(new_dir)



if~exist(new_dir, 'dir')

    mkdir(new_dir)

    disp('~~~~~')
    disp("Created new dir")
    disp(new_dir)
    disp('~~~~~')

end



end