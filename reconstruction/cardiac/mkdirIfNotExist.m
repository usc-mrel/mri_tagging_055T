function mkdirIfNotExist(dir)
    if ~exist(dir, 'dir')
        mkdir(dir);
    end
end