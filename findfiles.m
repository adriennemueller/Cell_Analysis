%FINDFILES Find all files in a folder, matching the pattern.
function [ files ] = findfiles( folder, pattern, ret_direc )
    files = [];
    matches = dir(fullfile(folder, pattern));
    
    % Don't return folders
    if (ret_direc == 0)
        matches = matches([matches.isdir] == 0);
    end

    % Return full file names
    for i = 1: length(matches)
        match = matches(i);
        match.name = fullfile(folder, match.name);
        files = [files ; match];
    end
    
    folders = dir(folder);
    folders = folders([folders.isdir] == 1);

    for i = 1: length(folders)
        child = folders(i);
        if ~strcmp(child.name, '.') && ~strcmp(child.name, '..')
            deeperFiles = findfiles( fullfile(folder, child.name), pattern, ret_direc );
            files = [files ; deeperFiles];
        else
            %disp(['Skipping ' child.name]);
        end
    end

end

