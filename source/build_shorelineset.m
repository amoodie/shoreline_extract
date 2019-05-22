function [] = build_shorelineset()
    %build_shorelineset processes the images and builds a set of csv shorelines from the images
    %
    % authored by Andrew J. Moodie and Brandee Carlson
    % 2015--2018
    % MIT License

    clear variables; 
    close all;

    %%%% SELECT THE RUNTIME PARAMETERS %%%%
    %
    % directory of the *raw* data folder (relative or absolute path to the folder)
    %     alternatively use cd <path to location> to execute elsewhere
    meta.directory = fullfile('..', 'data');
    %
    % do you want to make and save the processing image of the thresholding
    meta.make_thresh = true;
    %
    % do you want to make and save an RGB image too?
    meta.make_RGB = true;
    %
    % what are the coordinates for cropping to the delta extent
    meta.deltaULcoord = [633497, 4236964];
    meta.deltaLRcoord = [713380, 4162227];
    meta.crop_pts = [meta.deltaULcoord, meta.deltaLRcoord];
    meta.deltaLLcoord = [meta.deltaULcoord(1), meta.deltaLRcoord(2)]; % apex of the delta for cropping
    %
    %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
    
    % create the list of folders to loop through
    listing = dir(meta.directory); % all items in directory
    listing(ismember( {listing.name}, {'.', '..'})) = [];  % dont want . or ..
    direcoriesBool = [listing.isdir]; % logical of directories only
    folders = cellstr(  vertcat( listing(direcoriesBool).name )  ); % list of folder names only
    countfolders = length(folders);

    % loop through all the folder first to generate some metadata for sorting
    %     sort is important for making a movie if desired
    [sortidx, meta_sort] = get_sortorder(meta, folders, countfolders, meta.directory);
    
    % loop to process image --> shoreline 
    for i = 1:countfolders
        disp( ['operating on image folder: ', num2str(i), ' of ', num2str(countfolders)] )
        
        clear shoreline
        
        % grab the metadata
        meta = meta_sort{i};
        [meta.bandset] = set_bandset(meta.mission); % return bands to use as [thresh, R, G, B]
        
        % open the image that will be used for thresholding
        thresh_imgname = strcat(meta.name, '_B', meta.bandset(1), '.TIF');
        thresh_img = imread(char( fullfile(meta.imagefolder, thresh_imgname) ));
        
        % crop the image
        %    this section is Yellow River delta specific and would need to
        %    be rewritten for any other use.
        [meta.cropDim] = get_cropDim(meta.deltaULcoord, meta.deltaLRcoord, meta.res); % calculate the dimensions to crop at
        [thresh_crop] = crop_image(thresh_img, meta.ULcoord, meta.deltaULcoord, meta.cropDim, meta.res); % do the crop
        
        % threshold the image to a binary
        thresh_crop_adj = imadjust(thresh_crop, stretchlim(thresh_crop), [0 1], 1); % increase image contrast
        [thresh_val] = get_threshold(thresh_crop_adj); % determine the threshold value to use in binarization
        
        % find the shoreline
        [crop_edge] = find_shoreline(thresh_crop_adj, thresh_val, meta); % convert to binary and identify delta edge
        
        % concatenate edge into shoreline points
        [row, col] = find(crop_edge); % all points on shoreline
        shoreline_pts = horzcat(col, row); % shoreline point list, index ordered in array
        
        % convert the shoreline pts from img coords to geo-coords
        pivot_pt = max(shoreline_pts(:,2));
        shoreline_pts(:, 1) = (shoreline_pts(:, 1).*meta.res) + repmat(meta.deltaLLcoord(1), size(shoreline_pts, 1), 1);
        shoreline_pts(:, 2) = ((pivot_pt - shoreline_pts(:, 2)).*meta.res) + repmat(meta.deltaLLcoord(2), size(shoreline_pts, 1), 1);
        
        % sort the list into sequential shoreline trace
        [shoreline] = get_ordered(shoreline_pts);
        
        % make a RGB image
        if meta.make_RGB
            R_imgname = strcat(meta.name, '_B', meta.bandset(2), '.TIF');
            G_imgname = strcat(meta.name, '_B', meta.bandset(3), '.TIF');
            B_imgname = strcat(meta.name, '_B', meta.bandset(4), '.TIF');
            [RGB_fig] = plot_RGB(R_imgname, G_imgname, B_imgname, meta);
        end
        
        % write out shoreline data to .mat file
        save(fullfile( '..', 'output', strcat('shoreline_', datestr(meta.date, 'YYYY-mm-dd'), '.mat') ), 'shoreline')
        
        % write out shoreline data to a csv for use with other programs if needed
        outputname = fullfile( '..', 'output', strcat('shoreline_', datestr(meta.date, 'YYYY-mm-dd'), '.csv') );
        fid = fopen(outputname, 'w');
        fprintf(fid, 'X, Y\n');
        fclose(fid);
        dlmwrite(outputname, ...
            shoreline, '-append', 'precision', '%f');
        
        % write out metadata to .mat file
        save(fullfile( '..', 'output', strcat('meta_', datestr(meta.date, 'YYYY-mm-dd'), '.mat') ), 'meta')
        
        % end of loop
    end
    
end


function [sortidx, metadata_sort] = get_sortorder(meta, folders, countfolders, directory)
    % sort the folders for processing in order,
    % this is needed for making a movie while processing if desired
    
    % initialize
    dates = NaN(countfolders, 1); 
    
    % loop through to get metadata
    for i = 1:countfolders
        % concatenate to build the present folder and metadata file info
        imagefolder = fullfile(directory, folders(i)); 
        imagemetafile = fullfile(imagefolder, strcat(folders(i), '_MTL.txt'));
        [fidmetadata, message] = fopen(char(imagemetafile)); % fid of metadata file
        if fidmetadata < 0; disp(message); end
        
        % process the file to extract the metadata
        [metadata{i, 1}] = get_metadata(meta,  fidmetadata, imagefolder);
        
        % date is what we're after for sorting
        dates(i) = metadata{i}.date;
    end
    
    % perform the sort
    [~, sortidx] = sort(dates);
    [metadata_sort] = metadata(sortidx);
end


function [cropDim] = get_cropDim(ULcoord, LRcoord, res)
    % part of the cropping calculation, this is where I want to take 
    % a list of points and find the bounding box
    xDim = (LRcoord(1) - ULcoord(1)) / res; % x
    yDim = (ULcoord(2) - LRcoord(2)) / res; % y
    cropDim = [xDim yDim];
end


function [crop_img] = crop_image(image, ULcoord, cropULcoord, cropDim, resolution)
    ULidx = [(cropULcoord(1)-ULcoord(1)), (ULcoord(2)-cropULcoord(2))] ./ resolution; % WILL NOT WORK BECAUSE SPACING CHANGES???
    crop_img = imcrop(image, [ULidx(1) ULidx(2) cropDim(1) cropDim(2)]);
end


function [img_edge] = find_shoreline(img, thresh, meta)
    %find_shoreline performs the thresholding and extraction operations
    %
    % [img_edge] = find_shoreline(img, thresh, meta) is the main shoreline 
    %    extraction routine descibed in Moodie et al.
    %
    % inputs:
    %    img = the matrix of image intensity values to threshold with
    %    thresh = numeric intensity to binarize with
    %    meta = auxilary metadata for whether to make the plot 
    % outputs:
    %    img_edge = matrix same size as img with only the shoreline edge as true
    
    img_bw = im2bw(img, thresh);                        % threshold image
    img_fill = imfill(img_bw, 'holes');                 % fill it from the outside
    img_rms = ~bwareaopen(~img_fill, 30000);            % remove small isolated water-on-land objects
    img_rms2 = bwareaopen(img_rms, 500);                % remove small isolated land-in-water objects %%%%%%%% 2000
    img_structel = strel('disk', 5);                    % build structural object (something like a filter) of ('shape', size)
    img_open = imopen(img_rms2, img_structel);          % morphological closure with structure
    img_structel2 = strel('disk', 50);                  % build structural object (something like a filter) of ('shape', size)
    img_close = imclose(img_open, img_structel2);       % morphological closure with structure
    img_rms3 = bwareaopen(img_close, 10000);            % remove small isolated land-in-water objects
    img_pad = padarray(img_rms3, [1 0], 1, 'post');     % add row to end
    img_pad = padarray(img_pad, [0 1], 1, 'pre');       % add col to front
    img_fill2 = imfill(img_pad, 'holes');               % fill it from the outside
    img_unpad = img_fill2(1:end-1, 2:end);              % remove padding from outside
    img_fill3 = bwareafilt(img_unpad, 1, 'largest');    % retain only largest 
    img_edge = edge(img_fill3, 'sobel');                % find edge
    
    if meta.make_thresh
        fig = figure();
        subplot(2,3,1)
        name = 'raw image';
            imshow(img)
            title(name)
        subplot(2,3,2)
        name = 'apply threshold';
            imshow(img_bw)
            title(name)
        subplot(2,3,3)
        name = 'flood pixels';
            imshowpair(img_bw, img_fill)
            title(name)
        subplot(2,3,4)
        name = 'remove small objects';
            imshowpair(img_fill, img_rms2)
            title(name)
        subplot(2,3,5)
        name = 'morphological open and close';
            imshowpair(img_rms2, img_close)
            title(name)
        subplot(2,3,6)
        name = 'extract shoreline';
            imshow(img)
            hold on
            [row, col] = find(img_edge);
            shoreline = horzcat(col, row);
            plot(shoreline(:,1), shoreline(:,2), 'r.')
            title(name)
    end
    
end


function [bandset] = set_bandset(mission)
    %set_bandset gives the bands to use for a given mission
    % 
    % [bandset] = set_bandset(mission)
    % inputs:
    %    mission: 'LANDSAT_X' where X is the number. This should be
    %    stripped directly from the metadata
    % outputs:
    %    bandset: the set to use in the code later, 4 character string, [thresh, R, G, B]
    
    switch mission
        case {'LANDSAT_1', 'LANDSAT_2', 'LANDSAT_3'}
            bandset = ['7' '6' '5' '4'];
        case {'LANDSAT_4', 'LANDSAT_5'}
            bandset = ['7' '3' '2' '1'];
        case 'LANDSAT_7'
            bandset = ['7' '3' '2' '1'];
        case 'LANDSAT_8'
            bandset = ['7' '4' '3' '2'];
    end
    
end


function [meta] = get_metadata(meta, fidmetadata, imagefolder)
    
    % read the raw text file into a cell array of strings
    [metadata] = textscan(fidmetadata, '%s','delimiter', '\n');
    
    % strip out the desired info
    % date aquired
    [date] = strip_from_meta(metadata, 'DATE_ACQUIRED', 'str');
    [meta.date] = datenum(date);
    
    % cloud cover from file
    [meta.clouds] = strip_from_meta(metadata, 'CLOUD_COVER', 'num');
    
    % upper left geo-coordinate of image
    [ULXcoord] = strip_from_meta(metadata, 'CORNER_UL_PROJECTION_X_PRODUCT', 'num');
    [ULYcoord] = strip_from_meta(metadata, 'CORNER_UL_PROJECTION_Y_PRODUCT', 'num');
    [meta.ULcoord] = [ULXcoord, ULYcoord];

    % what mission shot the image
    [meta.mission] = strip_from_meta(metadata, 'SPACECRAFT_ID', 'str');
    
    % image resolution
     [meta.res] = strip_from_meta(metadata, 'GRID_CELL_SIZE_REFLECTIVE', 'num');

    % what is the scene ID (name)
    [meta.name] = strip_from_meta(metadata, 'LANDSAT_PRODUCT_ID', 'str');

    % image folder path
    [meta.imagefolder] = imagefolder;
    
end


function [value] = strip_from_meta(metadata, keystring, valuetype)
    %strip_from_meta strips out the value from metadata for a given string
    %
    % [value] = strip_from_meta(metadata, keystring, valuetype)
    % takes:
    %     metadata = array of metadata from file
    %     keystring = the string to search for and strip out the value
    %     valuetype = the type to return ('str' or 'num')
    % returns:
    %     value = corresponding entry in metadata formatted as valuetype
    
    % find the index and grab that line
    idx = find(~cellfun(@isempty,strfind(metadata{1,1}, keystring)) == 1);
    idx = idx(1);
    str = metadata{1,1}(idx);
    
    % split the string at the equals sign and strip quotes and whitespace
    splt = strsplit(char(str),'=');
    noquotes = strrep(splt(2), '"', '');
    nowhite = strtrim(noquotes{:});
   
    % process to the desired type
    if strcmp(valuetype, 'str')    
        value = nowhite;
    elseif strcmp(valuetype, 'num')
        value = str2double(nowhite);
    else
        error('invalid valuetype')
    end

end


function [thresh] = get_threshold(img)
    %get_threshold gets the threshold used to binarize the image
    %
    % works following the method described in Moodie et al.
    %
    
    img = im2double(img); % convert to double type for math manipulation
    img_long = reshape(img, size(img,1) * size(img,2), 1); % make long for histogram
    [Hcount, Hbin] = histcounts(img_long, 40); % get histogram count of intensity
    dx = Hbin(2) - Hbin(1); % histogram spacing
    dcdx = (Hcount(2:end) - Hcount(1:end-1)) ./ dx; % spatial gradient in count
    [maxy, maxidx] = max(Hcount); % value and index of max count
    range = [maxidx find(dcdx(maxidx:end)>=0, 1, 'first')+maxidx-2]; % bin range from max to first upslope
    m = mean(dcdx(range(1):range(2))); % mean slope across the range
    maxx = maxidx*dx; % max x location to start from
    yint = maxy - (m * (maxx)); % solve for y intercept of line
    xint = (-1*yint)/m; % project down to x axis
    thresh = xint; % that's the intensity threshold
end


function [line] = get_ordered(pointlist)
    %get_ordered sorts the points in the list into a sequential order along the shore
    %
    
    used = false(size(pointlist, 1), 1);
    cnt = 1;
    last = cnt;
    used(last) = true;
    line(cnt, :) = pointlist(cnt, :);
    lastdist = 0;
    while lastdist <= 60
        used(last) = true;
        
        distances = sqrt(((pointlist(last, 1) - pointlist(:, 1)).^2 + (pointlist(last, 2) - pointlist(:, 2)).^2));
        distances(used) = Inf;
        [lastdist, last] = min(distances);
        
        line(cnt+1, :) = pointlist(last, :);
        cnt = cnt + 1;
        if cnt > 15000
            error
        end
    end
    line = line(1:end-1, :);
end


function [thefig] = plot_RGB(R_imgname, G_imgname, B_imgname, meta)
    %plot_RGB manipulates and plots the RGB image
    [R_img] = imread(char(fullfile(meta.imagefolder, R_imgname))); % open red image
    [G_img] = imread(char(fullfile(meta.imagefolder, G_imgname)));
    [B_img] = imread(char(fullfile(meta.imagefolder, B_imgname)));
    [R_crop] = crop_image(R_img, meta.ULcoord, meta.deltaULcoord, meta.cropDim, meta.res); % do the crop
    [G_crop] = crop_image(G_img, meta.ULcoord, meta.deltaULcoord, meta.cropDim, meta.res); 
    [B_crop] = crop_image(B_img, meta.ULcoord, meta.deltaULcoord, meta.cropDim, meta.res); 
    clip = [0.2, 0.8];
    [R_adj] = imadjust(R_crop, stretchlim(R_crop, clip), [0.2 0.8], 1); % increase image contrast
    [G_adj] = imadjust(G_crop, stretchlim(G_crop, clip), [0.2 0.8], 1);
    [B_adj] = imadjust(B_crop, stretchlim(B_crop, clip), [0.2 0.8], 1);
    RGB = cat(3, R_adj, G_adj, B_adj);
        
    thefig = figure();
    imshow(RGB)
    title(datestr(meta.date))
    
end

