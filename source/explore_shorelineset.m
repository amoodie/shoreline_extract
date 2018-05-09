function explore_shorelineset()
    %explore_shorelineset explores the shoreline data set
    %
    % authored by Andrew J. Moodie and Brandee Carlson
    % 2015--2018
    % MIT License

    clear variables;
    close all;

    %%%% SELECT THE RUNTIME PARAMETERS %%%%
    %
    % directory of the *processed* data folder (relative or absolute path to the folder)
    %     alternatively use cd <path to location> to execute elsewhere
    meta.directory = fullfile('..', 'output');
    %
    %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%

    % load data
    [data] = load_data(meta.directory);
%     [auto.data] = load_data('auto');

%     T = load('./shorelines/deltaLLcoord.mat');
%     deltaCROPLLcoord = T.deltaLLcoord; % I have no idea what this is for
    deltaLLcoord = [611066, 4149751]; % this is Lijin, used for calculating the delta radius
    
    % prepare the channel centerline
    %    this below routine could be repeated for multiple lines and then
    %    the multiple lines could be processed in the loop too
    linename = fullfile('..', 'data', 'qinshuigou_channelline.csv');
    qing.raw = csvread(linename, 1, 0); % read xy points from csv file 
    qing.xy = [qing.raw(:,1), qing.raw(:,2)];
    qing.between = [0 sqrt((qing.xy(2:end, 1) - qing.xy(1:end-1, 1)).^2 + ...
        (qing.xy(2:end, 2) - qing.xy(1:end-1, 2)).^2)'];
    qing.along = cumsum(qing.between);
   
    % loop through all the shorelines
    for i = 1:size(data,1)
        % make subset to reduce the size of the data to check for intersection
        [data(i).lobe] = get_subset(data(i).shoreline, qing.xy); % subset of shoreline only in lobe area
                
        % calculate the mean radius and intersection
        [data(i).radius] = get_radius(data(i).shoreline, deltaLLcoord); % radius of the entire delta
        [data(i).qingint] = get_intersection(data(i).lobe, qing, data(i).meta.res); % intersection with Qingshuigou
    
    end
    
    % convert data into a single table for more convenient manipulation
    %   this step requires manually adding items to the table with the 
    %   data you want to use in analysis
    table_vars = {'date', 'deltaradius', 'qingint'};
    all.tab = array2table(horzcat(arrayfun(@(x) x.meta.date, data), arrayfun(@(x) x.radius.mean, data), arrayfun(@(x) x.qingint.distalong, data)), ...
        'VariableNames', table_vars);    
    
    % convert data into tables spanning time ranges for convenience of building models
    qing.tab = make_rangetable(all, [datenum('1976','YYYY'), datenum('07/31/1997','mm/dd/YYYY')]); % only data from qingshuigou lobe development
        % note that in this example case, the data in all.tab is identical to qing.tab
    
    % make the models
    delta.model = make_rate(all.tab, 'deltaradius ~ date');
    qing.model = make_rate(qing.tab, 'qingint ~ date');
    
    % demonstration of intersection extraction (last in loop)
    figure();
    subplot(2, 4, [1:2 5:6]); hold on;
        plot(data(i).shoreline(:,1), data(i).shoreline(:,2), 'k-', 'LineWidth', 1.5) % plot the shoreline
        plot(qing.xy(:,1), qing.xy(:,2), 'g-', 'LineWidth', 1.5) % plot the intersection line
        plot(data(i).qingint.shorelinecoords(1), data(i).qingint.shorelinecoords(2), 'r*', 'MarkerSize', 8) % plot the shoreline determined intersection
        plot(data(i).qingint.linecoords(1), data(i).qingint.linecoords(2), 'r*', 'MarkerSize', 8) % plot the centerline determined intersection
        xl = xlim; yl = ylim;
        xlim([min(data(i).shoreline(:,1)), xl(2)])
        ylim([min(data(i).shoreline(:,2)), yl(2)])
        set(gca, 'xTickLabels', xticks/1000)
        set(gca, 'yTickLabels', yticks/1000)
        xlabel('Easting UTM Zone 50S (km)')
        ylabel('Northing UTM Zone 50S (km)')
        title('find intersection')
        axis equal
        axis tight
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10)
    subplot(2, 4, [3:4 7:8]); hold on;
        plot(data(i).shoreline(:,1), data(i).shoreline(:,2), 'k-', 'LineWidth', 1.5) % plot the shoreline
        plot(qing.xy(:,1), qing.xy(:,2), 'g-', 'LineWidth', 1.5) % plot the intersection line
        plot(data(i).qingint.shorelinecoords(1), data(i).qingint.shorelinecoords(2), 'r*', 'MarkerSize', 8) % plot the shoreline determined intersection
        plot(data(i).qingint.linecoords(1), data(i).qingint.linecoords(2), 'r*', 'MarkerSize', 8) % plot the centerline determined intersection
        axis equal
        xlim([data(i).qingint.linecoords(1)-60, data(i).qingint.linecoords(1)+60])
        ylim([data(i).qingint.linecoords(2)-60, data(i).qingint.linecoords(2)+60])
        title('zoomed')
        set(gca, 'xTickLabel', [], 'yTickLabel', [])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10)
        
    % demonstration of the rate plot
    figure()
    s1 = subplot(1, 2, 1); hold on;
        plot(repmat(datenum('1976','YYYY'), 1, 2), [20 90], ':', 'Color', [0 0 0]) % lower bound of qingshuigou development
        plot(repmat(datenum('1997','YYYY'), 1, 2), [20 90], ':', 'Color', [0 0 0]) % upper bound of qingshuigou development
        plot(all.tab.date, all.tab.deltaradius ./ 1000, 'r.', 'MarkerSize', 10) % plot the data
        plot(delta.model.xs, delta.model.ys ./ 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0]) % plot the model
        datetick('x', 'YYYY')
        xlim([datenum('1970','YYYY') datenum('2000','YYYY')])
        params1 = {['rate = ', num2str(round(delta.model.m * 365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(delta.model.r2, 2))]};
        format1 = sprintf('%s\n', params1{:});
        text(0.1, 0.1, format1(1:end-1), ...
            'Color', [0 0 0], 'Parent', s1, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('YR delta development')
        xlabel('year')
        ylabel('mean delta radius (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10)
    s2 = subplot(1, 2, 2); hold on;
        plot(repmat(datenum('1976','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0]) % lower bound of qingshuigou development
        plot(repmat(datenum('1997','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0]) % upper bound of qingshuigou development
        plot(qing.tab.date, qing.tab.qingint ./ 1000, 'r.', 'MarkerSize', 10) % plot the data
        plot(qing.model.xs, qing.model.ys ./ 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0]) % plot the model
        datetick('x', 'YYYY')
        xlim([datenum('1970','YYYY') datenum('2000','YYYY')])
        params2 = {['rate = ', num2str(round(qing.model.m * 365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(qing.model.r2, 2))]};
        format2 = sprintf('%s\n', params2{:});
        text(0.5, 0.1, format2(1:end-1), ...
            'Color', [0 0 0], 'Parent', s2, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('Qingshuigou lobe development')
        xlabel('year')
        ylabel('distance downstream from datum (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10)
    set(gcf, 'Pos', [50 100 800 400], 'PaperPositionMode', 'auto')
    
    % output the data to a folder
    datatable = all.tab; %#ok<NASGU>
    save(fullfile( '..', 'output', 'processed_data.mat'), 'datatable');
    
end


function [radius] = get_radius(shoreline, deltaLLcoord)
    %get_radius returns a radius object
    %
    %
    
    % initialize
    npts = size(shoreline, 1);
    obs = NaN(npts, 1);
    
    % loop through all the points along the shoreline
    for i = 1:npts
        % calculate distance from apex to point
        xdiff = shoreline(i, 1) - deltaLLcoord(1);
        ydiff = shoreline(i, 2) - deltaLLcoord(2);
        obs(i) = sqrt(xdiff^2 + ydiff^2);
    end
    
    % calculate and store data into object
    radius.obs = obs;
    radius.mean = nanmean(obs);
    radius.std = nanstd(obs);
    
end


function [intersection] = get_intersection(shoreline, line, res)
    %get_intersection returns an intersection object
    %

    % number of points in the shoreline
    npts = size(shoreline, 1);
    
    % preallocate
    mindists = NaN(npts, 1);
    mindistsidx = NaN(npts, 1);
    
    % loop through every point in shoreline to determine if it intersects with line
    for i = 1:npts
        % calculate the distance from this shoreline pt to every point in line
        distances = sqrt(((shoreline(i, 1) - line.xy(:, 1)).^2 + (shoreline(i, 2) - line.xy(:, 2)).^2));
        
        % record the minimum distance for that shoreline point
        [mindists(i), mindistsidx(i)] = nanmin(distances);
        
    end
    % find the minimum distance of any of the points
    [mindist, mindistidx] = nanmin(mindists);
    
    % validate it is an intersection based on the grid spacing of the raster
    if sqrt(res^2) < mindist
        warning('intersection identified, but may not be real, returning NaN')
        mindist = NaN;
    end

    intersection.mindist = mindist; % minimum distance used to find the intersection
    intersection.shorelinecoords = shoreline(mindistidx, :); % the point of the intersection in the shoreline dataset
    intersection.linecoords = line.xy(mindistsidx(mindistidx), :);
    intersection.distalong = line.along(mindistsidx(mindistidx));
                
end


function [model] = make_rate(table, formula) 
    %make_rate makes a linear rate model for given formula
    %
    
    model.mdl = fitlm(table, formula);
    model.b = model.mdl.Coefficients.Estimate(1);
    model.m = model.mdl.Coefficients.Estimate(2);
    model.bserr = bootstrp(35, @mean, model.mdl.Residuals.Raw);
    model.CI = coefCI(model.mdl);
    model.err = mean(abs(model.CI(2,:)-model.m));
    model.r2 = model.mdl.Rsquared.ordinary;
    model.xs = linspace(min(table.date), max(table.date), 10);
    model.ys = ((model.m .* model.xs) + model.b);
end


function [table] = make_rangetable(all, range)
    %make_rangetable subsets the total dataset, yielding only data within date range
    %
    % inputs:
    %   all = the total datatable
    %   range = a two element vector with datenum cooresponding to lower 
    %       and upper bound of data you want retained
    % outputs:
    %   table = table containing all data in range
    %
    
    keep = and(all.tab.date >= range(1), all.tab.date <= range(2));
    table = all.tab(keep, :);
    
end


function [subset] = get_subset(shoreline, line)
    %get_subset finds the minimum number of points in data that could possibly intersect with line
    %
    % useful for reducing compute time in the intersection code
    %
    
    % make logicals of what's beyond the line limits
    trash_left = (shoreline(:,1) < min(line(:,1))); % left
    trash_right = (shoreline(:,1) > max(line(:,1))); % right
    trash_bottom = (shoreline(:,2) < min(line(:,2))); % bottom
    trash_top = (shoreline(:,2) > max(line(:,2))); % top
    
    % retain what is never outside limits
    keep = ~ or(  or(trash_left, trash_right), or(trash_bottom, trash_top)  );
    subset = shoreline(keep, :); 
    
end


function [data] = load_data(directory)

    % get the directory listing and identify the shoreline_.mat files
    listing = dir(directory); % all items in directory
    listing(ismember( {listing.name}, {'.', '..'})) = [];  % dont want . or ..
    filesBool = ~[listing.isdir]; % logical of files only
    listing = listing(filesBool);
    listingnames = arrayfun(@(x) x.name, listing, 'Unif', 0);
    listingshoreline = cellfun(@(x) and(contains(x, '.mat'), contains(x, 'shoreline_')), listingnames);
    shorelinefiles = listing(listingshoreline);
    listingmeta = cellfun(@(x) and(contains(x, '.mat'), contains(x, 'meta_')), listingnames);
    metafiles = listing(listingmeta);
    
    if size(shorelinefiles) ~= size(metafiles)
        error('number of shorelines and metadata files does not match')
    else
        countfiles = length(shorelinefiles);
    end
    
    % preallocate
    data = cell(countfiles, 2);
    
    % loop through and load data into a cell array
    for i = 1:countfiles
        ishoreline = load(fullfile(shorelinefiles(i).folder, shorelinefiles(i).name), 'shoreline');
        data{i, 1} = ishoreline.shoreline;
        imeta = load(fullfile(metafiles(i).folder, metafiles(i).name), 'meta');
        data{i, 2} = imeta.meta;
    end
    
    % convert to a table structure
    data = cell2struct(data, {'shoreline', 'meta'}, 2);
    
    % ensure the data are sorted in time sequence
    dates = cellfun(@(x) x.date, {data.meta});
    [~, datessortidx] = sort(dates);
    data = data(datessortidx);
    
end

