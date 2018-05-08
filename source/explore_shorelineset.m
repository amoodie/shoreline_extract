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
    delta.model = make_rate(qing.tab, 'qingint ~ date');
    
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
        set(gca, 'xTickLabel', [], 'yTickLabel', [])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10)
        
    
    % demonstration of the rate plot
    figure()
    subplot(1, 2, 1); hold on;
        plot(repmat(datenum('1976','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0]) % lower bound of qingshuigou development
        plot(repmat(datenum('1997','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0]) % upper bound of qingshuigou development
        plot(arrayfun(@(x) x.meta.date, data), arrayfun(@(x) x.qingint.distalong, data), 'r.', 'MarkerSize', 10) % plot the data
        plot(preQ8.deltaeval.xs, preQ8.deltaeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0]) % plot the model
        datetick('x', 'YYYY')
        xlim([datenum('1970','YYYY') datenum('2000','YYYY')])
        params1 = {['rate = ', num2str(round(preQ8.deltaeval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(preQ8.deltaeval.r2, 2))]};
        format1 = sprintf('%s\n', params1{:});
        annot1 = text(0.5, 0.1, format1(1:end-1), ...
            'Color', [0 0 0], 'Parent', s1, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean delta radius')
        xlabel('year')
        ylabel('mean delta radius (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
    
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

    npts = size(shoreline, 1);
%     intersection = cell(size(shoreline, 1), 4);
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


function [fig] = all_shorelines(manu, auto, deltaLLcoord, deltaCROPLLcoord, mask, fig)
    all.data = vertcat(manu.data, auto.data);
    filename = '../maps/intersection_channelline.csv';
    AC1channelline = csvread(filename, 1, 0);
    filename = '../maps/intersection_AC2.csv';
    AC2channelline = csvread(filename, 1, 0);
    filename = '../maps/intersection_2016.csv';
    modernchannelline = csvread(filename, 1, 0);

    [cmap] = colormap_fun(size(all.data, 1), 1); % second input is start idx
    [cmap] = parula(size(all.data, 1));
    figure(fig)
        cla
        hold on
            for p = 1:size(all.data, 1)
                plot(all.data{p, 2}(:, 1), all.data{p, 2}(:, 2), '.', 'Color', cmap(p,:))
            end
            plot(AC1channelline(:, 1), AC1channelline(:, 2), '.', 'Color', [0 0 0])
            plot(AC2channelline(:, 1), AC2channelline(:, 2), '.', 'Color', [0 0 0])
            plot(modernchannelline(:, 1), modernchannelline(:, 2), '.', 'Color', [0 0 0])
%             plot([mask.acoord(1) mask.maxh(1)], [mask.acoord(2), mask.maxh(2)], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.2)
%             plot([mask.acoord(1) mask.maxl(1)], [mask.acoord(2), mask.maxl(2)], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.2)
        text(all.data{1, 2}(1350, 1), all.data{1, 2}(970, 2), datestr(all.data{1, 1}, 'YYYY'))
        text(698657, 4195500, datestr(all.data{end, 1}, 'YYYY'))
%         title('shorelines')
        axis('equal')
        xlim([deltaCROPLLcoord(1) 720000])
        ylim([deltaCROPLLcoord(2) 4233000])
        xTicks = get(gca, 'XTick');
        xTicks = xTicks(1:2:end);
        set(gca, 'XTick', (xTicks), 'XTickLabel', (xTicks))
        yTicks = get(gca, 'YTick')';
        yTicks = yTicks(1:1:end);
        box on
        set(gca, 'YTick', (yTicks), 'YTickLabel', num2str(yTicks)) 
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
        colorbar('Ticks', [0 1], 'TickLabels',{'1855','2016'})
    print('-dpng', '-r300', './figs/XOM_shoreline_map.png');
    
end


function [movieFIG] = movie_fig(manu, auto, deltaLLcoord, deltaCROPLLcoord, mask, movieFIG)
    all.data = vertcat(manu.data, auto.data);
%     [cmap] = parula(size(all.data, 1));
    [cmap] = parula(77);
    filename = '../maps/intersection_channelline.csv';
    AC1channelline = csvread(filename, 1, 0);

    figure(movieFIG)
    axis('equal')
    xlim([deltaCROPLLcoord(1) 720000])
    ylim([deltaCROPLLcoord(2) 4233000])
    xTicks = get(gca, 'XTick');
    xTicks = xTicks(1:2:end);
    set(gca, 'XTick', (xTicks), 'XTickLabel', (xTicks))
    yTicks = get(gca, 'YTick')';
    yTicks = yTicks(1:1:end);
    box on
    set(gca, 'YTick', (yTicks), 'YTickLabel', num2str(yTicks)) 
    set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(movieFIG,'Visible','on');
    set(movieFIG, 'Pos', [100 100 700 600]);
    set(movieFIG, 'PaperPositionMode', 'auto')
    colorbar('Ticks', [0 1], 'TickLabels',{'1973','1996'}, 'position', [0.75 0.75 0.05 0.1]);
    cla
    hold on
        for p = 3:81
            plot(all.data{p, 3}(:, 1), all.data{p, 3}(:, 2), '-', 'Color', cmap(p,:))
            trace = plot(AC1channelline(:, 1), AC1channelline(:, 2), '-', 'Color', [0 0 0]);
            drawnow 
%             saveas(movieFIG, sprintf(strcat('./figs/shorelinemovie/%03d.png'), p-2))
            delete(trace)
        end
            
%             plot(AC2channelline(:, 1), AC2channelline(:, 2), '.', 'Color', [0 0 0])
%             plot(modernchannelline(:, 1), modernchannelline(:, 2), '.', 'Color', [0 0 0])
%             plot([mask.acoord(1) mask.maxh(1)], [mask.acoord(2), mask.maxh(2)], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.2)
%             plot([mask.acoord(1) mask.maxl(1)], [mask.acoord(2), mask.maxl(2)], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.2)


end


function [fig] = JGR_meanrad(manu, auto, qing, preQ8, fig)
    figure(fig)
    [colorOrder] = get(gca, 'ColorOrder');
    barColor = [0.8 0.8 0.8];
    slct = auto.tab.date < datenum('1996', 'YYYY');
    s1 = subplot(1, 2, 1);
        cla
        hold on
%             plot(repmat(datenum('1854','YYYY'), 1, 2), [0 100], ':', 'Color', [0 0 0])
%             plot(repmat(datenum('1997','YYYY'), 1, 2), [0 100], ':', 'Color', [0 0 0])
            e1 = errorbar([manu.tab.date; auto.tab.date(slct)], [manu.tab.meandeltarad; auto.tab.meandeltarad(slct)] / 1000, ...
                [manu.tab.meandeltaraderr; auto.tab.meandeltaraderr(slct)] / 1000, 'LineStyle', 'none', 'Color', barColor);
%             m = plot(manu.tab.date, manu.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
%             a = plot(auto.tab.date(slct), auto.tab.meandeltarad(slct) / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            f = plot(preQ8.tab.date, preQ8.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
            plot(preQ8.deltaeval.xs, preQ8.deltaeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1840','YYYY') datenum('2010','YYYY')])
        ylim([40 90])
        params1 = {['rate = ', num2str(round(preQ8.deltaeval.m*365.25)/1000), ' \pm', num2str(round(preQ8.deltaeval.err*365.25,-1)/1000), ' km/yr'], ...
            ['r^2 = ', num2str(round(preQ8.deltaeval.r2, 2))]};
        format1 = sprintf('%s\n', params1{:});
        annot1 = text(0.3, 0.1, format1(1:end-1), ...
            'Color', [0 0 0], 'Parent', s1, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean delta radius')
        xlabel('year')
        ylabel('mean delta radius (km)')
        box on
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    s3 = subplot(1, 2, 2);
    manu.tab.qingint(1) = 0;
    manu.tab.qingint(3) = NaN;
        cla
        hold on
%             plot(repmat(datenum('1976','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
%             plot(repmat(datenum('1997','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            e3 = errorbar([datenum('1855','YYYY'); manu.tab.date; auto.tab.date(slct)], [0; manu.tab.qingint; auto.tab.qingint(slct)] / 1000, ...
                [2000; manu.tab.qinginterr; auto.tab.qinginterr(slct)] / 1000, 'LineStyle', 'none', 'Color', barColor);
            m = plot(manu.tab.date, manu.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            a = plot(auto.tab.date(slct), auto.tab.qingint(slct) / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            f = plot(qing.tab.date, qing.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
            plot(qing.inteval.xs, qing.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1970', 'YYYY') datenum('2000', 'YYYY')])
        set(gca, 'XTick', datenum({'1975', '1980', '1985', '1990', '1995', '2000'}, 'YYYY'), 'XTickLabel', {'1975', '', '1985', '', '', '2000'})
        ylim([0 40])
        params3 = {['rate = ', num2str(round(qing.inteval.m*365.25,-1)/1000), ' \pm', num2str(round(qing.inteval.err*365.25,-1)/1000), ' km/yr'], ...
            ['r^2 = ', num2str(round(qing.inteval.r2, 2))]};
        format3 = sprintf('%s\n', params3{:});
        annot3 = text(0.35, 0.1, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
        box on
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(fig,'Visible', 'on');
    set(fig, 'Pos', [100 100 800 400]);
    set(fig, 'PaperPositionMode', 'auto')
    print('-dpng', '-r300', './figs/JGR_shoreline_data.png');
    print('-depsc','-r300','-painters', './figs/JGR_shoreline_data.eps');
end


function [fig] = all_data(manu, auto, qing, preQ8, Q8, mod, fig)
    figure(fig)
    [colorOrder] = get(gca,'ColorOrder');
    
    % mean delta radius
    s1 = subplot(2, 3, 1);
        cla
        hold on
            plot(repmat(datenum('1854','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0])
            m = plot(manu.tab.date, manu.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :));
            a = plot(auto.tab.date, auto.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :));
            plot(preQ8.deltaeval.xs, preQ8.deltaeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        legend([m a], {'manual', 'auto'}, 'Location', 'NorthWest')
        datetick('x', 'YYYY')
        xlim([datenum('1850','YYYY') datenum('2020','YYYY')])
        % model parameters writing
        params1 = {['rate = ', num2str(round(preQ8.deltaeval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(preQ8.deltaeval.r2, 2))]};
        format1 = sprintf('%s\n', params1{:});
        annot1 = text(0.5, 0.1, format1(1:end-1), ...
            'Color', [0 0 0], 'Parent', s1, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean delta radius')
        xlabel('year')
        ylabel('mean delta radius (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
        
    % qinshuigou lobe radius
    s2 = subplot(2, 3, 2);
        cla
        hold on
            plot(repmat(datenum('1976','YYYY'), 1, 2), [0 30], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 30], ':', 'Color', [0 0 0])
            plot(manu.tab.date, manu.tab.meanloberad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
            plot(auto.tab.date, auto.tab.meanloberad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
            plot(qing.radeval.xs, qing.radeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1950','YYYY') datenum('2020','YYYY')])
        % model parameters writing
        params2 = {['rate = ', num2str(round(qing.radeval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(qing.radeval.r2, 2))]};
        format2 = sprintf('%s\n', params2{:});
        annot2 = text(0.5, 0.1, format2(1:end-1), ...
            'Color', [0 0 0], 'Parent', s2, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean Qingshuigou lobe radius')
        xlabel('year')
        ylabel('mean lobe radius (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
        
    % qingshuigou lobe lengths
    s3 = subplot(2, 3, 3);
        cla
        hold on
            plot(repmat(datenum('1976','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(repmat(max(mod.tab.date), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(manu.tab.date, manu.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
            plot(auto.tab.date, auto.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
            plot(qing.inteval.xs, qing.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
            plot(qing.retreateval.xs, qing.retreateval.ys / 1000, 'LineStyle', ':', 'LineWidth', 2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1950','YYYY') datenum('2020','YYYY')])
        % model (progradation) parameters writing
        params = {['rate = ', num2str(round(qing.inteval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(qing.inteval.r2, 2))]};
        format3 = sprintf('%s\n', params{:});
        annot3 = text(0.5, 0.1, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        % model (retreat) parameters writing
        params = {['rate = ', num2str(round(qing.retreateval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(qing.retreateval.r2, 2))]};
        format3 = sprintf('%s\n', params{:});
        annot3 = text(0.5, 0.25, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('Qingshuiguo lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
        
    % Q8 lobe length
    s4 = subplot(2, 3, 4);
        cla
        hold on
        plot(repmat(min(Q8.tab.date), 1, 2), [0 35], ':', 'Color', [0 0 0])
        plot(repmat(max(Q8.tab.date), 1, 2), [0 35], ':', 'Color', [0 0 0])
        plot(repmat(max(mod.tab.date), 1, 2), [0 35], ':', 'Color', [0 0 0])
        plot(manu.tab.date, manu.tab.Q8int / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
        plot(auto.tab.date, auto.tab.Q8int / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
        plot(Q8.inteval.xs, Q8.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        plot(Q8.retreateval.xs, Q8.retreateval.ys / 1000, 'LineStyle', ':', 'LineWidth', 2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1990','YYYY') datenum('2020','YYYY')])
        % model (progradation) parameters writing
        params = {['rate = ', num2str(round(Q8.inteval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(Q8.inteval.r2, 2))]};
        format3 = sprintf('%s\n', params{:});
        annot3 = text(0.5, 0.1, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s4, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        % model (retreat) parameters writing
        params = {['rate = ', num2str(round(Q8.retreateval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(Q8.retreateval.r2, 2))]};
        format3 = sprintf('%s\n', params{:});
        annot3 = text(0.5, 0.3, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s4, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('Q8 lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
        
    % modern lobe lengths
    s5 = subplot(2, 3, 5);
        cla
        hold on
        plot(repmat(max(Q8.tab.date), 1, 2), [0 35], ':', 'Color', [0 0 0])
        plot(repmat(max(mod.tab.date), 1, 2), [0 35], ':', 'Color', [0 0 0])
        plot(manu.tab.date, manu.tab.modint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
        plot(auto.tab.date, auto.tab.modint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
        plot(mod.inteval.xs, mod.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1990','YYYY') datenum('2020','YYYY')])
        params = {['rate = ', num2str(round(mod.inteval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(mod.inteval.r2, 2))]};
        format3 = sprintf('%s\n', params{:});
        annot3 = text(0.5, 0.1, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s5, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('modern lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
        box on
        set(gca, 'LineWidth', 1.1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
        
    set(fig,'Visible', 'on');
    set(fig, 'Pos', [100 100 1200 800]);
    set(fig, 'PaperPositionMode', 'auto')
    
    print('-dpng', '-r300', './all_shoreline_data.png');
    print('-depsc', '-painters', '-r300', './all_shoreline_data.eps');
    
    set(fig, 'PaperUnits', 'Points', 'PaperSize', [fig.Position(3), fig.Position(4)])
    print('-dpdf', '-r300', './all_shoreline_data.pdf');
    
end


function [fig] = XOM_meanrad(manu, auto, lobe, modelset, fig)
    figure(fig)
    [colorOrder] = get(gca, 'ColorOrder');
    barColor = [0.8 0.8 0.8];
    s1 = subplot(1, 2, 1);
        cla
        hold on
            plot(repmat(datenum('1854','YYYY'), 1, 2), [0 100], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 100], ':', 'Color', [0 0 0])
            e1 = errorbar([manu.tab.date; auto.tab.date], [manu.tab.meandeltarad; auto.tab.meandeltarad] / 1000, ...
                [manu.tab.meandeltaraderr; auto.tab.meandeltaraderr] / 1000, 'LineStyle', 'none', 'Color', barColor);
            m = plot(manu.tab.date, manu.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            a = plot(auto.tab.date, auto.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            f = plot(modelset.tab.date, modelset.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
            plot(modelset.deltaeval.xs, modelset.deltaeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1850','YYYY') datenum('2020','YYYY')])
        params1 = {['rate = ', num2str(round(modelset.deltaeval.m*365.25)/1000), ' \pm', num2str(round(modelset.deltaeval.err*365.25,-1)/1000), ' km/yr'], ...
            ['r^2 = ', num2str(round(modelset.deltaeval.r2, 2))]};
        format1 = sprintf('%s\n', params1{:});
        annot1 = text(0.3, 0.1, format1(1:end-1), ...
            'Color', [0 0 0], 'Parent', s1, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean delta radius')
        xlabel('year')
        ylabel('mean delta radius (km)')
        box on
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    s3 = subplot(1, 2, 2);
    manu.tab.qingint(1) = 0;
    manu.tab.qingint(3) = NaN;
        cla
        hold on
            plot(repmat(datenum('1976','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            e3 = errorbar([manu.tab.date; auto.tab.date], [manu.tab.qingint; auto.tab.qingint] / 1000, ...
                [manu.tab.qinginterr; auto.tab.qinginterr] / 1000, 'LineStyle', 'none', 'Color', barColor);
            m = plot(manu.tab.date, manu.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            a = plot(auto.tab.date, auto.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            f = plot(lobe.tab.date, lobe.tab.qingint / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
            plot(lobe.inteval.xs, lobe.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1960', 'YYYY') datenum('2020', 'YYYY')])
        set(gca, 'XTick', [datenum('1975', 'YYYY') datenum('2000', 'YYYY')], 'XTickLabel', {'1975', '2000'})
        ylim([ 0 35])
        params3 = {['rate = ', num2str(round(lobe.inteval.m*365.25,-1)/1000), ' \pm', num2str(round(lobe.inteval.err*365.25,-1)/1000), ' km/yr'], ...
            ['r^2 = ', num2str(round(lobe.inteval.r2, 2))]};
        format3 = sprintf('%s\n', params3{:});
        annot3 = text(0.35, 0.1, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
        box on
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(fig,'Visible', 'on');
    set(fig, 'Pos', [100 100 800 400]);
    set(fig, 'PaperPositionMode', 'auto')
    print('-dpng', '-r300', './figs/XOM_shoreline_data.png');
    print('-depsc','-r300','-painters', './figs/XOM_shoreline_data.eps');
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

