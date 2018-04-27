function make_shoreline()
    
    % load data
    [manu.data] = load_data('manual');
    [auto.data] = load_data('auto');
    T = load('./shorelines/deltaLLcoord.mat');
    deltaCROPLLcoord = T.deltaLLcoord;
    deltaLLcoord = [611066 4149751]; % this is Lijin
    
    % make subsets
    [lobe_mask] = make_mask([6.77e5 4.177e6]);
    [manu.lobe] = get_subsets(manu.data, lobe_mask);
    [auto.lobe] = get_subsets(auto.data, lobe_mask);
    
    % maths
    obsN = size(manu, 1) + size(auto, 1); % number of total observations
    [manu.deltarad] = get_radius(manu.data, deltaLLcoord);
    [auto.deltarad] = get_radius(auto.data, deltaLLcoord);
    [manu.loberad] = get_radius(manu.lobe, lobe_mask.acoord);
    [auto.loberad] = get_radius(auto.lobe, lobe_mask.acoord);
    [manu.lobeint] = get_intersection(manu.lobe, 'AC1');
    [auto.lobeint] = get_intersection(auto.lobe, 'AC1');
    [manu.AC2int] = get_intersection(manu.lobe, 'AC2');
    [auto.AC2int] = get_intersection(auto.lobe, 'AC2');
    [manu.modint] = get_intersection(manu.lobe, 'modern');
    [auto.modint] = get_intersection(auto.lobe, 'modern');
    
    % make tables
    manu.tab = make_table(manu);
    auto.tab = make_table(auto);
    lobe.tab = make_rangetable(manu, auto, [datenum('1976','YYYY') datenum('07/31/1997','mm/dd/YYYY')]);
    all.tab = make_rangetable(manu, auto, [datenum('1850','YYYY') datenum('2020','YYYY')]);
    preAC2.tab = make_rangetable(manu, auto, [datenum('1850','YYYY') datenum('1997','YYYY')]);
    
    % make models
    preAC2.deltamodel = fitlm(preAC2.tab, 'meandeltarad ~ date');
    preAC2.deltaeval.b = preAC2.deltamodel.Coefficients.Estimate(1);
    preAC2.deltaeval.m = preAC2.deltamodel.Coefficients.Estimate(2);
    
    preAC2.deltaeval.bserr = bootstrp(35, @mean, preAC2.deltamodel.Residuals.Raw);
    preAC2.deltaeval.CI = coefCI(preAC2.deltamodel);
    preAC2.deltaeval.err = mean(abs(preAC2.deltaeval.CI(2,:)-preAC2.deltaeval.m));
    preAC2.deltaeval.r2 = preAC2.deltamodel.Rsquared.ordinary;
    preAC2.deltaeval.xs = linspace(min(preAC2.tab.date), max(preAC2.tab.date), 10);
    preAC2.deltaeval.ys = ((preAC2.deltaeval.m .* preAC2.deltaeval.xs) + preAC2.deltaeval.b);
    
    lobe.radmodel = fitlm(lobe.tab, 'meanloberad ~ date');
    lobe.radeval.b = lobe.radmodel.Coefficients.Estimate(1);
    lobe.radeval.m = lobe.radmodel.Coefficients.Estimate(2);
    lobe.radeval.CI = coefCI(lobe.radmodel);
    lobe.radeval.err = mean(abs(lobe.radeval.CI(2,:)-lobe.radeval.m));
    lobe.radeval.r2 = lobe.radmodel.Rsquared.ordinary;
    lobe.radeval.xs = linspace(min(lobe.tab.date), max(lobe.tab.date), 10);
    lobe.radeval.ys = ((lobe.radeval.m .* lobe.radeval.xs) + lobe.radeval.b);
    
    lobe.intmodel = fitlm(lobe.tab, 'lobeint ~ date');
    lobe.inteval.b = lobe.intmodel.Coefficients.Estimate(1);
    lobe.inteval.m = lobe.intmodel.Coefficients.Estimate(2);
    lobe.inteval.CI = coefCI(lobe.intmodel);
    lobe.inteval.err = mean(abs(lobe.inteval.CI(2,:)-lobe.inteval.m));
    lobe.inteval.r2 = lobe.intmodel.Rsquared.ordinary;
    lobe.inteval.xs = linspace(min(lobe.tab.date), max(lobe.tab.date), 10);
    lobe.inteval.ys = ((lobe.inteval.m .* lobe.inteval.xs) + lobe.inteval.b);
    
    % convert to lines for publication
%     manu.data(:,3) = manu.data(:,2);
%     for i = 1:size(auto.data,1)
%         input = auto.data{i,2};
%         [auto.data(i, 3)] = get_ordered(input);
%     end
    
    % make plots
    fig_alldata = figure('Visible', 'off');
    fig_shorelines = figure('Visible', 'off');
    fig_XOMmeanrad = figure('Visible', 'off');
    fig_JGRmeanrad = figure('Visible', 'off');
    [fig_alldata] = all_data(manu, auto, lobe, preAC2, fig_alldata);
%     [fig_shorelines] = all_shorelines(manu, auto, deltaLLcoord, deltaCROPLLcoord, lobe_mask, fig_shorelines);
%     [fig_XOMmeanrad] = XOM_meanrad(manu, auto, lobe, preAC2, fig_XOMmeanrad);
    [fig_JGRmeanrad] = JGR_meanrad(manu, auto, lobe, preAC2, fig_JGRmeanrad);
    
    [fig_movie] = movie_fig(manu, auto, deltaLLcoord, deltaCROPLLcoord, lobe_mask, fig_shorelines);
    
    save('./shorelines/shoreline_out.mat', 'all')
    save('./shorelines/meandeltarad_out.mat', 'preAC2')
end

function [rad] = get_radius(data, deltaLLcoord)
    Nobs = size(data, 1);
    rad = cell(Nobs, 4);
    for i = 1:Nobs
        pts = data{i, 2};
        Npts = size(pts, 1);
        if cell2mat(data(i, 1)) < datenum('0/01/1972'); err = 10000; else err = 2000; end
        for j = 1:Npts
            xdiff = pts(j, 1) - deltaLLcoord(1);
            ydiff = pts(j, 2) - deltaLLcoord(2);
            radobs(j) = sqrt(xdiff^2 + ydiff^2);
        end
        if ~exist('radobs','var')
            radobs = NaN;
        end
        rad(i, 1) = num2cell(data{i, 1});
        rad(i, 2) = {radobs};
        rad(i, 3) = {mean(radobs)};
%         rad(i, 4) = {sqrt(  sum( repmat(err/sqrt(Npts), 1, Npts))  )}; % error on the mean
        rad(i, 4) = {err};
%         try
%             rad(i, 4) = {mean(bootstrp(35, @mean, radobs))};
%         catch
%             rad(i, 4) = {NaN};
%         end
    end
end

function [intersection] = get_intersection(data, set)
    switch set
        case 'AC1'
            % define the line at some grid interval
            filename = '../../maps/intersection_channelline.csv';
            channelline = csvread(filename, 1, 0);
        case 'AC2'
            filename = '../../maps/intersection_AC2.csv';
            channelline = csvread(filename, 1, 0);
        case 'modern'
            filename = '../../maps/intersection_2016.csv';
            channelline = csvread(filename, 1, 0);
    end
    betweenintline = [0 sqrt((channelline(2:end, 1) - channelline(1:end-1, 1)).^2 + (channelline(2:end, 2) - channelline(1:end-1, 2)).^2)'];
    alongintline = cumsum(betweenintline);
    Nobs = size(data, 1);
    intersection = cell(size(data, 1), 4);
    for i = 1:Nobs
        if cell2mat(data(i, 1)) < datenum('0/01/1972'); err = 10000; else err = 2000; end
        if ~isempty(data{i, 2})
            clear intline
            % reduce the channelline data to only those within x AND y range of the shoreline
            keepx = and(channelline(:, 1) <= max(data{i, 2}(:, 1)), channelline(:, 1) >= min(data{i, 2}(:, 1)));
            keepy = and(channelline(:, 2) <= max(data{i, 2}(:, 2)), channelline(:, 2) >= min(data{i, 2}(:, 2)));
            keep = and(keepx, keepy);
            intline(:, 1) = channelline(keep, 1); % cut the dataset down to only below upperlim
            intline(:, 2) = channelline(keep, 2); % cut the dataset down
            % reduce the shorline data to only those within x AND y range of the new intline
            % NEED TO WRITE THIS CODE!
            Nintpts = size(intline, 1);
            mindists = NaN(Nintpts, 1);
            for j = 1:Nintpts
                distances = sqrt(((data{i, 2}(:, 1) - intline(j, 1)).^2 + (data{i, 2}(:, 2) - intline(j, 2)).^2));
                [mindists(j), ~] = min(distances);
            end
            if ~isempty(intline)
                [~, intidx] = min(mindists);
                intcoords = intline(intidx, :);
                intdist = alongintline(intidx);
            else
                intcoords = [NaN NaN];
                intdist = NaN;
            end
        else
            intcoords = [NaN NaN];
            intdist = NaN;
        end
        intersection(i, 1) = num2cell(data{i, 1});
        intersection(i, 2) = {intcoords};
        intersection(i, 3) = {intdist};
        intersection(i, 4) = {sqrt(err^2 + err^2)};
    end
end

function [table] = make_table(data)
    datamat = horzcat( cell2mat(data.data(:,1)), cell2mat(data.deltarad(:,3)), cell2mat(data.deltarad(:,4)), ...
        cell2mat(data.loberad(:,3)), cell2mat(data.loberad(:,4)), cell2mat(data.lobeint(:,3)), cell2mat(data.lobeint(:,4)), ...
        cell2mat(data.AC2int(:,3)), cell2mat(data.AC2int(:,4)), cell2mat(data.modint(:,3)), cell2mat(data.modint(:,4)) );
    table = array2table(datamat, ...
        'VariableNames', {'date', 'meandeltarad', 'meandeltaraderr', 'meanloberad', 'meanloberaderr', ...
        'lobeint', 'lobeinterr', 'AC2int', 'AC2interr', 'modint', 'modinterr'});
end

function [table] = make_rangetable(manu, auto, range)
    keep.manu = and(manu.tab.date >= range(1), manu.tab.date <= range(2));
    keep.auto = and(auto.tab.date >= range(1), auto.tab.date <= range(2));
    table = vertcat(manu.tab(keep.manu, :), auto.tab(keep.auto, :));
end

function [mask] = make_mask(acoord)
    mask.acoord = acoord; % apex coord
    mask.thet = 60; % opening angle
    mask.xdist = 1100 * (30); % max distance from apex (don't change 30)
    mask.maxh = [mask.acoord(1)+mask.xdist mask.acoord(2)+tand(mask.thet)*mask.xdist];
    mask.maxl = [mask.acoord(1)+mask.xdist mask.acoord(2)-tand(mask.thet)*mask.xdist];
end

function [subset] = get_subsets(data, mask)
    Nobs = size(data, 1);
    subset = cell(size(data));
    for i = 1:Nobs
        clear pts
        pts(:,2) = data{i, 2}((data{i, 2}(:,2) <= mask.maxh(2)), 2); % cut the dataset down to only below upperlim
        pts(:,1) = data{i, 2}((data{i, 2}(:,2) <= mask.maxh(2)), 1); % cut the dataset down
        Npts = size(pts, 1);
        keeps = zeros(Npts, 1);
        for j = 1:Npts
            xdiff = pts(j, 1) - mask.acoord(1);
            ydiff = pts(j, 2) - mask.acoord(2);
            hlim = [mask.acoord(1)+xdiff mask.acoord(2)+tand(mask.thet)*xdiff];
            llim = [mask.acoord(1)+xdiff mask.acoord(2)-tand(mask.thet)*xdiff];
            if pts(j, 2) < hlim(2) && pts(j, 2) > llim(2)
                keeps(j) = 1;
            end
        end
        subset(i, 1) = num2cell(data{i, 1});
        subset(i, 2) = {pts(logical(keeps), :)};
    end
end

function [line] = get_ordered(pointlist)
    used = false(size(pointlist, 1), 1);
    cnt = 1;
    last = cnt;
    used(last) = true;
    line(cnt, :) = pointlist(cnt, :);
    lastdist = 0;
    while lastdist <= 60;
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
    line = {line(1:end-1, :)};
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

function [fig] = JGR_meanrad(manu, auto, lobe, modelset, fig)
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
            f = plot(modelset.tab.date, modelset.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
            plot(modelset.deltaeval.xs, modelset.deltaeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1840','YYYY') datenum('2010','YYYY')])
        ylim([40 90])
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
    manu.tab.lobeint(1) = 0;
    manu.tab.lobeint(3) = NaN;
        cla
        hold on
%             plot(repmat(datenum('1976','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
%             plot(repmat(datenum('1997','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            e3 = errorbar([datenum('1855','YYYY'); manu.tab.date; auto.tab.date(slct)], [0; manu.tab.lobeint; auto.tab.lobeint(slct)] / 1000, ...
                [2000; manu.tab.lobeinterr; auto.tab.lobeinterr(slct)] / 1000, 'LineStyle', 'none', 'Color', barColor);
            m = plot(manu.tab.date, manu.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            a = plot(auto.tab.date(slct), auto.tab.lobeint(slct) / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            f = plot(lobe.tab.date, lobe.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
            plot(lobe.inteval.xs, lobe.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1970', 'YYYY') datenum('2000', 'YYYY')])
        set(gca, 'XTick', datenum({'1975', '1980', '1985', '1990', '1995', '2000'}, 'YYYY'), 'XTickLabel', {'1975', '', '1985', '', '', '2000'})
        ylim([0 40])
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
    print('-dpng', '-r300', './figs/JGR_shoreline_data.png');
    print('-depsc','-r300','-painters', './figs/JGR_shoreline_data.eps');
end



function [fig] = all_data(manu, auto, lobe, modelset, fig)
    figure(fig)
    [colorOrder] = get(gca,'ColorOrder');
    s1 = subplot(2, 3, 1);
        cla
        hold on
            plot(repmat(datenum('1854','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 70], ':', 'Color', [0 0 0])
            m = plot(manu.tab.date, manu.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :));
            a = plot(auto.tab.date, auto.tab.meandeltarad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :));
            plot(modelset.deltaeval.xs, modelset.deltaeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        legend([m a], {'manual', 'auto'}, 'Location', 'NorthWest')
        datetick('x', 'YYYY')
        xlim([datenum('1850','YYYY') datenum('2020','YYYY')])
        params1 = {['rate = ', num2str(round(modelset.deltaeval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(modelset.deltaeval.r2, 2))]};
        format1 = sprintf('%s\n', params1{:});
        annot1 = text(0.5, 0.1, format1(1:end-1), ...
            'Color', [0 0 0], 'Parent', s1, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean delta radius')
        xlabel('year')
        ylabel('mean delta radius (km)')
    s2 = subplot(2, 3, 2);
        cla
        hold on
            plot(repmat(datenum('1976','YYYY'), 1, 2), [0 30], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 30], ':', 'Color', [0 0 0])
            plot(manu.tab.date, manu.tab.meanloberad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
            plot(auto.tab.date, auto.tab.meanloberad / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
            plot(lobe.radeval.xs, lobe.radeval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1950','YYYY') datenum('2020','YYYY')])
        params2 = {['rate = ', num2str(round(lobe.radeval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(lobe.radeval.r2, 2))]};
        format2 = sprintf('%s\n', params2{:});
        annot2 = text(0.5, 0.1, format2(1:end-1), ...
            'Color', [0 0 0], 'Parent', s2, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('mean lobe radius')
        xlabel('year')
        ylabel('mean lobe radius (km)')
    s3 = subplot(2, 3, 3);
        cla
        hold on
            plot(repmat(datenum('1976','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(manu.tab.date, manu.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
            plot(auto.tab.date, auto.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
            plot(lobe.inteval.xs, lobe.inteval.ys / 1000, 'LineStyle', '--', 'LineWidth', 1.2, 'Color', [0 0 0])
        datetick('x', 'YYYY')
        xlim([datenum('1950','YYYY') datenum('2020','YYYY')])
        params = {['rate = ', num2str(round(lobe.inteval.m*365.25)), ' m/yr'], ...
            ['r^2 = ', num2str(round(lobe.inteval.r2, 2))]};
        format3 = sprintf('%s\n', params{:});
        annot3 = text(0.5, 0.1, format3(1:end-1), ...
            'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('Qingshuiguo lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
    s4 = subplot(2, 3, 4);
        cla
        hold on
        plot(manu.tab.date, manu.tab.AC2int / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
        plot(auto.tab.date, auto.tab.AC2int / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
        datetick('x', 'YYYY')
        xlim([datenum('1990','YYYY') datenum('2020','YYYY')])
%         params = {['rate = ', num2str(round(lobe.inteval.m*365.25)), ' m/yr'], ...
%             ['r^2 = ', num2str(round(lobe.inteval.r2, 2))]};
%         format3 = sprintf('%s\n', params{:});
%         annot3 = text(0.5, 0.1, format3(1:end-1), ...
%             'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('AC2 lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
    s5 = subplot(2, 3, 5);
        cla
        hold on
        plot(manu.tab.date, manu.tab.modint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(1, :))
        plot(auto.tab.date, auto.tab.modint / 1000, 'o', 'MarkerSize', 4, 'Color', colorOrder(2, :))
        datetick('x', 'YYYY')
        xlim([datenum('1990','YYYY') datenum('2020','YYYY')])
%         params = {['rate = ', num2str(round(lobe.inteval.m*365.25)), ' m/yr'], ...
%             ['r^2 = ', num2str(round(lobe.inteval.r2, 2))]};
%         format3 = sprintf('%s\n', params{:});
%         annot3 = text(0.5, 0.1, format3(1:end-1), ...
%             'Color', [0 0 0], 'Parent', s3, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        title('modern lobe length')
        xlabel('year')
        ylabel('intersection distance from datum (km)')
        
    set(fig,'Visible', 'on');
    set(fig, 'Pos', [100 100 1200 800]);
    set(fig, 'PaperPositionMode', 'auto')
%     print('-dpng', '-r300', './figs/shoreline_data.png');
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
    manu.tab.lobeint(1) = 0;
    manu.tab.lobeint(3) = NaN;
        cla
        hold on
            plot(repmat(datenum('1976','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            plot(repmat(datenum('1997','YYYY'), 1, 2), [0 35], ':', 'Color', [0 0 0])
            e3 = errorbar([manu.tab.date; auto.tab.date], [manu.tab.lobeint; auto.tab.lobeint] / 1000, ...
                [manu.tab.lobeinterr; auto.tab.lobeinterr] / 1000, 'LineStyle', 'none', 'Color', barColor);
            m = plot(manu.tab.date, manu.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            a = plot(auto.tab.date, auto.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'Color', [0.3 0.3 0.3]);
            f = plot(lobe.tab.date, lobe.tab.lobeint / 1000, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r');
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

%% stable
function [cmap] = colormap_fun(n, s)
%     RRR = linspace(0  ,255,n+1)/255;
%     GGG = linspace(0  ,0  ,n+1)/255;
%     BBB = linspace(255,0  ,n+1)/255;
    tfig = figure('Visible', 'off');
    [colorOrder] = get(gca, 'ColorOrder');
    ends = colorOrder(s:s+1, :) * 255;
    close(tfig)

    RRR = linspace(ends(1, 1), ends(2, 1), n+1) / 255;
    GGG = linspace(ends(1, 2), ends(2, 2), n+1) / 255;
    BBB = linspace(ends(1, 3), ends(2, 3), n+1) / 255;
    
    [cmap] = colormap([RRR',GGG',BBB']);
end

function [data] = load_data(srcstr)
    folders = strsplit(ls(strcat('./shorelines/', srcstr, '_shoreline/')));
    countfolders = length(folders)-1;
    data = cell(countfolders, 2);
    for i = 1:countfolders
        filename = strcat('./shorelines/', srcstr, '_shoreline/', char(folders(i)));
        points = csvread(filename, 1, 0);
        namesplit = strsplit(char(folders(i)), {'_', '.'});
        date = datenum(namesplit(2));
        data(i, 1) = num2cell(date);
        data(i, 2) = {points(:,[1 2])};
    end
    dates = arrayfun(@cell2mat, data(:,1));
    [~, idx] = sort(dates);
    data = data(idx, :);
end

