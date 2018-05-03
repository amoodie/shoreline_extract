function make_shoreline()
    
    % load data
    [manu.data] = load_data('manual');
    [auto.data] = load_data('auto');
    T = load('./shorelines/deltaLLcoord.mat');
    deltaCROPLLcoord = T.deltaLLcoord;
    deltaLLcoord = [611066, 4149751]; % this is Lijin
    
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
    [manu.qingint] = get_intersection(manu.lobe, 'AC1'); % intersection with Qingshuigou
    [auto.qingint] = get_intersection(auto.lobe, 'AC1'); 
    [manu.Q8int] = get_intersection(manu.lobe, 'AC2'); % intersection with Q8 lobe
    [auto.Q8int] = get_intersection(auto.lobe, 'AC2');
    [manu.modint] = get_intersection(manu.lobe, 'modern'); % intersection with modern lobe
    [auto.modint] = get_intersection(auto.lobe, 'modern');
    
    % make tables
    manu.tab = make_table(manu);
    auto.tab = make_table(auto);
    all.tab = make_rangetable(manu, auto, [datenum('1850','YYYY'), datenum('2020','YYYY')]); % all data from all time
    qing.tab = make_rangetable(manu, auto, [datenum('1976','YYYY'), datenum('07/31/1997','mm/dd/YYYY')]); % only data from qingshuigou lobe development
    postqing.tab = make_rangetable(manu, auto, [datenum('07/31/1997','mm/dd/YYYY'), datenum('2020', 'YYYY')]); % only data from *after* qingshuigou lobe development
    preQ8.tab = make_rangetable(manu, auto, [datenum('1850','YYYY'), datenum('1997','YYYY')]); % all data before abandonment of qingshuigou
    Q8.tab = make_rangetable(manu, auto, [datenum('07/31/1997','mm/dd/YYYY'), datenum('06/31/2007','mm/dd/YYYY')]); % only data from Q8 lobe development
    mod.tab = make_rangetable(manu, auto, [datenum('06/31/2007','mm/dd/YYYY'), datenum('2020','YYYY')]); % only data from modern lobe development
    
    % make models
    % entire delta radius model
    preQ8.deltamodel = fitlm(preQ8.tab, 'meandeltarad ~ date');
    preQ8.deltaeval.b = preQ8.deltamodel.Coefficients.Estimate(1);
    preQ8.deltaeval.m = preQ8.deltamodel.Coefficients.Estimate(2);
    preQ8.deltaeval.bserr = bootstrp(35, @mean, preQ8.deltamodel.Residuals.Raw);
    preQ8.deltaeval.CI = coefCI(preQ8.deltamodel);
    preQ8.deltaeval.err = mean(abs(preQ8.deltaeval.CI(2,:)-preQ8.deltaeval.m));
    preQ8.deltaeval.r2 = preQ8.deltamodel.Rsquared.ordinary;
    preQ8.deltaeval.xs = linspace(min(preQ8.tab.date), max(preQ8.tab.date), 10);
    preQ8.deltaeval.ys = ((preQ8.deltaeval.m .* preQ8.deltaeval.xs) + preQ8.deltaeval.b);
    
    % quingshuigou lobe radius model
    qing.radmodel = fitlm(qing.tab, 'meanloberad ~ date');
    qing.radeval.b = qing.radmodel.Coefficients.Estimate(1);
    qing.radeval.m = qing.radmodel.Coefficients.Estimate(2);
    qing.radeval.CI = coefCI(qing.radmodel);
    qing.radeval.err = mean(abs(qing.radeval.CI(2,:)-qing.radeval.m));
    qing.radeval.r2 = qing.radmodel.Rsquared.ordinary;
    qing.radeval.xs = linspace(min(qing.tab.date), max(qing.tab.date), 10);
    qing.radeval.ys = ((qing.radeval.m .* qing.radeval.xs) + qing.radeval.b);
    
    % quingshuigou lobe progradation rate model
    qing.intmodel = fitlm(qing.tab, 'qingint ~ date');
    qing.inteval.b = qing.intmodel.Coefficients.Estimate(1);
    qing.inteval.m = qing.intmodel.Coefficients.Estimate(2);
    qing.inteval.CI = coefCI(qing.intmodel);
    qing.inteval.err = mean(abs(qing.inteval.CI(2,:)-qing.inteval.m));
    qing.inteval.r2 = qing.intmodel.Rsquared.ordinary;
    qing.inteval.xs = linspace(min(qing.tab.date), max(qing.tab.date), 10);
    qing.inteval.ys = ((qing.inteval.m .* qing.inteval.xs) + qing.inteval.b);
    
    % qingshuigou lobe retreat after abandonment model
    qing.retreatmodel = fitlm(postqing.tab, 'qingint ~ date');
    qing.retreateval.b = qing.retreatmodel.Coefficients.Estimate(1);
    qing.retreateval.m = qing.retreatmodel.Coefficients.Estimate(2);
    qing.retreateval.CI = coefCI(qing.retreatmodel);
    qing.retreateval.err = mean(abs(qing.retreateval.CI(2,:)-qing.retreateval.m));
    qing.retreateval.r2 = qing.retreatmodel.Rsquared.ordinary;
    qing.retreateval.xs = linspace(min(postqing.tab.date), max(postqing.tab.date), 10);
    qing.retreateval.ys = ((qing.retreateval.m .* qing.retreateval.xs) + qing.retreateval.b);
    
    % Q8 lobe progradation model
    Q8.intmodel = fitlm(Q8.tab, 'Q8int ~ date');
    Q8.inteval.b = Q8.intmodel.Coefficients.Estimate(1);
    Q8.inteval.m = Q8.intmodel.Coefficients.Estimate(2);
    Q8.inteval.CI = coefCI(Q8.intmodel);
    Q8.inteval.err = mean(abs(Q8.inteval.CI(2,:)-Q8.inteval.m));
    Q8.inteval.r2 = Q8.intmodel.Rsquared.ordinary;
    Q8.inteval.xs = linspace(min(Q8.tab.date), max(Q8.tab.date), 10);
    Q8.inteval.ys = ((Q8.inteval.m .* Q8.inteval.xs) + Q8.inteval.b);
    
    % Q8 lobe retreat after abandonment model
    Q8.retreatmodel = fitlm(mod.tab, 'Q8int ~ date');
    Q8.retreateval.b = Q8.retreatmodel.Coefficients.Estimate(1);
    Q8.retreateval.m = Q8.retreatmodel.Coefficients.Estimate(2);
    Q8.retreateval.CI = coefCI(Q8.retreatmodel);
    Q8.retreateval.err = mean(abs(Q8.retreateval.CI(2,:)-Q8.retreateval.m));
    Q8.retreateval.r2 = Q8.retreatmodel.Rsquared.ordinary;
    Q8.retreateval.xs = linspace(min(mod.tab.date), max(mod.tab.date), 10);
    Q8.retreateval.ys = ((Q8.retreateval.m .* Q8.retreateval.xs) + Q8.retreateval.b);
    
    % modern lobe progradation model
    mod.intmodel = fitlm(mod.tab, 'modint ~ date');
    mod.inteval.b = mod.intmodel.Coefficients.Estimate(1);
    mod.inteval.m = mod.intmodel.Coefficients.Estimate(2);
    mod.inteval.CI = coefCI(mod.intmodel);
    mod.inteval.err = mean(abs(mod.inteval.CI(2,:)-mod.inteval.m));
    mod.inteval.r2 = mod.intmodel.Rsquared.ordinary;
    mod.inteval.xs = linspace(min(mod.tab.date), max(mod.tab.date), 10);
    mod.inteval.ys = ((mod.inteval.m .* mod.inteval.xs) + mod.inteval.b);
    
    % make plots
    fig_alldata = figure('Visible', 'off');
    fig_shorelines = figure('Visible', 'off');
    fig_XOMmeanrad = figure('Visible', 'off');
    fig_JGRmeanrad = figure('Visible', 'off');
    [fig_alldata] = all_data(manu, auto, qing, preQ8, Q8, mod, fig_alldata);
%     [fig_shorelines] = all_shorelines(manu, auto, deltaLLcoord, deltaCROPLLcoord, lobe_mask, fig_shorelines);
%     [fig_XOMmeanrad] = XOM_meanrad(manu, auto, lobe, preAC2, fig_XOMmeanrad);
    [fig_JGRmeanrad] = JGR_meanrad(manu, auto, qing, preQ8, fig_JGRmeanrad);
    
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
        cell2mat(data.loberad(:,3)), cell2mat(data.loberad(:,4)), cell2mat(data.qingint(:,3)), cell2mat(data.qingint(:,4)), ...
        cell2mat(data.Q8int(:,3)), cell2mat(data.Q8int(:,4)), cell2mat(data.modint(:,3)), cell2mat(data.modint(:,4)) );
    table = array2table(datamat, ...
        'VariableNames', {'date', 'meandeltarad', 'meandeltaraderr', 'meanloberad', 'meanloberaderr', ...
        'qingint', 'qinginterr', 'Q8int', 'Q8interr', 'modint', 'modinterr'});
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

%% stable
function [cmap] = colormap_fun(n, s)

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

