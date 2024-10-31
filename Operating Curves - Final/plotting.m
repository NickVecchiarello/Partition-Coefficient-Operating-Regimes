function y = plotting(input, b, Load_top, Load_bottom)
    LW = 4;                                                                 %Line Width
    limy = [0, 2*ceil(max(Load_top)/2)];
    limx = [input.min_Vr_range-0.05, input.max_Vr_range];
    linec = [0.3, 0.3, 0.3];
    shadec = [0.8, 0.8, 0.8];
    xlab = 'Phase Ratio (Vl/Vr)';
    ylab = 'Loading (mg/mL)';
    xtickinterval = 5;
    ytickinterval = 1;
    fs = 24;
    width = 1200;
    height = 800;
    figure
    hold on
    
    %% Determines if Top and Bottom Operating Curves Intersect
    [~, mi] = find((Load_top-Load_bottom <0 ));                             %Finds top and bottom corss point if it exists    
    if isempty(mi) ~= 1
       b = b(max(mi):end);
       Load_top = Load_top(max(mi):end);
       Load_bottom = Load_bottom(max(mi):end);
    else 
    end
   
    %% Plot top operating curve
    %plot(b, smooth(smooth(Load_top)), '--', 'LineWidth', 3, 'Color', [0.35, 0.35, 0.35])
    plot(b, smooth(smooth(Load_top)), 'LineWidth', LW, 'Color', linec)
    ylim(limy);
    xlim(limx);

    %% Plot bottom operating line
    %plot(b, smooth(Load_bottom), '--', 'LineWidth', 3, 'Color', [0.35, 0.35, 0.35])
    plot(b, smooth(Load_bottom), 'LineWidth', LW, 'Color', linec)

    %% Shade Operating Region
    b2 = [b, fliplr(b)];
    inBetween = [smooth(smooth(Load_top)); flip(smooth(smooth(Load_bottom)))];
    fill(b2, inBetween, shadec);
    xlabel(xlab);
    ylabel(ylab);
    Pix_SS = get(0,'screensize');
    set(gca,'FontSize',fs);
    set(gcf, 'position', [(Pix_SS(3)-width)/2, (Pix_SS(4)-height)/2, width, height]);
    xticks([linspace(floor(input.min_Vr_range), input.max_Vr_range, (input.max_Vr_range)/xtickinterval +1)]);
    yticks([linspace(0, max(limy), max(limy)/ytickinterval +1)]);
    h = gca; % Get axis to modify
    h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
    h.XAxis.MinorTickValues;
    h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
    h.YAxis.MinorTickValues;
    box on
end
