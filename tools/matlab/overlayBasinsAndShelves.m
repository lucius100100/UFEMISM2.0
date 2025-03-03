%function to overlay basins and shelves
function overlayBasinsAndShelves(ax, basinFile, shelvesFile)

    %basins
    basinShp = shaperead(basinFile);

    %subregion names
    subregionNames = {basinShp.Subregions};
    uniqueSubs     = unique(subregionNames);

    for s = 1:length(uniqueSubs)
        thisSubregion = uniqueSubs{s};
        if isempty(thisSubregion)
            continue; 
        end
        
        idx = strcmp(subregionNames, thisSubregion);

        %union polygon
        unionPoly = polyshape();
        isInit    = false;
        for i = find(idx)
            x = basinShp(i).X;
            y = basinShp(i).Y;
            if ~isempty(x) && ~isempty(y)
                try
                    p = polyshape(x, y, 'Simplify', true);
                    if ~isInit
                        unionPoly = p;
                        isInit = true;
                    else
                        unionPoly = union(unionPoly, p);
                    end
                catch ME
                    warning("Error building union for '%s': %s", thisSubregion, ME.message);
                end
            end
        end
        
        %plotting merged subregion
        if isInit
            plot(ax, unionPoly, ...
                'FaceColor',[0.7, 0.7, 0.7], ...   
                'EdgeColor','k', ...            
                'LineWidth',1, ...
                'FaceAlpha',1.0, ...          
                'HandleVisibility','off');
            
            %label subregions
            [cx, cy] = centroid(unionPoly);
            text(ax, cx, cy, thisSubregion, ...
                'Color','k', 'FontSize',8, 'FontWeight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end

        %plotting basins
        for i = find(idx)
            x = basinShp(i).X;
            y = basinShp(i).Y;
            plot(ax, x, y, ...
                'Color',[0.3, 0.3, 0.3], ... 
                'LineWidth',0.5, ...
                'HandleVisibility','off');
        end
    end

    %shelves
    shelfShp = shaperead(shelvesFile);

    for i = 1:length(shelfShp)
        x = shelfShp(i).X;
        y = shelfShp(i).Y;
        if ~isempty(x) && ~isempty(y)
            try
                pShelf = polyshape(x, y, 'Simplify', true);
                plot(ax, pShelf, ...
                    'FaceColor',[0.8, 0.8, 1], ... 
                    'FaceAlpha',0.3, ...           
                    'EdgeColor','k', ...
                    'LineStyle','--', ...
                    'LineWidth',1, ...
                    'DisplayName','Ice shelf');    
            catch ME
                warning("Error with ice shelf %d: %s", i, ME.message);
            end
        end
    end
    
end
