function output_txt = data_cursor_cross(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
    % Create the data cursor text string:
    pos = get(event_obj,'Position');
    output_txt = {['X: ',num2str(pos(1),4)],...
                  ['Y: ',num2str(pos(2),4)]};
    % If there is a Z-coordinate in the position, display it as well:
    if length(pos) > 2
      output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
    end
    % Find the hggroup object that goes with the event:
    posArgs = {'XData',pos(1),'YData',pos(2)};
    if length(pos) > 2
      posArgs(5:6) = {'ZData',pos(3)};
    end
    hMarker = findall(0,'Tag','DataTipMarker','Type','line',posArgs{:});
    hGroup = get(hMarker,'Parent');
    % If hGroup is a cell array, it means more than one data cursor is at the
    %   given position. Since there is no clear way to distinguish them, the
    %   crosshairs for both should be initialized/updated:
    if iscell(hGroup)
      cellfun(@update_crosshairs,hGroup);
    else
      update_crosshairs(hGroup);
    end
%--- Begin nested functions -----------------------------------------------
    function update_crosshairs(hObject)
      if isappdata(hObject,'DataCursorCrosshairs')
        % Update the crosshairs:
        hCrosshairs = getappdata(hObject,'DataCursorCrosshairs');
        if length(hCrosshairs) > 2
          set(hCrosshairs(1),'YData',pos([2 2]),'ZData',pos([3 3]));
          set(hCrosshairs(2),'XData',pos([1 1]),'ZData',pos([3 3]));
          set(hCrosshairs(3),'XData',pos([1 1]),'YData',pos([2 2]));
        else
          set(hCrosshairs(1),'YData',pos([2 2]));
          set(hCrosshairs(2),'XData',pos([1 1]));
        end
      else
        % Get the axes limits, create the crosshairs, set the crosshair
        %   handles as application data for the hggroup object, and set the
        %   DeleteFcn of the hggroup object so it removes the crosshairs:
        hAxes = get(hObject,'Parent');
        xLimits = xlim(hAxes);
        yLimits = ylim(hAxes);
        if length(pos) > 2
          zLimits = zlim(hAxes);
          hX = line('XData',xLimits,'YData',pos([2 2]),'ZData',pos([3 3]));
          hY = line('XData',pos([1 1]),'YData',yLimits,'ZData',pos([3 3]));
          hZ = line('XData',pos([1 1]),'YData',pos([2 2]),'ZData',zLimits);
          hCrosshairs = [hX hY hZ];
        else
          hX = line('XData',xLimits,'YData',pos([2 2]));
          hY = line('XData',pos([1 1]),'YData',yLimits);
          hCrosshairs = [hX hY];
        end
        setappdata(hObject,'DataCursorCrosshairs',hCrosshairs);
        set(hObject,'DeleteFcn',...
            @(src,event) delete(hCrosshairs(ishandle(hCrosshairs))));
      end
    end
end