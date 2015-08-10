classdef PlotInfo
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StartX
        EndX
        StateId
        LineStyle
        Color
    end
    
    methods
        function obj = PlotInfo(sx, ex, sid, linetype)
            obj.StartX = sx;
            obj.EndX = ex;
            obj.StateId = sid;
            obj.LineStyle = PickLineStyle(linetype);
            obj.Color = PickColor(sid);
        end
        
        function obj = MixState(info, sid)
            c = PickColor(sid);
            obj = info;
            obj.Color = (c + obj.Color)/2;
        end
    end
end


function f = PickLineStyle(lineType)
    if lineType == 3
        f = ':';
    elseif lineType == 2
        f = '-.';
    else
        f = '-';
    end
end

function c = PickColor(n)
    index = mod(n, 10);
    if index == 0
        index = 10;
    end
    if index == 1
        c = [1 0 0]; % 'r' red
    elseif index == 2
        c = [0 0 0]; % 'k' black
    elseif index == 3
        c =[0 1 0]; %  'g' green
    elseif index == 4
        c = [0 0 1]; % 'b' blue
    elseif index == 5
        c = [1 0 1]; % 'm' megenta
    elseif index == 6
        c = [0 1 1]; % 'c' cyan
    elseif index == 7
        c = [0.7 0.7 0.7];
    elseif index == 8
        c = [0.7 0.4 1];
    elseif index == 9
        c = [0.6 0 0.3];
    else
        c = [1 0.5 0];
    end
end