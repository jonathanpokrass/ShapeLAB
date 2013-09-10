function [] = log(varargin)
global IS_LOGGING_ENABLED
if IS_LOGGING_ENABLED
    disp([datestr(now,13) ' - ' sprintf(varargin{:})])
end
end
