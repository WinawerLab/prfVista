function rootPath = prfRootPath()
% Determine path to root of the PRFmodel directory
%
%        rootPath = prfRootPath;
%
% This function MUST reside in the directory at the base of the
% PRFmodel directory structure 
%
rootPath = which('prfRootPath');

rootPath = fileparts(rootPath);

return
