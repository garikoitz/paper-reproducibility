function rootPath = paperReprPath()
% Determine path to root of the code directory
%
%        rootPath = paperReprPath;
%
% This function MUST reside in the directory at the base of the code structure 
%
% Copyright Stanford team, GLU, 2019

rootPath = which('paperReprPath');

rootPath = fileparts(rootPath);

return
