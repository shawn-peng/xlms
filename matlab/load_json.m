function data = load_json(jsonfile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(jsonfile); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data = jsondecode(str); % Using the jsondecode function to parse JSON from string

end

