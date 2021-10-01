% Gneneric Address
% If pc replace / with \ otherwise replaces \ with /
function out = gen_addr(x)
if ispc
    out = strrep(x, '/', '\');
elseif ismac
    out = strrep(x, '\', '/');
end