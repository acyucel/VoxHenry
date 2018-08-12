% Octave work-around for 'fcnchk' MatLab function
% Reference: http://octave.1599824.n4.nabble.com/rootfinding-newtzero-td4631423.html
function [f,msg]=fcnchk_octave(x, n)
    f = x;
    msg.message = '';
    msg.identifier = '';
    msg = msg(zeros(0,1));
end
