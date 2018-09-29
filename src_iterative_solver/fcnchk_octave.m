% Octave work-around for 'fcnchk' MatLab function
% Reference: http://octave.1599824.n4.nabble.com/rootfinding-newtzero-td4631423.html
function [f,msg]=fcnchk_octave(x, n)
    f = x;
    % 'msg' contains any error message from fcnchk
    % so if 'x' is actually a matrix, let's put it in a message
    if isa(x,'float')
      msg.message = 'matrix';
      msg.identifier = '';
    else
      msg.message = '';
      msg.identifier = '';
      msg = msg(zeros(0,1));
    end

end
