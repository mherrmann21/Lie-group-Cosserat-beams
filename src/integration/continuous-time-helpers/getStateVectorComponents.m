function [R, x, eta] = getStateVectorComponents(xState)
    % Get the components of the full state vector, i.e.:
    %  * the array of node rotation matrix elements
    %  * the array of node position vectors
    %  * the array of node velocity vectors
    %
    % The state vector has the form
    % [ vec(R0), ..., vec(Rn), x0, ..., xn, eta0, ..., etan ]
    %

    % Get number of discrete nodes ( = Nr. of segments +1)
    nNodes = length(xState) / 18;

    R   = xState( 1 : 9*nNodes );
    x   = xState( 9*nNodes+1 : 12*nNodes );
    eta = xState( end-6*nNodes+1 : end );
end