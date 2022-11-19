function Dq = Dq(tau,q)

 % Dq - generalized fractal dimention
 %
 % tau is the scaling exponent. q is the index variable.

 %   [Dq = Dq(tau,q) returns the spetrum of Generalized fractal Dimention
 %   tau(hoeld(q),q), hoeld(q) and q.
 %   tau is a vector or 2D or 3D matrix. If tau is a matrix,
 %   spect returns the spetrum from each column of the matrix.
 %   See also mdfa,hoelder,tau
 %
 % References
 % 1.  Muzy, J., Bacry, E., Arneodo, A.: Multifractal formalism revisited with wavelets. International Journal of Bifurcation and Chaos in Applied Sciences and Engineering 04 (April 1994)

if nargin == 1, q = [-2.^(4:-1:-3),0,2.^(-3:4)]; end%(2)
sizx = size(tau);%(2)
Dq = tau./(repmat(q',[1,sizx(2:end)]) - 1);