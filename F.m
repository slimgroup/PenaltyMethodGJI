function D = F(m,Q,model)
% Modeling operator, returns D = PA(m)^{-1}Q,
% where A(m) is the Helmholtz operator and P is a sampling operator.
%
%
% use:
%   D = F(m,Q,model)
%
% input:
%   mt - model [s^2/m^2]
%   Q  - sources
%   D  - data
%   model.h - [dz,dx] gridspacing in z and x direction [m]
%   model.n - [nz,nx] number of gridpoints in z and x direction
%   model.f - frequency
%   model.Is - source location indices
%   model.Ir - receiver location indices
%
% This program is part of the paper
% "Mitigating local minima in full-waveform inversion by expanding the search space",
% T. van Leeuwen and F.J. Herrmann, 2013 (submitted to GJI).
%
% Copyright (C) 2013 Tristan van Leeuwen (tleeuwen@eos.ubc.ca)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% get sampling operators
Ps = getP(model.n,model.Is{:});
Pr = getP(model.n,model.Ir{:});

% get Helmholtz matrix
A = getA(model.f,m(:),model.h,model.n);

% solve for wavefield and sample data
U = A\Ps'*Q;
D = Pr*U;