function cycle_indices = findCycle(u, fs, initPos, ncycle) %#ok<INUSL> fs kept for interface compatibility
%FINDCYCLE Indices spanning NCYCLE cycles (fractional allowed) of a signal.
%   cycle_indices = findCycle(u, fs, initPos, ncycle) locates the
%   zero-crossing nearest the relative position INITPOS (in (0,1)) and
%   returns sample indices covering NCYCLE cycles from that point.
%   NCYCLE may be fractional (e.g. 1.5 for one and a half cycles);
%   half-cycles are resolved exactly at zero-crossings, finer fractions
%   are interpolated within the last half-cycle.

% ---- Orientation and input checks -----------------------------------
u = u(:)';                                   % force row vector
if nargin < 4 || isempty(ncycle), ncycle = 1; end
validateattributes(initPos, {'numeric'}, {'scalar', '>', 0, '<', 1})
validateattributes(ncycle,  {'numeric'}, {'scalar', '>=', 0.5})

% ---- Zero crossings ---------------------------------------------------
zeroIdx = find(u(1:end-1) .* u(2:end) < 0);  % sign changes
nz = numel(zeroIdx);

% Required crossing intervals: 2 per cycle -> may be fractional
nHalf     = 2 * ncycle;          % e.g. 1.5 cycles -> 3 half-cycles
nHalfFull = floor(nHalf);        % whole half-cycles at crossings
frac      = nHalf - nHalfFull;   % remainder inside the next half-cycle

nNeeded = nHalfFull + (frac > 0);   % crossing intervals we must fit
if nz < nNeeded + 1
    error("findCycle:notEnoughCycles", ...
        "Signal has %d zero-crossings; need at least %d for %.2f cycle(s).", ...
        nz, nNeeded + 1, ncycle);
end

% ---- Starting crossing near relative position initPos ----------------
k = round(nz * initPos);
k = max(k, 1);
k = min(k, nz - nNeeded);        % leave room for the requested span

% ---- End index: whole half-cycles + fractional remainder --------------
idxStart = zeroIdx(k);
idxEnd   = zeroIdx(k + nHalfFull);
if frac > 0
    % advance a fraction of the next half-cycle (in samples)
    idxNext = zeroIdx(k + nHalfFull + 1);
    idxEnd  = idxEnd + round(frac * (idxNext - idxEnd));
end

cycle_indices = idxStart:idxEnd;
end