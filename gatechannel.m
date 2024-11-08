function [x, y, loc, pk, xgate] = gatechannel(ch, xrange, binsize, sm, minpeakh, minpeakp)

% GATECHANNEL find peaks in a transformed FCS file channel
%
%    [x, y, loc, pk, gate] = GATECHANNEL(c, xr, bs, sm, mph, mpp) processes
%    a channel of transformed cytometry MFI data and identifies peaks in
%    the corresponding histogram.  The arguments are:
%        1. transformed (e.g. arcsinh) cytometry channel MFI data
%        2. the range of MFI (2-element vector)
%        3. bin size for constructing a histogram of MFI data
%        4. factor for loess-smoothing of histogram data
%        5. minimum acceptable peak height
%        6. minimum acceptable peak prominence
%    The return values are:
%        1. the x and y coordinates for the smoothed MFI histogram
%        2. the location of the peak (n-element vector for n peaks)
%        3. x-coordinates of peaks (n-element vector)
%        4. y-coordinates of peaks (n-element vector)
%        5. boundaries for gates (nx2 matrix for n peaks)
%    The boundaries of a peak are determined by identifying the y
%    coordinate at the midpoint of a peak's prominence, and noting where a
%    horizontal line positioned at that y coordinate intersects the sides
%    of the peak.
%
%    Usage:
%
%        [x, y, loc, pk g] = GATECHANNEL(c, [1 8], .05, .03, .1, .1);
%
%    Requires Matlab Signal Processing Toolbox and Curve Fitting Toolbox
%
% -------------------------------------------------------------------------
%
%    Authors: Alexander F. Rosenberg (afr@uab.edu) and Rodney G. King
%        with John T. Killian, Fen Zhou, Davide Botta, Todd J. Green,
%        Jobaida Akther, M. Emon Hossain, Shihong Qiu, Guang Yang,
%        Troy D. Randall and Frances E. Lund
%
%    University of Alabama at Birmingham
%    Department of Microbiology
%    April 11, 2023
%    Copyright (C) 2024 UAB Research Foundation
%    This software is offered with no guarantees of any kind.
%
%    see: "A high-throughput multiplex array for antigen-specific serology
%    with automated analysis", bioRxiv 2023 April.
%    doi: 10.1101/2023.03.29.534777
%
%    This file (part of the "CBA Toolbox") is free software: you can
%    redistribute it and/or modify it under the terms of the GNU General
%    Public License as published by the Free Software Foundation, version 3
%    of the License.  This file is distributed in the hope that it will be
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    General Public License for more details.  You should have received a
%    copy of the GNU General Public License along with this program.  If  
%    not, see https://www.gnu.org/licenses/gpl-3.0.en.html.


    % histogram of (transformed) channel data
    binedge = xrange(1):binsize:xrange(2);
    x = (binedge(2:end) - (binsize / 2))';
    y = histcounts(ch, binedge);
    y = y ./ max(y);
    y = smooth(x, y, sm, 'loess');

    % find peaks
    [pk, loc, ~, pr] = findpeaks(y, x,...
        'MinPeakHeight', minpeakh,...
        'MinPeakProminence', minpeakp,...
        'annotate', 'extents',...
        'WidthReference', 'halfprom');

    % number of peaks
    n = length(pk);
        
    % initialize transformed gate coords
    xgate = nan(n, 2);

    % compute peak width based on prominences
    for i = 1:n

        % midpoint for prominence
        ym = (pk(i) - pr(i)) + (pr(i) / 2);

        % distance to left edge of peak
        ileft = max(intersect(find(y < ym), find(x < loc(i))));
        m = (y(ileft + 1) - y(ileft)) / (x(ileft + 1) - x(ileft));
        b = y(ileft) - (m * x(ileft));
        xleft = (ym - b) / m;    

        % distance to right edge of peak
        iright = min(intersect(find(y < ym), find(x > loc(i))));
        m = (y(iright) - y(iright - 1)) / (x(iright) - x(iright - 1));
        b = y(iright) - (m * x(iright));
        xright = (ym - b) / m;     

        % gate
        xgate(i, :) = [xleft xright];

    end
              
end
