function runqc(A, cutoff)

% RUNQC Generates visualization to assess gating consistency, bead counts
%
%    RUNQC(X, threshold) takes as input a structure array created by
%    RUNEXPT and generates QC figures including plots showing gate
%    positions per sample (one plot for each bead size) and plots showing
%    the events per peak and singlets per bead size for each sample.
%
%    *** NOTE *** this is currently hardcoded for assays w/ two bead sizes
% I THOUGHT I FIXED THIS
%    
%    Usage:
%
%        RUNQC(X, 30);
%
%    Requires Matlab Statistics Toolbox
%
% -------------------------------------------------------------------------
%
%    Authors: Alexander F. Rosenberg (afr@uab.edu) and Rodney G. King
%        with John T. Killian, Todd J. Green, J. Akther, M. Emon Hossain,
%        Shihong Qiu, Guang Yang, Troy D. Randall and Frances E. Lund
%
%    University of Alabama at Birmingham
%    Department of Microbiology
%    April 11, 2023
%    Copyright (C) 2023 UAB Research Foundation
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


    % number of samples processed by RUNEXPT
    n = length(A);
    
    % number of bead sizes
    ns = length(unique(A(1).results.out.size));
    
    % number of axes (gate plots for each size and then two more axes)
    na = ns + 2;

    % plot dimension settings    
    bot = 65;
    pwg = 300;
    ph = 850;
    sp = 30;
    top = 10;
    left = 80;
    right = 20;
    xrange = [-1 9];
    
    % colors for different bead sizes (may need to add more)
    col = {[1 0 0], [0 0 1], [0 .6 0]};
    
    % figure dimensions    
    fw = left + (pwg * na) + (sp * (na - 1)) + right;
    fh = bot + ph + top;

    % render figure window    
    figure('position', [10 10 fw fh]);
    
    % axes positions
    xoff = (left + ((0:(na - 1)) * (pwg + sp))) / fw;
    axd = [bot / fh, pwg / fw, ph / fh];
    ax = nan(na, 1);
    for i = 1:ns
        ax(i, 1) = subplot('position', [xoff(i) axd], 'nextplot', 'add');
    end
    ax(ns + 1, 1) = subplot('position', [xoff(ns + 1) axd], 'nextplot', 'add');
    ax(ns + 2, 1) = subplot('position', [xoff(ns + 2) axd], 'nextplot', 'add');

    % initialize event count variables    
    total = nan(n, 1);
    B = struct;
    for i = 1:ns
        B(i).nbead = nan(n, 1);
        B(i).nsing = nan(n, 1);
    end

    % iterate through samples (elements in input stucture array)    
    for i = 1:n    

        % step through bead sizes
        for s = 1:ns
            
            % rows of results output table
            k = A(i).results.out.size == s;
            
            % get coords of peaks
            x = [A(i).results.out.left(k) A(i).results.out.right(k)];
            
            % event count for bead peaks
            e = A(i).results.out.events(k);
            
            % plot bead gate peaks for each size
            plot(ax(s), x', [i i], 'color', col{s}, 'linewidth', 1);
            
            % plot events per peak
            plot(ax(ns + 1), e, i, '.', 'color', col{s});
                    
            % annotate per-sample, per-peak event counts below threshold
            if ~isempty(find(e < cutoff, 1))
                plot(ax(ns + 1), e(e < cutoff), i, 'o',...
                    'color', col{s},...
                    'markersize', 12,...
                    'linewidth', 2);
            end
        
            % get counts for current sample                    
            total(i)      = A(i).results.counts.total_events;
            B(s).nbead(i) = A(i).results.counts.bead_size_events(s);            
            B(s).nsing(i) = A(i).results.counts.bead_size_singlets(s);
            
        end            

    end
    
    % plot # singlets per sample for each bead size
    for s = 1:ns
        plot(ax(ns + 2), B(s).nsing, 1:n, '-', 'color', col{s}, 'linewidth', 1.5);
    end

    % show plate divisions on plot (assume samples sorted by plate)
    platediv = find(diff(grp2idx({A.plate}'))) + .5;
    for s = 1:ns
        plot(ax(s), xrange, (platediv * [1 1])', '-', 'color', [1 1 1] * .3, 'linewidth', .5);
    end

    % clean up plots    
    set(ax(1:ns),...
        'xlim', xrange,...
        'ylim', [0 n + 1],...
        'linewidth', 1,...
        'tickdir', 'out',...
        'fontname', 'arial',...
        'fontsize', 12);
    if ns > 1, set(ax(ns), 'yticklabel', {}); end      
    set(ax(ns + 1),...
        'xlim', [0 max(get(ax(ns + 1), 'xlim'))],...
        'ylim', [0 n + 1],...
        'yticklabel', {},...
        'linewidth', 1,...
        'tickdir', 'out',...
        'fontname', 'arial',...
        'fontsize', 12);
    set(ax(ns + 2),...
        'xlim', [0 max(get(ax(ns + 2), 'xlim'))],...
        'ylim', [0 n + 1],...
        'yticklabel', {},...
        'linewidth', 1,...
        'tickdir', 'out',...
        'fontname', 'arial',...
        'fontsize', 12);

    % show plate divisions on plot (assume samples sorted by plate)
    plot(ax(ns + 1), get(ax(ns + 1), 'xlim'), (platediv * [1 1])', '-', 'color', [1 1 1] * .3, 'linewidth', .5);
    plot(ax(ns + 2), get(ax(ns + 2), 'xlim'), (platediv * [1 1])', '-', 'color', [1 1 1] * .3, 'linewidth', .5);
    
    %plot labels
    hy = ylabel(ax(1), 'Sample in Experiment', 'fontname', 'arial', 'fontsize', 18);
    hy.Position(1) = -2.8;
    
    for s = 1:ns
        hx = xlabel(ax(s), ['Gates: arcsinh MFI (Bead Size ' num2str(s) ')'], 'fontname', 'arial', 'fontsize', 18);
        hx.Position(2) = -1 * (n * .04);
    end
    hx = xlabel(ax(ns + 1), 'Events per Peak', 'fontname', 'arial', 'fontsize', 18);
    hx.Position(2) = -1 * (n * .04);
    hx = xlabel(ax(ns + 2), 'Singlets per Bead Size', 'fontname', 'arial', 'fontsize', 18);
    hx.Position(2) = -1 * (n * .04);

    % make figure background white    
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');



