function F = cbapeak(F)

% CBAPEAK find peaks for each bead size bucket
%
%    A = CBAPEAK(A) reads a data structure generated by CBASINGLET and
%    identifies peaks associated with beads in each bead-size series.
%
%    Usage:
%
%        A = CBAPEAK(A);
%
%    Requires Matlab Signal Processing Toolbox
%    Requires Matlab Curve Fitting Toolbox
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


    % structure for each bead size
    C = struct;
    
    % init error message
    errmsg = '';
    
    % augment data structure to annotate events as belonging to peak
    F.event.peak = nan(size(F.event, 1), 1);
    
    % get expected number peaks and number of bins per bead size (sorted)
    B = sortrows(F.beadconf, {'size', 'bin'});
    np = nan(F.nbeadsize, 1);
    nb = nan(F.nbeadsize, 1);
    for i = 1:F.nbeadsize
        k = unique(B.npeaks(B.size == i));
        if length(k) ~= 1
            error(['split beads ' num2str(i) 'have different numbers of features']);
        end
        np(i) = k;
        nb(i) = length(find(B.size == i));
    end
    
    % step through bead sizes
    for i = 1:F.nbeadsize
        
        % gates per bead size
        [C(i).x, C(i).y, C(i).loc, C(i).pk, C(i).gate] = gatechannel(...
            F.event.bead(F.event.size == i),...
            F.settings.beadpeak_range,...
            F.settings.beadpeak_bin_size,...
            F.settings.beadpeak_loess_smooth,...
            F.settings.beadpeak_min_peak_height,...
            F.settings.beadpeak_min_peak_prominence);
        
        % if more peaks than expected, accept top n (assume others noise)
        [~, isort] = sort(C(i).pk, 'descend');
        rnk = (1:length(C(i).pk))';
        rnk(isort) = rnk;
        ikeep = find(rnk <= np(i));        
        C(i).loc = C(i).loc(ikeep);
        C(i).pk = C(i).pk(ikeep);
        C(i).gate = C(i).gate(ikeep, :);
        
        % test for expected number of peaks
        if np(i) ~= size(C(i).gate, 1)
            errmsg = [errmsg ' bead size ' num2str(i)];
        end
        
        % get associated ag names per peak order for current bead size
        C(i).beadid = ones(1, np(i)) * i;
        C(i).peakid = 1:np(i);
        C(i).events = nan(1, np(i));
        
    end
    
    % if render plot
    if F.display == 1
        
        % pixel dimensions of components
        pw = 600;
        ph = 120;
        vsp = 35;
        bot = 60;
        left = 60;
        right = 10;
        top = 10;
        fw = left + pw + right;
        fh = bot + (ph * F.nbeadsize) + (vsp * (F.nbeadsize - 1)) + top;
        yoff = (((0:(F.nbeadsize - 1)) * (ph + vsp)) + bot) / fh;
                
        % draw figure
        FIG4 = figure('position', [25 10 fw fh]);   
    
        % peak plot for each bead size
        for k = 1:F.nbeadsize        
            subplot('position', [left / fw, yoff(k), pw / fw, ph / fh],...
                'nextplot', 'add');
            plot(C(k).x, C(k).y, '-',...
                'linewidth', 1.5,...
                'color', F.pref.colororder{k});
            for i = 1:size(C(k).gate, 1)            
                plot(C(k).gate(i, :), [1 1] * C(k).pk(i) + .03, '-',...
                    'linewidth', 2,...
                    'color', F.pref.colororder{k});        
                text(C(k).loc(i), C(k).pk(i) + .1, num2str(i),...
                    'fontname', 'arial',...
                    'fontsize', 12,...
                    'color', F.pref.colororder{k},...
                    'horizontalalignment', 'center');
            end
            set(gca,...
                'xlim', F.settings.beadpeak_range,...
                'ylim', [0 1.2],...
                'box', 'off',...
                'linewidth', 1,...
                'tickdir', 'out',...
                'fontname', 'arial',...
                'fontsize', 12);
            if k == 1
                xlabel('Bead Separation (arcsinh-xform)',...
                    'fontname', 'arial',...
                    'fontsize', 18);
            end
        end
    
        % make figure background white
        set(gcf, 'inverthardcopy', 'off', 'color', 'w');     

        % add to figure handles
        F.figures(4) = FIG4;        

    end
    
    % throw error for incorrect # peaks (after figure rendered)    
    if ~isempty(errmsg), error(['unexpected number of peaks:' errmsg]); end  
    
    % bin sizes for output table (run length encoded)
    a = cumsum([1 B.npeaks']);
    b = zeros(1, sum(B.npeaks));
    b(a(1:end - 1)) = 1;
    bins = B.bin(cumsum(b));
    
    % gates for output table
    gt = cell2mat({C(B.size).gate}');
    
    % is there a split channel
    if isnan(F.splitchan), split = 0; else, split = 1; end
    
    % prepare output table
    O = table;
    O.size  = [C(B.size).beadid]';
    O.bin   = bins(:);
    O.peak  = [C(B.size).peakid]';
    O.left  = gt(:, 1);
    O.right = gt(:, 2);
    
    % step through bead sizes, bins, peaks; assign relative peak orders
    E = struct;   
    for i = 1:F.nbeadbin        
        sz = B.size(i);        
        if split == 1
            bn = B.bin(i);
            sgate = F.counts.bead_size_split_gates(sz);
            if isnan(sgate)
                po = ['peakorder' num2str(sz)];
            else
                po = ['peakorder' num2str(sz) 'p' num2str(bn)];
            end
        else
            po = ['peakorder' num2str(sz)];
        end
        E(i).feat = strsplit(F.(po), ',');        
        for j = 1:B.npeaks(i)
            g = F.event.bead >= C(sz).gate(j, 1) & ...
                F.event.bead <= C(sz).gate(j, 2) & ...
                F.event.size == sz;            
            if split == 1 && ~isnan(sgate), g = g & F.event.bin == bn; end
            F.event.peak(g) = j;                        
            E(i).events(1, j) = sum(g);                
        end
    end        
    O.name = string([E.feat]');
    O.events = [E.events]';
    
    % if no split channel, remove "bin" field for backwards compatibility
    if split == 0, O = removevars(O, 'bin'); end
    
    % assign output structure
    F.out = O;
    
return
