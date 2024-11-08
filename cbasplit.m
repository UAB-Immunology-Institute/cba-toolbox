function F = cbasplit(F)

% CBASPLIT read cba flow structure and annotate beads by split-bin
%
%    A = CBASPLIT(A) reads a data structure generated by CBASINGLET and
%    annotates events based on their corresponding bead size split bins, if
%    a split channel was specified. This function returns an amended
%    version of the input structure.  If figures were requested (set
%    display flag = 1 in cbaload which sets a field in the input structure)
%    then a plot showing the split gate(s) is displayed.
%
%    Usage:
%
%        A = CBASPLIT(A);
%
%    Requires Matlab Signal Processing Toolbox
%    Requires Matlab Curve Fitting Toolbox
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


    % if there is no split channel, we can ignore this function
    if ~isnan(F.splitchan)

        % structure for each bead size
        C = struct;

        % init error message
        errmsg = '';

        % augment data structure to annotate events as belonging to split bin
        F.event.bin = nan(size(F.event, 1), 1);

        % split gate for each bead size (stays as NaN if no split used)    
        gate = nan(F.nbeadsize, 1);

        % step through bead sizes
        for i = 1:F.nbeadsize

            % gates per bead size
            [C(i).x, C(i).y, C(i).loc, C(i).pk, C(i).gate] = gatechannel(...
                F.event.split(F.event.size == i),...
                F.settings.beadpeak_range,...
                F.settings.beadpeak_bin_size,...
                F.settings.beadpeak_loess_smooth,...
                F.settings.beadpeak_min_peak_height,...
                F.settings.beadpeak_min_peak_prominence);

            % if two peaks
            if size(C(i).gate, 1) > 1
                gate(i) = mean([C(i).gate(1, 2) C(i).gate(2, 1)]);
            end

            % test for expected number of peaks
            np = size(C(i).gate, 1);
            if isnan(gate(i)) && np ~= 1
                errmsg = [errmsg ' bead size ' num2str(i) ': expecting 1 peak, found ' num2str(np)];
            end
            if ~isnan(gate(i)) && np ~= 2
                errmsg = [errmsg ' bead size ' num2str(i) ': expecting 2 peaks, found ' num2str(np)];
            end

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
            FIG3 = figure('position', [25 10 fw fh]);   

            % peak plot for each bead size
            for k = 1:F.nbeadsize        
                subplot('position', [left / fw, yoff(k), pw / fw, ph / fh],...
                    'nextplot', 'add');
                plot(C(k).x, C(k).y, '-',...
                    'linewidth', 1.5,...
                    'color', F.pref.colororder{k, 1});
                if ~isnan(gate(k))
                    plot([1 1] * gate(k), [0 1.2], '--',...
                        'linewidth', 1,...
                        'color', F.pref.colororder{k, 1});
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
            F.figures(3) = FIG3;        

        end

        % throw error for incorrect # peaks (after figure rendered)    
        if ~isempty(errmsg), error(['unexpected number of peaks:' errmsg]); end  

        % step through bead sizes, assign events to bin(s)
        q = F.beadconf(:, {'size', 'bin'});
        q.events = nan(size(q, 1), 1);
        for i = 1:F.nbeadsize        
            if isnan(gate(i))
                idx = F.event.size == i;
                F.event.bin(idx) = 1;
                q.events(q.size == i) = length(find(idx));
            else
                idx1 = F.event.size == i & F.event.split < gate(i);
                idx2 = F.event.size == i & F.event.split > gate(i);
                F.event.bin(idx1) = 1;
                F.event.bin(idx2) = 2;
                q.events(q.size == i & q.bin == 1) = length(find(idx1));
                q.events(q.size == i & q.bin == 2) = length(find(idx2));
            end
        end
        
        % add counts and gate(s) to stats structure
        F.counts.bead_size_bins = q;
        F.counts.bead_size_split_gates = gate;
        
    end
    
return
