function Q = runstd(A)

% RUNSTD Computes standard curves for each isotype/bead size/dilution
%
%    RUNSTD(X) takes as input a structure array created by RUNEXPT and 
%    extracts the structure array elements derived from standard samples in
%    order to create a standard curve by fitting to a 4-parameter logistic.
%    This function generates plots for calibration curves (one plot per
%    flow configuration/dilution) for each isotype/bead size.  It outputs a
%    data structure array (one structure per flow configuration/dilution)
%    containing the fit parameters for each isotype/bead size combination
%    (curve fits are computed separately for different bead sizes even for
%    the same isotype secondary).
%    
%    Usage:
%
%        Q = RUNSTD(X);
%
%    Uses the function L4P from Matlab file exchange:
%    Cardillo G. (2012) Four parameters logistic regression
%    https://it.mathworks.com/matlabcentral/fileexchange/38122
%
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


    % extract the standard curve samples from the RUNEXPT structure array 
    C = A([A.iscontrol]' == 1);

    % numbers of control samples
    nc = length(C);

    % get info from ctrl data
    D = table;
    D.fconf = cell(nc, 1);
    for i = 1:nc, D.fconf(i) = cellstr(C(i).results.flowconf); end
    D.dil = [C.dilution]';
    D.conc = [C.concentration]';

    % unique flow configs to look at - this is how many structure elements
% NOTE: this probably should be flow config AND dilution factor because you
% could have identical flow configs but with different dilution factors
    ufc = unique(D.fconf);

    % prepare output structure
    Q = struct;

    % step through flow configs
    for i = 1:length(ufc)
    
        % this flow config
        Q(i).flow = string(ufc(i));
    
        % this dilution
        udc = unique(D.dil(strcmp(D.fconf, ufc(i))));   % needs to be fixed - there are no dilutions for stds
        if ~isscalar(udc), error('Multiple dilutions for flow config'); end    % does this make sense
        Q(i).dilution = udc;  
    
        % get indices of standard samples for this flow config
        v = intersect(find(D.fconf == Q(i).flow), find(D.dil == Q(i).dilution));  % take out dilution
    
        %get secondaries
        temp = strsplit(Q(i).flow, {':', ','});
        Q(i).secondaries = temp(1:2:end); 
    
        % set up calib curve table
        Q(i).curve = table;
        Q(i).curve.fcs = {C(v).fcs}';
        Q(i).curve.conc = D.conc(v);
    
        % step through concentrations
        for j = 1:length(v)
        
            % get data from corresponding standard sample
            idx = intersect(v, find(Q(i).curve.conc(j) == D.conc));
            tab = C(idx).results.out;
        
            % number bead sizes
            nb = length(unique(tab.size));
        
            % step through secondaries
            for k = 1:length(Q(i).secondaries)
            
                % this isotype             
                iso = Q(i).secondaries(k);

                % step through bead sizes
                for b = 1:nb

                    % header names
                    nm = iso + "_" + num2str(b);

                    % values to average
                    values = tab.(iso)(tab.name == iso & tab.size == b);
                    Q(i).curve.(nm)(j) = mean(values);

                end
            end
        end
    
        % sort concentration table by concentration
        [Q(i).curve, Q(i).isort] = sortrows(Q(i).curve, 'conc');

    end

    % color key
    I = table;
    I.isotype = ["IgM" "IgG" "IgA"]';
    I.color = {[0 .6 0]; [.8 0 0]; [.3 .4 .8]};

    % symbols for bead sizes
    sym = {'o', '^'};

    % plot layout
    pw = 450;
    ph = 400;
    left = 70;
    bot = 80;
    top = 50;
    sp = 60;
    right = 40;
    fw = left + (length(Q) * pw) + ((length(Q) - 1) * sp) + right;
    fh = bot + ph + top;
    xoff = (left + ((0:(length(Q) - 1)) * (pw + sp))) / fw;
    figure('position', [10 10 fw fh]);

    % draw plots
    ymax = 0;
    for i = 1:length(Q)

        % flow config title for plot
        subplot('position', [xoff(i), (bot + ph) / fh, pw / fw, top / fh]);
        text(.5, .75, ['Flow Configuration: ' char(Q(i).flow)],...
            'fontname', 'arial',...
            'fontsize', 18,...
            'horizontalalignment', 'center');
        text(.5, .25, ['Dilution: ' num2str(Q(i).dilution)],...
            'fontname', 'arial',...
            'fontsize', 18,...
            'horizontalalignment', 'center');
        set(gca,...
            'xlim', [0 1],...
            'ylim', [0 1],...
            'xcolor', 'w',...
            'ycolor', 'w',...
            'xtick', [],...
            'ytick', [],...
            'box', 'off');

        % data axes
        ax(i) = subplot('position', [xoff(i), bot / fh, pw / fw, ph / fh],...
            'nextplot', 'add');

        % isotype bead combos
        B = table;
%         iso = repmat(Q(i).secondaries, 2, 1); % <- this was hardcoded???
        iso = repmat(Q(i).secondaries, nb, 1);
        B.beadsize = repmat((1:nb)', length(Q(i).secondaries), 1);
        B.beadsizelabel = repmat("", length(B.beadsize), 1);
        for q = 1:length(B.beadsize)
            w = find(A(1).results.beadsize.size == B.beadsize(q));
            B.beadsizelabel(q) = A(1).results.beadsize.label(w);
        end
        bead = repmat(strcat("_", string(num2str((1:nb)'))), 1, length(Q(i).secondaries));
        isobead = iso + bead;
        B.iso = iso(:);
        B.isobead = isobead(:);
        B.legend = B.beadsizelabel + " " + B.iso;

        % add fit param columns
        B.A = nan(length(B.iso), 1);
        B.B = nan(length(B.iso), 1);
        B.C = nan(length(B.iso), 1);
        B.D = nan(length(B.iso), 1);

        % plot handles
        h = nan(length(B.iso), 1);

        % step though isotype/bead combos compute fits, make plots
        for j = 1:length(B.iso)

            % plot standards data points
            x = log10(Q(i).curve.conc);
            y = Q(i).curve.(B.isobead(j));
            c = I.color{I.isotype == B.iso(j)};
            h(j) = plot(x, y, sym{B.beadsize(j)}, 'color', c, 'markerfacecolor', c);

            % populate fit struct
            L = struct;
            L.isotype = B.iso(j);
            L.beadsize = B.beadsize(j);
            L.fit = L4P(x, y);
            B.A(j) = L.fit.A;
            B.B(j) = L.fit.B;
            B.C(j) = L.fit.C;
            B.D(j) = L.fit.D;

            % fit values
            yf = L.fit.D + (L.fit.A - L.fit.D) ./ (1 + (x / L.fit.C).^L.fit.B);        
            Q(i).lp4(j) = L;

            % plot fit curve
            plot(x, yf, '-', 'color', c, 'linewidth', 1); 

        end

        % new fit parameter table
        Q(i).fitparam = removevars(B, {'isobead', 'legend'});
        Q(i).fitparam.Properties.VariableNames(3) = {'isotype'};

        % determine x-range for plot
        xrange = [...
            min(log10(Q(i).curve.conc)),...
            max(log10(Q(i).curve.conc))];
        
        % plot legend
        hleg = legend(h, B.legend, 'location', 'northwest');
        set(hleg, 'fontname', 'arial', 'fontsize', 18, 'box', 'off');

        % ensure same ylim for all plots
        if ymax < max(ylim), ymax = max(ylim); end

        % xlabel
        xlab = strip(cellstr(num2str(Q(i).curve.conc)));
        set(gca,...
            'xtick', x,...
            'xticklabel', xlab,...
            'xticklabelrotation', -60,...
            'xlim', xrange,...
            'tickdir', 'out',...
            'linewidth', 1,...
            'fontname', 'arial',...
            'fontsize', 12);
        xlabel('Arbitrary Concentration Units', 'fontname', 'arial', 'fontsize', 18);

    end

    % axes preferences
    set(ax, 'ylim', [0 ymax]);
    hy = ylabel(ax(1), 'arcsinh MFI', 'fontname', 'arial', 'fontsize', 18);
    hy.Position(1) = min(xlim) - (diff(xlim) * .07);

    % clean up plot
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');
    
return

