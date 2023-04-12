function Y = runcalc(A, Q, expt, varargin)

% RUNCALC Computes standard curves for each isotype/bead size/dilution
%
%    RUNCALC(X, Q, ex) takes as input a structure array created by RUNEXPT
%    as the first argument, a standard curve structure array created by 
%    RUNSTD as the second argument, and an experiment name as the third 
%    argument and returns a structure array of tables (one per flow 
%    config/dilution) ofback-calculated concentrations.  Additionally, 
%    output tables are written to TXT file(s) - one per flow
%    configuration/dilution.
%    
%    RUNCALC(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%    parameter/value pair values:
%        'outdir'       specify a directory for the PDF and/or the TXT
%                       file to be written to - default is the directory
%                       from which this function is run
%
%    Usage:
%
%        Y = RUNCALC(A, Q, 'myexpt');
%
%        Y = RUNCALC(A, Q, 'myexpt', 'outdir', '/Users/afr/output/');
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


    % check for optional arg
    pa = inputParser;
    addParameter(pa, 'outdir',  './',...
        @(x) validateattributes(x, {'string', 'char'}, {'nonempty'}));
    parse(pa, varargin{:});    
    outdir = pa.Results.outdir;

    % flowconfigs in calib stuct
    flowconfigs = [Q.flow]';

    % number of samples
    n = length(A);

    % get flowconfigs of all samples
    fc = nan(n, 1);
    for i = 1:n, fc(i) = find(A(i).results.flowconf == flowconfigs); end

    % flag for control samples
    c = [A.iscontrol]';

    % output structure array
    Y = struct;

    % feature names (use based on first sample, assume all the same
% NOTE: first sample is not guaranteed to be a non-control!!    
    aconf = num2str(A(1).results.out.size) + "_" + ...
        num2str(A(1).results.out.peak) + "_" + ...
        regexprep(A(1).results.out.name, ' ', '_');
    sz = A(1).results.out.size;

    % step through flowconfigs to build tables
    for i = 1:length(flowconfigs)

        % indices of samples with current flowconfig
        idx = find(fc == i & c == 0);
        nidx = length(idx);
        B = A(idx);
    
        % dilution
% NOTE: first sample is not guaranteed to be a non-control!!        
        df = B(1).dilution;
    
        % initialize table
        T = table;
        T.experiment = repmat(string(expt), nidx, 1);
        T.donor = string({B.donor}');
        T.day = [B.days]';
        T.fcs = string({B.fcs}');
        T.plate = string({B.plate}');
        T.well = string({B.well}');
        T.dilution = [B.dilution]';
        T.iscontrol = [B.iscontrol]';
        T.concentration = [B.concentration]';
        T.flowconfig = repmat(flowconfigs(i), nidx, 1);
    
        % events data
        events = nan(nidx, length(aconf));
        for k = 1:nidx, events(k, :) = B(k).results.out.events'; end
        eventsheader = "events_" + aconf;
                
        % MFI data per isotype
        iso = strsplit(flowconfigs(i), {',', ':'});
        iso = iso(1:2:end);
        mfi = [];
        mfiheader = [];
        for q = 1:length(iso)
            mat = nan(nidx, length(aconf));
            for k = 1:nidx, mat(k, :) = B(k).results.out.(iso(q))'; end
            mfi = [mfi mat];
            mfiheader = [mfiheader; iso(q) + "_" + aconf];
        end
        
        % get isotypes and bead sizes for each column of mfi data
        H = table;
        H.iso = reshape(repmat(iso, length(aconf), 1), length(iso) * length(aconf), 1);
        H.sz = repmat(sz, length(iso), 1);

        % get standard curve fit params
        F = Q(i).lp4;
        K = table;
        K.isotype = [F.isotype]';
        K.beadsize = [F.beadsize]';
        for w = 1:size(K, 1)
            K.A(w, 1) = F(w).fit.A;
            K.B(w, 1) = F(w).fit.B;
            K.C(w, 1) = F(w).fit.C;
            K.D(w, 1) = F(w).fit.D;
        end
        
        % make blank concentration table
        concheader = "conc_" + mfiheader;
        conc = nan(size(mfi));

        % iterate through isotypes and bead sizes, compute concentrations
        for iiso = 1:length(iso)
            for isz = 1:length(unique(K.beadsize))
                iparam = find(K.isotype == iso(iiso) & K.beadsize == isz);
                pa = K.A(iparam);
                pb = K.B(iparam);
                pc = K.C(iparam);
                pd = K.D(iparam);
                icol = find(H.iso == iso(iiso) & H.sz == isz);
                m = mfi(:, icol);
                temp = (10.^(pc * ((((pa - pd) ./ (m - pd)) - 1) .^ (1 / pb)))) * df / 1000; % ug/ml
                temp(m < pa) = -9999999;
                temp(m > pd) = 9999999;
                conc(:, icol) = temp;            
            end
        end

        % make output table into table
        T = [T,...
            array2table(events, 'VariableNames', eventsheader),...
            array2table(mfi,    'VariableNames', mfiheader),...
            array2table(conc,   'VariableNames', concheader)];

        % add to output
        Y(i).table = T;
        
        % name for text output file
        aname = strsplit(B(1).results.arrayconf, '/');
        aname = regexprep(aname{end}, '\.txt', '');
        fname = ...
            outdir + string(expt) + "_" + ...
            aname + "_dilution_" + num2str(df) + "_" + ...
            regexprep(char(Q(i).flow), {':', ','}, {'_', '_'}) + ".txt";
        
        % write to output file
        writetable(Y(i).table, fname, 'FileType', 'text', 'Delimiter', '\t');

    end
    
return

    





