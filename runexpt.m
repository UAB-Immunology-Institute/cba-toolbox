function X = runexpt(pt, aconf, prof, manifest)

% RUNEXPT Wrapper to run CBAFCS on a directory of FCS files
%
%    X = RUNEXPT(fcsdir, aconf, prof, manifest) processes a batch of FCS 
%    files in a specified directory (first argument) based on a specified 
%    bead array design text file (second argument), a peak-finding  
%    parameter profile text file (third file) and a sample manifest listing
%    all the FCS files to process (fourth argument).  The FCS files in the
%    manifest file must be present in the specified directory.  The output
%    is a structure array (one element per FCS file processed).  Each
%    structure element contains the fields present in the manifest for the
%    given FCS file as well as a sub-structure (the output of CBAFCS)
%    containing the analysis results.
%
%    Usage:
%
%        A = RUNEXPT(...
%               '/path/to/my/FCSfiles/',...
%               '/Users/afr/config/myarrayconfig.txt',...
%               '/Users/afr/config/mypeakparams.txt',...
%               '/Users/afr/beadarray/manifests/expt100.txt');
%
%    Uses the function "fca_readfcs.m" to read FCS files
%    Copyright (c) 2020, Laszlo Balkay, All rights reserved.
%
%    Uses the function "dscatter.m" to render 2-D scatter/density plots
%    Copyright (c) 2016, MathWorks Inc. All rights reserved.
%
%    Requires Matlab Signal Processing Toolbox and Curve Fitting Toolbox
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


    % read in header lines (beginning with "#")
    h = {};
    fid = fopen(manifest);
    while 1
        a = fgetl(fid);
        a = regexprep(a, '"', '');
        if isempty(regexp(a, '^#', 'once')), break; end
        a = regexprep(a, '\s*$', '');
        h = [h; cellstr(a)];
    end
    fclose(fid);
    nh = length(h);
    
    % parse attribute-value pairs and make a property-value table
    av = cell(nh, 2);
    for i = 1:nh, av(i, :) = strsplit(regexprep(h{i}, '^#', ''), '='); end
    K = cell2table(av(:, 2), 'VariableNames', {'value'});
    K.Properties.RowNames = string(av(:, 1));
    K.value = string(K.value);
    
    % get bead channel
    bead = K.value("beadchannel");
    
    % get timepoint column
    tp = K.value("timepoint");
    
    % now, knowing how many header lines, read in table
    T = readtable(manifest, 'NumHeaderLines', nh);
    
    % number of samples to process
    n = size(T, 1);
    
    % output structure for quantification results
    X = struct;
    
    % samples to keep (those without errors);
    keep = ones(n, 1);
    
    % add field to table for log output
    T.status = repmat({'ok'}, n, 1);
    T.number = (1:n)';
    
    % iterate through samples in manifest
    for i = 1:n
        X(i).donor         = T.donor{i};
        X(i).days          = T.(tp)(i);
        X(i).fcs           = T.fcs{i};
        X(i).plate         = T.plate{i};
        X(i).well          = T.well{i};        
        X(i).dilution      = T.dilution(i);
        X(i).iscontrol     = T.control(i);
        X(i).concentration = T.conc(i);
        try
            X(i).results = cbafcs(...
                [pt T.fcs{i}],...
                T.control(i),...
                aconf,...
                bead,...
                K.value(T.flow{i}),...
                prof,...
                'display', 'off');
            disp(['file ' num2str(i) ' of ' num2str(n) ': ok']);
        catch e
            disp(['file ' num2str(i) ' of ' num2str(n) ': ' e.message]);
            keep(i) = 0;
            T.status(i) = cellstr(e.message);
        end
    end
    
    % summary message
    disp([num2str(length(find(keep == 1))) ' / ' num2str(n) ' samples processed']);
    disp([num2str(length(find(keep == 0))) ' / ' num2str(n) ' samples omitted']);
                
    % error logs
    EL = T(keep == 0, {'number', 'fcs', 'status'});
    
    % export logfile
    if ~isempty(find(keep == 0, 1))
        dd = datetime;
        dd.Format = 'uuuuMMdd_HHmmss';
        logname = "cba_errlog_" + string(dd) + ".txt";
        writetable(EL, logname, 'Delimiter', '\t', 'WriteVariableNames', false);
        disp("see error log: " + logname);
    end
    
    % only return samples that were successful
    X = X(keep == 1);
    
return


