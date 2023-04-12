function F = cbafcs(fname, isstd, aconf, bchan, fconf, prof, varargin)

% CBAFCS Process serum antibody cytometric assay FCS file
%
%    A = CBAFCS(fcs, isstd, aconf, bchan, fconf, prof) processes an FCS
%    file to quantify antibody reactivities and returns a data structure
%    containing analysis parameters used, raw flow data, transformed flow
%    data and summarized gating and reactivity data.  The required
%    arguments are:
%        1. FCS file with complete path
%        2. "0" for a regular sample, "1" for a standard control sample
%        3. a bead array configuration file
%        4. channel name in FCS file to use for bead separation
%        5. flow configuration string for secondary isotypes
%        6. a peak finding parameter settings file
%
%    CBAFCS(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%    parameter/value pair values:
%        'pdf'          create a PDF report containing parameters, figures
%                       and summarized data (the display option must be set
%                       to "yes")
%        'txt'          create a tab delimited text file, with analysis
%                       parameters as header rows followed by reactivity
%                       table
%        'outdir'       specify a directory for the PDF and/or the TXT
%                       file to be written to - default is the directory
%                       from which this function is run
%        'display'      whether to render plots - default is 'yes'
%
%    Usage:
%
%        A = CBAFCS(...
%               'FCS/myflowfile.fcs',...
%               0,...
%               'config/mybeadarrayconf.txt',...
%               'FL4-A',...
%               'IgM:FL1-A,IgA:FL2-A,IgG:FL5-A',...
%               'config/mypeakparams.txt');
%
%        A = CBAFCS(...
%               'FCS/myflowfile.fcs',...
%               0,...
%               'config/mybeadarrayconf.txt',...
%               'FL4-A',...
%               'IgM:FL1-A,IgA:FL2-A,IgG:FL5-A',...
%               'config/mypeakparams.txt',...
%               'pdf', 'yes',...
%               'outdir', '~/Desktop');
%
%    Uses the function "fca_readfcs.m" to read FCS files
%    Copyright (c) 2020, Laszlo Balkay, All rights reserved.
%
%    Uses the function "dscatter.m" to render 2-D scatter/density plots
%    Copyright (c) 2016, MathWorks Inc. All rights reserved.
%
%    Requires Matlab Signal Processing Toolbox
%    Requires Matlab Curve Fitting Toolbox
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


    % version
    ver = '1.8';

    % validate input args
    ty = {'string', 'char'};
    ne = {'nonempty'};
    pa = inputParser;
    addRequired(pa,  'fname',          @(x) validatefile(x));
    addRequired(pa,  'isstd',          @(x) ismember(x, [0 1]));
    addRequired(pa,  'aconf',          @(x) validatefile(x));
    addRequired(pa,  'bchan',          @(x) validateattributes(x, ty, ne));
    addRequired(pa,  'fconf',          @(x) validateattributes(x, ty, ne));
    addRequired(pa,  'prof',           @(x) validatefile(x));
    addParameter(pa, 'pdf',     'no',  @(x) validateattributes(x, ty, ne));
    addParameter(pa, 'txt',     'no',  @(x) validateattributes(x, ty, ne));
    addParameter(pa, 'display', 'yes', @(x) validateattributes(x, ty, ne));
    addParameter(pa, 'outdir',  './',  @(x) validateattributes(x, ty, ne));
    parse(pa, fname, isstd, aconf, bchan, fconf, prof, varargin{:});    
    if lower(string(pa.Results.pdf)) == "yes",     rp =  1; else, rp =  0; end    
    if lower(string(pa.Results.txt)) == "yes",     rpx = 1; else, rpx = 0; end 
    if lower(string(pa.Results.display)) == "yes", ds =  1; else, ds =  0; end
    outdir = pa.Results.outdir;
    
    % load FCS file and initialize data structure
    F1 = cbaload(fname, isstd, aconf, bchan, fconf, prof, ds, ver);
    
    % separate beads by size
    F2 = cbasize(F1);
    
    % retain singlets
    F3 = cbasinglet(F2);
    
    % split bead sizes if necessary
    F4 = cbasplit(F3);
    
    % find bead peaks for each bead size series
    F5 = cbapeak(F4);
    
    % find MFI for each feature and isotype
    F = cbamfi(F5);
    
    % if text file output
    if rpx == 1, cbatxt(F, outdir); end
    
    % if PDF file output
    if rp == 1, cbapdf(F, outdir); end
    
return
    
    
