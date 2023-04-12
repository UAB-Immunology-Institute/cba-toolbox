function validatefile(arg)

% VALIDATEFILE Throw error if not a valid file
%
%    VALIDATEFILE(arg) throws an error if the input argument is not a valid
%    file.  This function is expecting scalar input (a filename which is
%    either a character array or a scalar string
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


    if ~(ischar(arg) || (isstring(arg) && isscalar(arg)))
        error('expecting one item, either a character array or string');
    end

    if ~isfile(arg)
        error(['file not found: ' char(arg)]);
    end

return
    