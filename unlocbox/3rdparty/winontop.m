function WasOnTop = winontop( FigureHandle, IsOnTop )
%WINONTOP allows to trigger figure's "Always On Top" state
% INPUT ARGUMENTS:
%  FigureHandle - Matlab's figure handle, scalar
%  IsOnTop      - logical scalar or empty array
%
% USAGE:
%  WinOnTop( hfigure, bool );
%  WinOnTop( hfigure );            - equal to WinOnTop( hfigure,true);
%  WinOnTop();                     - equal to WinOnTop( gcf, true);
%  WasOnTop = WinOnTop(...);       - returns boolean value "if figure WAS on top"
%  IsOnTop  = WinOnTop(hfigure,[]) - gets "if figure is on top" property
% 
% LIMITATIONS:
%  java enabled
%  figure must be visible
%  figure's "WindowStyle" should be "normal"
%
%
% Written by Igor
% i3v@mail.ru
%
% 16 June 2013 - Initial version
% 27 June 2013 - removed custom "ishandle_scalar" function call
%
%
% Copyright (c) 2013, Igor
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%      Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%      Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/3rdparty/winontop.php

% Copyright (C) 2012-2013 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.5.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Parse Inputs

if ~exist('FigureHandle','var');FigureHandle = gcf; end
assert(...
          isscalar( FigureHandle ) && ishandle( FigureHandle ) &&  strcmp(get(FigureHandle,'Type'),'figure'),...
          'WinOnTop:Bad_FigureHandle_input',...
          '%s','Provided FigureHandle input is not a figure handle'...
       );

assert(...
            strcmp('on',get(FigureHandle,'Visible')),...
            'WinOnTop:FigInisible',...
            '%s','Figure Must be Visible'...
       );

assert(...
            strcmp('normal',get(FigureHandle,'WindowStyle')),...
            'WinOnTop:FigWrongWindowStyle',...
            '%s','WindowStyle Must be Normal'...
       );
   
if ~exist('IsOnTop','var');IsOnTop=true;end
assert(...
          islogical( IsOnTop ) &&  isscalar( IsOnTop) || isempty( IsOnTop ), ...
          'WinOnTop:Bad_IsOnTop_input',...
          '%s','Provided IsOnTop input is neither boolean, nor empty'...
      );
%% Pre-checks

error(javachk('swing',mfilename)) % Swing components must be available.
  
  
%% Action

% Flush the Event Queue of Graphic Objects and Update the Figure Window.
drawnow expose

jFrame = get(handle(FigureHandle),'JavaFrame');

drawnow

WasOnTop = jFrame.fHG1Client.getWindow.isAlwaysOnTop;

if ~isempty(IsOnTop)
    jFrame.fHG1Client.getWindow.setAlwaysOnTop(IsOnTop);
end

end


