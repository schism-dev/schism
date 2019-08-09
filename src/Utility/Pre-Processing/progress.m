%Author: Ivica Janekovic
function progress(size,current,text_mode)

% progress(size,current,text_mode)

global progress_last_pix;
global progress_text;
global progress_handle;

length = 0;
cur = 0;

if nargin >= 1
  length = size;
end
if nargin >= 2
  cur = current;
end
if nargin == 3
  progress_text = text_mode;
end
if ( exist('progress_text') )
  if ( isempty(progress_text) )
    progress_text = 0;
  end
end

% If you want text
if ( progress_text == 1 )
  scale = 100/60;
  if ( length == 0 & cur == 0 )
    progress_last_pix = 0;
    fprintf('0');
    i = 1;
    while i<=60
      if ( mod(ceil(i*scale),25) <= .5 )
        fprintf('%03d',ceil(i*scale));
        i = i + 2;
      else
        fprintf('.')
      end
      i = i + 1;
    end
    fprintf('\n');
    tic;
  % Draw the current position
  elseif ( length > 0 & cur <= length ) 
    % Map to our scale
    pix = ceil((cur/length * 100)/scale);
    if ( pix == 60 & progress_last_pix ~= 60 )
      fprintf('==== %6.2f sec\n',toc);
    elseif ( pix > progress_last_pix )
      for i=progress_last_pix:pix-1,
        fprintf('=');
      end
    end
    progress_last_pix = pix;
  end
else
  % If you have graphics
  if ( length == 0 & cur == 0 )
    if ( ishandle(progress_handle) )
      close(progress_handle);
    end
    progress_handle=waitbar(0, 'Progress...');
    h_menu=uimenu(progress_handle,'Label','Control');
    h1=uimenu(h_menu,'Label','Pause','Callback','display(''paused'');pause;display(''running'');','Accelerator','s');
    h2=uimenu(h_menu,'Label','Suspend','Callback','keyboard','Accelerator','z');
    tic;
  else
    waitbar(cur/length);
    if ( cur/length >= 1 )
      if ( ishandle(progress_handle) )
        close(progress_handle);
      end    
      toc
    end
  end
end
