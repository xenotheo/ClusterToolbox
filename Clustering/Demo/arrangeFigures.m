function [] = arrangeFigures(n,p,q,pp,mm)
% function [] = arrangeFigures(n,p,q,pp,mm)
% position n figures as an pxq array in full sceen,
% (optional) pp the permutation order
% (optional) mm the screen number (for multiple monitors)
%
% NPP
%

  if ~exist( 'pp', 'var' ) || isempty( pp )
    pp = 1:n;
  end 
  if ~exist( 'mm', 'var' ) || isempty( mm )
    mm = 1;
  end 
  
  % select monitor
  wh = get(0,'MonitorPositions');
  
  mm = min( size( wh, 1 ), mm );
  wh = wh( mm, : );
  
  off = wh(1:2);
  dx  = wh(3) / q;
  dy  = wh(4) / p;

  k = 0;
  for i = 1:p
    for j = 1:q
      k = k + 1;
      if k <= n
        if pp(k) == 0
          continue
        end
        fh = figure(pp(k));
        set(fh,'Position', [(j-1)*dx+off(1),(i-1)*dy+off(2) dx dy/1.125]);
      end
    end
  end
  
end