function struct_push_back(sname,s)
% Add s to struct array sname. 
% yes it's orrible, but I prefer to hide the horror in here
% than to require callers to state the number of fields twice.

if nargin == 0
  % test
  clear mystruct
  tmp.a = 1;
  tmp.b = {'d', 2};
  
  struct_push_back mystruct tmp
  struct_push_back mystruct tmp
  mystruct
  return
end

if ~evalin('caller', ['exist(''' sname ''',''var'')'])
  evalin('caller', [sname '=' s ';']);
else
  evalin('caller', [sname '(end+1)=' s ';']);
end
