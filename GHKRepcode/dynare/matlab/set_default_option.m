function options=set_default_option(options,field,default)
  if ~isfield(options,field)
    if sscanf(version('-release'),'%d') < 13
      options = setfield(options,field,default);
    else
      eval('options.(field) = default;');
    end
  end
  
  % 06/07/03 MJ added ; to eval expression