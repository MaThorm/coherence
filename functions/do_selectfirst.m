function [ft_str] = do_selectfirst(ft_str)
% Takes first channel of ft_Struct
for ii = 1:length(ft_str)
    cfg = [];
    cfg.channel = ft_str(ii).label{1,1}
    ft_str(ii) = ft_selectdata(cfg,ft_str(ii))
end 
end

