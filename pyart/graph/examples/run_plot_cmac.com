rm cmac_rhi.dir/figure_multi_panel_hsrhi_20110524.015604.png 
./plot_cmac.py sgpxsaprrhicmacI5.c0.20110524.015604_NC4.nc cmac_rhi.dir figure_
display cmac_rhi.dir/figure_multi_panel_hsrhi_20110524.015604.png &
display cmac_rhi.dir/reference.png &

