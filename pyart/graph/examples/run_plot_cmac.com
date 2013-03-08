rm cmac_rhi.dir/figure_multi_panel_hsrhi_20110524.015604.png 
rm test_cmac_rhi.png
./plot_cmac.py sgpxsaprrhicmacI5.c0.20110524.015604_NC4.nc cmac_rhi.dir figure_

./plot_cmac_radar.py sgpxsaprrhicmacI5.c0.20110524.015604_NC4.nc test_cmac_rhi.png

display cmac_rhi.dir/figure_multi_panel_hsrhi_20110524.015604.png &
display cmac_rhi.dir/reference.png &
display test_cmac_rhi.png &


