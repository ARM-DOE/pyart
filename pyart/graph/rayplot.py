"""
Colored ray plotter
needs more documentation
"""

from pylab import figure, plot, title, xlabel, ylabel, subplot, xlim, ylim
import numpy as np


class rayplot:

    def __init__(self, radar, ncp_min=0.5, rhv_min=0.8, **kwargs):
        self.radar = radar
        self.ncp_min = ncp_min
        self.rhv_min = rhv_min
        self.is_coh = (radar.fields['norm_coherent_power']['data'] >
                       self.ncp_min)
        self.is_cor = radar.fields['copol_coeff']['data'] > self.ncp_min
        self.is_good = np.logical_and(self.is_coh, self.is_cor)
        self.good_idx = np.where(self.is_good)
        self.bad_idx = np.where(np.logical_not(self.is_good))
        self.distance = radar.range['data']
        self.figsize = kwargs.get('figsize', [15, 10])

    def __call__(self, ray, vars, limits=[0, -1]):
        nfig = len(vars)
        self.fig = figure(figsize=self.figsize)
        ax_list = []
        plot_list = {}
        i = 1
        ray_number = self.radar.sweep_start_ray_index['data'][ray[0]] + ray[1]
        for var in vars:
            cur_ax = subplot(nfig, 1, i)
            cur_good = plot(
                self.distance[limits[0]:limits[1]],
                np.ma.masked_where(np.logical_not(self.is_good),
                                   self.radar.fields[var]['data'])
                [ray_number, limits[0]:limits[1]])

            cur_bad = plot(
                self.distance[limits[0]:limits[1]],
                np.ma.masked_where(self.is_good,
                                   self.radar.fields[var]['data'])
                [ray_number, limits[0]:limits[1]], 'r-')

            ylabel(self.radar.fields[var]['units'])
            title(self.radar.fields[var]['standard_name'].replace('_', ' '))
            cur_ax.yaxis.grid(color='gray', linestyle='dashed')
            if i != nfig:
                cur_ax.get_xaxis().set_ticks([])
            else:
                xlabel(self.radar.range['standard_name'] + ' (' +
                       self.radar.range['units'] + ')')
            if var == 'copol_coeff':
                ylim([.5, 1.0])
            plot_list.update({var + '_good': cur_good})
            plot_list.update({var + '_bad': cur_bad})
            ax_list.append(cur_ax)
            i = i + 1
