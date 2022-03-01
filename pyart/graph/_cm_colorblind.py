"""
Data for colorblind friendly radar colormaps

"""
import numpy as np
import os


def yuv_rainbow_24(nc):
    path1 = np.linspace(0.8*np.pi, 1.8*np.pi, nc)
    path2 = np.linspace(-0.33*np.pi, 0.33*np.pi, nc)

    y = np.concatenate([np.linspace(0.3, 0.85, nc*2//5),
                        np.linspace(0.9, 0.0, nc - nc*2//5)])
    u = 0.40*np.sin(path1)
    v = 0.55*np.sin(path2) + 0.1

    rgb_from_yuv = np.array([[1, 0, 1.13983],
                             [1, -0.39465, -0.58060],
                             [1, 2.03211, 0]])
    cmap_dict = {'blue': [], 'green': [], 'red': []}
    for i in range(len(y)):
        yuv = np.array([y[i], u[i], v[i]])
        rgb = rgb_from_yuv.dot(yuv)
        red_tuple = (i/(len(y)-1.0), rgb[0], rgb[0])
        green_tuple = (i/(len(y)-1.0), rgb[1], rgb[1])
        blue_tuple = (i/(len(y)-1.0), rgb[2], rgb[2])
        cmap_dict['blue'].append(blue_tuple)
        cmap_dict['red'].append(red_tuple)
        cmap_dict['green'].append(green_tuple)

    return cmap_dict


# balance colormap from cmocean project
# courtesy Kristen Thyng

# Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M., & DiMarco, S. F. (2016).
# True colors of oceanography. Oceanography, 29(3), 10.

# HomeyerRainbow developed by Cameron Homeyer with assistance from Bobby Jackson


data_dir = os.path.split(__file__)[0]
bal_rgb_vals = np.genfromtxt(os.path.join(data_dir, 'balance-rgb.txt'))

datad = {'HomeyerRainbow': yuv_rainbow_24(15),
         'balance': bal_rgb_vals}
