import matplotlib.pyplot as plt
import numpy as np


def plot_convstrat_scheme(
    always_core_thres,
    use_cosine,
    max_diff=None,
    zero_diff_cos_val=None,
    use_addition=False,
    scalar_diff=None,
):
    """
    Plots the scheme used in the convective stratiform classification

    Parameters
    ----------
    always_core_thres : float
        All values above this threshold considered to be convective
    use_cosine : bool
        Boolean used to determine if cosine scheme should be used for identifying convective cores (True) or a scalar
        scheme (False)
    max_diff : float, optional
        Maximum difference between background average and reflectivity in order to be classified as convective.
        "a" value in Eqn. B1 in Yuter and Houze (1997)
    zero_diff_cos_val : float, optional
        Value where difference between background average and reflectivity is zero in the cosine function
        "b" value in Eqn. B1 in Yuter and Houze (1997)
    use_addition : bool, optional
        Determines if a multiplier (False) or addition (True) in the scalar difference scheme should be used
    scalar_diff : float, optional
        If using a scalar difference scheme, this value is the multiplier or addition to the background average

    """

    # create array of background values
    bkg_vals = np.linspace(0, 60, 100)

    # create difference array
    if use_cosine:
        # cosine scheme
        diff = max_diff * np.cos(np.pi * bkg_vals / (2 * zero_diff_cos_val))
    else:
        if use_addition:
            # scalar addition scheme
            diff = (bkg_vals + scalar_diff) - bkg_vals
        else:
            # scalar multiplier scheme
            diff = (bkg_vals * scalar_diff) - bkg_vals

    # if values are less than zero, set to zero
    diff[diff < 0] = 0
    # background values greater than always core thres, set to zero
    diff[bkg_vals > always_core_thres] = 0

    # Now plot
    ax = plt.gca()
    # plot difference line
    ax.plot(bkg_vals, diff, lw=2, color="black")
    # plot always core thres
    ax.axvline(x=always_core_thres, lw=1, ls="--", color="red")
    ax.text(always_core_thres + 2, 1, "Always Core Thres.", color="red")
    if use_cosine:
        # plot zero difference cosine value
        ax.axvline(x=zero_diff_cos_val, lw=1, ls="--", color="green")
        ax.text(zero_diff_cos_val + 2, 0.75, "Zero Diff. Cosine Val.", color="green")
        # plot max difference
        ax.axhline(y=max_diff, lw=1, ls="--", color="blue")
        ax.text(10, max_diff + 0.05, "Max. Diff.", color="blue")
    elif use_addition:
        # plot scalar
        ax.axhline(y=scalar_diff, lw=1, ls="--", color="orange")
        ax.text(10, scalar_diff + 0.05, "Scalar Diff.", color="orange")
    # add grid
    ax.grid()
    # set axis limits
    ax.set_ylim([0, np.max(diff) + 0.2])
    ax.set_xlim([np.min(bkg_vals), np.max(bkg_vals)])
    # set axis labels and title
    ax.set_ylabel("Difference (dBZ - dBZ$_{background}$)")
    ax.set_xlabel("Background Value (dBZ$_{background}$)")
    ax.set_title("Convective Stratiform Equation")
