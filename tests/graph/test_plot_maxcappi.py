import os

import matplotlib.pyplot as plt

import pyart


def test_plot_maxcappi_simple(outfile=None):
    """
    Test the basic functionality of plot_maxcappi.
    """
    # Create a test grid using Py-ART's testing utility
    grid = pyart.testing.make_target_grid()
    grid.z["data"] = grid.z["data"] * 10 + 100
    grid.metadata["instrument_name"] = "GRID"

    # Use plot_maxcappi with the generated grid
    pyart.graph.max_cappi.plot_maxcappi(
        grid=grid,
        field="reflectivity",
        savedir=None,  # Do not save the plot
        show_figure=False,  # Do not show the plot
    )

    if outfile:
        plt.savefig(outfile)
    plt.close()


def test_plot_maxcappi_with_save(outfile=None):
    """
    Test plot_maxcappi and save the output to a file.
    """
    # Create a test grid using Py-ART's testing utility
    grid = pyart.testing.make_target_grid()
    grid.z["data"] = grid.z["data"] * 10 + 100
    grid.metadata["instrument_name"] = "GRID"

    # Define the output file path
    outfile = outfile or "test_plot_maxcappi_output.png"

    # Use plot_maxcappi with the generated grid
    pyart.graph.max_cappi.plot_maxcappi(
        grid=grid,
        field="reflectivity",
        savedir=None,  # Handle saving manually below
        show_figure=False,  # Do not show the plot
    )

    # Save the figure to a file
    plt.savefig(outfile)
    plt.close()

    # Check if the file was created
    assert os.path.exists(outfile), "The plot was not saved as expected."


def test_plot_maxcappi_with_all_options(outfile=None):
    """
    Test plot_maxcappi with all options enabled.
    """
    import cartopy.crs as ccrs

    # Create a test grid using Py-ART's testing utility
    grid = pyart.testing.make_target_grid()
    grid.z["data"] = grid.z["data"] * 10 + 100
    grid.metadata["instrument_name"] = "GRID"

    # Use a custom projection for testing
    projection = ccrs.Mercator()

    # Use plot_maxcappi with additional options
    pyart.graph.max_cappi.plot_maxcappi(
        grid=grid,
        field="reflectivity",
        title="Test Max-CAPPI",
        lat_lines=None,
        lon_lines=None,
        add_map=True,
        projection=projection,
        colorbar=True,
        range_rings=True,
        dpi=150,
        savedir=None,
        show_figure=False,
    )

    if outfile:
        plt.savefig(outfile)
    plt.close()


if __name__ == "__main__":
    test_plot_maxcappi_simple("figure_plot_maxcappi_simple.png")
    test_plot_maxcappi_with_save("figure_plot_maxcappi_output.png")
    test_plot_maxcappi_with_all_options("figure_plot_maxcappi_all_options.png")
