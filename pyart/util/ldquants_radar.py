#!/usr/bin/env python
"""
 NAME:
   ldquants_radar.py

 PURPOSE:
   To read in the NWS Houston (KHGX) NEXRAD reflectivity archived data
   along with the LDQUANTs laser disdrometer data for comparison.

 SYNTAX:
   python ldquants_radar.py houldquantsM1.c1.YYYYMMDD.HHMMSS.nc
          houldquantsS1.c1.YYYYMMDD.HHMMSS.nc
          houvdisquantsM1.c1.20220109.000000.nc

 INPUT:
   houldquantsM1.c1.YYYYMMDD.HHMMSS.nc
       LDQUANTS   (M1 Site) netCDF file.
   houldquantsS1.c1.YYYYMMDD.HHMMSS.nc
       LDQUANTS   (S1 Site) netCDF file.
   houvdisquantsM1.c1.20220109.000000.nc
       VDISQUANTS (M1 Site) netCDF file.

 KEYWORDS:

 EXECUTION EXAMPLE:
   Linux example: python ldquants_radar.py
                  houldquantsM1.c1.20220109.000000.nc
                  houldquantsS1.c1.20220109.000000.nc

 MODIFICATION HISTORY:
   2022/02/18 - Joe O'Brien <obrienj@anl.gov> :
                Created using examples from the Py-Art quick start
                guide and LDQUANTS documentation
   2022/03/08 - Joe O'Brien <obrienj@anl.gov> :
                Updated to transition to pulling KHGX data from
                AWS server instead of command line input.
                Added 'pull_allscans' function, moved plotting
                into 'plot_ppirhi' function.
                Updated metadata for functions
   2022/03/14 - Joe O'Brien <obrienj@anl.gov> :
                Added 'plot_his' function for plotting
                comparison of equivalent reflectivity factor for
                the whole day.
   2022/03/16 - Joe O'Brien <obrienj@anl.gov> :
                Added handling/plotting for LDQUANTS S1 site and
                VDISQUANTS M1 site
   2022/03/25 - Joe O'Brien <obrienj@anl.gov> :
                Updating to conform to PEP8 and pylint testing

 NOTES:
   1) PyART quick start guide:
       https://arm-doe.github.io/pyart/
   2) NOAA NEXRAD ARCHIVE:
       ncdc.noaa.gov/nexradinv/
   3) PyART Github
       github.com/ARM-DOE/pyart
   4) LDQUANTS documentation
       - https://www.arm.gov/capabilities/vaps/ldquants
   5) Github used to store this analysis
       - https://github.com/ARM-Development/ldquants-radar
   6) Notebook used for nexradaws data pull/display
       - https://github.com/scollis/notebooks/blob/master/urban
            /Boston%20Snow%20Retrieval.ipynb
   7) Python nexradaws documentation
       nexradaws.readthedocs.io/_/downloads/en/latest/pdf/
"""

import sys
import time
import datetime
#import tempfile

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
##import cartopy.crs as ccrs
##import cartopy.feature as cfeature
##from cartopy.feature import NaturalEarthFeature

##from netCDF4 import Dataset
import netCDF4
import pyart
import nexradaws

#-------------------------
# I) Define Functions
#-------------------------

# Command line syntax statement.
def help_message():
    """
    To display to the command line the expected syntax for this
    program including any/all keywords or options

    Returns
    -------
    Print statements to the command line and closing of the program.
    """
    print("\n")
    print("Syntax: ldquantsRADAR.py <-h> LDQUANTS_M1_file"
          + "LDQUANTS_S1_file\n")
    print("  INPUT:  ")
    print("    LDQUANTS_M1_file:"
          + "- M1 site disdrometer LDQUANTS VAP file")
    print("    LDQUANTS_S1_file:"
          + "- S1 site disdrometer LDQUANTS VAP file\n")
    print("  OPTIONS:  ")
    print("    -h:"
          + "- Help statement. Print Syntax\n")
    print("  KEYWORDS: \n")
    print("EXAMPLE: python ldquantsRADAR.py"
          + " houldquantsM1.c1.20220109.000000.nc"
          + " houldquantsS1.c1.20220109.000000.nc\n")

# Distance between radar and target
# Distance between two points on a sphere.
def sphere_distance(rad_lat,tar_lat,rad_lon,tar_lon):
    """
    Calculation of the great circle distance between radar and target

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculated for a smooth sphere
    Radar and Target are at the same altitude (need to check)

    Parameter
    ---------
    rad_lat   : float, [degrees]
               latitude of the radar in degrees
    tar_lat   : float, [degrees]
               latidude of the target in degrees
    rad_lon   : float, [degrees]
               longitude of the radar in degrees
    tar_lon   : float, [degrees]
               longitude of the target in degress

    Returns
    -------
    distance : float, [meters]
               Great-Circle Distance between radar and target in meters
    """
    # convert latitude/longitudes to radians
    rad_lat = rad_lat * (np.pi/180.)
    tar_lat = tar_lat * (np.pi/180.)
    rad_lon = rad_lon * (np.pi/180.)
    tar_lon = tar_lon * (np.pi/180.)
    # difference in latitude  - convert from degrees to radians
    d_lat = (tar_lat - rad_lat)
    # difference in longitude - convert from degress to radians
    d_lon = (tar_lon - rad_lon)
    # Haversince formula
    numerator = ((np.sin(d_lat/2.0)**2.0)
                + np.cos(rad_lat) * np.cos(tar_lat) * (np.sin(d_lon/2.0)**2.0))
    distance = 2 * 6371000 * np.arcsin(np.sqrt(numerator))

    # return the output
    return distance

# Great Circle Bearing Calculation - Forward Azimuth Angle
def for_azimuth(rad_lat, tar_lat, rad_lon, tar_lon):
    """
    Calculation of inital bearing along a great-circle arc
    Known as Forward Azimuth Angle.

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculated for a smooth sphere
    Radar and Target are at the same altitude (need to check)

    Parameters
    ----------
    rad_lat  : float, [degrees]
              latitude of the radar in degrees
    tar_lat  : float, [degrees]
              latidude of the target in degrees
    rad_lon  : float, [degrees]
              longitude of the radar in degrees
    tar_lon  : float, [degrees]
              longitude of the target in degress

    Returns
    -------
    azimuth : float, [degrees]
        azimuth angle from the radar where
        target is located within the scan.
        output is in degrees.
    """
    # convert latitude/longitudes to radians
    rad_lat = rad_lat * (np.pi/180.)
    tar_lat = tar_lat * (np.pi/180.)
    rad_lon = rad_lon * (np.pi/180.)
    tar_lon = tar_lon * (np.pi/180.)
    # Differnce in longitudes
    d_lon = tar_lon - rad_lon
    # Determine x,y coordinates for arc tangent function
    corr_x = np.sin(d_lon) * np.cos(tar_lat)
    corr_y = ((np.cos(rad_lat) *np.sin(tar_lat))
             - (np.sin(rad_lat) * np.cos(tar_lat) * np.cos(d_lon)))
    # Determine forward azimuth angle
    azimuth = np.arctan2(corr_y,corr_x) * (180./np.pi)

    # Return the output
    return azimuth

# download all valid scans for a given date
def pull_allscans(connection, site, ndatetime, outdir=None):
    """
    To search for all available scans for a given NEXRAD
    radar site for a given date and download the data.

    Parameters
    ----------
    connection : nexradaws
                 nexradaws connection interface object
                 for the given date and site.
    site       : str
                 NWS NEXRAD Radar site identifier to
                 download data for.
    ndatetime  : Datetime
                 datetime object containing date to download
                 data for. Time zone needs to be set in UTC.
    outdir     : str, optional
                 directory path to store downloaded radar scans
                 to. Maintains AWS directory structure.

    Returns
    -------
    downloads  : nexradaws
                 nexradaws object containing downloaded files and paths
                 to these files
    """
    # Create a temporary directory to hold the data.
    ##tlocation   = tempfile.mkdtemp()
    # Check to see if data is available this date.
    if str(ndatetime.year) not in connection.get_avail_years():
        print('ERROR: Date Not Found for Radar Site\n')
        sys.exit(1)
    else:
        if (f'{ndatetime.month:02}' not in
                connection.get_avail_months(ndatetime.year)):
            print('ERROR: Date Not Found For Radar Site\n')
            sys.exit(1)
        else:
            if (f'{ndatetime.day:02}' in
                   connection.get_avail_days(ndatetime.year,ndatetime.month)):
                print('ERROR: Date Not Found For Radar Site\n')
                sys.exit(1)
    # Check to see if data is available for the site.
    site_check = connection.get_avail_radars(ndatetime.year, ndatetime.month,
                                            ndatetime.day)
    if site not in site_check:
        print('ERROR: Radar Site Not Found\n')
        help_message()
        sys.exit(1)
    # Grab all the scans for the date.
    # Use AWS server to grab the available scans based on the datetime
    # object and site.
    these_scans = connection.get_avail_scans(ndate.year,ndate.month,
                                             ndate.day,site)
    # Need to clean
    these_good_scans = []
    these_good_times = []
    for i in enumerate(these_scans):
        if i[1] is not None:
            these_good_times.append(i[1])
            these_good_scans.append(i[1])
    # download the data. search if output directory is defined.
    if outdir is not None:
        ndownloads = connection.download(these_good_scans,outdir,
                                        keep_aws_folders=True)
    else:
        ndownloads = connection.download(these_good_scans,site,
                                        keep_aws_folders=True)

    # return the downloaded files
    return ndownloads

def plot_ppirhi(display,nradar,ldquants,tsync,nazimuth,ndistance,ncolumn):
    """
    To display the PPI and peusdo-RHI scan for a given PyART
    RadarDisplay object for comparison with a ground site
    or target.

    Parameters
    ----------
    display   : RadarDisplay
                PyART RadarDisplay object from which reflectivity factor
                is displayed
    nradar    : Radar
                PyART Radar object form which metadata/time are used for
                filename creation
    ldquants  : netCDF4
                netCDF object of the VAP LDQUANTs data from which
                lat/lon and rain rate are displayed.
    tsync     : masked array
                numpy masked array of indices where time between
                ldquants/radar match
    nazimuth  : float, [degrees]
                azimuth value to display with peusdo-RHI scan
    ndistance : float, [km]
                distance from the radar to observed column
    ncolumn   : array
                radar derived rain rates for column above observation site

    Returns
    -------
        plot : *png
               Saves *png file for the plot.
    """
    # create the figure
    fig, axarr = plt.subplots(1,2,figsize=(12.80,6.45))
    # Plot a PPI of the lowest sweep.
    # For whatever reason, pyart.graph.RadarMapDisplay.plot_ppi_map
    # was not playing well with subplots
    display.plot_ppi('reflectivity',sweep=0,axarr=axarr[0])
    # Determine limits for the PPI plot, force to be square
    axarr[0].set_ylim(-100,100)
    axarr[0].set_xlim(-100,100)
    # Define a point on the graph for the LDQUANTS site.
    display.plot_label('M1',(ldquants.variables['lat'][:],
                       ldquants.variables['lon'][:]), ax=axarr[0], symbol='ro',
                       text_color='k')
    # Plot the psendo-RHI using the forward azimuth angle calculated
    # from for_azimuth function.
    display.plot_azimuth_to_rhi('reflectivity',nazimuth, ax=axarr[1])
    # Since this is disdrometer-radar comparison, shortening up plot limits
    axarr[1].set_xlim(0,50)
    axarr[1].set_ylim(0,20)
    # Highlight the locaiton of the disdrometer data using the distance
    # calculated in the sphere_distance function
    axarr[1].axvline((ndistance/1000.),0,1, color='k')
    # Plot text of the estimated rainfall rate from the radar gates over
    # the disdrometer site.
    axarr[1].text((ndistance/1000.) + 1,19,'Radar RR: '
             + str(np.around(np.sum(ncolumn),2))
             +' [mm/hr]',fontsize=9)
    axarr[1].text((ndistance/1000.) + 1,18,'     M1 RR: '
            + str(np.around(np.sum(ldquants.variables['rain_rate'][tsync]),2))
            +' [mm/hr]',fontsize=9)
    # Save the figure.
    nout = ('ldquantsRadar_' + nradar.metadata['instrument_name'] + '_'
            + nradar.time['units'].split(' ')[-1].split('T')[0] + '_'
            + begin_time[0:2] + begin_time[3:5] + begin_time[6:8]+'.png')
    plt.savefig(nout, dpi=100)
    # Close the figure.
    plt.close(fig)

def plot_his(data_m1,data_s1,brange,nbin):
    """
    To calculate the bi-dimensional histogram of two data samples
    and display

    Parameters
    ----------
    data_m1 : dict
             Dictionary containing int/float values of radar field
             parameters from the Houston NEXRAD radar and M1-site
             LDQUANTS data that are desired to be displayed
    data_s1 : dict
             Dictionary containing int/float values of radar
             field parameters from the Houston NEXRAD radar and
             M1-site LDQUANTS data that are desired to be displayed
    brange : list
             List in the format [xmin, xmax, ymin, ymax] for values
             to be binned within a 2D histogram
    nbin   : int,float
             value for the number of bins to apply to the 2D histogram

    Returns
    -------
    plot   : *png
             Saves *png file for the plot.
    """
    # create the subplots
    fig, axarr = plt.subplots(2,2,figsize=(12.80,6.45))
    # Plot the data.
    axarr[0,0].hist2d(data_m1['LD_Z'],data_m1['rhi_data'],bins=nbin,
                   range=[[brange[0],brange[1]],[brange[2],brange[3]]],
                   norm=mpl.colors.LogNorm(),cmap=pyart.graph.cm_colorblind.HomeyerRainbow)
    # Add a 1:1 ratio line
    axarr[0,0].plot(np.arange(brange[0]-10,brange[1]+10,1),
                 np.arange(brange[0]-10,brange[1]+10,1),'k')
    # Define axe titles for the figure.
    axarr[0,0].set_ylabel(r'KHGX $Z_e$ [dBZ]')
    axarr[0,0].set_xlabel(r'LDQUANTS M1-Site Derived $Z_e$ [dBZ]')
    # Plot the data.
    axarr[0,1].hist2d(data_m1['VD_Z'],data_m1['rhi_data'],bins=nbin,
                   range=[[brange[0],brange[1]],[brange[2],brange[3]]],
                   norm=mpl.colors.LogNorm(),cmap=mpl.cm.jet)
    axarr[0,1].plot(np.arange(brange[0] - 10,brange[1] + 10,1),
                 np.arange(brange[0] - 10,brange[1] + 10,1),'k')
    axarr[0,1].set_xlabel(r'VDISQUANTS M1-Site Derived $Z_e$ [dBZ]')
    axarr[0,1].set_ylabel(r'KHGX $Z_e$ [dBZ]')
    # Plot the data.
    axarr[1,0].hist2d(data_s1['LD_Z'],data_s1['rhi_data'],bins=nbin,
                   range=[[brange[0],brange[1]],[brange[2],brange[3]]],
                   norm=mpl.colors.LogNorm(),cmap=mpl.cm.jet)
    axarr[1,0].plot(np.arange(brange[0] - 10,brange[1] + 10,1),
                 np.arange(brange[0] - 10,brange[1] + 10,1),'k')
    axarr[1,0].set_xlabel(r'LDQUANTS S1-Site Derived $Z_e$ [dBZ]')
    axarr[1,0].set_ylabel(r'KHGX $Z_e$ [dBZ]')
    # Hide the blank subplot
    axarr[1,1].set_visible(False)
    # plot the title
    plt.title(data_m1['date'][0]
              + ' Equivalent Radar Reflectivity Factor Comparison')
    # Save the figure.
    nout = ('ldquantsRadar_' + radar.metadata['instrument_name'] + '_'
            + radar.time['units'].split(' ')[-1].split('T')[0]
            + 'fullComparison.png')
    plt.savefig(nout, dpi=100)
    # Close the figure.
    plt.close(fig)

#-------------------------
# II) Input
#-------------------------

# Define the starting time of the code.
t0 = time.time()

# Check for options within the system arguments
for param in sys.argv:
    if param.startswith('-h'):
        help_message()
        sys.exit()

# Check to make sure there were three input files.
if len(sys.argv) < 4:
    help_message()
    sys.exit()
else:
    m1file  = sys.argv[-3]
    s1file  = sys.argv[-2]
    vidfile = sys.argv[-1]

# add checks to make sure files are correct
if not m1file.split('.')[0].endswith('M1'):
    print('ERROR: LDQUANTS M1 Site file incorrectly entered')
    sys.exit(1)
if m1file.split('.')[0][3:-2] != 'ldquants':
    print('ERROR: LDQUANTS M1 Site file incorrectly entered')
    sys.exit(1)
if not s1file.split('.')[0].endswith('S1'):
    print('ERROR: LDQUANTS S1 Site file incorrectly entered')
    sys.exit(1)
if s1file.split('.')[0][3:-2] != 'ldquants':
    print("ERROR: LDQUANTS (S1) Site file incorrectly entered")
    sys.exit(1)
if not vidfile.split('.')[0].endswith('M1'):
    print('ERROR: VDISQUANTS S1 Site file incorrectly entered')
    sys.exit(1)
if vidfile.split('.')[0][3:-2] != 'vdisquants':
    print('ERROR: VDISQUANTS S1 Site file incorrectly entered')
    sys.exit(1)

# Define the date from the LDQUANTS file, create datetime object
m1_date = m1file.split('.')[2]
ndate = datetime.datetime(int(m1_date[0:4]),int(m1_date[4:6]),
                          int(m1_date[6:8]))

# Read the file. Define a dataset variable.
nc_M = netCDF4.Dataset(m1file,'r')
nc_S = netCDF4.Dataset(s1file,'r')
nc_V = netCDF4.Dataset(vidfile,'r')

#-------------------------------------------------------------------
# III) Establish nexradaws connection and pull radar data
#-------------------------------------------------------------------

# Connect to the Cloud - nexradaws is pup installable from anaconda
conn = nexradaws.NexradAwsInterface()

# Define a NEXRAD site for this comparison.
# As this is designed for TRACER, currently just hard-coding in
#    KHGX site.
NSITE = 'KHGX'

# Call the pull_allscans function to download the data.
downloads = pull_allscans(conn,NSITE,ndate)

#-------------------------------------------------------------------
# IV) Iterate over successful downloads.
#     Calculate distance and azimuth to target
#-------------------------------------------------------------------

# initiate a dictionary to hold the striped out data
# from each scan.
ndata_m1 = {'date':[],'time':[], 'xGate':[],'yGate':[],'zGate':[],
            'tarDis':[],'tarAzi':[],'rhi_data':[],'rain_data':[],
            'LD_rain':[],'LD_Z':[], 'VD_rain': [], 'VD_Z': [],
           }

ndata_s1 = {'date':[],'time':[], 'xGate':[],'yGate':[],'zGate':[],
            'tarDis':[],'tarAzi':[],'rhi_data':[],'rain_data':[],
            'LD_rain':[],'LD_Z':[]
           }

# iterate over the sucessfull downloads
for loc_file in downloads.iter_success():
    # check to make sure not to use the MDM files
    if loc_file.filepath[-3:] != 'MDM':
        # print the file
        print(loc_file)
        # Read in the file.
        radar = pyart.io.read(loc_file.filepath)
        # Create the Radar Display
        display = pyart.graph.RadarDisplay(radar)
        # Calculate estimated rainfall rate from reflectivty
        rain = pyart.retrieve.qpe.est_rain_rate_z(radar)
        # Add the estimated rainfall rate back into the radar object
        radar.add_field('est_rainfall_rate',rain)
        # call the sphere_distance function
        dis_M1 = sphere_distance(radar.latitude['data'][:],
                                    nc_M.variables['lat'][:],
                                    radar.longitude['data'],
                                    nc_M.variables['lon'][:])

        dis_S1 = sphere_distance(radar.latitude['data'][:],
                                    nc_S.variables['lat'][:],
                                    radar.longitude['data'],
                                    nc_S.variables['lon'][:])
        # call the for_azimuth function
        azim_M1  = for_azimuth(radar.latitude['data'][:],
                                nc_M.variables['lat'][:],
                                radar.longitude['data'],
                                nc_M.variables['lon'][:])

        azim_S1  = for_azimuth(radar.latitude['data'][:],
                                nc_S.variables['lat'][:],
                                radar.longitude['data'],
                                nc_S.variables['lon'][:])
        # Find the reflectivty data for the azimuth over the target,
        # find gates for that location
        # Note: x,y,z gates are in km
        rhi_data_m1,rhix_m1,rhiy_m1,rhiz_m1 = (pyart.graph.RadarDisplay(radar)
                                             ._get_azimuth_rhi_data_x_y_z(
                                             'reflectivity',azim_M1,
                                             edges=True,mask_tuple=None,
                                             gatefilter=None,
                                             filter_transitions=True))

        rhi_data_s1,rhix_gci,rhiy_s1,rhiz_s1 = (pyart.graph.RadarDisplay(radar)
                                             ._get_azimuth_rhi_data_x_y_z(
                                             'reflectivity',azim_S1,
                                             edges=True,mask_tuple=None,
                                             gatefilter=None,
                                             filter_transitions=True))
        # Find the estimated rainfall data for the azimuth over the target,
        # find gates for that location
        # Note: x,y,z gates are in km
        rain_data_m1,rainx_m1,rainy_m1,rainz_m1 = (pyart.graph.RadarDisplay(
                                                 radar).
                                                 _get_azimuth_rhi_data_x_y_z(
                                                 'est_rainfall_rate',
                                                 azim_M1,edges=True,
                                                 mask_tuple=None,
                                                 gatefilter=None,
                                                 filter_transitions=True))

        rain_data_s1,rainx_s1,rainy_s1,rainz_s1 = (pyart.graph.RadarDisplay(
                                                 radar).
                                                 _get_azimuth_rhi_data_x_y_z(
                                                 'est_rainfall_rate',
                                                 azim_S1,
                                                 edges=True,mask_tuple=None,
                                                 gatefilter=None,
                                                 filter_transitions=True))
        # Calculate distance from the x,y coordinates to target
        rhidis_m1  = np.sqrt((rhix_m1**2) + (rhiy_m1**2))*np.sign(rhiy_m1)
        rhidis_s1  = np.sqrt((rhix_gci**2) + (rhiy_s1**2))*np.sign(rhiy_s1)
        # calculate the gate distance to the target
        # mask all data in RHI besides vertical column above the target
        # need to iterate over every sweep
        colrain_m1 = np.ma.zeros(rhi_data_m1.shape[0])
        for i in range(rhi_data_m1.shape[0]):
	    # Find closest gate to the target
            tar_gate = np.where(abs(rhidis_m1[i,:] - (dis_M1/1000.))
                               == min(abs(rhidis_m1[i,:]
                               - (dis_M1/1000.))))
            # Create a new bool array to hold the new mask
            mask = np.ones(rhi_data_m1[i,:].shape,bool)
            mask[:] = True
            mask[tar_gate[0]] = False
            # Append each gates estimated rain rate to the zeros array.
            if rain_data_m1[i,tar_gate[0]].mask is False:
                colrain_m1[i] = rain_data_m1[i,tar_gate[0]]
            else:
                colrain_m1[i] = np.ma.masked
            # Mask every value but the gate over the target
            rhi_data_m1[i,:].mask  = mask
            rain_data_m1[i,:].mask = mask

        # calculate the gate distance to the target
        # mask all data in RHI besides vertical column above the target
        # need to iterate over every sweep
        colrain_s1 = np.ma.zeros(rhi_data_s1.shape[0])
        for i in range(rhi_data_s1.shape[0]):
	    # Find closest gate to the target
            tar_gate = np.where(abs(rhidis_s1[i,:] - (dis_S1/1000.))
                               == min(abs(rhidis_s1[i,:]
                               - (dis_S1/1000.))))
            # Create a new bool array to hold the new mask
            mask = np.ones(rhi_data_s1[i,:].shape,bool)
            mask[:] = True
            mask[tar_gate[0]] = False
            # Append each gates estimated rain rate to the zeros array.
            if rain_data_s1[i,tar_gate[0]].mask is False:
                colrain_s1[i] = rain_data_s1[i,tar_gate[0]]
            else:
                colrain_s1[i] = np.ma.masked
            # Mask every value but the gate over the target
            rhi_data_s1[i,:].mask  = mask
            rain_data_s1[i,:].mask = mask

        # Note:
        #       Time sync against daily file of the ldquants data.
        #       Time sync is going to take more thought out approach.
        #       Just trying to get a quick comparison with the disdrometer

        # I'm sure there's a better way to do this with netCDF4
        # datetime objects
        # Convert radar time to seconds from midnight
        # note: radar time units are seconds from start of the scan.
        begin_time     = radar.time['units'].split(' ')[2].split('T')[-1][0:-1]
        sfm_begin_time = (int(begin_time.split(':')[0])*3600.
                       + int(begin_time.split(':')[1])*60.
                       + int(begin_time.split(':')[2]))
        sfm_end_time   = sfm_begin_time+radar.time['data'][-1]

        # find the ldqaunt times that match the radar sfm.
        # note: ldquant times are already in sfm from the start of the day
        sync_x      = np.where(nc_M.variables['time'][:] > sfm_begin_time)
        sync_y      = np.where(nc_M.variables['time'][:] < sfm_end_time)
        time_sync_M1 = np.intersect1d(sync_x,sync_y)

        # find the ldqaunt times that match the radar sfm.
        # note: ldquant times are already in sfm from the start of the day
        sync_x_S1    = np.where(nc_S.variables['time'][:] > sfm_begin_time)
        sync_y_S1    = np.where(nc_S.variables['time'][:] < sfm_end_time)
        time_sync_S1 = np.intersect1d(sync_x_S1,sync_y_S1)

        # find the ldqaunt times that match the radar sfm.
        # note: ldquant times are already in sfm from the start of the day
        sync_x_V1    = np.where(nc_V.variables['time'][:] > sfm_begin_time)
        sync_y_V1    = np.where(nc_V.variables['time'][:] < sfm_end_time)
        time_sync_V1 = np.intersect1d(sync_x_V1,sync_y_S1)

        # Determine where x,y,z gates and rhi/rain data is valid
        mask_x,mask_y = np.where(rhi_data_m1.mask is False)

        # fill the ndata dictionary for further analysis
        ##for i in range(len(mask_x)):
        ndata_m1['date'].append(radar.time['units'].split(' ')[2][0:10])
        ndata_m1['time'].append(sfm_begin_time)
        ndata_m1['xGate'].append(rhix_m1[mask_x,mask_y][0])
        ndata_m1['yGate'].append(rhiy_m1[mask_x,mask_y][0])
        ndata_m1['zGate'].append(rhiz_m1[mask_x,mask_y][0])
        ndata_m1['tarDis'].append(dis_M1[0]/1000.)
        ndata_m1['tarAzi'].append(azim_M1[0])
        ndata_m1['rhi_data'].append(rhi_data_m1[mask_x,mask_y][0])
        ndata_m1['rain_data'].append(rain_data_m1[mask_x,mask_y][0])
        if np.ma.sum(nc_M.variables['rain_rate'][time_sync_M1] is not
                np.ma.masked):
            ndata_m1['LD_rain'].append(np.ma.around(
                                       np.ma.sum(nc_M.variables[
                                       'rain_rate'][time_sync_M1]),2))
        else:
            ndata_m1['LD_rain'].append(0.0)

        if (np.ma.average(nc_M.variables['reflectivity_factor_sband20c']
                [time_sync_M1]) is not np.ma.masked):
            ndata_m1['LD_Z'].append(np.ma.around(
                                    np.ma.average(nc_M.variables[
                                    'reflectivity_factor_sband20c']
                                    [time_sync_M1]),2))
        else:
            ndata_m1['LD_Z'].append(-75)

        if (np.ma.sum(nc_V.variables['rain_rate'][time_sync_V1]) is not
                np.ma.masked):
            ndata_m1['VD_rain'].append(np.ma.around(np.ma.sum(
                                       nc_V.variables['rain_rate']
                                       [time_sync_V1]),2))
        else:
            ndata_m1['VD_rain'].append(0.0)

        if (np.ma.average(nc_V.variables['reflectivity_factor_sband20c']
                [time_sync_V1]) is not np.ma.masked):
            (ndata_m1['VD_Z'].append(np.ma.around(np.ma.average(
                                     nc_V.variables[
                                     'reflectivity_factor_sband20c']
                                     [time_sync_V1]),2)))
        else:
            ndata_m1['VD_Z'].append(-75)

        ndata_s1['date'].append(radar.time['units'].split(' ')[2][0:10])
        ndata_s1['time'].append(sfm_begin_time)
        ndata_s1['xGate'].append(rhix_gci[mask_x,mask_y][0])
        ndata_s1['yGate'].append(rhiy_s1[mask_x,mask_y][0])
        ndata_s1['zGate'].append(rhiz_s1[mask_x,mask_y][0])
        ndata_s1['tarDis'].append(dis_S1[0]/1000.)
        ndata_s1['tarAzi'].append(azim_S1[0])
        ndata_s1['rhi_data'].append(rhi_data_s1[mask_x,mask_y][0])
        ndata_s1['rain_data'].append(rain_data_s1[mask_x,mask_y][0])
        if (np.ma.sum(nc_S.variables['rain_rate'][time_sync_S1]) is not
                np.ma.masked):
            (ndata_s1['LD_rain'].append(
                                        np.ma.around(np.ma.sum(
                                        nc_S.variables['rain_rate']
                                        [time_sync_S1]),2)))
        else:
            ndata_s1['LD_rain'].append(0.0)

        if (np.ma.average(nc_S.variables['reflectivity_factor_sband20c']
                [time_sync_S1]) is not np.ma.masked):
            (ndata_s1['LD_Z'].append(np.ma.around(np.ma.average(nc_S.
                                     variables['reflectivity_factor_sband20c']
                                     [time_sync_S1]),2)))
        else:
            ndata_s1['LD_Z'].append(-75)

        # Call the plotting function
        plot_ppirhi(display,radar,nc_M,time_sync_M1,azim_M1,dis_M1,colrain_m1)

#-----------------------------------------------
# V) Plot the 2D-Histrogram for the entire date
#----------------------------------------------

# Define the range of values to be binned
his_range = [-20, 70, -20, 70]
NBINS = 90

# call the plot_his function
plot_his(ndata_m1,ndata_s1,his_range,NBINS)

#--------------------
# VI) END OF PROGRAM
#--------------------

# Close the netCDF file.
nc_M.close()
nc_S.close()
nc_V.close()

# define the ending time of the program for testing purposes
t1 = time.time()

# print out run time
t3 = (t1-t0)/60.
print("Run Time: ", t3, ' min')
