__author__ = 'jwely'

import os
import shutil
import textio
import arcpy
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True


def copyfile(path, dest_dir, workspace = ""):
    """
    path      the full filepath to a file
    dest_dir  destination for copy
    returns   the full filepath of the new destination

    removes the workspace from the filepath to give a
    workspace relative filepath.
    """

    if os.path.isfile(path):

        head, tail = os.path.split(path)
        destination = os.path.join(workspace, dest_dir, tail)

        if not os.path.isfile(destination):
            shutil.copy(path, destination)
            print("Added {0}".format(destination))
        else:
            print("Found {0}".format(destination))

        return destination.replace(workspace + "\\", "")

    else:
        print("{0} is an invalid filepath!".format(path))
        return None


def clip_inputs(raster, shapefile, outpath):
    """
    uses arcpy to clip a raster to an input shapefile, used for reducing the calculation area.
    """

    print("Clipped and added {0}".format(outpath))
    arcpy.Clip_management(raster, "#", outpath, shapefile, "ClippingGeometry")
    out = arcpy.sa.Float(arcpy.sa.ExtractByMask(outpath, shapefile))
    out.save(outpath)
    return outpath


def prepare_metric_env(workspace, landsat_band2, landsat_band3, landsat_band4, landsat_band5,
                       landsat_band6, landsat_band7, landsat_band10, landsat_band11, landsat_metapath,
                       dem_path, hot_shape_path, cold_shape_path, wx_filepath, testflag = None,
                       recalc = None, crop_type = None, wx_elev = None,
                       wx_zom = None, LE_cold_cal_factor = None, mountains = None, L_green_fac = None,
                       clip_extent = None):
    """
    Saves a config file with all required attributes of this metric model, also copies
    the source data files into the template structure for good record keeping. all stored
    filepaths are relative to "workspace"

    workspace:              fresh folder to populate with the many metric model files and parameters
    landsat_filepath_list:  list of landsat band tiff filepaths, MUST be in order [2,3,4,5,6,7,10,11]
    landsat_metapath:       filepath to the landsat metadata
    dem_path:               filepath to the digital elevation model
    hot_shape_path:         filepath to the shapefile outlining the hot pixels
    cold_shape_path:        filepath to the shapefile outlining the cold pixels
    wx_filepath:            filepath to the weather data. downloaded by following NCDC download tutorial
    testflag:               parameter string to regulate intermediate output file generation
                                possible values include "ALL" , "LIMITED", "ET_ONLY"
    recalc:                 forces all calculations to be re-performed
    crop_type:              reference crop type, either "alfalfa" or "corn"
    wx_elev                 the elevation of the weather station used
    wx_zom                  estimate of the "zom" term at weather station (temporary)
    LE_cold_cal_factor      used to calibrate LE terms. This should probably always be = 1.05
    mountains               will be either True or False, defaults to False.
    L_green_fac             reference L factor for calculating SAVI (0.1 for Idaho, 0.5 for NorthCarolina)
    """

    # set default values if they are still None
    timezone = 0    # timezone value should always be zero, since all data is already in UTC

    if testflag is None:
        testflag = "ALL"

    if recalc is None:
        recalc = True

    if crop_type is None:
        crop_type = "alfalfa"

    if wx_elev is None:
        wx_elev = 1

    if wx_zom is None:
        wx_zom = 0.010

    if LE_cold_cal_factor is None:
        LE_cold_cal_factor = 1.05

    if mountains is None:
        mountains = False

    if L_green_fac is None:
        L_green_fac = 0.5


    # first copy the template "empty metric model" structure
    if not os.path.isdir(workspace):
        live_path = os.path.realpath(__file__)
        head, tail = os.path.split(live_path)
        template_path = os.path.join(head, "Empty_Metric_Model")
        shutil.copytree(template_path, workspace)
    else:
        raise Exception("input workspace must be a directory that does not already exist! one will be created here!")


    # set other inferred attributes of the working directory structure
    out_dir        = "output"
    middle_dir     = "intermediate_calculations"
    ref_pixel_dir  = "input_ref_pixels"
    landsat_dir    = "input_landsat"
    weather_dir    = "input_weather"
    dem_dir        = "input_dem"
    geodatabase    = "scratch"


    # move the landsat data and set the new path attributes
    landsat_metapath = copyfile(landsat_metapath, landsat_dir, workspace)

    bands = [landsat_band2, landsat_band3, landsat_band4, landsat_band5,
                 landsat_band6, landsat_band7, landsat_band10, landsat_band11]

    for i,band in enumerate(bands):

        if clip_extent is None:
            if os.path.exists(band + ".ovr"):
                copyfile(band + ".ovr", landsat_dir, workspace)
            bands[i] = copyfile(band, landsat_dir, workspace)

        else:
            copyfile(band + ".ovr", landsat_dir, workspace)
            bands[i]= os.path.join(workspace, landsat_dir, os.path.basename(band))
            clip_inputs(band, clip_extent, bands[i])


    # move the DEM and associated files
    for demfile in [dem_path + ext for ext in [".ovr", ".aux.xml", ".xml"]]:
        if os.path.exists(demfile):
            copyfile(dem_path + ".ovr", dem_dir, workspace)
    if os.path.exists(dem_path.replace(".tif", ".tfw")):
        copyfile(dem_path.replace(".tif", ".tfw"), dem_dir, workspace)

    if clip_extent is None:
        dem_path = copyfile(dem_path, dem_dir, workspace)
    else:
        outpath = os.path.join(workspace, dem_dir, os.path.basename(dem_path))
        dem_path = clip_inputs(dem_path, clip_extent, outpath)


    # moves the shapefiles for hot and cold pixels, and clip extent.
    extensions = [".cpg", ".dbf", ".prj", ".sbn", ".sbx", ".shx", ".shp"]

    for extension in extensions:
        copyfile(hot_shape_path.replace(".shp", extension), ref_pixel_dir, workspace)
        copyfile(cold_shape_path.replace(".shp", extension), ref_pixel_dir, workspace)
        if clip_extent is not None:
            copyfile(clip_extent.replace(".shp", extension), dem_dir, workspace)

    #hot_shape_path = copyfile(hot_shape_path, ref_pixel_dir, workspace)
    #cold_shape_path = copyfile(cold_shape_path, ref_pixel_dir, workspace)
    #clip_extent = copyfile(clip_extent, dem_dir, workspace)

    # move the weather data
    wx_filepath = copyfile(wx_filepath, weather_dir, workspace)


    # create the config file
    config_filepath = "config.txt"

    cdict =  {"config_filepath" : config_filepath,
              "metric_workspace": workspace,
              "out_dir"         : out_dir,
              "middle_dir"      : middle_dir,
              "ref_pixel_dir"   : ref_pixel_dir,
              "landsat_dir"     : landsat_dir,
              "weather_dir"     : weather_dir,
              "dem_dir"         : dem_dir,
              "geodatabase"     : geodatabase,

              "landsat_meta"    : landsat_metapath,
              "landsat_band2"   : bands[0],
              "landsat_band3"   : bands[1],
              "landsat_band4"   : bands[2],
              "landsat_band5"   : bands[3],
              "landsat_band6"   : bands[4],
              "landsat_band7"   : bands[5],
              "landsat_band10"  : bands[6],
              "landsat_band11"  : bands[7],

              "dem_path"        : dem_path,
              "hot_shp_path"    : hot_shape_path,
              "cold_shp_path"   : cold_shape_path,
              "weather_path"    : wx_filepath,

              "crop_type"       : crop_type,
              "timezone"        : timezone,
              "wx_elev"         : wx_elev,
              "wx_zom"          : wx_zom,
              "LE_ref"          : LE_cold_cal_factor,
              "mountains"       : mountains,
              "L_green_fac"     : L_green_fac,
              "testflag"        : testflag,
              "recalc"          : recalc,
              "clip_extent"     : clip_extent}

    config = textio.ioconfig()
    config.add_param(cdict)
    config.write(os.path.join(workspace, config_filepath))
    config.read(os.path.join(workspace, config_filepath))
    print("configuration file has been saved at {0}".format(os.path.join(workspace, config_filepath)))

    return config_filepath


# testing area
if __name__ == "__main__":
    pass
    # quick dirty way to build config for a bunch of stuff
