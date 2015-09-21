# arcpy imports
import arcpy
arcpy.env.overwriteOutput = True

from METRIC import metric_py
from METRIC import prepare_metric_env


__author__ = ["Kent Sparrow",
              "Jamie Vanderheiden",
              "Nathan Quian",
              "Jeffry Ely, jeff.ely.08@gmail.com"]


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""

        self.label = "DEVELOP_METRIC" 
        self.alias = "DEVELOP_METRIC"

        self.tools = [METRIC_configure, METRIC_run]
        return


class METRIC_configure(object):
    
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "METRIC_configure"
        self.description = "To configure the metric model"
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""

        
        workspace = arcpy.Parameter(
            displayName =   "Workspace to store intermediates and outputs (full path to a NEW directory)",
            name =          "workspace",
            datatype =      "GPString",
            parameterType = "Required",
            direction =     "Input")
        
        l8_meta = arcpy.Parameter(
            displayName =   "Landsat 8 Metadata filepath",
            name =          "l8_meta",
            datatype =      "DEFile",
            parameterType = "Required",
            direction =     "Input")
        l8_meta.filter.list = ['txt']

        band2 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 2 filepath (tif)",
            name =          "l8_b2",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band3 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 3 filepath (tif)",
            name =          "l8_b3",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band4 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 4 filepath (tif)",
            name =          "l8_b4",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band5 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 5 filepath (tif)",
            name =          "l8_b5",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band6 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 6 filepath (tif)",
            name =          "l8_b6",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band7 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 7 filepath (tif)",
            name =          "l8_b7",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band10 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 10 filepath (tif)",
            name =          "l8_b10",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        band11 = arcpy.Parameter(
            displayName =   "Landsat 8 Band 11 filepath (tif)",
            name =          "l8_b11",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input",
            category =      "Landsat bands")

        dem = arcpy.Parameter(
            displayName =   "Digital Elevation Model raster (tif)",
            name =          "DEM",
            datatype =      "GPRasterLayer",
            parameterType = "Required",
            direction =     "Input")

        hot_pix = arcpy.Parameter(
            displayName =   "shapefile of hot pixels (non-irrigated)",
            name =          "hot_pix",
            datatype =      "DEShapefile",
            parameterType = "Required",
            direction =     "Input")
        
        cold_pix = arcpy.Parameter(
            displayName =   "shapefile of cold pixels (irrigated)",
            name =          "cold_pix",
            datatype =      "DEShapefile",
            parameterType = "Required",
            direction =     "Input")

        wx_file = arcpy.Parameter(
            displayName =   "Weather data text file (txt)",
            name =          "wx_file",
            datatype =      "DEFile",
            parameterType = "Required",
            direction =     "Input")
        wx_file.filter.list = ['txt']
        
        crop_type = arcpy.Parameter(
            displayName =   "Reference crop type",
            name =          "crop",
            datatype =      "GPString",
            parameterType = "Required",
            direction =     "Input")
        crop_type.filter.list = ["alfalfa", "grass"]

        wx_elev = arcpy.Parameter(
            displayName =   "manual elevation of weather station (in meters)",
            name =          "wx_elevation",
            datatype =      "GPDouble",
            parameterType = "Optional",
            direction =     "Input")
            
        wx_zom = arcpy.Parameter(
            displayName =   "manual roughness length (zom) estimate at weather station location ",
            name =          "wx_zom",
            datatype =      "GPDouble",
            parameterType = "Optional",
            direction =     "Input")
            
        LE_ref = arcpy.Parameter(
            displayName =   "manual LE_cold calibration factor Usually 1.05",
            name =          "LE_cold_cal_factor",
            datatype =      "GPDouble",
            parameterType = "Optional",
            direction =     "Input")
        
        mountains = arcpy.Parameter(
            displayName =   "Mountainous terrain",
            name =          "mounts",
            datatype =      "GPBoolean",
            parameterType = "Optional",
            direction =     "Input")
        
        saveflag = arcpy.Parameter(
            displayName =   "Save_flag (decides which variables to save to hard disk)",
            name =          "saveflag",
            datatype =      "GPString",
            parameterType = "Required",
            direction =     "Input")
        saveflag.filter.list = ["ALL", "LIMITED", "ET-ONLY"]

        recalc = arcpy.Parameter(
            displayName =   "Force Recalculation",
            name =          "recalc",
            datatype =      "GPBoolean",
            parameterType = "Required",
            direction =     "Input")


        l_green_fac = arcpy.Parameter(
            displayName =   "reference SAVI adjust",
            name =          "l_green_fac",
            datatype =      "GPDouble",
            parameterType = "Optional",
            direction =     "Input")
        
        params = [l8_meta, band2, band3, band4, band5, band6, band7, band10, band11, dem,
                  hot_pix, cold_pix, wx_file, crop_type, workspace, wx_elev, wx_zom, LE_ref,
                  mountains, saveflag, recalc, l_green_fac]

        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True
    

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return
    

    def updateMessages(self, parameters): 
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return
    

    def execute(self, parameters, messages):
        """ Run tool sourcecode by calling metric.py """

        p = parameters
        ps = []
        for param in parameters:
            typestring = str(type(param.value))
            if "geop" in typestring:
                ps.append(str(param.value))
            else:
                ps.append(param.value)

        arcpy.AddMessage(ps)

        # build a config file with these inputs
        config_filepath = prepare_metric_env(ps[14], ps[2], ps[2], ps[3], ps[4], ps[5], ps[6], ps[7], ps[8], ps[0],
                                             ps[9], ps[10], ps[11], ps[12], ps[19], ps[20], ps[13], ps[15],
                                             ps[16], ps[17], ps[18], ps[21])

        # print some information to the screen
        for param in p:
            print("{0} : {1}".format(str(param.name).ljust(12," "), param.value))
            arcpy.AddMessage("{0} : {1}".format(str(param.name).ljust(12," "), param.value))

        print("Configuration file created at {0}".format(config_filepath))
        arcpy.AddMessage("Configuration file created at {0}".format(config_filepath))

        return



class METRIC_run(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "METRIC_run"
        self.description = "To run the METRIC model from existing configuration file"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        configpath = arcpy.Parameter(
            displayName =   "METRIC config filepath",
            name =          "configpath",
            datatype =      "DEFile",
            parameterType = "Required",
            direction =     "Input")
        configpath.filter.list = ['txt']

        params = [configpath]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """ Run tool sourcecode by calling metric.py """

        metric_py.run(str(parameters[0].value))
        return

# testing area
