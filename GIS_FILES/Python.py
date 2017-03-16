#=============================================================
# Import the system modules
import arcpy
import os
from arcpy import env

# Set current working directory (cwd)
cwd = os.getcwd()
# ===============================================================

# Input data is in the NAD 1983 coordinate system
# Create a spatial referene object for the output coordinate system: NAD 83
# California Zone 3
out_coordinate_system = arcpy.SpatialReference('NAD 1983 (2011) StatePlane California III FIPS 0403 (Meters)')

# Project the counties shapefile
arcpy.Project_management(os.path.join(cwd, "Data/tl_2015_us_county/tl_2015_us_county.shp"), os.path.join(cwd, "Pre_analyzed_files/counties_Nad83z3.shp"), out_coordinate_system)

# Project the states shapefile
arcpy.Project_management(os.path.join(cwd, "Data/tl_2015_us_state/tl_2015_us_state.shp"), os.path.join(cwd, "Pre_analyzed_files/state_Nad83z3.shp"), out_coordinate_system)

# Project the Sub-counties shapefile
arcpy.Project_management(os.path.join(cwd, "Data/tl_2015_06_cousub/tl_2015_06_cousub.shp"), os.path.join(cwd, "Pre_analyzed_files/Subcounties_Nad83z3.shp"), out_coordinate_system)

# Merge the water shapefiles
arcpy.Merge_management([os.path.join(cwd, "Data/tl_2015_06001_areawater/tl_2015_06001_areawater.shp"),os.path.join(cwd, "Data/tl_2015_06013_areawater/tl_2015_06013_areawater.shp"),
			os.path.join(cwd, "Data/tl_2015_06041_areawater/tl_2015_06041_areawater.shp"),os.path.join(cwd, "Data/tl_2015_06075_areawater/tl_2015_06075_areawater.shp"),
			os.path.join(cwd, "Data/tl_2015_06081_areawater/tl_2015_06081_areawater.shp")], os.path.join(cwd, "Pre_analyzed_files/Water_area.shp"))

# Project the merged water shapefile
arcpy.Project_management(os.path.join(cwd, "Pre_analyzed_files/Water_area.shp"), os.path.join(cwd, "Pre_analyzed_files/Water_Nad83z3.shp"), out_coordinate_system)

# ===============================================================

# SEE TUTORIAL ON ANALYZING THE RASTER FILES

# ===============================================================
# Project the new shapefiles for the output coordinate system: NAD 83
# California Zone 3
# ===============================================================

arcpy.Project_management(os.path.join(cwd, "Pre_analyzed_files/NED_Water.shp"), os.path.join(cwd, "Pre_analyzed_files/NED_Water_nad83z3.shp"), out_coordinate_system)

arcpy.Project_management(os.path.join(cwd, "Pre_analyzed_files/NED_stormsurge.shp"), os.path.join(cwd, "Pre_analyzed_files/NED_ss_Nad83z3.shp"), out_coordinate_system)

arcpy.Project_management(os.path.join(cwd, "Pre_analyzed_files/NED_mslr.shp"), os.path.join(cwd, "Pre_analyzed_files/NED_mslr_Nad83z3.shp"), out_coordinate_system)

arcpy.Project_management(os.path.join(cwd, "Pre_analyzed_files/NED_h09.shp"), os.path.join(cwd, "Pre_analyzed_files/NED_h09_Nad83z3.shp"), out_coordinate_system)

arcpy.Project_management(os.path.join(cwd, "Pre_analyzed_files/NED_uslr.shp"), os.path.join(cwd, "Pre_analyzed_files/NED_uslr_Nad83z3.shp"), out_coordinate_system)

# ===============================================================
# Erase the known water area from the flood risk areas:
# ===============================================================

# Baseline storm surge
arcpy.Erase_analysis(os.path.join(cwd, "Pre_analyzed_files/NED_ss_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Water_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/ss_Nad83z3.shp"), "")

# Mean sea-level rise
arcpy.Erase_analysis(os.path.join(cwd, "Pre_analyzed_files/NED_mslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Water_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/mslr_Nad83z3.shp"), "")

# Estimate from Heberger et al. (2009)
arcpy.Erase_analysis(os.path.join(cwd, "Pre_analyzed_files/NED_h09_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Water_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/h09_Nad83z3.shp"), "")

# Accounting for sea-level rise uncertainty
arcpy.Erase_analysis(os.path.join(cwd, "Pre_analyzed_files/NED_uslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Water_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/uslr_Nad83z3.shp"), "")

# ===============================================================
# Select and write the selected features to a new feature class
# for San Francisco #635, Oakland #319, and Alameda #320
# ===============================================================
# San Francisco #635
arcpy.SelectLayerByAttribute_management("counties_Nad83z3", "NEW_SELECTION", '"FID" = 635 ')
arcpy.CopyFeatures_management("counties_Nad83z3", "SanFrancisco")

# Clear the selection
arcpy.SelectLayerByAttribute_management("counties_Nad83z3", "CLEAR_SELECTION")

# Oakland #319
arcpy.SelectLayerByAttribute_management("Subcounties_Nad83z3", "NEW_SELECTION", '"FID" = 319 ')
arcpy.CopyFeatures_management("Subcounties_Nad83z3" , "Oakland")

# Clear the selection
arcpy.SelectLayerByAttribute_management("Subcounties_Nad83z3", "CLEAR_SELECTION")

# Alameda #320
arcpy.SelectLayerByAttribute_management("Subcounties_Nad83z3", "NEW_SELECTION", '"FID" = 320 ')
arcpy.CopyFeatures_management("Subcounties_Nad83z3", "Alameda")

# Clear the selection
arcpy.SelectLayerByAttribute_management("Subcounties_Nad83z3", "CLEAR_SELECTION")

# ===============================================================
# Clip the flood area polygons with the San Francisco, Oakland, and Alameda
# counties/ sub-counties
# ===============================================================

arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/ss_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/SanFrancisco.shp"), os.path.join(cwd, "Pre_analyzed_files/ss_SanFrancisco.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/mslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/SanFrancisco.shp"), os.path.join(cwd, "Pre_analyzed_files/mslr_SanFrancisco.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/h09_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/SanFrancisco.shp"), os.path.join(cwd, "Pre_analyzed_files/h09_SanFrancisco.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/uslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/SanFrancisco.shp"), os.path.join(cwd, "Pre_analyzed_files/uslr_SanFrancisco.shp"))

arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/ss_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Oakland.shp"), os.path.join(cwd, "Pre_analyzed_files/ss_Oakland.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/mslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Oakland.shp"), os.path.join(cwd, "Pre_analyzed_files/mslr_Oakland.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/h09_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Oakland.shp"), os.path.join(cwd, "Pre_analyzed_files/h09_Oakland.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/uslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Oakland.shp"), os.path.join(cwd, "Pre_analyzed_files/uslr_Oakland.shp"))

arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/ss_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Alameda.shp"), os.path.join(cwd, "Pre_analyzed_files/ss_Alameda.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/mslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Alameda.shp"), os.path.join(cwd, "Pre_analyzed_files/mslr_Alameda.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/h09_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Alameda.shp"), os.path.join(cwd, "Pre_analyzed_files/h09_Alameda.shp"))
arcpy.Clip_analysis(os.path.join(cwd, "Pre_analyzed_files/uslr_Nad83z3.shp"), os.path.join(cwd, "Pre_analyzed_files/Alameda.shp"), os.path.join(cwd, "Pre_analyzed_files/uslr_Alameda.shp"))

# ===============================================================
# Create an additional field to quantify the area in km^2 flooded
# in each county/ sub-county
# ===============================================================

arcpy.AddField_management("ss_SanFrancisco", "AREASQKM", "DOUBLE")
arcpy.AddField_management("mslr_SanFrancisco", "AREASQKM", "DOUBLE")
arcpy.AddField_management("h09_SanFrancisco", "AREASQKM", "DOUBLE")
arcpy.AddField_management("uslr_SanFrancisco", "AREASQKM", "DOUBLE")

arcpy.AddField_management("ss_Oakland", "AREASQKM", "DOUBLE")
arcpy.AddField_management("mslr_Oakland", "AREASQKM", "DOUBLE")
arcpy.AddField_management("h09_Oakland", "AREASQKM", "DOUBLE")
arcpy.AddField_management("uslr_Oakland", "AREASQKM", "DOUBLE")

arcpy.AddField_management("ss_Alameda", "AREASQKM", "DOUBLE")
arcpy.AddField_management("mslr_Alameda", "AREASQKM", "DOUBLE")
arcpy.AddField_management("h09_Alameda", "AREASQKM", "DOUBLE")
arcpy.AddField_management("uslr_Alameda", "AREASQKM", "DOUBLE")

# ===============================================================
