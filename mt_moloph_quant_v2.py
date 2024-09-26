#// Mitochodria molphology quantification, jython version
#// author: Masak Takaine
#
#// This FIJI macro allows automatic detection and analysis of mitochondria labeled by Cox4-mNG.
#// As input, two channel image file that contains a pair of fluorescence (Ch #1) and phase contrast (Ch #2) microscopic images is required.
#// This macro is optimized for images accuired by using a 100x objective lens.

#@ String (label="Comments, e.g., 2024-09-21") edate1
#@ File (label="Choose source Folder", style="directory") dirS0 
#@ File (label="Choose destination Folder", style="directory") dirD0
#@ Double (min=1, max=50, value=10, stepSize = 0.5, style="format:##.#",label="Radius for BackgroundSubtraction") rb_radius
#@ Double (min=0, max=10, value=0.2, stepSize = 0.1, style="format:#.#",label="Minimum size for particle analysis (in µm^2)") min_size

from __future__ import division
from ij import IJ, ImagePlus, Prefs
from ij.process import ImageStatistics as IS
options = IS.MEAN | IS.AREA | IS.STD_DEV  # many others
from ij.plugin import Duplicator
from ij.plugin import ZProjector as ZP
from ij.gui import Roi
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin import ChannelSplitter as CS
from ij.plugin import RGBStackMerge, RGBStackConverter
from ij.plugin.filter import BackgroundSubtracter as BS
createBackground = False
lightBackground = False
useParaboloid = False
doPresmooth = False
correctCorners = False
#from ij.plugin.filter import ParticleAnalyzer as PA
#options_pa = PA.SHOW_NONE + PA.CLEAR_WORKSHEET+PA.ADD_TO_MANAGER# PA.EXCLUDE_EDGE_PARTICLES + PA.SHOW_RESULTS+PA.SHOW_OUTLINES

from java.util import Date
from java.text import SimpleDateFormat

import os
from os import path
from org.apache.commons.math3.stat.descriptive import DescriptiveStatistics as DSS

timeStamp = SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(Date())
print timeStamp

# Save the Result Table in csv format
def save_result_table(directory, filename, result_table):
    resultfile = os.path.join(directory, filename + ".csv") 
    result_table.saveAs(resultfile)

# Save the image file in tiff format
def save_image_as_tif(directory, filename, image):
    outputfile = os.path.join(directory, filename + ".tif")
    IJ.saveAs(image, "TIFF", outputfile) # 保存する画像、tiff形式、パスを指定
    
# Analyze the image of fluorescent condensates
def analyse_mito_image(current_file_path):
	imp0 = IJ.openImage(current_file_path)
	filename = imp0.getTitle().split(".")[0]
	
	# Max intensity projection
	imp0_mip = ZP.run(imp0, "max all")
	
	cs = CS()
	image_list = cs.split(imp0_mip)
	ch1 = image_list[0] # Channel #1：fluorescence image
	ch2 = image_list[1] # Channel #2：bright field/phase-contrast image
	
	ch1_mask = ch1.duplicate()
	ch1_mask.setTitle("mask")
	
	# BackGround subtraction,  Important step to detect fine fluorescent condensates
	bs = BS()
	bs.rollingBallBackground(ch1_mask.getProcessor(), rb_radius, createBackground, lightBackground, useParaboloid, doPresmooth, correctCorners)
	
	ch1_bgsub = ch1_mask.duplicate()
	ch1_bgsub.setTitle("BG_subtracted")
	
	# Convert to binary image
	IJ.setAutoThreshold(ch1_mask, "Otsu dark")
	Prefs.blackBackground = True
	IJ.run(ch1_mask, "Convert to Mask", "")
	
	# Detection of mitochondria
	IJ.run("Set Measurements...", "area mean perimeter bounding fit shape feret's integrated redirect=None decimal=3")
	IJ.run(ch1_mask, "Analyze Particles...", "size="+str(min_size)+"-500 exclude clear add")  # size unit is µm, min_sizeはstr（）で文字列に変換してから連結
	rt = ResultsTable.getResultsTable()
	curr_pnums = rt.size()
	
	#rt.show("Results")
	rm = RoiManager.getRoiManager()
	if curr_pnums != 0: # if particles were detected
		rt.reset()
		rm.runCommand(ch1,"Show None")
		rm.runCommand(ch1,"Show All")
		rm.runCommand(ch1,"Measure")
		for j in range(0, rt.size()):
			rt.setValue("date", j, edate)
			rt.setValue("file", j, filename)
			rt.setValue("rb_radius", j, rb_radius)
			rt.setValue("min_size", j, min_size)
		rt.show("Results")
	else: # if no particels were detected
		rt.reset()
		rt.setValue("Area", 0, "NaN")
		rt.setValue("Mean", 0, "NaN")
		rt.setValue("Min", 0, "NaN")
		rt.setValue("Max", 0, "NaN")
		rt.setValue("X", 0, "NaN")
		rt.setValue("Y", 0, "NaN")
		rt.setValue("Circ.", 0, "NaN")
		rt.setValue("IntDen", 0, "NaN")
		rt.setValue("RawIntDen", 0, "NaN")
		rt.setValue("AR", 0, "NaN")
		rt.setValue("Round", 0, "NaN")
		rt.setValue("Solidity", 0, "NaN")
		rt.setValue("Perim.", 0, "NaN")
		rt.setValue("Feret", 0, "NaN")
		rt.setValue("date", 0, edate)
		rt.setValue("file", j, filename)
		rt.setValue("rb_radius", 0, rb_radius)
		rt.setValue("min_size", 0, min_size)
		rt.show("Results")
	
	# Extract each variable as a list of values from the table to calculate descriptive statistics
	area = rt.getColumnAsDoubles(rt.getColumnIndex("Area"))
	meanint = rt.getColumnAsDoubles(rt.getColumnIndex("Mean"))
	intden = rt.getColumnAsDoubles(rt.getColumnIndex("IntDen"))
	rounds = rt.getColumnAsDoubles(rt.getColumnIndex("Round"))
	ar = rt.getColumnAsDoubles(rt.getColumnIndex("AR"))
	circ = rt.getColumnAsDoubles(rt.getColumnIndex("Circ."))
	solidity = rt.getColumnAsDoubles(rt.getColumnIndex("Solidity"))
	perimeter = rt.getColumnAsDoubles(rt.getColumnIndex("Perim."))
	feret = rt.getColumnAsDoubles(rt.getColumnIndex("Feret"))
	
	# Collect into a dictionary
	params = {"area":area, "meanint":meanint, "intden":intden, "rounds":rounds, "ar":ar, "circ":circ, "solidity":solidity, "perimeter":perimeter, "feret":feret}

	# Generate DSS instance, input values of the variable, collect the instances into a dictionary
	dstats = {}
	for k,v in params.items():
		dss = DSS()
		for j in range(0, len(v)):
			dss.addValue(v[j])
		dstats[k] = dss
	
	# ch1_bgsub: BG-subtracted fluorescence image, ch2: phase contrast image, ch1_mask: mask image, rt: ResultsTable
	return ch1_bgsub,ch2,ch1_mask,rt,dstats,filename 
    # These return values are gathered into a tapple. 

#### Main code
# Insert a blank to prevent automatic modification on Excel.
edate = " "+edate1
# Make directories
dirD = os.path.join(str(dirD0), edate1 + "_output")
if not os.path.exists(dirD):
	os.mkdir(dirD)
dirBF = os.path.join(str(dirD), "BF")
if not os.path.exists(dirBF):
	os.mkdir(dirBF)                           
#Create a folder for mask images and ROI data
dirDR = os.path.join(str(dirD), "segmented")
if not os.path.exists(dirDR):
	os.mkdir(dirDR)
dirGreen = os.path.join(str(dirD), "green")
if not os.path.exists(dirGreen):
	os.mkdir(dirGreen)
dirCSV = os.path.join(str(dirD), edate1 + "_csv")
if not os.path.exists(dirCSV):
	os.mkdir(dirCSV)

# Acquire a list of files in the directory
filelist = os.listdir(str(dirS0))

# List comprehension, extract nd2 files.
nd2_files = [f for f in filelist if f.split(".")[-1] == "nd2"]
#filenames = [f.split(".")[0] for f in filelist]
nd2_files = sorted(nd2_files)

# Create a table that summarises the averages of particle parameters in a image file
particle_stat = ResultsTable()

for nd2_file in reversed(nd2_files):  # reversed() generates a revered iterator
    current_file_path = os.path.join(str(dirS0), nd2_file) 
    results = analyse_mito_image(str(current_file_path))
    ch1_bgsub = results[0] 
    ch2 = results[1]
    ch1_mask = results[2]   
    rt = results[3]    
    dstats = results[4]
    filename = results[5]
    
    save_result_table(str(dirCSV), filename, rt)
    save_image_as_tif(str(dirGreen), filename, ch1_bgsub)
    save_image_as_tif(str(dirBF), filename, ch2)
    save_image_as_tif(str(dirDR), filename, ch1_mask)
    
    particle_stat.addValue("date", edate)
    particle_stat.addValue("file_name", filename)
    particle_stat.addValue("rb_radi", rb_radius)
    particle_stat.addValue("min_sizes", min_size)
    particle_stat.addValue("particle_nums", dstats["area"].getN())
    particle_stat.addValue("mean_areas", dstats["area"].getMean())
    particle_stat.addValue("mean_perim.", dstats["perimeter"].getMean())
    particle_stat.addValue("mean_meanints", dstats["meanint"].getMean())
    particle_stat.addValue("mean_intdens", dstats["intden"].getMean())
    particle_stat.addValue("mean_feret", dstats["feret"].getMean())
    particle_stat.addValue("mean_roundness", dstats["rounds"].getMean())
    particle_stat.addValue("mean_AR", dstats["ar"].getMean())
    particle_stat.addValue("mean_circ", dstats["circ"].getMean())
    particle_stat.addValue("mean_solidity", dstats["solidity"].getMean())
    particle_stat.addRow()

#　Remove the last empty row
particle_stat.deleteRow(particle_stat.size() - 1)
save_result_table(str(dirD), edate1+"_particle_stat", particle_stat)

print "Done. \n"
IJ.run("Clear Results")
rm = RoiManager.getRoiManager()
rm.reset()
IJ.run("Close All")

timeStamp = SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(Date())
print timeStamp