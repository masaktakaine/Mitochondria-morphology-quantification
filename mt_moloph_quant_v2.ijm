// Mitochodria molphology quantification
// author: Masak Takaine

// This FIJI macro allows automatic detection and analysis of mitochondria labeled by Cox4-mNG.
// As input, two channel image file that contains a pair of fluorescence (Ch #1) and phase contrast (Ch #2) microscopic images is required.
// This macro is optimized for images accuired by using a 100x objective lens.

macro "Mt_molphology_quantification" {

#@ String (label="Comments, e.g., 2024-09-21") edate1
#@ File (label="Choose source Folder", style="directory") dirS1 
#@ File (label="Choose destination Folder", style="directory") dirD1
#@ Double (min=1, max=50, value=10, stepSize = 0.5, style="format:##.#",label="Radius for BackgroundSubtraction") rb_radius
#@ Double (min=0, max=10, value=0.2, stepSize = 0.1, style="format:#.#",label="Minimum size for particle analysis (in µm^2)") min_size
#@ String (label="Hide/Show the active image? The Show slows the analysis.", choices={"hide","show"}, style="radioButtonHorizontal") sbm

setBatchMode(sbm); // hides the active image, required ImageJ 1.48h or later
dirS = dirS1 + File.separator; // "File.separator" returns the file name separator character depending on the OS system used.
dirD1 = dirD1 + File.separator;

//setOption("ExpandableArrays", true);  // Enables/disables support for auto-expanding arrays, In ImageJ 1.53g or newer, arrays automatically expand in size as needed.
edate = " "+edate1; // Insert a blank to prevent automatic modification on Excel.

imagefilelist = getFileList(dirS);
date =newArray;
file_name = newArray; 	// An array to store filenames
rb_radi = newArray;
min_sizes = newArray;
particle_nums = newArray;

mean_areas = newArray();
mean_meanints = newArray();
mean_intdens = newArray();
sum_intdens = newArray();
mean_roundness = newArray();
mean_AR = newArray();
mean_circ = newArray();
mean_solidity = newArray();
mean_perim = newArray();
mean_feret = newArray();

dirD = dirD1 +"/"+ edate1 + "_output";
File.makeDirectory(dirD);

dirBF = dirD + "/phase/";		//Create a folder for phase-contrast or biright-field images
File.makeDirectory(dirBF);				

dirGreen = dirD + "/green/";  //Create a folder for fluorescent images
File.makeDirectory(dirGreen);				

dirDR = dirD + "/segmented/";  //Create a folder for mask images and ROI data
File.makeDirectory(dirDR);

dirCSV = dirD +"/" + edate1 + "_csv/";
File.makeDirectory(dirCSV);	  // Create a folder for csv files

for (i = 0; i < imagefilelist.length; i++) {
    currFile = dirS + imagefilelist[i];
    if((endsWith(currFile, ".nd2"))||(endsWith(currFile, ".oib"))||(endsWith(currFile, ".zvi"))) { // process if files ending with .oib or .nd2, or .zvi
		run("Bio-Formats Macro Extensions"); 
		Ext.openImagePlus(currFile)}
	else if ((endsWith(currFile, ".tif"))||(endsWith(currFile, ".tiff"))) {// process if files ending with .tif or .tiff (hyperstack files)
			open(currFile); 
		}
run("Clear Results");						// Reset Results window
print("\\Clear"); 
title = getTitle();
// Remove the extension from the filename
title_s = replace(title, "\\.nd2", ""); 
title_s = replace(title_s, "\\.tif", "");
title_s = replace(title_s, "\\.tiff", "");

run("Z Project...", "projection=[Max Intensity]");
title2 = getTitle();
run("Split Channels");
c1 = "C1-" + title2; // Channel #1：fluorescence image
c2 = "C2-" + title2; // Channel #2：bright field/phase-contrast image

selectWindow(c2);
BFID = getImageID();
selectImage(BFID);
saveAs("Tiff", dirBF + title_s);
close();

selectWindow(c1);						
GRID = getImageID();
run("Grays");
run("Duplicate...", "title=Temp"); //Duplicate and rename
selectWindow("Temp");
run("Subtract Background...", "rolling=rb_radius disable");
run("Duplicate...", "title=BG_subtracted"); //Duplicate and rename
wait(200); // 200 ms待つ、待たないと描画させないで画像処理した場合に処理がおかしくなる

// The intensity threshold is determined by the method of Otsu.
selectWindow("Temp");
setAutoThreshold("Otsu dark");
setOption("BlackBackground", true);
run("Convert to Mask");

// Detection of mitochodria
run("Set Measurements...", "area mean perimeter bounding fit shape feret's integrated redirect=None decimal=3"); 
run("Analyze Particles...", "size=min_size-500 circularity=0-1.00 show=Masks exclude clear add");
particle_nums[i] = nResults;

if (nResults !=0) {  // if particles were detected
run("Clear Results");
//selectImage(GRID);; // Superimpose ROIs on the original fluorescence image
selectWindow("BG_subtracted");

roiManager("Show None");
roiManager("Show All");
roiManager("Measure");
for(k=0; k<nResults; k++) {
 setResult("date",k,edate);
 setResult("file",k,title_s);
 setResult("rb_radius", k, rb_radius);
 setResult("min_size", k, min_size);
}

// Because Table.getColumn() does not work in batch-mode ("hide"),
//the specified column in the Results table is obtained as an array by iterating getResult().
colArea = newArray();
colMean = newArray();
colIntDen = newArray();
colRound = newArray();
colAR = newArray();
colCirc = newArray();
colSolid = newArray();
colPerim = newArray();
colFeret = newArray();

for (p=0; p<nResults; p++){ 
	colArea[p] =getResult("Area", p);
	colMean[p] =getResult("Mean",p);
	colIntDen[p] = getResult("IntDen",p);
	colRound[p] = getResult("Round",p);
	colAR[p] = getResult("AR",p);
	colCirc[p] = getResult("Circ.",p);
	colSolid[p] = getResult("Solidity",p);
	colPerim[p] = getResult("Perim.",p);
	colFeret[p] = getResult("Feret", p);
}
Array.getStatistics(colArea, min1, max1, mean1, stdDev1);
// Return min, max, mean and stdDev of the Area column to min1, max1, mean1 and stdDev1, respectively
Array.getStatistics(colMean, min2, max2, mean2, stdDev2);
Array.getStatistics(colIntDen, min3, max3, mean3, stdDev3);
Array.getStatistics(colRound, min4, max4, mean4, stdDev4);
Array.getStatistics(colAR, min5, max5, mean5, stdDev5);
Array.getStatistics(colCirc, min6, max6, mean6, stdDev6);
Array.getStatistics(colSolid, min7, max7, mean7, stdDev7);
Array.getStatistics(colPerim, min8, max8, mean8, stdDev8);
Array.getStatistics(colFeret, min9, max9, mean9, stdDev9);
saveAs("Results", dirCSV + title_s + ".csv");
} else{ // if no particels were detected
	run("Clear Results");
	setResult("Area", 0, "NaN");
	setResult("Mean", 0, "NaN");
	setResult("Min", 0, "NaN");
	setResult("Max", 0, "NaN");
	setResult("X", 0, "NaN");
	setResult("Y", 0, "NaN");
	setResult("Circ.", 0, "NaN");
	setResult("IntDen", 0, "NaN");
	setResult("RawIntDen", 0, "NaN");
	setResult("AR", 0, "NaN");
	setResult("Round", 0, "NaN");
	setResult("Solidity", 0, "NaN");
	setResult("Perim.", 0, "NaN");
	setResult("Feret", 0, "NaN");
	setResult("date",0,edate);
 	setResult("file",0,title_s);
 	setResult("rb_radius", 0, rb_radius);
	setResult("min_size", 0, min_size);
 saveAs("Results", dirCSV + title_s + ".csv");
 wait(100);
}

selectWindow("Mask of Temp");
roiManager("Show None");
roiManager("Show All"); 	// Show all ROIs to save the ROIs as overlays
saveAs("Tiff", dirDR +title_s);
close();

//selectImage(GRID);; // Superimpose ROIs on the original fluorescence image
selectWindow("BG_subtracted");
roiManager("Show None");
roiManager("Show All");
saveAs("Tiff", dirGreen + title_s);
close();

run("Close");	
run("Close All");
run("Clear Results");
print("\\Clear");
roiManager("reset");

date[i] = edate;
file_name[i] = title_s;
rb_radi[i] = rb_radius;
min_sizes[i] = min_size;

if (particle_nums[i] !=0) {
mean_areas[i] = mean1;
mean_meanints[i] = mean2;
mean_intdens[i] = mean3;
sum_intdens[i] = particle_nums[i]*mean_intdens[i];
mean_roundness[i] = mean4;
mean_AR[i] = mean5;
mean_circ[i] = mean6;
mean_solidity[i] = mean7;
mean_perim[i] = mean8;
mean_feret[i] = mean9;
} else{
	mean_areas[i] = NaN;
	mean_meanints[i] = NaN;
	mean_intdens[i] = NaN;
	sum_intdens[i] = NaN;
	mean_roundness[i] = NaN;
	mean_AR[i] = NaN;
	mean_circ[i] = NaN;
	mean_solidity[i] = NaN;
	mean_perim[i] = NaN;
	mean_feret[i] = NaN;
}
}
 
// Show statistics of the particles in a new table.	
Array.show("particle_stat(row numbers)", date, file_name, rb_radi, min_sizes, particle_nums, mean_areas, mean_perim,
mean_meanints, mean_intdens, sum_intdens, mean_feret, mean_roundness, mean_AR, mean_circ, mean_solidity);
	selectWindow("particle_stat");
    saveAs("Results", dirD +"/"+ edate1 + "_particle_stat.csv");
	run("Close");	
    run("Clear Results");						// Reset Results window
	print("\\Clear"); 							// Reset Log window
	roiManager("reset");						// Reset ROI manager
	run("Close All");
	run("Close All");
	showMessage(" ", "<html>"+"<font size=+2>Process completed<br>");
}