root = "/Volumes/HooToo/2021 05-17 MDCK OFF ICAM/";

function filename(ln,time,pos){
	return "LN"+LN+"_"+time+"/LN"+ln+"_"+time+"_MMStack_Pos"+pos+".ome.tif";
}



// Here's the plan: I have quantified 2 fields from lane 2; will concatenate those and treat as one data point; 
// next need to quantify lane 1 from this same dataset and use the two measurements to make the ligand-free part of the figure

LN=4;
pos=2;
nframes = 9;
imnames = newArray(nframes);
concatstring = " title=Concat open ";

ch_blue = "C1";
ch_green = "C2";
ch_red = "C3";
ch_farred = "C4";
ch_gray = "C5";

for (i=1; i<nframes+1; i++){
	open(root+filename(LN,i,pos));
	imnames[i-1]=getTitle();
	concatstring = concatstring+"image"+toString(i)+"="+imnames[i-1]+" ";
}

run("Concatenate...", concatstring);

//run("Concatenate...", "all_open title=sequence open");
run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1 frame=[60 sec]");
run("Split Channels");

colors = newArray(ch_blue,ch_red,ch_farred);

selectWindow(ch_gray+"-Concat");
run("Translate...", "x=3 y=0 interpolation=None stack");

for (i=0; i<3; i++){
	selectWindow(colors[i]+"-Concat");
	run("Subtract Background...", "rolling=25 stack");
	run("Gamma...", "value=0.6 stack");
	run("Gaussian Blur...", "sigma=1 stack");
	normalize();
}

run("Merge Channels...", "c1="+ch_farred+"-Concat c2="+ch_red+"-Concat c3="+ch_blue+"-Concat c4="+ch_gray+"-Concat create keep ignore");
run("Stack to RGB", "frames");

function normalize(){
	getHistogram(values,counts,65536);
	cumsum = 0;
	thr_lo = 512*512*0.8;
	thr_hi = 512*512*0.99;
	thr_lo_val = 0;
	thr_hi_val = 0;
	unsurpassed = true;
	for (i=0; i<65536; i++){
		cumsum = cumsum+counts[i];
		if (cumsum>=thr_lo){
			if (thr_lo_val==0){
				thr_lo_val=i+1;
			}
		}
		if (cumsum>=thr_hi){
			if (thr_hi_val==0){
				thr_hi_val=i+1;
			}
		}
	}
	setMinAndMax(thr_lo_val,thr_hi_val);
	run("Apply LUT", "stack");
}
