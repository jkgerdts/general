
cycle = 24;

run("Close All");

//open("/Volumes/HooToo/2021 04-22 HuG/LN5_crawling_with_fluor_q20/LN5_crawling_with_fluor_q20_MMStack_Pos0.ome.tif");

open("/Volumes/HooToo/2021 04-22 HuG/LN5_crawling_with_fluor_q20/LN5_crawling_with_fluor_q20_MMStack_Pos1.ome.tif");

im = getTitle();
channels = newArray("blue","red","farred");

for (i=0; i<3; i++){
	selectWindow(im);
	len = nSlices;
	run("Make Substack...", "delete slices=1-"+len+"-"+cycle-i);
	rename(channels[i]);
}

for (i=0; i<3; i++){
	selectWindow(channels[i]);
	run("Subtract Background...", "rolling=25 stack");
	run("Gamma...", "value=0.6 stack");
	run("Gaussian Blur...", "sigma=1 stack");
	normalize();
	run("Make Substack...", "  slices=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4");
	rename(channels[i]+"_expanded");
	selectWindow(channels[i]);
	//close();
}


selectWindow(im);
setSlice(nSlices);
run("Add Slice");

run("Merge Channels...", "c1=farred_expanded c2=red_expanded c3=blue_expanded c4="+im+" create keep ignore");

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
