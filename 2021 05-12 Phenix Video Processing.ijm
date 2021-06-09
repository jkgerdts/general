root = "/Volumes/New Volume/2021 05 immune synapse pilot__2021-05-10T00_18_13-Measurement 1/Images/";

r=7;
c=9;
f=1;
red_ch = 1;
green_ch = 2;
p_len = 5;
sk_len=30;


for (f=8; f<9; f++){
	openmovie(r,c,f,p_len,green_ch,"green",sk_len);
	openmovie(r,c,f,p_len,red_ch,"red",sk_len);
	run("Merge Channels...", "c1=red_movie c2=green_movie create");
	rename(f);
}

function filename(r,c,f,p,ch,sk){
	return "r"+zeropad(r)+"c"+zeropad(c)+"f"+zeropad(f)+"p"+zeropad(p)+"-ch"+ch+"sk"+sk+"fk1fl1.tiff";	
}

function zeropad(integer){
	st = d2s(integer,0);
	if (st.length==1){
		return "0"+st;
	}
	return st;
}

function openimage(r,c,f,p_len,ch,color,sk){
	for (p=1; p<p_len+1; p++){
		open(root+filename(r,c,f,p,ch,sk));
		rename(color+"_slice_"+p);
	}
	run("Images to Stack", "name=Stack title=slice use");
	rename(color);
	run("Subtract Background...", "rolling=100 stack");
	run("Gamma...", "value=0.6 stack");
	run("Z Project...", "projection=[Max Intensity]");
	rename(color+"_project");
	selectWindow(color);
	close();
	selectWindow(color+"_project");
}

function openmovie(r,c,f,p_len,ch,color,sk_len){
	for (sk=1; sk<sk_len+1; sk++){
		openimage(r,c,f,p_len,ch,color,sk);
	}
	run("Images to Stack", "name=Stack title=project use");
	rename(color+"_movie");
}


//run("Concatenate...", "all_open title=sequence open");
//run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1 frame=[60 sec]");
//run("Split Channels");

//stacklabels = newArray("C1-sequence","C2-sequence","C3-sequence");

//for (i=0; i<3; i++){
//	selectWindow(stacklabels[i]);
//	run("Subtract Background...", "rolling=25 stack");
//	run("Gamma...", "value=0.6 stack");
//	run("Gaussian Blur...", "sigma=1 stack");
//	normalize();
//}

function normalize(){
	getHistogram(values,counts,65536);
	cumsum = 0;
	thr_lo = 512*512*0.8;
	thr_hi = 512*512*0.995;
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
