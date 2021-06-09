x=newArray(46,58,58,82,74,104,120,120,133,80,40,40,11,6,221,240,252,230,361,216,119,120,77,82,80,76,72,70,180,183,176,177,160,158,161,149,149,149,148,148,132,93,88,115,125,168,224,221,223,288,331,330,424,420,402,403,417,441,455,464,451,456,456,457,467,465,483,374,382,328,378,398,382,375,397,396,398,398,245,151,105,59,288,330,354,350,348,486,458,464,462,459,459,454,453,321,262,44,46,43,42,28,24,94,96,106,95,88,82,220,219,221,221);
y=newArray(110,30,36,25,30,101,120,128,164,184,250,241,360,362,269,362,393,423,82,115,121,127,246,255,263,262,258,258,311,320,324,324,303,303,306,255,257,259,262,266,230,272,286,288,338,232,177,189,194,139,152,160,212,228,331,336,317,321,308,320,324,346,345,340,376,381,340,407,450,451,451,143,72,100,18,19,17,26,59,28,44,96,188,222,226,230,234,272,216,214,92,102,115,142,149,166,129,407,454,460,464,481,491,438,443,391,372,418,410,492,495,498,504);
t=newArray(0,0,2,0,4,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0,0,1,0,1,2,0,1,2,0,1,0,1,0,1,2,0,1,2,3,4,0,0,1,0,1,0,0,0,1,0,0,2,0,3,0,1,0,0,0,0,0,0,1,2,0,1,0,0,0,0,0,0,0,1,0,1,2,3,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,2,3);

L = x.length;

red=measure("C3-sequence",x,y,t);
green=measure("C2-sequence",x,y,t);
blue=measure("C1-sequence",x,y,t);

print("red,green,blue");
for (i=0; i<L; i++){
	print(red[i]+","+green[i]+","+blue[i]);
}

function measure(imname,xarr,yarr,tarr){
	selectWindow(imname);
	arr = newArray(xarr.length);
	for (i=0; i<t.length; i++){
		setSlice(t[i]+1);
		makeRectangle(x[i]-1,y[i]-1,3,3);
		run("Measure");
		arr[i]=getResult("Mean");
	}
	return arr;
}