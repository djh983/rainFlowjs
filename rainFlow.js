////////////////////////////////////////
/////Rain Flow JS Library
/////2022/02/17
/////Daniel Hill
////////////////////////////////////////


////////////////////////////////////////
//Function 1: "getPks"
////Input 1: "dset" = dataset 1D array
////Input 2: "tm" = time data array for Input 1
////Output 1: "pkSet" = 1D array of each peak and valley, defined as a change in trajectory
////Output 2: "tmSet" = time data for Output 1

function getPks(dset, tm){
var i = 1;
var traj = 0;
var pkSet = [];
var tmSet = [];
dset.forEach(function(cv){
	if(cv < dset[i]){
		if(traj != -1){
			pkSet.push(cv);
			tmSet.push(tm[i-1]);				
		}
		traj = -1;
	}
	else if(cv > dset[i]){
		if(traj != 1){
			pkSet.push(cv);
			tmSet.push(tm[i-1]);
		}
		traj = 1;
	}else{
		traj = 0;
	}
	i++;
})

return[pkSet, tmSet];
}
/////////////////////////////////////////


/////////////////////////////////////////
//Function 2: "hystGate"
///Input 1: "dset0" = dataset 1D array
///Input 2: "gt" = hysteresis gate size 
///Input 3: "tm0" = time data array for Input 1 
///Output 1: "filtSet" = filtered data where peak/valleys less than the gate are removed
///Output 2: "tmSet" = time data array for Output 1

function hystGate(dset0, gt, tm0){
var dsetarr = getPks(dset0, tm0);
var i = 0;
var j = 0;
var cvtmp = 0;
var filtSet = [];
var tmSet = [];
var dset = dsetarr[0];
var tm = dsetarr[1];
var tmp = 0;
var traj = 0;
dset.forEach(function(cv){
if(j == 0){
filtSet.push(cv);
tmSet.push(tm[0]);
cvtmp = cv;
j++;
}
if(Math.abs(cv - cvtmp) > gt){
	if(cv > cvtmp){ //upward trajectory
		if(traj == 1){
			filtSet.pop();
			tmSet.pop();
			j--;
		}
	traj = 1;
	}else{
		if(traj == -1){
			filtSet.pop();
			tmSet.pop();
			j--;
		}
	traj = -1; //downward trajectory
	}
	filtSet.push(cv);
	tmSet.push(tm[i]);
	if(j == 1){ //special case for 1st point
		for(var k=0; k<i;k++){
			if((filtSet[0] - dset[k])*traj > 0){
				filtSet[0] = dset[k];
				tmSet[0] = tm[k];
			}
		}
	}
	j++;
	
}else if((cv - cvtmp)*traj > 0){
filtSet.pop();
tmSet.pop();
filtSet.push(cv); //replace previous
tmSet.push(tm[i]);
}

cvtmp = filtSet[j-1];
i++;
})

return [filtSet, tmSet];
}
////////////////////////////////////////////////

////////////////////////////////////////////////
//Function 3: "rainFlowBasic" = basic discrete rainflow count using 4 point method, no average correction
///Input 1: "dset" = dataset array of peaks and valleys; recommended to be filtered through hystGate 1st
///Input 2: "numBin" = desired number of bins for the output
///Input 3: "gt" = the gate used in the hysteresis filter, sets the minimum for the bins
///Output 1: "bins" = array of bin values
///Output 2: "cnts" = array of counts for each bin

function rainFlowBasic(dset,numBin, gt){ //basic discrete rainflow count using 4 point method, no average
gt*=1;
var dsetr = [];
var amp = [];
dset.forEach(function(cv){
	dsetr.push(cv);
});

var bins = [];
var cnts = [];

var a = 0;

while(a == 0){
for(var i=0; i<dsetr.length; i++){
	if((i+4) >= dsetr.length){
		a = 1;
		break;
	}
	var p1 = dsetr[i];
	var p2 = dsetr[i+1];
	var p3 = dsetr[i+2];
	var p4 = dsetr[i+3];
	if(Math.abs(p1-p4) > Math.abs(p2-p3)){ //cycle is bound
		var rv = Math.abs(p2-p3);
		amp.push(rv);
		dsetr.splice(i+1, 2);	
		break;
	}	
} ///for
} ///while

var mxamp = amp.reduce(function(a, b){
	return Math.max(a,b);
}, -Infinity);

var stp = Math.round((mxamp - gt)/numBin);
for(var i=0; i<numBin; i++){
	bins.push(i*stp + gt);
	cnts.push(0);
}

amp.forEach(function(cv){
	for(var k = bins.length -1; k >= 0; k--){
		if(cv >= bins[k]){
			cnts[k]++;
				break;
		}
	}

});
return [bins, cnts];
}
///////////////////////////////////////////////////////////////////

//////////////////////////////////
//Function 4: "fromTo"  from-to rainflow count using 4 point method, includes mean effects
///Input 1: "dset" = peak-valley dataset, suggested to be filtered through hysteresis gate
///Input 2: "numBin" = desired number of bins for the output
///Output 1: "bins" = array of bin values
///Output 2: "cnts" = bin# x bin# array of counts for each bin
///Output 3: "xvars" = array of 'to' values
///Output 4: "yvars" = array of 'fom' values

function fromTo(dset, numBin){
var dsetr = [];
var amp = [];
dset.forEach(function(cv){
	dsetr.push(cv);
});

var mx = dsetr.reduce(function(a, b){
	return Math.max(a,b);
}, -Infinity);
var mn = dsetr.reduce(function(a,b){
	return Math.min(a,b);
}, Infinity);

var stp = Math.round((mx - mn)/numBin);

var bins = [];
var cnts = [[]]; //2D array
var xvars = [];
var yvars = [];

for (var i=0; i<numBin; i++){
	bins.push(mn + i*stp);
	cnts[i] = [];
	for(var j=0; j<numBin; j++){
		cnts[i][j] = 0; //all matrix elements are initialized to zero
	}
}
var a = 0;
while(a == 0){
for(var i=0; i<dsetr.length; i++){
	if((i+4) >= dsetr.length){
		a = 1;
		break;
	}
	var p1 = dsetr[i];
	var p2 = dsetr[i+1];
	var p3 = dsetr[i+2];
	var p4 = dsetr[i+3];
	if(Math.abs(p1-p4) > Math.abs(p2-p3)){
		for(var j=bins.length - 1; j>=0; j--){
			if(p2 >= bins[j]){
				break;
			}
		}
		for(var k=bins.length - 1; k>=0; k--){
			if(p3 >= bins[k]){
				break;
			}
		}
		cnts[k][j]++;
		//xvars.push(bins[k]);
		//yvars.push(bins[j]);
		xvars.push(p2);
		yvars.push(p3);

		dsetr.splice(i+1, 2);
		break;
	}
} ///for
} ///while
return [bins, cnts, xvars, yvars];
}
/////////////////////////////////////////////////////////////


/////////////////////////////////////
//Function 5: "ampAvg"  2D rainflow using 4 point method
///Input 1: "dset" = filtered peak-valley dataset
///Input 2: "numBin" = desired number of bins
///Input 3: "gt" = gate size used in filtering
///Output 1: "amp" = array of amplitudes
///Output 2: "avgs" = array of averages
///Output 3: "ampbins" = array of amplitude bins
///Output 4: "avgbins" = array of average bins
///Output 5: "cnts" = ampbins# x avgbins# array of counts for each bin
///Output 6: "damp" = array of discrete amplitudes (bins) to simplify 3dplots
///Output 7: "davg" = array of discrete averages (bins) to simplify 3dplots

function ampAvg(dset, numBin, gt){ //Rainflow with averaging
var dsetr = [];
var amp = [];
var avgs = [];
var ampbins = [];
var avgbins = [];
var cnts = [[]];
var damp = [];
var davg = [];
dset.forEach(function(cv){
	dsetr.push(cv);
});


var a = 0;

while(a == 0){
for(var i=0; i<dsetr.length; i++){
	if((i+4) >= dsetr.length){
		a = 1;
		break;
	}
	var p1 = dsetr[i];
	var p2 = dsetr[i+1];
	var p3 = dsetr[i+2];
	var p4 = dsetr[i+3];
	if(Math.abs(p1-p4) > Math.abs(p2-p3)){ //cycle is bound
		var rv = Math.abs(p2-p3);
		var av = (p2 + p3)/2;
		amp.push(rv);
		avgs.push(av);
		dsetr.splice(i+1, 2);	
		break;
	}	
} ///for
} ///while

var mx = avgs.reduce(function(a, b){
	return Math.max(a,b);
}, -Infinity);
var mn = avgs.reduce(function(a,b){
	return Math.min(a,b);
}, Infinity);
var mxamp = amp.reduce(function(a, b){
	return Math.max(a,b);
}, -Infinity);

var ampstp = Math.round(Math.abs((mxamp)/numBin));
var avgstp = Math.round(Math.abs((mx-mn)/numBin)); 
for (var i=0; i<numBin; i++){
	avgbins.push(mn + i*avgstp);
	ampbins.push(i*ampstp);
	cnts[i] = [];
	for(var j=0; j<numBin; j++){
		cnts[i][j] = 0; //all matrix elements are initialized to zero
		
	}
}
var thresh = 0.4;
for(var i=0; i<amp.length;i++){
	for(var j=ampbins.length -1; j>=0; j--){
		if(amp[i]>=ampbins[j]){
			var tmpamp = amp[i];
			break;
		}
	}
	for(var k=avgbins.length -1; k>=0; k--){
		if(avgs[i]>=avgbins[k]){
			var tmpavg = avgs[i];
			break;
		}
	}
	cnts[j][k]++;
	if((j>numBin*thresh) && ((k<numBin*thresh/2) || (k>numBin*thresh/2))){
		damp.push(tmpamp);
		davg.push(tmpavg);
	}
}
return [amp, avgs, ampbins, avgbins, cnts, damp, davg];
}
//////////////////////////////////////////////////////////


//////////////////////////////////////////////
//Function 6:"meanCorrectedRf"  uses the 4 point method for rainflow and Walker method (SWT in general form) to correct for mean stresses; should only be used with strain data
///Input 1: "ampset" = array of amplitudes from ampAvg function
///Input 2: "meanset" = array of means from ampAvg function
///Input 3: "numBin" = desired number of bins for the rainflow output
///Output 1: "rfcbins" = mean corrected rainflow bins
///Output 2: "rfccnts" = mean corrected rainflow counts

function meanCorrectedRf(ampset, meanset, numBin){ //Rainflow with mean strain correction, using Walker Method (SWT in general form);
var dsetr = [];
var gam = 0.651;  //walker coefficient, estimated for steels with Tensile strength 300-350MPa
var rfcbins = [];
var rfccnts = [];
var coramp = [];

var i = 0;
ampset.forEach(function(cv){
	var mn = meanset[i];
	var sr = (mn - (cv/2))/(mn + (cv/2));  //stress ratio
	
	if(sr < 0.98){     // be sure the mean correction is not excessive
	var tmp = (2/(1-sr))^gam; 	
	}else{
		tmp = 1;
	}
	coramp.push(tmp*cv);
	i++;
});	

var mx = coramp.reduce(function(a, b){
	return Math.max(a,b);
}, -Infinity);

var stp = Math.round(mx/numBin);

for(var j = 0; j<numBin; j++){
	rfcbins.push(j*stp);
	rfccnts.push(0);
}

coramp.forEach(function(cv){
for(var k = numBin-1; k >=0; k--){
	if(cv > rfcbins[k]){
		rfccnts[k]++;
		break;
	}	
}
});
return [rfcbins, rfccnts];
}
////////////////////////////////////////////////////////////

////////////////////////////////////////
//Function 7: "makeRFTable"  creates an HTML table of rainflow data
///Input 1: "bn" = array of bins
///Input 2: "cnt" = array of counts for each bin
///No outputs: appends html page with a table

function makeRFTable(bn, cnt){
	var i = 0;
	const newTab = document.createElement("table");
	var th1 = document.createElement('tr');
	var thd1 = document.createElement('th');
	var thd2 = document.createElement('th');
	thd1.appendChild(document.createTextNode("Bins"));
	thd2.appendChild(document.createTextNode("Counts"));
	th1.appendChild(thd1);
	th1.appendChild(thd2);
	newTab.appendChild(th1);
	bn.forEach(function(cv){
		var newRow = document.createElement('tr');
		var Td1 = document.createElement('td');
		var Td2 = document.createElement('td');
		Td1.appendChild(document.createTextNode(bn[i]));
		Td2.appendChild(document.createTextNode(cnt[i]));
		newRow.appendChild(Td1);
		newRow.appendChild(Td2);
		newTab.appendChild(newRow);
		i++;
	});
	newTab.className = "tbl";
	document.body.appendChild(newTab);
}
//////////////////////////////////////////

/////////////////////////////////////////////////
//Function 8: "getNfs" calculates damage based on material data and strain value
///Input 1: "ymod" Young's modulus
///Input 2: "sigfp" Fatigue strength coefficient
///Input 3: "bconst" Fatiuge strength exponent
///Input 4: "epsfp" Fatigue ductility coefficient
///Input 5: "cconst" Fatigue ductility exponent
///Input 6: "strnval" strain value in microstrain

function getNfs(ymod, sigfp, bconst, epsfp, cconst, strnval){
	var numf = 1;
	var numx = 10**18;
	var i = 1;
	var epst0 = strlf(ymod, sigfp, bconst, epsfp, cconst, numf);
	var epst1 = strlf(ymod, sigfp, bconst, epsfp, cconst, numx);
	var stp = 10;	
	if(strnval < epst1){
		num2 = numx;
	}else{
		var num2 = (numx + numf)/2;
		epst1 = strlf(ymod, sigfp, bconst, epsfp, cconst, num2);
		while(Math.abs(epst1 - strnval) > 10){
			if(epst1 > strnval){
				numf = num2; //guess is higher strain; move left bound
			}else{
				numx = num2; //guess is lower strain; move right bound
			}
			num2 = (numx + numf)/2;		
			epst1 = strlf(ymod, sigfp, bconst, epsfp, cconst, num2);
		}
	}	
	return num2;
}


function strlf(ymod, sigfp, bconst, epsfp, cconst, numf){
var x = 2*numf;
var a = sigfp/ymod;
var epse = a*(Math.pow(x,bconst));
var epsp = epsfp*(Math.pow(x,cconst));
var tmp = (epse + epsp)*1000000;
return tmp;

}

function getStrLfPlt(ymod, sigfp, bconst, epsfp, cconst){
var nfcs = [];
var nfstr = [];
for(var i=0; i<16; i++){
	nfcs.push(10**i);
}

nfcs.forEach(function(cv){
var tmpsta = (sigfp/ymod)*((2*cv)**bconst);
var tmpstb = epsfp*((2*cv)**cconst);
var tmpst = (tmpsta + tmpstb)*1000000; //microstrain
nfstr.push(tmpst);
});

return [nfcs, nfstr];  //(x,y)

}
