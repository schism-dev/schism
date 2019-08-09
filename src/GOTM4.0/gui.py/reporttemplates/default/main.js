var showdefaults = false;
function showHideDefaults() {
	newdisplaystyle = '';
	newweight = 'bold';
	if (showdefaults) {
		newdisplaystyle = 'none';
		newweight = 'normal';
	}
	
	rows = document.getElementById('tableScenario').rows;
	for (i=0; i<rows.length; i++) {
		tr = rows[i];
		if (tr.getAttribute('default')!=null)
			tr.style.display = newdisplaystyle;
		else
			tr.cells[tr.cells.length-1].style.fontWeight = newweight;
	}
	
	showdefaults = (!showdefaults);
	aShowHide = document.getElementById('aShowHideDefaults');
	if (showdefaults) {
		aShowHide.innerHTML = 'Hide default values';
	} else {
		aShowHide.innerHTML = 'Show default values';
	}
}
