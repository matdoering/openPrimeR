// Returns a function, that, as long as it continues to be invoked, will not
// be triggered. The function will be called after it stops being called for
// N milliseconds. If `immediate` is passed, trigger the function on the
// leading edge, instead of the trailing.
function debounce(func, wait, immediate) {
	var timeout;
	return function() {
		var context = this, args = arguments;
		var later = function() {
			timeout = null;
			if (!immediate) func.apply(context, args);
		};
		var callNow = immediate && !timeout;
		clearTimeout(timeout);
		timeout = setTimeout(later, wait);
		if (callNow) func.apply(context, args);
	};
}; 
function showStuff(id) {
    document.getElementById(id).style.display = "block";
}
function hideStuff(id) {
    document.getElementById(id).style.display = "none";
}

function my_debounce(func, wait) {
    // make sure we don't show the spinner everytime we're set to busy
    // let's wait <wait> time before calling spinner
    var timeout = null;
    if (timeout) {
        clearTimeout(timer); //cancel the previous timer.
        timeout = null;
    }
    timeout = setTimeout(func, wait);
};
var showSpinnerLazy = function() {
    if ($('html').hasClass('shiny-busy')) {
        showStuff("progressContainer");
        showStuff("progressBlock");
        showStuff("progressText");
        $('html').addClass('my-shiny-busy');
    }
}
var hideSpinner = function() {
    if ($('html').hasClass('my-shiny-busy')) {
        $('html').removeClass('my-shiny-busy');
        hideStuff("progressContainer");
        hideStuff("progressBlock");
        hideStuff("progressText");
    }
}
$(document).on('shiny:busy', function(event) {
    my_debounce(showSpinnerLazy, 2000);
});

//$(document).on('shiny:busy', function(event) {
 //   hideSpinner();
//});    

//'html').attr('class')=='shiny-busy'
$(document).on('shiny:idle', function(event) {
    hideSpinner();
});    

