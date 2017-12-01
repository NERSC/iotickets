/* global myUsername, chartHead, chartTail, verifyData, chart_iotimeperc, chart_rwsize, chart_iorate, chart_iosizehist, chart_iosizehistweighted */
/* exported myDarshanCharts */

var queryUser = "";
var originalWidth;
var inJobID;
var lastJob;
var globalResizeTimer;
var months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"];

$( document ).ready(function($) {
  "use strict";


});

$( window ).resize(function() {
  "use strict";
	// the charts need to be rebuilt on page width change or they get messed up
	//		this SHOULD be done using a highcharts method, but this will have
	//		to wait until the whole program is rebuilt using objects and better
	//		scoping practices

  if(globalResizeTimer !== null) { window.clearTimeout(globalResizeTimer); }

  globalResizeTimer = window.setTimeout(function() {

		var width = document.documentElement.clientWidth;

		if( width !== originalWidth ) {

			if( lastJob !== null ) {

				build(lastJob);

			} else {

				buildChain();

			}

			originalWidth = width;

		}

  }, 1000);

});

var divArray = [

	// this is mostly for populating the page with ajax-loader gifs

	"#darshan-charts-super-header",
	"#darshan-charts-range-slider",
	"#darshan-charts-job-selector",
	"#darshan-charts-left1",
	"#darshan-charts-left2",
	"#darshan-charts-right1",
	"#darshan-charts-right2",
	"#darshan-charts-bottom1",
	"#darshan-charts-bottom2",
	"#darshan-charts-bottom3"

];

function getNavPlaceholderString() {
  "use strict";

  var myString = '<div id="date-ranger" style="display:none;">Date Range:<div class="head-room"></div>';
  myString += '<center><div id="darshan-charts-range-slider"></div></center>';
  myString += '<div class="elbow-room"></div></div>';

  myString += 'Individual Task or Aggregate:<div '/*class="col-lg-10"*/+' id="darshan-charts-job-selector" style="padding-top:8px;"></div>';

  return myString;

}

function reportError(str) {
  "use strict";

	console.log(str);

	for(var i = 2; i < divArray.length; i++) {
		$(divArray[i]).html('');
	}

	var myString = "";
	myString += chartHead(' ERROR');
	myString += str;
	myString += chartTail();

	$("#darshan-charts-error").html(myString);

}

function reportNotice(str) {
  "use strict";
	console.log(str);

	for(var i = 3; i < divArray.length; i++) {
		$(divArray[i]).html('');
	}

	var myString = "";
	myString += chartHead('NOTICE');
	myString += str;
	myString += chartTail();

	$("#darshan-charts-notice").html(myString);

	//alert(str);

}

function displayLoader() {
  "use strict";
	var myString = "";
	myString += chartHead('WORKING');
	myString += '<center><img src="images/blue-ajax-loader.gif"></center>';
	myString += chartTail();

	$("#darshan-charts-loader").html(myString);

}

function clearLoader() {
  //alert("clear loader")
  "use strict";
  $("#darshan-charts-loader").html('');
  //alert("clear loader 2")

}

function clearAllMessages() {
  "use strict";
  $("#darshan-charts-error").html('');
  $("#darshan-charts-notice").html('');

}

function myDarshanCharts() {

  "use strict";

  lastJob = null;
  globalResizeTimer = null;
  $("#darshan-charts-nav-desktop").html(getNavPlaceholderString());

  clearAllMessages();

  queryUser = myUsername;

  if (myUsername.indexOf("invalid") !== -1) {

      clearAllMessages();

      lastJob = null;
      globalResizeTimer = null;
      reportError("<a class=\"btn btn-primary\" href=\"login.php\">Please Login</a>");

  } else {

    inJobID = getUrlVars()["jobid"];
    console.log("Starting Darshan Charts "+inJobID);
    if (inJobID) {
      $("#myheader").html("IO Statistics for "+inJobID)
      $('#date-ranger').hide()
    } else {
      $('#date-ranger').show()
    }
    $('#job-well').show();
    testMongoConnection();

  }

}

function testMongoConnection() {

  "use strict";

	displayLoader();

	$.ajax({

		url: "non-newt-queries/mongo-test.php",
		type: "GET",
		dataType: "json",

		success: function(res) {

			var connected = res.connected;

			// this needs to be fixed up
			if( !connected ) {

				reportError("MongoDB error: cannot establish a connection");

			} else {

				buildChain();

			}

		},

		error: function(xhr, status, error) {

			reportError(error);

		}

	});

}

function buildChain() {

  "use strict";
  buildRangeSlider();

}

function dateToString(dateObj) {

	var month = parseInt(dateObj.getUTCMonth()) * 1 + 1;
	var dateString = dateObj.getUTCFullYear() + "-" + ("0" + month).slice(-2) + "-" + ("0" + dateObj.getUTCDate()).slice(-2);
	dateString += " " + ("0" + dateObj.getUTCHours()).slice(-2) + ":" + ("0" + dateObj.getUTCMinutes()).slice(-2);
	return dateString;

}

function buildRangeSlider() {

  "use strict";
	clearAllMessages();

	$.ajax({

		url: "non-newt-queries/time-range.php?username="+queryUser,
		type: "GET",
		dataType: "json",

		success: function(res) {

			if(res['oldest'] === null && res['newest'] === null) {

				reportError("No jobs found for user: "+queryUser);

			}

			var myString = "<div id=\"slider\"></div>";

			$("#darshan-charts-range-slider").html(myString);

                        console.log("About to define oldest "+res['oldest']+" "+res['newest'])

			var oldest = new Date(1356998400000);
                        if ( res['oldest'] !== null) { 
                          oldest = new Date(res['oldest']*1000);
                        }
			oldest.setUTCHours(0);
			oldest.setUTCMinutes(0);
			oldest.setUTCSeconds(0);

			var newest = new Date(res['newest']*1000);

			$("#slider").dateRangeSlider(
                                { scales: [{
                                  first: function(value){ return value; },
                                  end: function(value) {return value; },
                                  next: function(value){
                                    var next = new Date(value);
                                    return new Date(next.setMonth(value.getMonth() + 1));
                                  },
                                  label: function(value){
                                    return months[value.getMonth()];
                                  },
                                  format: function(tickContainer, tickStart, tickEnd){
                                    tickContainer.addClass("myCustomClass");
                                   }
                                 }]
                                },
				{ arrows: false },
				{ bounds: {min: oldest, max: newest}},
				{ defaultValues: {min : oldest, max: newest}},
				{ step: {days: 1}}
			);

			populateJobSelector(res['oldest'],res['newest']);

			$("#slider").bind("valuesChanged", function(e, data){

				displayLoader();

				$("#darshan-charts-job-selector").html('');

				var minTime = data.values.min;
				minTime.setUTCHours(0);
				minTime.setUTCMinutes(0);
				minTime.setUTCSeconds(0);
				var maxTime = data.values.max;
				maxTime.setUTCHours(23);
				maxTime.setUTCMinutes(59);
				maxTime.setUTCSeconds(59);

				minTime = minTime.getTime()/1000;
				maxTime = maxTime.getTime()/1000;
				populateJobSelector(minTime, maxTime);

			});

		},

		error: function(xhr, status, error) {

			reportError(error);

		}

	});

}

function populateJobSelector(oldest, newest) {

  "use strict";

	var queryString = "non-newt-queries/jobs.php?username="+queryUser;

        if (inJobID) {
		queryString += "&jobid=" + inJobID;
	} else if(oldest && newest) {
		queryString += "&oldest=" + oldest + "&newest=" + newest;
	}

	$.ajax({

		url: queryString,
		type: "GET",
		dataType: "json",

		success: function(res) {

		  if (res.length < 1) {

                        //alert("In populate")
                        //$('#job-well').hide();
                        clearLoader();
		        reportError("No darshan data found.");

	          } else {

			var myString = "";
			myString += "<select class=\"selectpicker show-menu-arrow show-tick\" data-width=\"100%\" id=\"job-selector\">";

                        if (inJobID) {
			  myString += "<option value=\"" + "!range! " + oldest + " " + newest + "\">" + "aggregate from from all "+inJobID+" tasks</option>";
                        } else {
			  myString += "<option value=\"" + "!range! " + oldest + " " + newest + "\">" + "aggregate from " + dateToString(new Date(oldest*1000)) + " to " + dateToString(new Date(newest*1000)) + "</option>";
                        } 

			for(var i = 0; i < res.length; i++) {

				var myDate = dateToString(new Date(res[i]['end_time']*1000));
				myString += "<option value=\"" + res[i]['_id'] + "\">" + res[i]['host'] + ", " + res[i]['jobid'] + ", " + res[i]['exe'].replace(/^.*[\\\/]/, '') + ", " + myDate + "</option>";

			}

			myString += "</select>";
			$("#darshan-charts-job-selector").html(myString);

                        $('#job-selector').selectpicker();

			queryAggregate(build, oldest, newest);

			$('#job-selector').change(function() {

				var selectValue = $(this).val().split(" ");
				console.log(selectValue);

				setTimeout(function() {
					
					if( selectValue[0] === '!range!' ) {
						queryAggregate(build, selectValue[1], selectValue[2]);
					} else {
						query(build, selectValue[0]);
					}

				},250);

			});
                  }
		},

		error: function(xhr, status, error) {

			reportError(error);

		}

	});

}

function query(handle, job) {

  console.log('In query '+job);

  "use strict";
	$.ajax({

		url: "non-newt-queries/get-job.php?_id=" + job,
		type: "GET",
		dataType: "json",

		success: function(res) {

			lastJob = res;
			handle(res);

		},

		error: function(xhr, status, error) {

			reportError(error);

		}

	});

}

function queryAggregate(handle, from, until) {

  "use strict";

  var myURL = "non-newt-queries/get-aggregate.php?username="+queryUser+"&oldest="+from+"&newest="+until;
  if (inJobID) {
    myURL = "non-newt-queries/get-aggregate.php?username="+queryUser+"&jobid="+inJobID;
  } 
 
  $.ajax({

		url: myURL,
		type: "GET",
		dataType: "json",

		success: function(res) {

			lastJob = res;
			handle(res);

		},

		error: function(xhr, status, error) {

			reportError(error);

		}

  });

}

var build = function buildCharts(json) {

  "use strict";
	clearAllMessages();

	if( verifyData(json) ) {

		var firstChartDiv = 3;

		/*for(var i = firstChartDiv; i < divArray.length; i++) {

			$(divArray[i]).html('<center><img src="images/blue-ajax-loader.gif"></center>');

		}*/

		displayLoader();

		setTimeout(function() {

			clearLoader();
			buildLeft1(firstChartDiv, json);
			//buildLeft2(firstChartDiv, json);
			buildRight1(firstChartDiv, json);
			buildRight2(firstChartDiv, json);
			//buildBottom1(firstChartDiv, json);
			buildBottom2(firstChartDiv, json);
			buildBottom3(firstChartDiv, json);

		}, 250);

	} else {

		reportNotice("Unviable job data.");

	}

};

function buildLeft1(firstDiv, json) {

  "use strict";
	var element = divArray[firstDiv];

	var myString = "";
	myString += chartHead(' Percentage Time Spent (I/O vs Computation)');
	myString += "<div id=\"chart-iotimeperc\" class=\"center\"></div>";
	myString += chartTail();

	$(element).html(myString);

	chart_iotimeperc(json, 'chart-iotimeperc');

}

function buildRight1(firstDiv, json) {

  "use strict";
	var element = divArray[firstDiv+2];

	var myString = "";
	myString += chartHead(' Average POSIX I/O Amount (Read/Written GB)');
	myString += "<div id=\"chart-rwsize\"></div>";
	myString += chartTail();

	$(element).html(myString);

	chart_rwsize(json, '#chart-rwsize');

}

function buildRight2(firstDiv, json) {

  "use strict";
	var element = divArray[firstDiv+3];

	var myString = "";
	myString += chartHead(' Estimated IO Rate');
	myString += "<div id=\"chart-iorate\"></div>";
	myString += chartTail();

  $(element).html(myString);

	chart_iorate(json, '#chart-iorate');

}

function buildBottom2(firstDiv, json) {

  "use strict";
	var element = divArray[firstDiv+5];

	var myString = "";
	myString += chartHead(' Number of Transactions in Each Size Range');
	myString += "<div id=\"chart-iosizehist\"></div>";
	myString += chartTail();

	$(element).html(myString);

	chart_iosizehist(json, '#chart-iosizehist');
	//chart_iosizehist(json, element);

}

function buildBottom3(firstDiv, json) {

  "use strict";
	var element = divArray[firstDiv+6];

	var myString = "";
	myString += chartHead(' Amount Written With Different Transaction Sizes');
	myString += "<div id=\"chart-iosizehistweighted\"></div>";
	myString += chartTail();

	$(element).html(myString);

	chart_iosizehistweighted(json, '#chart-iosizehistweighted');
	//chart_iosizehistweighted(json, element);

}

