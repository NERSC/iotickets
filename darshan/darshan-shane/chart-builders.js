// Courtesy of: Shane Pearlman

/* global Morris, Highcharts */
/* exported chartHead, chartTail, verifyData, chart_iotimeperc, chart_iorate, chart_rwsize, chart_iosizehist, chart_iosizehistweighted */

function toArray(json) {

	"use strict";
	var arr = [];

    arr.push(json['_id']);
    arr.push(json['jobid']);
    arr.push(json['start_time']);
    arr.push(json['end_time']);
    arr.push(json['nprocs']);
    arr.push(json['host']);
    arr.push(json['exe']);
    arr.push(json['uid']);
    arr.push(json['total_CP_BYTES_READ']);
    arr.push(json['total_CP_BYTES_WRITTEN']);
    arr.push(json['total_CP_SIZE_READ_0_100']);
    arr.push(json['total_CP_SIZE_READ_100_1K']);
    arr.push(json['total_CP_SIZE_READ_1K_10K']);
    arr.push(json['total_CP_SIZE_READ_10K_100K']);
    arr.push(json['total_CP_SIZE_READ_100K_1M']);
    arr.push(json['total_CP_SIZE_READ_1M_4M']);
    arr.push(json['total_CP_SIZE_READ_4M_10M']);
    arr.push(json['total_CP_SIZE_READ_10M_100M']);
    arr.push(json['total_CP_SIZE_READ_100M_1G']);
    arr.push(json['total_CP_SIZE_READ_1G_PLUS']);
    arr.push(json['total_CP_SIZE_WRITE_0_100']);
    arr.push(json['total_CP_SIZE_WRITE_100_1K']);
    arr.push(json['total_CP_SIZE_WRITE_1K_10K']);
    arr.push(json['total_CP_SIZE_WRITE_10K_100K']);
    arr.push(json['total_CP_SIZE_WRITE_100K_1M']);
    arr.push(json['total_CP_SIZE_WRITE_1M_4M']);
    arr.push(json['total_CP_SIZE_WRITE_4M_10M']);
    arr.push(json['total_CP_SIZE_WRITE_10M_100M']);
    arr.push(json['total_CP_SIZE_WRITE_100M_1G']);
    arr.push(json['total_CP_SIZE_WRITE_1G_PLUS']);
    arr.push(json['total_CP_F_POSIX_READ_TIME']);
    arr.push(json['total_CP_F_POSIX_WRITE_TIME']);
    arr.push(json['total_CP_F_POSIX_META_TIME']);
    arr.push(json['agg_perf_by_cumul']);
    arr.push(json['agg_perf_by_open']);
    arr.push(json['agg_perf_by_open_lastio']);
    arr.push(json['agg_perf_by_slowest']);
    arr.push(json['total_bytes']);

	return arr;

}

function verifyData(json) {

  "use strict";
	var data = toArray(json);

	for(var i = 8; i < data.length; i++) {
		if( typeof data[i] !== "number" ) {
			return false;
    }
	}

	return true;

}

// function d3test(element,dataset) {

// 	d3.select(element).selectAll("div")
// 	    .data(dataset)
// 	    .enter()
// 	    .append("div")
// 		.style("width", function(d) { return d * 10 + "px"; })
// 		.text(function(d) { return d; });

// }

function chart_iotimeperc(json, element) {

  "use strict";
	var readTime = json['total_CP_F_POSIX_READ_TIME'];
	var writeTime = json['total_CP_F_POSIX_WRITE_TIME'];
	var metaTime = json['total_CP_F_POSIX_META_TIME'];
        var ioTime = json['total_bytes']/1024/1024/json['agg_perf_by_slowest']
	//var ioTime = readTime+writeTime;
        if (json['io_time']) {
          ioTime = json['io_time'];
        }

        var myStartTime = json['start_time'];
        var myEndTime = json['end_time'];
        var myTotTime = myEndTime - myStartTime + 1;
        if (json['tot_time']) {
          myTotTime = json['tot_time']+1;
        }

        console.log("Drawing Pie start "+myStartTime+" end "+myEndTime+" tot "+myTotTime+" meta "+metaTime+" write "+writeTime+" read "+readTime+" tot "+ioTime);

        var myComputeTime = myTotTime - ioTime;

	var percentIOTime = Math.floor(ioTime/(myTotTime)*100000)/1000;
	//var percentMetaTime = Math.floor(metaTime/(myTotTime)*100000)/1000;
	var percentComputeTime = Math.floor(myComputeTime/(myTotTime)*100000)/1000;

	var metaColor = "#256781";
	var ioColor = "#5977B3";
	var computeColor = "#CC8352";

	//if(!percentIOTime) {

		Morris.Donut({

			element: element,
			colors: [ioColor, computeColor],
			formatter: function(y) {
				return (y + "%");
			},
			data: [
				{label: "IO", value: percentIOTime},
				{label: "Compute", value: percentComputeTime}
			]

		});

	//}

}

function chart_iosizehist(json, element) {

  "use strict";
	var jsonArr = toArray(json);

	var read = [];
	for(var i = 10; i < 20; i++) {
		read.push(jsonArr[i]);
	}

	var written = [];
	for(i = 20; i < 30; i++) {
		written.push(jsonArr[i]);
	}
	
	$(element).highcharts({

		chart: {
	        type: 'column',
			marginLeft: 100
        },

        title: {
            text: ''
        },

        xAxis: {
            categories: ['0-100B', '100B-1KB', '1KB-10KB', '10-100KB', '100KB-1MB',
				'1MB-4MB', '4MB-10MB', '10MB-100MB', '100MB-1GB', '1GB  or More']
        },

        yAxis: {
            allowDecimals: true,
            min: 0,
            title: {
                text: 'number of transactions'
            }
        },

        tooltip: {
            //enabled: false
        },

        plotOptions: {
            column: {
                stacking: null
            }
        },

        series: [{
            name: 'Read',
            data: read
        }, {
            name: 'Written',
            data: written
        }],

		credits: {

			enabled: false

		}

    });

}

function chart_iosizehistweighted(json, element) {

  "use strict";
	var read = [];
    read.push(json['total_CP_SIZE_READ_0_100']/1024/1024/1024*50);
    read.push(json['total_CP_SIZE_READ_100_1K']/1024/1024/1024*(1024+100)/2);
    read.push(json['total_CP_SIZE_READ_1K_10K']/1024/1024/1024*(10*1024+1024)/2);
    read.push(json['total_CP_SIZE_READ_10K_100K']/1024/1024/1024*(100*1024+10*1024)/2);
    read.push(json['total_CP_SIZE_READ_100K_1M']/1024/1024/1024*(1024*1024+100*1024)/2);
    read.push(json['total_CP_SIZE_READ_1M_4M']/1024/1024/1024*(4*1024*1024+1*1024*1024)/2);
    read.push(json['total_CP_SIZE_READ_4M_10M']/1024/1024/1024*(10*1024*1024+4*1024*1024)/2);
    read.push(json['total_CP_SIZE_READ_10M_100M']/1024/1024/1024*(100*1024*1024+10*1024*1024));
    read.push(json['total_CP_SIZE_READ_100M_1G']/1024/1024/1024*(1024*1024*1024+100*1024*1024)/2);
    read.push(json['total_CP_SIZE_READ_1G_PLUS']/1024/1024/1024*(1024*1024*1024*2));

	for(var i = 0; i < read.length; i++) read[i] = Math.floor(read[i]*1000)/1000;

	var written = [];
    written.push(json['total_CP_SIZE_WRITE_0_100']/1024/1024/1024*50);
    written.push(json['total_CP_SIZE_WRITE_100_1K']/1024/1024/1024*(1024+100)/2);
    written.push(json['total_CP_SIZE_WRITE_1K_10K']/1024/1024/1024*(10*1024+1024)/2);
    written.push(json['total_CP_SIZE_WRITE_10K_100K']/1024/1024/1024*(100*1024+10*1024)/2);
    written.push(json['total_CP_SIZE_WRITE_100K_1M']/1024/1024/1024*(1024*1024+100*1024)/2);
    written.push(json['total_CP_SIZE_WRITE_1M_4M']/1024/1024/1024*(4*1024*1024+1*1024*1024)/2);
    written.push(json['total_CP_SIZE_WRITE_4M_10M']/1024/1024/1024*(10*1024*1024+4*1024*1024)/2);
    written.push(json['total_CP_SIZE_WRITE_10M_100M']/1024/1024/1024*(100*1024*1024+10*1024*1024));
    written.push(json['total_CP_SIZE_WRITE_100M_1G']/1024/1024/1024*(1024*1024*1024+100*1024*1024)/2);
    written.push(json['total_CP_SIZE_WRITE_1G_PLUS']/1024/1024/1024*(1024*1024*1024*2));

	for(i = 0; i < written.length; i++) written[i] = Math.floor(written[i]*1000)/1000;
	
	$(element).highcharts({

		chart: {
	        type: 'column',
			marginLeft: 100
        },

        title: {
            text: ''
        },

        xAxis: {
            categories: ['0-100B', '100B-1KB', '1KB-10KB', '10-100KB', '100KB-1MB',
				'1MB-4MB', '4MB-10MB', '10MB-100MB', '100MB-1GB', '1GB  or More']
        },

        yAxis: {
            allowDecimals: true,
            min: 0,
            title: {
                text: 'amount of I/O (GB)'
            }
        },

        tooltip: {
            //enabled: false
        },

        plotOptions: {
            column: {
                stacking: null
            }
        },

        series: [{
            name: 'Read',
            data: read
        }, {
            name: 'Written',
            data: written
        }],

		credits: {

			enabled: false

		}

    });

}

function chart_rwsize(json, element) {

  "use strict";
	var read = Math.ceil(json['total_CP_BYTES_READ']/1024/1024/1024*100)/100;
	var written = Math.ceil(json['total_CP_BYTES_WRITTEN']/1024/1024/1024*100)/100;

    $(element).highcharts({
        chart: {
            type: 'bar',
            height: 300
        },
        title: {
            text: ''
        },
        xAxis: {
			lineWidth: 3,
			labels: {
				enabled: false
			},
            categories: [],
            title: {
                text: null
            }
        },
        yAxis: {
			gridLineWidth: 0,
            min: 0,
            title: {
                text: 'Size (GB)',
                align: 'high'
            },
            labels: {
                overflow: 'justify'
            }
        },
        tooltip: {
            //valueSuffix: ' GB'
			formatter: function() {
				return this.y + ' GB';
			}
        },
        plotOptions: {
			series: {
				animation : true
			},
            bar: {
                dataLabels: {
                    enabled: true
                }
            }
        },
        legend: {
            layout: 'horizontal',
            align: 'right',
            verticalAlign: 'top',
            x: -30,
            y: 25,
            floating: true,
            borderWidth: 0,
            backgroundColor: (Highcharts.theme && Highcharts.theme.legendBackgroundColor || '#FFFFFF'),
            shadow: true
        },
        credits: {
            enabled: false
        },
        series: [{
            name: 'Read',
            data: [read]
        }, {
            name: 'Written',
            data: [written]
        }]
    }); 

}

function chart_iorate(json, element) {

  "use strict";
	var idealrate = 30000;

  var rate = Math.round(json['agg_perf_by_cumul']*100)/100;
  idealrate = Math.round(idealrate*100)/100;

  $(element).highcharts({

    chart: {
      type: 'bar',
      height: 300
    },
        title: {
            text: ''
        },
    xAxis: {
      lineWidth: 3,
      labels: {
        enabled: false
      },
      categories: [],
      title: {
        text: null
      }
    },

    yAxis: {
      gridLineWidth: 0,
      min: 0,
      title: {
        text: 'Rate (MB/Second)',
        align: 'high'
      },
      labels: {
        overflow: 'justify',
      }
    },

    tooltip: {
      //enabled: false
      formatter: function() {
        return this.y * 1 + ' MB/s';
      }
    },

    plotOptions: {
      series: {
        animation : true
      },
      bar: {
        dataLabels: {
          enabled: true,
          formatter: function() {
            return this.y * 1;
          }
        }
      }
    },

    legend: {
      layout: 'horizontal',
      align: 'right',
      verticalAlign: 'bottom',
      x: -30,
      y: -40,
      floating: true,
      borderWidth: 0,
      backgroundColor: (Highcharts.theme && Highcharts.theme.legendBackgroundColor || '#FFFFFF'),
      shadow: true
    },

    credits: {
      enabled: false,
    },

    series: [{
      name: 'Measured',
      color:'#700070',
      data: [(rate)]
    }, {
      name: 'Theoratical Peak',
      color:'#D9D9D9',
      data: [(idealrate)]
    }]

  });

}

function chartHead(title) {

  "use strict";
	return '<div class="panel panel-default"><div class="panel-heading" onClick=""><i class="fa fa-bar-chart-o fa-fw"></i>' + title + '</div><div class="panel-collapse collapse in"><div class="panel-body">';

}

function chartTail() {

  "use strict";
	return '</div></div></div>';

}

