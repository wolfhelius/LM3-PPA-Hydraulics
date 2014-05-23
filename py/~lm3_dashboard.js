<!DOCTYPE html>
<meta charset="utf-8">

<style>
	h1 { 
	  background-color: #2CA25F;
      font: 24px sans-serif;
      margin-top: 0;
      margin-bottom: 0;
	}

 	h2 {
 	  background-color: #99D8C9;
	  font: 16px sans-serif;
      margin-top: 0;
      margin-bottom: 0;
	}
	
 	body {
	  background-color: #E5F5F9;
	  font: 12px sans-serif;
      margin-top: 0;
      margin-bottom: 0;
	}

	sidebar {
 	    background-color: #E5F5F9;
		font: 12px sans-serif;
		position: fixed;
		left: 0px
		top: 0px
		width: 100px;
	}
		
	.axis path,
	.axis line {
	  fill: none;
	  stroke: #000;
	  shape-rendering: crispEdges;
	}

// 	.x.axis path {
// 	  display: none;
// 	}

	svg {
	  font: 12px sans-serif;
	}
	
	div #hoverText {
		font: 14px sans-serif;
		position: absolute;
		width: auto;
		height: auto;
		padding: 10px;
		background-color: white;
		border-radius: 10px;
		box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
		pointer-events: none;
		opacity: 0;
	}


	.line {
	  fill: none;
	  stroke-width: 2.0px;
	 // hover behavior
	  -moz-transition: all 0.3s;
      -o-transition: all 0.3s;
      -webkit-transition: all 0.3s;
      transition: all 0.3s;
	}

	.line:hover {
	  stroke-width: 3.0px;
	  filter: invert(100%);
	  box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
	}	
	
	div #BxN_L_t {
// 		clear: both;
// 		width:100%;
// 		height:300px;
	}

</style>


<body>
	<script src="http://mbostock.github.com/d3/d3.js?2.7.2"></script>
	<script src="http://cdn.craig.is/js/mousetrap/mousetrap.min.js"></script>
	<script src="http://underscorejs.org/underscore-min.js"></script>
	<script type="text/javascript" src="http://code.jquery.com/jquery-1.6.2.min.js"></script> 
	
	<!-- START DOM -->

	<div id="header" >
		<h1>LM3 Visualization</h1>
	</div>
	
	<div id="sidebar">
		<h2>LM3 Runs</h2>
	</div>
	
	<div id="body">
		<h2>LM3 Plots</h2>
		
		<div id="BxN_L_time">
		<!-- Biomass*N lines versus time --!>
		</div>

		<div id="bliving_L_time">
		<!-- Biomass*N stacked versus time --!>
		</div>

		<div id="BxN_S_time">
		<!-- Biomass*N line versus time --!>
		</div>

		<div id="BA_L_time">
		<!-- Basal area stacked versus time --!>
		</div>

		<div id="nindivs_L_time">
		<!-- Number of individuals versus time --!>
		</div>

		<div id="age_L_time">
		<!-- Age lines versus time --!>
		</div>

		<div id="height_B_CAxN">
		<!-- Height versus Crown Area * N --!>
		</div>

	</div>

	<div id="footer">	
	</div>

	<!-- END DOM -->
	

	<script type="text/javascript">	
		function output(inp) {
			document.body.appendChild(document.createElement('pre')).innerHTML = inp;
		}		
		
		// Globals
		var tileJSON;
		var BxN = [], BA = [], CAxN = [], time = [], bliving = [], age = [], nindivs = [], height = [];
		
		var color = d3.scale.category10();
		
		d3.json("/tile.json", function(json) {
		  tileJSON = json;
		  
		  //these go in the call to ensure they run only after data has been loaded
		  preprocess();
		  addLinLinePlot("time", "BxN", "Biomass per cohort (kg C/m2)");
		  addLinLinePlot("time", "bliving", "Biomass per tree (kg C)");
		  addLinLinePlot("time","age", "Age per cohort");
		  addLinLinePlot("time","nindivs", "Pop density (#/m2)");
		  addLinLinePlot("time","BA", "Basal Area (m2/ha)");
		  //addBarPlot("CAxN","height","Total Crown Area");
		});

		
		function preprocess(){		
			//delete tileJSON["_id"];
						
			// pre-process whole-cohort variables
			for (var cohort in tileJSON){ 
				if (cohort != "_id" & cohort != "species_meta"){
					var nY = tileJSON[cohort].age.length;
					var sY = tileJSON[cohort].startyear;
					var timeVec = d3.range(sY, sY+nY);
					tileJSON[cohort].time = d3.range(sY, sY+nY);
					tileJSON[cohort].BxN = [];
					tileJSON[cohort].BA = [];
					tileJSON[cohort].CAxN = [];
					for (var iY = 0; iY<nY; iY++){
						tileJSON[cohort].BxN[iY] = tileJSON[cohort].bliving[iY]*tileJSON[cohort].nindivs[iY];
						tileJSON[cohort].BA[iY] = Math.pow((tileJSON[cohort].dbh[iY])/2,2)*tileJSON[cohort].nindivs[iY]*10000;
						tileJSON[cohort].CAxN[iY] = tileJSON[cohort].crownarea[iY]*tileJSON[cohort].nindivs[iY];	
					}		
				}		
			}	
			for (var cohort in tileJSON){
				if (cohort != "_id" & cohort != "species_meta"){
					BxN = BxN.concat(tileJSON[cohort].BxN);
					BA = BA.concat(tileJSON[cohort].BA);
					CAxN = CAxN.concat(tileJSON[cohort].CAxN);
					bliving = bliving.concat(tileJSON[cohort].bliving);
					age = age.concat(tileJSON[cohort].age);
					nindivs = nindivs.concat(tileJSON[cohort].nindivs);
					time = time.concat(tileJSON[cohort].time);
					height = height.concat(tileJSON[cohort].height);
				}
			}
			time = d3.range(d3.min(time), d3.max(time)+1).map(function(d){return d*1});;
		}
		
		function addLogLinePlot(xName, yName, title){
				
			var margin = {top: 20, right: 20, bottom: 30, left: 50},
				width = 640 - margin.left - margin.right,
				height = 300 - margin.top - margin.bottom;

			//var parseDate = d3.time.format("%d-%b-%y").parse;

			var x = d3.time.scale()
				.range([0, width])
				.domain(d3.extent(window[xName])).nice();

			var y = d3.scale.log()
				.range([height, 0])
				.domain(d3.extent(window[yName])).nice();

			var xAxis = d3.svg.axis()
				.scale(x)
				.orient("bottom");

			var yAxis = d3.svg.axis()
				.scale(y)
				.orient("left");

			var tag = '#'.concat(yName,'_L_', xName);

			var svg = d3.select(tag).append("svg")
				.attr("width", width + margin.left + margin.right)
				.attr("height", height + margin.top + margin.bottom)
				.append("g")
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

			  svg.append("g")
				  .attr("class", "x axis")
				  .attr("transform", "translate(0," + height + ")")
				  .call(xAxis);

			  svg.append("g")
				  .attr("class", "y axis")
				  .call(yAxis)
				  .append("text")
				  .attr("transform", "rotate(-90)")
				  .attr("y", 6)
				  .attr("dy", ".71em")
				  .style("text-anchor", "end")
				  .text(title);


				var line = d3.svg.line()
					.x(function(d) { return x(d.x); })
					.y(function(d) { return y(d.y); });


				// note there is a distinction between .data (dynamic) and .datum (static)
				// https://github.com/mbostock/d3/wiki/Selections#wiki-datum
				for (var cohort in tileJSON){ 
					if (cohort != "_id" & cohort != "species_meta"){
						var cc_data = []; 
						
						for (i = 0; i<tileJSON[cohort].age.length; i++){
							var d = {};
							d.x = tileJSON[cohort][xName][i];
							d.y = tileJSON[cohort][yName][i];
							cc_data.push(d);
						}
						
						var hoverText = "Species: ".concat(tileJSON[cohort].species.toString());
						var hover = d3.select(tag)
							.append("div")
							.style("position", "absolute")
							.style("z-index", "10")
							.style("visibility", "hidden")
							.text(hoverText).attr("id", "hoverText");							
							
						svg.append("path")
							  .data([cc_data])
							  .attr("class", "line")
							  .attr("id", cohort)
							  .attr("d", line)
							  .attr("stroke", color(tileJSON[cohort].species))
							  .on("mouseover", function(){return hover.style("visibility", "visible");})
							  .on("mousemove", function(){return hover.style("top", (event.pageY-10)+"px").style("left",(event.pageX+10)+"px");})
							  .on("mouseout", function(){return hover.style("visibility", "hidden");});

						}
				}

		}
		function addLinLinePlot(xName, yName, title){
				
			var margin = {top: 20, right: 20, bottom: 30, left: 50},
				width = 640 - margin.left - margin.right,
				height = 300 - margin.top - margin.bottom;

			//var parseDate = d3.time.format("%d-%b-%y").parse;

			var x = d3.time.scale()
				.range([0, width])
				.domain(d3.extent(window[xName])).nice();

			var y = d3.scale.linear()
				.range([height, 0])
				.domain(d3.extent(window[yName])).nice();

			var xAxis = d3.svg.axis()
				.scale(x)
				.orient("bottom");

			var yAxis = d3.svg.axis()
				.scale(y)
				.orient("left");

			var tag = '#'.concat(yName,'_L_', xName);

			var svg = d3.select(tag).append("svg")
				.attr("width", width + margin.left + margin.right)
				.attr("height", height + margin.top + margin.bottom)
				.append("g")
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

			  svg.append("g")
				  .attr("class", "x axis")
				  .attr("transform", "translate(0," + height + ")")
				  .call(xAxis);

			  svg.append("g")
				  .attr("class", "y axis")
				  .call(yAxis)
				  .append("text")
				  .attr("transform", "rotate(-90)")
				  .attr("y", 6)
				  .attr("dy", ".71em")
				  .style("text-anchor", "end")
				  .text(title);


				var line = d3.svg.line()
					.x(function(d) { return x(d.x); })
					.y(function(d) { return y(d.y); });


				// note there is a distinction between .data (dynamic) and .datum (static)
				// https://github.com/mbostock/d3/wiki/Selections#wiki-datum
				for (var cohort in tileJSON){ 
					if (cohort != "_id" & cohort != "species_meta"){
						var cc_data = []; 
						
						for (i = 0; i<tileJSON[cohort].age.length; i++){
							var d = {};
							d.x = tileJSON[cohort][xName][i];
							d.y = tileJSON[cohort][yName][i];
							cc_data.push(d);
						}
						
						var hoverText = "Species: ".concat(tileJSON[cohort].species.toString());
						var hover = d3.select(tag)
							.append("div")
							.style("position", "absolute")
							.style("z-index", "10")
							.style("visibility", "hidden")
							.text(hoverText).attr("id", "hoverText");							
							
						svg.append("path")
							  .data([cc_data])
							  .attr("class", "line")
							  .attr("id", cohort)
							  .attr("d", line)
							  .attr("stroke", color(tileJSON[cohort].species))
							  .on("mouseover", function(){return hover.style("visibility", "visible");})
							  .on("mousemove", function(){return hover.style("top", (event.pageY-10)+"px").style("left",(event.pageX+10)+"px");})
							  .on("mouseout", function(){return hover.style("visibility", "hidden");})

						}
				}

		}
		function addBarPlot(xName, yName, title){
				
			var margin = {top: 20, right: 20, bottom: 30, left: 50},
				width = 640 - margin.left - margin.right,
				height = 300 - margin.top - margin.bottom;

			//var parseDate = d3.time.format("%d-%b-%y").parse;

			var x = d3.time.scale()
				.range([0, width])
				.domain(d3.extent(window[xName])).nice();

			var y = d3.scale.linear()
				.range([height, 0])
				.domain(d3.extent(window[yName])).nice();

			var xAxis = d3.svg.axis()
				.scale(x)
				.orient("bottom");

			var yAxis = d3.svg.axis()
				.scale(y)
				.orient("left");

			var tag = '#'.concat(yName,'_B_', xName);

			var svg = d3.select(tag).append("svg")
				.attr("width", width + margin.left + margin.right)
				.attr("height", height + margin.top + margin.bottom)
				.append("g")
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

			  svg.append("g")
				  .attr("class", "x axis")
				  .attr("transform", "translate(0," + height + ")")
				  .call(xAxis);

			  svg.append("g")
				  .attr("class", "y axis")
				  .call(yAxis)
				  .append("text")
				  .attr("transform", "rotate(-90)")
				  .attr("y", 6)
				  .attr("dy", ".71em")
				  .style("text-anchor", "end")
				  .text(title);


				var line = d3.svg.line()
					.x(function(d) { return x(d.x); })
					.y(function(d) { return y(d.y); });


				// note there is a distinction between .data (dynamic) and .datum (static)
				// https://github.com/mbostock/d3/wiki/Selections#wiki-datum
				for (var cohort in tileJSON){ 
					var leftPos = 0;
					if (cohort != "_id" & cohort != "species_meta" ){
						var cc_data = []; 
						
						for (i = 0; i<tileJSON[cohort].age.length; i++){
							var d = {};
							d.w = tileJSON[cohort][xName][i];
							d.y = tileJSON[cohort][yName][i];
							if (i==0){						
								d.x = 0;
							} else {
								d.x = cc_data[i-1].x+cc_data[i-1].w;
							}
							cc_data.push(d);
						}
						
						var hoverText = "Species: "+tileJSON[cohort].species.toString()+"\n"
										"Cohort: "+ cohort.toString();
						var hover = d3.select(tag)
							.append("div")
							.style("position", "absolute")
							.style("z-index", "10")
							.style("visibility", "hidden")
							.text(hoverText).attr("id", "hoverText");							
							
						svg.append("rect")
							  .data([cc_data])
							  .attr("class", "line")
							  .attr("id", cohort)
							  .attr("x", function(d) {return x(leftPos);})
							  .attr("width", function(d) {return x(d.x);})
							  .attr("y", function(d) {return y(d.y);})
							  .attr("height", function(d) {return height - y(d.y);})							  
							  .attr("fill", color(tileJSON[cohort].species))
							  .attr("stroke", "black")
							  .attr("stroke-width", "1px")
							  .on("mouseover", function(){return hover.style("visibility", "visible");})
							  .on("mousemove", function(){return hover.style("top", (event.pageY-10)+"px").style("left",(event.pageX+10)+"px");})
							  .on("mouseout", function(){return hover.style("visibility", "hidden");});

						}
						leftPos = leftPos+cc_data						
				}

		}
		
	</script>
</body>		 