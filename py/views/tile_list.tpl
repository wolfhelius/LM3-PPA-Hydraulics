<html>
<style>
 	body {
	  font: 12px sans-serif;
      margin-top: 0;
      margin-bottom: 0;
	}
</style>

<head>
	<title>LM3 Runs</title>
</head>
<body>
% for key in tiles.keys():
	<li><a href=http://localhost:8080/lm3_dashboard/{{key}}>{{tiles[key]}}: {{key}}</a>
	<a href=http://localhost:8080/delete/{{key}}>delete</a>	
	</li>
% end
</body>

</html>