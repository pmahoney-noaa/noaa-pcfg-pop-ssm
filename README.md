# Projecting PCFG Gray Whale Abundance


<script src="README_files/libs/htmlwidgets-1.6.4/htmlwidgets.js"></script>
<link href="README_files/libs/datatables-css-0.0.0/datatables-crosstalk.css" rel="stylesheet" />
<script src="README_files/libs/datatables-binding-0.33/datatables.js"></script>
<script src="README_files/libs/jquery-3.6.0/jquery-3.6.0.min.js"></script>
<link href="README_files/libs/dt-core-1.13.6/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="README_files/libs/dt-core-1.13.6/css/jquery.dataTables.extra.css" rel="stylesheet" />
<script src="README_files/libs/dt-core-1.13.6/js/jquery.dataTables.min.js"></script>
<link href="README_files/libs/crosstalk-1.2.1/css/crosstalk.min.css" rel="stylesheet" />
<script src="README_files/libs/crosstalk-1.2.1/js/crosstalk.min.js"></script>


<div class="cell-output-display">

<div class="datatables html-widget html-fill-item" id="htmlwidget-11050a4084b48a7bc9fd" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-11050a4084b48a7bc9fd">{"x":{"filter":"none","vertical":false,"data":[["Base","Base","Base","AR1v1","AR1v1","AR1v1","AR1v2","AR1v2","AR1v2","ENP Calves","ENP Calves","ENP Calves","Calves/Strandings","Calves/Strandings","Calves/Strandings","Calves only","Calves only","Calves only","Strandings only","Strandings only","Strandings only"],[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.4650022255284227,0.5112246517637251,0.3335885008226741,0.3499817215475712,0.3279981815654215,0.3087872011108891,0.5368877267319965,0.5308277648321691,0.5124110423003756,0.6496370370936217,0.5872171327923503,0.684939989791537,0.3269771940532265,0.3750481729183305,0.4603830349977082,0.4822992988924115,0.5007867455542321,0.5166107326644891,0.3509242897107823,0.3535176516209442,0.3763357961568637]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Model<\/th>\n      <th>chain<\/th>\n      <th>num_divergent<\/th>\n      <th>num_max_treedepth<\/th>\n      <th>ebfmi<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"className":"dt-center","targets":[1,2,3,4]},{"name":"Model","targets":0},{"name":"chain","targets":1},{"name":"num_divergent","targets":2},{"name":"num_max_treedepth","targets":3},{"name":"ebfmi","targets":4}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render"],"jsHooks":[]}</script>

</div>

## Leave-one-out Cross Validation (LOO)

## Retrospective analysis

![](README_files/figure-commonmark/retroFig-1.png)

## Model-specific trends and projections

Note, the number of PCFG calves used in projected years (models
Calves/Strandings and Calves only) are not accurate and likely represent
underestimates of reality.

![](README_files/figure-commonmark/trendFig-1.png)

## Model-specific predictions for Y<sub>final</sub> + 2

![](README_files/figure-commonmark/fig-final-proj-1.png)
