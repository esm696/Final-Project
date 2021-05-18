# Voronoi bar plot
In this work, we need to use OVITO GUI to calculate Voronoi index for the Cu-Zr particle. Here are the discriptions of the files.
    -**ovito.voro.py** This is the python script you need to write in OVITO "Python Script" modifier, and run it after run a "Voronoi analysis" modifier (One example).  <br>
	     -Ovito input file: dump.melt86500000.cfg.gz<br>
		 -Output file: Voronoi_text<br>
	-**voronoianalysis.py** This is the python script used to generate bar plot for 4 different particles.<br>
	     - input files: Voronoi_IGC4, Voronoi_IGC6, Voronoi_IGC10, Voronoi_MQ<br>

**IF NECESSARY, PLEASE CHANGE THE PATH WRITING IN EVERY CODE**
