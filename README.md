A transcriptome browser for Ribo-Seq and RNA-Seq data

Trips-viz is a transcriptome browser designed to visualize Ribosome profiling and RNA-seq data at the level of a single gene/transcript isoform as opposed to at the genome level. Trips-viz also provides you the ability to vizualize data from a gene under different conditions, get meta-information at an individual dataset level such as read length distribution or triplet periodicity, and provides the ability to find differentially expressed or translated genes.

To run Trips-Viz locally follow these steps:
1. rename trips_empty_db.sqlite to trips.sqlite, 
2. rename config_template.py to config.py. 
3. Run the command python __init__.py true 5000 and navigate to 0.0.0.0:5000/ in a browser (if running on a server replace 0.0.0.0 with the server IP). 
