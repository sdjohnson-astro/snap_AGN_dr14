The following contains the project code catalogs, and data for SDSS quasar-quasar pairs

Python programs are in code/
Catalogs are in catalogs/
spectra data are in spectra/fg and spectra/bg for foreground and background quasars.

Task 1: check foreground quasar spectra
	The first task in the project is to check the classification and redshift
	measurements for quasars from the SDSS dr14 quasar catalog (Paris+2014)
	
	To do so, cd into the code directory "cd code" and then start the spectrum
	checking user interface for foreground quasars with the command in the terminal:
		"python checkForeground.py"
	
	This will open up a user interface for checking foreground quasar spectra.
	It plots each quasar spectrum with flux in white, a model in red, and error
	array in blue.
	
	Vertical lines mark the expected positions of emission lines at the quasar redshift.
	
	The goal of this plot is to determine whether the quasar is acceptable or not.
	To be acceptable both:
		(1) The red model should be a good fit to the white data.
		and
		(2) The vertical lines should be aligned with observed emission lines in
		    the white data.
			 
	
	To smooth or unsmooth the spectra by typing '+' or '-'.
	
	To zoom in/out or pan, use the mouse.
	
	The title of each plot shows the index, plate-MJD-fiberID, and quality class.
	
	To mark a quasar as acceptable, type 'a', to mark it as rejected, type 'r'.
	If you are uncertain, type 'u'.
	
	Make sure that that the quality classification has updated in the plot title.
	
	After marking the classification, go to the next spectrum by typing 'n'.
	If you need to go back, type 'b'.
	
	When you go to the next spectrum, the interface automatically saves the new
	classification to ../catalogs/pairs.fits
	
	
	There are a total of almost 12,000 of these to check so getting through all
	of them will take a while. Once you get the hang of it, it should move pretty
	quickly though. Something like 3 seconds per quasar. That means a total of
	about 10 hours. It is a tedious and boring process so my advice is to break
	it up into ~10 minute chunks and continue going through programming lesson
	between chunks. 
	
	You can close out the program and then restart at the index you left off at.
	For example if you went through index 1005 then you can start the interface
	at that index by typing "python checkForeground.py -i 1005".
	
	To convert the updated pair catalog with classification into a human readable
	text file type "python convertPairTable.py". This is useful if you lose track
	of the index you left off at which you can find out by running the code and
	then looking at the "pairs.txt" in the catalog directory. If you haven't
	checked a spectrum, it will have a -1 in the "QUALITY_FG" column.
	
Task 2: Repeat the process for background quasar spectra using
	"python checkBackground.py"
	
	Everything else is the same. Except now the quality flag is "QUALITY_BG" in
	the pair table.
	
	Going through another ~12000
	
	