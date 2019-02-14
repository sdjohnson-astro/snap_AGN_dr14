Start the measurement routine by typing "python measureMgII.py" in the terminal

If you are starting off from some previous index, say 512, type "python measureMgII.py -i 512" into the terminal.

This will show you a zoomed in region in the background quasar spectrum around
the expected location of an absorption feature at the foreground quasar
redshift. You need then to mark the location of the absorption features if
present, mark regions that are free of features to fit the continuum, and mark
the boundary of the two absorption features to enable the measurements. Then
finally you evaluate the quality of the measurement as "accepted", "rejected",
"uncertaint", or "multiple" where the last category is for systems that show 
many absorption features and need further evaluations.

Some of the key commands are as before (e.g. "n" or "b" moves you to the next or
previous index). But there are also a bunch of new ones. Here is a list of all
the commands and then a list of step-by-step instructions.

"b": move back one index (and save the catalog)
"n": move forward one index (and save the catalog)
"m": mark the left-most absorption line pair (MgII 2796). This will cause two
     yellow vertical lines to move with the left-most one where your mouse was
     when "m" was clicked. If the absorption lines are real, they should line up
     in pairs with marked by these lines.
"z": Set the left continuum region lower boundary
"x": Set the left continuum region upper boundary
"c": Set the right continuum lower boundary
"v": Set the right continuum upper boundary. After setting all the boundaries
     the code will fit the continuum and overplot the fit in green
"1": Set the lower boundary for the first absorption feature
"2": Set the upper boundary for the first absorption feature. After setting
     both the code will fit the absorption profile and over-plot the result in
     green.
"3": Set the lower boundary for the first absorption feature
"4": Set the upper boundary for the second absorption feature. After setting
     both the code will fit the absorption profile and over-plot the result in
     green.
"a": accept the continuum fit and absorption measurements
"r": reject due to problems with the spectrum or the continuum fit and  
      absorption measurements
"u": mark the absorption feature as uncertain (e.g. if you only see one line)
"l": if there are more than two absorption features. This is a flag to give it
     a closer look
     
Step-by-step instructions after starting the interface:
(1) Mark the position of the left-most absorption feature by putting your mouse
    over it and typing "m". The two vertical yellow lines will be drawn with
    the left-most one at the location of your mouse. These should mark the
    absorption lines.
(2) Mark the left lower, left upper, right lower, and right upper continuum
    regions by typing "z", "x", "c", "v" with your mouse at the desired
    locations. Each location will be marked by a thing vertical green line.
    Once all four boundaries are marked a thicker green line will be drawn
    showing the continuum fit. This should go through the data and pass
    above the absorption features.
(3) Mark the left absorption feature boundaries by typing "1" and "2". Once both
    boundaries are marked the code will fit the line and plot the fit in green.
(4) Mark the left absorption feature boundaries by typing "1" and "2". Once both
    boundaries are marked the code will fit the line and plot the fit in
    green.
(5): Mark the quality as accepted ("a"), rejected ("r"), uncertain ("u"), or 
     many lines ("l")
     
     
If you do not see the two absorption features then perform only the continuum
step.