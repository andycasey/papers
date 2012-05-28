#!/usr/bin/python
# *-* coding:utf-8 -*-
# vim: set fileencoding=utf-8 :
"""
Generalized histogram python code

Andy Casey - June 2010

Create a generalized histogram given a TOPCAT ASCII file

	
Inputs:

	None.
	
	Simply run 'python gen-histogram.py' or your python executable weapon of choice
	then it will ask you what's needed.
	
	Tab completion will intuitively show you the expressions you can use.

Outputs:

	kernel.png
	cdf.png
	

Notes:

	Use typical python-like (engrish) expression for the logic expression.
	
	Examples of logic expression:
	
	1. 	(( EWMgH > 0 and sensible ) or ( EWMgH > 0 and HGHT1fx >= 0.7)) and interpMgH != 0
	
	2. 	my-giant-subset or my-giant-subset-2
	
	3. 	5 > EWMgH >= 4
	
	If your expression is invalid, it will tell you - it cannot be broken. If you give
	it an expression that returns zero data items, it will completely ignore your
	expression and use all rows in the file.
	
	Using logic expressions like this WON'T work:
	
	1. 	EWMgH>0 and l>32
	
	2. 	VGSR1 > 0 && EWCaT > 0
	
	Because the header items need to be words (separated by spaces), so 'EWMgH>0' won't
	be recognized but 'EWMgH >0' is okay. And the ||, && logic isn't python-like, so
	you'll have to use and/or instead.

"""


import os, glob, readline, atexit
import random as rand
from itertools import count, izip
from numpy import *
from matplotlib.pyplot import *
from pylab import *
#from matplotlib import rcParams
#matplotlib.rcParams['text.usetex'] = True

options = []

# Tab completion
def complete(text, state):
	for option in options:
		if option.startswith(text):
			if not state:
				return option
			else:
				state -= 1

readline.parse_and_bind('tab: complete')
readline.set_completer(complete)
#' \t\n`~!@#$%^&*()-=+[{]}\\|;:\'",<>/?'
readline.set_completer_delims(' \t\n`~!@#$%^&*=+[{]}\\|;:\'"<>/?')
histfile = os.path.join(os.environ['HOME'], '.pyhist')
try:
	readline.read_history_file(histfile)
except IOError:
	pass
atexit.register(readline.write_history_file, histfile)


# Grab the data file we wish to use

fileName = ''

while not fileName:
	options = glob.glob('*')
	fileName = raw_input('Data file: ').strip()

	if not fileName:
		exit()
	
	if not os.path.exists(fileName):
		print 'File not found'
		fileName = ''
		break
	
	try:
		data = open(fileName, 'r')
	except Exception:
		fileName = ''

	# File exists, and is open. Check for #headers
	headers = data.readline()
	lenHeaders = headers
	headers = headers.split()

	if ((headers[0] != '#') or (3 > len(headers))):
		print 'No headers found in ' + fileName + ',.. must be a valid TOPCAT ASCII file'
		fileName = ''

	else:
		# Remove the beginning '#'	
		headers.pop(0)
	
# We have the file open in data.* and we have the headers
options = headers

histoVar, histoErr = ['', '']

while not histoVar:
	histoVar = raw_input('Variable header of interest: ').strip()

	if not histoVar:
		exit()
	
	try:
		varIdx = headers.index(histoVar)
	except Exception:
		print 'Header not found'
		histoVar = ''
			

# Get the errors associated with that variable

while not histoErr:
	histoErr = raw_input('Variable error header of interest: ').strip()
	
	if not histoErr:
		exit()
		
	try:
		errIdx = headers.index(histoErr)
	except Exception:
		print 'Header not found'
		histoErr = ''


# Ask for any additional logic required

logic = ''
while not logic:
	logic = raw_input('Logical expression: ').strip()

	# Logic isn't required	
	if not logic:
		break

	# If it's given though..
	try:
		# Use a test line
			
		testLine = data.readline().split()
		data.seek(len(lenHeaders))
		
		# Create the expression
		logicActual = ''
		
		for logicExp in logic.split():
			try:
				idx = headers.index(logicExp)
				term = testLine[idx]
				
				# For 'true' and 'false' they must be 'True' and 'False' for python to handle them well
				if (type(term) == type(str())):
					term = term.replace('true', 'True')
					term = term.replace('false', 'False')
					
				logicActual += ' ' + term
			except Exception:
				logicActual += ' ' + logicExp
		
		# logicActual is our expression
		
		if eval(logicActual):
			readline.add_history(logic)
		else:
			readline.add_history(logic)
	
	except Exception:
		print "Invalid expression"		
		logic = ''
		
# Bandwidth
bandwidth = ''

while not bandwidth:
	bandwidth = raw_input('Kernel bandwidth calculation (e.g. sample, all, or any float; default is all): ')
	
	if not bandwidth or (bandwidth == 'all'):
		bandwidth = 'histoErrsAll'
		break
		
	elif (bandwidth == 'sample'):
		bandwidth = 'histoErrs'
		break
			
	else:
		try:
			bandwidth = float(bandwidth)
		
		except Exception:
			print 'Kernel bandwidth must be "sample", "all", or some float'
			bandwidth = ''


# Select a comparison function
observedHeaders = options

options =	[
						'halo',
						'sgr-prolate',
						'sgr-oblate',
						'sgr-triaxial',
						'sgr-spherical'
					]
					
comparisonFunction = ''

while not comparisonFunction:
	comparisonFunction = raw_input('Comparison function (e.g. halo, sgr-prolate): ')
	
	if (comparisonFunction == 'halo'):
		# good
		break
	
	comparisonFiles = {
											'sgr-prolate' 	: 'prolate.dat',
											'sgr-oblate'		:	'oblate.dat',
											'sgr-triaxial'	:	'SgrTriax_DYN.dat',
											'sgr-spherical'	:	'spherical.dat'
										}
	
	try:
		comparisonFile = comparisonFiles[comparisonFunction]
			
	except Exception:
		print 'Unknown comparison function, double-tab for options'
		comparisonFunction = ''

	if comparisonFunction:
		if not os.path.exists(comparisonFile):
			print 'Cannot find file ' + comparisonFile + ' for ' + comparisonFunction + ' comparison function'
			comparisonFunction = ''
	
		else:
				
			# Get the model column headers
	
			model = open(comparisonFile, 'r')
			modelHeaders = model.readline().split()
		
			if ((modelHeaders[0] != '#') or (3 > len(modelHeaders))):
				print 'No headers found in ' + comparisonFile + ',.. must be a valid TOPCAT ASCII file'
				comparisonFunction = ''
			
			else:			
				# Remove the beginning '#'
				modelHeaders.pop(0)
		
	

#COMPARISON FUNCTION IS A HALO GAUSSIAN

if (comparisonFunction == 'halo'):

	options = []

	# Get a sample number

	sampleNumber = ''

	while not sampleNumber:

		sampleNumber = raw_input('Sample size number: ')
	
		if not sampleNumber:
			sampleNumber = 180
			break
		
		else:
			try:
				sampleNumber = int(sampleNumber)
			except Exception:
				print 'Sample size must be a positive integer'
				sampleNumber = ''
		
			if (1 > sampleNumber):
				print 'Sample size must be a positive integer'
				sampleNumber = ''

			
	# Mean value
	
	gaussianMean = ''
	defaultGaussianMean = 0.00 # Sirko et al 2004

	while not gaussianMean:
		#σμ
		gaussianMean = raw_input('Gaussian mean distribution μ (e.g. Vgsr,μ = 0.0 km/s Sirko et al. 2004): ')
	
		if not gaussianMean and (gaussianMean != 0):
			gaussianMean = defaultGaussianMean
			break
	
		else:
			try:
				gaussianMean = float(gaussianMean)
			except Exception:
				print 'Mean distribution must be a real float'
				gaussianMean = ''
		

	# Sigma value

	gaussianSigma = ''
	defaultGaussianSigma = 100.0 # Sirko et al 2004

	while not gaussianSigma:

		gaussianSigma = raw_input('Gaussian standard deviation σ (e.g. Vgsr,σ = 100.0 km/s Sirko et al. 2004): ')
	
		if not gaussianSigma and (gaussianSigma != 0):
			gaussianSigma = defaultGaussianSigma
			break
		
		else:
			try:
				gaussianSigma = float(gaussianSigma)
			except Exception:
				print 'Standard deviation must be a real float'
				gaussianSigma = ''
				
else:
	
	#COMPARISON FUNCTION IS SOME SGR MODEL
	
	# We will need to match up positions in the observed and model data
	
	# Ask for position limits from observedHeaders options
	# Ask for the equivalent limiting column headers
	
	if logic:
		print '	Note: Logical expression will not be applied to model points'
	
	options = observedHeaders
	
	obsPositionLimits = ''
	
	while not obsPositionLimits:
		
		obsPositionLimits = raw_input('Position column header(s) from observed data (e.g. "lambdaSgr" or "l, b"): ')
		
		if (type(obsPositionLimits) != type(str())):
			obsPositionLimits = ''
			break
			
		if obsPositionLimits:
			obsPositionLimits = obsPositionLimits.split(',')
						
			for i in obsPositionLimits:
				i = i.strip()
				
				try:
					j = observedHeaders.index(i)
				
				except Exception:
					print 'The option ' + i + ' does not appear to be a column header'
					obsPositionLimits = ''
					break

			if obsPositionLimits:
				if (len(obsPositionLimits) == 1):
					lSgrObsVar = obsPositionLimits[0].strip()
					lSgrObsIdx = observedHeaders.index(lSgrObsVar)
				
				elif (len(obsPositionLimits) == 2):
					xObsVar = obsPositionLimits[0].strip()
					yObsVar = obsPositionLimits[1].strip()
			
					xObsIdx = observedHeaders.index(xObsVar)
					yObsIdx = observedHeaders.index(yObsVar)
				
				else:
					print "Don't know how to handle more than 2 position inputs"
					obsPositionLimits = ''
			
	
		# Check to see if the x and y columns exist in the model headers
		
		try:
			if (len(obsPositionLimits) == 1):
				lSgrModelVar = lSgrObsVar
				lSgrModelIdx = modelHeaders.index(lSgrModelVar)
			
			else:
				xModelVar = xObsVar
				yModelVar = yObsVar
		
				xModelIdx = modelHeaders.index(xModelVar)
				yModelIdx = modelHeaders.index(yModelVar)
				
		except Exception:
	
			options = modelHeaders
			modelPositionLimits = ''
		
			while not modelPositionLimits:
				modelPositionLimits = raw_input('Equivalent position column header(s) from model data (e.g. "lambdaSgr" or "ra, dec"): ')
			
				if (type(modelPositionLimits) != type(str())):
					modelPositionLimits = ''
					break
			
				if modelPositionLimits:
					modelPositionLimits = modelPositionLimits.split(',')
					
					for i in modelPositionLimits:
						i = i.strip()
					
						try:
							j = modelHeaders.index(i)
						
						except Exception:
							print 'The option ' + i + ' does not appear to be a column header'
							modelPositionLimits = ''
							break
				
					if modelPositionLimits:
						if (len(modelPositionLimits) != len(obsPositionLimits)):
							print 'Degrees of model position limits and observed position limits does not match'
							modelPositionLimits = ''
							
						else:
							if (len(modelPositionLimits) == 1):
								lSgrModelVar = modelPositionLimits[0].strip()
								lSgrModelIdx = modelHeaders.index(lSgrModelVar)
						
							else:
								xModelVar = modelPositionLimits[0].strip()
								yModelVar = modelPositionLimits[1].strip()
				
								xModelIdx = modelHeaders.index(xModelVar)
								yModelIdx = modelHeaders.index(yModelVar)
						
				

			# Get our model var index
		modelVar = histoVar

		try:
			modelIdx = modelHeaders.index(modelVar)

		except Exception:
			print '	Could not find "' + histoVar + '" in ' + comparisonFile + ' headers'

			modelVar = ''	
			options = modelHeaders

			while not modelVar:
				modelVar = raw_input('Equivalent model variable header of interest: ')
	
				if (type(modelVar) == type(str())):
					try:
						modelVar = modelVar.strip()
						modelIdx = modelHeaders.index(modelVar)
		
					except Exception:
						print 'Header not found'
						modelVar = ''
		
				else:
					modelVar = ''

#INPUTS OVER

# READ IN EVERYTHING WE NEED
# Key variables: filename, histoVar, histoErr, varIdx, errIdx

histoVars = []
histoErrs = []
histoErrsAll = []
#comparisonFunction = 'halo'
"""
if (comparisonFunction != 'halo'):
	xPosLimits = [999, -999] # some dummy values that are [high, low] for [min, max] overwriting
	yPosLimits = xPosLimits
	lSgrPosLimits = xPosLimits
"""
while not len(histoVars):
	# Load up our variables and errors

	for line in data:
		logicActual = ''
		line = line.split()
	
		# Check for logic
		if logic:
			for logicExp in logic.split():
				try:
					idx = headers.index(logicExp)
					term = line[idx]
				
					# For 'true' and 'false' they must be 'True' and 'False' for python to handle them well
					if (type(term) == type(str())):
						term = term.replace('true', 'True')
						term = term.replace('false', 'False')
					
					logicActual += ' ' + term
				except Exception:
					logicActual += ' ' + logicExp
		
			histoErrsAll.append(float(line[errIdx]))		

			if (comparisonFunction != 'halo'):
				# Lambda limits
				
				if (len(obsPositionLimits) == 1):
					lSgrPosLimits = [min([lSgrPosLimits[0], float(line[lSgrObsIdx])]),
													 max([lSgrPosLimits[1], float(line[lSgrObsIdx])])]
			
				else:
					# Position limits
					
					# Minimum and maximum 'x' values
					xPosLimits = [min([xPosLimits[0], float(line[xObsIdx])]),
												max([xPosLimits[1], float(line[xObsIdx])])]
								
					# Minimum and maximum 'y' values
					yPosLimits = [min([yPosLimits[0], float(line[yObsIdx])]),
												max([yPosLimits[1], float(line[yObsIdx])])]
			
														
			# logicActual is our expression
			if eval(logicActual):
				histoVars.append(float(line[varIdx]))
				histoErrs.append(float(line[errIdx]))

			else:
				continue
			
		else:
			histoVars.append(float(line[varIdx]))
			histoErrs.append(float(line[errIdx]))
			histoErrsAll.append(float(line[errIdx]))
			
			if comparisonFile:
				# Lambda limits
				if (len(obsPositionLimits) == 1):
					lSgrPosLimits = [min([lSgrPosLimits[0], float(line[lSgrObsIdx])]),
										  		 max([lSgrPosLimits[1], float(line[lSgrObsIdx])])]
			
				else:
					# Position limits
					
					# Minimum and maximum 'x' values
					xPosLimits = [min([xPosLimits[0], float(line[xObsIdx])]),
												max([xPosLimits[1], float(line[xObsIdx])])]
								
					# Minimum and maximum 'y' values
					yPosLimits = [min([yPosLimits[0], float(line[yObsIdx])]),
												max([yPosLimits[1], float(line[yObsIdx])])]
			
								
								
	if not len(histoVars):
		print "Ignoring the logical expression as it returned zero data rows"
		logic = ''
		data = open(fileName, 'r')
		data.readline()

# Generate bandwidth

if (type(bandwidth) == type(str())):
	if (bandwidth == 'histoErrs'):
		bandwidth = mean(histoErrs)
		
	elif (bandwidth == 'histoErrsAll'):
		bandwidth = mean(histoErrsAll)
	
	else:
		print 'WARNING: Bandwidth was not calculated properly, using sample selection bandwidth'
		
		q = raw_input('Push enter to continue')
		bandwidth = mean(histoErrs)

if (comparisonFunction != 'halo'):

	# COMPARISON FILE
	
	# By now we have xObsIdx, yObsIdx, xModelIdx, yModelIdx, xPosLimits
	# We can generate a sample from our model data

	comparVars = []
	for line in model:
	
		line = line.split()
		
		# Check if this is within our xPosLimits, yPosLimits
		if (len(obsPositionLimits) == 1):
			if (lSgrPosLimits[1] > float(line[lSgrModelIdx]) > lSgrPosLimits[0]):
				comparVars.append(float(line[modelIdx]))
		else:
			if ((xPosLimits[1] > float(line[yModelIdx]) > xPosLimits[0]) and (yPosLimits[1] > float(line[xModelIdx]) > yPosLimits[0])): # wtf?
				comparVars.append(float(line[modelIdx]))



# It's GISTOGRAM time
# Create our x-variable range
if (comparisonFunction != 'halo'):
	totalVars = comparVars + histoVars
else:
	totalVars = histoVars

totalVars = histoVars
minVar, minIdx = min(izip(totalVars, count()))
maxVar, maxIdx = max(izip(totalVars, count()))

x = arange(minVar -2, maxVar +1, 0.01)
y = zeros(len(x))

# CLEAR!
numOnScreen = 50.0
#os.system('clear')

# Generate the generalized histogram
i, j = [0.0, len(histoVars)]
print bandwidth
for mean, sigma in zip(histoVars, histoErrs):
	y += exp((-(x-mean)**2)/(2*bandwidth**2))/sqrt(2*pi)
	i += 1
	
	completed = int(floor(100*i/j)/(100/numOnScreen))
	completedStr = "#"*completed
	incompleteStr = "-"*int(numOnScreen-completed)
	
#	os.system('clear')
#	print "Generating observed functions\n["+completedStr+incompleteStr+"]	("+str(round(100*i/j, 2))+"%)"


# Normalize and smooth the kernel


# COMPARISON FUNCTIONS


if (comparisonFunction == 'halo'):
	# Generate the Gaussian distribution function
	g = exp(-((x-gaussianMean)**2)/(2*gaussianSigma**2))*(1/sqrt(2*pi*gaussianSigma**2))

	# Generate the N sigma lines
	sigmas = []
	#sigmas.append(1) # 67% confidence level
	#sigmas.append(2) # 95% confidence level
	sigmas.append(3) # 99.7% confidence level

	#gpSigmas = [(g*sampleNumber + k*sqrt(g*sampleNumber**2))/sampleNumber**2 for k in sigmas]
	gpSigmas = [g + k*sqrt(g)/sampleNumber for k in sigmas]
	gmSigmas = [g - k*sqrt(g)/sampleNumber for k in sigmas]

	# Some minus (unrealistic) values may result, let's minimize the floor to zero
	gmSigmas = [[max(0, gmSigmas[i][j]) for j in range(0, len(gmSigmas[i]))] for i in range(0, len(gmSigmas))]

	# Normalize and smooth the kernel
	y = y/(len(histoVars)*bandwidth)
	
else:	
	z = zeros(len(x))
	i, j = [0.0, len(comparVars)]
	
	for mean in comparVars:
	
		z += exp((-(x-mean)**2)/(2*bandwidth**2))/sqrt(2*pi)
		i += 1
		
		completed = int(floor(100*i/j)/(100/numOnScreen))
		completedStr = "#"*completed
		incompleteStr = "-"*int(numOnScreen-completed)
	
		os.system('clear')
		print "Generating comparison functions\n["+completedStr+incompleteStr+"]	("+str(round(100*i/j, 2))+"%)"
		
	# Normalize and smooth the kernels
	#y = y/((len(histoVars) + len(comparVars))*bandwidth)	
	#z = z/((len(histoVars) + len(comparVars))*bandwidth)
	y = y/(len(histoVars)*bandwidth)
	z = z/(len(comparVars)*bandwidth)

	


# Generate the cumulative distribution function
binLimits = range(int(floor(min(x))) - 1, int(ceil(max(x))) + 2)
cdf = zeros(len(binLimits))

# Bin it up into steps of 1
for (binNum, cdfValue) in zip(digitize(x, binLimits), y): 
    cdf[binNum] += cdfValue

# Get the CDF
cdf = [cdf[i] + sum(cdf[0:i]) for i in range(0, len(cdf))]


# Write the results to output files
#file = open('kernel.dat', 'w')
#file.writeline("#x	y")
#for xa, ya in x, y:
#	file.writeline(str(x)+"	"+str(y)+"\n")
#file.close()

# Plot the result

clf()
plot(x, y, 'r-')
xlabel(r'$V_{los,GC}$')
ylabel(r'$\sum{}P(V_{los,GC})$')

if (comparisonFunction == 'halo'):
#	plot(x, g, 'b-')
	
	# Calculate chi-squared
	dataFile = open(histoVar + '-' + comparisonFunction + '.txt', 'w')
	for i in range(0, len(g)):
		dataFile.writelines(str(x[i])+' '+str(g[i])+' '+str(y[i])+' '+str(gmSigmas[0][i])+' '+str(gpSigmas[0][i])+"\n")
	dataFile.close()
		
	
	chi = sum((y - g)**2)/len(y)
	title(r'$n_{O} = ' + str(len(histoVars)) +', n_{G} = '+ str(sampleNumber) + ', \sum\chi^{2} = ' + str(chi) + ', L = ' + logic + '$')
	
#	for gpSigma in gpSigmas:
#		plot(x, gpSigma, 'b-.')

#	for gmSigma in gmSigmas:
#		plot(x, gmSigma, 'b-.')

else:
	plot(x, z, 'b-')

	# Save data to txt
	
	dataFile = open(histoVar + '-' + comparisonFunction + '.txt', 'w')
	for i in range(0, len(z)):
		dataFile.writelines(str(x[i])+' '+str(z[i])+' '+str(y[i])+"\n")
	dataFile.close()
	
	# Calculate chi-squared
	
	chi = sum((y - z)**2)/len(y)
	title(r'$n_{O} = ' + str(len(histoVars)) + ', n_{M} = ' + str(len(comparVars)) + ', \sum\chi^{2} = ' + str(chi) + ', L = ' + logic + '$')
	

plotFileName = histoVar + '-' + comparisonFunction + '.png'
savefig(plotFileName)


clf()
plot(binLimits, cdf, 'b.')
xlabel(r'$V_{los,GC}$')
cdfFileName = 'cdf-' + histoVar + '-' + comparisonFunction + '.png'
savefig(cdfFileName)

print 'Generated ' + plotFileName + ' and ' + cdfFileName
print "	" + str(len(histoVars)) + " observed values used"
if (comparisonFunction == 'halo'):
	print "	" + str(sampleNumber) + " sample number used for gaussian with;"
	print "		μ : " + str(gaussianMean)
	print "		σ : " + str(gaussianSigma)
	if (len(sigmas) > 1):
		print '	Displaying ' + '-, '.join(map(str, sigmas)) + '-sigma lines (+/-)'
	else:
		print '	Displaying ' + str(sigmas[0]) + '-sigma line (+/-)'
		
else:
	print "	" + str(len(comparVars)) + " model values used from the bounding box;"
	
	if (len(obsPositionLimits) == 1):
		print '		' + lSgrObsVar + ' : [' + str(round(lSgrPosLimits[0], 2)) + ', ' + str(round(lSgrPosLimits[1], 2)) + ']'
		print '		Translations: ' + lSgrObsVar + ',observed = ' + lSgrModelVar + ',model'
		
	else:
		print "		" + xObsVar + " : [" + str(round(xPosLimits[0], 2)) + ", " + str(round(xPosLimits[1], 2)) + "]"
		print "		" + yObsVar + " : [" + str(round(yPosLimits[0], 2)) + ", " + str(round(yPosLimits[1], 2)) + "]"
		print "		Translations: " + xObsVar + ",observed = " + xModelVar + ",model  and  " + yObsVar + ",observed = " + yModelVar + ",model"

print "\n	Chi-squared fit: " + str(chi)
	
	
	

