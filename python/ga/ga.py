import random
import numpy as np
import scipy.io.wavfile
from scikits.talkbox.features import mfcc
import os
from scipy.spatial.distance import euclidean
from os import listdir
from os.path import isfile, join, basename, splitext
import os

## todo - put these settings in an ini file
#CONST_DX_PARAM_COUNT = 147
CONST_DX_PARAM_COUNT = 15
CONST_DIR = "../"
#CONST_DX_VST = CONST_DIR + "Dexed.vst"
CONST_DX_VST = CONST_DIR + "mda DX10.vst"
CONST_MIDI_FILE = CONST_DIR + "midi_export.mid"
CONST_MRS_WATSON = "\"" +  CONST_DIR + "MrsWatson-0.9.8/Mac OS X/mrswatson" + "\""
#CONST_OUT_FILE = CONST_DIR + "output.wav"
CONST_OUT_FILE = "./output.wav"

# renders random sounds until it gets wone with mean power per frame  > min_power
def get_loud_random_sound(filename = "loud.wav", params_file = "loud_params.txt", min_power = 500000.0):
	power = 0
	while power < min_power:
		params = render_random_dx_sound(filename)
		spec = wav_to_spec(filename)
		power = get_weighted_power_per_frame(spec)
	print "Achieved power of "+str(power)
	# now write the params
	string_to_file(str(params), params_file)

# renders several random sounds and returns the mean power per frame 
# that splits the population on ratio
# e.g. if ratio = 0.75, you get the power of the 75th percentile 
def get_calibrated_power(filename, ratio = 0.95, samples = 100):
	powers = []
	for sample in range(0, samples):
		render_random_dx_sound(filename)
		spec = wav_to_spec(filename)
		#print spec
		power = get_weighted_power_per_frame(spec)
		powers.append(power)
	powers.sort()
	print powers
	ind = int(ratio * len(powers))
	return powers[ind]

def get_power_per_frame(spec):
	power = 0
	for frame in spec:
		power = power + sum(frame)
	return power/len(spec)

def get_weighted_power_per_frame(spec):
	power = 0
	# scale the power by a linear downward ramp
	# -> high frequencies less favoured
	bins = len(spec[0])
	# inverse linear ramp
	ramp = np.array([(bins - float(i)) / bins  for i in range(0, bins)])
	for i in range(0, len(spec)):
		frame = spec[0]
		power = power + sum(frame * ramp)
	return power/len(spec)		

# compute_population_fitnesses should be a function
# that takes two arguments: the target's feature vector and a 2D array representing 
# several sets of parameter settings
# this function should return a list of objects:
# [{"fitness":123, "params":[1, 2, 3, 4]}, {...}]
def ga_optimise(compute_population_fitnesses, 
				synth = CONST_DX_VST, 
				target = [[0.1], [0.1]], 
				param_count = CONST_DX_PARAM_COUNT, 
				iterations = 3, 
				pop_size = 10, 
				crossovers = 3, 
				mutation_rate = 0.1, 
				log = True, 
				data_folder = "./data/"
				):
	fs = []
	pop = get_population(pop_size, param_count)
	string_to_file("iteration, fitness, params\n", data_folder+"log.txt", append = False)
	reps = 0
	for i in range(0, iterations):
		print i
		fs_and_params = compute_population_fitnesses(target, synth, pop)
		fmax = get_fittest(fs_and_params)
		print "Iteration "+str(i)+" fittest: "+str(fmax["fitness"])
		fs.append(fmax)
		# render out a wav file
		render_individual(fmax['params'], synth, data_folder + "best_"+str(i)+".wav")
		string_to_file(str(i) + ","+str(fmax["fitness"])+",\""+str(fmax["params"]) + "\"\n", data_folder+"log.txt", append = True)
		# now check if we have converged
		if len(fs) > 1:
			if fs[len(fs) - 2]["fitness"] == fmax["fitness"]:
				reps = reps + 1
			else:
				reps = 0
		if reps == 5:# die after 5 repeated fitnesses.
			break
		pop = generate_next_population(fs_and_params, crossovers, mutation_rate)
	return fs

def get_fittest(fs_and_params):
	fmax = fs_and_params[0]
	for f in fs_and_params:
		if f['fitness'] > fmax['fitness']:
			fmax = f
	return fmax

def string_to_file(s, file, append = False):
	if append:
		text_file = open(file, "a")
	else:
		text_file = open(file, "w")
	text_file.write(s.encode('utf-8'))
	text_file.close()

# generates the next population by selecting parents 
# and cross breeding
def generate_next_population(fs_and_params = [{"fitness":0.1, "params":[0.1, 0,1]}, {"fitness":0.01, "params":[0.2, 0,1]}], 
							 crossovers = 3, mutation_rate=0.1):
	new_pop = []
	# keep the fittest one...
	new_pop.append(get_fittest(fs_and_params)['params'])
	for i in range(1, len(fs_and_params)):
		fs = [item['fitness'] for item in fs_and_params]
		ind1 = select_index(fs)
		ind2 = select_index(fs)
		params = crossover(fs_and_params[ind1]['params'], 
							fs_and_params[ind2]['params'], crossovers)
		params = mutate(params, mutation_rate)
		# then set some values to random 
		if random.random > mutation_rate:
			for i in range(0, int(len(params) * mutation_rate)):
				ind = int(len(params) * random.random())
				params[ind] = random.random()
		new_pop.append(params)
	return new_pop

# selects an index from the sent list, treating each element as a probability of being chosen, proportionate 
# to the sum of the list elements
# in other words it samples from a probaility distribution
def select_index(fitnesses = [0.1, 0.1, 0.4]):
	total = sum(fitnesses)
	target = random.random() * total
	pos = 0
	for i in range(0, len(fitnesses)):
		pos = pos + fitnesses[i]
		if target <= pos:
			ind = i
			break
	return ind

#
#   perform crossover on the sent genome, aiming to get 'crosses' crossovers.
#   
def crossover (params1, params2, crosses):
	loci = get_crossover_points(params1, params2, crosses)
	params3 = crossover_genomes_at_loci(params1, params2, loci);
	return params3
    
#   
# returns an array of count crossover points for a genome with the sent length
#  
def get_crossover_points(params1, params2, count):
	loci = []
	if len(params1) > len(params2):
	    params = params2
	else:
	    params = params1	
	if count > len(params):
	    count =  len(params) / 2
	for i in range(0, count):	
	    locus = round(random.random() * (len(params)-1))
	    # TODO might want to check we don't already have this crossover locus..
	    loci.append(locus)
	return loci

# 
# perform a crossover between the two sent genomes, crossing over at the sent loci. returns a new genome
#
def crossover_genomes_at_loci(params1, params2, loci):
    params3 = [] #we'll store the new genome here
    src = params1#begin reading from genome1
    curr_src = 0
    loci_ind = 0 #begin at the first crossover point
    if len(params1) > len(params2):
        max = len(params2)
    else:
        max = len(params1)
  	for genome_ind in range(0, max):
  		if loci_ind < len(loci) and genome_ind == loci[loci_ind]: #{// time to crossover
			curr_src = 1 - curr_src;
			if curr_src == 0:
			    src = params1
			else:
			    src = params2
			loci_ind = loci_ind + 1
		params3.append(src[genome_ind])
    return params3;

def get_mutant_pop(params, size, rate = 0.1):
	pop = []
	for i in range(0, size):
		pop.append(mutate(params, rate))
	return pop

def mutate(params, rate = 0.1):
	p_new = []
	for i in range(0, len(params)):
		change = (random.random() - 0.5) * rate
		p_new.append(min(1.0, np.abs(params[i] + change)))
	return p_new

def render_and_find_fittest(target, pop):
	fs = compute_population_fitnesses(target, pop)
	fmax = get_fittest(fx)
	return fmax

# renders the sent population of parameter settings 
# into MFCCs, then computes the fitness for each one
def compute_population_fitnesses(target, synth, pop):
	candidates = population_to_mfccs(pop, synth)
	results = []
	for i in range(0, len(pop)):
		result = {"params":pop[i], "fitness":compute_fitness(target, candidates[i])}
		results.append(result)
	return results

def compute_fitness(target, candidate):
	# first truncate the longer array
	len_diff = np.abs(len(target) - len(candidate))
	# truncate the longest one
	a1 = target
	a2 = candidate
	if len(target) > len(candidate):# truncate a1, the target
		a1 = a1[0:len(candidate)]
	if len(candidate) > len(target):# truncate a2, the candidate
		a2 = a2[0:len(target)]
	dist = 0
	for i in range(0, len(a1)):
		dist = dist + euclidean(a1[i], a2[i])
	if dist == 0:
		return 1
	else:
		return 1/dist

def population_to_mfccs(pop, synth):
	files = render_population(pop, synth, CONST_DIR)
	mfccs = []
	for f in files:
		#print "Extracting mfccs "+f
		mfccs.append(wav_to_mfcc(f))
		os.unlink(f)
	return mfccs

# sets nans and other nasties to zero
def zero_nans(data):
	def f(x):
		if (np.isnan(x) or np.isinf(x)):
			return 0
		else:
			return x
	data = map(f, data)
	return data

def render_preset(synth, index, filename):
	cmd = get_mrs_watson_preset_command(CONST_MRS_WATSON, synth, CONST_MIDI_FILE, index, filename)
#	print cmd
	os.system(cmd)	

def render_individual(params, synth, filename):
	#print "Rendering wav to "+filename
	param_string = get_mrs_watson_param_string(params)
	cmd = get_mrs_watson_command(CONST_MRS_WATSON, synth, CONST_MIDI_FILE, param_string, filename)
	os.system(cmd)

def render_population(pop, synth, output_folder):
	files = []
	for i in range(0, len(pop)):
		#print "Rendering "+str(i)
		filename = output_folder + "output_"+str(i)+".wav" 	
		param_string = get_mrs_watson_param_string(pop[i])
		cmd = get_mrs_watson_command(CONST_MRS_WATSON, synth, CONST_MIDI_FILE, param_string, filename)
		os.system(cmd)
		files.append(filename)
	return files

def get_population(size, param_count):
	pop = []
	for i in range(0, size):
		pop.append(get_random_params(param_count))
	return pop 


# utlility functions for the dx7 programmer

# get the mel spectrum
def wav_to_spec(wavfile):
	sample_rate, X = scipy.io.wavfile.read(wavfile)
	ceps, mspec, spec = mfcc(X)
	return spec

# http://stackoverflow.com/questions/25988749/mfcc-feature-descriptors-for-audio-classification-using-librosa
def wav_to_mfcc(wavfile):
	sample_rate, X = scipy.io.wavfile.read(wavfile)
	ceps, mspec, spec = mfcc(X)
	# replace nans and infs with zero
	ceps = map(zero_nans, ceps)
	return ceps
# generates a randome parameter set for the DX7 synth and renders it out to disk 
# using that parameter set.
def render_random_dx_sound(outfile):
	params = get_random_params(CONST_DX_PARAM_COUNT)
	param_string = get_mrs_watson_param_string(params)
	cmd = get_mrs_watson_command(CONST_MRS_WATSON, CONST_DX_VST, CONST_MIDI_FILE, param_string, outfile)
	os.system(cmd)
	return params
# get a command that will render the sent synth's output
def get_mrs_watson_command(mrswatson, vsti, midifile, params, outfile):
	cmd = mrswatson + " --channels 1 --quiet --plugin \""  + vsti + "\""
	cmd = cmd + " --midi-file "+midifile
	cmd = cmd + params
	cmd = cmd + " --output \"" + outfile + "\""
	return cmd

def get_mrs_watson_preset_command(mrswatson, vsti, midifile, preset_index, outfile):
	cmd = mrswatson + " --channels 1 --quiet --plugin \""  + vsti + ","+str(preset_index)+"\""
	cmd = cmd + " --midi-file "+midifile
	cmd = cmd + " --output \"" + outfile + "\""
	return cmd

# get a parameter argument string for mrswatson for the sent params
# e.g. --parameter 1,0.3 --parameter 0,0.75
def get_mrs_watson_param_string(params):
	inds = [i for i in range(0, len(params))]
	p_and_is = []
	for i in range(0, len(params)):
		p_and_is.append([inds[i], params[i]])
	ps = "".join([" --parameter "+str(p[0])+","+str(p[1]) for p in p_and_is])
	return ps

# generate a list of random values between 0 and 1 of the specified length
def get_random_params(length):
	params = np.random.rand(length)
	return params


def get_files_in_dir(dir, filter = None):
	if filter != None:
		files = [f for f in listdir(dir) if f.find(filter) != -1 and isfile(join(dir, f))]
	else:
		files = [f for f in listdir(dir) if isfile(join(dir, f))]
	return files


if __name__ == '__main__':
	for i in range(0, 5):
		render_random_dx_sound("random.wav")

