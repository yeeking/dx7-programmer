import ga.ga as ga
import os
import datetime


def ga_optimise(synth, param_count, target, output_dir, iterations = 10, pop_size = 500):
	fs = ga.ga_optimise(compute_population_fitnesses = ga.compute_population_fitnesses, 
				target = target, 
				synth = synth, 
				param_count = param_count, 
				iterations = iterations, 
				pop_size = pop_size, 
				crossovers = param_count / 5, 
				mutation_rate = 0.5, 
				log = True, 
				data_folder = output_dir)
	return fs


if __name__ == '__main__':
	vst_synth = "../mda DX10.vst"
	vst_param_count = 15
	target_dir = "../runs/" + datetime.datetime.now().strftime("%Y%m%d%H%M%s") + "/"
	os.mkdir(target_dir)

	# first generate the target sounds
	# which are the 32 presets from the synth
	for i in range(0, 32):
		filename = target_dir + "preset_"+str(i)+".wav"
		ga.render_preset(vst_synth, i, filename)

	for i in range(0, 32):
		filename = target_dir + "preset_"+str(i)+".wav"
		print "Looking for "+filename
		target_mfccs = ga.wav_to_mfcc(filename)
		data_folder = target_dir + "_preset_"+str(i) + "/"		
		try:
			os.mkdir(data_folder)
		except:
			print "data folder already there."
		ga.string_to_file("synth: "+vst_synth + "\npreset: "+str(i), data_folder + "details.txt")
		ga_optimise(vst_synth, vst_param_count, target_mfccs, data_folder)
		

	# targets = ga.get_files_in_dir(test_dir, filter = "wav")
	# for target in targets:
	# 	print "Looking for "+target
	# 	target_mfccs = ga.wav_to_mfcc("test.wav")
	# 	data_folder = "data/data_"+target+"/"
	# 	try:
	# 		os.mkdir(data_folder)
	# 	except:
	# 		print "data folder already there."
	# 	ga_optimise(vst_synth, vst_param_count, target_mfccs, data_folder)


