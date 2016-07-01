import unittest
import ga.ga as ga
import os.path
import os
import numpy as np
from scipy.spatial.distance import euclidean

class TestMixerMethods(unittest.TestCase):
	# def test_something(self):
	# 	self.assertEquals(True, True)

	# def test_get_population(self):
	# 	pop = ga.get_population(3, ga.CONST_DX_PARAM_COUNT)
	# 	self.assertEquals(len(pop), 3)
	# 	self.assertEquals(len(pop[0]), ga.CONST_DX_PARAM_COUNT)

	# def test_render_population(self):
	# 	pop = ga.get_population(3, ga.CONST_DX_PARAM_COUNT)
	# 	files = ga.render_population(pop, "./")
	# 	for f in files:
	# 		self.assertTrue(os.path.isfile(f))
	# 	# now clean up
	# 	for f in files:
	# 		os.unlink(f)


	# def test_zero_nans(self):
	# 	a = [np.inf, np.nan, -np.inf, 24]
	# 	a = ga.zero_nans(a)
	# 	self.assertEquals(len(a), 4)
	# 	self.assertEquals(a[0], 0)
	# 	self.assertEquals(a[1], 0)
	# 	self.assertEquals(a[2], 0)
	# 	self.assertEquals(a[3], 24)

	# def test_population_to_mfccs(self):
	# 	pop = ga.get_population(3, ga.CONST_DX_PARAM_COUNT)
	# 	mfccs = ga.population_to_mfccs(pop)
	# 	self.assertEquals(len(pop), len(mfccs))
	# 	# this mfcc library generates frames of 13 coeffs
	# 	self.assertEquals(len(mfccs[0][0]), 13)
	# 	# check there are no nans anywhere


	# def test_compute_fitness(self):
	# 	a1 = [[1,2,3]]
	# 	a2 = [[1,2,3]]
	# 	f = ga.compute_fitness(a1, a2)
	# 	self.assertEquals(f, 1)
	# 	a1 = [[1,2,3]]
	# 	a2 = [[0,2,3]]
	# 	f = ga.compute_fitness(a1, a2)
	# 	self.assertEquals(f, 1/ euclidean(a1[0], a2[0]))
	# 	# test the different length code
	# 	a1 = [[1,2,3]]
	# 	a2 = [[0,2,3], [1,2,3]]
	# 	# should ignore the second frame in a2..
	# 	self.assertEquals(f, 1/ euclidean(a1[0], a2[0]))

	# def test_compute_population_fitnesses(self):
	# 	target = ga.wav_to_mfcc("test.wav")
	# 	pop = ga.get_population(3, ga.CONST_DX_PARAM_COUNT)
	# 	results = ga.compute_population_fitnesses(target, pop)
	# 	print results
	# def test_crossover(self):
	# 	p1 = [0.1, 0.1, 0.1, 0.1]
	# 	p2 = [0.2, 0.2, 0.2, 0.2]
	# 	p3 = ga.crossover(p1, p2, 5)
	# 	self.assertEquals(len(p3), len(p1))

	# def test_select_index(self):
	# 	fs = [0.1, 0.1, 0.4, 0.1, 0.6]
	# 	ind = ga.select_index(fitnesses = fs)
	# 	self.assertTrue(ind < len(fs))

	def test_mutate(self):
		p1 = [0.1, 0.1, 0.1]
		for i in range(0, 100):
			p1 = ga.mutate(p1, rate = 1.0)
			self.assertTrue(max(p1) <= 1.0)
			self.assertTrue(min(p1) >= 0.0)


	def test_ga_optimise(self):
		# simple population processer function 
		# that just assigns fitness based on the sum of
		# a set of parameters
		def cpf(target, pop):
			print "cpf called!"
			fs_and_ps = []
			for ind in range(0, len(pop)):
				# fitness is a sum of the parameter settings
				# so higher parameter values should be selected for...
				f_and_p = {"fitness":sum(pop[ind]), "params":pop[ind]}
				fs_and_ps.append(f_and_p)
			return fs_and_ps

		iterations = 10
		fs = ga.ga_optimise(compute_population_fitnesses = cpf, 
				target = [[0.1], [0.1]], 
				param_count = ga.CONST_DX_PARAM_COUNT, 
				iterations = iterations, 
				pop_size = 1000, 
				crossovers = 3, 
				mutation_rate = 0.5	
				)
		self.assertEquals(len(fs), iterations)
		self.assertTrue(max(fs) <= ga.CONST_DX_PARAM_COUNT)
#	def test_get_mutant_pop(self):
		

if __name__ == '__main__':
    unittest.main()
