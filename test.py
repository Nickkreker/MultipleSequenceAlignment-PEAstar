import utils
import numpy as np

pam250 = utils.load_pam250_matrix()
alignment = np.array([['A', 'T'], ['_', 'A']])
utils.calc_alignment_score(alignment, pam250)