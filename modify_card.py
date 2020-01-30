import os
import sys
import numpy as np


pt_min = 1.5
pt_max = 18.0
pt_step = 0.5

pt_array = np.arange(pt_min, pt_max, pt_step)

lines = []
remove = []

with open("./ncteq.card") as card:
    text = card.readlines()

text = [line for line in text if not ("--pTmin" in line)]


pt_list = [[pt, pt+pt_step] for pt in pt_array]


for pt_pair in pt_list:
    text.append("{:4.2f}".format(pt_pair[0]) + "\t" + "{:4.2f}".format(pt_pair[1]) + "\t--pTmin pTmax\n")


with open("./test.card", "w+") as card:
    card.writelines(text)