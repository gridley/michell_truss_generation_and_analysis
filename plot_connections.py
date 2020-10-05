#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.colors
from matplotlib import cm
import sys

input_file = sys.argv[1]
pins = []
joints = []

inserting_pins = False
inserting_joints = False
inserting_connections = False

with open(input_file) as fh:
    for line in fh:
        if 'joints' in line:
            inserting_pins = False
            inserting_connections = False
            inserting_joints = True
            continue
        elif 'pins' in line:
            inserting_pins = True
            inserting_connections = False
            inserting_joints = False
            continue
        elif 'connections' in line:
            inserting_pins = False
            inserting_connections = True
            inserting_joints = False
            continue
        elif 'forces' in line:
            break

        split_line = line.split()
        if inserting_pins:
            pins.append((float(split_line[0]), float(split_line[1])))
        elif inserting_joints:
            joints.append((float(split_line[0]), float(split_line[1])))
        elif inserting_connections:
            try:
                first_char = [c for c in split_line[0] if not c.isdigit()][0]
                second_char = [c for c in split_line[1] if not c.isdigit()][0]
            except:
                print(first_char)
                print(second_char)
                raise Exception("^^^")
            first_num = int(''.join([c for c in split_line[0] if c.isdigit()]))
            second_num = int(''.join([c for c in split_line[1] if c.isdigit()]))
            if first_char == 'j':
                point1 = joints[first_num]
            elif first_char == 'p':
                point1 = pins[first_num]
            else:
                raise Exception('hmm? {}'.format(first_char))
            if second_char == 'j':
                point2 = joints[second_num]
            elif second_char == 'p':
                point2 = pins[second_num]
            else:
                raise Exception('hmm? {}'.format(second_char))
            plt.plot((point1[0], point2[0]), (point1[1], point2[1]))
        else:
            pass

plt.show()
