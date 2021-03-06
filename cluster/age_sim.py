#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Aaron LI
# 2015/04/12
#

"""
Simulate the evolution time (Gyr) for each cluster,
generated by a *uniformly* random distribution
between 0 and 2.7.
"""

import sys
import random


def main():
    if len(sys.argv) != 3:
        print("Usage:")
        print("    %s <number> <out: t.dat>" % sys.argv[0])
        sys.exit(-1)

    number = int(sys.argv[1])
    outfile = open(sys.argv[2], "w")

    for i in range(number):
        # randomly (uniform distribution) assign cluster evolution time
        age = random.random() / 1.0 * 2.7
        outfile.write("%.6f\n" % age)

    outfile.close()


if __name__== "__main__":
    main()

