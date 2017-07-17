#!/usr/bin/env python3
#
# Generate the needed frequency values which are required to
# calculate the spectra of ICA extracted cluster.
#
# The extracted spectrum has data points of 65, 80, 95, ..., 185 MHz.
#
# However, ICA separation requires in total 5 input images of frequencies
# (f-4, f-2, f, f+2, f+4) MHz, and each images has frequency band width
# of 1 MHz, which combined from (f-0.5, f-0.25, f, f+0.25, f+0.5) MHz images.
#
# Weitian LI
# 2015/03/26
#

begin = 65
end = 185
step = 15
sep1 = 2 # separting frequency between each input images for ICA
bandwidth = 1
sep2 = 0.25

def gen_spec_freqs():
    # To include the end point
    freqs = []
    for i in range(begin, end+1, step):
        for j in [i-2*sep1, i-sep1, i, i+sep1, i+2*sep1]:
            for f in [j-2*sep2, j-sep2, j, j+sep2, j+2*sep2]:
                freqs.append(f)
    return freqs


def main():
    freqs = gen_spec_freqs()
    for f in freqs:
        print("%06.2f" % f) # format: xxx.xx


if __name__ == "__main__":
    main()

